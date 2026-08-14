"""MCMC posteriors and corner plots for the MRP fits.

Two things live here:

  Tier 0  ``corner_from_mc``   -- the existing Monte-Carlo refit draws, plotted
          jointly.  ``monte_carlo_fits`` (driver_recovery) and ``monte_carlo``
          (combined_hmf) already store all four parameters per draw; nothing has
          ever plotted their covariance.  Free, but it inherits the optimiser's
          non-convergence, so treat it as a diagnostic.

  Tier 1  ``run_emcee``        -- a proper posterior from the SAME objective.

Why the objective drops straight in.  ``make_massfn`` returns

    chi^2 = sum ((log10 phi_obs - log10 phi_mod)/sigma)^2  +  2 V sum phi_mod dlogM

and the second term is exactly the Poisson likelihood of observing ZERO groups
above the fitted range: -2 ln P(0|mu) = 2 mu, with mu = V sum phi dlogM.  So the
whole thing is already a valid -2 ln L and

    ln P = -0.5 * massfn(theta) + flat priors

needs no reinterpretation of Driver's method.  That is the point of doing it this
way rather than inventing a new likelihood: the MCMC and the published
Nelder-Mead fit are answering the same question, so the NM answer is a
regression test on the sampler.

PARAMETERISATION.  ``massfn`` takes phi LINEARLY.  We sample log10(phi) and pass
10**lp in, because a flat prior on a linear scale parameter is not scale
invariant and lets the walker wander to phi <= 0, where the model's log10 is NaN.
Reported medians are therefore of log10 phi*, which is what gets quoted anyway.

VALIDATION.  For a configuration with a genuine interior minimum (GSR with
Driver's Tempel SDSS) the best-lnP sample must reproduce the NM answer.  That
comparison is Jacobian-free, unlike comparing to the posterior median.  Where the
fit is ill-posed (GAMA-only) the two will NOT agree, and that disagreement is the
result -- see CLAUDE.md "The central finding".
"""

import numpy as np

PARAM_LABELS = [r"$\log_{10} M_*$", r"$\log_{10}\phi_*$", r"$\alpha$", r"$\beta$"]
PARAM_NAMES = ["logM*", "logphi*", "alpha", "beta"]

# Flat prior box.  The lower M* edge matters: the multi-start scan found the
# penalty term overflowing to chi^2 ~ 1e35 below logM* ~ 11, so walkers must be
# kept out of it rather than allowed to score a huge-but-finite value there.
BOUNDS = np.array([
    [11.0, 16.5],   # log10 M*
    [-8.0, 0.0],    # log10 phi*
    [-2.5, 0.5],    # alpha
    [0.05, 3.0],    # beta
])


def log_prob(theta, fn, bounds=BOUNDS):
    """ln P for emcee.  ``fn`` is a make_massfn objective taking LINEAR phi."""
    if np.any(theta < bounds[:, 0]) or np.any(theta > bounds[:, 1]):
        return -np.inf
    ms, lp, al, be = theta
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        chi2 = fn(np.array([ms, 10.0**lp, al, be]))
    # Guard the overflow region explicitly: a non-finite or absurd chi^2 must be
    # -inf, not a large finite number a walker can get stuck against.
    if not np.isfinite(chi2) or chi2 > 1e12:
        return -np.inf
    return -0.5 * chi2


def run_emcee(fn, start, nwalkers=64, nsteps=8000, burn=2000, thin=10,
              bounds=BOUNDS, seed=42, progress=True):
    """Sample the posterior.  ``start`` is (mstar, phi_LINEAR, alpha, beta), i.e.
    the Nelder-Mead answer -- walkers are scattered around it.

    Returns (chain, best, info) with chain (N,4) in log10-phi parameterisation.
    """
    import emcee

    rng = np.random.default_rng(seed)
    p0c = np.array([start[0], np.log10(abs(start[1])), start[2], start[3]])
    ndim = 4

    p0 = np.empty((nwalkers, ndim))
    scatter = np.array([0.10, 0.10, 0.05, 0.03])
    for i in range(nwalkers):
        while True:
            trial = p0c + scatter * rng.standard_normal(ndim)
            if np.isfinite(log_prob(trial, fn, bounds)):
                p0[i] = trial
                break

    # M* and phi* are anti-correlated at rho ~ -0.99 along a CURVED ridge.  The
    # default StretchMove crawls along it (tau ~ 215 at 3000 steps); the
    # differential-evolution pair is the standard fix for exactly this geometry
    # and cuts the autocorrelation by roughly an order of magnitude.
    moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(fn, bounds),
                                    moves=moves)
    sampler.run_mcmc(p0, nsteps, progress=progress)

    chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
    lnp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    best = chain[np.argmax(lnp)]

    try:
        tau = sampler.get_autocorr_time(quiet=True)
    except Exception:
        tau = np.full(ndim, np.nan)
    # n_eff per parameter: (nsteps - burn) * nwalkers / tau.  Below ~50
    # independent samples per parameter the credible intervals are not reliable.
    neff = (nsteps - burn) * nwalkers / tau
    info = dict(acceptance=float(np.mean(sampler.acceptance_fraction)),
                tau=tau, neff=neff, nsamples=len(chain),
                converged=bool(np.all(neff > 50)),
                lnp_max=float(np.max(lnp)), chi2_min=float(-2 * np.max(lnp)))
    return chain, best, info


def summarise(chain, label="posterior", nm=None):
    """Print 16/50/84 and, if given, the Nelder-Mead answer for comparison."""
    q = np.percentile(chain, [16, 50, 84], axis=0)
    print(f"\n  --- {label} ---")
    head = f"  {'param':10s}{'median':>9}{'-1sig':>8}{'+1sig':>8}"
    if nm is not None:
        head += f"{'Nelder-Mead':>13}{'(NM-med)/sig':>14}"
    print(head)
    for i, n in enumerate(PARAM_NAMES):
        lo, md, hi = q[0, i], q[1, i], q[2, i]
        sig = 0.5 * (hi - lo)
        line = f"  {n:10s}{md:9.3f}{md - lo:8.3f}{hi - md:8.3f}"
        if nm is not None:
            line += f"{nm[i]:13.3f}{(nm[i] - md) / sig if sig > 0 else np.nan:14.2f}"
        print(line)
    c = np.corrcoef(chain.T)
    print(f"  corr(logM*, logphi*) = {c[0, 1]:+.3f}   "
          f"corr(logM*, alpha) = {c[0, 2]:+.3f}")
    return q


def _to_log_phi(mc):
    """MC refit dicts store phi LINEARLY; convert to the sampler's coordinates."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.column_stack([mc["mstar"], np.log10(np.abs(mc["phistar"])),
                                mc["alphastar"], mc["betastar"]])


def corner_plot(datasets, out, truths=None, bounds=None, title=None):
    """Overlay one or more sample sets on shared corner axes.

    ``datasets`` is [(chain, label, colour), ...].  Overlaying is the whole point
    here: Tier 0 (bootstrap refits) and Tier 1 (MCMC) should agree where the fit
    converges and visibly disagree where it does not.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import corner

    if bounds is None:
        allc = np.vstack([d[0] for d in datasets])
        rng = [(np.percentile(allc[:, i], 0.5), np.percentile(allc[:, i], 99.5))
               for i in range(allc.shape[1])]
    else:
        rng = [tuple(b) for b in bounds]

    fig = None
    for chain, label, colour in datasets:
        fig = corner.corner(
            chain, labels=PARAM_LABELS, range=rng, fig=fig, color=colour,
            plot_datapoints=False, plot_density=False, fill_contours=True,
            levels=(0.393, 0.865), smooth=1.0,
            hist_kwargs=dict(density=True, lw=1.6),
            contour_kwargs=dict(linewidths=1.0),
            truths=truths, truth_color="k",
        )

    handles = [plt.Line2D([], [], color=c, lw=3, label=l)
               for _, l, c in datasets]
    if truths is not None:
        handles.append(plt.Line2D([], [], color="k", lw=1.2,
                                  label="Nelder-Mead (Driver's method)"))
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=11,
               bbox_to_anchor=(0.98, 0.98))
    if title:
        fig.suptitle(title, fontsize=13, y=1.005)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")
