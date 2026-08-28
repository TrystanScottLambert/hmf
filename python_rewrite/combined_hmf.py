#!/usr/bin/env python3
"""Python port of ``allhmf.r``: the combined multi-survey HMF.

Driver's combined figure -- GAMA + SDSS + REFLEX II fitted, with 2PIGG and
Tempel+14 shown but not fitted -- run with the **new Nessie GAMA DMU** in place
of his GAMA catalogue.  Driver's own GAMA points are overplotted so the effect of
swapping the group catalogue is visible directly.

Everything except the GAMA input is ``allhmf.r`` as written: the same
``../data`` comparison files, the same h-conversions (lines 300-334), the same
three-volume penalty, the same ``optim(maxit=500, parscale=c(1,1,1,0.1))``, the
same Monte-Carlo error loop and the same Omega_M inset.

The GAMA binned table is generated here rather than read from disk.  ``allhmf.r``
reads ``gamahmfGAMA5.csv``, which is not available, but it uses only three of its
columns -- V1 (bin centre), V4 (log10 phi_corr) and V8 (combined fractional
error) -- and those are exactly what ``driver_recovery.py`` produces and what was
verified against published table 1.  ``--write-tables`` dumps them in Driver's
column order if you want the CSVs themselves.

SDSS
----
``sdsshmf5.csv`` was never available, but ``sdsshmf.r`` is, so the SDSS leg is
rebuilt from the Tempel+14 DR10 catalogues by ``sdss_hmf.py`` and used directly.
``--sdss`` takes ``auto`` (default, Tempel+14), ``nessie`` (the Nessie SDSS
catalogue via ``nessie_sdss_hmf.py``), ``none``, or a path to a V1..V8 table if
you later get Driver's own file.

The two SDSS variants are kept as **separate figures**, not overlaid: they share
a footprint and a galaxy sample and differ only in the group finder, so the
honest comparison is two runs of the identical pipeline.  The default ``--out``
follows ``--sdss`` so they cannot overwrite each other:

    hmf_combined_nessie.pdf       GAMA = Nessie, SDSS = Tempel+14
    hmf_combined_nessie_sdss.pdf  GAMA = Nessie, SDSS = Nessie

``--sdss-max`` caps the SDSS bins entering the fit.  Nessie SDSS reaches logM
15.65 against Tempel's 15.05 and the penalty integrates over ``max(allx)``, so
the two runs do not see the same penalty range unless you set it.

Validation
----------
Run with Driver's GAMA, this port reproduces his published table 2 to within
0.035 in every parameter, for both GR and GSR -- which exercises the GAMA HMF,
the SDSS HMF, the REFLEX conversions, the combined objective, the penalty and
the optimiser together.  ``--myoption`` prints the comparison.

Usage
-----
    python combined_hmf.py                      # GAMA(Nessie) + SDSS + REFLEX
    python combined_hmf.py --omega-prior        # add Driver's Omega_M term
    python combined_hmf.py --myoption GR --sdss none
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import driver_recovery as dr
import nessie_sdss_hmf as nsd
import new_gama_hmf as ng
import sdss_hmf
from driver_recovery import (FULLSKY, G, MSOL, PARSEC, RRandom, bin_hmf, co_vol,
                             cosvar, lcdm_curve, r_optim_nm)

# ---------------------------------------------------------------------------
# Constants -- allhmf.r lines 246-268.  Note several differ from gamahmf.r.
# ---------------------------------------------------------------------------

HO = 67.37
OMEGAM = 0.3147
OMEGAMERR = 0.0074           # line 250, used by the Omega_M prior
OMEGAL = 1 - OMEGAM
RHOCRIT = 3 * (1000 * HO / (1e6 * PARSEC)) ** 2 / (8 * np.pi * G)
FITBINWID = 0.001            # line 256 -- gamahmf.r uses 0.01
LOGBIN = 0.2
ZLIMIT = 0.25
ZLIMITSDSS = 0.08            # line 262

# allhmf.r line 263 uses 175 deg^2 for GAMA where gamahmf.r line 228 uses 179.92.
# Driver's own two scripts disagree; his value is kept for his catalogue and the
# true area is used for the new one.
AREA_GAMA_DRIVER = 175.0
AREA_SDSS = 7221.0           # line 264
VOLUME_REFLEXII = 13000000.0  # line 265, hardcoded

ALPHA_FIX = -1.864908        # line 257, used by myoption="FIX"
COSVAR_REFLEX = 0.05         # line 348
DRIVER_FIT_C = (100 / 255, 149 / 255, 237 / 255)   # == ng.OLD_C, his points
OUR_FIT_C = "#d62728"        # == our red points; the fit must match its data
PUBLISHED_C = "#c8781e"      # Driver+22's *published* fit, drawn in the background

# Driver+22's headline result, quoted in his abstract and table 2 as the GSR
# row: log10(M*) = 14.13, log10(phi*) = -3.96, alpha = -1.68, beta = 0.63.
# This ONE line is drawn on every figure, whatever --myoption is set to, so the
# reference does not move from plot to plot.  (The validation printout still
# uses PUBLISHED_TABLE2[myoption], which is a different question -- there we are
# checking our port against his fit to that same sample combination.)
DRIVER_ABSTRACT_FIT = (14.13, -3.96, -1.68, 0.63)
# Driver+22 table 2 GAMA5 -- his GAMA-ONLY fit, for figure 1 only.
DRIVER_GAMA_ONLY_FIT = (13.51, -3.19, -1.27, 0.47)

# Which Driver+22 table 2 row a figure should draw behind it: the one fitted to
# the same sample combination the figure plots, not always his headline GSR.
PUB_FIT = {
    "gsr": (DRIVER_ABSTRACT_FIT, None),
    "gama": (DRIVER_GAMA_ONLY_FIT, "Driver+22 published GAMA-only fit"),
    "gs": ((14.35, -4.38, -1.96, 0.60), "Driver+22 published GAMA+SDSS fit"),
}

MLIMIT_GAMA = 12.7           # line 276
MLIMIT_SDSS = 12.9           # line 278

# Driver et al. (2022) table 2, for validation.  These are fits to *his* GAMA.
PUBLISHED_TABLE2 = {
    "G": (13.51, -3.19, -1.27, 0.47),     # GAMA5
    "S": (13.38, -3.00, -1.57, 0.47),     # SDSS5
    "R": (14.33, -4.30, -1.62, 0.79),     # REFLEX II
    "GS": (14.35, -4.38, -1.96, 0.60),
    "GR": (13.72, -3.44, -1.29, 0.55),
    "SR": (14.44, -4.52, -1.85, 0.79),
    "GSR": (14.13, -3.96, -1.68, 0.63),
}


def survey_volume(z, area):
    return area / FULLSKY * 1e9 * co_vol(np.atleast_1d(z))[0]


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def driver_table(b):
    """The eight columns ``gamahmf.r`` line 451 writes, in its reversed order.

    V1 bin centre, V2 raw N, V3 log10(weighted counts), V4 log10(phi_corr),
    V5 Poisson frac, V6 Monte-Carlo frac, V7 cosmic variance, V8 combined frac.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return pd.DataFrame({
            "V1": b.gamax[::-1],
            "V2": b.raw_counts[::-1],
            "V3": np.log10(b.wcounts)[::-1],
            "V4": np.log10(b.gamay)[::-1],
            "V5": b.rootnerr[::-1],
            "V6": b.mcerr[::-1],
            "V7": np.full(len(b.gamax), b.cosvariance),
            "V8": b.gamaf[::-1],
        })


def load_reflex(path="../data/reflex.csv"):
    """REFLEX II, allhmf.r lines 306-314.

    Note line 307 uses the *already converted* ``reflex$x`` on its right-hand
    side, because line 306 has run; the order is preserved here.
    """
    r = pd.read_csv(path)
    x = r["x"].values + np.log10(70 / HO)
    y = r["Curve1"].values + 4.0 * np.log10(HO / 70) - 14.0 + x + 1
    f = np.full(len(x), 1 / np.sqrt(20))
    f[0] = 1 / np.sqrt(3)          # R's reflexf[1]
    f[len(f) - 1] = 1 / np.sqrt(3)  # R's reflexf[43]
    return x, y, f


def load_tpigg(path="../data/tpigg.dat"):
    """2PIGG, allhmf.r lines 318-321.  Plotted, never fitted."""
    t = pd.read_csv(path, sep=r"\s+", header=None)
    return (t[0].values + np.log10(100 / HO),
            t[1].values + np.log10((HO / 100) ** 3), t[2].values, t[3].values)


def load_elmo(path="../data/elmo.csv"):
    """Tempel+14 SDSS DR10 band, allhmf.r lines 326-332.  Plotted, never fitted."""
    e = pd.read_csv(path, header=None)
    x = e[0].values + 10.0 + np.log10(100 / HO)
    s = (HO / 100) ** 3
    v2, v3, v4 = e[1].values * s, e[2].values * s, e[3].values * s
    with np.errstate(divide="ignore", invalid="ignore"):
        lo, hi = np.log10(v2 - v3), np.log10(v2 + v4)
    lo = np.nan_to_num(lo, nan=-6.5, neginf=-6.5)
    return x, np.log10(v2), lo, hi


# ---------------------------------------------------------------------------
# Objective -- allhmf.r lines 200-240
# ---------------------------------------------------------------------------

def make_massfn(allx, ally, allf, volumes, omega_prior=False, fix_alpha=None):
    """chi^2 plus the penalty, summed over every contributing survey volume.

    ``omega_prior`` adds line 238: the fitted MRP is integrated to a matter
    density and compared with OmegaM = 0.3147 +/- 0.0074.  That term, together
    with the high-mass anchor from SDSS and REFLEX, is what keeps the combined
    fit from running away to low M* the way the GAMA-only fit does.
    """
    allxxx = np.max(allx) + np.arange(1, 11) * LOGBIN
    ln10 = np.log(10.0)
    sig = allf / ln10
    xo = np.arange(0.0, 18 + 1e-9, FITBINWID)

    def fn(p):
        if fix_alpha is None:
            mstar, phi, alpha, beta = p
        else:
            mstar, phi, beta = p
            alpha = fix_alpha
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            shape = (ln10 * beta * np.exp(-10 ** (beta * (allxxx - mstar)))
                     * (phi * (10 ** allxxx / 10 ** mstar) ** (alpha + 1)))
            penalty = 2 * np.sum(shape) * LOGBIN * np.sum(volumes)
            model = np.log10(beta * ln10 * np.exp(-10 ** (beta * (allx - mstar)))
                             * (phi * (10 ** allx / 10 ** mstar) ** (alpha + 1)))
            chi2 = np.sum(((ally - model) / sig) ** 2) + penalty
            if omega_prior:
                y = (ln10 * phi * beta * (10 ** xo / 10 ** mstar) ** (alpha + 1)
                     * np.exp(-10 ** (beta * (xo - mstar))))
                om = np.sum(y * 10 ** xo) * FITBINWID * MSOL / (1e6 * PARSEC) ** 3 / RHOCRIT
                chi2 = chi2 + ((om - OMEGAM) ** 2) / OMEGAMERR ** 2
            return chi2

    return fn


def omega_matter(par, above=None):
    """Integrate a fitted MRP to a matter density (allhmf.r lines 409, 415)."""
    x = np.arange(0.0, 18 + 1e-9, FITBINWID)
    y = dr.mrp(x, *par)
    if above is not None:
        y, x = y[x > above], x[x > above]
    return np.sum(y * 10 ** x) * FITBINWID * MSOL / (1e6 * PARSEC) ** 3 / RHOCRIT


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

def assemble(sets, myoption):
    """allhmf.r lines 358-388: pick which surveys enter the fit."""
    want = {"GSR": "GSR", "Omega": "GSR", "FIX": "GSR", "GS": "GS", "GR": "GR",
            "SR": "SR", "R": "R", "G": "G"}[myoption]
    x, y, f, vols = [], [], [], []
    for key in want:
        if key not in sets:
            raise SystemExit(
                f"myoption={myoption} needs the '{key}' dataset, which is not "
                f"available.\n  SDSS in particular needs Driver's sdsshmf5.csv "
                f"(see --sdss); the raw\n  ../data/sdssdr10table*.fits are the "
                f"Tempel+14 catalogues, not a binned HMF.")
        xx, yy, ff, vv = sets[key]
        x.append(xx)
        y.append(yy)
        f.append(ff)
        vols.append(vv)
    return np.concatenate(x), np.concatenate(y), np.concatenate(f), np.array(vols)


def fit_combined(allx, ally, allf, vols, myoption, omega_prior, phimrp,
                 mstarmrp, alphamrp, betamrp, maxit=500):
    """allhmf.r lines 449-457.  parscale=c(1,1,1,0.1) here, not c(1,1,1,0.5)."""
    if myoption == "FIX":
        fn = make_massfn(allx, ally, allf, vols, omega_prior, fix_alpha=ALPHA_FIX)
        par, val, nfe, conv = r_optim_nm(np.array([mstarmrp, phimrp, betamrp]), fn,
                                         maxit=maxit, reltol=1e-8,
                                         parscale=np.array([1.0, 1.0, 0.1]))
        return np.array([par[0], par[1], ALPHA_FIX, par[2]]), val, nfe, conv
    fn = make_massfn(allx, ally, allf, vols, omega_prior)
    return r_optim_nm(np.array([mstarmrp, phimrp, alphamrp, betamrp]), fn,
                      maxit=maxit, reltol=1e-8,
                      parscale=np.array([1.0, 1.0, 1.0, 0.1]))


def monte_carlo(sets, myoption, omega_prior, phimrp, mstarmrp, alphamrp, betamrp,
                cosvars, iters, seed, verbose=True):
    """allhmf.r lines 352-416.

    Two quirks of the original are reproduced.  The cosmic-variance term is
    applied in *log* space (``V4`` is log10 phi, and ``cv = V4 * rnorm(...)``)
    while the measurement error is applied in linear space.  And line 354 draws
    ``rnorm(length(gama$V4))`` for REFLEX, not ``length(reflexy)``, so with 14
    GAMA bins and 43 REFLEX points R recycles the 14 deviates across the 43
    values; ``np.resize`` reproduces that.
    """
    rng = RRandom(seed)
    out = {k: np.full(iters, np.nan) for k in
           ("mstar", "phistar", "alphastar", "betastar", "omegam", "omegam2")}
    curves = []
    ngama = len(sets["G"][0])
    for i in range(iters):
        cv = {}
        # draw order as in R: reflex, sdss, gama -- all drawn regardless of use
        cv["R"] = sets["R"][1] * np.resize(rng.norm(ngama, COSVAR_REFLEX),
                                           len(sets["R"][1])) if "R" in sets else None
        if "S" in sets:
            cv["S"] = sets["S"][1] * rng.norm(len(sets["S"][1]), cosvars["S"])
        cv["G"] = sets["G"][1] * rng.norm(ngama, cosvars["G"])

        pert = {k: (v[0], v[1] + (cv[k] if cv.get(k) is not None else 0.0), v[2], v[3])
                for k, v in sets.items()}
        allx, ally_log, allf, vols = assemble(pert, myoption)
        allf = np.where(np.isnan(allf), 0.9999, allf)
        allf = np.where(np.isinf(allf), 0.0, allf)
        allf = np.where(allf >= 1.0, 0.9999, allf)

        yy = 10 ** ally_log
        mock = yy + yy * rng.norm(len(allf), allf)
        m = (mock > 0) & ~np.isnan(mock)
        if m.sum() < 5:
            continue
        with np.errstate(divide="ignore"):
            par, _, _, _ = fit_combined(allx[m], np.log10(mock[m]), allf[m], vols,
                                        myoption, omega_prior, phimrp, mstarmrp,
                                        alphamrp, betamrp)
        out["mstar"][i], out["phistar"][i] = par[0], par[1]
        out["alphastar"][i], out["betastar"][i] = par[2], par[3]
        out["omegam"][i] = omega_matter(par)
        out["omegam2"][i] = omega_matter(par, above=MLIMIT_GAMA)
        if i < 1000:
            curves.append(par)
        if verbose and (i + 1) % 200 == 0:
            print(f"    {i + 1}/{iters}", end="\r", flush=True)
    if verbose:
        print(" " * 30, end="\r")
    out["curves"] = curves
    return out


# ---------------------------------------------------------------------------
# Plot -- allhmf.r lines 293-512
# ---------------------------------------------------------------------------

def plot_combined(fit_par, mc, sets, driver_gama, extras, mrpx, mrpy, factor,
                  outfile, myoption, omega_prior, fit_par_driver=None,
                  sdss_label="SDSS DR10 (Tempel+14)", fit_method=None,
                  figsize=(7.87, 4.72), pub_fit=None, pub_label=None,
                  show_extras=True, show_driver_gama=True,
                  sdss_label_short=None, show_driver_fit=True):
    """The combined figure.

    ``show_extras=False`` drops the REFLEX, 2PIGG and Tempel-curve *points*
    (and their legend rows), leaving only the survey data actually being
    fitted; the LCDM line and Driver's published fit stay, as references.
    ``pub_fit`` overrides which published fit is drawn behind -- the GAMA-only
    figure needs his GAMA5 row, not the GSR one every combined figure uses.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import plotting                       # the paper's house style
    from matplotlib.collections import LineCollection

    # A 3.54 in panel cannot carry the full-page legend: the same 8 pt text is
    # proportionally twice as wide and runs across the data.  Shorten the rows
    # and shrink the type when the figure is narrow.
    compact = figsize[0] < 5.0
    lfs = 5.6 if compact else 8.0

    def _lab(long, short):
        return short if compact else long

    # Each fitted curve takes the colour of the data it was fitted to: red for
    # the Nessie GAMA points, cornflower for Driver's.  These used to be the
    # other way round -- our fit was drawn in exactly the cornflower of HIS
    # data points, which invites the reader to pair the wrong line with the
    # wrong sample.
    cf = OUR_FIT_C
    tp_x, tp_y, tp_up, tp_do = extras["tpigg"]
    el_x, el_y, el_lo, el_hi = extras["elmo"]

    fig = plt.figure(figsize=figsize, dpi=600)
    ax = fig.add_axes([0.105, 0.115, 0.885, 0.875])
    xfit = np.arange(0.0, 18 + 1e-9, 0.01)

    if mc is not None and mc["curves"]:
        segs = []
        for par in mc["curves"][:1000]:
            with np.errstate(divide="ignore", invalid="ignore"):
                yy = np.log10(dr.mrp(xfit, *par))
            m = np.isfinite(yy) & (yy > -9) & (xfit > 12.5) & (xfit < 16.2)
            if m.sum() > 2:
                segs.append(np.column_stack([xfit[m], yy[m]]))
        # rasterised: ~1000 curves make a vector PDF enormous and slow to
        # open, while the axes, text and headline curve stay vector.
        ax.add_collection(LineCollection(segs, colors=[cf], linewidths=0.6,
                                         alpha=0.01, rasterized=True))

    # Driver+22's *published* table 2 fit for this sample combination, drawn
    # first so it sits behind the data.  This is his printed answer, not our
    # refit of it -- the two are different lines and are labelled as such.
    pub = pub_fit if pub_fit is not None else DRIVER_ABSTRACT_FIT
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xfit, np.log10(dr.mrp(xfit, pub[0], 10 ** pub[1], pub[2],
                                      pub[3])),
                color=PUBLISHED_C, lw=3.2, alpha=0.9, zorder=2)

    if show_extras:
        ax.fill_between(el_x, el_lo, el_hi, color=(0.5, 0.5, 0.5, 0.10), lw=0)
        ax.plot(el_x, el_y, color="cyan", lw=1.2, zorder=2)
    with np.errstate(divide="ignore"):
        ax.plot(mrpx - 0.08, np.log10(mrpy) - np.log10(factor) + 0.08, ls="--",
                color="black", lw=2, zorder=3)

    if show_extras:
        rx, ry, rf, _ = sets["R"]
        ax.errorbar(rx, ry,
                    yerr=[np.abs(np.log10(1 - rf)), np.abs(np.log10(1 + rf))],
                    fmt="D", ms=4, color="forestgreen", elinewidth=0.9,
                    capsize=0, ls="none", zorder=5)
        ax.errorbar(tp_x, tp_y, yerr=[np.abs(tp_do), np.abs(tp_up)], fmt="o",
                    ms=4, color="grey", elinewidth=0.9, capsize=0, ls="none",
                    zorder=5)
    if "S" in sets:
        sx, sy, sf, _ = sets["S"]
        ax.errorbar(sx, sy, yerr=[np.abs(np.log10(1 - sf)), np.abs(np.log10(1 + sf))],
                    fmt="P", ms=5, color="purple", elinewidth=0.9, capsize=0,
                    ls="none", zorder=5)

    # Driver's own GAMA, for comparison
    dgx, dgy, dgf = driver_gama
    if not show_driver_gama:
        dgx, dgy, dgf = dgx[:0], dgy[:0], dgf[:0]
    ax.errorbar(dgx, dgy, yerr=[np.abs(np.log10(1 - dgf)), np.abs(np.log10(1 + dgf))],
                fmt="o", ms=5, mfc="none", mec=ng.OLD_C, ecolor=ng.OLD_C,
                elinewidth=0.9, capsize=0, ls="none", zorder=6)
    gx, gy, gf, _ = sets["G"]
    ax.errorbar(gx, gy, yerr=[np.abs(np.log10(1 - gf)), np.abs(np.log10(1 + gf))],
                fmt="o", ms=5, color="red", elinewidth=0.9, capsize=0, ls="none",
                zorder=7)

    # the identical fit run on Driver's GAMA, faint, so the effect of swapping
    # only the group catalogue can be read straight off the figure
    if fit_par_driver is not None and show_driver_gama and show_driver_fit:
        with np.errstate(divide="ignore", invalid="ignore"):
            # Solid and heavier, underneath our fit, which is dotted with wide
            # gaps.  The two curves nearly coincide by construction, so the only
            # way both stay visible is to let one show THROUGH the other; a
            # dashed line under a dotted line just fragments into invisibility.
            ax.plot(xfit, np.log10(dr.mrp(xfit, *fit_par_driver)),
                    color=DRIVER_FIT_C, lw=2.6, ls="-", zorder=4)
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xfit, np.log10(dr.mrp(xfit, *fit_par)), color=cf, lw=1.8,
                ls=(0, (1, 2.5)), zorder=8)

    ax.set_xlim(12.75, 16)
    ax.set_ylim(-8, -2)
    ax.set_xlabel(r"log$_{10}$(Halo Mass / M$_\odot$)",
                  fontsize=9 if compact else 10)
    ax.set_ylabel(r"log$_{10}$(number density) [Mpc$^{-3}$ dex$^{-1}$]",
                  fontsize=9 if compact else 10)
    ax.tick_params(direction="in", top=True, right=True, which="both",
                   labelsize=7.5 if compact else 9)
    ax.minorticks_on()

    keys = [("o", "red", _lab(f"GAMA (Nessie DMU) z<{ZLIMIT} and N>4",
                              "GAMA (Nessie)"), False),
            ("o", ng.OLD_C, _lab("GAMA (Driver+22) z<0.25 and N>4",
                                 "GAMA (Driver+22)"), True)]
    if not show_driver_gama:
        keys = keys[:1]
    if "S" in sets:
        keys.append(("P", "purple",
                     _lab(f"{sdss_label} z<{ZLIMITSDSS} and N>4",
                          sdss_label_short or sdss_label), False))
    if show_extras:
        keys += [("D", "forestgreen",
                  "REFLEX II, x-ray, z~0.1 (Bohringer et al. 2017)", False),
                 ("o", "grey", "2PIGG z < 0.12 (Eke et al. 2008)", False)]
    # The legend is hand-laid-out, so its spacing has to adapt to how many
    # entries this particular combination produces -- otherwise adding one
    # (as the published fit did) pushes the last row off the bottom of the
    # axes and into the tick labels.
    n_rows = len(keys) + 2 + show_extras + (fit_par_driver is not None
                                            and show_driver_gama
                                            and show_driver_fit)
    y0, y_floor = -4.95, -7.62
    dy = min(0.325, (y0 - y_floor) / max(n_rows - 1, 1))
    for k, (mk, c, lab, hollow) in enumerate(keys):
        yy = y0 - dy * k
        ax.plot([12.85], [yy], marker=mk, ms=4 if compact else 5, color=c,
                mfc=("none" if hollow else c), mec=c)
        ax.text(12.95, yy, " " + lab, color=c, va="center", fontsize=lfs)
    yy = y0 - dy * (len(keys) - 1)
    if show_extras:
        yy -= dy
        ax.plot([12.8, 12.9], [yy, yy], color="cyan")
        ax.text(12.95, yy, " SDSS DR10 mass function (Tempel et al. 2014)",
                color="cyan", va="center", fontsize=lfs)
    yy -= dy
    ax.plot([12.8, 12.9], [yy, yy], color=cf, lw=1.8, ls=(0, (1, 2.5)))
    # With --mcmc this line is the best-lnP posterior sample and the band is
    # posterior draws, not Driver's optimiser point and refits.  Say so.
    ax.text(12.95, yy, _lab(f" Best fit MRP function to {myoption}",
                            f" Best fit to {myoption}")
            + (" + $\\Omega_M$ prior" if omega_prior else "")
            + (f" ({fit_method})" if fit_method else ""), va="center",
            fontsize=lfs)
    if fit_par_driver is not None and show_driver_gama and show_driver_fit:
        yy -= dy
        ax.plot([12.8, 12.9], [yy, yy], color=DRIVER_FIT_C, lw=2.6, ls="-")
        ax.text(12.95, yy, _lab(" same fit using Driver+22 GAMA",
                                " same fit, Driver+22 GAMA"), va="center",
                fontsize=lfs, color=DRIVER_FIT_C)
    yy -= dy
    ax.plot([12.8, 12.9], [yy, yy], color=PUBLISHED_C, lw=3.2, alpha=0.9)
    ax.text(12.95, yy, " " + (pub_label or "Driver+22 published GSR fit"),
            va="center", fontsize=lfs, color=PUBLISHED_C)
    yy -= dy
    ax.plot([12.8, 12.9], [yy, yy], color="black", ls="--", lw=2)
    ax.text(12.95, yy, " $\\Lambda$CDM expectation from MRP", va="center",
            fontsize=lfs)

    # Omega_M inset, allhmf.r lines 494-510
    if mc is not None:
        iax = fig.add_axes([0.63, 0.60, 0.345, 0.375])
        om = mc["omegam"][np.isfinite(mc["omegam"])]
        iax.hist(om, bins=np.arange(0, 1.001, 0.01), color="0.55")
        iax.axvspan(np.quantile(om, 0.16), np.quantile(om, 0.84),
                    color=(1, 0, 0, 0.25), lw=0)
        iax.axvline(OMEGAM, color="black", lw=2, ls=":")
        iax.axvline(np.quantile(om, 0.50), color="red", lw=1)
        iax.set_xlim(0, 1)
        iax.set_xlabel(r"$\Omega_M$ (total)", fontsize=7 if compact else 9)
        iax.set_ylabel("Frequency", fontsize=8)
        iax.tick_params(labelsize=7, direction="in")
        lo, hi = np.quantile(om, 0.16), np.quantile(om, 0.84)
        med = np.quantile(om, 0.50)
        # .2g rendered 0.199 as "0.2", which is not a publication number.
        iax.text(0.97, 0.92, r"$\Omega_M$ = " + f"{med:.3f}"
                 + f"$^{{+{hi - med:.3f}}}_{{-{med - lo:.3f}}}$",
                 transform=iax.transAxes, ha="right", va="top", color="red",
                 fontsize=9)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  wrote {outfile}")


# ---------------------------------------------------------------------------

def gama_set(g3c, vlimit, area, seed, nmc_edb, nboot, mlimit, fit_max=None,
             fit_min=None):
    """Bin a GAMA catalogue and reduce it to allhmf.r's (V1, V4, V8) + volume."""
    b, _ = bin_hmf(g3c, vlimit, seed=seed, nmc=nmc_edb, verbose=False, nboot=nboot)
    t = driver_table(b)
    lo = mlimit if fit_min is None else max(mlimit, fit_min)
    keep = (t.V1 > lo) & np.isfinite(t.V4)              # allhmf.r line 276
    if fit_max is not None:
        keep &= t.V1 <= fit_max
    t = t[keep]
    return (t.V1.values, t.V4.values, t.V8.values,
            survey_volume(ZLIMIT, area)), b


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--myoption", default="GSR",
                   choices=["G", "R", "GR", "GS", "SR", "GSR", "Omega", "FIX"])
    p.add_argument("--omega-prior", action="store_true",
                   help="add the OmegaM term (implied by --myoption Omega)")
    p.add_argument("--sdss", default="auto",
                   help="'auto' rebuilds Tempel+14 with sdss_hmf.py, 'nessie' "
                        "uses the Nessie SDSS catalogue, 'none' disables, or a "
                        "path to a V1..V8 table")
    p.add_argument("--sdss-max", type=float, default=None,
                   help="drop SDSS bins above this logM before fitting.  The "
                        "penalty integrates over max(allx), and Nessie SDSS "
                        "reaches 15.65 against Tempel's 15.05, so this is the "
                        "range check group 300223 taught us to run.")
    p.add_argument("--mass-mode", default="tempel_eq8",
                   choices=["tempel_eq8", "tempel_nfw", "tempel_rms", "robotham", "shift"],
                   help="mass estimator for the NESSIE SDSS leg only (--sdss "
                        "nessie).  'tempel_eq8' is the catalogue default, a "
                        "Hernquist mass with Nessie's cbrt(3) bug. "
                        "'tempel_nfw' corrects both and puts it on the same "
                        "NFW scale as Tempel's published col15, which is what "
                        "Driver's SDSS leg uses.")
    p.add_argument("--mass-shift", type=float, default=0.0,
                   help="dex shift, with --mass-mode shift")
    p.add_argument("--ml-cut", type=float, default=None,
                   help="drop groups whose mass-to-light ratio exceeds the "
                        "running median at their mass by more than this many "
                        "dex (1.0 recommended).  Targets the mass being wrong "
                        "rather than capping the weight of a mass that is "
                        "believed, so it is preferred to --vmax-floor, which "
                        "is a biased estimator.  Applied to both GAMA legs.")
    p.add_argument("--vmax-floor", type=float, default=1e-3,
                   help="minimum Vmax as a fraction of the survey volume, i.e. "
                        "a cap on 1/Vmax.  Driver's own vlimitmin is 1e-3, "
                        "which never binds.  Raising it suppresses nearby "
                        "low-Vmax groups (e.g. GAMA 205509, 38%% of the 14.2 "
                        "bin) but is a BIASED estimator whose bias is "
                        "mass-dependent -- see CLAUDE.md before using it for "
                        "anything but a robustness check.")
    p.add_argument("--fit-max", type=float, default=None,
                   help="clamp BOTH survey legs (GAMA and SDSS) at this logM.  "
                        "The high-mass tails are single-group bins with "
                        "fractional errors of ~1, so they add no constraint but "
                        "do set max(allx) and hence the penalty range.  REFLEX "
                        "is deliberately left uncapped: it is the anchor, and "
                        "its high-mass points are the constraint, not noise.")
    p.add_argument("--name", default=None,
                   help="basename for this deliverable: writes hmf_<name>.pdf, "
                        "corner_<name>.pdf and chain_<name>.npz.  Overrides the "
                        "tag-derived filenames, which are unreadable once more "
                        "than two options are set.")
    p.add_argument("--profile-mstar", action="store_true",
                   help="print the chi2 profile in log10 M* (the other three "
                        "parameters minimised at each point) and exit.  Says "
                        "whether a railed posterior means the data prefer a "
                        "low M* or simply do not constrain it.")
    p.add_argument("--no-driver-gama", action="store_true",
                   help="drop Driver+22's GAMA comparison points and the "
                        "'same fit using Driver+22 GAMA' curve, leaving only "
                        "the Nessie data sets.")
    p.add_argument("--no-extras", action="store_true",
                   help="drop the REFLEX, 2PIGG and Tempel-curve points and "
                        "their legend rows, leaving only the survey data being "
                        "fitted.  Model curves are kept.")
    p.add_argument("--full-page", action="store_true",
                   help="size the figure to span both columns (7.1 in) rather "
                        "than the default 20x12 cm.")
    p.add_argument("--pub-fit", default="gsr", choices=["gsr", "gama", "gs"],
                   help="which Driver+22 table 2 fit to draw behind: his GSR "
                        "(default), his GAMA-only GAMA5 row, or his GAMA+SDSS "
                        "GS row -- each figure should show the row fitted to "
                        "the same sample combination it plots.")
    p.add_argument("--no-driver-fit", action="store_true",
                   help="drop the 'same fit using Driver+22 GAMA' curve while "
                        "keeping his data points.  On the combined figures his "
                        "published fit lies on top of it, so the extra line "
                        "adds nothing.")
    p.add_argument("--mass-err", default="driver",
                   help="mass-error curve for the Nessie GAMA leg: 'driver' "
                        "(his hardcoded NFOF_YY/MASSCORR, the default, which "
                        "the reproductions depend on) or a vuvuzela CSV such "
                        "as gama_masserr.csv, which also replaces MASSCORR "
                        "with the measured multiplicity debiasing.")
    p.add_argument("--sdss-mass-err", default="driver",
                   help="mass-error curve for the Nessie SDSS leg, e.g. "
                        "sdss_masserr.csv.  Scatter only: the debiasing is "
                        "used by the 'robotham' mass mode alone.  Driver's "
                        "Tempel leg always keeps his own curve.")
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--nmc-edb", type=int, default=1001)
    p.add_argument("--nboot", type=int, default=0,
                   help="group-bootstrap errors instead of Poisson, applied to "
                        "both the GAMA and SDSS legs (0 = as Driver)")
    p.add_argument("--iters", type=int, default=1001)
    p.add_argument("--maxit", type=int, default=500,
                   help="optim budget (allhmf.r uses 500).  GSR needs ~2000 to "
                        "actually converge; check whether a convergence=1 "
                        "result is a minimum or just a stopped simplex.")
    p.add_argument("--fit-min", type=float, default=None,
                   help="raise the LOWER mass cut on both survey legs.  Driver "
                        "uses 12.7 (GAMA) / 12.9 (SDSS), allhmf.r lines 276/278; "
                        "this only ever tightens them, never loosens.  REFLEX is "
                        "left alone, as it is with --fit-max.")
    p.add_argument("--mcmc", action="store_true",
                   help="sample the posterior with emcee on the SAME objective "
                        "the Nelder-Mead fit uses, and write a corner plot "
                        "overlaying it on the Monte-Carlo refit draws.  The "
                        "penalty term is the Poisson zero-count likelihood, so "
                        "chi^2 is already -2lnL and needs no reinterpretation.")
    p.add_argument("--mcmc-steps", type=int, default=8000)
    p.add_argument("--mcmc-walkers", type=int, default=64)
    p.add_argument("--out", default=None,
                   help="default depends on --sdss so the two SDSS variants "
                        "cannot overwrite each other")
    p.add_argument("--write-tables", action="store_true")
    args = p.parse_args()

    # Built unconditionally: --mcmc names its corner plot from it too, so it must
    # exist and be complete even when --out was given explicitly.  Only the
    # canonical GSR runs get the canonical names; anything else -- a different
    # myoption, a dropped SDSS leg, a capped range -- earns a suffix, so a quick
    # diagnostic cannot silently overwrite a deliverable.
    tag = "" if args.myoption == "GSR" else f"_{args.myoption}"
    if args.sdss == "nessie":
        tag += "_sdss"
    elif args.sdss == "none":
        tag += "_nosdss"
    elif args.sdss != "auto":
        tag += "_sdssfile"
    if args.ml_cut is not None:
        tag += f"_mlcut{args.ml_cut:g}"
    if args.mass_mode != "tempel_eq8":
        tag += f"_{args.mass_mode.replace('tempel_', '')}"
    if args.vmax_floor != 1e-3:
        tag += f"_vfloor{args.vmax_floor:g}"
    if args.fit_min is not None:
        tag += f"_min{args.fit_min:g}"
    if args.fit_max is not None:
        tag += f"_clamp{args.fit_max:g}"
    elif args.sdss_max is not None:
        tag += f"_max{args.sdss_max:g}"
    if args.omega_prior:
        tag += "_omega"
    if args.iters != 1001:
        tag += f"_it{args.iters}"
    if args.mcmc:
        tag += "_mcmc"   # never overwrite the Nelder-Mead deliverable
    if args.out is None:
        args.out = (f"hmf_{args.name}.pdf" if args.name
                    else f"hmf_combined_nessie{tag}.png")

    omega_prior = args.omega_prior or args.myoption == "Omega"
    mrpx, mrpy, factor, phimrp = lcdm_curve()
    mstarmrp, alphamrp, betamrp = dr.MSTARMRP, dr.ALPHAMRP, dr.BETAMRP

    print("Combined HMF (allhmf.r) with the new Nessie GAMA DMU")
    print(f"  myoption = {args.myoption}   OmegaM prior = {omega_prior}")

    # --- GAMA: new Nessie DMU (fitted) and Driver's (comparison only) --------
    g_new, v_new = ng.build_groups_new(verbose=False,
                                       vmax_floor_frac=args.vmax_floor,
                                       masserr=args.mass_err)
    g_new, _ = ng.apply_ml_cut(g_new, args.ml_cut, verbose=args.ml_cut is not None)
    gset, b_new = gama_set(g_new, v_new, ng.NEW_AREA, args.seed, args.nmc_edb,
                           args.nboot, MLIMIT_GAMA, fit_max=args.fit_max,
                           fit_min=args.fit_min)
    g_old, v_old = dr.build_groups("../data/G3CFoFGroupv10.fits",
                                   "../data/GAMAGalsInGroups.csv", verbose=False,
                                   vmax_floor_frac=args.vmax_floor)
    g_old, _ = ng.apply_ml_cut(g_old, args.ml_cut, verbose=False)
    oset, b_old = gama_set(g_old, v_old, AREA_GAMA_DRIVER, args.seed, args.nmc_edb,
                           args.nboot, MLIMIT_GAMA, fit_max=args.fit_max,
                           fit_min=args.fit_min)
    print(f"  GAMA (Nessie)  : {len(gset[0])} bins, area {ng.NEW_AREA} deg^2, "
          f"volume {gset[3]:.4e} Mpc^3")
    print(f"  GAMA (Driver)  : {len(oset[0])} bins, area {AREA_GAMA_DRIVER} deg^2 "
          f"(allhmf.r line 263)")

    rx, ry, rf = load_reflex()
    sets = {"G": gset, "R": (rx, ry, rf, VOLUME_REFLEXII)}
    print(f"  REFLEX II      : {len(rx)} points, volume {VOLUME_REFLEXII:.3e} Mpc^3")

    sdss_label = None
    if args.sdss != "none":
        if args.sdss == "auto":
            s, _, _ = sdss_hmf.build(seed=args.seed, nmc=args.nmc_edb,
                                     verbose=False, nboot=args.nboot,
                                     vmax_floor_frac=args.vmax_floor)
            sdss_label = "SDSS DR10 (Tempel+14)"
        elif args.sdss == "nessie":
            # Same footprint (7221 deg^2) and the same galaxies -- the Nessie
            # file is Tempel's table 1 -- so only the grouping differs.
            s, _, _ = nsd.build(seed=args.seed, nmc=args.nmc_edb,
                                verbose=False, nboot=args.nboot,
                                vmax_floor_frac=args.vmax_floor,
                                mass_mode=args.mass_mode,
                                mass_shift=args.mass_shift,
                                masserr=args.sdss_mass_err)
            sdss_label = {
                "tempel_eq8": "SDSS (Nessie)",
                "tempel_nfw": "SDSS (Nessie, NFW mass)",
                "tempel_rms": "SDSS (Nessie, NFW mass + Tempel rms $\\sigma$)",
                "robotham": "SDSS (Nessie, Robotham mass as GAMA)",
                "shift": f"SDSS (Nessie, {args.mass_shift:+.3f} dex)",
            }[args.mass_mode]
        else:
            s = pd.read_csv(args.sdss)
            sdss_label = "SDSS (file)"
        lo_s = MLIMIT_SDSS if args.fit_min is None else max(MLIMIT_SDSS, args.fit_min)
        s = s[(s.V1 > lo_s) & np.isfinite(s.V4)]
        cap = min([c for c in (args.sdss_max, args.fit_max) if c is not None],
                  default=None)
        if cap is not None:
            s = s[s.V1 <= cap]
        # allhmf.r line 264 does NOT subtract the zmin volume here, though
        # sdsshmf.r line 248 does.  His value is kept in his place.
        sets["S"] = (s.V1.values, s.V4.values, s.V8.values,
                     survey_volume(ZLIMITSDSS, AREA_SDSS))
        print(f"  {sdss_label:<16}: {len(s)} bins, area {AREA_SDSS} deg^2, "
              f"z<{ZLIMITSDSS}, logbin 0.1, max bin {s.V1.max():.2f}")

    extras = dict(tpigg=load_tpigg(), elmo=load_elmo())

    if args.write_tables:
        driver_table(b_new).to_csv("gamahmfNessie5.csv", index=False)
        driver_table(b_old).to_csv("gamahmfGAMA5.csv", index=False)
        print("  wrote gamahmfNessie5.csv and gamahmfGAMA5.csv "
              "(Driver's V1..V8 column order)")

    cosvars = {"G": cosvar(gset[3] / 3, 3), "S": cosvar(
        sets["S"][3], 1) if "S" in sets else 0.0}
    print(f"  cosvar GAMA    : {cosvars['G']:.5f}")

    # --- fit -----------------------------------------------------------------
    allx, ally, allf, vols = assemble(sets, args.myoption)
    allf = np.where(allf >= 1.0, 0.9999, allf)

    if args.profile_mstar:
        # chi2 profile in log10 M*: fix it, minimise over the other three.
        # This is what says whether a railed posterior means "the data prefer
        # M* = 11" or "the data do not constrain M* at all" -- very different
        # statements, and only the profile distinguishes them.
        from scipy.optimize import minimize
        fn = make_massfn(allx, ally, allf, vols, omega_prior)
        grid = np.arange(11.0, 15.01, 0.25)
        best_c, prof = np.inf, []
        for ms in grid:
            r = min((minimize(lambda p: fn(np.array([ms, p[0], p[1], p[2]])),
                              x0=np.array([p0, a0, b0]), method="Nelder-Mead",
                              options=dict(maxiter=4000, xatol=1e-6,
                                           fatol=1e-8))
                     for p0 in (10 ** -3.2, 10 ** -2.6)
                     for a0 in (-1.5, -0.5)
                     for b0 in (0.3, 0.6)), key=lambda r: r.fun)
            prof.append(r.fun)
            best_c = min(best_c, r.fun)
        print(f"\n  chi2 profile in log10 M*  ({args.myoption}, "
              f"{len(allx)} points)")
        print(f"  {'logM*':>7s} {'chi2':>10s} {'dchi2':>8s}")
        for ms, c in zip(grid, prof):
            print(f"  {ms:7.2f} {c:10.2f} {c - best_c:8.2f}")
        return

    par, val, nfe, conv = fit_combined(allx, ally, allf, vols, args.myoption,
                                       omega_prior, phimrp, mstarmrp, alphamrp,
                                       betamrp, maxit=args.maxit)
    print(f"\n  fitting {len(allx)} points from {args.myoption}")
    print(f"  Monte-Carlo: {args.iters} refits ...")
    mc = monte_carlo(sets, args.myoption, omega_prior, phimrp, mstarmrp, alphamrp,
                     betamrp, cosvars, args.iters, args.seed)

    # the same fit with Driver's GAMA, so the catalogue swap is quantified
    sets_drv = dict(sets)
    sets_drv["G"] = oset
    ax_d, ay_d, af_d, vol_d = assemble(sets_drv, args.myoption)
    af_d = np.where(af_d >= 1.0, 0.9999, af_d)
    par_d, val_d, _, conv_d = fit_combined(ax_d, ay_d, af_d, vol_d, args.myoption,
                                           omega_prior, phimrp, mstarmrp, alphamrp,
                                           betamrp, maxit=args.maxit)

    print("\n  --- combined MRP fit ---")
    print(f"  {'':<13}{'Nessie DMU':>10}{'Driver+22':>11}{'diff':>8}")
    for k, nm in enumerate(["log10(M*)", "log10(phi*)", "alpha", "beta"]):
        a = par[k] if k != 1 else np.log10(abs(par[1]))
        bq = par_d[k] if k != 1 else np.log10(abs(par_d[1]))
        print(f"  {nm:<13}{a:10.3f}{bq:11.3f}{a - bq:+8.3f}")
    print(f"  {'chi2':<13}{val:10.2f}{val_d:11.2f}"
          f"   (convergence {conv} / {conv_d})")
    par_plot, mc_plot, par_d_plot = par, mc, par_d   # --mcmc overrides these
    if args.mcmc:
        import mcmc_hmf as mh

        print(f"\n  --- MCMC on the same objective ({args.myoption}) ---")
        fn = make_massfn(allx, ally, allf, vols, omega_prior)
        chain, best, info = mh.run_emcee(
            fn, par, nwalkers=args.mcmc_walkers, nsteps=args.mcmc_steps,
            burn=args.mcmc_steps // 4, seed=args.seed)
        nm = np.array([par[0], np.log10(abs(par[1])), par[2], par[3]])
        print(f"  {info['nsamples']} samples, acceptance {info['acceptance']:.2f}, "
              f"tau {np.array2string(info['tau'], precision=0)}, "
              f"n_eff {np.array2string(info['neff'], precision=0)}")
        if not info["converged"]:
            print("  WARNING: n_eff < 50 for at least one parameter -- the "
                  "credible intervals below are not yet reliable, run longer.")
        mh.summarise(chain, f"posterior, {args.myoption}", nm=nm)
        print(f"\n  best-lnP sample vs Nelder-Mead (Jacobian-free comparison):")
        for i, n in enumerate(mh.PARAM_NAMES):
            print(f"  {n:10s}{best[i]:9.3f}{nm[i]:11.3f}{best[i] - nm[i]:+9.3f}")
        print(f"  {'chi2':10s}{info['chi2_min']:9.2f}{val:11.2f}"
              f"{info['chi2_min'] - val:+9.3f}")

        # Swap the posterior in EVERYWHERE the Monte-Carlo refits were used --
        # the spaghetti band, the Omega_M inset and the printed summary -- not
        # just the corner plot.  Leaving the figure on the Nelder-Mead curve
        # while the corner plot showed the posterior would put two different
        # answers in the same deliverable.  The refits are kept for the corner
        # overlay, which is where the comparison belongs.
        rng_pick = np.random.default_rng(args.seed)
        idx = rng_pick.choice(len(chain), size=min(1000, len(chain)),
                              replace=False)
        lin = np.column_stack([chain[:, 0], 10.0 ** chain[:, 1],
                               chain[:, 2], chain[:, 3]])
        mc_plot = dict(
            mstar=chain[:, 0], phistar=lin[:, 1],
            alphastar=chain[:, 2], betastar=chain[:, 3],
            curves=[lin[i] for i in idx],
            omegam=np.array([omega_matter(p) for p in lin[idx]]),
            omegam2=np.array([omega_matter(p, above=MLIMIT_GAMA)
                              for p in lin[idx]]),
        )
        par_plot = np.array([best[0], 10.0 ** best[1], best[2], best[3]])

        # The "same fit using Driver+22 GAMA" comparison line must come from the
        # SAME method, or the figure silently compares an MCMC curve against an
        # optimiser point and the catalogue difference is confounded with the
        # fitting method.  Second chain, on Driver's GAMA data.
        fn_d = make_massfn(ax_d, ay_d, af_d, vol_d, omega_prior)
        chain_d, best_d, info_d = mh.run_emcee(
            fn_d, par_d, nwalkers=args.mcmc_walkers, nsteps=args.mcmc_steps,
            burn=args.mcmc_steps // 4, seed=args.seed + 1)
        if not info_d["converged"]:
            print("  WARNING: Driver-GAMA comparison chain has n_eff < 50")
        par_d_plot = np.array([best_d[0], 10.0 ** best_d[1], best_d[2], best_d[3]])
        print(f"  Driver-GAMA comparison chain: chi2 {info_d['chi2_min']:.2f} "
              f"vs NM {val_d:.2f}, logM* {best_d[0]:.3f} vs {par_d[0]:.3f}")
        print("  figure, Omega_M inset and summary now use the POSTERIOR; "
              "the headline curve is the best-lnP sample")

        # Only our posterior is shown; the point estimates are stars.  The MC
        # refits and the Nelder-Mead crosshair are deliberately not drawn -- the
        # crosshair was our own refit of our own catalogue and read as if it were
        # Driver's published answer, which it is not.
        med = np.median(chain, axis=0)
        pub_marker = PUB_FIT[args.pub_fit][0]
        markers = [
            (np.array(pub_marker),
             "Driver+22 published "
             + {"gama": "GAMA-only", "gs": "GAMA+SDSS"}.get(args.pub_fit,
                                                            "GSR"), "white"),
            (med, "This work (posterior median)", "#d4af37"),
        ]
        sets_t0 = [(chain, f"This work, MCMC posterior ({args.myoption})",
                    "#4878a8")]
        np.savez(f"chain_{args.name}.npz" if args.name else f"chain{tag}.npz",
                 chain=chain, chain_driver=chain_d,
                 median=med, best=best, nm=nm, driver_published=pub_marker,
                 chi2_best=info["chi2_min"], chi2_driver=info_d["chi2_min"],
                 chi2_nm=val, chi2_nm_driver=val_d,
                 # the inset's own draws, so mrp_table quotes the identical
                 # numbers rather than re-subsampling and differing in the
                 # third decimal
                 omega_draws=mc_plot["omegam"],
                 omega_draws_above=mc_plot["omegam2"])
        cout = (f"corner_{args.name}.pdf" if args.name
                else "corner" + (tag if tag else "_GSR") + ".png")
        mh.corner_plot(sets_t0, cout, markers=markers,
                       title=f"{args.myoption}: MRP posterior, this work")

    pub = PUBLISHED_TABLE2.get(args.myoption)
    clamped = args.fit_max is not None or args.sdss_max is not None
    if pub is not None and not omega_prior and clamped:
        print("\n  (published table 2 comparison suppressed: Driver does not "
              "clamp the mass\n   range, so a clamped fit is not comparable "
              "with it.  Re-run without\n   --fit-max/--sdss-max to validate "
              "the port.)")
    elif pub is not None and not omega_prior:
        print(f"\n  Validation -- Driver's GAMA against his published table 2 "
              f"({args.myoption}):")
        print(f"  {'':<13}{'this port':>10}{'published':>11}{'diff':>8}")
        for k, nm in enumerate(["log10(M*)", "log10(phi*)", "alpha", "beta"]):
            v = par_d[k] if k != 1 else np.log10(abs(par_d[1]))
            print(f"  {nm:<13}{v:10.3f}{pub[k]:11.2f}{v - pub[k]:+8.3f}")

    print()
    src = "MCMC posterior" if args.mcmc else "Monte-Carlo refits"
    names = ["log10(M*)", "log10(phi*)", "alpha", "beta"]
    vals = [par_plot[0], np.log10(abs(par_plot[1])), par_plot[2], par_plot[3]]
    chains = [mc_plot["mstar"], np.log10(np.abs(mc_plot["phistar"])),
              mc_plot["alphastar"], mc_plot["betastar"]]
    # Intervals come from whatever the FIGURE uses, so the printed numbers and
    # the plotted band can never disagree.
    print(f"  {'param':<13}{'best':>9}{'16th':>9}{'84th':>9}   ({src})")
    for nm, v, c in zip(names, vals, chains):
        c = c[np.isfinite(c)]
        print(f"  {nm:<13}{v:9.3f}{np.quantile(c, .16):9.3f}"
              f"{np.quantile(c, .84):9.3f}")
    om2 = omega_matter(par_plot, above=MLIMIT_GAMA)
    o2 = mc_plot["omegam2"][np.isfinite(mc_plot["omegam2"])]
    print(f"  chi2 = {val:.3f}   fevals = {nfe}   convergence = {conv}"
          + ("  (budget exhausted)" if conv == 1 else ""))
    # The TOTAL is the headline: it is what the inset histogram plots and what
    # we quote.  Its systematic dwarfs its statistical error -- under a uniform
    # mass-scale shift the total scales as exactly 10**delta (no boundary term),
    # and this project's measured calibration spread is ~0.15 dex -- so print
    # the band, because otherwise nobody applies it.
    # Take the MEDIAN OF THE DRAWS, not omega_matter(par_plot).  The inset
    # annotates the median and mrp_table reads the same draws, so using the
    # point estimate here made the terminal disagree with its own figure in the
    # third decimal.
    o1 = mc_plot["omegam"][np.isfinite(mc_plot["omegam"])]
    omt = float(np.median(o1))
    SYS = 0.15
    print(f"\n  OmegaM (TOTAL, = the inset)  = {omt:.4f} "
          f"(+{np.quantile(o1, .84) - omt:.4f} / -{omt - np.quantile(o1, .16):.4f})"
          f"  stat")
    print(f"    x 10^(+/-{SYS} dex) mass scale -> {omt * 10 ** -SYS:.4f} .. "
          f"{omt * 10 ** SYS:.4f}   syst")
    print(f"    vs Planck OmegaM = {OMEGAM}: {omt / OMEGAM:.3f} of it")
    # Driver's constant, allhmf.r 415/462.  Our fit floor may now be higher than
    # 12.7 (--fit-min), in which case this carries some extrapolation below the
    # data -- another reason the total is the cleaner thing to quote.
    print(f"  OmegaM (logM > {MLIMIT_GAMA}, Driver's cut) = {om2:.4f} "
          f"(+{np.quantile(o2, .84) - om2:.4f} / -{om2 - np.quantile(o2, .16):.4f})")

    import plotting as _pl
    plot_combined(par_plot, mc_plot, sets, oset[:3], extras, mrpx, mrpy, factor,
                  args.out, args.myoption, omega_prior,
                  fit_par_driver=par_d_plot,
                  sdss_label=sdss_label or "SDSS DR10 (Tempel+14)",
                  sdss_label_short={"tempel_rms": "SDSS (Nessie, rms $\\sigma$)",
                                    "tempel_nfw": "SDSS (Nessie, NFW)",
                                    "tempel_eq8": "SDSS (Nessie)"}.get(
                                        args.mass_mode) if args.sdss == "nessie"
                                    else None,
                  fit_method="MCMC posterior" if args.mcmc else None,
                  figsize=(7.10, 4.26) if args.full_page else _pl.SINGLE_COLUMN,
                  pub_fit=PUB_FIT[args.pub_fit][0],
                  pub_label=PUB_FIT[args.pub_fit][1],
                  show_extras=not args.no_extras,
                  show_driver_gama=not args.no_driver_gama,
                  show_driver_fit=not args.no_driver_fit)


if __name__ == "__main__":
    main()
