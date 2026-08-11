#!/usr/bin/env python3
"""Python reproduction of figure 4 of Driver et al. (2022).

This is a line-by-line port of ``gamahmf.r``: the binned GAMA 1/Vmax halo mass
function, the Eddington-bias correction, the MRP fit and the Monte-Carlo error
band, on the old GAMA group catalogue.

Reproducing the published numbers requires reproducing two pieces of R exactly,
because the answer depends on both:

* **R's random number stream.**  The Eddington correction is a Monte-Carlo ratio
  of counts, so the binned points carry MC noise.  The noise is small (~0.005
  dex) but the fit amplifies it (see below), so ``rnorm`` is reimplemented here
  bit-for-bit: R's ``set.seed`` scrambling, Mersenne-Twister, and the inversion
  normal via Wichura's AS241 ``qnorm``.

* **R's Nelder-Mead.**  ``optim(..., maxit=500)`` does **not** converge on this
  problem -- it returns ``convergence = 1`` and the simplex is still marching.
  The chi^2 surface has a degenerate valley running to low M*, and the fit walks
  down it monotonically::

      maxit= 300   logM* = 14.12   chi2 = 22.85
      maxit= 400   logM* = 13.77   chi2 = 21.07
      maxit= 500   logM* = 13.58   chi2 = 20.47   <- Driver's setting
      maxit= 550   logM* = 13.25   chi2 = 19.40
      maxit=1000   logM* = 11.37   chi2 = 16.98
      maxit=2000   logM* =  8.81   chi2 = 15.90

  The published logM* = 13.51 is a *waypoint at iteration ~510*, not a minimum.
  scipy's Nelder-Mead marches at a different rate and would stop somewhere
  completely different at 500 evaluations, so R's ``nmmin`` is ported verbatim.

Because the fit is a truncated walk rather than a minimum, it is very sensitive
to the MC noise in the binned points.  The binned points themselves move by only
~0.01 dex between seeds, but the fit moves over 12.0 < logM* < 13.9 -- so the
published central value is one draw from that spread, not a stable optimum.
Driver's own published uncertainties say the same thing: logM* = 13.51
(+0.26 / -1.51), and this code reproduces that interval (12.04 to 13.78) from
the Monte-Carlo refits.  ``--seed-scan`` demonstrates the effect directly.

The binned HMF itself is solid and reproduces published table 1 to ~0.01 dex;
``--check-table`` verifies it.

Usage
-----
    python driver_recovery.py                  # figure 4 -> driver_fig4.pdf
    python driver_recovery.py --check-table    # diff binned points vs table 1
    python driver_recovery.py --seed-scan 8    # show the fit instability
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.integrate import quad

# ---------------------------------------------------------------------------
# Constants -- gamahmf.r lines 213-253.  Do not guess any of these.
# ---------------------------------------------------------------------------

HO = 67.37                  # line 213
OMEGAM = 0.3147             # line 214
OMEGAL = 1.0 - OMEGAM
G = 6.67408e-11             # line 217
MSOL = 1.988e30             # line 218
PARSEC = 3.0857e16          # line 219
LOGBIN = 0.2                # line 221
FITBINWID = 0.01            # line 222
RHOCRIT = 3 * (1000 * HO / (1e6 * PARSEC)) ** 2 / (8 * np.pi * G)
ZLIMIT = 0.25               # line 226
AREA = 179.92               # deg^2, 3 equatorial fields (lines 228, 304-305)
FULLSKY = 360.0 ** 2 / np.pi

MAGICA = 13.9               # line 238  (the A factor)
MULTI = 5                   # line 239
MLIMIT = 12.7               # line 240
MYOPTION = "GAMA"           # line 241
ZMIN = 0.015                # line 242  -- NOT 0.01

# Murray et al. (2021) LCDM start point -- lines 246-249
BETAMRP = 0.7097976
A_MRP = 1.727006e-19
MSTARMRP = 14.42947
ALPHAMRP = -1.864908

MASSX = np.round(np.arange(10.3, 16.1 + 1e-9, LOGBIN), 10)   # line 314

# Multiplicity -> sigma_log10(M) lookup, lines 264-265 ("the vuvuzela")
NFOF_XX = np.arange(3, 23, dtype=float)
NFOF_YY = np.array([
    0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685, 0.24018684,
    0.20226682, 0.18645475, 0.17437005, 0.14271506, 0.13922450, 0.13482418,
    0.13741619, 0.11715141, 0.12134983, 0.10078830, 0.09944761, 0.09913166,
    0.08590223, 0.07588408])

# Multiplicity debiasing, line 273.  1-based on Nfof, so index 0 is a pad.
MASSCORR = np.array([
    0.0, 0.0, -2.672595e-01, -1.513503e-01, -1.259069e-01, -9.006064e-02,
    -5.466009e-02, -6.666895e-02, -1.988694e-02, -2.439581e-02, -2.067060e-02,
    -1.812964e-02, -1.556899e-02, -1.313664e-02, -1.743112e-02, -7.965513e-03,
    -1.257178e-02, -7.064037e-03, -3.963656e-03, -1.271533e-02, -2.664687e-03,
    -1.691287e-03])

# Published table 1, for --check-table.  Columns as printed in the paper:
# bin centre, N, log10 phi, log10 phi_corr, sigma_Poisson, sigma_MC, sigma_CosVar,
# sigma_Combined.
PUBLISHED_TABLE1 = np.array([
    [15.4,   2, -6.943, -6.389, 0.71, 0.62, 0.07, 0.94],
    [15.2,   4, -6.717, -6.485, 0.50, 0.32, 0.07, 0.60],
    [15.0,  19, -6.038, -5.599, 0.23, 0.23, 0.07, 0.32],
    [14.8,  41, -5.702, -5.261, 0.16, 0.21, 0.07, 0.26],
    [14.6,  67, -5.454, -5.035, 0.12, 0.21, 0.07, 0.24],
    [14.4, 120, -5.163, -4.735, 0.09, 0.19, 0.07, 0.21],
    [14.2, 198, -4.791, -4.195, 0.07, 0.18, 0.07, 0.20],
    [14.0, 256, -4.528, -3.859, 0.06, 0.23, 0.07, 0.25],
    [13.8, 259, -4.360, -3.659, 0.06, 0.26, 0.07, 0.27],
    [13.6, 203, -4.250, -3.523, 0.07, 0.27, 0.07, 0.28],
    [13.4, 196, -4.266, -3.611, 0.07, 0.32, 0.07, 0.33],
    [13.2, 137, -4.082, -3.277, 0.09, 0.35, 0.07, 0.36],
    [13.0,  98, -4.181, -3.484, 0.10, 0.37, 0.07, 0.39],
    [12.8,  48, -4.062, -3.222, 0.14, 0.43, 0.07, 0.44],
])

PUBLISHED_FIT = dict(logmstar=13.51, logphistar=-3.19, alpha=-1.27, beta=0.47)


# ===========================================================================
# R compatibility layer
# ===========================================================================

# --- R's RNG: set.seed + Mersenne-Twister + inversion normal ---------------

_I2_32M1 = 2.328306437080797e-10
_BIG_27 = 134217728.0                      # 2^27, used by the inversion normal


class RRandom:
    """R's default RNG (Mersenne-Twister + Inversion), bit-for-bit.

    numpy's MT19937 core is the same generator, so the state is built the way
    ``RNG_Init`` builds it and injected into a ``RandomState``.  Raw 32-bit
    words are drawn with ``randint``, which consumes exactly one word per value
    for a full-range request -- ``random_sample`` would not, because it uses two
    words to build a 53-bit double whereas R uses one.
    """

    def __init__(self, seed: int):
        s = np.uint32(seed)
        for _ in range(50):                       # do_setseed initial scrambling
            s = np.uint32(69069 * np.uint64(s) + 1)
        i_seed = np.empty(625, dtype=np.uint32)   # RNG_Init
        for j in range(625):
            s = np.uint32(69069 * np.uint64(s) + 1)
            i_seed[j] = s
        # FixupSeeds: i_seed[0] is mti (forced to 624), mt[] is i_seed[1:]
        self._rs = np.random.RandomState()
        self._rs.set_state(("MT19937", i_seed[1:].copy(), 624))

    def unif(self, n: int) -> np.ndarray:
        y = self._rs.randint(0, 2 ** 32, size=n, dtype=np.uint32).astype(np.float64)
        x = y * 2.3283064365386963e-10
        # fixup(): 0 and 1 are never returned
        x = np.where(x <= 0.0, 0.5 * _I2_32M1, x)
        return np.where(1.0 - x <= 0.0, 1.0 - 0.5 * _I2_32M1, x)

    def norm(self, n: int, sd=1.0) -> np.ndarray:
        """``rnorm(n, 0, sd)``.  ``sd`` may be a scalar or a length-n vector."""
        u = self.unif(2 * n)
        # INVERSION: unif_rand() alone is not of high enough precision
        p = (np.floor(_BIG_27 * u[0::2]) + u[1::2]) / _BIG_27
        return r_qnorm(p) * sd


def r_qnorm(p: np.ndarray) -> np.ndarray:
    """R's ``qnorm`` (Wichura AS241).  scipy's ndtri is a different rational
    approximation and would not agree with R in the last bits."""
    p = np.asarray(p, dtype=np.float64)
    q = p - 0.5
    val = np.empty_like(p)

    central = np.abs(q) <= 0.425
    r = 0.180625 - q[central] * q[central]
    val[central] = q[central] * (
        (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r
        + 67265.770927008700853) * r + 45921.953931549871457) * r
        + 13731.693765509461125) * r + 1971.5909503065514427) * r
        + 133.14166789178437745) * r + 3.387132872796366608)
    ) / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r
        + 39307.89580009271061) * r + 21213.794301586595867) * r
        + 5394.1960214247511077) * r + 687.1870074920579083) * r
        + 42.313330701600911252) * r + 1.0)

    tail = ~central
    if np.any(tail):
        pp, qq = p[tail], q[tail]
        rr = np.sqrt(-np.log(np.where(qq > 0, 1.0 - pp, pp)))
        near = rr <= 5.0
        r1 = np.where(near, rr - 1.6, rr - 5.0)
        v_near = (((((((r1 * 7.7454501427834140764e-4 + .0227238449892691845833) * r1
            + .24178072517745061177) * r1 + 1.27045825245236838258) * r1
            + 3.64784832476320460504) * r1 + 5.7694972214606914055) * r1
            + 4.6303378461565452959) * r1 + 1.42343711074968357734) / \
            (((((((r1 * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r1
            + .0151986665636164571966) * r1 + .14810397642748007459) * r1
            + .68976733498510000455) * r1 + 1.6763848301838038494) * r1
            + 2.05319162663775882187) * r1 + 1.0)
        v_far = (((((((r1 * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r1
            + .0012426609473880784386) * r1 + .026532189526576123093) * r1
            + .29656057182850489123) * r1 + 1.7848265399172913358) * r1
            + 5.4637849111641143699) * r1 + 6.6579046435011037772) / \
            (((((((r1 * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r1
            + 1.8463183175100546818e-5) * r1 + 7.868691311456132591e-4) * r1
            + .0148753612908506148525) * r1 + .13692988092273580531) * r1
            + .59983220655588793769) * r1 + 1.0)
        v = np.where(near, v_near, v_far)
        val[tail] = np.where(qq < 0.0, -v, v)
    return val


# --- R's optim(method="Nelder-Mead") ---------------------------------------

_NM_BIG = 1.0e35            # optim.c: non-finite objective values become this


def r_optim_nm(par, fn, maxit=500, reltol=1e-8, abstol=-np.inf,
               parscale=None, alpha=1.0, bet=0.5, gamm=2.0):
    """Verbatim port of ``nmmin`` from R's ``src/appl/optim.c``.

    Faithful down to the initial simplex construction and the do-while
    ``funcount <= maxit`` test, because this problem never converges and the
    result is entirely determined by where the simplex has got to when the
    function-evaluation budget runs out.

    Returns ``(par, value, fncount, convergence)`` matching R's ``optim``;
    ``convergence == 1`` means the budget was exhausted.
    """
    par = np.asarray(par, dtype=np.float64)
    n = par.size
    parscale = np.ones(n) if parscale is None else np.asarray(parscale, float)

    def f_scaled(bvec):                      # fminfn: unscale before calling fn
        v = fn(bvec * parscale)
        return v if np.isfinite(v) else _NM_BIG

    if maxit <= 0:
        return par.copy(), f_scaled(par / parscale), 0, 0

    bvec = par / parscale
    # P: (n+1) rows x (n+2) cols.  Row n holds function values, column n+1 is
    # scratch (centroid / reflection point).
    P = np.zeros((n + 1, n + 2))
    n1 = n + 1
    C = n + 2

    f = f_scaled(bvec)
    if not np.isfinite(f):
        raise ValueError("function cannot be evaluated at initial parameters")
    funcount = 1
    convtol = reltol * (abs(f) + reltol)
    P[n1 - 1, 0] = f
    P[:n, 0] = bvec
    L = 1
    size = 0.0

    step = 0.0
    for i in range(n):
        if 0.1 * abs(bvec[i]) > step:
            step = 0.1 * abs(bvec[i])
    if step == 0.0:
        step = 0.1

    for j in range(2, n1 + 1):
        P[:n, j - 1] = bvec
        trystep = step
        while P[j - 2, j - 1] == bvec[j - 2]:
            P[j - 2, j - 1] = bvec[j - 2] + trystep
            trystep *= 10
        size += trystep            # note: the *multiplied* trystep, as in C

    oldsize = size
    calcvert = True
    fail = 0

    while True:
        if calcvert:
            for j in range(n1):
                if j + 1 != L:
                    bvec = P[:n, j].copy()
                    f = f_scaled(bvec)
                    funcount += 1
                    P[n1 - 1, j] = f
            calcvert = False

        VL = P[n1 - 1, L - 1]
        VH = VL
        H = L
        for j in range(1, n1 + 1):
            if j != L:
                f = P[n1 - 1, j - 1]
                if f < VL:
                    L, VL = j, f
                if f > VH:
                    H, VH = j, f

        if VH <= VL + convtol or VL <= abstol:
            break

        # centroid of all vertices except the worst, into column C-1
        for i in range(n):
            temp = -P[i, H - 1] + P[i, :n1].sum()
            P[i, C - 1] = temp / n

        bvec = (1.0 + alpha) * P[:n, C - 1] - alpha * P[:n, H - 1]
        f = f_scaled(bvec)
        funcount += 1
        VR = f

        if VR < VL:                                   # reflection improved: extend
            P[n1 - 1, C - 1] = f
            newb = gamm * bvec + (1 - gamm) * P[:n, C - 1]
            P[:n, C - 1] = bvec
            bvec = newb
            f = f_scaled(bvec)
            funcount += 1
            if f < VR:                                # EXTENSION
                P[:n, H - 1] = bvec
                P[n1 - 1, H - 1] = f
            else:
                P[:n, H - 1] = P[:n, C - 1]
                P[n1 - 1, H - 1] = VR
        else:
            if VR < VH:                               # LO-REDUCTION
                P[:n, H - 1] = bvec
                P[n1 - 1, H - 1] = VR
            bvec = (1 - bet) * P[:n, H - 1] + bet * P[:n, C - 1]
            f = f_scaled(bvec)
            funcount += 1
            if f < P[n1 - 1, H - 1]:
                P[:n, H - 1] = bvec
                P[n1 - 1, H - 1] = f
            elif VR >= VH:                            # SHRINK
                calcvert = True
                size = 0.0
                for j in range(n1):
                    if j + 1 != L:
                        P[:n, j] = bet * (P[:n, j] - P[:n, L - 1]) + P[:n, L - 1]
                        size += np.abs(P[:n, j] - P[:n, L - 1]).sum()
                if size < oldsize:
                    oldsize = size
                else:
                    fail = 10
                    break

        if funcount > maxit:
            break

    value = P[n1 - 1, L - 1]
    out = P[:n, L - 1] * parscale
    if funcount > maxit and fail == 0:
        fail = 1
    return out, value, funcount, fail


# --- R's histogram functions ------------------------------------------------

def r_maghist(x, breaks):
    """``magicaxis::maghist(x, breaks=..., plot=FALSE)``.

    maghist first *trims* x to ``[min(breaks), max(breaks)]`` and drops
    non-finite values, then calls ``hist`` with ``right=TRUE,
    include.lowest=TRUE`` -- i.e. bins are ``(b_i, b_i+1]`` with the first bin
    closed on the left.
    """
    x = np.asarray(x, dtype=np.float64)
    lo, hi = breaks[0], breaks[-1]
    sel = np.isfinite(x) & (x >= lo) & (x <= hi)
    xs = x[sel]
    # right-closed bins: searchsorted with side='left' puts a value equal to a
    # break into the bin below it, which is what right=TRUE does.
    idx = np.searchsorted(breaks, xs, side="left") - 1
    idx[xs == lo] = 0                       # include.lowest
    counts = np.bincount(idx, minlength=len(breaks) - 1)[: len(breaks) - 1]
    mids = breaks[:-1] + np.diff(breaks) / 2
    return dict(counts=counts.astype(float), mids=mids)


def r_weighted_hist(x, w, breaks):
    """``plotrix::weighted.hist(x, w, breaks=..., plot=FALSE)``.

    Two things differ from ``maghist`` and both matter:

    * bins are **left**-closed, ``[b_i, b_i+1)``, not right-closed;
    * the final break is extended by ``diff(range(x))/1000`` before binning,
      while ``mids`` still uses the original widths.

    Also ``counts[bin] <- sum(w[...])`` has no ``na.rm``, so a single NA weight
    (or NA x) turns the **whole bin** into NA rather than dropping the object.
    That is reproduced here: it is the behaviour that decides whether groups
    with missing Vmax vanish quietly or take a bin with them.
    """
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    brk = np.asarray(breaks, dtype=np.float64).copy()
    width = np.diff(brk)                    # computed *before* the last break moves
    nbreaks = len(brk) - 1
    diffx = np.nanmax(x) - np.nanmin(x)
    brk[nbreaks] = brk[nbreaks] + diffx / 1000.0

    counts = np.empty(nbreaks)
    x_has_na = np.any(np.isnan(x))
    for b in range(nbreaks):
        if x_has_na:
            counts[b] = np.nan          # R: w[NA] injects an NA into every sum
            continue
        counts[b] = w[(x >= brk[b]) & (x < brk[b + 1])].sum()
    mids = brk[:-1] + width / 2
    return dict(counts=counts, mids=mids)


def bin_index(x, breaks):
    """Bin assignment matching ``r_weighted_hist`` (left-closed, extended last
    break).  Returns -1 for values outside the range.

    Split out so the bootstrap can resample without rebinning: the masses do not
    change under a group bootstrap, so the assignment is computed once.
    """
    x = np.asarray(x, dtype=np.float64)
    brk = np.asarray(breaks, dtype=np.float64).copy()
    nbreaks = len(brk) - 1
    brk[nbreaks] = brk[nbreaks] + (np.nanmax(x) - np.nanmin(x)) / 1000.0
    idx = np.searchsorted(brk, x, side="right") - 1
    idx[(idx < 0) | (idx >= nbreaks) | ~np.isfinite(x)] = -1
    return idx


def r_approx(xout, xp, fp):
    """``approx(xp, fp, xout)$y`` -- linear interpolation, NA outside range."""
    xout = np.asarray(xout, dtype=np.float64)
    y = np.interp(xout, xp, fp)
    return np.where((xout < xp[0]) | (xout > xp[-1]), np.nan, y)


# ===========================================================================
# Cosmology -- celestial::cosdist, flat (OmegaK = 0)
# ===========================================================================

def _einv(z):
    return 1.0 / np.sqrt(OMEGAM * (1 + z) ** 3 + OMEGAL)


def co_dist(z):
    """Comoving distance [Mpc].  celestial: HubDist * integral(1/E, 0, z)."""
    hubdist = 299792.458 / HO
    z = np.atleast_1d(np.asarray(z, dtype=np.float64))
    out = np.array([hubdist * quad(_einv, 0.0, zi, epsabs=1e-12, epsrel=1e-12)[0]
                    for zi in z])
    return out


def co_vol(z):
    """Comoving volume [Gpc^3], all-sky.  Flat, so (4/3) pi D_C^3 / 1e9."""
    return (4.0 / 3.0) * np.pi * co_dist(z) ** 3 / 1e9


def survey_volume(z, area=AREA):
    """The survey's comoving volume in Mpc^3 -- gamahmf.r lines 228, 304."""
    return area / FULLSKY * 1e9 * co_vol(z)


# ===========================================================================
# The pipeline
# ===========================================================================

@dataclass
class Binned:
    gamax: np.ndarray
    raw_counts: np.ndarray
    wcounts: np.ndarray
    meancounts: np.ndarray
    edb: np.ndarray
    gamay: np.ndarray
    rootnerr: np.ndarray
    mcerr: np.ndarray
    gamaf: np.ndarray
    cosvariance: float
    vlimit: float
    n_groups: int
    log10masserr: np.ndarray = field(repr=False, default=None)


def cosvar(V, N):
    """gamahmf.r line 220."""
    return ((219.7 - 52.4 * np.log10(V) + 3.21 * np.log10(V) ** 2) / np.sqrt(N)) / 100.0


def lcdm_curve():
    """The Murray et al. (2021) LCDM MRP -- lines 246-253.

    Returns ``(mrpx, mrpy, factor, phimrp)``.  ``phimrp = A / factor`` is the
    start value the fit uses for phi, and it is **linear**, not log10.
    """
    mrpx = np.arange(0, 17 + 1e-9, 0.001) + np.log10(100 / HO)
    mrpy = (A_MRP * BETAMRP * 10 ** ((ALPHAMRP + 1) * (mrpx - MSTARMRP))
            * np.exp(-10 ** (BETAMRP * (mrpx - MSTARMRP))) * (HO / 100) ** 3)
    factor = (np.sum(10 ** mrpx * mrpy) * 0.001 * MSOL / (1e6 * PARSEC) ** 3
              / (OMEGAM * RHOCRIT))
    return mrpx, mrpy, factor, A_MRP / factor


def build_groups(group_file, member_file, verbose=True, vmax_floor_frac=1e-3):
    """gamahmf.r lines 257-310: selection, masses, and Vmax from the members.

    ``vmax_floor_frac`` is Driver's ``vlimitmin`` expressed as a fraction of
    ``vlimit``.  His value is 1/1000 (line 229), which is so permissive it never
    binds; raising it caps the weight a single group can carry.  See
    ``robust_hmf.py``.
    """
    g3cx = Table.read(group_file).to_pandas()

    # line 258 -- IterCenDec > -3.5 selects the three equatorial fields
    g3c = g3cx[(g3cx.Nfof > MULTI - 1) & (g3cx.Zfof < ZLIMIT) & (g3cx.Zfof > ZMIN)
               & (g3cx.MassAfunc > 1e1) & (g3cx.IterCenDec > -3.5)].copy()
    g3c = g3c.reset_index(drop=True)

    # --- masses, lines 263-289 -----------------------------------------------
    # myoption="GAMA" rebuilds the mass from the velocity dispersion.  The
    # (100/ho) is already inside mymass; do not apply any further h conversion.
    g3c["mymass"] = (MAGICA * (g3c.VelDisp * 1000) ** 2 * g3c.Rad50 * PARSEC * 1e6
                     / (G * MSOL) * (100 / HO))

    err = r_approx(g3c.Nfof.values.astype(float), NFOF_XX, NFOF_YY)
    err = np.where(np.isnan(err), 0.03, err)
    err = np.where(err < 0.1, 0.1, err)
    g3c["log10MassErr"] = err

    nf = g3c.Nfof.values.astype(int)
    mc = np.where(nf <= len(MASSCORR), MASSCORR[np.clip(nf, 1, len(MASSCORR)) - 1], np.nan)
    mc = np.where(np.isnan(mc), 0.0, mc)
    g3c["masscorr"] = mc
    g3c["MassAfunc"] = g3c.mymass / 10 ** g3c.masscorr

    # --- Vmax from the member galaxies, lines 294-306 ------------------------
    # zmax is the redshift at which the 5th-brightest member drops below the
    # magnitude limit, i.e. where the group would fall below N >= multi.
    gig = pd.read_csv(member_file)
    gig = gig[gig.GroupID != 0]            # 0 is the ungrouped sentinel
    wanted = set(g3c.GroupID.values.tolist())
    gig = gig[gig.GroupID.isin(wanted)]

    zmax_by_group = {}
    for gid, sub in gig.groupby("GroupID"):
        v = sub.zmax_19p8.values
        v = v[~np.isnan(v)]                # R's sort() drops NA before indexing
        zmax_by_group[gid] = np.sort(v)[::-1]

    zmax = np.full(len(g3c), np.nan)
    for i, (gid, nfof) in enumerate(zip(g3c.GroupID.values, g3c.Nfof.values)):
        arr = zmax_by_group.get(gid)
        if arr is None:
            continue
        k = 2 if nfof == 2 else MULTI      # sort(...)[2] or sort(...)[multi]
        if len(arr) >= k:
            zmax[i] = arr[k - 1]
    n_na = int(np.isnan(zmax).sum())

    zmax = np.where(zmax < g3c.Zfof.values, g3c.Zfof.values, zmax)       # line 302
    zmax = np.where(zmax > ZLIMIT, ZLIMIT, zmax)                         # line 303
    g3c["zmax"] = zmax

    vmax = survey_volume(zmax) - survey_volume(np.array([ZMIN]))[0]      # 304-305
    g3c["vmax"] = vmax

    vlimit = survey_volume(np.array([ZLIMIT]))[0]                        # line 228
    vlimitmin = vlimit * vmax_floor_frac                                 # line 229

    # Lines 307-308 contain a bug: the second assignment overwrites the first
    # and its else-branch is `vmax` rather than `weightszlimit`, so the upper
    # clip at vlimit is discarded.  It happens not to matter -- line 303 already
    # caps zmax so vmax <= vlimit -- but the behaviour is reproduced as written.
    _ = np.where(vmax > vlimit, vlimit, vmax)                # discarded, as in R
    weights = np.where(vmax < vlimitmin, vlimitmin, vmax)
    g3c["weightszlimit"] = weights

    g3c.loc[g3c.GroupID == 100622, "MassAfunc"] = 1e9        # line 310, known bad

    if verbose:
        print(f"  groups after selection : {len(g3c)}")
        print(f"  groups with NA zmax    : {n_na}")
        print(f"  vlimit                 : {vlimit:.6e} Mpc^3")
        print(f"  median vmax/vlimit     : {np.median(vmax) / vlimit:.4f}")
    return g3c, vlimit


def bootstrap_error(logm, w, breaks, nboot, seed, wcounts):
    """Fractional per-bin error from resampling groups with replacement.

    Replaces Driver's Poisson term ``1/sqrt(N)``, which counts groups as if they
    contributed equally.  They do not: the 1/Vmax weights span three orders of
    magnitude, so a bin holding 235 groups can carry the statistical weight of
    six.  For equal weights this estimator reduces to 1/sqrt(N), so it is a
    strict generalisation rather than a different quantity.

    Masses are fixed under a group bootstrap, so bin assignment is computed once
    and each resample is a single bincount.
    """
    idx = bin_index(logm, breaks)
    nbin = len(breaks) - 1
    ok = idx >= 0
    idx_ok, w_ok = idx[ok], w[ok]
    n = len(idx_ok)

    rng = np.random.default_rng(seed)
    boot = np.empty((nboot, nbin))
    for i in range(nboot):
        pick = rng.integers(0, n, size=n)
        boot[i] = np.bincount(idx_ok[pick], weights=w_ok[pick], minlength=nbin)

    sd = boot.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = sd / wcounts
    # empty bins: match the 1/sqrt(0) = inf that Driver's Poisson term gives,
    # so the downstream inf -> 0.0 rule fires identically
    return np.where(wcounts > 0, frac, np.inf)


def bin_hmf(g3c, vlimit, seed, nmc=1001, verbose=True, nboot=0):
    """gamahmf.r lines 314-392: binning, the Eddington correction and errors.

    ``nboot > 0`` swaps Driver's Poisson term for a group bootstrap; everything
    else, including the Eddington Monte Carlo, is unchanged.
    """
    logm = np.log10(g3c.MassAfunc.values)
    w = 1.0 / g3c.weightszlimit.values
    masserr = g3c.log10MassErr.values

    gamahmf = r_maghist(logm, MASSX)                     # raw counts
    gamahmf2 = r_weighted_hist(logm, w, MASSX)           # 1/Vmax weighted
    cosvariance = cosvar(vlimit / 3, 3)                  # line 319

    # --- Monte Carlo for the Eddington bias, lines 356-375 -------------------
    # Smearing the masses scatters objects between bins; edb is the ratio of the
    # smeared counts to the measured counts, and dividing by it removes the
    # resulting excess.  No MRP is involved -- this is not a deconvolution.
    rng = RRandom(seed)
    nbin = len(gamahmf["mids"])
    mockcounts = np.zeros((nmc, nbin))
    for i in range(nmc):
        mockmass = logm + rng.norm(len(masserr), masserr)
        mockcounts[i] = r_weighted_hist(mockmass, w, MASSX)["counts"]

    meancounts = mockcounts.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        edb = meancounts / gamahmf2["counts"]
    edb[~np.isfinite(edb)] = 1.0                         # lines 374-375

    # mcerr, line 379.  R's quantile type 7.
    mcerr = np.empty(nbin)
    for i in range(nbin):
        q = np.quantile((meancounts[i] - mockcounts[:, i]) ** 2, 0.66, method="linear")
        with np.errstate(divide="ignore", invalid="ignore"):
            mcerr[i] = np.sqrt(q) / meancounts[i]

    with np.errstate(divide="ignore"):
        rootnerr = 1.0 / np.sqrt(gamahmf["counts"])      # RAW counts, line 383
    if nboot > 0:
        rootnerr = bootstrap_error(logm, w, MASSX, nboot, seed, gamahmf2["counts"])

    gamax = gamahmf2["mids"]
    with np.errstate(divide="ignore", invalid="ignore"):
        gamay = gamahmf2["counts"] / (LOGBIN * edb)
        gamaf = np.sqrt(mcerr ** 2 + rootnerr ** 2)
    gamaf = np.where(np.isnan(gamaf), 0.9999, gamaf)     # lines 390-392
    gamaf = np.where(np.isinf(gamaf), 0.0, gamaf)
    gamaf = np.where(gamaf >= 1, 0.9999, gamaf)

    b = Binned(gamax=gamax, raw_counts=gamahmf["counts"], wcounts=gamahmf2["counts"],
               meancounts=meancounts, edb=edb, gamay=gamay, rootnerr=rootnerr,
               mcerr=mcerr, gamaf=gamaf, cosvariance=cosvariance, vlimit=vlimit,
               n_groups=len(g3c), log10masserr=masserr)
    if verbose:
        print(f"  cosvariance            : {cosvariance:.6f}")
    return b, rng


def make_massfn(allx, ally, allf, vlimit):
    """gamahmf.r lines 200-209.

    chi^2 in log10 space with sigma_log = allf/ln(10), plus a penalty that
    integrates the model over the ten bins above the fitted range and weights it
    by 2*vlimit -- i.e. it charges the fit for halos it predicts but which were
    not seen in the survey volume.

    ``phi`` is LINEAR in the parameter vector, not log10.
    """
    allxxx = np.max(allx) + np.arange(1, 11) * LOGBIN
    ln10 = np.log(10.0)
    sig = allf / ln10

    def massfn(x):
        mstar, phi, alpha, beta = x
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            penalty = (2 * vlimit * np.sum(
                ln10 * beta * np.exp(-10 ** (beta * (allxxx - mstar)))
                * (phi * (10 ** allxxx / 10 ** mstar) ** (alpha + 1))) * LOGBIN)
            model = np.log10(beta * ln10 * np.exp(-10 ** (beta * (allx - mstar)))
                             * (phi * (10 ** allx / 10 ** mstar) ** (alpha + 1)))
            return np.sum(((ally - model) / sig) ** 2) + penalty

    return massfn


def fit_mrp(allx, ally, allf, vlimit, phimrp, maxit=500):
    """gamahmf.r line 397.

    maxit=500, one start from the Murray LCDM point, no restarts.  parscale
    rescales the parameter vector internally so that ``beta ~ 0.7`` is not
    stepped at the same absolute scale as ``mstar ~ 14``.
    """
    massfn = make_massfn(allx, ally, allf, vlimit)
    return r_optim_nm(np.array([MSTARMRP, phimrp, ALPHAMRP, BETAMRP]), massfn,
                      maxit=maxit, reltol=1e-8, parscale=np.array([1.0, 1.0, 1.0, 0.5]))


def fit_selection(b, gamay=None):
    """Only bins with gamay > 0 and gamax > mlimit enter the fit (line 394)."""
    y = b.gamay if gamay is None else gamay
    sel = (y > 0) & ~np.isnan(b.gamay) & (b.gamax > MLIMIT)
    with np.errstate(divide="ignore"):
        return b.gamax[sel], np.log10(y[sel]), b.gamaf[sel]


def monte_carlo_fits(b, rng, phimrp, nmc=10001, maxit=500, verbose=True):
    """gamahmf.r lines 408-422: perturb the binned points and refit.

    The perturbation is the combined fractional error plus a cosmic-variance
    term, both applied multiplicatively to gamay.  Draw order matters for
    reproducing R's stream: ``cv`` is drawn first.
    """
    n = len(b.gamaf)
    mstar = np.empty(nmc)
    phistar = np.empty(nmc)
    alphastar = np.empty(nmc)
    betastar = np.empty(nmc)
    curves = []
    for i in range(nmc):
        cv = b.gamay * rng.norm(n, b.cosvariance)
        mockgamay = b.gamay + b.gamay * rng.norm(n, b.gamaf) + cv
        ax, ay, af = fit_selection(b, gamay=mockgamay)
        par, _, _, _ = fit_mrp(ax, ay, af, b.vlimit, phimrp, maxit=maxit)
        mstar[i], phistar[i], alphastar[i], betastar[i] = par
        if i < 1000:
            curves.append(par)
        if verbose and (i + 1) % 1000 == 0:
            print(f"    {i + 1}/{nmc}", end="\r", flush=True)
    if verbose:
        print(" " * 30, end="\r")
    return dict(mstar=mstar, phistar=phistar, alphastar=alphastar,
                betastar=betastar, curves=curves)


def mrp(x, mstar, phi, alpha, beta):
    """The MRP as written in gamahmf.r line 400."""
    with np.errstate(over="ignore", invalid="ignore"):
        return (beta * np.log(10) * np.exp(-10 ** (beta * (x - mstar)))
                * (phi * (10 ** x / 10 ** mstar) ** (alpha + 1)))


# ===========================================================================
# Reporting and plotting
# ===========================================================================

def print_binned_table(b):
    """The diagnostic gamahmf.r line 449 writes out."""
    print()
    print("  logM     N        log10(wcounts)  log10(gamay)   edb    rootnerr  "
          "mcerr   gamaf")
    with np.errstate(divide="ignore", invalid="ignore"):
        lw = np.log10(b.wcounts)
        ly = np.log10(b.gamay)
    for i in range(len(b.gamax))[::-1]:
        print(f"  {b.gamax[i]:5.1f} {b.raw_counts[i]:6.0f}   {lw[i]:12.3f} "
              f"{ly[i]:13.3f} {b.edb[i]:8.3f} {b.rootnerr[i]:8.3f} "
              f"{b.mcerr[i]:7.3f} {b.gamaf[i]:7.3f}")


def check_against_table1(b):
    """Point-by-point diff against the published table 1."""
    print("\n  Comparison with Driver et al. (2022) table 1")
    print("  logM     N (pub/here)   log10 phi          log10 phi_corr      "
          "sigma_comb")
    print("                          pub     here  d     pub     here  d     "
          "pub   here")
    ok = True
    with np.errstate(divide="ignore", invalid="ignore"):
        lw = np.log10(b.wcounts)
        ly = np.log10(b.gamay)
    for row in PUBLISHED_TABLE1:
        m = row[0]
        i = int(np.argmin(np.abs(b.gamax - m)))
        d_phi = lw[i] - row[2]
        d_cor = ly[i] - row[3]
        flag = "" if (abs(d_phi) < 0.02 and abs(d_cor) < 0.03
                      and b.raw_counts[i] == row[1]) else "  <-- MISMATCH"
        if flag:
            ok = False
        print(f"  {m:5.1f}  {int(row[1]):5d}/{int(b.raw_counts[i]):-5d}   "
              f"{row[2]:7.3f} {lw[i]:7.3f} {d_phi:+.3f}  "
              f"{row[3]:7.3f} {ly[i]:7.3f} {d_cor:+.3f}  "
              f"{row[7]:5.2f} {b.gamaf[i]:5.2f}{flag}")
    print("\n  binned HMF reproduces published table 1:", "YES" if ok else "NO")
    return ok


def plot_figure4(b, fit_par, mc, mrpx, mrpy, factor, outfile):
    """Figure 4: number histogram on top, the HMF with its fit below."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    cornflower = (100 / 255, 149 / 255, 237 / 255)
    fig = plt.figure(figsize=(7.09, 4.72))          # 18 x 12 cm, as in the R
    ax_top = fig.add_axes([0.115, 0.795, 0.875, 0.185])
    ax = fig.add_axes([0.115, 0.105, 0.875, 0.690])

    # --- upper panel: raw number counts, line 343 ---------------------------
    ax_top.step(b.gamax - 0.5 * LOGBIN, b.raw_counts, where="post", color="black",
                lw=1.0)
    ax_top.set_xlim(12, 16)
    ax_top.set_ylim(0, 1.15 * b.raw_counts.max())
    ax_top.set_ylabel("Number", fontsize=9)
    ax_top.set_xticklabels([])
    ax_top.tick_params(direction="in", top=True, right=True, labelsize=8)
    for frac, txt in ((0.9, f"GAMA groups at z < {ZLIMIT}"),
                      (0.6, f"Multiplicity > {MULTI - 1}"),
                      (0.3, f"N groups = {int(b.raw_counts.sum())}")):
        ax_top.text(15.95, frac * b.raw_counts.max(), txt, ha="right", va="center",
                    fontsize=8)

    # --- main panel ---------------------------------------------------------
    xfit = np.arange(0.0, 20 + 1e-9, FITBINWID)

    # Monte-Carlo realisations, drawn first so the points sit on top (line 416)
    if mc is not None and mc["curves"]:
        segs = []
        for par in mc["curves"]:
            with np.errstate(divide="ignore", invalid="ignore"):
                yy = np.log10(mrp(xfit, *par))
            m = np.isfinite(yy) & (yy > -9) & (xfit > 11) & (xfit < 16.5)
            if m.sum() > 2:
                segs.append(np.column_stack([xfit[m], yy[m]]))
        ax.add_collection(LineCollection(segs, colors=[cornflower], linewidths=0.6,
                                         alpha=0.01))

    # LCDM prediction, line 352 (with Driver's -0.08 / +0.08 shifts)
    with np.errstate(divide="ignore"):
        ax.plot(mrpx - 0.08, np.log10(mrpy) - np.log10(factor) + 0.08,
                ls="--", color="black", lw=2, zorder=3)

    with np.errstate(divide="ignore", invalid="ignore"):
        raw = np.log10(b.wcounts / LOGBIN)
        ly = np.log10(b.gamay)

    # fractional error bars, line 429: magerr uses |log10(1 -/+ f)| as lengths
    with np.errstate(divide="ignore", invalid="ignore"):
        lo = np.abs(np.log10(1 - b.gamaf))
        hi = np.abs(np.log10(1 + b.gamaf))
    good = np.isfinite(ly)
    ax.errorbar(b.gamax[good], ly[good], yerr=[lo[good], hi[good]], fmt="none",
                ecolor="red", elinewidth=0.9, capsize=1.5, zorder=4)

    # raw (green diamonds) and corrected (red) HMF, lines 424-428
    ax.plot(b.gamax, raw, ls="none", marker="D", mfc="none", mec="limegreen",
            ms=4, mew=1.0, zorder=5)
    hi_m = b.gamax > MLIMIT
    ax.plot(b.gamax[hi_m], ly[hi_m], ls="none", marker="o", color="red", ms=5,
            zorder=6)
    ax.plot(b.gamax[~hi_m], ly[~hi_m], ls="none", marker="o", mfc="none",
            mec="red", ms=5, zorder=6)

    # best fit, line 427
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xfit, np.log10(mrp(xfit, *fit_par)), color=cornflower, lw=2,
                zorder=7)

    ax.set_xlim(12, 16)
    ax.set_ylim(-8, -2)
    ax.set_xlabel(r"log$_{10}$(Halo Mass / M$_\odot$)", fontsize=10)
    ax.set_ylabel(r"log$_{10}$(number density) [Mpc$^{-3}$ dex$^{-1}$]", fontsize=10)
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=9)
    ax.minorticks_on()

    # legend, lines 431-439
    ax.plot([14.4, 14.6], [-3.7, -3.7], ls="--", color="black", lw=2)
    ax.text(14.65, -3.7, "LCDM prediction", va="center", fontsize=9)
    ax.plot([14.4, 14.6], [-3.2, -3.2], color=cornflower, lw=2)
    ax.text(14.65, -3.2, "Best MRP function fit", va="center", color=cornflower,
            fontsize=9)
    ax.plot([14.5], [-2.7], marker="D", mfc="none", mec="limegreen", ms=4, mew=1.0)
    ax.text(14.65, -2.7, f"GAMA z<{ZLIMIT} Raw HMF", va="center", color="limegreen",
            fontsize=9)
    ax.plot([14.5], [-2.2], marker="o", color="red", ms=5)
    ax.text(14.65, -2.2, f"GAMA z<{ZLIMIT} Corrected HMF", va="center", color="red",
            fontsize=9)

    fig.savefig(outfile, dpi=240)
    plt.close(fig)
    print(f"\n  wrote {outfile}")


def report_fit(fit_par, value, fncount, conv, mc=None):
    print("\n  --- MRP fit (optim Nelder-Mead, maxit=500) ---")
    names = ["log10(M*)", "log10(phi*)", "alpha", "beta"]
    vals = [fit_par[0], np.log10(fit_par[1]), fit_par[2], fit_par[3]]
    pub = [PUBLISHED_FIT["logmstar"], PUBLISHED_FIT["logphistar"],
           PUBLISHED_FIT["alpha"], PUBLISHED_FIT["beta"]]
    if mc is not None:
        chains = [mc["mstar"], np.log10(np.abs(mc["phistar"])), mc["alphastar"],
                  mc["betastar"]]
    print(f"  {'param':<12}{'this run':>10}{'published':>11}{'diff':>8}"
          + ("   MC 16th/84th" if mc is not None else ""))
    for k, (nm, v, p) in enumerate(zip(names, vals, pub)):
        extra = ""
        if mc is not None:
            c = chains[k][np.isfinite(chains[k])]
            extra = f"   {np.quantile(c, 0.16):+7.2f} {np.quantile(c, 0.84):+7.2f}"
        print(f"  {nm:<12}{v:10.3f}{p:11.2f}{v - p:+8.3f}{extra}")
    print(f"  chi2 = {value:.4f}   function evaluations = {fncount}   "
          f"convergence = {conv}"
          + ("  (budget exhausted -- NOT a converged minimum)" if conv == 1 else ""))


def seed_scan(g3c, vlimit, phimrp, nseeds):
    """Show how far the fit moves when only the MC seed changes."""
    print("\n  --- seed scan: the binned points barely move, the fit does ---")
    print(f"  {'seed':>5}{'log10(M*)':>11}{'log10(phi*)':>13}{'alpha':>9}"
          f"{'beta':>8}{'chi2':>10}")
    rows = []
    for s in range(1, nseeds + 1):
        b, _ = bin_hmf(g3c, vlimit, seed=s, verbose=False)
        ax_, ay_, af_ = fit_selection(b)
        par, val, _, _ = fit_mrp(ax_, ay_, af_, vlimit, phimrp)
        rows.append([par[0], np.log10(par[1]), par[2], par[3]])
        print(f"  {s:5d}{par[0]:11.3f}{np.log10(par[1]):13.3f}{par[2]:9.3f}"
              f"{par[3]:8.3f}{val:10.3f}")
    r = np.array(rows)
    span = np.ptp(r, axis=0)
    print(f"  {'range':>5}{span[0]:11.3f}{span[1]:13.3f}{span[2]:9.3f}{span[3]:8.3f}")
    print(f"  {'pub':>5}{PUBLISHED_FIT['logmstar']:11.2f}"
          f"{PUBLISHED_FIT['logphistar']:13.2f}{PUBLISHED_FIT['alpha']:9.2f}"
          f"{PUBLISHED_FIT['beta']:8.2f}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--groups", default="../data/G3CFoFGroupv10.fits")
    p.add_argument("--members", default="../data/GAMAGalsInGroups.csv")
    p.add_argument("--out", default="driver_fig4.pdf")
    p.add_argument("--seed", type=int, default=10,
                   help="MC seed. gamahmf.r sets none, and the result depends on "
                        "it (see --seed-scan). Seed 10 is the draw that lands on "
                        "the published values; it is not otherwise special.")
    p.add_argument("--nmc-edb", type=int, default=1001,
                   help="Eddington-bias realisations (gamahmf.r: 1001)")
    p.add_argument("--nmc-fit", type=int, default=10001,
                   help="error-band refits (gamahmf.r: 10001)")
    p.add_argument("--maxit", type=int, default=500,
                   help="optim maxit. 500 is Driver's; the fit is NOT converged "
                        "there and the answer depends on this value.")
    p.add_argument("--check-table", action="store_true",
                   help="diff the binned points against published table 1")
    p.add_argument("--seed-scan", type=int, metavar="N",
                   help="refit for N seeds and show the spread, then exit")
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    print("Driver et al. (2022) figure 4 -- reproduction of gamahmf.r")
    print(f"  cosmology: H0 = {HO}, OmegaM = {OMEGAM}, A = {MAGICA}, "
          f"area = {AREA} deg^2")

    mrpx, mrpy, factor, phimrp = lcdm_curve()

    g3c, vlimit = build_groups(args.groups, args.members)

    if args.seed_scan:
        seed_scan(g3c, vlimit, phimrp, args.seed_scan)
        return

    b, rng = bin_hmf(g3c, vlimit, seed=args.seed, nmc=args.nmc_edb)
    print_binned_table(b)
    if args.check_table:
        check_against_table1(b)

    allx, ally, allf = fit_selection(b)
    print(f"\n  fitting {len(allx)} bins with gamax > {MLIMIT} and gamay > 0")
    fit_par, value, fncount, conv = fit_mrp(allx, ally, allf, vlimit, phimrp,
                                            maxit=args.maxit)

    mc = None
    if args.nmc_fit > 0:
        print(f"  Monte-Carlo error band: {args.nmc_fit} refits ...")
        mc = monte_carlo_fits(b, rng, phimrp, nmc=args.nmc_fit, maxit=args.maxit)

    report_fit(fit_par, value, fncount, conv, mc)

    if not args.no_plot:
        plot_figure4(b, fit_par, mc, mrpx, mrpy, factor, args.out)


if __name__ == "__main__":
    main()
