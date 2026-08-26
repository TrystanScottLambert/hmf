#!/usr/bin/env python3
"""Survey figures for the group catalogues: GAMA cones, the SDSS wedge, the sky.

Three deliverables, sized for where they go in the paper:

    --cones    one single-column figure per GAMA region  (cone_g09.png, ...)
    --puck     the SDSS wedge, across both columns        (sdss_puck.png)
    --sky      GAMA and SDSS in one 360-degree polar plot (sky_polar.png)

All three are the same picture in different guises: right ascension as the
polar angle, comoving distance as the radius, declination flattened away.
Galaxies are black dots; groups are red circles **sized by multiplicity**.

Every wedge carries two radial scales -- comoving distance up the left edge and
redshift up the right -- so no in-panel arc labels are needed and there is only
one family of arcs, all drawn identically.

**The GAMA cones stretch the angular axis.**  A region spans ~12 deg of RA and
is drawn as a 65 deg wedge, otherwise the wedge is a needle.  The radial axis is
true comoving distance in every panel, so structure is never distorted radially.
The SDSS wedge and the 360-degree plot need no stretch and get none.

Groups are the **full released catalogues** at multiplicity >= 3 -- not the
multiplicity >= 5, z > 0.015 subset the mass function is measured from.  These
figures accompany the data release, so they show everything in it apart from
the HMF cuts, which are applied later.  Pairs (N = 2) are excluded: they are a
separate product (``get_pair_dmu`` in ``make_gama_dmu``), they are 65 per cent
of the GAMA catalogue and 61 per cent of the SDSS one, and drawing them buries
everything else.  ``--nmin`` changes the floor.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

import driver_recovery as dr
import plotting

COSMO = FlatLambdaCDM(H0=dr.HO, Om0=0.3147)

GAMA_GALS = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CGal.fits"
GAMA_GROUPS = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits"
SDSS_DIR = "/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS"
SDSS_GALS = f"{SDSS_DIR}/sdss_galaxies.parquet"
SDSS_GROUPS = f"{SDSS_DIR}/sdss_groups.parquet"

GAL_C = "k"
GRP_C = "#E03B2F"
GAMA_C = "#E03B2F"
SDSS_C = "#1F5FA8"

# Marker areas for the multiplicity key.  N**0.7 keeps a 758-member cluster
# from being 150 times the area of a 5-member group while still reading as a
# clear size sequence.
MULT_KEY = [3, 5, 10, 20, 50, 100]


def mult_size(n, scale=2.6, cap=150.0):
    return scale * np.minimum(np.asarray(n, float), cap) ** 0.7


def _fits_frame(path):
    t = Table.read(path)
    out = {}
    for c in t.colnames:
        a = np.asarray(t[c])
        if a.dtype.byteorder == ">":
            a = a.astype(a.dtype.newbyteorder("="))
        out[c] = a
    return pd.DataFrame(out)


def _wedge_xy(ra, dist, ra_mid, stretch, pa=90.0):
    """RA -> polar angle about the apex; ``pa`` is where the wedge points."""
    theta = np.radians((np.asarray(ra) - ra_mid) * stretch + pa)
    return dist * np.cos(theta), dist * np.sin(theta)


def frame_box(ax, x0, x1, y0, y1):
    """Pin the tight-bbox crop to an exact rectangle.

    ``axis("off")`` leaves no axes frame, so ``bbox_inches="tight"`` crops to
    the drawn artists and the output aspect ratio becomes whatever the labels
    happen to reach.  An invisible patch spanning the intended box makes the
    crop deterministic.
    """
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec="none",
                           zorder=0))


def mult_legend(ax, loc="upper left", title=r"$N_{\mathrm{fof}}$", fontsize=6.5):
    import matplotlib.pyplot as plt
    # mew matches the tick and outline weight so the key reads as part of the
    # same drawing rather than a lighter afterthought.
    handles = [plt.Line2D([], [], ls="", marker="o", mfc="none", mec=GRP_C,
                          mew=1.0, ms=np.sqrt(mult_size(n)), label=f"{n}")
               for n in MULT_KEY]
    handles[-1].set_label(f"$\\geq${MULT_KEY[-1]}")
    return ax.legend(handles=handles, loc=loc, frameon=False, title=title,
                     prop={"size": fontsize, "weight": "bold"},
                     title_fontproperties={"size": fontsize + 1.0,
                                           "weight": "bold"},
                     handletextpad=0.6, labelspacing=0.8, borderpad=0.3,
                     scatterpoints=1)


# ---------------------------------------------------------------------------
# The shared wedge frame: outline, arcs, and the two radial scales
# ---------------------------------------------------------------------------

def _edge_geometry(edge_dir, outward_sign):
    """Outward normal and text rotation for one straight edge of a wedge.

    A tick on a radial edge has to be drawn along the edge *normal*: stepping
    radially would just slide along the edge itself and draw nothing.  The
    label is then offset along the same normal so it sits clear of the data,
    and rotated to lie along the edge when the edge is nearer horizontal than
    vertical -- otherwise it stays upright, which is how a y axis reads.
    """
    normal = edge_dir + outward_sign * np.pi / 2.0
    rot = 0.0
    if abs(np.cos(edge_dir)) > abs(np.sin(edge_dir)):     # a shallow edge
        rot = np.degrees(edge_dir)
        rot -= 180.0 * np.round(rot / 180.0)
    return normal, rot


def draw_frame(ax, ra_lo, ra_hi, ra_mid, stretch, dmax, dticks, zticks,
               fontsize=7.5, ticklen=0.018, labelpad=7.0, lw=0.9):
    """Outline the wedge and hang the two radial scales off its straight edges.

    Left edge  -> comoving distance, at round values in Mpc, each with an arc.
    Right edge -> redshift, at round values, tick marks only.

    Two independent scales rather than one set of arcs labelled twice: that way
    both axes carry round numbers.  Every arc is drawn identically -- the only
    heavier lines are the wedge outline and the ticks.
    """
    half = np.radians((ra_hi - ra_lo) * stretch / 2.0)
    arc = np.linspace(ra_lo, ra_hi, 600)
    for r in (ra_lo, ra_hi):
        ax.plot(*_wedge_xy(np.array([r, r]), np.array([0.0, dmax]), ra_mid,
                           stretch), color="k", lw=lw, zorder=5)
    ax.plot(*_wedge_xy(arc, np.full_like(arc, dmax), ra_mid, stretch),
            color="k", lw=lw, zorder=5)

    for edge_ra, edge_dir, sign, values, arcs in (
            (ra_hi, np.pi / 2.0 + half, +1.0, dticks, True),
            (ra_lo, np.pi / 2.0 - half, -1.0, zticks, False)):
        normal, rot = _edge_geometry(edge_dir, sign)
        nx, ny = np.cos(normal), np.sin(normal)
        for v in values:
            d = float(COSMO.comoving_distance(v).value) if not arcs else float(v)
            if d > dmax * 1.001 or d <= 0:
                continue
            if arcs and d < dmax:
                ax.plot(*_wedge_xy(arc, np.full_like(arc, d), ra_mid, stretch),
                        color="k", lw=0.5, ls=(0, (4, 3)), alpha=0.45, zorder=5)
            x0, y0 = _wedge_xy(np.array([edge_ra]), np.array([d]), ra_mid, stretch)
            t = ticklen * dmax
            ax.plot([x0[0], x0[0] - nx * t], [y0[0], y0[0] - ny * t],
                    color="k", lw=lw, zorder=6)                  # inward tick
            ax.annotate(f"{v:g}", (x0[0], y0[0]), fontsize=fontsize,
                        ha="center", va="center", rotation=rot,
                        rotation_mode="anchor", fontweight="bold", zorder=7,
                        xytext=(labelpad * nx, labelpad * ny),
                        textcoords="offset points")


def _nice_step(span, target=3):
    """A round tick step giving roughly ``target`` ticks across ``span``."""
    cands = np.array([0.5, 1, 2, 2.5, 5, 10, 15, 20, 30, 45, 60])
    return float(cands[np.argmin(np.abs(span / cands - target))])


def ra_ticks(ax, ra_lo, ra_hi, ra_mid, stretch, dmax, step=None, target=3,
             fontsize=8.5, ticklen=0.018, labelpad=0.040, lw=0.9):
    """Right ascension ticked and labelled around the outer arc, pointing in."""
    step = step or _nice_step(ra_hi - ra_lo, target)
    first = np.ceil((ra_lo + 0.02 * (ra_hi - ra_lo)) / step) * step
    for ra_t in np.arange(first, ra_hi - 0.02 * (ra_hi - ra_lo), step):
        x0, y0 = _wedge_xy(np.array([ra_t]), np.array([dmax]), ra_mid, stretch)
        x1, y1 = _wedge_xy(np.array([ra_t]), np.array([dmax * (1 - ticklen)]),
                           ra_mid, stretch)
        ax.plot([x0[0], x1[0]], [y0[0], y1[0]], color="k", lw=lw, zorder=6)
        x2, y2 = _wedge_xy(np.array([ra_t]), np.array([dmax * (1 + labelpad)]),
                           ra_mid, stretch)
        ax.annotate(f"{ra_t:g}$^\\circ$", (x2[0], y2[0]), fontsize=fontsize,
                    ha="center", va="center", zorder=7, fontweight="bold")


def edge_titles(ax, ra_lo, ra_hi, ra_mid, stretch, dmax, fontsize=9.0,
                pad=0.10, mode="end", opening=None):
    """Name the two radial scales.

    ``mode="end"`` puts each title above the far end of its straight edge;
    ``mode="along"`` rotates it to lie along the edge, offset outwards, which
    is the only thing that works when the edge ends are the busiest corners of
    the figure.
    """
    if mode == "end":
        xl, yl = _wedge_xy(np.array([ra_hi]), np.array([dmax * (1 + pad)]),
                           ra_mid, stretch)
        ax.annotate("comoving distance\n(Mpc)", (xl[0], yl[0]),
                    fontsize=fontsize, ha="center", va="bottom",
                    linespacing=1.3, fontweight="bold")
        xr, yr = _wedge_xy(np.array([ra_lo]), np.array([dmax * (1 + pad)]),
                           ra_mid, stretch)
        ax.annotate("redshift", (xr[0], yr[0]), fontsize=fontsize,
                    ha="center", va="bottom", fontweight="bold")
        return

    half = np.radians(opening / 2.0) if opening is not None \
        else np.radians((ra_hi - ra_lo) * stretch / 2.0)
    for label, sign in (("comoving distance (Mpc)", +1), ("redshift", -1)):
        edge = np.pi / 2.0 + sign * half
        norm = edge + sign * np.pi / 2.0
        r = dmax * 0.58
        x = r * np.cos(edge) + dmax * pad * np.cos(norm)
        y = r * np.sin(edge) + dmax * pad * np.sin(norm)
        rot = np.degrees(edge)
        if rot > 135.0:
            rot -= 180.0
        ax.annotate(label, (x, y), fontsize=fontsize, ha="center", va="center",
                    rotation=rot, rotation_mode="anchor", fontweight="bold")


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_gama(verbose=True, nmin=3):
    """Every GAMA galaxy, and every GAMA group above the multiplicity floor."""
    gal = _fits_frame(GAMA_GALS)
    gal["GAMARegion"] = gal.GAMARegion.str.decode("utf-8")
    g = _fits_frame(GAMA_GROUPS)
    g["GAMARegion"] = g.GAMARegion.str.decode("utf-8")
    grp = pd.DataFrame({"ra": g.IterCenRA.values, "dec": g.IterCenDec.values,
                        "z": g.Zfof.values, "n": g.Nfof.values.astype(int),
                        "region": g.GAMARegion.values})
    npair = int((grp.n < nmin).sum())
    grp = grp[grp.n >= nmin]
    if verbose:
        print(f"  GAMA galaxies            : {len(gal)}")
        print(f"  GAMA groups (N >= {nmin})     : {len(grp)}"
              f"   ({npair} below the floor, dropped)")
    return gal, grp


def load_sdss(verbose=True, nmin=3):
    """Every Nessie SDSS galaxy, and every group above the multiplicity floor."""
    gal = pd.read_parquet(SDSS_GALS)
    g = pd.read_parquet(SDSS_GROUPS)
    grp = pd.DataFrame({"ra": g.iter_ra.values, "dec": g.iter_dec.values,
                        "z": g.median_redshift.values,
                        "n": g.multiplicity.values.astype(int)})
    npair = int((grp.n < nmin).sum())
    grp = grp[grp.n >= nmin]
    if verbose:
        print(f"  SDSS galaxies            : {len(gal)}")
        print(f"  SDSS groups (N >= {nmin})     : {len(grp)}"
              f"   ({npair} below the floor, dropped)")
    return gal, grp


# ---------------------------------------------------------------------------
# The GAMA cones
# ---------------------------------------------------------------------------

def cone_region(gal, g3c, region, outfile, zmax=0.25, opening=30.0,
                dticks=(200, 400, 600, 800, 1000),
                zticks=(0.05, 0.10, 0.15, 0.20, 0.25), group_scale=1.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gg = gal[(gal.GAMARegion == region) & (gal.Z < zmax) & (gal.Z > 0.0)]
    gs = g3c[(g3c.region == region) & (g3c.z > 0) & (g3c.z < zmax)]
    ra_lo, ra_hi = gg.RAcen.min(), gg.RAcen.max()
    ra_mid = 0.5 * (ra_lo + ra_hi)
    stretch = opening / (ra_hi - ra_lo)
    half = np.radians(opening / 2.0)
    dmax = COSMO.comoving_distance(zmax).value

    fig = plt.figure(figsize=(3.4, 6.8), dpi=600)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

    d_gal = COSMO.comoving_distance(gg.Z.values).value
    ax.scatter(*_wedge_xy(gg.RAcen.values, d_gal, ra_mid, stretch),
               s=0.40, c=GAL_C, lw=0, alpha=0.8, rasterized=True, zorder=2)

    d_grp = COSMO.comoving_distance(gs.z.values).value
    ax.scatter(*_wedge_xy(gs.ra.values, d_grp, ra_mid, stretch),
               s=mult_size(gs.n.values, scale=3.64 * group_scale),
               facecolors="none",
               edgecolors=GRP_C, lw=0.5, zorder=3)

    draw_frame(ax, ra_lo, ra_hi, ra_mid, stretch, dmax, dticks, zticks,
               fontsize=8.0, labelpad=9.5)
    ra_ticks(ax, ra_lo, ra_hi, ra_mid, stretch, dmax, target=3, fontsize=8.5)
    edge_titles(ax, ra_lo, ra_hi, ra_mid, stretch, dmax, fontsize=9.5,
                pad=0.128, mode="along", opening=opening)

    ax.set_aspect("equal")
    ax.axis("off")
    # Frame pinned to twice as tall as it is wide.
    xhalf, ylo = dmax * 0.345, -dmax * 0.20
    frame_box(ax, -xhalf, xhalf, ylo, ylo + 4.0 * xhalf)
    ax.set_xlim(-xhalf * 1.02, xhalf * 1.02)
    ax.set_ylim(ylo * 1.02, ylo + 4.0 * xhalf * 1.02)

    mult_legend(ax, loc="lower left", fontsize=8)
    ax.text(0.97, 0.015, region.upper(), transform=ax.transAxes, fontsize=11,
            ha="right", va="bottom", fontweight="bold")

    plotting.end_plot(outfile)
    plt.close(fig)
    print(f"  wrote {outfile}   ({len(gs)} groups, {len(gg)} galaxies)")


# ---------------------------------------------------------------------------
# The SDSS wedge
# ---------------------------------------------------------------------------

def sdss_puck(gal, grp, outfile, zmax=0.08, width=8.3,
              dticks=(100, 200, 300), zticks=(0.02, 0.04, 0.06, 0.08),
              group_scale=1.0):
    """The SDSS wedge, drawn in Cartesian coordinates.

    A polar axes keeps the bounding box of the whole circle even when
    ``set_thetamin`` trims it to a sector, so the sector ends up small in the
    middle of a mostly empty canvas.  Building the wedge by hand instead means
    the figure is exactly the size of the thing being drawn -- which matters
    here because this one runs across both columns.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = gal[(gal.zobs > 0) & (gal.zobs < zmax)]
    ra_lo, ra_hi = g.RAJ2000.min(), g.RAJ2000.max()
    ra_mid = 0.5 * (ra_lo + ra_hi)
    half = np.radians((ra_hi - ra_lo) / 2.0)
    dmax = COSMO.comoving_distance(zmax).value

    fig = plt.figure(figsize=(width, width / 1.94), dpi=600)
    ax = fig.add_axes([0.005, 0.02, 0.99, 0.96])

    d = COSMO.comoving_distance(g.zobs.values).value
    ax.scatter(*_wedge_xy(g.RAJ2000.values, d, ra_mid, 1.0),
               s=0.15, c=GAL_C, lw=0, alpha=0.65, rasterized=True, zorder=2)

    gp = grp[(grp.z > 0) & (grp.z < zmax)]
    dg = COSMO.comoving_distance(gp.z.values).value
    ax.scatter(*_wedge_xy(gp.ra.values, dg, ra_mid, 1.0),
               s=mult_size(gp.n.values, scale=0.45 * group_scale),
               facecolors="none",
               edgecolors=GRP_C, lw=0.35, zorder=3)

    draw_frame(ax, ra_lo, ra_hi, ra_mid, 1.0, dmax, dticks, zticks,
               fontsize=8.0, ticklen=0.022, labelpad=8.0)
    edge_titles(ax, ra_lo, ra_hi, ra_mid, 1.0, dmax, fontsize=9.5, pad=0.150,
                mode="along")

    ra_ticks(ax, ra_lo, ra_hi, ra_mid, 1.0, dmax, step=20.0, fontsize=8.5,
             labelpad=0.032)

    ax.set_aspect("equal")
    ax.axis("off")
    xm = dmax * np.sin(half) * 1.10
    ax.set_xlim(-xm, xm)
    ax.set_ylim(-dmax * 0.20, dmax * 1.16)

    mult_legend(ax, loc="upper left", fontsize=8.5)
    ax.text(0.985, 0.97, f"Nessie SDSS, $z < {zmax}$",
            transform=ax.transAxes, fontsize=11, ha="right", va="top",
            fontweight="bold")

    plotting.end_plot(outfile)
    plt.close(fig)
    print(f"  wrote {outfile}   ({len(gp)} groups, {len(g)} galaxies)")


# ---------------------------------------------------------------------------
# Both surveys, all the way round
# ---------------------------------------------------------------------------

def sky_polar(gama_gal, sdss_gal, outfile, zmax=0.25, width=7.0, zstep=0.05):
    """RA all the way round, declination flattened, both surveys on one scale.

    This is the figure that shows what the two catalogues actually are: SDSS a
    wide shallow cap, GAMA four narrow spokes reaching several times as deep.

    Rings are equally spaced in redshift, so the outermost ring *is* the edge of
    the plot and nothing is labelled twice.  Note this makes the rings unevenly
    spaced on the page, since comoving distance is not linear in z.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zrings = np.arange(zstep, zmax + 1e-9, zstep)
    dmax = COSMO.comoving_distance(zmax).value

    fig = plt.figure(figsize=(width, width), dpi=600)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

    def xy(ra, dist):
        t = np.radians(np.asarray(ra, float))
        return dist * np.cos(t), dist * np.sin(t)

    s = sdss_gal[(sdss_gal.zobs > 0) & (sdss_gal.zobs < zmax)]
    ax.scatter(*xy(s.RAJ2000.values,
                   COSMO.comoving_distance(s.zobs.values).value),
               s=0.08, c=SDSS_C, lw=0, alpha=0.42, rasterized=True, zorder=2)
    gm = gama_gal[(gama_gal.Z > 0) & (gama_gal.Z < zmax)]
    ax.scatter(*xy(gm.RAcen.values,
                   COSMO.comoving_distance(gm.Z.values).value),
               s=0.08, c=GAMA_C, lw=0, alpha=0.5, rasterized=True, zorder=3)

    ring = np.linspace(0, 360, 1200)
    for z in zrings:
        d = float(COSMO.comoving_distance(z).value)
        last = z >= zmax - 1e-9
        ax.plot(*xy(ring, np.full_like(ring, d)), color="k",
                lw=0.7 if last else 0.5, ls="-" if last else (0, (4, 3)),
                alpha=1.0 if last else 0.35, zorder=4)
        ax.annotate(f"{z:g}", (0.0, d), fontsize=8, ha="center", va="center",
                    zorder=7, fontweight="bold",
                    bbox=dict(fc="w", ec="none", alpha=0.9, pad=0.6))

    for ra_t in np.arange(0, 360, 30):
        x0, y0 = xy(np.array([ra_t]), np.array([dmax]))
        x1, y1 = xy(np.array([ra_t]), np.array([dmax * 0.982]))
        ax.plot([x0[0], x1[0]], [y0[0], y1[0]], color="k", lw=0.8, zorder=6)
        x2, y2 = xy(np.array([ra_t]), np.array([dmax * 1.055]))
        ax.annotate(f"{ra_t:.0f}$^\\circ$", (x2[0], y2[0]), fontsize=8.5,
                    ha="center", va="center", zorder=7, fontweight="bold")

    ax.set_aspect("equal")
    ax.axis("off")
    lim = dmax * 1.10
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    handles = [plt.Line2D([], [], ls="", marker="o", ms=4, mfc=c, mec="none",
                          label=l)
               for c, l in [(GAMA_C, f"GAMA ({len(gm):,} galaxies)"
                             .replace(",", "\u2009")),
                            (SDSS_C, f"SDSS ({len(s):,} galaxies)"
                             .replace(",", "\u2009"))]]
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, frameon=False,
              bbox_to_anchor=(0.0, 1.0))

    plotting.end_plot(outfile)
    plt.close(fig)
    print(f"  wrote {outfile}   ({len(gm)} GAMA + {len(s)} SDSS galaxies)")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cones", action="store_true", help="one per GAMA region")
    p.add_argument("--puck", action="store_true", help="the SDSS wedge")
    p.add_argument("--sky", action="store_true", help="both surveys, 360 deg")
    p.add_argument("--all", action="store_true")
    p.add_argument("--opening", type=float, default=30.0,
                   help="drawn wedge opening angle in degrees (default 30)")
    p.add_argument("--sdss-zmax", type=float, default=0.08)
    p.add_argument("--nmin", type=int, default=3,
                   help="multiplicity floor for the plotted groups (default 3; "
                        "N = 2 pairs are a separate catalogue product)")
    p.add_argument("--group-scale", type=float, default=1.0,
                   help="multiply every group marker area by this.  The full "
                        "catalogue is ~9x the HMF sample, so values below 1 "
                        "are what let the galaxies show through again.")
    args = p.parse_args()
    if args.all:
        args.cones = args.puck = args.sky = True
    if not (args.cones or args.puck or args.sky):
        p.error("pick at least one of --cones / --puck / --sky / --all")

    print("Survey figures")
    gama_gal = g3c = sdss_gal = sdss_grp = None
    if args.cones or args.sky:
        gama_gal, g3c = load_gama(nmin=args.nmin)
    if args.puck or args.sky:
        sdss_gal, sdss_grp = load_sdss(nmin=args.nmin)

    if args.cones:
        for region in ["g09", "g12", "g15", "g23"]:
            cone_region(gama_gal, g3c, region, f"cone_{region}.png",
                        opening=args.opening, group_scale=args.group_scale)
    if args.puck:
        sdss_puck(sdss_gal, sdss_grp, "sdss_puck.png", zmax=args.sdss_zmax,
                  group_scale=args.group_scale)
    if args.sky:
        sky_polar(gama_gal, sdss_gal, "sky_polar.png")


if __name__ == "__main__":
    main()
