"""One publication style for every figure in this project.

Call ``plotstyle.apply()`` before creating a figure.  Everything that draws
should go through this rather than setting rcParams locally, so the figure set
stays internally consistent.

House rules, applied here:

* larger fonts throughout (labels 15, ticks 13) -- the previous defaults were
  legible on screen and unreadable at journal column width;
* tick marks pointing **inwards**, on all four sides;
* **minor ticks on**, also inwards;
* PNG output at 200 dpi rather than PDF.

PNG is deliberate: these figures carry thousands of Monte-Carlo spaghetti lines
and filled contours, which make vector PDFs enormous and slow to open, and every
consumer of them (referee reports, slides, the web) rasterises anyway.
"""

import matplotlib

DPI = 200
EXT = ".png"


def apply(scale=1.0):
    """Set the house style.  ``scale`` multiplies every font size."""
    matplotlib.rcParams.update({
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "font.size": 13 * scale,
        "axes.titlesize": 15 * scale,
        "axes.labelsize": 15 * scale,
        "xtick.labelsize": 13 * scale,
        "ytick.labelsize": 13 * scale,
        "legend.fontsize": 12 * scale,
        # inward ticks on all four sides, with minor ticks
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 7 * scale,
        "ytick.major.size": 7 * scale,
        "xtick.minor.size": 3.5 * scale,
        "ytick.minor.size": 3.5 * scale,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.minor.width": 0.9,
        "ytick.minor.width": 0.9,
        "axes.linewidth": 1.1,
        "legend.frameon": False,
    })


def as_png(path):
    """Rewrite an output path to PNG.  Used so existing --out defaults and any
    user-supplied .pdf name land on the house format without special-casing."""
    if path is None:
        return None
    for e in (".pdf", ".eps", ".svg", ".ps"):
        if path.lower().endswith(e):
            return path[: -len(e)] + EXT
    return path if path.lower().endswith(EXT) else path + EXT
