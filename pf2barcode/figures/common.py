"""
This file contains functions that are used in multiple figures.
"""

import sys
import time
from string import ascii_letters

import matplotlib
import numpy as np
import seaborn as sns
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import kruskal

matplotlib.use("AGG")

matplotlib.rcParams["legend.labelspacing"] = 0.2
matplotlib.rcParams["legend.fontsize"] = 8
matplotlib.rcParams["xtick.major.pad"] = 1.0
matplotlib.rcParams["ytick.major.pad"] = 1.0
matplotlib.rcParams["xtick.minor.pad"] = 0.9
matplotlib.rcParams["ytick.minor.pad"] = 0.9
matplotlib.rcParams["legend.handletextpad"] = 0.5
matplotlib.rcParams["legend.handlelength"] = 0.5
matplotlib.rcParams["legend.framealpha"] = 0.5
matplotlib.rcParams["legend.markerscale"] = 0.7
matplotlib.rcParams["legend.borderpad"] = 0.35
matplotlib.rcParams["svg.fonttype"] = "none"

def getSetup(
    figsize: tuple[float, float], gridd: tuple[int, int]
) -> tuple[list[Axes], Figure]:
    """Establish figure set-up with subplots."""
    sns.set_theme(
        style="whitegrid",
        font_scale=0.7,
        color_codes=True,
        palette="colorblind",
        rc={"grid.linestyle": "dotted", "axes.linewidth": 0.6},
    )

    # Setup plotting space and grid
    f = plt.figure(figsize=figsize, layout="constrained")
    gs1 = gridspec.GridSpec(gridd[0], gridd[1], figure=f)

    # Get list of axis objects
    ax = [f.add_subplot(gs1[x]) for x in range(gridd[0] * gridd[1])]

    return ax, f


def subplotLabel(axs: list[Axes]):
    """Place subplot labels on figure."""
    for ii, ax in enumerate(axs):
        ax.text(
            -0.2,
            1.2,
            ascii_letters[ii],
            transform=ax.transAxes,
            fontweight="bold",
            va="top",
        )


def genFigure():
    """Main figure generation function."""
    start = time.time()
    nameOut = "figure" + sys.argv[1]

    exec(f"from pf2barcode.figures.{nameOut} import makeFigure", globals())
    ff = makeFigure()  # type: ignore # noqa: F821

    if ff is not None:
        ff.savefig(
            f"./output/{nameOut}.svg", dpi=300, bbox_inches="tight", pad_inches=0
        )

    print(f"Figure {sys.argv[1]} is done after {time.time() - start} seconds.\n")

def kruskal_pvalues(X):

    # Initialize an array of zeros to store p-values for each PC
    pvalues = np.zeros(X.obsm["X_pca"].shape[1])

    # Loop through each principal component (PC) to compute p-values
    for jj in range(X.obsm["X_pca"].shape[1]):
        cells = [] 
        # Loop through each unique barcode in the "SW" observation column
        for barcodes in X.obs["SW"].unique():
            # Select data corresponding to current iteration of barcode and PC
            cells_selected = X.obsm["X_pca"][X.obs["SW"] == barcodes, jj]
            # Flatten data and add to list
            cells.append(cells_selected.flatten())
        # Perform the Kruskal-Wallis H-test
        pvalues[jj] = kruskal(*cells).pvalue

    return(pvalues)

