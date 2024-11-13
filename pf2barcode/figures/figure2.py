"""
Figure 2
"""

import pandas as pd
import scanpy as sc
import seaborn as sns

from pf2barcode.imports import import_CCLE

from .common import (
    getSetup,
    subplotLabel,
)

def makeFigure():

    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((6, 6), (1, 1))
    subplotLabel(ax)

    X = import_CCLE()

    Xsel = X[:, X.var.index.str.contains('AXL')]
    Xsel.uns

    sc.pl.pca(X, color="SW", components=['1, 2'], ax=ax[0], size=10)

    return f
