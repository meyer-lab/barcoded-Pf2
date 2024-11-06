import seaborn as sns
import numpy as np
from scipy.stats import kruskal
import matplotlib.pyplot as plt

from pf2barcode.imports import import_CCLE

from .common import (
    getSetup,
    subplotLabel,
)                                        

def makeFigure():

    X = import_CCLE()

    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((10, 6), (1, 1))
    subplotLabel(ax)

    idx_selected = 0
    cells = []
    pvalues = np.zeros(X.obsm["X_pca"].shape[1])

    for jj in range(X.obsm["X_pca"].shape[1]):
        for barcodes in X.obs["SW"].unique():
            cells_selected = X.obsm["X_pca"][X.obs["SW"] == barcodes, idx_selected]
            cells.append(cells_selected.flatten())

        pvalues[jj] = kruskal(*cells).pvalue

    sns.barplot(x=np.arange(pvalues.shape[0]), y=-np.log10(pvalues))
    plt.xlabel("PC")
    plt.ylabel("-log10(p-value)")
    
    return f