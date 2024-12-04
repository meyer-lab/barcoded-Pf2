import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import kruskal

from pf2barcode.imports import import_CCLE

from .common import (
    getSetup,
    subplotLabel,
)                                        

def makeFigure():

    X = import_CCLE()

    ax, f = getSetup((10, 6), (1, 1))
    subplotLabel(ax)

    pvalues = np.zeros(X.obsm["X_pca"].shape[1])

    for jj in range(X.obsm["X_pca"].shape[1]):
        cells = []
        for barcodes in X.obs["SW"].unique():
            cells_selected = X.obsm["X_pca"][X.obs["SW"] == barcodes, jj]
            cells.append(cells_selected.flatten())

        pvalues[jj] = kruskal(*cells).pvalue
        

    sns.barplot(x=np.arange(pvalues.shape[0]), y=-np.log10(pvalues))
    plt.xlabel("PC")
    plt.ylabel("-log10(p-value)")
    
    return f