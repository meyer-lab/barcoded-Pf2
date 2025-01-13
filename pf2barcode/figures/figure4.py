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
    
    # Get a list of the axis objects and create a figure.
    ax, f = getSetup((10, 6), (1, 1))
    subplotLabel(ax)

    #  Initialize an array of zeros to store p-values for each PC
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
        
    # Barplot setup
    sns.barplot(x=np.arange(pvalues.shape[0]), y=-np.log10(pvalues))
    plt.xlabel("PC")
    plt.ylabel("-log10(p-value)")

    return f
