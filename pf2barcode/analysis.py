"""
This file contains the kruskal_pvalues function that is used in figure 4 and 5.
"""

import numpy as np
from scipy.stats import kruskal


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

    return pvalues
