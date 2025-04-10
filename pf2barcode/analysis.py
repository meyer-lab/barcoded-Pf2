"""
This file contains the kruskal_pvalues function that is used in figure 4 and 5.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
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


def anova_pvalues(X):
    # get the number of PCs used
    n_pcs = X.obsm["X_pca"].shape[1]
    # initialize an array of empty P-values for each PC
    pvalues = np.zeros(n_pcs)

    # loop over every PC
    for jj in range(n_pcs):
        # construct a pandas dataframe with the PC and the associated group label
        df = pd.DataFrame({"PC": X.obsm["X_pca"][:, jj], "Group": X.obs["SW"].values})

        # apply ordinary least squares to model the PC values as a function of group membership
        model = smf.ols("PC ~ Group", data=df).fit()
        # apply anova to the model
        anova_table = sm.stats.anova_lm(model, typ=2)
        # assign the p values from the anova table
        pvalues[jj] = anova_table["PR(>F)"]["Group"]

    return pvalues
