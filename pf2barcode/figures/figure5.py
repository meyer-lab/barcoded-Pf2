"""
Visualizes a scatter plot, showing the relationship between % variance explained by each PC and the statistical significance of the PC (-log10(p-value))

Specifically, % of variance explained measures how much of the variability each PC captures

Helps identify which PCs are both highly explanatatory (captures a large portion of the variance) and statistically significiant, with each point containing annotations of its corresponding PC number
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal

from pf2barcode.imports import import_CCLE

from .common import (
    getSetup,
    subplotLabel,
)

from .analysis import kruskal_pvalues


def makeFigure():
    X = import_CCLE()

    # Get a list of the axis objects and create a figure
    ax, f = getSetup((10, 6), (1, 1))
    subplotLabel(ax)

    # Extract PCA results
    X_pca = X.obsm["X_pca"]  
   
    # Compute % variance explained
    total_variance = np.sum(np.var(X_pca, axis = 0))
    variance_explained = (np.var(X_pca, axis = 0) / total_variance) * 100  

    # Implement kruskal_pvalues function
    pvalues = kruskal_pvalues(X)

    neg_log_pvalues = -np.log10(pvalues)

    # Create scatter plot
    sns.scatterplot(x=variance_explained, y=neg_log_pvalues, color="blue")

    plt.xlabel("Variance Explained (%)")
    plt.ylabel("-log10(p-value)")
    plt.title("PCA Significance vs Variance Explained")

    for i in range(len(variance_explained)):
        plt.annotate(
            f"PC{i + 1}",  
            xy = (variance_explained[i], neg_log_pvalues[i]), 
            xytext = (5, 5), 
            textcoords = "offset points", 
            ha = 'center',  
            fontsize = 6,  
)

    return f
