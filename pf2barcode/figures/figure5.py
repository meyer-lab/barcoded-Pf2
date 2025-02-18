"""
Visualizes a scatter plot, showing the relationship between % variance and -log10(pvalue)
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
    X_pca = X.obsm["X_pca"]  # PCA-transformed data
   
    # Compute % variance explained
    total_variance = np.sum(np.var(X_pca, axis = 0))
    variance_explained = (np.var(X_pca, axis = 0) / total_variance) * 100  # Convert to %

    # Implement kruskal_pvalues function
    pvalues = kruskal_pvalues(X)

   
    # Convert p-values to -log10 scale
    neg_log_pvalues = -np.log10(pvalues)

    # Create scatter plot
    sns.scatterplot(x = variance_explained, y = neg_log_pvalues, color = "blue")

    # Label axes
    plt.xlabel("Variance Explained (%)")
    plt.ylabel("-log10(p-value)")
    plt.title("PCA Significance vs Variance Explained")

    for i in range(len(variance_explained)):
        plt.annotate(
            f"PC{i + 1}",  # Text to display
            xy = (variance_explained[i], neg_log_pvalues[i]),  # Point to annotate (as a tuple)
            xytext = (5, 5),  # Offset for the text
            textcoords = "offset points",  # Coordinate system for xytext
            ha = 'center',  # Horizontal alignment of the text
            fontsize = 6,  # Font size
)

    return f