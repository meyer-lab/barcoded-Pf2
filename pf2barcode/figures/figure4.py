"""
Generates a bar plot visualizing the relationship of PCs and computed
negative log10 p-values from the Kruskal-Wallis H-test

Computed p-values determines if distributions of PCs are statistically
significantly across different groups, and the negative log10 transformation
of the p-values allows for easier identification and interpretation of signficant PCs
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from pf2barcode.analysis import anova_pvalues, kruskal_pvalues
from pf2barcode.figures.common import (
    getSetup,
    subplotLabel,
)
from pf2barcode.imports import import_CCLE


def makeFigure():
    X = import_CCLE()

    # Get a list of the axis objects and create a figure.
    ax, f = getSetup((10, 6), (1, 1))
    subplotLabel(ax)

    # Implement kruskal_pvalues function
    pvalues = kruskal_pvalues(X)

    # Barplot setup
    sns.barplot(x=np.arange(pvalues.shape[0]), y=-np.log10(pvalues))
    plt.xlabel("PC")
    plt.ylabel("-log10(p-value)")

    return f


def makeFigureAnova():
    X = import_CCLE()

    # Get a list of the axis objects and create a figure.
    ax, f = getSetup((10, 6), (1, 1))
    subplotLabel(ax)

    # Implement anova_pvalues function
    pvalues = anova_pvalues(X)

    # Barplot setup
    sns.barplot(x=np.arange(pvalues.shape[0]), y=-np.log10(pvalues))
    plt.xlabel("PC")
    plt.ylabel("-log10(p-value)")

    return f
