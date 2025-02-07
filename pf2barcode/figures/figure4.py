import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import kruskal

from pf2barcode.imports import import_CCLE

from .common import (
    getSetup,
    subplotLabel,
    kruskal_pvalues, 
)

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


