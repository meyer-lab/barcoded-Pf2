"""
Creates a plot of the gene associations with the principal components
"""

import os
import scanpy as sc
from pf2barcode.imports import import_CCLE

def makeFigure():
    X = import_CCLE("pca")

    # Get rid of cells with no barcode
    X = X[X.obs["SW"] != "unknown"]

    # remove cells with barcodes having less than 10 cells
    good_SW = X.obs["SW"].value_counts().index[X.obs["SW"].value_counts() > 10]
    X = X[X.obs["SW"].isin(good_SW)]

    # Output directory
    sc.settings.figdir = "./output/"

    # Plot PCA loadings from X.varm["PCs"]
    sc.pl.pca_loadings(
        X,
        components=[1, 2, 3, 4, 5],  
        show=False,
        save="figure6.svg"
    )

    os.rename("./output/pca_loadingsfigure6.svg", "./output/figure6.svg")

    return None