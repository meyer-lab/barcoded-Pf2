"""
Figure 3
"""

import os

import pandas as pd
import scanpy as sc
from anndata import AnnData

from pf2barcode.imports import import_CCLE


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    X = import_CCLE("glm_pca")

    # Get rid of cells with no barcode
    X = X[X.obs["SW"] != "unknown"]

    # remove cells with barcodes having less than 10 cells
    good_SW = X.obs["SW"].value_counts().index[X.obs["SW"].value_counts() > 10]
    X = X[X.obs["SW"].isin(good_SW)]

    pcadata = AnnData(
        X.obsm["X_pca"],
        obs=X.obs,
        var=pd.DataFrame(
            index=[f"PC{i}" for i in range(1, X.obsm["X_pca"].shape[1] + 1)]
        ),
    )
    pcadata.X /= pcadata.X.mean(axis=1, keepdims=True)

    sc.settings.figdir = "./output/"
    sc.pl.heatmap(
        pcadata,
        pcadata.var_names,
        groupby="SW",
        dendrogram=True,
        vmax=20.0,
        vmin=-20.0,
        cmap="magma",
        show=False,
        save="figure3.svg",
    )
    os.rename("./output/heatmapfigure3.svg", "./output/figure3.svg")

    return None
