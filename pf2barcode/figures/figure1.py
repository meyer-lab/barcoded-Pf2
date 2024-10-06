"""
Figure 3
"""

import pandas as pd
import seaborn as sns

from pf2barcode.imports import import_CCLE

from .common import (
    getSetup,
    subplotLabel,
)


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((6, 6), (2, 2))
    subplotLabel(ax)

    X = import_CCLE()

    # Get rid of cells with no barcode
    X = X[X.obs["SW"] != "unknown"]  # type: ignore

    # remove cells with barcodes having less than 10 cells
    good_SW = X.obs["SW"].value_counts().index[X.obs["SW"].value_counts() > 10]
    X = X[X.obs["SW"].isin(good_SW)]  # type: ignore

    df = pd.DataFrame({"reads": X.obs["n_counts"], "SW": X.obs["SW"]})

    sns.boxplot(data=df, x="reads", y="SW", ax=ax[0])

    return f
