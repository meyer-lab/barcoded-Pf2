"""
Test the cross validation accuracy.
"""

import numpy as np

from ..imports import import_CCLE, import_GSE150949


def test_imports():
    """Test import functions."""
    X = import_CCLE()
    print(f"Data shape: {X.shape}")

    X = import_GSE150949(
        "/opt/extra-storage/GSE150949/GSE150949_pooled_watermelon.data.matrix.csv.gz"
    )
    Y = import_GSE150949(
        "/opt/extra-storage/GSE150949/GSE150949_pooled_watermelon.count.matrix.csv.gz"
    )
    assert isinstance(X, np.ndarray), "X is not a numpy array"
    assert isinstance(Y, np.ndarray), "Y is not a numpy array"
    assert X.size == Y.size
