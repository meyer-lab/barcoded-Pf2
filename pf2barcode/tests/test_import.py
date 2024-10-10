"""
Test the cross validation accuracy.
"""

import numpy as np

from ..imports import (
    import_CCLE,
    DataToScores
)


def test_imports():
    """Test import functions."""
    X = import_CCLE()
    print(f"Data shape: {X.shape}")
    
    X = DataToScores('/Users/quinnmccall/Desktop/GSE150949_pooled_watermelon.data.matrix.csv.gz', '/Users/quinnmccall/Desktop/GSE150949_metaData_with_lineage.txt.gz')
    Y = DataToScores('/Users/quinnmccall/Desktop/GSE150949_pooled_watermelon.count.matrix.csv.gz', '/Users/quinnmccall/Desktop/GSE150949_metaData_with_lineage.txt.gz')
    assert isinstance(X, np.ndarray), "X is not a numpy array"
    assert isinstance(Y, np.ndarray), "Y is not a numpy array"
    assert X.size == Y.size

test_imports()
