"""
Test the cross validation accuracy.
"""

from ..imports import import_CCLE, import_GSE150949


def test_imports():
    """Test import functions."""
    X = import_CCLE()
    print(f"Data shape: {X.shape}")

    data, count = import_GSE150949()
