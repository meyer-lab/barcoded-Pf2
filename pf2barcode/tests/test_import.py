"""
Test the cross validation accuracy.
"""

from ..imports import import_CCLE


def test_imports():
    """Test import functions."""
    X = import_CCLE()
    print(f"Data shape: {X.shape}")
