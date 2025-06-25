from pathlib import Path

import anndata
import hdf5plugin  # noqa: F401
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_array, csr_matrix
from scipy.special import xlogy
from sklearn.utils.sparsefuncs import (
    inplace_column_scale,
    mean_variance_axis,
)


def prepare_dataset(X: anndata.AnnData, geneThreshold: float) -> anndata.AnnData:
    assert isinstance(X.X, csr_matrix)
    assert np.amin(X.X.data) >= 0.0

    # Filter out genes with too few reads
    readmean, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    X = X[:, readmean > geneThreshold]

    # Normalize read depth
    sc.pp.normalize_total(
        X, exclude_highly_expressed=False, inplace=True, key_added="n_counts"
    )

    # Scale genes by sum
    readmean, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    readsum = X.shape[0] * readmean
    inplace_column_scale(X.X, 1.0 / readsum)  # type: ignore

    # Transform values
    X.X.data = np.log10((1000.0 * X.X.data) + 1.0)  # type: ignore

    return X


def prepare_dataset_dev(X: anndata.AnnData) -> anndata.AnnData:
    X.X = csr_array(X.X)  # type: ignore
    assert np.amin(X.X.data) >= 0.0

    # Remove cells and genes with fewer than 30 reads
    X = X[X.X.sum(axis=1) > 5, X.X.sum(axis=0) > 5]

    # Copy so that the subsetting is preserved
    X._init_as_actual(X.copy())

    # deviance transform
    y_ij = X.X.toarray()  # type: ignore

    # counts per cell
    n_i_col = y_ij.sum(axis=1).reshape(-1, 1)

    # MLE of gene expression
    pi_j = y_ij.sum(axis=0) / np.sum(n_i_col)
    mu_ij = n_i_col * pi_j[None, :]

    # --- Calculate Deviance Terms using numerically stable xlogy ---
    # D = 2 * [ y*log(y/mu) + (n-y)*log((n-y)/(n-mu)) ]
    # D = 2 * [ (xlogy(y, y) - xlogy(y, mu)) + (xlogy(n-y, n-y) - xlogy(n-y, n-mu)) ]

    n_minus_y = n_i_col - y_ij
    n_minus_mu = n_i_col - mu_ij

    # Term 1: y * log(y / mu) = xlogy(y, y) - xlogy(y, mu)
    # xlogy handles y=0 case correctly returning 0.
    term1 = xlogy(y_ij, y_ij) - xlogy(y_ij, mu_ij)

    # Term 2: (n-y) * log((n-y) / (n-mu)) = xlogy(n-y, n-y) - xlogy(n-y, n-mu)
    # xlogy handles n-y=0 case correctly returning 0.
    term2 = xlogy(n_minus_y, n_minus_y) - xlogy(n_minus_y, n_minus_mu)

    # Calculate full deviance: D = 2 * (term1 + term2)
    # Handle potential floating point inaccuracies leading to small negatives
    deviance = 2 * (term1 + term2)
    deviance = np.maximum(deviance, 0.0)  # Ensure non-negative before sqrt

    # Calculate signed square root residuals: sign(y - mu) * sqrt(D)
    signs = np.sign(y_ij - mu_ij)

    # Reset sign to 0 if deviance is exactly 0 (e.g. y=0, mu=0 or y=n, mu=n)
    # Avoids sign(-0.0) sometimes being -1
    signs[deviance == 0] = 0

    X.X = signs * np.sqrt(deviance)
    return X


def import_CCLE(pca_option="dev_pca") -> anndata.AnnData:
    # pca option should be passed as either pca or glm_pca
    """Imports barcoded cell data."""
    adatas = {}
    barcode_dfs = []

    # Get the directory containing this file
    current_dir = Path(__file__).parent

    for name in (
        # "HCT116_tracing_T1",
        # "T2_HCT116",
        # "T1_MDAMB231",
        "T2_MDAMB231",
    ):
        data = anndata.read_text(current_dir / "data" / f"{name}_count_mtx.tsv.bz2").T
        barcodes = pd.read_csv(
            current_dir / "data" / f"{name}_SW.txt", sep="\t", index_col=0, header=0
        )
        barcodes["sample"] = name
        barcodes = barcodes[~barcodes.index.duplicated(keep="first")]

        data.obs = data.obs.join(barcodes, how="left")

        barcode_dfs.append(barcodes)
        adatas[name] = data

    X = anndata.concat(adatas, label="sample", index_unique="-")
    X.X = csr_matrix(X.X)

    counts = X.obs["SW"].value_counts()
    counts = counts[counts > 5]

    X = X[X.obs["SW"].isin(counts.index), :]

    # Copy so that the subsetting is preserved
    X._init_as_actual(X.copy())

    # Counts per cell
    X.obs["n_counts"] = X.X.sum(axis=1)

    # conditional statement for either dev_pca or pca
    if pca_option == "dev_pca":
        X = prepare_dataset_dev(X)
        sc.pp.pca(X, n_comps=20, svd_solver="arpack")
    else:
        X = prepare_dataset(X, geneThreshold=0.001)
        sc.pp.pca(X, n_comps=20, svd_solver="arpack")

    return X
