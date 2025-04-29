from pathlib import Path

import anndata
import hdf5plugin  # noqa: F401
import numpy as np
import pandas as pd
import scanpy as sc
from anndata.io import read_text
from glmpca.glmpca import glmpca
from scipy.sparse import csr_array, csr_matrix
from sklearn.preprocessing import scale
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
    n_i = y_ij.sum(axis=1)

    # MLE of gene expression
    pi_j = y_ij.sum(axis=0) / np.sum(n_i)

    non_y_ij = n_i[:, None] - y_ij
    mu_ij = n_i[:, None] * pi_j[None, :]
    signs = np.sign(y_ij - pi_j[None, :])

    first_term = 2 * y_ij * np.log(np.maximum(y_ij, 1.0) / mu_ij)
    second_term = 2 * non_y_ij * np.log(non_y_ij / (n_i[:, None] - mu_ij))

    X.X = signs * np.sqrt(np.maximum(first_term + second_term, 0.0))

    X.X = scale(X.X)

    assert np.all(np.isfinite(X.X))  # type: ignore
    return X


def import_CCLE(pca_option="other") -> anndata.AnnData:
    # pca option should be passed as either pca or glm_pca
    """Imports barcoded cell data."""
    adatas = {}
    barcode_dfs = []

    for name in (
        # "HCT116_tracing_T1",
        # "T2_HCT116",
        "T1_MDAMB231",
        "T2_MDAMB231",
    ):
        data = read_text(Path("./pf2barcode/data/" + name + "_count_mtx.tsv.bz2")).T
        barcodes = pd.read_csv(
            "./pf2barcode/data/" + name + "_SW.txt", sep="\t", index_col=0, header=0
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

    # conditional statement for either glm_pca or pca
    if pca_option == "glm_pca":
        # convert from sparse to dense matrix
        # matrix must be transposed for this implementation of glm_pca to give the same
        # shape output matrix as scanpy PCA
        X_dense = X.X.toarray().T

        # run glm_pca with the dense matrix and 20 components
        glmpca_result = glmpca(
            X_dense, L=20, verbose=True, ctl={"maxIter": 2, "eps": 0.0001}
        )

        X.obsm["X_pca"] = glmpca_result["factors"]
        X.varm["PCs"] = glmpca_result["loadings"]
    if pca_option == "dev_pca":
        X = prepare_dataset_dev(X)
        sc.pp.pca(X, n_comps=20, svd_solver="arpack")
    else:
        X = prepare_dataset(X, geneThreshold=0.001)
        sc.pp.pca(X, n_comps=20, svd_solver="arpack")
    return X
