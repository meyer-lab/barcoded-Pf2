from pathlib import Path

import anndata
import hdf5plugin  # noqa: F401
import numpy as np
import pandas as pd
import scanpy as sc
from glmpca.glmpca import glmpca
from scipy.sparse import csr_matrix
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


def import_CCLE(pca_option) -> anndata.AnnData:
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
        data = anndata.io.read_text(
            Path("./pf2barcode/data/" + name + "_count_mtx.tsv.bz2")
        ).T
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
    # conditional statement for either glm_pca or pca
    if pca_option == "glm_pca":
        # convert from sparse to dense matrix
        # matrix must be transposed for this implementation of glm_pca to give the same
        # shape output matrix as scanpy PCA
        X_dense = X.X.toarray().T
        # run glm_pca with the dense matrix and 30 components
        glmpca_result = glmpca(X_dense, L=30, verbose=True)
        # construct a new anndata object using the output dict from glm_pca
        X_new = anndata.AnnData(X=glmpca_result["factors"])
        # copy over the data.obs from the original data
        X_new.obs = X.obs.copy()
        X = X_new
        # manually add an x_pca column to the anndata object so that it works with the figure files
        X.obsm["X_pca"] = X.X
    else:
        X = prepare_dataset(X, geneThreshold=0.001)
        sc.pp.pca(X, n_comps=30, svd_solver="arpack")
    return X


def process_GSE150949(data_file):
    """Preprocessing for GSE150949. This is performed to generate the annData files
    that we will use more regularly."""
    # read in the meta data using pandas
    metadata = pd.read_csv(
        "/opt/extra-storage/GSE150949/GSE150949_metaData_with_lineage.txt.gz",
        delimiter="\t",
        engine="python",
    )
    # separate the columns to merge with the data
    columns = metadata[["full_cell_barcode", "lineage_barcode"]]
    # create anndata object of the data file
    data = anndata.io.read_csv(data_file, delimiter=",")
    # merge the data file object with the metadata coluns
    data.obs = data.obs.join(columns, how="left")
    data.obs["full_cell_barcode"] = data.obs["full_cell_barcode"].astype(str)
    data.obs["lineage_barcode"] = data.obs["lineage_barcode"].astype(str)
    data.X = csr_matrix(data.X)
    return data


def import_GSE150949():
    # read in the meta data using anndata
    data = anndata.io.read_h5ad(
        "/opt/extra-storage/GSE150949/GSE150949_pooled_watermelon.data.h5"
    )
    count = anndata.io.read_h5ad(
        "/opt/extra-storage/GSE150949/GSE150949_pooled_watermelon.count.h5"
    )

    return data, count
