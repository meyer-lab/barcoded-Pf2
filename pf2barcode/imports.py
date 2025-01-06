from pathlib import Path

import anndata
import hdf5plugin  # noqa: F401
import numpy as np
import pandas as pd
import scanpy as sc
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


def import_CCLE() -> anndata.AnnData:
    """Imports barcoded cell data."""
    adatas = {}
    barcode_dfs = []

    for name in (
        # "HCT116_tracing_T1",
        # "T2_HCT116",
        "T1_MDAMB231",
        "T2_MDAMB231",
    ):
        data = anndata.read_text(
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

    X = prepare_dataset(X, geneThreshold=0.001)

    sc.pp.pca(X, n_comps=30, svd_solver="arpack")
    return X



def process_GSE150949(data_file):
    """Preprocessing for GSE150949. This is performed to generate the annData files
    that we will use more regularly."""
    # read in the meta data using pandas
    metadata = pd.read_csv(
        "pf2barcode/data/GSE150949_pooled_watermelon.metadata.matrix.csv.gz",
        delimiter=",",
        engine="python",
    )
    # Checked for duplicates in the cell column to determine if cell is lineage barcodes or full cell barcodes 
    # result: No duplicates found in the 'cell' column.
    # no duplicates indicates that the cell column is the full cell barcode
    
    # set metadata index to the cell barcodes 
    metadata.set_index('cell', inplace=True)
  
    # create anndata object of the data file
    data = anndata.read_csv(data_file, delimiter=",", dtype=str)
    data = data.T
    data.obs = data.obs.join(metadata, how="left")
    # the cell column is now the data.obs.index (can be added manually if needed)
    print(data.obs.index)
    print(data.obs["orig.ident"])
    return data


def import_GSE150949():
    # read in the meta data using anndata
    data = anndata.read_h5ad(
        "/opt/extra-storage/GSE150949/GSE150949_pooled_watermelon.data.h5"
    )
    count = anndata.read_h5ad(
        "/opt/extra-storage/GSE150949/GSE150949_pooled_watermelon.count.h5"
    )

    return data, count
