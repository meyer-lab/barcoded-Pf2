import scanpy as sc
from pf2barcode.imports import import_GSE150949
import seaborn as sns

def makeFigures():
    X, Y = import_GSE150949()
    sc.pp.pca(X, n_comps=2, svd_solver="arpack")
    sc.pp.pca(Y, n_comps=2, svd_solver="arpack")
    print(X.var_names)
    print(X.obs.columns)
    sc.pl.pca(X, color="lineage_barcode", components="1,2", size=10, show=False, save="pca_X2.png")
    sc.pl.pca(Y, color="lineage_barcode", components="1,2", size=10, show=False, save="pca_Y2.png")

makeFigures()
