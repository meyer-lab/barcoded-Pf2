#%%
import scanpy as sc
from pf2barcode.imports import import_CCLE
import seaborn as sns

X = import_CCLE()

Xsel = X[:, X.var.index.str.contains('AXL')]
Xsel.uns

sc.pl.pca(X, color="SW", components=['1, 2'], size=10)

# %%
