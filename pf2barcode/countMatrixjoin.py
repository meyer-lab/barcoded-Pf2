import pandas as pd
import anndata
import scanpy as sc

#read in the metadata
metadataLineage = pd.read_csv('/Users/quinnmccall/Desktop/GSE150949_metaData_with_lineage.txt.gz', compression = 'gzip', delimiter = '\t', engine = 'python')
#select the columns to merge
columns = metadataLineage[['full_cell_barcode', 'lineage_barcode']]
#read in the watermelon count matrix
countMatrix = anndata.read_csv('/Users/quinnmccall/Desktop/GSE150949_pooled_watermelon.count.matrix.csv.gz', delimiter=',')
#merge the count matrix with the selected metadata columns
countMatrix.obs = countMatrix.obs.join(columns, how ='left')
#perform PCA
sc.tl.pca(countMatrix)
#print the scores matrix
print(countMatrix.obsm['X_pca'])

print('PCA computation done and scores added to AnnData object.')