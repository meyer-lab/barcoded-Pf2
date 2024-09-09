import pandas as pd
import anndata
import scanpy as sc

#read in the metadata
metadataLineage = pd.read_csv('/Users/quinnmccall/Desktop/GSE150949_metaData_with_lineage.txt.gz', compression = 'gzip', delimiter = '\t', engine = 'python')
#select the columns to merge
columns = metadataLineage[['full_cell_barcode', 'lineage_barcode']]
#read in the watermelon data matrix
watermelonData = anndata.read_csv('/Users/quinnmccall/Desktop/GSE150949_pooled_watermelon.data.matrix.csv.gz', delimiter=',')
#merge the watermelon data matrix with the selected columns from the metadata
watermelonData.obs = watermelonData.obs.join(columns, how ='left')
#perform PCA
sc.tl.pca(watermelonData)
#print the scores matrix
print(watermelonData.obsm['X_pca'])

print('PCA computation done and scores added to AnnData object.')