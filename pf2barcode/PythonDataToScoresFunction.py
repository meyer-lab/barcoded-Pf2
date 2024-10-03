import pandas as pd
import anndata
import scanpy as sc

def DataToScores(data_file, metadata_file):
    #read in the meta data using pandas
    metadata = pd.read_csv(metadata_file, delimiter = '\t', engine = 'python')
    #separate the columns to merge with the data
    columns = metadata[['full_cell_barcode', 'lineage_barcode']]
    #create anndata object of the data file
    data = anndata.read_csv(data_file, delimiter=',')
    #merge the data file object with the metadata coluns
    data.obs = data.obs.join(columns, how ='left')
    #run PCA
    sc.tl.pca(data)
    #print the scores plot 
    print(data.obsm['X_pca'])
    print('PCA computation done and scores added to AnnData object.')
#run the function
DataToScores('/Users/quinnmccall/Desktop/GSE150949_pooled_watermelon.data.matrix.csv.gz', '/Users/quinnmccall/Desktop/GSE150949_metaData_with_lineage.txt.gz')
DataToScores('/Users/quinnmccall/Desktop/GSE150949_pooled_watermelon.count.matrix.csv.gz', '/Users/quinnmccall/Desktop/GSE150949_metaData_with_lineage.txt.gz')