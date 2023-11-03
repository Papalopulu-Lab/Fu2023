library(Seurat)
library(dplyr)
library(ggplot2)

# If your expression matrix contains a mix of different pre-determined cell clusters or labels, it may be wise to run OscoNet seperately for each cell type.
# To do so specify the cell label you would like to isolate and run below. Remember to modify your output file names (counts, clusters) to reflect this.

cellType <- 'cell label' 

recluster<-function (ds){
  
  ## If you are running OscoNet on only a subpopulation then uncomment and use the code below, otherwise proceed as normal.
  
  # cols <- colnames(ds)
  # cellsList <- c()
  # for (I in 1:length(cols)) {if (startsWith(cols[I], 'MCF7_RM_CD44L') == T) {cellsList <- c(cellsList, cols[I])}}
  # df = subset(ds, select = cellsList)
  
  df = ds
  
  pbmc_MGH <- CreateSeuratObject(counts = df)
  pbmc_MGH <- ScaleData(object = pbmc_MGH)
  pbmc_MGH <- NormalizeData(object = pbmc_MGH)
  pbmc_MGH <- FindVariableFeatures(object = pbmc_MGH)
  pbmc_MGH <- RunPCA(object = pbmc_MGH)
  pbmc_MGH <- FindNeighbors(object = pbmc_MGH)
  pbmc_MGH <- FindClusters(object = pbmc_MGH)
  pbmc_MGH <- RunTSNE(object = pbmc_MGH)
  pbmc_MGH <- RunUMAP(pbmc_MGH, dims = 1:10)
  DimPlot(pbmc_MGH, reduction = "umap")
  DimPlot(pbmc_MGH, reduction = "tsne")
  
  counts <- pbmc_MGH@assays$RNA@counts
  clusters <- pbmc_MGH$seurat_clusters
  names(clusters) <- colnames(pbmc_MGH)
  out <- list("counts" = counts, "clusters" = clusters, "shape" = "round")
  return(out) 
}

# Select appropriate line to match file extension of your matrix 

ds=read.csv("/OscoNet/YourDataset/Data/expressionMatrix.csv", sep=',')
# ds=read.csv("/OscoNet/YourDataset/Data/expressionMatrix.txt", sep='\t')
# ds=read.csv("/OscoNet/YourDataset/Data/expressionMatrix.tsv", sep='\t')

out=recluster(ds)

# Modify these file names below if only running on a subset of this data.
write.csv(out$counts, file="/OscoNet/YourDataset/Data/counts.csv")
write.csv(out$clusters, file="/OscoNet/YourDataset/Data/clusters.csv")
rm(ds)








