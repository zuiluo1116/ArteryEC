# chooseBioCmirror () 

# BiocManager::install(c("scry","irlba","dendextend","BiocNeighbors","matrixStats","RSpectra","sfsmisc","fastcluster","data.tree"))
# BiocManager::install(c("data.table","tidyverse","stringr","COSG","plyr","httr","jsonlite"))
# BiocManager::install(c("hdf5r"))
# BiocManager::install(c("stirngr"))

# devtools::install_github("hhoeflin/hdf5r")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("tidyverse")
# devtools::install_github("igrabski/sc-SHC")

library(scSHC)
library(Seurat)
library(SeuratDisk)
library(stringr)
setwd("./sc-SHC")


files <- list.files()
files <- files[str_detect(files,pattern = "h5seurat")]

# [1] "GSE155468_EC.h5seurat"             "GSE155514_EC.h5seurat"             "GSE159677_EC.h5seurat"            
# [4] "GSE213740_EC.h5seurat"             "GSE216860_EC.h5seurat"             "NOte"                             
# [7] "rstudio-server-1.4.1103-amd64.deb"


GSE155468_EC <- LoadH5Seurat(files[1]) #919
GSE155514_EC <- LoadH5Seurat(files[2]) #1223
GSE159677_EC <- LoadH5Seurat(files[3]) #5852
GSE213740_EC <- LoadH5Seurat(files[4]) #9403
GSE216860_EC <- LoadH5Seurat(files[5]) #8662


library(scSHC)
# sc <- GSE155468_EC
# sc <- GSE155514_EC
# sc <- GSE159677_EC
# sc <- GSE213740_EC
sc <- GSE216860_EC


data <- sc@assays[["RNA"]]@counts
batch_ident <- sc$orig.ident
clusters <- scSHC(data,batch = batch_ident,parallel = T)
# new_clusters <- testClusters(data, as.character(batch_ident))
sc$new_cluster <- clusters[[1]]

pdf(file = "860.dim.pdf")
DimPlot(sc,group.by = "new_cluster")
dev.off()

length(clusters[[1]]) #919
length(SHC_468[[1]]) #919
length(SHC_514[[1]]) #1223
length(SHC_677[[1]]) #5852
length(SHC_740[[1]]) #9403
length(SHC_860[[1]]) #8662


# SHC_468 <- clusters
# SHC_514 <- clusters
# SHC_677 <- clusters
# SHC_740 <- clusters
SHC_860 <- clusters


save(list = c("SHC_468","SHC_514","SHC_677","SHC_740","SHC_860"),file = "clusters.RData")

rm(list = ls())

gc()


sc_list <- list(GSE155468_EC,GSE155514_EC,GSE159677_EC,GSE213740_EC,GSE216860_EC)
gene_list <- lapply(sc_list, function(x){
  x <- x@assays[["RNA"]]@counts@Dimnames[[1]]
})

GSE155468_EC
gene_sel <- gene_list[[1]]

for (j in 2:length(gene_list)) {
  gene_sel <- intersect(gene_sel,gene_list[[j]])
}


sc_sel_list <- lapply(sc_list, function(x){
  x <- x[gene_sel,]
  tmp <- x
  DefaultAssay(tmp) <- "RNA"
  tmp@assays$SCT <- NULL
  x <- tmp
})
GSEs <- c("GSE155468","GSE155514","GSE159677","GSE213740","GSE216860")

for (i in 1:length(sc_sel_list)) {
  sc_sel_list[[i]]$GSE <- GSEs[i]
  sc_sel_list[[i]]$id <- paste0(sc_sel_list[[i]]$GSE,"_",sc_sel_list[[i]]$orig.ident)
}

sc_sel <- sc_sel_list[[1]]

for (j in 2:length(gene_list)) {
  sc_sel <- merge(sc_sel,sc_sel_list[[j]])
}
levels(factor(sc_sel$id))


data <- sc_sel@assays[["RNA"]]@counts
batch_ident <- sc_sel$id
clusters <- scSHC(data,batch = batch_ident,parallel = T)
# new_clusters <- testClusters(data, as.character(batch_ident))
sc$new_cluster <- clusters[[1]]

save(list = c("SHC_468","SHC_514","SHC_677","SHC_740","SHC_860","clusters"),file = "clusters.RData")





