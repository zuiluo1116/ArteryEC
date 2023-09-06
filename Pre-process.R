library(Seurat)
library(harmony)
library(SeuratDisk)
library(data.table)
library(tidyverse)
library(stringr)
library(COSG)
library(plyr)
library(httr)
library(jsonlite)
setwd("E:/works/INSC-revision")

## function 
Topp.chr_list <- function(a)
{
  tmp <- a
  url = "https://toppgene.cchmc.org/API/lookup"
  encode = "json"
  payload = list(Symbols = tmp)
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  data = fromJSON(rawToChar(response$content))
  results = fromJSON(json)
  IDs = results[["Genes"]][["Entrez"]]
  url = "https://toppgene.cchmc.org/API/enrich"
  encode = "json"
  payload = list(Genes = IDs)
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  results = fromJSON(json)
  res <- results[["Annotations"]]
  res <- res[which(res$Category %in% "ToppCell"),]
  return(res)
}

## setseed
set.seed(1116)

## Pre-process
# Create Seurat objects from the original data shared in the GEO database
GSEs <- list.files()
GSEs <- GSEs[stringr::str_detect(GSEs,pattern = "GSE")]
# [1] "GSE155468" "GSE155514" "GSE159677" "GSE213740" "GSE216860"

dir(recursive = T)
# [1] "GSE155468/Con4/Con4.txt"            "GSE155468/Con6/Con6.txt"           
# [3] "GSE155468/Con9/Con9.txt"            "GSE155468/TAA1/TAA1.txt"           
# [5] "GSE155468/TAA2/TAA2.txt"            "GSE155468/TAA3/TAA3.txt"           
# [7] "GSE155468/TAA4/TAA4.txt"            "GSE155468/TAA5/TAA5.txt"           
# [9] "GSE155468/TAA6/TAA6.txt"            "GSE155468/TAA7/TAA7.txt"           
# [11] "GSE155468/TAA8/TAA8.txt"            "GSE155514/RPE004/RPE004.txt"       
# [13] "GSE155514/RPE005/RPE005.txt"        "GSE155514/RPE006/RPE006.txt"       
# [15] "GSE159677/P1_AC.h5"                 "GSE159677/P1_PA.h5"                
# [17] "GSE159677/P2_AC.h5"                 "GSE159677/P2_PA.h5"                
# [19] "GSE159677/P3_AC.h5"                 "GSE159677/P3_PA.h5"                
# [21] "GSE213740/AD_1/barcodes.tsv.gz"     "GSE213740/AD_1/features.tsv.gz"    
# [23] "GSE213740/AD_1/matrix.mtx.gz"       "GSE213740/AD_2/barcodes.tsv.gz"    
# [25] "GSE213740/AD_2/features.tsv.gz"     "GSE213740/AD_2/matrix.mtx.gz"      
# [27] "GSE213740/AD_3/barcodes.tsv.gz"     "GSE213740/AD_3/features.tsv.gz"    
# [29] "GSE213740/AD_3/matrix.mtx.gz"       "GSE213740/AD_4/barcodes.tsv.gz"    
# [31] "GSE213740/AD_4/features.tsv.gz"     "GSE213740/AD_4/matrix.mtx.gz"      
# [33] "GSE213740/AD_5/barcodes.tsv.gz"     "GSE213740/AD_5/features.tsv.gz"    
# [35] "GSE213740/AD_5/matrix.mtx.gz"       "GSE213740/AD_6/barcodes.tsv.gz"    
# [37] "GSE213740/AD_6/features.tsv.gz"     "GSE213740/AD_6/matrix.mtx.gz"      
# [39] "GSE213740/Normal_1/barcodes.tsv.gz" "GSE213740/Normal_1/features.tsv.gz"
# [41] "GSE213740/Normal_1/matrix.mtx.gz"   "GSE213740/Normal_2/barcodes.tsv.gz"
# [43] "GSE213740/Normal_2/features.tsv.gz" "GSE213740/Normal_2/matrix.mtx.gz"  
# [45] "GSE213740/Normal_3/barcodes.tsv.gz" "GSE213740/Normal_3/features.tsv.gz"
# [47] "GSE213740/Normal_3/matrix.mtx.gz"   "GSE216860/Normal_1/barcodes.tsv.gz"
# [49] "GSE216860/Normal_1/features.tsv.gz" "GSE216860/Normal_1/matrix.mtx.gz"  
# [51] "GSE216860/Normal_2/barcodes.tsv.gz" "GSE216860/Normal_2/features.tsv.gz"
# [53] "GSE216860/Normal_2/matrix.mtx.gz"   "GSE216860/Normal_3/barcodes.tsv.gz"
# [55] "GSE216860/Normal_3/features.tsv.gz" "GSE216860/Normal_3/matrix.mtx.gz"  
# [57] "GSE216860/Normal_4/barcodes.tsv.gz" "GSE216860/Normal_4/features.tsv.gz"
# [59] "GSE216860/Normal_4/matrix.mtx.gz"   "GSE216860/Normal_5/barcodes.tsv.gz"
# [61] "GSE216860/Normal_5/features.tsv.gz" "GSE216860/Normal_5/matrix.mtx.gz"  
# [63] "GSE216860/Normal_6/barcodes.tsv.gz" "GSE216860/Normal_6/features.tsv.gz"
# [65] "GSE216860/Normal_6/matrix.mtx.gz" 
mdir <- getwd()

# GSE155468
rm(list = setdiff(ls(),c("GSEs","Topp.chr_list","mdir")))

i <- 1
setwd(paste0(mdir,"/",GSEs[i]))
files <- dir(recursive = T)
read <- list()

for (j in 1:length(files)) {
  tmp <- fread(files[j])
  rnms <- tmp[,1]
  tmp <- as.matrix(tmp[,-1])
  rownames(tmp) <- rnms$V1
  read[[j]] <- CreateSeuratObject(counts = tmp,
                                  project = unlist(str_split(files[j],pattern = "/",n = 2))[1])
}
read <- lapply(read, function(x){
  tmp <- x
  tmp[["pt.mt"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^MT-")
  tmp[["pt.rb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^RP[SL]")
  tmp[["pt.hb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^HB[ABDEGMZ]")
  x <- tmp
})

reads <- read[[1]]
for (k in 2:length(read)) {
  reads <- merge(reads,read[[k]])
}

reads <- RenameCells(reads,add.cell.id = paste0(GSEs[i],"_"))
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")
reads <- subset(reads,pt.rb<40 & pt.mt<10 & pt.hb<2)
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")


# Method SCTransform -> Harmony
sce <- reads %>% 
  SCTransform() %>% 
  RunPCA(npcs=50, verbose=FALSE) %>% 
  RunHarmony(group.by.vars="orig.ident", 
             assay.use="SCT", max.iter.harmony = 20) %>%
  RunUMAP(reduction="harmony", dims=1:30) %>%
  FindNeighbors(reduction = "umap",dims = 1:2) %>% FindClusters(resolution = c(0.1))
sce <- FindClusters(sce,resolution = c(0.1))
FeaturePlot(sce,c("CDH5","CLDN5"),blend = T,label = T)

GSE155468_EC <- subset(sce,seurat_clusters == 12)
SaveH5Seurat(reads,"../obj/GSE155468_raw",overwrite = T)
SaveH5Seurat(sce,"../obj/GSE155468_SCT",overwrite = T)
SaveH5Seurat(GSE155468_EC,"../obj/GSE155468_EC",overwrite = T)

# GSE155514
rm(list = setdiff(ls(),c("GSEs","Topp.chr_list","mdir")))

i <- 2
setwd(paste0(mdir,"/",GSEs[i]))
files <- dir(recursive = T)
read <- list()
# rownames(tmp)
for (j in 1:length(files)) {
  # j <- 1
  tmp <- as.data.frame(fread(files[j]))
  rnms <- tmp[,1]
  tmp <- as.matrix(tmp[,-1])
  rownames(tmp) <- c(rnms)
  read[[j]] <- CreateSeuratObject(counts = tmp,
                                  project = unlist(str_split(files[j],pattern = "/",n = 2))[1])
}
read <- lapply(read, function(x){
  tmp <- x
  tmp[["pt.mt"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^MT-")
  tmp[["pt.rb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^RP[SL]")
  tmp[["pt.hb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^HB[ABDEGMZ]")
  x <- tmp
})

reads <- read[[1]]
for (k in 2:length(read)) {
  reads <- merge(reads,read[[k]])
}

reads <- RenameCells(reads,add.cell.id = paste0(GSEs[i],"_"))
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")
reads <- subset(reads,pt.rb<40 & pt.mt<10 & pt.hb<2)
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")


# Method SCTransform -> Harmony
sce <- reads %>% 
  SCTransform() %>% 
  RunPCA(npcs=50, verbose=FALSE) %>% 
  RunHarmony(group.by.vars="orig.ident", 
             assay.use="SCT", max.iter.harmony = 20) %>%
  RunUMAP(reduction="harmony", dims=1:30) %>%
  FindNeighbors(reduction = "umap",dims = 1:2) %>% FindClusters(resolution = c(0.1))

FeaturePlot(sce,c("CDH5","CLDN5"),blend = T,label = T)

GSE155514_EC <- subset(sce,seurat_clusters %in% c(3,7))

SaveH5Seurat(reads,"../obj/GSE155514_raw",overwrite = T)
SaveH5Seurat(sce,"../obj/GSE155514_SCT",overwrite = T)
SaveH5Seurat(GSE155514_EC,"../obj/GSE155514_EC",overwrite = T)


# GSE159677
rm(list = setdiff(ls(),c("GSEs","Topp.chr_list","mdir")))

library(hdf5r)
i <- 3
setwd(paste0(mdir,"/",GSEs[i]))
files <- dir(recursive = T)
newdir <- as.data.frame(stringr::str_split_fixed(string = files,pattern = "/",n=4))
newdir$V1 <- mapvalues(newdir$V1,
                       from = levels(factor(newdir$V1)),
                       to = c("P1AC","P1PA","P2AC","P2PA","P3AC","P3PA"))#paste0(paste0("P",1:3),c("AC","PA"))

dirs <- c("P1AC","P1PA","P2AC","P2PA","P3AC","P3PA")
for (dir in 1:length(dirs)) {
  dir.create(dirs[dir])
}

newdirs <- paste0(newdir$V1,"/",newdir$V4)
file.rename(from = files,to = newdirs)
files <- list.files()[-1]

read <- list()
for (j in 1:length(files)) {
  tmp <- Read10X(files[j])
  read[[j]] <- CreateSeuratObject(counts = tmp,
                                  project = files[j])
}
read <- lapply(read, function(x){
  tmp <- x
  tmp[["pt.mt"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^MT-")
  tmp[["pt.rb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^RP[SL]")
  tmp[["pt.hb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^HB[ABDEGMZ]")
  x <- tmp
})

reads <- read[[1]]
for (k in 2:length(read)) {
  reads <- merge(reads,read[[k]])
}

reads <- RenameCells(reads,add.cell.id = paste0(GSEs[i],"_"))
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")
reads <- subset(reads,pt.rb<40 & pt.mt<10 & pt.hb<2)
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")


# Method2 SCTransform -> Harmony
sce <- reads %>% 
  SCTransform() %>% 
  RunPCA(npcs=50, verbose=FALSE) %>% 
  RunHarmony(group.by.vars="orig.ident", 
             assay.use="SCT", max.iter.harmony = 20) %>%
  RunUMAP(reduction="harmony", dims=1:30) %>%
  FindNeighbors(reduction = "umap",dims = 1:2) %>% FindClusters(resolution = c(0.1))
# sce <- FindClusters(sce,resolution = c(0.1))
FeaturePlot(sce,c("CDH5","CLDN5"),blend = T,label = T)

GSE159677_EC <- subset(sce,seurat_clusters %in% c(3,10))

SaveH5Seurat(reads,"../obj/GSE159677_raw",overwrite = T)
SaveH5Seurat(sce,"../obj/GSE159677_SCT",overwrite = T)
SaveH5Seurat(GSE159677_EC,"../obj/GSE159677_EC",overwrite = T)

# GSE213740
rm(list = setdiff(ls(),c("GSEs","Topp.chr_list","mdir")))

i <- 4
setwd(paste0(mdir,"/",GSEs[i]))
files <- list.files()
read <- list()
# rownames(tmp)
for (j in 1:length(files)) {
  # j <- 1
  tmp <- Read10X(files[j])
  read[[j]] <- CreateSeuratObject(counts = tmp,
                                  project = files[j])
}

read <- lapply(read, function(x){
  tmp <- x
  tmp[["pt.mt"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^MT-")
  tmp[["pt.rb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^RP[SL]")
  tmp[["pt.hb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^HB[ABDEGMZ]")
  x <- tmp
})

reads <- read[[1]]
for (k in 2:length(read)) {
  reads <- merge(reads,read[[k]])
}

reads <- RenameCells(reads,add.cell.id = paste0(GSEs[i],"_"))
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")
readsa <- subset(reads,pt.rb<40 & pt.mt<10 & pt.hb<2)
VlnPlot(readsa,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")


# Method SCTransform -> Harmony
sce <- readsa %>% 
  SCTransform() %>% 
  RunPCA(npcs=50, verbose=FALSE) %>% 
  RunHarmony(group.by.vars="orig.ident", 
             assay.use="SCT", max.iter.harmony = 20) %>%
  RunUMAP(reduction="harmony", dims=1:30) %>%
  FindNeighbors(reduction = "umap",dims = 1:2) %>% FindClusters(resolution = c(0.1))

FeaturePlot(sce,c("CDH5","CLDN5"),blend = T,label = T)
VlnPlot(sce,c("CDH5","CLDN5","TOP2A"),stack = T,flip = T)
FeaturePlot(sce,c("CDH5"),blend = T,label = T)


cosgs <- COSG::cosg(sce,assay = "RNA")
cosg_list <- as.list(cosgs$names)
Topp_list <- list()
Topp_list <- lapply(cosg_list, Topp.chr_list)

Topp_anno <- lapply(Topp_list, function(x){
  x <- as.data.frame(stringr::str_split_fixed(x[["Name"]],pattern = "\\|",n = 2))[1:20,]
})


GSE213740_EC <- subset(sce,seurat_clusters %in% c(3,12,15))
SaveH5Seurat(reads,"../obj/GSE213740_raw",overwrite = T)
SaveH5Seurat(sce,"../obj/GSE213740_SCT",overwrite = T)
SaveH5Seurat(GSE213740_EC,"../obj/GSE213740_EC",overwrite = T)

# GSE216860
rm(list = setdiff(ls(),c("GSEs","Topp.chr_list","mdir")))

i <- 5
setwd(paste0(mdir,"/",GSEs[i]))
files <- list.files()
read <- list()
# rownames(tmp)
for (j in 1:length(files)) {
  # j <- 1
  tmp <- Read10X(files[j])
  read[[j]] <- CreateSeuratObject(counts = tmp,
                                  project = files[j])
}

read <- lapply(read, function(x){
  tmp <- x
  tmp[["pt.mt"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^MT-")
  tmp[["pt.rb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^RP[SL]")
  tmp[["pt.hb"]] <- PercentageFeatureSet(tmp, assay = "RNA", pattern = "^HB[ABDEGMZ]")
  x <- tmp
})

reads <- read[[1]]
for (k in 2:length(read)) {
  reads <- merge(reads,read[[k]])
}

reads <- RenameCells(reads,add.cell.id = paste0(GSEs[i],"_"))
VlnPlot(reads,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")
readsa <- subset(reads,pt.rb<40 & pt.mt<10 & pt.hb<2)
VlnPlot(readsa,features = c("pt.mt","pt.rb","pt.hb"),group.by = "orig.ident")


# Method SCTransform -> Harmony
sce <- readsa %>% 
  SCTransform() %>% 
  RunPCA(npcs=50, verbose=FALSE) %>% 
  RunHarmony(group.by.vars="orig.ident", 
             assay.use="SCT", max.iter.harmony = 20) %>%
  RunUMAP(reduction="harmony", dims=1:30) %>%
  FindNeighbors(reduction = "umap",dims = 1:2) %>% FindClusters(resolution = c(0.1))

FeaturePlot(sce,c("CDH5","CLDN5"),blend = T,label = T)
VlnPlot(sce,c("CDH5","CLDN5","TOP2A"),stack = T,flip = T)
FeaturePlot(sce,c("CDH5"),blend = T,label = T)


cosgs <- COSG::cosg(sce,assay = "RNA")
cosg_list <- as.list(cosgs$names)
Topp_list <- list()
Topp_list <- lapply(cosg_list, Topp.chr_list)

Topp_anno <- lapply(Topp_list, function(x){
  x <- as.data.frame(stringr::str_split_fixed(x[["Name"]],pattern = "\\|",n = 2))[1:20,]
})


GSE216860_EC <- subset(sce,seurat_clusters %in% c(1,8,13))
SaveH5Seurat(reads,"../obj/GSE216860_raw",overwrite = T)
SaveH5Seurat(sce,"../obj/GSE216860_SCT",overwrite = T)
SaveH5Seurat(GSE216860_EC,"../obj/GSE216860_EC",overwrite = T)



































# # Method1 Scale -> Normalize -> Harmony
# scRNA <- reads %>% ScaleData() %>% NormalizeData()  %>% FindVariableFeatures() %>% RunPCA()
# scRNA <- scRNA %>% RunHarmony(group.by.vars = "orig.ident",max.iter.harmony = 20) %>% 
#   RunUMAP(reduction = "harmony",dims = 1:30)
# scRNA <- scRNA %>% FindNeighbors(reduction = "umap",dims = 1:2) %>% FindClusters(resolution = c(0.1,0.4))
# DimPlot(scRNA,group.by = "orig.ident",label = T)
# 
# DimPlot(sce,group.by = "orig.ident",label = T)
# 
# FeaturePlot(sce,features = c("MYH11","CDH5","DCN","PTPRC"))
# FeaturePlot(scRNA,features = c("MYH11","CDH5","DCN","PTPRC"))
# FeaturePlot(scRNA,features = c("CXCR3","CD3D"))
# 
# DimPlot(scRNA,group.by = "RNA_snn_res.0.1",label = T)
# DimPlot(sce,group.by = "SCT_snn_res.0.1",label = T)
# 
# Idents(scRNA) <- "RNA_snn_res.0.1"
# 
# tmp <- cosg(scRNA,assay = "RNA")
# cosg_list <- as.list(tmp$names)
# Topp_list <- list()
# Topp_list <- lapply(cosg_list, Topp.chr_list)
# Topp_anno <- lapply(Topp_list, function(x){
#   x <- as.data.frame(stringr::str_split_fixed(x[["Name"]],pattern = "\\|",n = 2))[1:20,]
# })
# 
# scRNA@assays[["RNA"]]@scale.data <- as.matrix(0)
# SaveH5Seurat(scRNA,"GSE155468",overwrite = T)
# SaveH5Seurat(sce,"GSE155468_SCT",overwrite = T)
# 
# gene <- as.data.frame(scRNA@assays[["RNA"]]@counts@Dimnames[[1]])
# grep("^HLA-",gene[,1],value = T)
# 
# VlnPlot(scRNA,c("LYZ","S100A8","CD3G","NKG7","CD79A","CD79B","FLT3","GATA3",
#                 "CD34","ACTC1","DCN","RYR2","TOP2A","PTPRC","CDH5","MYH11"),stack = T,flip=T)
# 
# FeaturePlot(scRNA,c("MS4A1","MS4A2"),blend = T)
# FeaturePlot(scRNA,c("DCN","PDGFRB"),blend = T)
# FeaturePlot(scRNA,c("NKG7","KLRB1"),blend = T)
# 
# 
# scRNA[["pt.hla"]] <- PercentageFeatureSet(scRNA, assay = "RNA", pattern = "^HLA-")
# Nebulosa::plot_density(scRNA,c("pt.hla"))
# df <- scRNA@meta.data
# df$newtype <- df$RNA_snn_res.0.1
# df$newtype <- mapvalues(df$newtype,
#                         from = 0:16,
#                         to = c("T"), # 0
#                         c("Mye"), # 1 
#                         c("SMC"), # 2 
#                         c("T"), # 3 
#                         c("Mye"), # 4 
#                         c("T"), # 5 
#                         c("NK"), # 6 
#                         c("Fib"), # 7 
#                         c("Mye"), # 8 
#                         c("SMC"), # 9 
#                         c(""), # 10
#                         c("SMC"), # 11 
#                         c(""), # 12 
#                         c(""), # 13 
#                         c(""), # 14 
#                         c("Neu"), # 15 
#                         c("")  # 16
#                         )
# df$newtype <- mapvalues(df$newtype,
#                         from = 0:16,
#                         to = c(""), # 0
#                         c(""), # 1 
#                         c(""), # 2 
#                         c(""), # 3 
#                         c(""), # 4 
#                         c(""), # 5 
#                         c(""), # 6 
#                         c(""), # 7 
#                         c(""), # 8 
#                         c(""), # 9 
#                         c(""), # 10
#                         c(""), # 11 
#                         c(""), # 12 
#                         c(""), # 13 
#                         c(""), # 14 
#                         c(""), # 15 
#                         c("")  # 16
# )


