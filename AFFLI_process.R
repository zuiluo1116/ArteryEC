#scFEA

GSE155468_FLUX <- read.csv('./scFEA/GSE155468_flux.csv', header = T, row.names = 1)
GSE155468_FLUX <- t(data.matrix(GSE155468_FLUX))
dim(GSE155468_FLUX)

obj <- subset(scw,GSE=="GSE155468")
obj[["FLUX"]] <- CreateAssayObject(counts = GSE155468_FLUX)
DefaultAssay(obj) <- 'FLUX'
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = F)
obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = F)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 10, reduction.name = 'pca.flux', verbose = F)
obj <- FindNeighbors(obj, dims = 1:2, verbose = F)
obj <- FindClusters(object = obj,graph.name = "SCT_snn")
obj <- RunTSNE(obj, dims = 1:2, assay = 'FLUX', reduction.name = "tsne.flux", verbose = F)
DimPlot(obj, reduction = "tsne.flux",group.by = "scissor_LSS",label = T) + ggtitle('tSNE of Flux')
feamarkers <- FindAllMarkers(object = obj,only.pos = T,logfc.threshold = 0)

obj@assays$SCT <- NULL

GSEs <- levels(factor(GSEs))
scFEA <- list()
getwd()
setwd("E:/works/artery ECs/scFEA")
read.csv("GSE155468_flux.csv", header = T, row.names = 1)

rm(flux_tmp)
scFEA <- t(matrix(rownames(GSE155468_FLUX)))

for (i in 1:length(levels(factor(GSEs)))) {
  t <- read.csv(paste0(GSEs[i],"_flux.csv"), header = T, row.names = 1)
  t <- data.matrix(t)
  scFEA<- rbind(scFEA,t)
  # print(paste0(GSEs[i],"_flux.csv"))
}
scFEA <- scFEA[-1,]

scFEA1 <- apply(scFEA, MARGIN = c(1,2),FUN = function(x){return(as.numeric(x))})
# scl <- scw
rm(scl)
scw[["scFEA"]] <- NULL
scw[["scFEA"]] <- CreateAssayObject(counts = t(scFEA1))
head(scw@assays$RNA@counts)[1:4,1:4]

obj <- DietSeurat(scw)
obj[["FLUX"]] <- CreateAssayObject(counts = t(scFEA1))

# %%R
DefaultAssay(obj) <- 'FLUX'
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 20, verbose = F)
obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = F)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 10, reduction.name = 'pca.flux', verbose = F)
obj <- FindNeighbors(obj, dims = 1:2, verbose = F,reduction = "pca.flux")
obj <- FindClusters(obj, resolution = 0.5, verbose = F)
obj <- RunTSNE(obj, dims = 1:2, assay = 'FLUX', reduction.name = "tsne.flux",reduction = "pca.flux" ,verbose = F)
DimPlot(obj, reduction = "tsne.flux",group.by = "X3d.clu",split.by = "GSE")
Idents(obj) <- "X3d.clu"
fams_fea <- FindAllMarkers(obj,assay = "FLUX",only.pos = F,slot = "scale.data",
                           logfc.threshold = 0,min.diff.pct = 0,
                           min.pct = 0,min.cells.feature = 0,min.cells.group = 0)
avg_fea <- AverageExpression(obj,assays = "FLUX",group.by = "X3d.clu",slot = "scale.data")[[1]]

library(scRNAtoolVis)
library(scales)
library(ComplexHeatmap)
markerVocalno(markers = fams_fea,labelCol = hue_pal()(10))
ComplexHeatmap::Heatmap(avg_fea)


library(SeuratDisk)
library(Seurat)
rm("GSE216860.list")
GSE <- LoadH5Seurat("GSE216860.h5seurat")

DimPlot(GSE)
FeaturePlot(GSE,features = "CDH5",label = T,label.color = "red",label.size = 12)
FeaturePlot(GSE,features = "TAGLN",label = T,label.color = "red",label.size = 12)
FeaturePlot(GSE,features = "PTPRC",label = T,label.color = "red",label.size = 12)
FeaturePlot(GSE,features = "DCN",label = T,label.color = "red",label.size = 12)

GSE$seurat_clusters <- paste0(scw$GSE,"_",scw$orig.ident)
df.newstate <- data.frame(
  clu = GSE$seurat_clusters
)
table(df.newstate$clu)

df.newstate[which(df.newstate$clu %in% c(6,9,13,18,20,22,26,29,32,43)),2] <- "EC"
df.newstate[which(df.newstate$clu %in% c(0,7,10,11,14,15,16,23,25,34,35,39)),2] <- "SMC"
df.newstate[which(df.newstate$clu %in% c(3,4,5,17,19,21,27,30,31,33,38,40,42)),2] <- "IMM"
df.newstate[which(df.newstate$clu %in% c(1,2,8,12,24,28)),2] <- "Fib"

GSE$celltype <- df.newstate$V2
DimPlot(GSE,group.by = "celltype")
prop.table(table(GSE$celltype,GSE$orig.ident),margin = 2)

levels(factor(df.newstate$clu[is.na(df.newstate$V2)]))
rm(list=ls())

#################################################################################
# > load("E:/works/artery ECs/monocle.Rdata")
# > load("E:/works/artery ECs/vas_EC.Rdata")
cdst <- cdsk_3d[,select.cells]
tmp <- ls()
tmp <- tmp[-which(tmp=="cdst")]
rm(list = tmp)

cdsk_3d <- learn_graph(cdsk_3d)
cdst = cluster_cells(cdst)

cdst <- learn_graph(cdst)
plot_cells(cdst)
cdst <- order_cells(cdst)
Track_genes_NEW <- graph_test(cdst, neighbor_graph="principal_graph")
library(tidyverse)
save(cds,file = "cds1222.RData")
load("NEW.clu.RData")
df.3d <- df.3d[scw@assays[["RNA"]]@counts@Dimnames[[2]],]
scw$NEW.clu <- df.3d$clusters


df.3d[which(df.3d$clusters %in% c(1,2,4,6,9,10)),3] <- "ACKR1+ EC"
df.3d[which(df.3d$clusters %in% c(5,7)),3] <- "TCIM+ EC"
df.3d[which(df.3d$clusters %in% c(3)),3] <- "ITLN1+ EC"
df.3d[which(df.3d$clusters %in% c(8)),3] <- "Lymphatic EC"
names(df.3d)[3] <- "celltype"
scw$celltype <- df.3d$celltype
scw$celltype <- factor(df.3d$celltype,levels = c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))

DimPlot(scw,group.by = "celltype",reduction = "scphere")

library(monocle3)
data<-GetAssayData(scw,assay ='RNA',slot ='counts')
cell_metadata <-scw@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(data))
rownames(gene_annotation)<-rownames(data)
cds <-new_cell_data_set(data,cell_metadata =cell_metadata,gene_metadata =gene_annotation)
# 主成分
cds <- preprocess_cds(cds, num_dim = 50)
# plot_pc_variance_explained(cds)
# 降维 默认UMAP
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_cells(cds,show_trajectory_graph = F,color_cells_by = "celltype")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scw, reduction = "scphere")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds,reduction_method = "UMAP")
# plot_cells(cds,show_trajectory_graph = F,color_cells_by = "cluster")

cds <- learn_graph(cds)
# 给定起点绘制轨迹
cds <- order_cells(cds)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph")

# Track_genes.fli1 <- Track_genes2[,c(5,2,3,4,1,6)] %>% 
#   filter(q_value < 1e-3) %>% 
#   arrange(desc(morans_I))
Track_genes.fli <- Track_genes[,c(5,2,3,4,1,6)] %>% 
  filter(q_value < 1e-3) %>% 
  arrange(desc(morans_I))
genelist <- pull(Track_genes.fli, gene_short_name) %>% as.character()
set.seed(1116)
gene_module <- find_gene_modules(cds[genelist,], resolution=5e-3, cores = 1)
gene_module.fli <- gene_module[1:50,]
gene_module.fli <- gene_module[1:200,]

cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=cds@colData@listData[["celltype"]])
cell_group2 <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=cds@colData@listData[["NEW.clu"]])
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
agg_mat2 <- aggregate_gene_expression(cds, gene_module, cell_group2)
# write.csv(Track_genes.fli, "Trajectory_genes.fli.csv", row.names = F)
genelist <- pull(Track_genes.fli, gene_short_name) %>% as.character()
set.seed(1116)

gene_module.fli <- gene_module[1:200,]
### by celltype
gene_module.df <- melt(table(gene_module.fli$module))
names(gene_module.df) <- c("clusters","count")
gene_module.df <- gene_module.df[order(gene_module.df$count,decreasing = T),]
agg_mat <- agg_mat[as.character(gene_module.df$clusters),]
# agg_mat2 <- agg_mat2[as.character(gene_module.df$clusters),]
gene_module.df <- cbind(gene_module.df,agg_mat)

for (i in 1:nrow(gene_module.df)) {
  # gene_module.df[i,7] <- max(gene_module.df[i,3:6])
  gene_module.df[i,7] <- names(which.max(gene_module.df[i,3:6]))
}

for (i in 1:nrow(gene_module.df)) {
  l <- as.numeric(which.max(gene_module.df[i,3:6]))+2
  k <- setdiff(3:6,l)
  gene_module.df[i,8] <- (gene_module.df[i,l]-gene_module.df[i,k[1]])+
                          (gene_module.df[i,l]-gene_module.df[i,k[2]])+
                          (gene_module.df[i,l]-gene_module.df[i,k[3]])

}


names(gene_module.df)[7] <- "Enriched.Cluster"
names(gene_module.df)[8] <- "all.diff"
gene_module.df$Enriched.Cluster <- factor(gene_module.df$Enriched.Cluster,
                                          levels = c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))



# rownames(agg_mat)
# pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",
#                    cellwidth = 8, cellheight = 5,border_color = "white",
#                    cluster_cols = F)

agg_sel <- agg_mat[rownames(gene_module.df)[which(gene_module.df$count!=0)],]
# agg_sel2 <- agg_mat2[rownames(gene_module.df)[which(gene_module.df$count!=0)],]
row_anno <- gene_module.df[rownames(agg_sel),c(2,7,8)]
col_anno <- as.data.frame(c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))
rownames(col_anno) <- col_anno[,1]
names(col_anno) <- "Cell.Type"
names(row_anno) <- c("Count","Enriched.Cluster","All.Diff")
ann_colors = list(
  Count = c("white", "firebrick"),
  Enriched.Cluster = c(`ITLN1+ EC` = "#F8766D",
                       `TCIM+ EC` = "#7CAE00",
                       `ACKR1+ EC` = "#00BFC4",
                       `Lymphatic EC` = "#C77CFF"),
  Cell.Type = c(`ITLN1+ EC` = "#F8766D",
                       `TCIM+ EC` = "#7CAE00",
                       `ACKR1+ EC` = "#00BFC4",
                       `Lymphatic EC` = "#C77CFF"),
  All.Diff = c("#e9e9e9","#7570B3")
)
row_anno$Enriched.Cluster <- factor(row_anno$Enriched.Cluster,
                                    levels = c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))
pdf( 
  file = "./Figures/Fig.6A.pdf", # 文件名称
  width = 6,           # 宽
  height = 6)
ComplexHeatmap::pheatmap(agg_sel, scale="none", clustering_method="ward.D2",
                   cellwidth = 18, cellheight = 10,border_color = "white",
                   cluster_cols = T,annotation_row = row_anno,annotation_col = col_anno,
                   annotation_colors = ann_colors)
dev.off()


#### 
gene_module.df <- melt(table(gene_module.fli$module))
names(gene_module.df) <- c("clusters","count")
gene_module.df <- gene_module.df[order(gene_module.df$count,decreasing = T),]
agg_mat2 <- agg_mat2[as.character(gene_module.df$clusters),]
gene_module.df <- cbind(gene_module.df,agg_mat2)

for (i in 1:nrow(gene_module.df)) {
  # gene_module.df[i,7] <- max(gene_module.df[i,3:6])
  gene_module.df[i,13] <- names(which.max(gene_module.df[i,3:12]))
}

for (i in 1:nrow(gene_module.df)) {
  l <- as.numeric(which.max(gene_module.df[i,3:12]))+2
  k <- setdiff(3:12,l)
  gene_module.df[i,14] <- (gene_module.df[i,l]-gene_module.df[i,k[1]])+
    (gene_module.df[i,l]-gene_module.df[i,k[2]])+
    (gene_module.df[i,l]-gene_module.df[i,k[3]])+
    (gene_module.df[i,l]-gene_module.df[i,k[4]])+
    (gene_module.df[i,l]-gene_module.df[i,k[5]])+
    (gene_module.df[i,l]-gene_module.df[i,k[6]])+
    (gene_module.df[i,l]-gene_module.df[i,k[7]])+
    (gene_module.df[i,l]-gene_module.df[i,k[8]])+
    (gene_module.df[i,l]-gene_module.df[i,k[9]])
}


names(gene_module.df)[13] <- "Enriched.Cluster"
names(gene_module.df)[14] <- "all.diff"
# gene_module.df$Enriched.Cluster <- ifelse(gene_module.df$Enriched.Cluster %in% c(1,2,4,6,9,10),"ACKR1+ EC",
#                                           ifelse(gene_module.df$Enriched.Cluster %in% c(5,7),"TCIM+ EC",
#                                                  ifelse(gene_module.df$Enriched.Cluster %in% c(3),"ITLN1+ EC","Lymphatic EC")))



# rownames(agg_mat)
# pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",
#                    cellwidth = 8, cellheight = 5,border_color = "white",
#                    cluster_cols = F)

agg_sel2 <- agg_mat2[rownames(gene_module.df)[which(gene_module.df$count!=0)],]
# agg_sel2 <- agg_mat2[rownames(gene_module.df)[which(gene_module.df$count!=0)],]
row_anno <- gene_module.df[rownames(agg_sel2),c(2,13,14)]
col_anno <- melt(table(scw$celltype,scw$NEW.clu))
col_anno <- col_anno[which(col_anno$value!=0),]
rownames(col_anno) <- col_anno$Var2
col_anno <- as.data.frame(col_anno[,1])
names(col_anno) <- "Cell.Type"
names(row_anno) <- c("Count","Enriched.Cluster","All.Diff")
ann_colors = list(
  Count = c("white", "firebrick"),
  Cell.Type= c(`ITLN1+ EC` = "#F8766D",
                       `TCIM+ EC` = "#7CAE00",
                       `ACKR1+ EC` = "#00BFC4",
                       `Lymphatic EC` = "#C77CFF"),
  Enriched.Cluster= c(`3` = "#F8766D",
                `5` = "#7CAE00",
                `7` = "#7CAE00",
                `1` = "#00BFC4",
                `2` = "#00BFC4",
                `4` = "#00BFC4",
                `6` = "#00BFC4",
                `9` = "#00BFC4",
                `10` = "#00BFC4",
                `8` = "#C77CFF"),
  All.Diff = c("#e9e9e9","#7570B3")
)
# row_anno$Enriched.Cluster <- factor(row_anno$Enriched.Cluster,
#                                     levels = c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))
pdf( 
  file = "./Figures/sFig.4A.pdf", # 文件名称
  width = 6,           # 宽
  height = 6)
ComplexHeatmap::pheatmap(agg_sel2, scale="none", clustering_method="ward.D2",
                         cellwidth = 18, cellheight = 10,border_color = "white",
                         cluster_cols = T,annotation_row = row_anno,annotation_col = col_anno,
                         annotation_colors = ann_colors)
dev.off()



save(gene_module,file = "gene_module1223.RData")
gene_module.list <- split(gene_module$id,gene_module$module)
gene_module.sel <- gene_module[which(gene_module$module %in% c(23,17,34,25)),1:2]

plot_cells(cds,genes = gene_module.sel,show_trajectory_graph = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)
library(RobustRankAggreg)

module.23 <- unlist(c(gene_module[which(gene_module$module == 23),1]))
module.23.con <- bitr(module.23,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.17 <- unlist(c(gene_module[which(gene_module$module == 17),1]))
module.17.con <- bitr(module.17,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.25 <- unlist(c(gene_module[which(gene_module$module == 25),1]))
module.25.con <- bitr(module.25,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.34 <- unlist(c(gene_module[which(gene_module$module == 34),1]))
module.34.con <- bitr(module.34,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.24 <- unlist(c(gene_module[which(gene_module$module == 24),1]))
module.24.con <- bitr(module.24,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.55 <- unlist(c(gene_module[which(gene_module$module == 55),1]))
module.55.con <- bitr(module.55,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.102 <- unlist(c(gene_module[which(gene_module$module == 102),1]))
module.102.con <- bitr(module.102,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.4 <- unlist(c(gene_module[which(gene_module$module == 4),1]))
module.4.con <- bitr(module.4,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)

ChEA.TFs.RRA <- function(x){
  url = "https://maayanlab.cloud/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = "myQuery", gene_set = x)
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  TFs = fromJSON(json)
  TFs.only <- TFs
  for (i in 1:length(TFs)) {
    TFs.only[i] <- dplyr::select(TFs[[i]],"TF")
  }
  freq=as.data.frame(table(unlist(TFs.only)))
  ag=aggregateRanks(TFs.only)
  ag$Freq=freq[match(ag$Name,freq$Var1),2]
  return(ag)
}

module4.TFs.RRA <- ChEA.TFs.RRA(module.4.con$SYMBOL)
module102.TFs.RRA <- ChEA.TFs.RRA(module.102.con$SYMBOL)
module55.TFs.RRA <- ChEA.TFs.RRA(module.55.con$SYMBOL)
module17.TFs.RRA <- ChEA.TFs.RRA(module.17.con$SYMBOL)
module23.TFs.RRA <- ChEA.TFs.RRA(module.23.con$SYMBOL)
module25.TFs.RRA <- ChEA.TFs.RRA(module.25.con$SYMBOL)
module34.TFs.RRA <- ChEA.TFs.RRA(module.34.con$SYMBOL)
module24.TFs.RRA <- ChEA.TFs.RRA(module.24.con$SYMBOL)



library(UpSetR)


load("FAMS_type.RData")

FAMs_ITLN1 <- FAMs_type %>% 
  filter(cluster == "ITLN1+ EC") %>% pull( gene) %>% as.character()
FAMs_ITLN1.con <- bitr(FAMs_ITLN1,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
FAMs_TCIM <- FAMs_type %>% 
  filter(cluster == "TCIM+ EC") %>% pull( gene) %>% as.character()
FAMs_TCIM.con <- bitr(FAMs_TCIM,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
FAMs_ACKR1 <- FAMs_type %>% 
  filter(cluster == "ACKR1+ EC") %>% pull( gene) %>% as.character()
FAMs_ACKR1.con <- bitr(FAMs_ACKR1,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
FAMs_Lymphatic <- FAMs_type %>% 
  filter(cluster == "Lymphatic EC") %>% pull( gene) %>% as.character()
FAMs_Lymphatic.con <- bitr(FAMs_Lymphatic,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)


FAMs_ITLN1.TFs.RRA <- ChEA.TFs.RRA(FAMs_ITLN1.con$SYMBOL)
FAMs_TCIM.TFs.RRA <- ChEA.TFs.RRA(FAMs_TCIM.con$SYMBOL)
FAMs_ACKR1.TFs.RRA <- ChEA.TFs.RRA(FAMs_ACKR1.con$SYMBOL)
FAMs_Lymphatic.TFs.RRA <- ChEA.TFs.RRA(FAMs_Lymphatic.con$SYMBOL)


TFs <- list(
  # rownames(module4.TFs.RRA)[1:50],
  # rownames(module102.TFs.RRA)[1:50],
  rownames(module23.TFs.RRA)[1:50],
  # rownames(module24.TFs.RRA)[1:50],
  rownames(module17.TFs.RRA)[1:50],
  rownames(module34.TFs.RRA)[1:50],
  # rownames(module55.TFs.RRA)[1:50],
  rownames(module25.TFs.RRA)[1:50],
  rownames(FAMs_ITLN1.TFs.RRA)[1:50],
  rownames(FAMs_TCIM.TFs.RRA)[1:50],
  rownames(FAMs_ACKR1.TFs.RRA)[1:50],
  rownames(FAMs_Lymphatic.TFs.RRA)[1:50]
)
# names(TFs) <- c("Module 4","Module 102","Module 23","Module 24","Module 17","Module 34","Module 55","Module 25",
#                 "ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC")
names(TFs) <- c("Module 23","Module 17","Module 34","Module 25",
                "ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC")
upset(fromList(TFs),keep.order = T,nsets = 10,order.by = "freq")


datd_tmp <- fromList(TFs)
rownames(datd_tmp) <- unique(unlist(TFs))

upset(datd_tmp,  
      queries = list(list(query = intersects,   #UpSetR 内置的intersects query
                          params = list("Module 34"), ##指定作用的交集
                          color = "red", ##设置颜色，未设置会调用默认调色板
                          active = F,   # TRUE：条形图被颜色覆盖，FALSE：条形图顶端显示三角形
                          query.name = "ITLN1+")))#, # 添加query图例
                     # list(query = intersects,  params = list("Action", "Drama"), active = T,query.name = "Emotional action"), 
                     # list(query = intersects,  params = list("Drama", "Comedy", "Action"), color = "orange", active = T)),query.legend = "top")

# intersect(TFs$`Module 102`,TFs$`ITLN1+ EC`)
# intersect(TFs$`Module 4`,TFs$`ITLN1+ EC`)
intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`)
intersect(TFs$`Module 17`,TFs$`TCIM+ EC`)
intersect(TFs$`Module 24`,TFs$`TCIM+ EC`)
intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`)
# intersect(TFs$`Module 55`,TFs$`ACKR1+ EC`)
# intersect(TFs$`Module 25`,TFs$`Lymphatic EC`)

upset(fromList(TFs),keep.order = T,nsets = 8,order.by = "degree")

rm(ITLN1_EC_TFs)











expr_data <- AverageExpression(scw,assays = "RNA")
expr_data <- expr_data[[1]]
# 1
ITLN1_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`)
)
ITLN1_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`),
  "Modules" = module23.TFs.RRA[match(ITLN1_EC_TFs$TFs,module23.TFs.RRA$Name),2],
  "Uniques" = FAMs_ITLN1.TFs.RRA[match(ITLN1_EC_TFs$TFs,FAMs_ITLN1.TFs.RRA$Name),2]
)
ITLN1_EC_TFs$int <- ITLN1_EC_TFs$Modules*ITLN1_EC_TFs$Uniques
ITLN1_EC_TFs <- arrange_(ITLN1_EC_TFs,"int")

ITLN1_EC_TFs <- ITLN1_EC_TFs[1:10,]
ITLN1_EC_TFs <- ITLN1_EC_TFs[order(ITLN1_EC_TFs$int,decreasing = T),]
ITLN1_EC_TFs$TFs <- factor(ITLN1_EC_TFs$TFs,
                           levels = ITLN1_EC_TFs$TFs)

ITLN1_EC_TFs$Modules <- -log10(ITLN1_EC_TFs$Modules)
ITLN1_EC_TFs$Uniques <- log10(ITLN1_EC_TFs$Uniques)

tmp <- melt(ITLN1_EC_TFs[,1:3])

ggplot(tmp, aes(x = TFs, y = value, fill = variable)) +
  geom_bar(stat="identity", color = 'white', alpha = 1, width = 0.95)+
  # scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0),limits = c(-14,14)) +
  theme_classic()+coord_flip()+NoLegend()
  # +
  # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
Nebulosa::plot_density(scw,features = ITLN1_EC_TFs$TFs,reduction = "scphere")

expr_tf <- expr_data[as.character(ITLN1_EC_TFs$TFs),]
expr_tf <- t(scale(t(expr_tf)))
ComplexHeatmap::pheatmap(expr_tf, scale="none", clustering_method="ward.D2",
                         cellwidth = 18, cellheight = 10,border_color = "white",
                         cluster_cols = F,#annotation_row = row_anno,annotation_col = col_anno,
                         annotation_colors = ann_colors)




# 2
# intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`)
intersect(TFs$`Module 17`,TFs$`TCIM+ EC`)
# intersect(TFs$`Module 24`,TFs$`TCIM+ EC`)
# intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`)
TCIM_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 17`,TFs$`TCIM+ EC`)
)
TCIM_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 17`,TFs$`TCIM+ EC`),
  "Modules" = module17.TFs.RRA[match(TCIM_EC_TFs$TFs,module17.TFs.RRA$Name),2],
  "Uniques" = FAMs_TCIM.TFs.RRA[match(TCIM_EC_TFs$TFs,FAMs_TCIM.TFs.RRA$Name),2]
)
TCIM_EC_TFs$int <- TCIM_EC_TFs$Modules*TCIM_EC_TFs$Uniques
TCIM_EC_TFs <- arrange_(TCIM_EC_TFs,"int")

TCIM_EC_TFs <- TCIM_EC_TFs[1:10,]
TCIM_EC_TFs <- TCIM_EC_TFs[order(TCIM_EC_TFs$int,decreasing = T),]
TCIM_EC_TFs$TFs <- factor(TCIM_EC_TFs$TFs,
                           levels = TCIM_EC_TFs$TFs)

TCIM_EC_TFs$Modules <- -log10(TCIM_EC_TFs$Modules)
TCIM_EC_TFs$Uniques <- log10(TCIM_EC_TFs$Uniques)

tmp <- melt(TCIM_EC_TFs[,1:3])

ggplot(tmp, aes(x = TFs, y = value, fill = variable)) +
  geom_bar(stat="identity", color = 'white', alpha = 1, width = 0.95)+
  # scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0),limits = c(-14,14)) +
  theme_classic()+coord_flip()+NoLegend()
# +
# theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
Nebulosa::plot_density(scw,features = TCIM_EC_TFs$TFs,reduction = "scphere")

Idents(scw) <- "celltype"
scRNAtoolVis::AverageHeatmap(scw,markerGene = TCIM_EC_TFs$TFs)

expr_tf <- expr_data[as.character(TCIM_EC_TFs$TFs),]
expr_tf <- t(scale(t(expr_tf)))
ComplexHeatmap::pheatmap(expr_tf, scale="none", clustering_method="ward.D2",
                         cellwidth = 18, cellheight = 10,border_color = "white",
                         cluster_cols = F,#annotation_row = row_anno,annotation_col = col_anno,
                         annotation_colors = ann_colors)

  
# 3
# intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`)
# intersect(TFs$`Module 17`,TFs$`ACKR1+ EC`)
intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`)
# intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`)
ACKR1_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`)
)
ACKR1_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`),
  "Modules" = module34.TFs.RRA[match(ACKR1_EC_TFs$TFs,module34.TFs.RRA$Name),2],
  "Uniques" = FAMs_ACKR1.TFs.RRA[match(ACKR1_EC_TFs$TFs,FAMs_ACKR1.TFs.RRA$Name),2]
)
ACKR1_EC_TFs$int <- ACKR1_EC_TFs$Modules*ACKR1_EC_TFs$Uniques
ACKR1_EC_TFs <- arrange_(ACKR1_EC_TFs,"int")

ACKR1_EC_TFs <- ACKR1_EC_TFs[1:10,]
ACKR1_EC_TFs <- ACKR1_EC_TFs[order(ACKR1_EC_TFs$int,decreasing = T),]
ACKR1_EC_TFs$TFs <- factor(ACKR1_EC_TFs$TFs,
                          levels = ACKR1_EC_TFs$TFs)

ACKR1_EC_TFs$Modules <- -log10(ACKR1_EC_TFs$Modules)
ACKR1_EC_TFs$Uniques <- log10(ACKR1_EC_TFs$Uniques)

tmp <- melt(ACKR1_EC_TFs[,1:3])

ggplot(tmp, aes(x = TFs, y = value, fill = variable)) +
  geom_bar(stat="identity", color = 'white', alpha = 1, width = 0.95)+
  # scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0),limits = c(-14,14)) +
  theme_classic()+coord_flip()+NoLegend()
# +
# theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
Nebulosa::plot_density(scw,features = ACKR1_EC_TFs$TFs,reduction = "scphere")

Idents(scw) <- "celltype"
scRNAtoolVis::AverageHeatmap(scw,markerGene = ACKR1_EC_TFs$TFs)
expr_tf <- expr_data[as.character(ACKR1_EC_TFs$TFs),]
expr_tf <- t(scale(t(expr_tf)))
ComplexHeatmap::pheatmap(expr_tf, scale="none", clustering_method="ward.D2",
                         cellwidth = 18, cellheight = 10,border_color = "white",
                         cluster_cols = F,#annotation_row = row_anno,annotation_col = col_anno,
                         annotation_colors = ann_colors)

# 4
# intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`)
# intersect(TFs$`Module 17`,TFs$`Lymphatic+ EC`)
intersect(TFs$`Module 25`,TFs$`Lymphatic EC`)
# intersect(TFs$`Module 34`,TFs$`Lymphatic+ EC`)
Lymphatic_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 34`,TFs$`Lymphatic EC`)
)
Lymphatic_EC_TFs <- data.frame(
  "TFs" = intersect(TFs$`Module 34`,TFs$`Lymphatic EC`),
  "Modules" = module34.TFs.RRA[match(Lymphatic_EC_TFs$TFs,module34.TFs.RRA$Name),2],
  "Uniques" = FAMs_Lymphatic.TFs.RRA[match(Lymphatic_EC_TFs$TFs,FAMs_Lymphatic.TFs.RRA$Name),2]
)
Lymphatic_EC_TFs$int <- Lymphatic_EC_TFs$Modules*Lymphatic_EC_TFs$Uniques
Lymphatic_EC_TFs <- arrange_(Lymphatic_EC_TFs,"int")

Lymphatic_EC_TFs <- Lymphatic_EC_TFs[1:10,]
Lymphatic_EC_TFs <- Lymphatic_EC_TFs[order(Lymphatic_EC_TFs$int,decreasing = T),]
Lymphatic_EC_TFs$TFs <- factor(Lymphatic_EC_TFs$TFs,
                           levels = Lymphatic_EC_TFs$TFs)

Lymphatic_EC_TFs$Modules <- -log10(Lymphatic_EC_TFs$Modules)
Lymphatic_EC_TFs$Uniques <- log10(Lymphatic_EC_TFs$Uniques)
Lymphatic_EC_TFs$Uniques <- -Lymphatic_EC_TFs$Uniques

tmp <- melt(Lymphatic_EC_TFs[,1:3])

ggplot(tmp, aes(x = TFs, y = value, fill = variable)) +
  geom_bar(stat="identity", position = "dodge",color = 'white', alpha = 0.8, width = 0.95)+#coord_polar(start = 0)+
  # scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0),limits = c(-20,14)) +
  ggpubr::theme_pubr()+NoLegend()+coord_polar(start = 0)
  # annotate("segment",x=0.5,xend =0.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=1.5,xend =1.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=2.5,xend =2.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=3.5,xend =3.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=4.5,xend =4.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=5.5,xend =5.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=6.5,xend =6.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=7.5,xend =7.5,y=-1,yend=15,colour="black")+
  # annotate("segment",x=8.5,xend =8.5,y=-1,yend=15,colour="black")+
  # # annotate("segment",x=9.5,xend =9.5,y=-1,yend=15,colour="black")+
  # geom_segment(aes(x=0, y=15,xend=9.5,yend =15),size=1.5,color="#3B9AB2",
  #              arrow = arrow(length = unit(0, "npc"),type="closed"))
  # +
# theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

df_tf_expr <- as.data.frame(expr_tf)
for (i in 1:nrow(df_tf_expr)) {
  # gene_module.df[i,7] <- max(gene_module.df[i,3:6])
  df_tf_expr[i,5] <- names(which.max(df_tf_expr[i,1:4]))
}
names(df_tf_expr)[5] <- "Cluster"
df_tf_expr$Cluster <- factor(df_tf_expr$Cluster,levels = c(levels(as.factor(scw$celltype))))
df_tf_expr <- df_tf_expr[order(df_tf_expr$Cluster),]
df_tf_expr <- as.matrix(df_tf_expr[,1:4])
range(df_tf_expr)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
mycol <- colorRamp2(c(-1.5, 0, 1.5),c("#0da9ce", "white", "#e74a32"))

df_tf_expr <- df_tf_expr[unique(tmp$TFs),]
#环形热图绘制：
circos.heatmap(df_tf_expr,
               col = mycol)
circos.clear()#绘制完成后，需要重置循环布局再开启下一次绘图，不然会遇到报错或轨道继续在同一图中叠加；

circos.par(gap.after = c(0)) #调整圆环首位间距；
circos.heatmap(df_tf_expr,
               col = mycol,
               cluster = F, #是否对行聚类
               # dend.side = "inside", #聚类树方向：inside显示在圆环内圈，inside为圆环外圈；
               rownames.side ="outside", #行名方向；
               rownames.col = "black",
               track.height = 0.6, #轨道高度，即圆环/格子的粗细
               rownames.cex = 1)
circos.clear()
lg <- Legend(title = "Express",
             col_fun = mycol,
             direction = c("vertical"),
             title_position = c('topcenter'))
draw(lg, x = unit(0.9, "npc"), y = unit(0.5, "npc"), just = c("right", "center"))
df_tf_expr2 <- as.data.frame(expr_tf)
df_tf_expr2 <- df_tf_expr2[rownames(df_tf_expr),]
for (i in 1:nrow(df_tf_expr2)) {
  # gene_module.df[i,7] <- max(gene_module.df[i,3:6])
  df_tf_expr2[i,5] <- names(which.max(df_tf_expr2[i,1:4]))
}
df_tf_expr2$TFs <- rownames(df_tf_expr2)
df_tf_expr2 <- split(df_tf_expr2$TFs,df_tf_expr2$V5)

tmp1 <- melt(ITLN1_EC_TFs[,1:3])
tmp1 <- tmp1[which(tmp1$TFs %in% df_tf_expr2$`ITLN1+ EC`),]

tmp2 <- melt(TCIM_EC_TFs[,1:3])
tmp2 <- tmp2[which(tmp2$TFs %in% df_tf_expr2$`TCIM+ EC`),]

tmp3 <- melt(ACKR1_EC_TFs[,1:3])
tmp3 <- tmp3[which(tmp3$TFs %in% df_tf_expr2$`ACKR1+ EC`),]

tmp4 <- melt(Lymphatic_EC_TFs[,1:3])
tmp4 <- tmp4[which(tmp4$TFs %in% df_tf_expr2$`Lymphatic EC`),]

tmp <- rbind(tmp1,tmp2)
tmp <- rbind(tmp,tmp3)
tmp <- rbind(tmp,tmp4)

tmp$TFs <- factor(tmp$TFs,levels = rownames(df_tf_expr))
tmp$value <- ifelse(tmp$value<0,-tmp$value,tmp$value)
ggplot(tmp, aes(x = TFs, y = value, fill = variable)) +
  geom_bar(stat="identity", position = "dodge",
           color = 'white', alpha = 0.8, width = 0.95)+#coord_polar(start = 0)+
  # scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0),limits = c(-20,14)) +
  ggpubr::theme_pubr()+coord_polar(start = 0)

df_tf_expr2

scales::hue_pal()(4)
"#F8766D" "#7CAE00" "#00BFC4" "#C77CFF"



Nebulosa::plot_density(scw,features = Lymphatic_EC_TFs$TFs,reduction = "scphere")

Idents(scw) <- "celltype"
scRNAtoolVis::AverageHeatmap(scw,markerGene = Lymphatic_EC_TFs$TFs)

TFs_sel <- as.character(unique(unlist(list(ACKR1_EC_TFs$TFs,
                                ITLN1_EC_TFs$TFs,
                                TCIM_EC_TFs$TFs,
                                Lymphatic_EC_TFs$TFs))))
expr_tf <- expr_data[TFs_sel,]
expr_tf <- t(scale(t(expr_tf)))
ComplexHeatmap::pheatmap(expr_tf, scale="none", clustering_method="ward.D2",
                         cellwidth = 18, cellheight = 10,border_color = "white",
                         cluster_cols = F,#annotation_row = row_anno,annotation_col = col_anno,
                         annotation_colors = ann_colors)



  
  
  
  
  
  
  
library(ggplot2)
library(tidyverse)

ggplot(table,aes(Cluster,Cell,fill=Sample))+
  geom_col()+
  coord_flip()+
  scale_y_continuous(breaks = seq(-40, 40, 10), labels = as.character(abs(seq(-40, 40, 10))), limits = c(-40, 40))+  #设置y轴合适范围，并将负数转化为正数
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(),axis.text = element_text(color="black",size=12),axis.title = element_text(color = "black",size=15)) +geom_hline(yintercept = 0, size = 0.4) + #添加中间的线，或者不添加也行
  annotate('text',label = 'Y_SYF', 1, 40)+
  annotate('text',label = 'O_SYF', 1, -40)  #添加注释信息




# write.csv(Track_genes.fli, "Trajectory_genes.fli.csv", row.names = F)

Track_genes_sig10 <- Track_genes_NEW.fli %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()


Track_genes_sig100 <- Track_genes.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

cdst_genes_sig100 <- Track_genes_NEW.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Track_genes_NEW.fli <- Track_genes_NEW.fli %>% arrange(desc("morans_I"))
Track_genes.fli <- Track_genes.fli %>% arrange(desc("morans_I"))

df <- data.frame(
  old = Track_genes_sig100,
  new = cdst_genes_sig100
)

setdiff(Track_genes_sig100,cdst_genes_sig100)
setdiff(cdst_genes_sig100,Track_genes_sig100)








devtools::install_github("AckerDWM/gg3D")

library("gg3D")
qplot(x=0, y=0, z=0, geom="blank") + 
  theme_void() +
  axes_3D()
data <- cdst@int_colData@listData[["reducedDims"]]@listData[["UMAP"]]

cdst_track <- graph_test(cdst, neighbor_graph="principal_graph")
save(cdst_track,file = "Track_gene_vasEC.RData")
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# Matrix::rBind() to rbind()






# data <- cdst@principal_graph_aux@listData[["UMAP"]][["dp_mst"]]
# data <- cdst@principal_graph_aux@listData[["UMAP"]][["dp_mst"]]
# data <- as.data.frame(t(data))
names(data) <- c("X","Y","Z")

theta=135 #方位角的度数
phi=0 # 渐近线
ggplot(data, aes(x=X, y=Y, z=Z)) + 
  axes_3D(theta=theta, phi=phi) + stat_3D(theta=theta, phi=phi) +
  axis_labs_3D(theta=theta, phi=phi, size=3, hjust=c(1,1,1.2,1.2,1.2,1.2), vjust=c(-.5,-.5,-.2,-.2,1.2,1.2)) +
  labs_3D(theta=theta, phi=phi, hjust=c(1,0,0), vjust=c(1.5,1,-.2),labs=c("X", "Y", "Z")) +theme_void()


library(ggplot2)

ggplot(data, aes(x = "X", y = "Y", z = "Z"), padding = "padding:250px 500px 100px 100px;")


write.csv(cdst_track,"cdst_track.csv")


head(colData(cds))[1:10]


cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=cds@clusters@listData[["UMAP"]][["clusters"]])







#scFEA

GSEs <- levels(factor(scw$GSE))
flux <- t(matrix(rep(0,168)))
for (i in 1:length(GSEs)) {
  flux_tmp <- read.csv(paste0("./scFEA/",GSEs[i],"_flux.csv"), header = T, row.names = 1)
  flux_tmp <- data.matrix(flux_tmp)
  flux <- rbind(flux,flux_tmp)
}
colnames(flux) <- colnames(flux_tmp)
flux <- flux[scw@assays[["RNA"]]@counts@Dimnames[[2]],]

scw@assays$FLUX <- CreateAssayObject(counts = t(flux))
scw@assays[["FLUX"]]@scale.data <- scale(t(flux))

balance <- t(matrix(rep(0,70)))
for (i in 1:length(GSEs)) {
  balance_tmp <- read.csv(paste0("./scFEA/",GSEs[i],"_balance.csv"), header = T, row.names = 1)
  balance_tmp <- data.matrix(balance_tmp)
  balance <- rbind(balance,balance_tmp)
}
colnames(balance) <- colnames(balance_tmp)
balance <- balance[scw@assays[["RNA"]]@counts@Dimnames[[2]],]

scw@assays$BALANCE <- CreateAssayObject(counts = t(balance))
scw@assays[["BALANCE"]]@scale.data <- scale(t(balance))




DefaultAssay(scw) <- "BALANCE"
DimPlot(scw,reduction = "scphere",label=T,group.by = "NEW.clu")
Idents(scw) <- "NEW.clu"
levels(Idents(scw))
fams_balance <- FindAllMarkers(scw,logfc.threshold = 0,min.pct = 0,only.pos = T)
fams_balance_fli <- fams_balance %>% dplyr::filter(avg_log2FC>log2(1) & p_val <0.05 & pct.1-pct.2>0.05)

DefaultAssay(scw) <- "FLUX"
fams_flux <- FindAllMarkers(scw,logfc.threshold = 0,min.pct = 0,only.pos = T)
fams_flux_fli <- fams_flux %>% dplyr::filter(avg_log2FC>log2(1) & p_val <0.05 & pct.1-pct.2>0.05)
log2(1.2)

fams_scfea_fli <- split(fams_scfea_fli$gene,fams_scfea_fli$cluster)
scw@assays$FLUX@key <- "flux_"
scw@assays$BALANCE@key <- "balance_"

DefaultAssay(scw) <- "BALANCE"
FeaturePlot(scw,"M-58",reduction = "scphere")
Nebulosa::plot_density(scw,features = "M-58",reduction = "scphere")
VlnPlot(scw,features = "M-123",group.by = "NEW.clu")



flux_5_7 <- FindConservedMarkers(scw,ident.1 = 5,ident.2 = 7,logfc.threshold = 0,min.pct = 0,
                                 grouping.var = "GSE",assay = "FLUX")

scFEA_dict <- read.csv("scFEA/fea.module.csv",row.names = 1)
scFEA_dict$X <- rownames(scFEA_dict)
setdiff(1:171,scFEA_dict$Module_id)
scFEA_dict$X <- stringr::str_replace(scFEA_dict$X,pattern = "_",replacement = "-")


FEA <- AverageExpression(scw,assays = "scFEA",slot = "scale.data")[[1]]
for (variable in vector) {
  
}
Nebulosa::plot_density(scw,features = "lysine",reduction = "scphere")

scw@assays$scFEA@key <- "scfea_"

scFEA_genes <- read.csv("scFEA/genes.csv",header = F,row.names = 1)

scFEA_count <- scFEA_genes[,1]
scFEA_genes <- scFEA_genes[,-1]
scFEA_genes <- as.data.frame(t(scFEA_genes))
scFEA_genes <- as.list(scFEA_genes)
for (i in 1:length(scFEA_genes)) {
  scFEA_genes[[i]] <- scFEA_genes[[i]][1:scFEA_count[i]]
}

fea_genes <- as.data.frame(unique(unlist(scFEA_genes)))
names(fea_genes) <- "Gene.Symbol"

fea_genes$module <- unlist(gene_module[match(fea_genes$Gene.Symbol,gene_module$id),2])
tmp_fea <- melt(table(fea_genes$module))

SM_dict <- read.csv("scFEA/su.module.csv")
FEA_Genes_dict <- plyr::ldply(scFEA_genes, data.frame)
names(FEA_Genes_dict) <- c("Module_ID","Gene_Symbol")
FEA_Genes_dict[,3:4] <- scFEA_dict[match(FEA_Genes_dict$Module_ID,scFEA_dict$old),6:7]
FEA_Genes_dict[,5:7] <- SM_dict[match(FEA_Genes_dict$Supermodule_id,SM_dict$SM_ID),2:4]

FEA_Genes_dict[,8] <- gene_module[match(FEA_Genes_dict$Gene_Symbol,gene_module$id),2]
names(FEA_Genes_dict)[8] <- "Monocle_Module"

FAMs_type.fli <- FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$avg_log2FC>0),]
FAMs_type.fli <- FAMs_type.fli[!duplicated(FAMs_type.fli$gene),]

FEA_Genes_dict[,9] <- FAMs_type[match(FEA_Genes_dict$Gene_Symbol,FAMs_type.fli$gene),6]
names(FEA_Genes_dict)[9] <- "Cell_type"

tmp_fea_sel <- FEA_Genes_dict[!duplicated(FEA_Genes_dict$Gene_Symbol),]

tmp_fea_sel <- melt(table(tmp_fea_sel$Supermodule_id,tmp_fea_sel$Monocle_Module))
names(tmp_fea_sel) <- c("FEA_SU_module","Monocle_module","value")

DefaultAssay(scw) 
VlnPlot(scw,features = "M-34",group.by = "celltype")
VlnPlot(scw,features = "M-35",group.by = "celltype")

FeatureScatter(scw,feature1 = "M-34",feature2 = "M-35",group.by = "celltype",jitter = T)
DefaultAssay(scw) <- "BALANCE"
FeatureScatter(scw,feature1 = "Acetyl.CoA",feature2 = "Fatty.Acid",
               group.by = "celltype",jitter = T,shuffle = T)

RidgePlot(scw,features = c(paste0("M-",1:14)),group.by = "celltype")


scw$celltype
FLUX_sel <- as.data.frame(flux[,c("M_34","M_35")])
FLUX_sel <- cbind(FLUX_sel,scw$celltype)
names(FLUX_sel)[3] <- "Celltype"

FLUX_sel <- as.data.frame(flux[,c("M_34","M_35")])
FLUX_sel <- cbind(FLUX_sel,scw$celltype)
names(FLUX_sel)[3] <- "Celltype"
library()
ggplot(data = FLUX_sel,aes(x=M_34,y=M_35,color=Celltype))+
  geom_jitter(size=1, shape = 19) + geom_smooth(method="lm") +
  ggpubr::theme_pubr()
  # geom_abline(slope = 0.33, intercept = -4.69,size=2,color="darkRED") #coloc2







tmp_fea_sel <- FEA_Genes_dict[!duplicated(FEA_Genes_dict$Gene_Symbol),]
  
tmp_fea_sel <- melt(table(tmp_fea_sel$Monocle_Module,tmp_fea_sel$Cell_type))

tmp_fea_sel <- tmp_fea_sel %>% filter(value!=0) %>% arrange(desc(value))














################################################################################


# scenic

data("pbmc3k") 
sct <- subset(scw,GSE=="GSE159677")
exprMat <-as.matrix(sct@assays$RNA@data)
dim(exprMat) 
exprMat[l:4,1:4] 

names(sct@meta.data)
cellinfo <- sct@meta.data[,c(42,43)] 
table(cellinfo$celltype) 


write.csv(t(exprMat),file = "EC_159677.csv") 

