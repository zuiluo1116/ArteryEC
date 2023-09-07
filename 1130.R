#Merge
library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(harmony)
# atherosclerosis and aneurysm
# rm(list = ls())
getwd()
scRNA <- LoadH5Seurat("D:/ZuiLuo/works/aorta_human/scHuman_raw.h5seurat")
DimPlot(scint,group.by = "seurat_clusters",label = T,split.by = "GSE")
DefaultAssay(scRNA)
FeaturePlot(scRNA,features = c("TAGLN","PTPRC","CDH5","VCAM1"))

VlnPlot(scRNA,features = "CDH5") #4
int_EC <- subset(scRNA,seurat_clusters %in% c(4))
SaveH5Seurat(int_EC,"int_EC")


# load
scint <- LoadH5Seurat("int_EC.h5seurat")
# GSE179159 <- LoadH5Seurat("GSE179159_EC.h5seurat") # exclusion
GSE189179 <- LoadH5Seurat("GSE189179_EC.h5seurat")
GSE213740 <- LoadH5Seurat("GSE213740_EC.h5seurat")
GSE216860 <- LoadH5Seurat("GSE216860_EC.h5seurat")

mdata_scint <- scint@meta.data
mdata_scint <- select(mdata_scint,c("orig.ident","GSE","state"))
scint <- CreateSeuratObject(scint@assays$RNA@counts, meta.data = mdata_scint)

# mdata_GSE179159 <- GSE179159@meta.data
# mdata_GSE179159 <- select(mdata_GSE179159,c("orig.ident","GSE"))
# levels(factor(mdata_GSE179159$orig.ident))
# DimPlot(GSE179159,label = T,split.by = "orig.ident")
# which(mdata_GSE179159$orig.ident %in% c(Coronary) )
# rm(list = c("GSE179159","mdata_GSE179159"))

mdata_GSE189179 <- GSE189179@meta.data
mdata_GSE189179 <- select(mdata_GSE189179,c("orig.ident","GSE"))
levels(factor(mdata_GSE189179$orig.ident))
mdata_GSE189179[which(mdata_GSE189179$orig.ident %in% c("MT1")),3] <- "4kPa"
mdata_GSE189179[which(mdata_GSE189179$orig.ident %in% c("MT2")),3] <- "4kPa+TGFB"
mdata_GSE189179[which(mdata_GSE189179$orig.ident %in% c("MT3")),3] <- "TC"
mdata_GSE189179[which(mdata_GSE189179$orig.ident %in% c("MT4")),3] <- "TC+TGFB"
names(mdata_GSE189179)[3] <- "state"
GSE189179@meta.data <- mdata_GSE189179
GSE189179 <- CreateSeuratObject(GSE189179@assays$RNA@counts, meta.data = mdata_GSE189179)


mdata_GSE213740 <- GSE213740@meta.data
mdata_GSE213740 <- select(mdata_GSE213740,c("orig.ident","GSE"))
levels(factor(mdata_GSE213740$orig.ident))
mdata_GSE213740[which(mdata_GSE213740$orig.ident %in% c("AD_2","AD_3","AD_4","AD_5","AD_6")),3] <- "AD_lesion"
mdata_GSE213740[which(mdata_GSE213740$orig.ident %in% c("Normal_1","Normal_2","Normal_3")),3] <- "AD_control"
names(mdata_GSE213740)[3] <- "state"
GSE213740@meta.data <- mdata_GSE213740
GSE213740 <- CreateSeuratObject(GSE213740@assays$RNA@counts, meta.data = mdata_GSE213740)

mdata_GSE216860 <- GSE216860@meta.data
mdata_GSE216860 <- select(mdata_GSE216860,c("orig.ident","GSE"))
levels(factor(mdata_GSE216860$orig.ident))
mdata_GSE216860[,3] <- "age_aorta"
names(mdata_GSE216860)[3] <- "state"
GSE216860@meta.data <- mdata_GSE216860
DimPlot(GSE216860,group.by = "state")
GSE216860 <- CreateSeuratObject(GSE216860@assays$RNA@counts, meta.data = mdata_GSE216860)

scRNA <- merge(x=scint,y=c(GSE189179,GSE213740,GSE216860))
mdata <- scRNA@meta.data
levels(factor(mdata$state))
mdata[which(mdata$state %in% c("4kPa","4kPa+TGFB","TC","TC+TGFB")),5] <- paste0("MT_",mdata[which(mdata$state %in% c("4kPa","4kPa+TGFB","TC","TC+TGFB")),5] )
mdata[which(mdata$state %in% c("TAA")),5] <- "TAA_lesion"
mdata[which(mdata$state %in% c("TAC")),5] <- "TAA_control"
mdata[which(mdata$state %in% c("TAC")),5] <- "TAA_control"
mdata[which(mdata$state %in% c("AAA")),5] <- "AAA_lesion"
mdata[which(mdata$state %in% c("AAC")),5] <- "AAA_control"
mdata[which(mdata$state %in% c("CAAC")),5] <- "CAA_lesion"
mdata[which(mdata$state %in% c("CAPA")),5] <- "CAA_control"
mdata[which(mdata$state %in% c("CAAS")),5] <- "CAs_lesion"

mdata[which(mdata$orig.ident %in% c("AD_2","AD_3","AD_4","AD_5","AD_6")),5] <- "AD_lesion"
mdata[which(mdata$orig.ident %in% c("Normal_1","Normal_2","Normal_3")),5] <- "AD_control"
mdata[which(mdata$GSE %in% c("GSE216860")),5] <- "Aged_aorta"

table(mdata$GSE,mdata$state)
levels(factor(mdata$state))
mdata$state <- factor(mdata$state,levels = c("TAA_control","TAA_lesion","CAs_lesion","CAA_control","CAA_lesion",
                                             "AAA_control","AAA_lesion","MT_TC","MT_TC+TGFB","MT_4kPa","MT_4kPa+TGFB",
                                             "AD_control","AD_lesion","Aged_aorta"))


levels(factor(mdata$order))
table(names$order,names$state)

mdata$order <- paste0(mdata$GSE,"_",mdata$orig.ident)
mdata$order <- factor(mdata$order,
                           levels = c("GSE155468_Con4","GSE155468_Con6","GSE155468_Con9",
                                      "GSE155468_TAA1","GSE155468_TAA2","GSE155468_TAA3","GSE155468_TAA4","GSE155468_TAA5","GSE155468_TAA6","GSE155468_TAA7","GSE155468_TAA8",
                                      "GSE155514_RPE004","GSE155514_RPE005","GSE155514_RPE006",
                                      "GSE159677_Patient1.PA","GSE159677_Patient2.PA","GSE159677_Patient3.PA",
                                      "GSE159677_Patient1.CA","GSE159677_Patient2.CA","GSE159677_Patient3.CA",
                                      "GSE166676_0703","GSE166676_1127","GSE166676_0226","GSE166676_0524","GSE166676_05241","GSE166676_0617",
                                      "GSE189179_MT3","GSE189179_MT4","GSE189179_MT1","GSE189179_MT2",
                                      "GSE213740_Normal_1","GSE213740_Normal_2","GSE213740_Normal_3",
                                      "GSE213740_AD_2","GSE213740_AD_3","GSE213740_AD_4","GSE213740_AD_5","GSE213740_AD_6",
                                      "GSE216860_Normal_1","GSE216860_Normal_2","GSE216860_Normal_3","GSE216860_Normal_4","GSE216860_Normal_5","GSE216860_Normal_6"))

names <- mdata
names$old.name <- rownames(names)
names$barcode <- stringr::str_extract(names$old.name,pattern = "^[AGCT]{16}")

order <- 1:44
names(order) <- levels(factor(names$order))
names$no <- order[match(names$order,names(order))]
table(names$GSE,names$no)

names$new.name <- paste0(names$barcode,"_",names$no)
scbak <- scRNA

mdata$new.name <- names[match(rownames(mdata),names$old.name),10]
rownames(mdata) <- mdata$new.name
mdata$new.name <- NULL
mdata$state <- NULL
mdata$order <- NULL

scRNA <- RenameCells(scRNA,new.names = rownames(mdata))
scRNA@meta.data <- mdata


#########################################################################################
#filter out
scRNA <- subset(scRNA.sct.int,GSE!="GSE189179")
Idents(scRNA.sct.int) <- "integrated_snn_res.0.4"
sets <- 0:17 #8 9 16 12 14 16
sets <- sets[-c(9,10,17,13,15,17)]
DefaultAssay(scRNA.sct.int)
VlnPlot(scRNA.sct.int,features = c("CDH5","TAGLN","PTPRC","DCN"),stack = T,flip = T)
scRNA <- subset(scRNA.sct.int,integrated_snn_res.0.4 %in% sets)
scRNA <- subset(scRNA,GSE!="GSE189179")
DimPlot(scRNA,label = T,split.by = "GSE",ncol = 3)
VlnPlot(scRNA,features = c("CDH5","TAGLN","PTPRC","DCN"),stack = T,flip = T)

scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = scRNA@meta.data)
#END of MERGE
############################################################################################

scRNA <- SCTransform(scRNA)

scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 

ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
scRNA <- scRNA %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = c(seq(0,1,0.1)))
Idents(scRNA) 
DimPlot(scRNA,group.by = "orig.ident", split.by = "GSE",shuffle = T,ncol = 3)
DimPlot(scRNA,label = T,split.by = "GSE")

sets2 <- 0:20
sets2 <- sets2[-c(20,21)]
scRNA <- subset(scRNA, SCT_snn_res.1 %in% sets2)
DimPlot(scRNA,label = T)

DefaultAssay(scRNA) <- "RNA"

VlnPlot(scRNA,features = c("CDH5","TAGLN","PTPRC","DCN"),stack = T,flip = T)


#########手动筛选
expr <- scRNA@assays$RNA@data
gene_expression <- expr %>% .[c("DCN","PTPRC","TAGLN"),] %>% as.data.frame()
gene_expression <- as.data.frame(t(gene_expression))
gene_expression$cell <- rownames(gene_expression)

#选择阳性细胞
gene_expression_sel <- gene_expression[which(gene_expression$DCN>5|gene_expression$PTPRC>1|gene_expression$TAGLN>50),]
exclu <- scRNA[,rownames(gene_expression_sel)]
cells <- colnames(scRNA)


out <- rownames(gene_expression_sel)
select.cells <- CellSelector(plot = plot)
out <- union(out,select.cells)

exclu <- setdiff(cells,out)
scRNA.fil <- scRNA[,exclu]
# FeaturePlot(scRNA,features = c())
plot <- DimPlot(scRNA.fil,label = T)

SaveH5Seurat(scRNA.fil,"scRNA_filtered_harmony")








FeaturePlot(scRNA.fil,features = c("CDH5","TAGLN","PTPRC","DCN"),ncol = 2)


FeaturePlot(scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt",
                               "percent.rp","percent.hb"),ncol = 3)
table(scRNA@meta.data$GSE)
SaveH5Seurat(scRNA,"scRNA_harmony")
##########################################################################################

scRNAlist <- SplitObject(scRNA.fil,split.by = "GSE")
listk <- ls()
listk <- listk[-15]
rm(list = listk)
scRNAlist$GSE155468 <- merge(scRNAlist$GSE166676,scRNAlist$GSE155468)
scRNAlist$GSE166676 <- NULL
scRNAlist <- parallel::mclapply(scRNAlist, FUN=function(x) SCTransform(x), mc.cores = 1) 

scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000)

scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features) 

# plan("multisession", workers = 10)

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
scRNA.sct.int <- IntegrateData(scRNA.anchors, normalization.method="SCT") #
DefaultAssay(scRNA.sct.int) <- "integrated"
scRNA.sct.int <- ScaleData(scRNA.sct.int, verbose = FALSE)
scRNA.sct.int <- RunPCA(scRNA.sct.int, npcs = 30, verbose = FALSE)
scRNA.sct.int <- RunUMAP(scRNA.sct.int, reduction = "pca", dims = 1:30)
scRNA.sct.int <- FindNeighbors(scRNA.sct.int, reduction = "pca", dims = 1:30)
scRNA.sct.int <- FindClusters(scRNA.sct.int, resolution = c(seq(0,1,0.1)))
DimPlot(scRNA.sct.int,label = T)
scRNA.sct.int$order <- paste0(scRNA.sct.int$GSE,"_",scRNA.sct.int$orig.ident)
DimPlot(scRNA.sct.int,label = T,split.by = "order",ncol = 7)

DefaultAssay(scRNA.sct.int) <- "RNA"

FeaturePlot(scRNA.sct.int,features = "DCN")

################################################################################
scRNA.sct.int.fil <- subset(scRNA.sct.int,GSE!= "GSE166676")

DimPlot(scRNA.sct.int.fil,label = T,split.by = "order",ncol = 7)
FeaturePlot(scRNA.sct.int.fil,features = "DKK2")

Idents(scRNA.sct.int.fil) <- "integrated_snn_res.0.4"
table(scRNA.sct.int.fil@meta.data$orig.ident,scRNA.sct.int.fil@meta.data$SCT_snn_res.0.4)



##################################### try2
scRNA <- CreateSeuratObject(scRNA.sct.int.fil@assays$RNA@counts, meta.data = scRNA.sct.int.fil@meta.data)

scRNA <- SCTransform(scRNA)

scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 

ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
scRNA <- scRNA %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = c(seq(0,1,0.1)))
Idents(scRNA) 
DimPlot(scRNA,group.by = "orig.ident", split.by = "GSE",shuffle = T,ncol = 3)
DimPlot(scRNA,label = T,split.by = "order",ncol = 6)

cosg2 <- COSG::cosg(scRNA,groups = "all",assay = "RNA",slot = "count",n_genes_user = 50)


















cosg <- COSG::cosg(scRNA.sct.int.fil,groups = "all",assay = "RNA",slot = "count",n_genes_user = 50)

DimPlot(scRNA.sct.int.fil,label = T)

library(plotly)
# Re-run UMAPs that you have accurate calculations for all UMAP(s)
sc.sample <- scRNA.sct.int.fil
DefaultAssay(sc.sample) <- "integrated"


sc.sample <- RunUMAP(sc.sample,
                     dims = 1:30,
                     n.components = 3,reduction = "pca")
head(sc.sample[["umap"]]@cell.embeddings)

umap_1 <- sc.sample[["umap"]]@cell.embeddings[,1]
umap_2 <- sc.sample[["umap"]]@cell.embeddings[,2]
umap_3 <- sc.sample[["umap"]]@cell.embeddings[,3]
# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = sc.sample, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "integrated_snn_res.0.4"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data),plot.data$integrated_snn_res.0.4)

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", 
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")
plot_ly(data = plot.data, size = 0.4,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~integrated_snn_res.0.4, 
        colors = cb_palette,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


###################################################################


expr <- scRNA@assays$RNA@data


MT <- grep(x = scRNA@assays[["RNA"]]@counts@Dimnames[[1]],pattern = "^MT-")
MT.genes <- scRNA@assays[["RNA"]]@counts@Dimnames[[1]][MT]

RP <- grep(x = scRNA@assays[["RNA"]]@counts@Dimnames[[1]],pattern = "^RP[SL]")
RP.genes <- scRNA@assays[["RNA"]]@counts@Dimnames[[1]][RP]

HB <- grep(x = scRNA@assays[["RNA"]]@counts@Dimnames[[1]],pattern = "^HB")
HB.genes <- scRNA@assays[["RNA"]]@counts@Dimnames[[1]][HB]
HB.genes <- HB.genes[4:13]

outgenes <- union(MT.genes,RP.genes)
outgenes <- union(outgenes,HB.genes)

remain_genes <- setdiff(expr@Dimnames[[1]],outgenes)

exprw <- expr[remain_genes,]

scw <- CreateSeuratObject(exprw, meta.data = scRNA.sct.int.fil@meta.data)

scw <- SCTransform(scw)

scw <- RunPCA(scw, npcs=50, verbose=FALSE)
sce <- RunPCA(sce, npcs = 30, verbose = FALSE)
scw <- RunHarmony(scw, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 

ElbowPlot(scw, ndims = 50)
pc.num=1:50
scw <- RunTSNE(scw, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
scw <- scw %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = c(seq(0,1,0.1)))

scw <- FindClusters(scw,resolution = 0.4)
Idents(scw) 
DimPlot(scw,group.by = "orig.ident", split.by = "GSE",shuffle = T,ncol = 3)
DimPlot(scw,label = T,split.by = "order",ncol = 6)
DimPlot(scw,label = T)
cosgw <- COSG::cosg(scw,groups = "all",assay = "RNA",slot = "count",n_genes_user = 100)
cosgwk <- as.data.frame(cosgw[1])
genelist <- character()
cluster <- character()
for (i in 1:ncol(cosgwk)) {
  genelist <- c(genelist,cosgwk[,i])
  cluster <- c(cluster,rep((i-1),100))
}
genelist <- as.data.frame(genelist)
genelist$clusters <- cluster
genelist$clusters <- factor(genelist$clusters,levels = c(0:13))
# genelist <- split(genelist$genelist,genelist$clusters)

a <- genelist
b <- paste0(i)
c <- "ToppCell"

Topp.gene <- function(a,b,c)
{
  tmp <- a[which(a$cluster== b ),1]
  url = "https://toppgene.cchmc.org/API/lookup"
  encode = "json"
  payload = list(Symbols = tmp)
  response = POST(url = url, body = payload, encode = encode)
  json = httr::content(response, "text")
  data = fromJSON(rawToChar(response$content))
  results = fromJSON(json)
  IDs = results[["Genes"]][["Entrez"]]
  url = "https://toppgene.cchmc.org/API/enrich"
  encode = "json"
  payload = list(Genes = IDs)
  response = POST(url = url, body = payload, encode = encode)
  json = httr::content(response, "text")
  results = fromJSON(json)
  res <- results[["Annotations"]]
  res <- res[which(res$Category %in% c),]
  return(res)
}
toppgene_scw <- list()

cis <- levels(factor(results[["Annotations"]][["Category"]]))
toppgene_clu2 <- Topp.gene(genelist,2,cis)

for (i in 0:13) {
  # i=2
  x <- Topp.gene(genelist,paste0(i),c("ToppCell"))
  toppgene_scw[[i+1]] <- x
}

names(toppgene_scw) <- 0:13

###################################################################


library(monocle3)
data<-GetAssayData(scw,assay ='RNA',slot ='counts')
cell_metadata <-scw@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(data))
rownames(gene_annotation)<-rownames(data)
cds <-new_cell_data_set(data,cell_metadata =cell_metadata,gene_metadata =gene_annotation)
# 主成分
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
# 降维 默认UMAP
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_cells(cds,show_trajectory_graph = F)

# 替换UMAP坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scw, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
# 识别分群（基于UMAP二维）
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds@colData@listData[["mono.clu"]] <- cds@clusters@listData[["UMAP"]][["clusters"]]

plot_cells(cds, reduction_method="UMAP", color_cells_by="mono.clu")

df.seu <- select(scw@meta.data,"seurat_clusters")
df.seu$names <- rownames(df.seu)
df.mon <- data.frame("clusters" =  cds@clusters@listData[["UMAP"]][["clusters"]],
                     "names" = cds@colData@rownames)
df.seu$mono.clu <- df.mon[match(df.seu$names,df.mon$names),1]
scw@meta.data$mono.clu <- df.seu$mono.clu
DimPlot(scw,group.by = "mono.clu",split.by = "order",ncol = 6)

cds <- learn_graph(cds)
# 给定起点绘制轨迹
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "mono.clu", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=T)
# debug monocle3 line 93
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# Matrix::rBind() to rbind()
Track_genes <- graph_test(cds, neighbor_graph="principal_graph")

write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)

Track_genes.fli <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes.fli, "Trajectory_genes.fli.csv", row.names = F)

Track_genes_sig10 <- Track_genes.fli %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
Track_genes_sig100 <- Track_genes.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

FeaturePlot(scw,features = "SCG3")
# Track_genes_sig100 <- Track_genes.fli %>% top_n(n=100, morans_I) %>%
#   pull(gene_short_name) %>% as.character()

plot_genes_in_pseudotime(cds[Track_genes_sig10,], color_cells_by="mono.clu", 
                         min_expr=0.5, ncol = 5)
ggsave(filename = "2_8_1.pdf",height = 6,width = 8)
VlnPlot(scRNA,features = Track_genes_sig10,group.by = "mono.clu",stack = T,flip = T)+NoLegend()
ggsave(filename = "2_8_2.pdf",height = 6,width = 8)

#回传
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(scw@meta.data)]
scw$pseudotime <- pseudotime

# 识别模块基因
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
set.seed(1116)
gene_module <- find_gene_modules(cds[genelist,], resolution=5e-3, cores = 1)

# gene_module <- read.csv("./figures/2_9_5.Genes_Module.csv")
write.csv(gene_module, "Genes_Module.csv", row.names = F)

cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=cds@clusters@listData[["UMAP"]][["clusters"]])
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pseudo.order <- scw@meta.data %>%
  dplyr::select(c("pseudotime","mono.clu")) %>%
  group_by(mono.clu) %>%
  summarise(pseudo.mean = mean(pseudotime), pseudo.median = median(pseudotime), n = n()) %>%
  arrange(pseudo.mean)
scRNA@meta.data$mono.clu <- factor(scRNA@meta.data$mono.clu,levels = pseudo.order$mono.clu)


gene_module.list <- gene_module
gene_module.list$name <- paste("Module",gene_module.list$module,sep = " ")
write.csv(gene_module.list,"gene_module.list.csv")



cds_subset <- choose_cells(cds)
plot_cells(cds_subset, color_cells_by = "mono.clu", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph")

cds_subset_LSS <- cds_subset
LSS_trajectory_gene <- subset_pr_test_res
LSS_trajectory_gene <- LSS_trajectory_gene[order(LSS_trajectory_gene$morans_test_statistic,decreasing = T),]

cds_subset_Unstable <- cds_subset
US_trajectory_gene <- subset_pr_test_res
US_trajectory_gene <- US_trajectory_gene[order(US_trajectory_gene$morans_test_statistic,decreasing = T),]

intersect(US_trajectory_gene$gene_short_name[1:50],LSS_trajectory_gene$gene_short_name[1:50])

head(agg_mat)[1:4,1:4]
agg_mat <- as.data.frame(agg_mat[,pseudo.order$mono.clu])
agg_sel <- agg_mat[,c(2,4,8)]

pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",
                   cellwidth = 10, cellheight = 4,border_color = "white",
                   cluster_cols = F)
###################################################################
# 
# library(ReactomeGSA)
# gsva_result <- analyse_sc_clusters(scw, verbose = TRUE)
# pathway_expression <- pathways(gsva_result)
# colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
# pathway_expression[1:3,]
# gsva_diff <- abs(pathway_expression$X4-pathway_expression$X1)
# names(gsva_diff) <- pathway_expression$Name
# gsva_diff[which(gsva_diff==max(gsva_diff))]
# 
# max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
#   values <- as.numeric(row[2:length(row)])
#   return(data.frame(name = row[1], min = min(values), max = max(values)))
# }))
# 
# max_difference$diff <- max_difference$max - max_difference$min
# # sort based on the difference
# max_difference <- max_difference[order(max_difference$diff, decreasing = T), ] 
# head(max_difference)
# dir.create("./ReactomeGSA")
# 
# Idents(scw) <- "mono.clu"
# for (i in 1:20) {
#   p1 <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[i]) 
#   p1$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> df
#   # ?element_blank
#   df$cluster_id <- factor(stringr::str_split_fixed(df$cluster_id,"X",n = 2)[,2],
#                           levels = pseudo.order$mono.clu)
#   plot.tmp <- df %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+ 
#     geom_bar(stat='identity') + theme_bw()+ 
#     theme(axis.text.x = element_text(angle =45,hjust = .9,size = 10,vjust = 0.9))+ 
#     ggtitle(max_difference$name[i]) + theme(legend.position="none")+
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank(), 
#           axis.line = element_line(colour = "black"),
#           plot.title = element_text(size=15,hjust = 0.5),axis.title = element_text(size=6))
#   ggsave(plot = plot.tmp,filename = paste0("./ReactomeGSA/",i,".",max_difference$name[i],".pdf"),width = 8,height = 6)
# }
# DimPlot(scw,group.by = "mono.clu",label = T,label.size = 18)
# 


expr_monoclu <- AverageExpression(scw,slot = "count",group.by = "mono.clu")[[1]]
#### KEGG
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(homologene)
ls("package:KEGG.db")
# [1] "KEGG"             "KEGG_dbconn"      "KEGG_dbfile"      "KEGG_dbInfo"      "KEGG_dbschema"   
# [6] "KEGGENZYMEID2GO"  "KEGGEXTID2PATHID" "KEGGGO2ENZYMEID"  "KEGGMAPCOUNTS"    "KEGGPATHID2EXTID"
# [11] "KEGGPATHID2NAME"  "KEGGPATHNAME2ID" 
kegg.ent <- as.list(KEGGPATHID2EXTID)
names(kegg.ent) <- stringr::str_replace_all(names(kegg.ent),pattern = "mmu",replacement = "hsa")

# 匹配ID 和 name
tmp <- data.frame("Pathway"=unlist(as.list(KEGGPATHID2NAME)))
tmp$ID <- paste0("hsa",rownames(tmp))
kegg.dict <- tmp

tmp4 <- as.data.frame(unlist(kegg.ent))
tmp4 <- tmp4[!duplicated(unlist(kegg.ent)),]
tmp4 <- bitr(tmp4,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
kegg.genelist <- tmp4
dict_tmp <- homologene(kegg.genelist$SYMBOL,inTax = 10090, outTax = 9606)
kegg.genelist$Human <- dict_tmp[match(kegg.genelist$SYMBOL,dict_tmp$`10090`),2]
kegg.genesets <- lapply(kegg.ent,FUN = match.gen, ref = kegg.genelist, que_col = 1, ref_col=3)
kegg.genesets.ID <- kegg.genesets
names(kegg.genesets) <- tmp[match(names(kegg.genesets),tmp$ID),1]

kegg.manul <- GSVA::gsva(expr = expr_monoclu,gset.idx.list = kegg.genesets)


# 
# que %in% ref[,que_col]
# kegg.ent[[1]][1] %in% tmp4[,1]
# que <- tmp4[match(kegg.ent[[1]][1],tmp4[,1]),2]
# kegg.genesets <- gmt2list("./KEGG/kegg.gmt")

# 重命名list
kegg.genesets.ID <- kegg.genesets
names(kegg.genesets) <- tmp[match(names(kegg.genesets),tmp$ID),1]

kegg.manul <- GSVA::gsva(expr = expr_monoclu,gset.idx.list = kegg.genesets)
kegg.manul.ext <- cbind(rownames(kegg.manul),kegg.manul)

max_diff_kegg <- do.call(rbind, apply(kegg.manul.ext, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values),
                      clu_4 = row[5], clu_1 = row[2], clu_2 = row[3]))
  }))
max_diff_kegg$diff <- max_diff_kegg$max - max_diff_kegg$min
max_diff_kegg$diff_4v1 <- abs(as.numeric(max_diff_kegg$clu_4) - as.numeric(max_diff_kegg$clu_1))
# sort based on the difference
max_diff_kegg <- max_diff_kegg[order(max_diff_kegg$diff_4v1, decreasing = T), ]








FeaturePlot(scw,features = "IGFBP7")
DimPlot(scw,label = T)


#####################################################################
#借助Scissor评估起点
# devtools::install_github('sunduanchen/Scissor')

library(Scissor)
# GSE120521
bulk_dataset <- as.matrix(bulk)
rownames(bulk_dataset) <- bulk_dict$`bulk_ori[, 1]`
colnames(bulk_dataset) <- c(paste0("Stable_",1:4),paste0("Unstable_",1:4))
phenotype <- c(rep(0,4),rep(1,4))
names(phenotype) <- c(rep("stable",4),rep("unstable",4))
class(phenotype)
getwd()
DefaultAssay(scw) <- "RNA"
##Execute Scissor to select the informative cells
#use Scissor to select the phenotype-associated cell subpopulations
infos1 <- Scissor_SCT(bulk_dataset, scw, phenotype, alpha = 0.05, tag = tag,
                  family = "binomial", #二分类
                  Save_file = 'Scissor_stable.RData')
infos2 <- Scissor_SCT(bulk_dataset, scw, phenotype, alpha = NULL, cutoff = 0.03, 
                  family = "binomial", Load_file = 'Scissor_stable.RData')



# common <- intersect(rownames(bulk_dataset), rownames(scw))
# alpha = 0.05
# sc_exprs <- as.matrix(scw@assays$RNA@data)
# network <- as.matrix(scw@graphs$SCT_snn)
#####################################################################
Scissor_SCT <- function (bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                         cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                         Save_file = "Scissor_inputs.RData", Load_file = NULL) 
{
  library(Seurat)
  library(Matrix)
  library(preprocessCore)
  if (is.null(Load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    if (class(sc_dataset) == "Seurat") {
      sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
      network <- as.matrix(sc_dataset@graphs$SCT_snn)
    }
    else {
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", 
                                         verbose = F)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), 
                           verbose = F)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, 
                                  verbose = F)
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    diag(network) <- 0
    network[which(network != 0)] <- 1
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, 
    ])
    dataset1 <- normalize.quantiles(dataset0)
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 
                                     1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
                      z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      tag <- c(1,0)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", 
                     paste(tmp[1:(length(z) - 1)], collapse = ", "), 
                     ", and ", tmp[length(z)], " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      }
      else {
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    save(X, Y, network, Expression_bulk, Expression_cell, 
         file = Save_file)
  }
  else {
    load(Load_file)
  }
  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
               0.6, 0.7, 0.8, 0.9)
  }
  for (i in 1:length(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, nlambda = 100, 
                  nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    }
    else {
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
                  length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", 
                  formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
    cat("\n")
  }
  print("|**************************************************|")
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
                          family = family), Coefs = Coefs, Scissor_pos = Cell1, 
              Scissor_neg = Cell2))
}
#####################################################################

# 
# diag(network) <- 0
# network[which(network != 0)] <- 1
# x1 <- bulk_dataset[common, ]
# x2 <- sc_exprs[common, ]
# rm(sc_exprs)
# dataset0 <- cbind(x1, x2)
# 
# dataset1 <- normalize.quantiles(dataset0)
# rownames(dataset1) <- rownames(dataset0)
# colnames(dataset1) <- colnames(dataset0)
# rm(dataset0)
# 
# Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
# Expression_cell <- dataset1[, (ncol(bulk_dataset) + 
#                                  1):ncol(dataset1)]
# X <- cor(Expression_bulk, Expression_cell)
# quality_check <- quantile(X)
# if (quality_check[3] < 0.01) {
#   warning("The median correlation between the single-cell and bulk samples is relatively low.")
# }
# Y <- as.numeric(phenotype$Stable)
# z <- table(Y)
# if (length(z) != length(tag)) {
#   stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
# }
# else {
#   print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
#                 z[1], tag[1], z[2], tag[2]))
#   print("Perform logistic regression on the given phenotypes:")
# }
# 
# save(X, Y, network, Expression_bulk, Expression_cell, 
#      file = "Scissor_stable.RData")
# family = "binomial"
# cutoff = 0.2
# for (i in 1:length(alpha)) {
#   set.seed(123)
#   fit0 <- APML1(X, Y, family = family, penalty = "Net", 
#                 alpha = alpha[i], Omega = network, nlambda = 100, 
#                 nfolds = min(10, nrow(X)))
#   fit1 <- APML1(X, Y, family = family, penalty = "Net", 
#                 alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
#   if (family == "binomial") {
#     Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
#   }
#   else {
#     Coefs <- as.numeric(fit1$Beta)
#   }
#   Cell1 <- colnames(X)[which(Coefs > 0)]
#   Cell2 <- colnames(X)[which(Coefs < 0)]
#   percentage <- (length(Cell1) + length(Cell2))/ncol(X)
#   print(sprintf("alpha = %s", alpha[i]))
#   print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
#                 length(Cell1), length(Cell2)))
#   print(sprintf("The percentage of selected cell is: %s%%", 
#                 formatC(percentage * 100, format = "f", digits = 3)))
#   if (percentage < cutoff) {
#     break
#   }
#   cat("\n")
# }
# print("|**************************************************|")
# infos1 <- list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
#                         family = family), Coefs = Coefs, Scissor_pos = Cell1, 
#             Scissor_neg = Cell2)

#visualize the Scissor selected cells by adding a new annotation in the Seurat object
Scissor_select <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_select) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_select[infos2$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_select[infos2$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_select, col.name = "scissor")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'umap', split.by = 'scissor', pt.size = 1.2, order = c(2,1))#可视化
DimPlot(scw,group.by = "mono.clu",label = T)+NoLegend()
table(scw$scissor,scw$mono.clu)
table(scw$scissor_LSS,scw$mono.clu)

##############################################################

#GSE199709 
read <- list.files("./GSE199709_RAW")
read <- read[-10]
kk <- list()
for (i in 1:length(read)) {
  kk[[i]] <- read.table(paste0("./GSE199709_RAW/",read[i]),header = T)
}

kk <- as.data.frame(kk)
kk <- kk[!duplicated(kk$gene_id),]
rownames(kk) <- kk$gene_id
kk <- kk[,seq(2,18,2)]

GSE199709 <- as.matrix(kk)[,4:9]
colnames(GSE199709) <- c(paste0("LSS_",1:3),paste0("Control_",1:3))
phenotype2 <- c(rep(0,3),rep(1,3))
names(phenotype2) <- c(rep("Control",3),rep("LSS",3))

infos3 <- Scissor_SCT(GSE199709, scw, phenotype2, alpha = 0.05, tag = tag,
                      family = "binomial", #二分类
                      Save_file = 'Scissor_LSS.RData')
infos4 <- Scissor_SCT(GSE199709, scw, phenotype2, alpha = NULL, cutoff = 0.03, 
                      family = "binomial", Load_file = 'Scissor_LSS.RData')

Scissor_select2 <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_select2) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_select2[infos4$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_select2[infos4$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_select2, col.name = "scissor_LSS")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'umap', split.by = 'scissor_LSS', pt.size = 1.2, order = c(2,1))#可视化
FeaturePlot(scw,features = "SREBF2")



SaveH5Seurat(scw,"20221201")
save(list = c("cdsk_3d","cds","cdsk","Track_genes","Track_genes_scphere","gene_module"),
     file = "mono_object")

GSE120521 <- bulk_dataset
save(list = c("GSE120521","GSE199709","GSE210522","GSE211662"),file = "bulk.RData")
save(list = c("scRNA","sce","scRNA.sct.int","scRNA.sct.int.fil"),file = "scs.RData")

#########################################
# GSE210522
library(Scissor)
GSE210522 <- as.matrix(read.table("GSE210522/raw.txt",header = T,row.names = 1))
samples <- as.data.frame(names(GSE210522))

GSE210522_sel <- as.matrix(na.omit(GSE210522[,17:28]))
dim(GSE210522_sel)

GSE210522_IL1b <- GSE210522_sel[,c(10:12,4:6)]
phenotype4 <- c(rep(0,3),rep(1,3))
names(phenotype4) <- c(rep("VEGF",3),rep("IL1b",3))
tag2 <- c("VEGF","IL1b")

infos7 <- Scissor_SCT(GSE210522_IL1b, scw, phenotype4, alpha = 0.05, tag = tag2,
                      family = "binomial", #二分类
                      Save_file = 'Scissor/GSE210522_IL1b.RData')

Scissor_Il1b <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_Il1b) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_Il1b[infos7$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_Il1b[infos7$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_Il1b, col.name = "scissor_Il1b")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'scphere', group.by = 'scissor_Il1b', pt.size = 1.2, order = c(2,1))#可视化


GSE210522_TNFa <- GSE210522_sel[,c(1:3,7:9)]
phenotype5 <- c(rep(0,3),rep(1,3))
names(phenotype5) <- c(rep("Ctrl",3),rep("TNFa",3))
tag3 <- c("Ctrl","TNFa")

infos8 <- Scissor_SCT(GSE210522_TNFa, scw, phenotype5, alpha = 0.05, tag = tag3,
                      family = "binomial", #二分类
                      Save_file = 'Scissor/GSE210522_TNFa.RData')
Scissor_TNFa <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_TNFa) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_TNFa[infos8$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_TNFa[infos8$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_TNFa, col.name = "Scissor_TNFa")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'scphere', group.by = 'Scissor_TNFa', pt.size = 1.2, order = c(2,1))#可视化

GSE210522_VEGF<- GSE210522_sel[,c(1:3,10:12)]
phenotype6 <- c(rep(0,3),rep(1,3))
names(phenotype5) <- c(rep("Ctrl",3),rep("VEGF",3))
tag4 <- c("Ctrl","VEGF")

infos9 <- Scissor_SCT(GSE210522_VEGF, scw, phenotype6, alpha = 0.05, tag = tag4,
                      family = "binomial", #二分类
                      Save_file = 'Scissor/GSE210522_VEGF.RData')
Scissor_VEGF <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_VEGF) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_VEGF[infos8$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_VEGF[infos8$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_VEGF, col.name = "Scissor_VEGF")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'scphere', group.by = 'Scissor_VEGF', pt.size = 1.2, order = c(2,1))#可视化

GSE210522_Il1b_VEGF<- GSE210522_sel[,c(4:6,10:12)]
phenotype7 <- c(rep(0,3),rep(1,3))
names(phenotype7) <- c(rep("Il1b",3),rep("VEGF",3))
tag5 <- c("Il1b","VEGF")

infos9 <- Scissor_SCT(GSE210522_Il1b_VEGF, scw, phenotype7, alpha = 0.05, tag = tag5,
                      family = "binomial", #二分类
                      Save_file = 'Scissor/GSE210522_Il1b_VEGF.RData')
Scissor_Il1b_VEGF <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_Il1b_VEGF) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_Il1b_VEGF[infos9$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_Il1b_VEGF[infos9$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_Il1b_VEGF, col.name = "Scissor_Il1b_VEGF")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'scphere', group.by = 'Scissor_Il1b_VEGF', pt.size = 1.2, order = c(2,1))#可视化


#########################################
# GSE211162 NOT found
GSE211662 <- read.table("GSE211662/GSE211662.csv",sep = ";",header = T,row.names = 1)

GSE211662 <- as.matrix(GSE211662)[,1:6]
phenotype3 <- c(rep(0,3),rep(1,3))
names(phenotype3) <- c(rep("LSS",3),rep("HSS",3))
gc()
infos5 <- Scissor_SCT(GSE211662, scw, phenotype3, alpha = 0.05, tag = tag,
                      family = "binomial", #二分类
                      Save_file = 'Scissor_LH.RData')
infos6 <- Scissor_SCT(GSE211662, scw, phenotype3, alpha = NULL, cutoff = 0.03, 
                      family = "binomial", Load_file = 'Scissor_LH.RData')

Scissor_select3 <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_select3) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_select3[infos5$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_select3[infos5$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_select3, col.name = "scissor_LH")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'umap', split.by = 'scissor_LH', pt.size = 1.2, order = c(2,1))#可视化
FeaturePlot(scw,features = "SREBF2")


DimPlot(scw,group.by = "scissor_LSS",label = T,cols = c("Grey","RED","Blue"),pt.size = 1.2)
DimPlot(scw,group.by = "scissor",label = T,cols = c("Grey","RED","Blue"),pt.size = 1.2)



DimPlot(scw,split.by = "GSE",ncol = 3,group.by = "scissor_LSS",
        cols = c("Grey","RED","Blue"),pt.size = 1.2)

GSEs <- levels(factor(scw$GSE))

# GSE199709 LSS split GSE
LSS_lists_A <- list()
LSS_lists_B <- list()
LSS_lists_pos <- character()
LSS_lists_neg <- character()
total <- scw@graphs[["SCT_snn"]]

for (i in 1:length(GSEs)) {
  tmp.sc <- subset(scw,GSE==GSEs[i])
  cellnames <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
  matrixsel <-total[cellnames,cellnames] 
  tmp.sc@graphs[["SCT_snn"]] <- matrixsel
  infosx <- Scissor_SCT(GSE199709, tmp.sc, phenotype2, alpha = 0.05, tag = tag,
                        family = "binomial", #二分类
                        Save_file = paste0("./Scissor/",GSEs[i],"_.Rdata"))
  infosy <- Scissor_SCT(GSE199709, tmp.sc, phenotype2, alpha = NULL, cutoff = 0.03, 
                        family = "binomial", Load_file = paste0("./Scissor/",GSEs[i],"_.Rdata"))
  LSS_lists_A[[i]] <- infosx
  LSS_lists_B[[i]] <- infosy
  LSS_lists_pos <- c(LSS_lists_pos,LSS_lists_B$Scissor_pos)
  LSS_lists_neg <- c(LSS_lists_neg,LSS_lists_B$Scissor_neg)
}

Scissor_select_LSS <- rep(0, ncol(scw))#创建一个列表，用来表示4102个细胞，
names(Scissor_select_LSS) <- colnames(scw)#给列表中每一个数赋予细胞编号
Scissor_select_LSS[LSS_lists_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_select_LSS[LSS_lists_neg] <- 2#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = Scissor_select_LSS, col.name = "scissor_LSS_split")#将表示4102个细胞分类的列表添加到sc_dataset这个Seurat对象中
DimPlot(scw, reduction = 'umap',split.by = "GSE", group.by = 'scissor_LSS_split',
        pt.size = 1.2, order = c(2,1),label = T,cols = c("Grey","RED","Blue"))#可视化


# GSE120521 LSS split GSE
US_lists_A <- list()
US_lists_B <- list()
US_lists_pos <- character()
US_lists_neg <- character()
total <- scw@graphs[["SCT_snn"]]

for (i in 1:length(GSEs)) {
  tmp.sc <- subset(scw,GSE==GSEs[i])
  cellnames <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
  matrixsel <-total[cellnames,cellnames] 
  tmp.sc@graphs[["SCT_snn"]] <- matrixsel
  infosx <- Scissor_SCT(bulk_dataset, tmp.sc, phenotype, alpha = 0.05, tag = tag,
                        family = "binomial", #二分类
                        Save_file = paste0("./Scissor/",GSEs[i],"_.Rdata"))
  infosy <- Scissor_SCT(bulk_dataset, tmp.sc, phenotype, alpha = NULL, cutoff = 0.03, 
                        family = "binomial", Load_file = paste0("./Scissor/",GSEs[i],"_.Rdata"))
  US_lists_A[[i]] <- infosx
  US_lists_B[[i]] <- infosy
  US_lists_pos <- c(US_lists_pos,infosy$Scissor_pos)
  US_lists_neg <- c(US_lists_neg,infosy$Scissor_neg)
}









##########################################################################
#StringDB
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
# BiocManager::install("STRINGdb")
library(STRINGdb)
library(igraph)
# detach("package:igraph", unload = TRUE)
library(ggraph)


LSS_genes.fli <- LSS_trajectory_gene[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

LSS_genes_sig100 <- LSS_genes.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

US_genes.fli <- US_trajectory_gene[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

US_genes_sig100 <- US_genes.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# 需要分析的基因

LSS_gene.con <- LSS_genes_sig100 %>% bitr(fromType = "SYMBOL",
                                  toType = c("ENSEMBL","ENTREZID"),
                                  OrgDb = "org.Hs.eg.db",
                                  drop = T)
LSS_gene.con.GO <- enrichGO(gene = LSS_gene.con$ENSEMBL,
                          OrgDb = org.Hs.eg.db,keyType = "ENSEMBL",pvalueCutoff = 0.05)

LSS_gene.con.GO <- clusterProfiler::simplify(LSS_gene.con.GO)
LSS_gene.con.GO <- pairwise_termsim(LSS_gene.con.GO)
LSS_gene.con.GO <- setReadable(LSS_gene.con.GO,OrgDb = org.Hs.eg.db)



US_gene.con <- US_genes_sig100 %>% bitr(fromType = "SYMBOL",
                                          toType = c("ENSEMBL","ENTREZID"),
                                          OrgDb = "org.Hs.eg.db",
                                          drop = T)
US_gene.con.GO <- enrichGO(gene = US_gene.con$ENSEMBL,
                            OrgDb = org.Hs.eg.db,keyType = "ENSEMBL",pvalueCutoff = 0.05)

US_gene.con.GO <- clusterProfiler::simplify(US_gene.con.GO)
US_gene.con.GO <- pairwise_termsim(US_gene.con.GO)
US_gene.con.GO <- setReadable(US_gene.con.GO,OrgDb = org.Hs.eg.db)

gene.pot <- intersect(US_genes_sig100,LSS_genes_sig100)




###############################################################
# 创建STRINGdb对象
string_db <- STRINGdb$new( version="11.5", species=9606,
                           score_threshold=400, network_type="physical")
gene <- gene.pot %>% bitr(fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db",
                          drop = T)
data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",
                                      removeUnmappedRows = FALSE)
string_db$plot_network(data_mapped$STRING_id)
data_links <- data_mapped$STRING_id[!is.na(data_mapped$STRING_id)] %>% string_db$get_interactions()

# 转换stringID为Symbol，只取前两列和最后一列
links <- data_links %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)

# 节点数据
nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
# 创建网络图
# 根据links和nodes创建
net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
# 添加一些参数信息用于后续绘图
# V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
igraph::V(net)$size <- igraph::degree(net)/5 #
igraph::E(net)$width <- igraph::E(net)$weight/10
# 使用ggraph绘图
# ggraph是基于ggplot2的包，语法和常规ggplot2类似
ggraph(net,layout = "kk")+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>0, label=name), size = 5, repel = T)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size="none")+
  theme_graph()





################################################################################
load("E:/works/artery ECs/GSE208194_genes.RData")
intersect(fli$Gene.name,Track_genes_sig100)

table(scw$GSE,scw$mono.clu)


FeaturePlot(scw,features = "PROX1",label = T)
FeaturePlot(scw,features = "MMRN1")

################################################################################

for (i in 1:length(GSEs)) {
  sctmp <- subset(scw,GSE==GSEs[i])
  test <- GetAssayData(sctmp,slot = "count",assay = "RNA")
  write.csv(test,file = paste0("./exprs/",GSEs[i],"_expr.csv"))
}

write.csv(ExportGeneExpr, file='./input/Seurat_geneExpr.csv', row.names = T)
test=as.data.frame(test)

################################################################################
#scFEA

GSE155468_FLUX <- read.csv('./scFEA/GSE155468_flux.csv', header = T, row.names = 1)
GSE155468_FLUX <- t(data.matrix(GSE155468_FLUX))
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

scFEA <- list()
for (i in 1:length(GSEs)) {
  flux_tmp <- read.csv("./scFEA/",GSEs[i],"_flux.csv", header = T, row.names = 1)
  flux_tmp <- t(data.matrix(flux_tmp))
  obj <- subset(scw,GSE==GSEs[i])
  obj[["FLUX"]] <- CreateAssayObject(counts = flux_tmp)
  obj@assays$SCT <- NULL
  DefaultAssay(obj) <- 'FLUX'
  scFEA[[i]] <- obj
}




rm(list=ls())
gc()

