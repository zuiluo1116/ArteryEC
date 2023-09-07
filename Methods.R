################################################################################
# Figures

##  Library
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
library(SeuratDisk)
library(Seurat)
library(SeuratObject)
library(tidyverse)

library(harmony)
library(sphereplot)
library(monocle3,lib.loc = "D:/Program Files/R/R-4.2.1/library")
library(densitycut)
library(robustbase)
library(RColorBrewer)
library(SCEnt)
library(ggridges)
library(Nebulosa)

library(scRNAtoolVis)
library(ggplotify)
library(corrplot)
library(scales)
library(entropy)

library(irGSEA, lib.loc = "D:/Program Files/R/R-4.2.1/library")
library(AUCell)
library(UCell)
library(GSVA)
library(reshape2)

library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)

library(ggpubr)
library(ggplot2)

library(UpSetR)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOplot)

library(UpSetR)
library(ggupset)
library(httr)
library(jsonlite)
library(RobustRankAggreg)

library(Scissor)

library(ggforce)
library(tidyr)
library(ggalluvial)

max(scw$nCount_RNA)
min(scw$nCount_RNA)
min(scw$nFeature_RNA)
max(scw$nFeature_RNA)
VlnPlot(scw,features = c("nCount_RNA","nFeature_RNA"),group.by = "GSE")
VlnPlot(scRNA,features = c("nCount_RNA","nFeature_RNA"),group.by = "GSE")

minGene=1000
maxGene=4000
maxUMI=15000
pctMT=10
pctHB=1
scRNA <- subset(scw,nCount_RNA < maxUMI 
                & nFeature_RNA > minGene 
                & nFeature_RNA < maxGene)
DimPlot(scRNA,reduction = "scphere",group.by = "NEW.clu")
DimPlot(scw,reduction = "scphere",group.by = "NEW.clu")
FeaturePlot(scRNA,reduction = "scphere",features = c("nCount_RNA","nFeature_RNA"))

getwd()
load("E:/works/artery ECs/1213.Rdata")
scRNA <- LoadH5Seurat("scRNA_filtered_harmony.h5seurat")

VlnPlot(scRNA,features = c("MYH11","CDH5","DCN","PTPRC"),stack = T)

## GSE155468
load("E:/works/GSE155468/0701.RData")

sce <- merge(x=sceList[[1]],y=(c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                                 sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]],sceList[[10]])))
sce <- SCTransform(sce)

sce <- RunPCA(sce, npcs=50, verbose=FALSE)

### 整合方法1：单个样本间进行整合（推荐，效果更好）
sce <- RunHarmony(sce, group.by.vars="orig.ident", 
                  assay.use="SCT", max.iter.harmony = 20)

ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
sce <- RunUMAP(sce,reduction="harmony", dims=pc.num)
DimPlot(sce,reduction = "umap")

FeaturePlot(sce,features = c("MYH11","CDH5","DCN","PTPRC"))

getwd()
ggsave(filename = "GSE155468.pdf",width = 8,height = 6)

## GSE155514
load("E:/works/GSE155514/Carotid human/human.Rdata.RData")
FeaturePlot(sce,features = c("MYH11","CDH5","DCN","PTPRC"))
DefaultAssay(sce) <- "RNA"
ggsave(filename = "GSE155514.pdf",width = 8,height = 6)


## GSE159677
dir <- list.files("E:/works/GSE159677/raw")
sample_names = stringr::str_replace_all(dir," ",".")

scRNAlist1 <- list()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = paste0("E:/works/GSE159677/raw/",dir[i]))
  #不设置min.cells过滤基因会导致CellCycleScoring报错：
  #Insufficient data values to produce 24 bins.  
  scRNAlist1[[i]] <- CreateSeuratObject(counts, project=sample_names[i],)
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist1[[i]] <- RenameCells(scRNAlist1[[i]], add.cell.id = sample_names[i])   
  #计算线粒体基因比例
  if(T){    
    scRNAlist1[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist1[[i]], pattern = "^mt-") 
  }
  #计算核糖体基因比例
  if(T){
    scRNAlist1[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist1[[i]], pattern = "^Rp[sl]")
  }
  #计算红细胞基因比例
  if(T){
    HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb-bt","Hbb-bs","Hbb-bh1","Hbb-bh2","Hbb-y","Hbq1a","Hbq1b")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist1[[i]]))
    scRNAlist1[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist1[[i]], features=HB.genes) 
  }
}

scRNA <- merge(scRNAlist1[[1]], y=c(scRNAlist1[[2]], scRNAlist1[[3]], scRNAlist1[[4]], scRNAlist1[[5]], scRNAlist1[[6]]))

scRNA$orig.ident

gc()

scRNA <- SCTransform(scRNA)
##2.3 使用harmony整合数据</li>
### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

### 整合方法1：单个样本间进行整合（推荐，效果更好）
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
# group.by.vars参数是设置按哪个分组来整合
# max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
#??????RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）

#   ###整合方法2：按其他分类如不同分组来校正（实测效果不如方法1）
# if(F){
#   scRNA$batches <- scRNA$GSE
#   scRNA2 <- RunHarmony(scRNA, group.by.vars="GSE", assay.use="SCT")
# }

## 2.4 降维聚类
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunUMAP(scRNA,reduction="harmony", dims=pc.num)
FeaturePlot(scRNA,features = c("MYH11","CDH5","DCN","PTPRC"))
ggsave(filename = "GSE159677.pdf",width = 8,height = 6)

## GSE213740
load("E:/works/artery ECs/GSE213740.RData")

scRNAlist1 <- GSE213740.list[c(1,2,3,7,8,9)]
scRNA <- merge(scRNAlist1[[1]], y=c(scRNAlist1[[2]], scRNAlist1[[3]], scRNAlist1[[4]], scRNAlist1[[5]], 
                                    scRNAlist1[[6]]))
rm(GSE213740.list)
rm(scRNAlist1)
gc()
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunUMAP(scRNA,reduction="harmony", dims=pc.num)
FeaturePlot(scRNA,features = c("MYH11","CDH5","DCN","PTPRC"))
library(ggplot2)
ggsave(filename = "GSE213740.pdf",width = 8,height = 6)


## GSE216860
load("E:/works/artery ECs/GSE216860.RData")
scRNAlist1 <- GSE216860.list
scRNA <- merge(scRNAlist1[[1]], y=c(scRNAlist1[[2]], scRNAlist1[[3]], scRNAlist1[[4]], scRNAlist1[[5]], 
                                    scRNAlist1[[6]]))
rm(GSE216860.list)
rm(scRNAlist1)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunUMAP(scRNA,reduction="harmony", dims=pc.num)
FeaturePlot(scRNA,features = c("MYH11","CDH5","DCN","PTPRC"))
ggsave(filename = "GSE216860.pdf",width = 8,height = 6)



















scw <- SeuratDisk::LoadH5Seurat("20221201.h5seurat")
DimPlot(scw,group.by = "GSE",shuffle = T)
ggsave("UMAP_GSE.pdf",height = 6,width = 8)

DimPlot(sct,split.by = "orig.ident")
FeaturePlot(scw,features = "EPCAM",reduction = "scphere")
VlnPlot(scw,features = "PLAUR",group.by = "NEW.clu",assay = "RNA",slot = "count")

VlnPlot(scw,features = "ICAM1",group.by = "celltype")

Nebulosa::plot_density(scw,features = "PLAUR",reduction = "scphere")

# 构建scphere RNA表达矩阵
test <- GetAssayData(scw,slot = "count",assay = "RNA")
# 输出 NOT RUN AGAIN
# Matrix::writeMM(test,"scw.mtx")
# write.table(test@Dimnames[[1]],file = 'genes.tsv',
#             quote = F,sep = '\t',
#             col.names = F,row.names = F)
# write.table(test@Dimnames[[2]],file = 'barcodes.tsv',quote = F,
#             col.names = F,row.names = F)
#### Python scPhere model

#### R visualization

###sFig.1 3D下降维展示
# 3D可视化 定义功能 基于rgl
PlotSphere = function(x, cluster, col, density=FALSE, legend=FALSE) {
  if (missing(col)) {
    col = distinct.col
  }
  if (density == FALSE) {
    col.point = densitycut::AssignLabelColor(col, cluster, 
                                             uniq.label = sort(unique(cluster)))
  } else {
    colours = colorRampPalette((brewer.pal(7, "YlOrRd")))(10)
    FUN = colorRamp(colours)
    cluster = (cluster - min(cluster)) / diff(range(cluster))
    col.point = rgb(FUN(cluster), maxColorValue=256)
  }
  plot3d(x[, 1:3], col = col.point, 
         xlim=c(-1, 1), ylim=c(-1, 1), zlim=c(-1, 1), 
         box=FALSE, axes=FALSE, xlab='', ylab='', zlab='')
  
  arrow3d(c(0, -1.35, 0), c(0, 1.35, 0), 
          col = 'gray', s=0.04, type='extrusion', lit=FALSE)
  
  spheres3d(0, 0, 0, lit=FALSE, color='#dbd7d2', 
            alpha=0.6, radius=0.99)
  spheres3d(0, 0, 0, radius=0.9999, lit=FALSE, color='gray', 
            front='lines', alpha=0.6)
  
  if (density == FALSE) {
    id = !duplicated(cluster)
    col.leg = AssignLabelColor(col, cluster)[id]
    leg = cluster[id]
    names(col.leg) = leg
    
    if (legend == TRUE) {
      legend3d("topright", legend = leg, 
               pch = 16, col = col.leg, cex=1, inset=c(0.02)) 
    }
    
    cluster.srt = sort(unique(cluster))
    k.centers = sapply(cluster.srt, function(zz) {
      cat(zz, '\t')
      id = cluster == zz
      center = colMedians(as.matrix(x[id, , drop=FALSE]))
    })
    
    k.centers = t(k.centers)
    
    cluster.size = table(cluster)[as.character(cluster.srt)]
    id = which(cluster.size > 0)
    
    if (length(id) > 0) {
      k.centers = k.centers[id, , drop=FALSE]
      cluster.srt = cluster.srt[id]
    }
    
    k.centers = k.centers / sqrt(rowSums(k.centers^2)) * 1.15
    text3d(cex=6,k.centers, texts=cluster.srt, col='black')
  }
}

# 读取3D坐标
x = read.delim('./scphere/scw_latent_250epoch.tsv', 
               sep=' ', header=FALSE)
# 转换2D坐标
y <- car2sph(x)

# 坐标命名为细胞barcode 注意对应
rownames(x) <- test@Dimnames[[2]]
rownames(y) <- test@Dimnames[[2]]

# 构建monocle对象
load("E:/works/artery ECs/monocle.Rdata")
data <- GetAssayData(scw,assay ='RNA',slot ='counts')
cell_metadata <- scw@meta.data
gene_annotation <- data.frame(gene_short_name =rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)

# 主成分
cds <- preprocess_cds(cds, num_dim = 50)
# plot_pc_variance_explained(cds)
# 降维 默认UMAP
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA

# 替换2D UMAP坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
scphere.embed <- y
# 传入monocle 替换UMAP坐标
scphere.embed <- scphere.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- scphere.embed
# plot_cells(cds,show_trajectory_graph = F)

# 传入seurat 创建scphere降维
seu.embed <- Embeddings(scw,reduction = "umap")
scphere.embed <- scphere.embed[rownames(seu.embed),]

scw@reductions[["scphere"]] <- CreateDimReducObject(embeddings = scphere.embed,
                                                    assay = "SCT",global = T,key = "scPhere_")

# 构建3D Monocle3对象
cdsk_3d <- cds
# 3D trajectories（4步走）
cdsk_3d = reduce_dimension(cdsk_3d, max_components = 3)
cdsk_3d.embed <- cdsk_3d@int_colData$reducedDims$UMAP
scphere_3d.embed <- x
scphere_3d.embed <- scphere_3d.embed[rownames(cdsk_3d.embed),1:3] #排序
colnames(scphere_3d.embed) <- colnames(cdsk_3d.embed)
cdsk_3d@int_colData$reducedDims$UMAP <- scphere_3d.embed
plot_cells_3d(cdsk_3d,color_cells_by = "celltype",
              color_palette = scales::hue_pal()(10))
cds <- order_cells(cds)
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')


order_cells(cdsk_3d,reduction_method = "UMAP",root_cells = )

pseudotime(cdsk_3d) <- pseudotime
str(cdsk_3d)

# 识别分群
cdsk_3d = cluster_cells(cdsk_3d)
cdsk_3d = cluster_cells(cdsk_3d,resolution=0.0002)
plot_cells_3d(cdsk_3d)
# 2 monocle_cds <- cluster_cells(monocle_cds, resolution = 0.0001, reduction_method = "umap")
cdsk_3d <- learn_graph(cdsk_3d)
cdsk_3d <- order_cells(cdsk_3d)

sc.df <- as.data.frame(cdsk_3d@colData@listData[["scissor"]])
rownames(sc.df) <- cdsk_3d@assays@data@listData[["counts"]]@Dimnames[[2]]

time.df <- as.data.frame(cdsk_3d@colData@listData[["pseudotime"]])
time.df <- time.df[rownames(sc.df),]
cdsk_3d@colData@listData[["pt"]] <- time.df


plot_cells_3d(cdsk_3d,color_cells_by = "pt")

## TOP 50 genes
pt.sel <- Track_genes.fli[1:50,]
exp <- AverageExpression(scw,assays = "RNA",slot = "counts",group.by = "celltype")[[1]]
expsel <- exp[pt.sel$gene_short_name,]
pt.sel <- cbind(pt.sel,expsel)
write.csv(pt.sel,"pt.sel.csv")




min(pseudotime)
max(pseudotime)

FeaturePlot(scw,"PRCP",reduction = "scphere",min.cutoff = "q30")
VlnPlot(scw,"VCAM1",group.by = "celltype")
VlnPlot(scw,"FGF18",group.by = "celltype")

library(viridis)
c1=viridis(25653, alpha = 1, begin = 0, end = 1, direction = 1, option = "A")
length(unique(pseudotime))

df.3d <- data.frame("clusters" = cdsk_3d@clusters@listData[["UMAP"]][["clusters"]],
                    "names" = cdsk_3d@colData@rownames)

# df.3d <- df.3d[cdsk_3d@colData@rownames,]
# cdsk_3d@colData@listData[["3d.clu"]] <- df.3d$clusters
# df.3d <- df.3d[cds@colData@rownames,]
# cds@colData@listData[["3d.clu"]] <- df.3d$clusters
# plot_cells(cdsk_3d,color_cells_by = "X3d.clu")

# 分群回传seurat
df.3d <- df.3d[scw@assays[["RNA"]]@counts@Dimnames[[2]],]
save(df.3d,file = "NEW.clu.RData")
load("E:/works/artery ECs/NEW.clu.RData")

scw@meta.data$`NEW.clu` <- df.3d$clusters
DimPlot(scw,reduction = "scphere",group.by = "NEW.clu")
clus_num <- as.integer(scw$NEW.clu)


# 3D 可视化
PlotSphere(x, clus_num,scales::hue_pal()(10))


# 调整视角
view_start <- structure(c(0.8660254, -0, 0.5000000, 0,
                          
                          0, 1, 0, 0,
                          
                          -0.5, -0, 0.8660254 , 0,
                          
                          0, 0, 0, 1), .Dim = c(4L, 4L))

view_end <- structure(c(-0.8660254, -0, -0.5000000, 0,
                        
                        0, 1, 0, 0,
                        
                        0.5, -0, -0.8660254 , 0,
                        
                        0, 0, 0, 1), .Dim = c(4L, 4L))

par3d(userMatrix = view_start)
# 指定视角 出动图
movie3d(spin3d(axis=c(0,1,0),rpm=7.2),duration=10,fps=50,movie = "3dNEWplot",dir = "./Figures/")

# 指定视角 出图
par3d(userMatrix = view_start)
rgl.postscript('./Figures/3dplot.part1.pdf', fmt = 'pdf')
par3d(userMatrix = view_end)
rgl.postscript('./Figures/3dplot.part2.pdf', fmt = 'pdf')

# reduction dimension and cluster visualization
dmplot <- DimPlot(scw,group.by = "NEW.clu",reduction = "scphere",
                  label = T,label.size = 8,label.box = T) +
  NoLegend()+ggtitle(label = "")
ggsave(dmplot,filename = "Figures/Fig.1A.pdf",height = 6,width = 8)


# distribution among GSE datasets
sptplot <- DimPlot(scw,group.by = "NEW.clu",reduction = "scphere",split.by = "GSE",ncol = 3) +
  NoLegend()+ggtitle(label = "")
ggsave(sptplot,filename = "Figures/Fig.1B.pdf",height = 6,width = 8)


### 随机抽样 计算熵
scw

Idents(scw) <- "GSE"
entropy <- numeric()
sample_entropy <- data.frame()

for (j in 1:200) {
  set.seed(1116+j)
  sctmp <- scw[,sample(1:25869, 2000)]
  cellnum1 <- as.matrix(t(cbind(table(sctmp$GSE,sctmp$NEW.clu))))
  
  for (i in 1:nrow(cellnum1)){ 
    entropy[i]=gene_hom(as.numeric(cellnum1[i,]))
  }
  sample_entropy <- rbind(sample_entropy,entropy)
}

# 
# cellnum1
# entropy2 <- as.numeric()
# for (i in 1:ncol(cellnum1)){ 
#   entropy2[i]=gene_hom(as.numeric(cellnum1[,i]))
# }


names(sample_entropy) <- levels(as.factor(scw$NEW.clu))
# sample_entropy <- sample_entropy[-1,]
entropy_tmp <- as.matrix(sample_entropy)
# entropy_tmp <- (entropy_tmp)/(max(entropy_tmp))
entropy_tmp <- (entropy_tmp-min(entropy_tmp))/(max(entropy_tmp)-min(entropy_tmp))

# sample_test <- melt(sample_entropy)
sample_test <- melt(entropy_tmp)

table(sample_test$Var2)
# 
# class(sample_test$variable)
# class(sample_test$Var2)
sample_test$Var2 <- as.factor(sample_test$Var2)

group_heter <- ggplot(sample_test, aes(x = value , y = Var2 , fill = Var2)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none")
ggsave(group_heter,filename = "Figures/Fig.1C.pdf",height = 6,width = 2)

group_heter

# 
# shannon.entropy <-function(x,type='raw'){
#   if(type=='raw'){
#     myfreqs <- table(x)/length(x) 
#     myvec <- as.data.frame(myfreqs)[,2]
#   }else{
#     myvec=x
#   }
#   -sum(myvec * log2(myvec))
# }
# metric.entropy <-function(x,type='raw'){
#   if(type=='raw'){
#     myfreqs <- table(x)/length(x) 
#     myvec <- as.data.frame(myfreqs)[,2]
#   }else{
#     myvec=x
#   }
#   -sum(myvec * log(myvec,length(x)))
# }

# 
# # heatmap plot for unique expression gene
Idents(scw) <- "NEW.clu"

# FAMs_3d <- FindAllMarkers(scw,logfc.threshold = 0.4,slot = "count")

FAMs_1 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 1,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_2 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 2,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_3 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 3,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_4 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 4,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_5 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 5,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_6 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 6,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_7 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 7,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_8 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 8,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_9 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 9,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")
FAMs_10 <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = 10,only.pos = T,
                               logfc.threshold = 0.4,slot = "data")

act <- list(
  FAMs_1=rownames(FAMs_1),
  FAMs_2=rownames(FAMs_2),
  FAMs_3=rownames(FAMs_3),
  FAMs_4=rownames(FAMs_4),
  FAMs_5=rownames(FAMs_5),
  FAMs_6=rownames(FAMs_6),
  FAMs_7=rownames(FAMs_7),
  FAMs_8=rownames(FAMs_8),
  FAMs_9=rownames(FAMs_9),
  FAMs_10=rownames(FAMs_10)
)

length(unlist(act))
actallgenes <- unique(unlist(act))

FAMs_3d_scale_NEW <- FindAllMarkers(scw,logfc.threshold = 0.4,slot = "data")
FAMs_3d_scale_fli <- FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$gene %in% actallgenes),]
load("FAMS.RData")
save(list = c("act","FAMs_3d_scale_fli","FAMs_3d_scale_NEW"),file = "FAMS.RData")

vocaplot <- markerVocalno(markers = FAMs_3d_scale_fli,
                          labelCol = scales::hue_pal()(10))
vocaplot
ggsave(vocaplot,filename = "Figures/Fig.1D.pdf",height = 8,width = 16)

Idents(scw) <- "NEW.clu"
fea_lymEC <- FeaturePlot(scw,features = c("LYVE1","CCL21"),blend = T,reduction = "scphere")
ggsave(fea_lymEC,filename = "Figures/sFig.2A.pdf",height = 6,width = 32)

vln_lymEC <- VlnPlot(scw,features = c("LYVE1","CCL21"),combine = T)
ggsave(vln_lymEC,filename = "Figures/sFig.2B.pdf",height = 6,width = 16)


lym_featplot <- FeaturePlot(scw,features = c("LYVE1","CCL21"),blend = T,reduction = "scphere")

ggsave(lym_featplot,filename = "Figures/sFig.1C.pdf",height = 6,width = 22)

help(plot_density)




Idents(scw) <- "NEW.clu"
expr.clu <- AverageExpression(scw,slot = "data",assay = "RNA")[[1]]# 1, 2, 3-7
# expr.clu <- AverageExpression(scw,slot = "count",assay = "SCT")[[1]] #1, 2, 4&7,3&5&6
expr_feature <- expr.clu[actallgenes,]
cor_X3d.clu <- cor(expr_feature, method = c("spearman"))
# cor_X3d.clu <- cor(expr.clu, method = c("spearman"))
# cor_X3d.clu <- cor(expr.clu, method = )

# Saveplot
# 1. 画布
pdf( 
  file = "./Figures/Fig.2A.pdf", # 文件名称
  width = 6,           # 宽
  height = 6)              # 分辨率
# 2. 绘图
corrplot <- corrplot(cor_X3d.clu, type = "lower", order = "FPC", 
                     tl.col = "black", tl.srt = 45,addCoef.col = T)
# 3. 关闭画布
dev.off()



## cluster rename
df.newtype <- data.frame(
  clu = scw$NEW.clu
)

df.newtype[which(df.newtype$clu %in% c(1,2,4,6,9,10)),2] <- "ACKR1+ EC"
df.newtype[which(df.newtype$clu %in% c(5,7)),2] <- "TCIM+ EC"
df.newtype[which(df.newtype$clu %in% c(3)),2] <- "ITLN1+ EC"
df.newtype[which(df.newtype$clu %in% c(8)),2] <- "Lymphatic EC"

scw$celltype <- factor(df.newtype$V2,levels = c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))

dmplot_new <- DimPlot(scw,group.by = "celltype",reduction = "scphere")+NoLegend()
ggsave(dmplot_new,filename = "Figures/Fig.2B.pdf",height = 6,width = 8)

## proportion bar plot
cellnum1 <-table(scw$GSE,scw$celltype)
group<-rownames(cellnum1)
celltype <- colnames(cellnum1)
cell.prop<-as.data.frame(prop.table(t(cellnum1)))
colnames(cell.prop)<-c("Celltype","sample","proportion")
data4plot1 <-as.data.frame(cell.prop)

p1<-ggplot(cell.prop,aes(sample,proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=scales::hue_pal()(4))+#自定义fill的颜色
  ggtitle("cell proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+NoLegend()
p1
ggsave(p1,filename = "Figures/Fig.2B2.pdf",height = 8,width = 2)
VlnPlot(scw,features = "ACKR1",group.by = "NEW.clu")

# sck <- DietSeurat(scw)
sck <- subset(scw,celltype!="Lymphatic EC")
select.cells <- sck@assays[["RNA"]]@data@Dimnames[[2]]
save(select.cells,file = "vas_EC.Rdata")
# tmp <- ls()
# tmp <- tmp[-which(tmp %in% c("scw","cds","cdsk_3d"))]
# save(list = tmp,file = "1218.RData")


intersect(rownames(FAMs_5),rownames(FAMs_7))

den_vasEC <- Nebulosa::plot_density(scw,features = c("ITLN1"),
                                    reduction = "scphere",size = 0.4)+
  NoLegend()+ 
  Nebulosa::plot_density(scw,features = c("TCIM"),
                         reduction = "scphere",size = 0.4)+
  NoLegend()+
  Nebulosa::plot_density(scw,features = c("ACKR1"),
                         reduction = "scphere",size = 0.4)+
  NoLegend()+
  Nebulosa::plot_density(scw,features = c("CCL21"),
                         reduction = "scphere",size = 0.4)+
  NoLegend()

ggsave(den_vasEC,filename = "Figures/Fig.2C.pdf",height = 6,width = 12)
library(Seurat)
library(ggplot2)
rid_plot <- RidgePlot(scw,features = c("ACKR1","ITLN1","TCIM","CCL21"),
          same.y.lims = F,stack = T,slot = "data",group.by = "celltype",fill.by = "ident",log = T)+NoLegend()
ggsave(rid_plot,filename = "Figures/Fig.2C2N.pdf",height = 4,width = 8)

# Nebulosa::plot_density(scw,"RGCC",reduction = "scphere")
# FeaturePlot(scw,features = "RGCC",reduction = "scphere")

# expr_feature <- expr.clu[actallgenes,]
# cor_unique_gene <- cor(expr_feature, method = c("spearman"))
# corrplot <- corrplot(cor_X3d.clu, type = "lower", order = "FPC", 
#                      tl.col = "black", tl.srt = 45,addCoef.col = T)

Idents(scw) <- "celltype"
DotPlot(scw,features = c("ITLN1","TCIM","ACKR1","CCL21"))
FeaturePlot(scw,"ACKR1",reduction = "scphere",min.cutoff = "q30") #clus 1 2 4 6 9 10
FeaturePlot(scw,"ITLN1",reduction = "scphere",min.cutoff = "q30") #clus 3 
FeaturePlot(scw,"RGCC",reduction = "scphere",min.cutoff = "q30") # clus 5 7
dot_feaplot <- jjDotPlot(object = scw,
          gene = c("ITLN1","TCIM","ACKR1","CCL21"),id = "NEW.clu",ytree = T)
ggsave(dot_feaplot,filename = "Figures/sFig.1D.pdf",height = 6,width = 8)

# FeaturePlot(scw,features = c("ITLN1","ADGRF5"),blend = T,reduction = "scphere",min.cutoff = "q30")




# Circlize
# library(ComplexHeatmap)
# library(circlize)
# library(dendextend)
# library(dendsort)
# BiocManager::install("dendsort")
Idents(scw) <- "celltype"
FAMs_type <- FindAllMarkers(scw,logfc.threshold = 0.4,slot = "data")
FAMs_type$diff.pct <- FAMs_type$pct.1-FAMs_type$pct.2
FAMs_type.sel <- FAMs_type %>% 
  group_by_("cluster") %>% 
  filter(diff.pct>0.15 & avg_log2FC>0.5 &  pct.1 >0.3) %>% 
  arrange("diff.pct") %>% 
  slice_max(order_by = diff.pct, n = 20)
save(list = c("FAMs_type","FAMs_type.sel"),file = "FAMS_type.RData")

table(FAMs_type.sel$cluster)
expr.count <- AverageExpression(scw,slot = "scale.data",assay = "RNA")[[1]]#
hm.count <- expr.count[FAMs_type.sel$gene,]

# load("FAMS_type.RData")
#根据数据范围定义热图的颜色梯度：
range(hm.count)
mycol <- colorRamp2(c(-0.6740817, 0, 1.7099878),c("#0da9ce", "white", "#e74a32"))

#环形热图绘制：
circos.heatmap(hm.count,
               col = mycol)
circos.clear()#绘制完成后，需要重置循环布局再开启下一次绘图，不然会遇到报错或轨道继续在同一图中叠加；

circos.par(gap.after = c(50)) #调整圆环首位间距；
circos.heatmap(hm.count,
               col = mycol,
               cluster = TRUE, #是否对行聚类
               dend.side = "inside", #聚类树方向：inside显示在圆环内圈，inside为圆环外圈；
               rownames.side ="outside", #行名方向；
               rownames.col = "black",
               rownames.cex = 0.6)
circos.clear()
# 1. 画布
pdf( 
  file = "./Figures/Fig.2D.pdf", # 文件名称
  width = 6,           # 宽
  height = 6)              # 分辨率
# 2. 绘图
#聚类树美化：
circos.par(gap.after=c(0))
circos.heatmap(hm.count,
               col = mycol,
               cluster = T,
               dend.side = "inside",
               rownames.side = "outside",
               rownames.col = "black",
               rownames.cex = 0.6,
               track.height = 0.3, #轨道高度，即圆环/格子的粗细
               dend.track.height = 0.08,#聚类树高度调整
               dend.callback = function(dend, m, si) { #聚类树的回调
                 # color_branches(dend,k = 4,col = 1:4) #修改聚类树颜色
                 color_branches(dend,k = 4,col = scales::hue_pal()(4)) #修改聚类树颜色
               })

lg <- Legend(title = "Scaled Expression",
             col_fun = mycol,
             direction = c("horizontal"),
             title_position = c('topcenter'))
draw(lg, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))
circos.clear()
# 3. 关闭画布
dev.off()

## sample state rename
## cluster rename
scw$state <- paste0(scw$GSE,"_",scw$orig.ident)
df.newstate <- data.frame(
  clu = scw$state
)
table(df.newstate$clu)
orig_states <- levels(factor(scw$state))

df.newstate[which(df.newstate$clu %in% orig_states[c(1:3,16,18,20,26:32)]),2] <- "Control"
df.newstate[is.na(df.newstate$V2),2] <- "Burden"
scw$state <- factor(df.newstate$V2,levels = c("Control","Burden"))
DimPlot(scw,group.by = "NEW.clu",split.by = "state",reduction = "scphere")

## GSE159677 own-control
sct <- subset(scw,GSE == "GSE159677")

cellnum1 <-table(sct$orig.ident,sct$celltype)
group<-rownames(cellnum1)
celltype <- colnames(cellnum1)
cell.prop<-as.data.frame(prop.table(t(cellnum1)))
colnames(cell.prop)<-c("Celltype","sample","proportion")

cell.prop$sample <- str_sub(cell.prop$sample,start = 10,end = 11)

prop.plot <- ggplot(cell.prop, aes(x=sample, y=proportion, fill = sample)) + 
  geom_boxplot()+
  # geom_violin()+
  # geom_jitter(stroke = 0.01,width = 0.2,size=3)+
  geom_point(size=3)+
  # geom_line(data=filter(mydata,treatment=="v"),aes(group==ID))+
  facet_wrap(~Celltype, scale = "free_x") +
  scale_y_continuous(name = "Clusters Score") +
  scale_x_discrete(name = "Continent") +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  theme_bw()
ggsave(prop.plot,filename = "Figures/Fig.3A1.pdf",height = 6,width = 8)

## paired boxplot
# library(ggpubr)
# library(ggplot2)
compare_means(proportion ~ Celltype, data = cell.prop,paired = T)
p <- ggpaired(cell.prop, x = "sample", y = "proportion",linetype = 2,
         color = "Celltype", 
         line.color = "gray", line.size = 0.4,
         facet.by = "Celltype", short.panel.labs = FALSE)

# p + stat_compare_means(label = "p.format", paired = TRUE)
ggsave(p,filename = "Figures/Fig.3A.pdf",height = 6,width = 8)

## dimplot GSE159677
dmp_GSE159677 <- DimPlot(sct,split.by = "state",reduction = "scphere")+NoLegend()
ggsave(dmp_GSE159677,filename = "Figures/sFig.2A.pdf",height = 6,width = 8)

FeaturePlot(scw,features = "PROX1",reduction = "scphere")


####### Scissor
load("bulk.RData")

# library(Scissor)
Scissor_SCT <- function (bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                         cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                         Save_file = "Scissor_inputs.RData", Load_file = NULL) {
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


## GSE199709 
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

phe_199709 <- c(rep(1,3),rep(0,3))
tag_199709 <- c("LSS","Control")

sck <- subset(scw,GSE=="GSE159677")
sck <- sck@assays$RNA@scale.data
sck@graphs$SCT_snn <- scw@graphs$SCT_snn[sck@assays[["RNA"]]@counts@Dimnames[[2]],
                                         sck@assays[["RNA"]]@counts@Dimnames[[2]]]
# head(scw@graphs$SCT_snn)[1:4,1:4]
gc()
pri_199709 <- Scissor_SCT(GSE199709, sck, phe_199709,cutoff = 0.3, tag = tag_199709,
                      family = "binomial", #二分类
                      Save_file = './Scissor/Sci_199709.RData')
gc()
# 
# 
# pt3_199709 <- Scissor_SCT(GSE199709, sck, phe_199709, alpha = NULL, cutoff = 0.1, 
#                       family = "binomial", Load_file = './Scissor/Sci_199709.RData')

select_199709 <- rep(0, ncol(scw))#创建一个列表，用来表示scRNA中细胞，
names(select_199709) <- colnames(scw)#给列表中每一个数赋予细胞编号
select_199709[pri_199709$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
select_199709[pri_199709$Scissor_neg] <- -1#被选为Scissor-的细胞赋值为2
sck <- subset(scw,GSE=="GSE159677")
scw <- AddMetaData(scw, metadata = select_199709, col.name = "Sci_199709")

table(sck$celltype,sck$Sci_199709)

DimPlot(scw,group.by = "Sci_199709",split.by = "GSE",
        cols = c("#0000ff","#f0f0f0","#ff0000"),
        reduction = "scphere")

DotPlot(scw,features = "PARP2")
load("GSE199709_gene.RData")


plot_cells(cds,
           genes=TOPgenes_burden_high,
           show_trajectory_graph=FALSE)



## GSE120521

phe_120521 <- c(rep(1,4),rep(0,4))
tag_120521 <- c("Stable","Unstable")

sck <- subset(scw,GSE=="GSE159677")
sck <- GetAssayData(sck,assay = "RNA",slot = "count")
# gc()

# sctPA <- subset(scw, GSE=="GSE159677")
# rm(sctPA)
# sctPA <- GetAssayData(sctPA,assay = "RNA",slot = "count")
# 
# rm(sck)
# class(GSE120521)

pa_1205211 <- Scissor_SCT(GSE120521, sck, phe_120521,alpha = 0.05, tag = tag_120521,
                           family = "binomial", #二分类,
                           Save_file = './Scissor/Sci_120521_PA.RData')


# 
# 
# pt3_120521 <- Scissor_SCT(GSE120521, scw, phe_120521, alpha = NULL, cutoff = 0.1, 
#                           family = "binomial", Load_file = './Scissor/Sci_120521.RData')
gc()
sel_pa_120521 <- rep("NR", ncol(scw))#创建一个列表，用来表示scRNA中细胞，
names(sel_pa_120521) <- colnames(scw)#给列表中每一个数赋予细胞编号

sel_pa_120521[pa_1205211$Scissor_pos] <- "Pos"#被选为Scissor+的细胞赋值为1
sel_pa_120521[pa_1205211$Scissor_neg] <- "Neg"#被选为Scissor-的细胞赋值为2


scw <- AddMetaData(scw, metadata = sel_pa_120521, col.name = "Sci_120521")

table(scw$celltype,scw$Sci_120521)
# rm(sck)
sct <- subset(scw,GSE=="GSE159677")

DimPlot(sct,group.by = "Sci_120521",
        cols = c("#0000ff","#f0f0f0","#ff0000"),
        reduction = "scphere")
ggsave(sci120521,filename = "Figures/Fig.3B.pdf",height = 6,width = 8)

df_120521 <- table(sct$GSE,sct$celltype,sct$state,sct$Sci_120521)
data <- reshape2::melt(df_120521)
# data$Var3 <- paste0("Clu_",data$Var3)
head(data)
data <- data[which(data$Var4!="NR"),]
data <- data[which(data$Var2!="Lymphatic EC"),]
data <- data[which(data$Var1=="GSE159677"),]

library(ggforce)
library(ggalluvial)
library(ggpubr)
data <- gather_set_data(data, 1:4)

sankey_120521_R <- ggplot(data,
       aes(axis1 = Var3, axis2 = Var2, axis3 = Var4,
           y= value)) +
  scale_x_discrete(limits = c("Sample","Cell.type", "Scissor"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Var3)) +
  scale_fill_manual(values=c("#0000ff","#ff0000"))+
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() +
  ggtitle("GSE120521",
          "Scissor for Stable(-) or Unstable(+)")+NoLegend()
ggsave(sankey_120521_R,filename = "Figures/Fig.3C1.pdf",height = 6,width = 8)

sankey_120521_L <- ggplot(data,
                        aes(axis1 = Var3, axis2 = Var2, axis3 = Var4,
                            y= value)) +
  scale_x_discrete(limits = c("Sample","Cell.type", "Scissor"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Var2)) + 
  scale_fill_manual(values=scales::hue_pal()(4))+
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() +
  ggtitle("GSE120521",
          "Scissor for Stable(-) or Unstable(+)")+NoLegend()
ggsave(sankey_120521_L,filename = "Figures/Fig.3C2.pdf",height = 6,width = 8)


load('./Scissor/Sci_120521_PA.RData')
numbers <- length(pa_1205211$Scissor_pos) + length(pa_1205211$Scissor_neg)
result1 <- reliability.test(X, Y, network, alpha = 0.05, family = "binomial", 
                            cell_num = numbers, n = 10, nfold = 10)

evaluate_summary <- evaluate.cell('./Scissor/Sci_120521_PA.RData', pa_1205211, 
                                  FDR = 0.05, bootstrap_n = 100)





## GSE132651
load("E:/works/artery ECs/GSE132651.RData")
GSE132651 <- scale(GSE132651)

GSE132651 <- as.matrix(GSE132651)
phe_132651 <- c(rep(0,6),rep(1,13))
tag_132651 <- c("NL","ABNL")
pri_132651 <- Scissor_SCT(GSE132651, sct, phe_132651,cutoff = 0.3, tag = tag_132651,
                          family = "binomial", #二分类
                          Load_file = './Scissor/Sci_132651.RData',
                          Save_file = './Scissor/Sci_132651.RData')
# 
# pt3_132651 <- Scissor_SCT(GSE132651, scw, phe_132651, alpha = NULL, cutoff = 0.1, 
#                           family = "binomial", Load_file = './Scissor/Sci_132651.RData')

select_132651 <- rep("NR", ncol(sct))#创建一个列表，用来表示scRNA中细胞，
names(select_132651) <- colnames(sct)#给列表中每一个数赋予细胞编号
select_132651[pri_132651$Scissor_pos] <- "Pos"#被选为Scissor+的细胞赋值为1
select_132651[pri_132651$Scissor_neg] <- "Neg"#被选为Scissor-的细胞赋值为2
sct <- AddMetaData(sct, metadata = select_132651, col.name = "Sci_132651")

FeaturePlot(scw,"UBE2G2",reduction = "scphere")
table(scw$celltype,scw$Sci_132651)
rm(sck)
sck <- subset(scw,GSE=="GSE159677")
FeaturePlot(sck,"SLC16A3",reduction = "scphere")
VlnPlot(sck,features = "LRRFIP1")
DimPlot(scw,group.by = "Sci_132651",split.by = "GSE",
        cols = c("#f0f0f0","#ff0000","#0000ff"),
        reduction = "scphere")


## GSE211662

DefaultAssay(scw) <- "RNA"
phe_211662 <- c(rep(0,3),rep(1,3))
tag_211662 <- c("LSS","HSS")
pri_211662 <- Scissor_SCT(GSE211662, scw, phe_211662,alpha = 0.05, tag = tag_211662,
                          family = "binomial", #二分类
                          Load_file = './Scissor/Sci_211662.RData')
#                         Save_file = './Scissor/Sci_211662.RData')
# 
# pt3_211662 <- Scissor_SCT(GSE211662, scw, phe_211662, alpha = NULL, cutoff = 0.1, 
#                           family = "binomial", Load_file = './Scissor/Sci_211662.RData')

select_211662 <- rep("NR", ncol(scw))#创建一个列表，用来表示scRNA中细胞，
names(select_211662) <- colnames(scw)#给列表中每一个数赋予细胞编号
select_211662[pri_211662$Scissor_pos] <- "Pos"#被选为Scissor+的细胞赋值为1
select_211662[pri_211662$Scissor_neg] <- "Neg"#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = select_211662, col.name = "Sci_211662")

FeaturePlot(scw,"UBE2G2",reduction = "scphere")
table(scw$celltype,scw$Sci_211662)
rm(sck)
sck <- subset(scw,GSE=="GSE159677")
FeaturePlot(sck,"SLC16A3",reduction = "scphere")
VlnPlot(sck,features = "LRRFIP1")
sci211662_plot <- DimPlot(sck,group.by = "Sci_211662",#split.by = "GSE",
        cols = c("#0000ff","#f0f0f0","#ff0000"),
        reduction = "scphere")
ggsave(sci211662_plot,filename = "Figures/Fig.3D.pdf",height = 6,width = 8)

# 
# library(reshape2)
# library(ggforce)
# library(tidyr)
# library(ggalluvial)


df_211662 <- table(scw$GSE,scw$celltype,scw$state,scw$Sci_211662)
data <- reshape2::melt(df_211662)
# data$Var3 <- paste0("Clu_",data$Var3)
head(data)
data <- data[which(data$Var4!="NR"),]
data <- data[which(data$Var2!="Lymphatic EC"),]
data <- data[which(data$Var1=="GSE159677"),]

library(ggforce)
library(ggalluvial)
library(ggpubr)
data <- gather_set_data(data, 1:4)

sankey_211662_R <- ggplot(data,
                          aes(axis1 = Var3, axis2 = Var2, axis3 = Var4,
                              y= value)) +
  scale_x_discrete(limits = c("Sample","Cell.type", "Scissor"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Var3)) +
  scale_fill_manual(values=c("#0000ff","#ff0000"))+
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() +
  ggtitle("GSE211662",
          "Scissor for Stable(-) or Unstable(+)")+NoLegend()
ggsave(sankey_211662_R,filename = "Figures/Fig.3C1.pdf",height = 6,width = 8)

sankey_211662_L <- ggplot(data,
                          aes(axis1 = Var3, axis2 = Var2, axis3 = Var4,
                              y= value)) +
  scale_x_discrete(limits = c("Sample","Cell.type", "Scissor"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Var2)) + 
  scale_fill_manual(values=scales::hue_pal()(3))+
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() +
  ggtitle("GSE211662",
          "Scissor for Stable(-) or Unstable(+)")+NoLegend()
ggsave(sankey_211662_L,filename = "Figures/Fig.3C2.pdf",height = 6,width = 8)

## GSE173669
GSE173669 <- read.csv("GSE173669/GSE173669.csv",row.names = 1)
names(GSE173669) <- c(paste0("MaNl_",1:2),paste0("HgOl_",1:2))
GSE173669 <- as.matrix(GSE173669)
DefaultAssay(scw) <- "RNA"
phe_173669 <- c(rep(0,2),rep(1,2))
tag_173669 <- c("MaNl","HgOl")
pri_173669 <- Scissor_SCT(GSE173669, scw, phe_173669,cutoff = 0.3,alpha = 0.2, tag = tag_173669,
                          family = "binomial", #二分类
                          Load_file = './Scissor/Sci_173669.RData')
# 
# pt3_173669 <- Scissor_SCT(GSE173669, scw, phe_173669, alpha = NULL, cutoff = 0.1, 
#                           family = "binomial", Load_file = './Scissor/Sci_173669.RData')

select_173669 <- rep("NR", ncol(scw))#创建一个列表，用来表示scRNA中细胞，
names(select_173669) <- colnames(scw)#给列表中每一个数赋予细胞编号
select_173669[pri_173669$Scissor_pos] <- "Pos"#被选为Scissor+的细胞赋值为1
select_173669[pri_173669$Scissor_neg] <- "Neg"#被选为Scissor-的细胞赋值为2
scw <- AddMetaData(scw, metadata = select_173669, col.name = "Sci_173669")

# FeaturePlot(scw,"UBE2G2",reduction = "scphere")
# table(scw$celltype,scw$Sci_173669)
# rm(sck)
sck <- subset(scw,GSE=="GSE159677")
# FeaturePlot(sck,"SLC16A3",reduction = "scphere")
# VlnPlot(sck,features = "LRRFIP1")
DimPlot(scw,group.by = "Sci_173669",#split.by = "GSE",
                    cols = c("#0000ff","#f0f0f0","#ff0000"),
                    reduction = "scphere")
VlnPlot(sck,group.by = "Sci_173669",features = "SCD")


## GSE207919

GSE207919 <- read.csv("GSE207919/GSE207919.txt",sep = "\t")
rownames(GSE207919) <- GSE207919$Gene.Symbol
GSE207919 <- GSE207919[,c("CXA_137_056_S11","CXB_125_068_S12","CXC_113_080_S13","CTA_101_092_S14","CTB_184_009_S15","CTC_172_021_S16")]
names(GSE207919) <- c(paste0("T0h_",1:3),paste0("T16h_",1:3))
GSE207919 <- as.matrix(GSE207919)

## GSE206927
GSE206927 <- read.csv("GSE206927/GSE206927.txt",sep = "\t",row.names = 1)
GSE206927 <- as.matrix(GSE206927)



## GSE118446

GSE118446 <- read.csv("GSE118446/GSE118446.txt",sep = "\t")
GSE118446 <- GSE118446[!duplicated(GSE118446$name),]
rownames(GSE118446) <- GSE118446$name
GSE118446 <- GSE118446[,8:25]
GSE118446 <- as.matrix(GSE118446)

FeaturePlot(scw,features = c("CLU","EEF1G","RGS5"),reduction = "scphere",min.cutoff = "q30")


## GSE184512
GSE184512 <- read.csv("GSE184512/GSE184512.txt",sep = "\t")
GSE184512 <- GSE184512[!duplicated(GSE184512$gene_symbol),]
rownames(GSE184512) <- GSE184512$gene_symbol
GSE184512 <- GSE184512[,9:12]
names(GSE184512) <- c(paste0("NC_",1:2),paste0("TNF_",1:2))
GSE184512 <- as.matrix(GSE184512)

## GSE201466
GSE201466 <- read.csv("GSE201466/GSE201466.txt",sep = "\t")
GSE201466_matrix <- read.csv("GSE201466/GSE201466_series_matrix.txt",sep = "")

GSE201466 <- GSE201466[!duplicated(GSE201466$Gene.Symbol),]
rownames(GSE201466) <- GSE201466$Gene.Symbol
GSE201466 <- GSE201466[,7:24]

names(GSE201466) <- stringr::str_split_fixed(names(GSE201466),pattern = "_",n=2)[,1]
names(GSE201466) <- stringr::str_remove(names(GSE201466),"X")
GSE201466 <- GSE201466[,c(paste0(1,c("A","B","C")),
                          paste0(2,c("A","B","C")),
                          paste0(3,c("A","B","C")))]
names(GSE201466) <- c(paste0(c("0h_"),1:3),paste0(c("4h_"),1:3),paste0(c("10h_"),1:3))
GSE201466 <- as.matrix(GSE201466)

##201522
GSE210522 <- as.matrix(read.table("GSE210522/raw.txt",header = T,row.names = 1))
# samples <- as.data.frame(names(GSE210522))
GSE210522_sel2 <- as.matrix(na.omit(GSE210522[,5:28]))
dim(GSE210522_sel)

##  Fig.4
Idents(scw) <- "celltype"
FAMs_type <- FindAllMarkers(scw,logfc.threshold = 0.4,slot = "data")
FAMs_type$diff.pct <- FAMs_type$pct.1-FAMs_type$pct.2
FAMs_type.sel <- FAMs_type %>% 
  group_by_("cluster") %>% 
  filter(diff.pct>0.15 & avg_log2FC>0.5 &  pct.1 >0.3) %>% 
  arrange("diff.pct") %>% 
  slice_max(order_by = diff.pct, n = 20)
geneset <- split(FAMs_type.sel$gene,FAMs_type.sel$cluster)
geneset_df <- data.frame(
   Gene.Symbol = unlist(geneset),
   Cell.Type = factor(c(rep("ITLN1+ EC",20),
                 rep("TCIM+ EC",20),
                 rep("ACKR1+ EC",20),
                 rep("Lymphatic+ EC",20)),levels = c("ITLN1+ EC",
                                                     "TCIM+ EC",
                                                     "ACKR1+ EC",
                                                     "Lymphatic+ EC"))
   )
save(geneset,file = "geneset.1224.RData")
load("geneset.1224.RData")

pdf( 
  file = "./Figures/Fig.4A.pdf", # 文件名称
  width = 8,           # 宽
  height = 4)              # 分辨率
# 2. 绘图
library(monocle3)

plot_cells(cds = cds,genes = geneset_df,show_trajectory_graph = F)+NoLegend()

# 3. 关闭画布
dev.off()



# library(GSVA)
# library(reshape2)
es1 <- gsva(GSE120521, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
# es2 <- gsva(GSE199709, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es3 <- gsva(GSE210522_sel2, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es4 <- gsva(GSE211662, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
# es5 <- gsva(GSE167024, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es6 <- gsva(GSE104140, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es7 <- gsva(GSE173669, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
# es8 <- gsva(GSE211662, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es9 <- gsva(GSE207919, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es10 <- gsva(GSE118446, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es11 <- gsva(GSE184512, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es12 <- gsva(GSE201466, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es13 <- gsva(GSE206927, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)


es1.data <- melt(t(as.data.frame(es1)))
es1.data$Var1 <- str_split_fixed(es1.data$Var1,pattern = "_",n=2)[,1]
names(es1.data) <- c("Treatment","Cluster","Value")

es2.data <- melt(t(as.data.frame(es2)))
es2.data$Var1 <- str_split_fixed(es2.data$Var1,pattern = "_",n=2)[,1]
names(es2.data) <- c("Treatment","Cluster","Value")

es3.data <- melt(t(as.data.frame(es3)))
es3.data$Var1 <- paste0(str_split_fixed(es3.data$Var1,pattern = "_",n=3)[,1],"_",
                        str_split_fixed(es3.data$Var1,pattern = "_",n=3)[,2])
names(es3.data) <- c("Treatment","Cluster","Value")

es4.data <- melt(t(as.data.frame(es4)))
es4.data$Var1 <- str_split_fixed(es4.data$Var1,pattern = "_",n=2)[,1]
names(es4.data) <- c("Treatment","Cluster","Value")

es5.data <- melt(t(as.data.frame(es5)))
es5.data$Var1 <- str_split_fixed(es5.data$Var1,pattern = "_",n=2)[,1]
names(es5.data) <- c("Treatment","Cluster","Value")
library(reshape2)
es6.data <- melt(t(as.data.frame(es6)))
es6.data$Var1 <- str_split_fixed(es6.data$Var1,pattern = "_",n=2)[,1]
names(es6.data) <- c("Treatment","Cluster","Value")
es6.data <- es6.data[-which(es6.data$Treatment %in% c("intimal xanthoma","pathologic intimal thickening")),]
write.csv(es6.data,"es6.data.csv")

es7.data <- melt(t(as.data.frame(es7)))
es7.data$Var1 <- str_split_fixed(es7.data$Var1,pattern = "_",n=2)[,1]
names(es7.data) <- c("Treatment","Cluster","Value")

es8.data <- melt(t(as.data.frame(es8)))
es8.data$Var1 <- str_split_fixed(es8.data$Var1,pattern = "_",n=2)[,1]
names(es8.data) <- c("Treatment","Cluster","Value")

es9.data <- melt(t(as.data.frame(es9)))
es9.data$Var1 <- str_split_fixed(es9.data$Var1,pattern = "_",n=2)[,1]
names(es9.data) <- c("Treatment","Cluster","Value")

es10.data <- melt(t(as.data.frame(es10)))
es10.data$Var1 <- str_split_fixed(es10.data$Var1,pattern = "_",n=2)[,1]
names(es10.data) <- c("Treatment","Cluster","Value")
levels(factor(es10.data$Treatment))
es10.data <- es10.data[-which(es10.data$Treatment %in% c("HUVEC.IL1b","HUVEC.TGFb2")),]


es11.data <- melt(t(as.data.frame(es11)))
es11.data$Var1 <- str_split_fixed(es11.data$Var1,pattern = "_",n=2)[,1]
names(es11.data) <- c("Treatment","Cluster","Value")

es12.data <- melt(t(as.data.frame(es12)))
es12.data$Var1 <- str_split_fixed(es12.data$Var1,pattern = "_",n=2)[,1]
names(es12.data) <- c("Treatment","Cluster","Value")
es12.data$Treatment <- factor(es12.data$Treatment,levels = c("0h","4h","10h"))

es13.data <- melt(t(as.data.frame(es13)))
es13.data$Var1 <- str_split_fixed(es13.data$Var1,pattern = "_",n=2)[,1]
names(es13.data) <- c("Treatment","Cluster","Value")


esTNF.data <- rbind(es9.data,es11.data)
esTNF.data <- rbind(esTNF.data,es12.data)
levels(factor(esTNF.data$Treatment))
esTNF.data[which(esTNF.data$Treatment=="0h"),1] <- "GSE201466_TNF0h"
esTNF.data[which(esTNF.data$Treatment=="4h"),1] <- "GSE201466_TNF4h"
esTNF.data[which(esTNF.data$Treatment=="10h"),1] <- "GSE201466_TNF10h"
esTNF.data[which(esTNF.data$Treatment=="NC"),1] <- "GSE184512_NC"
esTNF.data[which(esTNF.data$Treatment=="TNF"),1] <- "GSE184512_TNF"
esTNF.data[which(esTNF.data$Treatment=="T0h"),1] <- "GSE207919_NC"
esTNF.data[which(esTNF.data$Treatment=="T16h"),1] <- "GSE207919_TNF"

ggplot(es3.data, aes(x=Treatment, y=Value, fill = Treatment)) + 
  geom_boxplot()+
  # geom_violin()+
  # geom_jitter(stroke = 0.1,width = 0.05)+
  geom_point()+
  facet_wrap(~Cluster, scale = "free_x") +
  scale_y_continuous(name = "Clusters Score") +
  theme_bw()+ theme(legend.position="top")+ 
  # theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_x_discrete("Samples", labels = NULL)

# unique()
#绘制boxplot
library(ggpubr)
scaleFUN <- function(x) sprintf("%.2f", x) 

sss_104140 <- ggboxplot(es6.data, x = "Treatment", y = "Value",linetype = 3,add = "jitter",
         color = "Treatment", label.rectangle = T,combine = T,point.size = 2,
         line.color = "gray", line.size = 0.4,point.color="black",
         facet.by = "Cluster")+
  scale_x_discrete("Samples", labels = NULL)+ 
  scale_y_continuous(labels=scaleFUN) +
  # theme_pubr() +
  ggtitle("GSE104140")
ggsave(sss_104140,filename = "Figures/Fig.4B.pdf",height = 6,width = 8)

sss_TNF <- ggboxplot(esTNF.data, x = "Treatment", y = "Value",linetype = 3,add = "jitter",
                        color = "Treatment", label.rectangle = T,combine = T,point.size = 2,
                        line.color = "gray", line.size = 0.4,point.color="black",
                        facet.by = "Cluster")+
  scale_x_discrete("Samples", labels = NULL)+
  scale_y_continuous(labels=scaleFUN) +
  # theme_pubr() +
  ggtitle("GSE184512 & GSE207919 & GSE201466")
ggsave(sss_TNF,filename = "Figures/Fig.4C.pdf",height = 6,width = 8)
# 
# sss_201466 <- ggboxplot(es12.data, x = "Treatment", y = "Value",linetype = 3,add = "jitter",
#                         color = "Treatment", label.rectangle = T,combine = T,point.size = 2,
#                         line.color = "gray", line.size = 0.4,point.color="black",
#                         facet.by = "Cluster")+
#   scale_x_discrete("Samples", labels = NULL)+
#   scale_y_continuous(labels=scaleFUN) +
#   # theme_pubr() +
#   ggtitle("GSE201466")
# ggsave(sss_201466,filename = "Figures/Fig.4D.pdf",height = 4,width = 8)

sss_118446 <- ggboxplot(es10.data, x = "Treatment", y = "Value",linetype = 3,add = "jitter",
                        color = "Treatment", label.rectangle = T,combine = T,point.size = 2,
                        line.color = "gray", line.size = 0.4,point.color="black",
                        facet.by = "Cluster")+
  scale_x_discrete("Samples", labels = NULL)+
  scale_y_continuous(labels=scaleFUN) +
  # theme_pubr() +
  ggtitle("GSE118446")
ggsave(sss_118446,filename = "Figures/Fig.4D.pdf",height = 6,width = 8)

sss_210522 <- ggboxplot(es3.data, x = "Treatment", y = "Value",linetype = 3,add = "jitter",
                        color = "Treatment", label.rectangle = T,combine = T,point.size = 2,
                        line.color = "gray", line.size = 0.4,point.color="black",
                        facet.by = "Cluster")+
  scale_x_discrete("Samples", labels = NULL)+
  scale_y_continuous(labels=scaleFUN) +
  # theme_pubr() +
  ggtitle("GSE210522")
ggsave(sss_210522,filename = "Figures/Fig.4E.pdf",height = 6,width = 8)
Nebulosa::plot_density(scw,
                       features = c("TGFBR1","NOTCH2","BMP2","ANKRD1","OSGIN2"),
                       reduction = "scphere")


DimPlot(scw,group.by = "Sci_173669",cols = c("grey","red","blue"),reduction = "scphere")


## Fig.5

df_track_genes <- data.frame(
  Gene.Symbol = Track_genes_sig50
)

FAMs_type.sel <- FAMs_type %>% filter(avg_log2FC > 0)
df_track_genes$Clu <- FAMs_type.sel[match(df_track_genes$Gene.Symbol,FAMs_type.sel$gene),6]
scales::hue_pal()(4)
group

FeaturePlot(scw,"FAU",reduction = "scphere",min.cutoff = "q30")

Nebulosa::plot_density(scw,features = "FAU",reduction = "scphere")


# cluster profiler compare cluster


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
test <- as.data.frame(scw@assays[["RNA"]]@data@Dimnames[[1]])
gene_dict <- bitr(scw@assays[["RNA"]]@data@Dimnames[[1]],fromType = "SYMBOL",
                  toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)

FAMs_type.sel$ENTREZID <- gene_dict[match(FAMs_type.sel$gene,gene_dict$SYMBOL),2]

FAMs_type.sel_KKK <- FAMs_type.sel

FAMs_type.sel_KKK <- FAMs_type.sel_KKK[!FAMs_type.sel_KKK$gene %in% c("HLA-DRA","HLA-DPA1","HLA-DRB1"),]

merge.clu_list <- split(FAMs_type.sel_KKK$ENTREZID,FAMs_type.sel_KKK$cluster)
merge.clu_list <- na.omit(merge.clu_list)
names(merge.clu_list)
com.KEGG <- compareCluster(merge.clu_list,fun = "enrichKEGG",organism="hsa", pvalueCutoff=0.05)
com.KEGG <- setReadable(com.KEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
com.KEGG_plot <- dotplot(com.KEGG,label_format=40) +   theme_bw()
ggsave(com.KEGG_plot,filename = "Figures/Fig.5B.pdf",height = 8,width =6)

com.GO <- compareCluster(merge.clu_list,fun = "enrichGO",OrgDb = org.Hs.eg.db, pvalueCutoff=0.05)
com.GO <- simplify(com.GO)
com.GO <- pairwise_termsim(com.GO)
com.GO <- setReadable(com.GO,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
com.GO_plot <- dotplot(com.GO,label_format=40) +   theme_bw()
ggsave(com.GO_plot,filename = "Figures/Fig.5C.pdf",height = 8,width =6)


com.Pathway <- compareCluster(merge.clu_list,fun = "enrichPathway",pvalueCutoff=0.05)
com.Pathway@readable <- F
com.Pathway <- setReadable(com.Pathway,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
com.Pathway_plot <- dotplot(com.Pathway,label_format=40) +   theme_bw() 
ggsave(com.Pathway_plot,filename = "Figures/Fig.5D.pdf",height = 8,width =6)

install.packages('GOplot')

# library(GOplot)
library(GOplot)
data(EC) 

EC$david 

com_GO <- com.GO@compareClusterResult
com_GO <- com_GO[,c(1,2,3,9,7)]
names(com_GO) <- names(EC$david)
com_GO$Genes <- stringr::str_replace_all(com_GO$Genes,pattern = "/",replacement = ",")

com_KEGG <- com.KEGG@compareClusterResult
com_KEGG <- com_KEGG[,c(1,2,3,9,7)]
names(com_KEGG) <- names(EC$david)
com_KEGG$Genes <- stringr::str_replace_all(com_KEGG$Genes,pattern = "/",replacement = ",")

com_Pathway <- com.Pathway@compareClusterResult
com_Pathway <- com_Pathway[,c(1,2,3,9,7)]
names(com_Pathway) <- names(EC$david)
com_Pathway$Genes <- stringr::str_replace_all(com_Pathway$Genes,pattern = "/",replacement = ",")

list_gene <- FAM
list_gene <- list_gene[,c(7,2,3,4,1,5,8)]
names(list_gene) <- names(tmp_genelist)

go_circ <- circle_dat(com_GO, list_gene)
go_circ$category <- factor(go_circ$category,levels = group)
# go_circ$term <- str_to_title(go_circ$term)

kegg_circ <- circle_dat(com_KEGG, list_gene)
kegg_circ$category <- factor(kegg_circ$category,levels = group)
# kegg_circ$term <- str_to_title(kegg_circ$term)

pathway_circ <- circle_dat(com_Pathway, list_gene)
pathway_circ$category <- factor(pathway_circ$category,levels = group)
# pathway_circ$term <- str_to_title(pathway_circ$term)
# circ <- circ[circ$category!="Lymphatic EC",]


## GO

data <- go_circ
display <- "multiple"
title <- "GO Enrichment"
colour <- c(scales::hue_pal()(4))
cols <- colour
labels <- 3
rang <- c(1, 10)
data$adj_pval <- -log(data$adj_pval, 10)
sub <- data[!duplicated(paste0(data$category,"_",data$term)), ]
sub2 <- sub %>% 
  group_by_("category") %>% 
  arrange("zscore") %>% 
  slice_max(order_by = zscore, n = labels)
lmts.x <- c(range(sub$zscore)[1],range(sub$zscore)[2]+0.5)
lmts.y <- c(1,range(sub$adj_pval)[2]+0.5)
terms <- sub[duplicated(sub$term),3]
sub[which(sub$term %in% terms),3]

g <- ggplot(sub, aes(zscore, adj_pval, fill = category, 
                     size = count)) + 
  labs(title = title, x = "z-score",
       y = "-log (adj p-value)") + 
  geom_point(shape = 21, col = "black",
             alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.05), col = "orange") + 
  scale_size(range = rang, guide = "none")

GO.com_Plot <- g + facet_grid(. ~ category, space = "free_x", 
                                scales = "fixed") + 
  scale_fill_manual(values = cols, 
                    guide = "none") + 
  scale_x_continuous(limits = lmts.x)+
  scale_y_continuous(limits = lmts.y,labels=scaleFUN)+
  theme_pubr()+
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[1]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[1],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[2]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[2],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[3]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[3],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[4]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[4],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  # ggrepel::geom_text_repel(data = sub[which(sub$term %in% terms),], aes(x = zscore, y = adj_pval,
  #                                                                             label = str_wrap(term,width = 20)),
  #                           force=20,color="black",size=3,point.padding = 0.5,hjust = 0.5,
  #                           arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
  #                           segment.color="red",segment.size=0.2,segment.alpha=1,nudge_y=-0.5) +
  # ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval,
  # label = str_wrap(term,width = 20)), size = 3,max.iter = 1) + 
  theme(axis.title = element_text(size = 14, 
                                  face = "bold"), axis.text = element_text(size = 14), 
        axis.line = element_line(colour = "grey80"), 
        axis.ticks = element_line(colour = "grey80"), 
        panel.border = element_rect(fill = "transparent", 
                                    colour = "grey80"), panel.background = element_blank(), 
        panel.grid = element_blank(), plot.background = element_blank())
ggsave(GO.com_Plot,filename = "Figures/Fig.5C1.pdf",height = 4,width = 8)

sub3 <- sub %>% 
  group_by_("category") %>% 
  arrange("zscore") %>% 
  slice_max(order_by = zscore, n = 5)
com.GO_dot <- ggplot(sub3,aes(x = category,y = str_wrap(term,width = 40)))+
  geom_point(aes(color = adj_pval,
                 size = count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("")+ylab("GO_Enrichment")+
  theme_bw()
ggsave(com.GO_dot,filename = "Figures/Fig.5C2.pdf",height = 6,width = 8)


## KEGG

data <- kegg_circ
display <- "multiple"
title <- "KEGG Enrichment"
colour <- c(scales::hue_pal()(4))
cols <- colour
labels <- 3
rang <- c(1, 10)
data$adj_pval <- -log(data$adj_pval, 10)
sub <- data[!duplicated(paste0(data$category,"_",data$term)), ]
sub2 <- sub %>% 
  group_by_("category") %>% 
  arrange("zscore") %>% 
  slice_max(order_by = zscore, n = 3)
lmts.x <- c(range(sub$zscore)[1],range(sub$zscore)[2]+0.5)
lmts.y <- c(1,range(sub$adj_pval)[2]+0.5)
terms <- sub[duplicated(sub$term),3]
sub[which(sub$term %in% terms),3]

g <- ggplot(sub, aes(zscore, adj_pval, fill = category, 
                     size = count)) + 
  labs(title = title, x = "z-score",
       y = "-log (adj p-value)") + 
  geom_point(shape = 21, col = "black",
             alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.05), col = "orange") + 
  scale_size(range = rang, guide = "none")

KEGG.com_Plot <- g + facet_grid(. ~ category, space = "free_x", 
                              scales = "fixed") + 
  scale_fill_manual(values = cols, 
                    guide = "none") + 
  scale_x_continuous(limits = lmts.x)+
  scale_y_continuous(limits = lmts.y,labels=scaleFUN)+
  theme_pubr()+
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[1]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[1],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[2]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[2],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[3]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[3],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[4]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[4],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  # ggrepel::geom_text_repel(data = sub[which(sub$term %in% terms),], aes(x = zscore, y = adj_pval,
  #                                                                             label = str_wrap(term,width = 20)),
  #                           force=20,color="black",size=3,point.padding = 0.5,hjust = 0.5,
  #                           arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
  #                           segment.color="red",segment.size=0.2,segment.alpha=1,nudge_y=-0.5) +
  # ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval,
  # label = str_wrap(term,width = 20)), size = 3,max.iter = 1) + 
  theme(axis.title = element_text(size = 14, 
                                  face = "bold"), axis.text = element_text(size = 14), 
        axis.line = element_line(colour = "grey80"), 
        axis.ticks = element_line(colour = "grey80"), 
        panel.border = element_rect(fill = "transparent", 
                                    colour = "grey80"), panel.background = element_blank(), 
        panel.grid = element_blank(), plot.background = element_blank())
ggsave(KEGG.com_Plot,filename = "Figures/Fig.5B1.pdf",height = 4,width = 8)

sub3 <- sub %>% 
  group_by_("category") %>% 
  arrange("zscore") %>% 
  slice_max(order_by = zscore, n = 5)
com.KEGG_dot <- ggplot(sub3,aes(x = category,y = str_wrap(term,width = 40)))+
  geom_point(aes(color = adj_pval,
                 size = count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("")+ylab("KEGG_Enrichment")+
  theme_bw()
ggsave(com.KEGG_dot,filename = "Figures/Fig.5B2.pdf",height = 6,width = 8)

## Pathway
data <- pathway_circ
display <- "multiple"
title <- "PATHWAY Enrichment"
colour <- c(scales::hue_pal()(4))
cols <- colour
labels <- 3
rang <- c(1, 10)
data$adj_pval <- -log(data$adj_pval, 10)
sub <- data[!duplicated(paste0(data$category,"_",data$term)), ]
sub2 <- sub %>% 
  group_by_("category") %>% 
  arrange("zscore") %>% 
  slice_max(order_by = zscore, n = 3)
lmts.x <- c(range(sub$zscore)[1],range(sub$zscore)[2]+0.5)
lmts.y <- c(1,range(sub$adj_pval)[2]+0.5)
terms <- sub[duplicated(sub$term),3]
sub[which(sub$term %in% terms),3]

g <- ggplot(sub, aes(zscore, adj_pval, fill = category, 
                     size = count)) + 
  labs(title = title, x = "z-score",
       y = "-log (adj p-value)") + 
  geom_point(shape = 21, col = "black",
             alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.05), col = "orange") + 
  scale_size(range = rang, guide = "none")

PATHWAY.com_Plot <- g + facet_grid(. ~ category, space = "free_x", 
                                scales = "fixed") + 
  scale_fill_manual(values = cols, 
                    guide = "none") + 
  scale_x_continuous(limits = lmts.x)+
  scale_y_continuous(limits = lmts.y,labels=scaleFUN)+
  theme_pubr()+
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[1]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[1],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[2]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[2],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[3]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[3],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  ggrepel::geom_label_repel(data = sub2[which(sub2$category==group[4]),], aes(x = zscore, y = adj_pval,
                                                                              label = str_wrap(term,width = 20)),
                            force=20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,fill = cols[4],alpha = 0.8, 
                            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                            segment.color="grey20",segment.size=0.2,segment.alpha=1,nudge_y=0.5) +
  # ggrepel::geom_text_repel(data = sub[which(sub$term %in% terms),], aes(x = zscore, y = adj_pval,
  #                                                                             label = str_wrap(term,width = 20)),
  #                           force=20,color="black",size=3,point.padding = 0.5,hjust = 0.5,
  #                           arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
  #                           segment.color="red",segment.size=0.2,segment.alpha=1,nudge_y=-0.5) +
  # ggrepel::geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval,
  # label = str_wrap(term,width = 20)), size = 3,max.iter = 1) + 
  theme(axis.title = element_text(size = 14, 
                                  face = "bold"), axis.text = element_text(size = 14), 
        axis.line = element_line(colour = "grey80"), 
        axis.ticks = element_line(colour = "grey80"), 
        panel.border = element_rect(fill = "transparent", 
                                    colour = "grey80"), panel.background = element_blank(), 
        panel.grid = element_blank(), plot.background = element_blank())
ggsave(PATHWAY.com_Plot,filename = "Figures/Fig.5D1.pdf",height = 4,width = 8)

sub3 <- sub %>% 
  group_by_("category") %>% 
  arrange("zscore") %>% 
  slice_max(order_by = zscore, n = 5)
com.PATHWAY_dot <- ggplot(sub3,aes(x = category,y = str_wrap(term,width = 40)))+
  geom_point(aes(color = adj_pval,
                 size = count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("")+ylab("PATHWAY_Enrichment")+
  theme_bw()
ggsave(com.PATHWAY_dot,filename = "Figures/Fig.5D2.pdf",height = 6,width = 8)


RidgePlot(scw,"PTPRC")

Nebulosa::plot_density(scw,"ITLN1",reduction = "scphere")
VlnPlot(scw,"ITLN1",group.by = "NEW.clu")
## Fig.6 
# load("cds1222.RData")




# 
# library(monocle3,lib.loc = "D:/Program Files/R/R-4.2.1/library")
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
## NEW 
Track_genes2 <- graph_test(cds, neighbor_graph="principal_graph")


save(Track_genes2,file = "trackgenes.new.RData")
save(cds,file = "cds.NEW.RData")


# Track_genes.fli1 <- Track_genes2[,c(5,2,3,4,1,6)] %>% 
#   filter(q_value < 1e-3) %>% 
#   arrange(desc(morans_I))
Track_genes.fli <- Track_genes[,c(5,2,3,4,1,6)] %>% 
  filter(q_value < 1e-3) %>% 
  arrange(desc(morans_I))
genelist <- pull(Track_genes.fli, gene_short_name) %>% as.character()
set.seed(1116)


gene_module <- find_gene_modules(cds[genelist,], resolution=5e-3, cores = 1)
# save(gene_module,file = "gene_module1223.RData")
load("gene_module1223.RData")
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


### by 3d clusters (NEW.clu) 
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


## TFs prediction
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

module.23 <- unlist(c(gene_module[which(gene_module$module == 23),1]))
module.23.con <- bitr(module.23,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
module.24 <- unlist(c(gene_module[which(gene_module$module == 24),1]))
module.24.con <- bitr(module.24,fromType = "SYMBOL",
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
module.35 <- unlist(c(gene_module[which(gene_module$module == 35),1]))
module.35.con <- bitr(module.35,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)

module17.TFs.RRA <- ChEA.TFs.RRA(module.17.con$SYMBOL)
module24.TFs.RRA <- ChEA.TFs.RRA(module.24.con$SYMBOL)
module23.TFs.RRA <- ChEA.TFs.RRA(module.23.con$SYMBOL)
module25.TFs.RRA <- ChEA.TFs.RRA(module.25.con$SYMBOL)
module34.TFs.RRA <- ChEA.TFs.RRA(module.34.con$SYMBOL)


FAMs_ITLN1 <- FAMs_type %>% filter(cluster == "ITLN1+ EC") %>% pull( gene) %>% as.character()
FAMs_ITLN1.con <- bitr(FAMs_ITLN1,fromType = "SYMBOL",
                       toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
FAMs_TCIM <- FAMs_type %>% filter(cluster == "TCIM+ EC") %>% pull( gene) %>% as.character()
FAMs_TCIM.con <- bitr(FAMs_TCIM,fromType = "SYMBOL",
                      toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
FAMs_ACKR1 <- FAMs_type %>% filter(cluster == "ACKR1+ EC") %>% pull( gene) %>% as.character()
FAMs_ACKR1.con <- bitr(FAMs_ACKR1,fromType = "SYMBOL",
                       toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db,drop = F)
FAMs_Lymphatic <- FAMs_type %>% filter(cluster == "Lymphatic EC") %>% pull( gene) %>% as.character()
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
upset(fromList(TFs),keep.order = T,nsets = 8,order.by = "freq")
intersect(TFs$`Module 23`,TFs$`ITLN1+ EC`)
intersect(TFs$`Module 17`,TFs$`TCIM+ EC`)
intersect(TFs$`Module 24`,TFs$`TCIM+ EC`)
intersect(TFs$`Module 34`,TFs$`ACKR1+ EC`)
intersect(TFs$`Module 25`,TFs$`Lymphatic EC`)



intersect(module.23.con$SYMBOL,FAMs_ITLN1.con$SYMBOL) #85
# intersect(module.17.con$SYMBOL,FAMs_TCIM.con$SYMBOL) #28
# intersect(module.24.con$SYMBOL,FAMs_TCIM.con$SYMBOL) #19
intersect(union(module.24.con$SYMBOL,module.17.con$SYMBOL),FAMs_TCIM.con$SYMBOL) #47

intersect(module.34.con$SYMBOL,FAMs_ACKR1.con$SYMBOL) #58
intersect(module.25.con$SYMBOL,FAMs_Lymphatic.con$SYMBOL) #70


ITLN1.TFs.RRA <- ChEA.TFs.RRA(intersect(module.23.con$SYMBOL,FAMs_ITLN1.con$SYMBOL))
TCIM.TFs.RRA <- ChEA.TFs.RRA(intersect(union(module.24.con$SYMBOL,module.17.con$SYMBOL),FAMs_TCIM.con$SYMBOL))
ACKR1.TFs.RRA <- ChEA.TFs.RRA(intersect(module.34.con$SYMBOL,FAMs_ACKR1.con$SYMBOL))
Lymphatic.TFs.RRA <- ChEA.TFs.RRA(intersect(module.25.con$SYMBOL,FAMs_Lymphatic.con$SYMBOL))


# TFs.type <- list(
#   rownames(ITLN1.TFs.RRA)[1:20],
#   rownames(TCIM.TFs.RRA)[1:20],
#   rownames(ACKR1.TFs.RRA)[1:20],
#   rownames(Lymphatic.TFs.RRA)[1:20]
# )
# names(TFs.type) <- c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC")
# upset(fromList(TFs.type),keep.order = T,nsets = 4,order.by = "freq")



expr_data <- AverageExpression(scw,assays = "RNA")
expr_data <- expr_data[[1]]

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
ITLN1_EC_TFs$Uniques <- -log10(ITLN1_EC_TFs$Uniques)
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
TCIM_EC_TFs$Uniques <- -log10(TCIM_EC_TFs$Uniques)
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
ACKR1_EC_TFs$Uniques <- -log10(ACKR1_EC_TFs$Uniques)
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



tmp1 <- melt(ITLN1_EC_TFs[,1:3])

tmp2 <- melt(TCIM_EC_TFs[,1:3])

tmp3 <- melt(ACKR1_EC_TFs[,1:3])

tmp4 <- melt(Lymphatic_EC_TFs[,1:3])

tmp <- rbind(tmp1,tmp2)
tmp <- rbind(tmp,tmp3)
tmp <- rbind(tmp,tmp4)

expr_tf <- expr_data[as.character(unique(tmp$TFs)),]
df_tf_expr2 <- as.data.frame(expr_tf)
# df_tf_expr2 <- df_tf_expr2[rownames(df_tf_expr),]
for (i in 1:nrow(df_tf_expr2)) {
  # gene_module.df[i,7] <- max(gene_module.df[i,3:6])
  df_tf_expr2[i,5] <- names(which.max(df_tf_expr2[i,1:4]))
}
df_tf_expr <- df_tf_expr2
df_tf_expr2$TFs <- rownames(df_tf_expr2)
df_tf_expr2 <- split(df_tf_expr2$TFs,df_tf_expr2$V5)

tmp1 <- tmp1[which(tmp1$TFs %in% df_tf_expr2$`ITLN1+ EC`),]
tmp2 <- tmp2[which(tmp2$TFs %in% df_tf_expr2$`TCIM+ EC`),]
tmp3 <- tmp3[which(tmp3$TFs %in% df_tf_expr2$`ACKR1+ EC`),]
tmp4 <- tmp4[which(tmp4$TFs %in% df_tf_expr2$`Lymphatic EC`),]
tmp <- rbind(tmp1,tmp2)
tmp <- rbind(tmp,tmp3)
tmp <- rbind(tmp,tmp4)

df_tf_expr$V5 <- factor(df_tf_expr$V5,levels = group)
df_tf_expr <- df_tf_expr[order(df_tf_expr$V5),]
df_tf_expr <- df_tf_expr[as.character(unique(tmp$TFs)),1:4]
df_tf_expr <- as.matrix(df_tf_expr)
df_tf_expr <- t(scale(t(df_tf_expr)))


tmp$TFs <- factor(tmp$TFs,levels = rownames(df_tf_expr))
tmp$value <- ifelse(tmp$value<0,-tmp$value,tmp$value)
ggplot(tmp, aes(x = TFs, y = value, fill = variable)) +
  geom_bar(stat="identity", position = "dodge",
           color = 'white', alpha = 0.8, width = 0.95)+#coord_polar(start = 0)+
  # scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0),limits = c(-20,14)) +
  ggpubr::theme_pubr()+coord_polar(start = 0)

library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
range(df_tf_expr)
mycol <- colorRamp2(c(-1.5, 0, 1.5),c("#0da9ce", "white", "#e74a32"))

circos.par(gap.after = c(0)) #调整圆环首位间距；
circos.heatmap(df_tf_expr,
               col = mycol,
               cluster = F, #是否对行聚类
               # dend.side = "inside", #聚类树方向：inside显示在圆环内圈，inside为圆环外圈；
               rownames.side ="outside", #行名方向；
               rownames.col = "black",
               track.height = 0.6, #轨道高度，即圆环/格子的粗细
               rownames.cex = 1)
lg <- Legend(title = "Express",
             col_fun = mycol,
             direction = c("vertical"),
             title_position = c('topcenter'))
draw(lg, x = unit(0.9, "npc"), y = unit(0.5, "npc"), just = c("right", "center"))
circos.clear()


### Fig.7 scFEA
#scFEA

GSEs <- levels(factor(scw$GSE))

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

Idents(scw) <- "celltype"
BAL_tmp <- AverageExpression(scw,assays = "BALANCE",slot = "count")[[1]]
Idents(scw) <- "NEW.clu"
BAL_tmp2 <- AverageExpression(scw,assays = "BALANCE",slot = "count")[[1]]

Metabolites <- c("Glucose","G6P","G3P","X3PD","Pyruvate","Acetyl.CoA","Oxaloacetate",
                 "Citrate","X2OG","Succinyl.CoA","Succinate","Fumarate","Malate",
                 "Lactate")
colnames(balance)
BAL_TCA <- BAL_tmp[Metabolites,]
BAL_TCA2 <- BAL_tmp2[Metabolites,]
pmea <- ComplexHeatmap::pheatmap(BAL_TCA, scale="row", clustering_method="ward.D2",
                                 cellwidth = 20, cellheight = 15,border_color = "white",
                                 #annotation_row = row_anno,#annotation_col = ann_colors,
                                 cluster_cols = F,cluster_rows = F,show_rownames = T,use_raster =F)
pmea <- ggplotify::as.ggplot(pmea)
pmea
ggsave(pmea,filename = "Figures/Fig.7A.pdf",height = 6,width = 8)


col_anno <- melt(table(scw$celltype,scw$NEW.clu))
col_anno <- col_anno[which(col_anno$value!=0),]
rownames(col_anno) <- col_anno$Var2
col_anno <- as.data.frame(col_anno[,1])
names(col_anno) <- "Cell.Type"
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
pmea2 <- ComplexHeatmap::pheatmap(BAL_TCA2, scale="row", clustering_method="ward.D2",
                                  cellwidth = 20, cellheight = 15,border_color = "white",
                                  annotation_col = col_anno,annotation_colors = ann_colors,
                                  cluster_cols = F,cluster_rows = F,show_rownames = T,use_raster =F)

manual_order = c(3,5,7,6,9,4,1,10,2,8)
pmea2@column_order <- manual_order
pmea2 <- ggplotify::as.ggplot(pmea2)
ggsave(pmea2,filename = "Figures/sFig.5A.pdf",height = 6,width = 8)


Idents(scw) <- "celltype"
DefaultAssay(scw) <- "FLUX"
FeatureScatter(scw,feature1 = "M-5",feature2 = "M-6",shuffle = T,pt.size = 1)
VlnPlot(scw,features = c("M-5","M-6"))

FLU_tmp <- AverageExpression(scw,assays = "FLUX",slot = "count")[[1]]



FLU_all <- data.frame(
  FLU_all = flux[,"M_6"]/flux[,"M_5"],
  scw$celltype
)

box_glyco <- ggplot(FLU_all, aes(x=scw.celltype, y=FLU_all, fill = scw.celltype)) + 
  geom_boxplot()+
  # geom_violin()+
  # geom_jitter(stroke = 0.01,width = 0.2,size=3)+
  # geom_point(size=0.4)+
  # # geom_line(data=filter(mydata,treatment=="v"),aes(group==ID))+
  # facet_wrap(~Celltype, scale = "free_x") +
  scale_y_continuous(name = "Glycolysis Score") +
  scale_x_discrete(name = "") +
  theme_bw()+NoLegend()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
ggsave(box_glyco,filename = "Figures/Fig.7B.pdf",height = 6,width = 8)



BAL_all2 <- data.frame(
  scw$celltype,
  balance[,"Fatty.Acid"],
  balance[,"X.E.E..Farnesyl.PP"],
  balance[,"Cholesterol"]
)
BAL_all2 <- melt(BAL_all2)

levels(factor(BAL_all2$variable))

BAL_all2$variable <- case_when(
  BAL_all2$variable == "balance....Fatty.Acid.." ~ "Fatty Acid",
  BAL_all2$variable == "balance....X.E.E..Farnesyl.PP.." ~ "(E,E)-Farnesyl-PP",
  BAL_all2$variable == "balance....Cholesterol.." ~ "Cholesterol")

box_FFA <- ggplot(BAL_all2, aes(x=scw.celltype, y=value, fill = scw.celltype)) + 
  geom_boxplot()+
  # geom_violin()+
  # geom_jitter(stroke = 0.01,width = 0.2,size=0.4)+
  # geom_point(size=0.4)+
  # # geom_line(data=filter(mydata,treatment=="v"),aes(group==ID))+
  facet_wrap(~variable, scale = "free_x") +
  scale_y_continuous(name = "Lipid Balance") +
  scale_x_discrete(name = "") +
  theme_bw()+NoLegend()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))

box_FFA
ggsave(box_FFA,filename = "Figures/Fig.7C.pdf",height = 6,width = 8)
















##########################################
# library(UpSetR)
# library(ggupset)
act <- act[c(paste0("FAMs_",c(1,2,4,6,9,10,5,7,8,3)))]
df.tmp <- melt(fromList(act)) 
tidy_movies

UpSetR::upset(fromList(act),nsets = 10,nintersects = 50,
      set_size.show = T,order.by = "freq",decreasing = F,
      group.by = "degree",keep.order = T)


genefeature <- split(FAMs_3d_scale_fli$gene,FAMs_3d_scale_fli$cluster)
genefeature_fli <- lapply(genefeature, function(x){return(x[c(1:5,(length(x)-4):length(x))])})

table(FAMs_3d_scale_fli$cluster)
table(FAMs_3d_scale_NEW$cluster)

upset(fromList(genefeature_fli),nsets = 10,nintersects = 50,
      set_size.show = T,order.by = "freq",decreasing = F)
DimPlot(scw,group.by = "NEW.clu",label = T,label.box = T,reduction = "scphere")

# 1,4,6,9,10
# 3





# FAMs_3d_roc <- FindAllMarkers(scw,slot = "scale.data",assay = "SCT",test.use = "roc")
# names(FAMs_3d_scale) <- names(FAMs_3d)
# names(FAMs_3d_roc) <- names(FAMs_3d)
# 
# # FAMs_3d_scale <- FAMs_3d
# FAMs_sel2 <- FAMs_3d %>%
#   group_by_("cluster") %>%
#   filter(avg_diff>0.2 & pct.2<0.5) %>%
#   arrange("diff.pct") %>%
#   slice_max(order_by = avg_diff, n = 20)
# table(FAMs_sel2$cluster)

# vocalno plot

FAMs_3d_scale.bak <- FAMs_3d_scale
save(FAMs_3d_scale_NEW,file = "FAMs_3d_scale_NEW.RData")
names(FAMs_3d_scale)[2] <- "avg_log2FC"
vocaplot <- markerVocalno(markers = FAMs_3d_scale_NEW,
                          labelCol = scales::hue_pal()(10))
vocaplot
FAMs_3d_scale.bak <- FAMs_3d_scale
names(FAMs_3d)
names(FAMs_3d_scale)

ggsave(vocaplot,filename = "Figures/Fig.1E(scale).pdf",height = 6,width = 8)
ggsave(vocaplot,filename = "Figures/Fig.1E(data).pdf",height = 6,width = 8)

FAMs_3d_scale_NEW <- FAMs_3d_scale_NEW[order(FAMs_3d_scale_NEW$avg_log2FC,decreasing = T),]

FAMs_3d_scale_NEW.UP <- FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$avg_log2FC>0),]
FAMs_3d_scale_NEW.DW <- FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$avg_log2FC<0),]


genefeature <- split(FAMs_3d_scale$gene,FAMs_3d_scale$cluster)
genefeature_UP <- split(FAMs_3d_scale_NEW.UP$gene,FAMs_3d_scale_NEW.UP$cluster)
genefeature_DW <- split(FAMs_3d_scale_NEW.UP$gene,FAMs_3d_scale_NEW.UP$cluster)



genefeature_scale.UP_simple <- lapply(genefeature_scale, function(x){return(x[c(1:20)])})
upset(fromList(genefeature_scale.UP_simple),nsets = 10,
      empty.intersections = "off",nintersects = 18,
      set_size.show = T,order.by = "freq")



table(FAMs_sel$cluster)
# Idents(scw)

# FAMs_SCT <- FindAllMarkers(scw,assay = "SCT",slot = "scale.data")
# FAMs_selk <- FAMs_SCT %>% 
#   group_by_("cluster") %>% 
#   filter(pct.1-pct.2>0.2 & avg_log2FC>0.4 & pct.2<0.5) %>% 
#   arrange("diff.pct") %>% 
#   slice_max(order_by = avg_log2FC, n = 25)
# markerVocalno(markers = FAMs_SCT,
#                           topn = 5,
#                           labelCol = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7"))

### Figure 2 Track genes
## …………………………………………………… ##
# 
# Track_genes <- read.csv("Trajectory_genes.fli.csv")
# Track_genes.fli <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

load("E:/works/artery ECs/track.RData")
Track_genes <- Track_genes %>% arrange(desc(morans_I))
Track_genes_sig50 <- Track_genes.fli %>% top_n(n=50, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Track_genes_sig100 <- Track_genes.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Track_genes_sig200 <- Track_genes.fli %>% top_n(n=200, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Track_genes_sig2000 <- Track_genes.fli %>% top_n(n=2000, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# track genes similarity
Idents(scw) <- "X3d.clu"
DefaultAssay(scw) <- "RNA"
scw <- ScaleData(scw)
scw <- NormalizeData(scw)
Idents(scw) <- "X3d.clu"
FAMs_X3d_op <- FindAllMarkers(scw,only.pos = T,logfc.threshold = 0.4)

# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
# 
# library(ggstatsplot)
# library(ggplot2,lib.loc = "D:/Program Files/R/library")
# getwd()
# BiocManager::install("ggstatsplot")
# BiocManager::install("pacman")
# detach("package:ggplot2", unload = TRUE)

# ggcorrmat(data = cor_X3d.clu, colors = c("#B2182B", "white", "#4D4D4D"), title = "Correlalogram Example of ggstatsplot charts makes",
#                 subtitle = "processed charts with ggcorrmat()", caption = "Visualization by DataCharm",
#                 ggtheme = hrbrthemes::theme_ipsum(base_family = "Roboto Condensed"), ) + theme(plot.title = element_text(hjust = 0.5,
#                                                                                                                          vjust = 0.5, color = "black", size = 10, margin = margin(t = 1, b = 12)), plot.subtitle = element_text(hjust = 0,
#                                                                                                                                                                                                                                 vjust = 0.5, size = 8), plot.caption = element_text(face = "bold", size = 10))
# 
cor_heat <- ComplexHeatmap::Heatmap(cor_X3d.clu)
cor_heat
cor_heat <- as.ggplot(cor_heat)
DimPlot(scw,group.by = "NEW.clu",reduction = "scphere",label = T)
ggsave(cor_heat,filename = "Figures/sFig.2A.pdf",height = 6,width = 8)

FeaturePlot(scw,features = "LYVE1",reduction = "scphere",min.cutoff = "q10")

VlnPlot(scw,features = "EDN1")
pattern_group <- data.frame("orig.clu"=scw$X3d.clu)
pattern_group$`merge.clu` <- ifelse(pattern_group$orig.clu == 1,"1",
                                    ifelse(pattern_group$orig.clu == 2,"2",
                                           ifelse(pattern_group$orig.clu %in% c(3,5,6),"3","4")))

scw$`merge_SCT.clu` <- factor(pattern_group$`merge.clu`,levels = 1:4)
Idents(scw) <- "merge_SCT.clu"
FAMs_merge2 <- FindAllMarkers(scw,only.pos = T,logfc.threshold = 0.4,
                              assay = "RNA",slot = "count")
FAMs_merge2$diff.pct <- FAMs_merge2$pct.1-FAMs_merge2$pct.2
FAMs_sel2 <- FAMs_merge2 %>% 
  group_by_("cluster") %>% 
  filter(diff.pct>0.15 & avg_log2FC>0.5 &  pct.1 >0.3) %>% 
  arrange("diff.pct") %>% 
  slice_max(order_by = diff.pct, n = 20)
table(FAMs_sel2$cluster)
Idents(scw) <- "merge_SCT.clu"
DefaultAssay(scw)
FeaturePlot(scw,"RGCC",reduction = "scphere",slot = "scale.data",label = T)+DimPlot(scw,reduction = "scphere",label = T)


# celltype distribution

cellnum1 <-table(scw$GSE,scw$NEW.clu)
group<-rownames(cellnum1)
celltype <- colnames(cellnum1)
cell.prop<-as.data.frame(prop.table(t(cellnum1)))
colnames(cell.prop)<-c("Celltype","sample","proportion")
data4plot1 <-as.data.frame(cell.prop)

p1<-ggplot(cell.prop,aes(sample,proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=hue_pal()(10))+#自定义fill的颜色
  ggtitle("cell proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+RotatedAxis()
p1
ggsave(p1,filename = "Figures/Fig.1D.pdf",height = 2,width = 6)


# original
hemapplot <- AverageHeatmap(object = scw,assays = "RNA",slot = "count",showRowNames = F,annoCol = TRUE,
                            myanCol = c("#ff6633","#9999ff","#99ccff","#66ff66"),
                            markerGene = FAMs_sel2$gene,markGenes = Track_genes_sig50,cluster_rows = F)
hemapplot
hemapplot <- as.ggplot(hemapplot)
ggsave(hemapplot,filename = "Figures/Fig.2B.pdf",height = 6,width = 8)

# Circlize
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
BiocManager::install("dendsort")

expr.count <- AverageExpression(scw,slot = "scale.data",assay = "RNA")[[1]]# 1, 2, 3-7
hm.count <- expr.count[FAMs_sel2$gene,]
FAMs_sel2$gene
#根据数据范围定义热图的颜色梯度：
range(hm.count)
mycol <- colorRamp2(c(-0.8828134, 0, 0.8828134),c("#0da9ce", "white", "#e74a32"))

#环形热图绘制：
circos.heatmap(hm.count,
               col = mycol)
circos.clear()#绘制完成后，需要重置循环布局再开启下一次绘图，不然会遇到报错或轨道继续在同一图中叠加；

circos.par(gap.after = c(50)) #调整圆环首位间距；
circos.heatmap(hm.count,
               col = mycol,
               cluster = TRUE, #是否对行聚类
               dend.side = "inside", #聚类树方向：inside显示在圆环内圈，inside为圆环外圈；
               rownames.side ="outside", #行名方向；
               rownames.col = "black",
               rownames.cex = 0.6)
circos.clear()

#聚类树美化：
circos.par(gap.after=c(90))
circos.heatmap(hm.count,
               col = mycol,
               cluster = FALSE,
               # dend.side = "outside",
               rownames.side = "outside",
               rownames.col = "black",
               rownames.cex = 0.6,
               track.height = 0.5, #轨道高度，即圆环/格子的粗细
               dend.track.height = 0.08, #聚类树高度调整
               dend.callback = function(dend, m, si) { #聚类树的回调
                 color_branches(dend,k = 4,col = 1:4) #修改聚类树颜色
               })

lg <- Legend(title = "Express",
             col_fun = mycol,
             direction = c("vertical"),
             title_position = c('topcenter'))
draw(lg, x = unit(0.9, "npc"), y = unit(0.5, "npc"), just = c("right", "center"))
circos.clear()

# @输出图形 Fig.2B_circ.pdf

## geneset remake
EndMT_genesets <- read.csv("ENDMT.geneset.csv")
EndMT_genesets <- as.list(read.csv("ENDMT.geneset.csv")[,1:4])
EndMT_genesets[["EC.fate"]] <- EndMT_genesets[["EC.fate"]][1:13]
EndMT_genesets[["EndMT"]] <- EndMT_genesets[["EndMT"]][1:34]
EndMT_genesets[["MSC.fate"]] <- EndMT_genesets[["MSC.fate"]][1:5]
EndMT_genesets[["Inflammation"]] <- EndMT_genesets[["Inflammation"]][1:33]


Mouse_genesets <- read.csv("genesets_mouse.csv")
Mouse_genesets <- split(Mouse_genesets$Human.Symbol,Mouse_genesets$Term)



Idents(scw) <- "merge_SCT.clu"
sck <- DietSeurat(scw)
rm(sctmp)


sctmp <- irGSEA.score(object = scw, assay = "RNA",
                      slot = "count", seeds = 1116, ncores = 1,
                      custom = T, geneset = Mouse_genesets,
                      kcdf = 'Gaussian')

# DefaultAssay(sctmp) <- 
result.dge <- irGSEA.integrate(object = sctmp, 
                               group.by = "GSE",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
#> Calculate differential gene set : AUCell
#> Calculate differential gene set : UCell
#> Calculate differential gene set : singscore
#> Calculate differential gene set : ssgsea
class(result.dge)
#> [1] "list"
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "AUCell",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
DefaultAssay(sctmp) <- "AUCell"
DefaultAssay(sctmp) <- "UCell"
DefaultAssay(sctmp) <- "singscore"
DefaultAssay(sctmp) <- "ssgsea"

table(scw$GSE)
sck <- subset(sctmp,GSE=="GSE213740")
Nebulosa::plot_density(sck,"Development",reduction = "scphere")+
  Nebulosa::plot_density(sck,"Differentiation",reduction = "scphere")+
  Nebulosa::plot_density(sck,"ECM-Organization",reduction = "scphere")+
  Nebulosa::plot_density(sck,"Integrin-Signaling",reduction = "scphere")+
  Nebulosa::plot_density(sck,"Lipoprotein-Handling",reduction = "scphere")

FeaturePlot(sctmp,levels(as.factor(result.dge$AUCell$Name)),reduction = "scphere",min.cutoff = "q50",pt.size = 0.4)

names(Mouse_genesets)

DefaultAssay(sctmp) <- "RNA"

Nebulosa::plot_density(sctmp,"ACKR1",reduction = "scphere",slot = "scale.data")
Nebulosa::plot_density(sctmp,"ACKR1",reduction = "scphere",slot = "count")
FeaturePlot(scw,features = "ACKR1",reduction = "scphere")

sck <- FindNeighbors(scw,reduction = "scphere",graph.name = "scp_graph",dims = 1:2)
sck <- FindClusters(sck,graph.name = "scp_graph",resolution = 0.3)
sck$NEW.clu <- scw$NEW.clu
DimPlot(sck,group.by = "NEW.clu",reduction = "scphere")

irGSEA.density.scatterplot(object = sctmp,
                           method = "UCell",
                           show.geneset = "Cluster 4",
                           reduction = "scphere")


densityheatmap <- irGSEA.densityheatmap(object = sctmp,
                                        method = "UCell",
                                        show.geneset = "ECM-Organization",cluster.levels = c(1,3,4,2))










# Gene Module
Track_genes$merge.clu <- FAMs_merge2[match(Track_genes$gene_short_name,FAMs_merge2$gene),6]

## module selection
Track_match <- data.frame(
  "Track_genes_sig50" = Track_genes_sig50
)

Track_match[["cluster"]] = FAMs_X3d_op[match(Track_genes_sig50,FAMs_X3d_op$gene),"cluster"]
Track_match[["merge.cluster"]] = FAMs_merge2[match(Track_match$Track_genes_sig50,FAMs_merge2$gene),"cluster"]
# Track_match[["module"]] = as.data.frame(gene_module[match(Track_genes_sig100,gene_module$id),"module"])[,1]
Track_match[["module2"]] = as.data.frame(gene_module[match(Track_match$Track_genes_sig50,gene_module$id),"module"])[,1]

tmps <- cbind(table(Track_match$merge.cluster,Track_match$module2))
tmps <- tmps[,order(colSums(tmps),decreasing = T)]
tmps <- tmps[,which(colSums(tmps)!=0)]
colnames(tmps)
# [1] "4"  "21" "91" "33" "77" "23" "30" "36" "54" "62" "64" "68" "1"  "8"  "51" "61" "69" "81"
plot_cells(cds,
           genes=gene_module %>% filter(module %in% c(4)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)
plot_cells(cds,
           genes=gene_module %>% filter(module %in% c(21)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)

genelist <- pull(Track_genes.fli, gene_short_name) %>% as.character()

gene_module <- find_gene_modules(cds[genelist,], resolution=5e-3, cores = 1)
gene_module.list <- gene_module
gene_module.list$name <- paste("Module",gene_module.list$module,sep = " ")
# table(gene_module$supermodule)

merge.clu <- as.data.frame(scw$merge_SCT.clu)
merge.clu <- merge.clu[row.names(colData(cds)),]
cds@colData@listData[["merge.clu"]] <- merge.clu

cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=cds@colData@listData[["merge.clu"]])


agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
agg_sel <- agg_mat[paste0("Module ",colnames(tmps)),]

pheatmap::pheatmap(agg_sel, scale="column", clustering_method="ward.D2",
                   cellwidth = 10, cellheight = 10,border_color = "white",
                   cluster_cols = F)











GSEs <- levels(factor(scw$GSE))
scFEA <- t(matrix(rep(0,168)))
for (i in 1:length(GSEs)) {
  flux_tmp <- read.csv(paste0("./scFEA/",GSEs[i],"_flux.csv"), header = T, row.names = 1)
  flux_tmp <- data.matrix(flux_tmp)
  scFEA <- rbind(scFEA,flux_tmp)
}
colnames(scFEA) <- colnames(flux_tmp)
scFEA <- scFEA[scw@assays[["RNA"]]@counts@Dimnames[[2]],]

scw@assays$scFEA <- CreateAssayObject(counts = t(scFEA))
scw@assays[["scFEA"]]@scale.data <- scale(t(scFEA))

DefaultAssay(scw) <- "scFEA"

DefaultAssay(scw) <- "RNA"
Idents(scw) <- "celltype"
levels(Idents(scw))[1]

fams_scfea <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = "ITLN1+ EC",
                                   assay = "scFEA")
fams_scfea <- FindAllMarkers(scw,logfc.threshold = 0,min.pct = 0)
fams_scfea_fli <- fams_scfea %>% dplyr::filter(pct.1-pct.2>0.05)
fams_scfea_fli <- split(fams_scfea_fli$gene,fams_scfea_fli$cluster)
UpSetR::upset(UpSetR::fromList(fams_scfea_fli))

scFEA_dict <- read.csv("scFEA/fea.module.csv",row.names = 1)
scFEA_dict$X <- rownames(scFEA_dict)
setdiff(1:171,scFEA_dict$Module_id)
scFEA_dict$X <- stringr::str_replace(scFEA_dict$X,pattern = "_",replacement = "-")
fams_scfea_fli <- lapply(fams_scfea_fli, function(x){return(scFEA_dict[match(x,scFEA_dict$X),6])})

scFEA_dict.sel <- scFEA_dict[which(scFEA_dict$X %in% unlist(fams_scfea_fli)),]
# ComplexHeatmap::Heatmap(scFEA[,unlist(fams_scfea_fli)])
fams_scfea_df <- as.data.frame(
  unlist(fams_scfea_fli)
  )
names(fams_scfea_df) <- "Moudles"
fams_scfea_df$celltype <- c(rep(names(fams_scfea_fli)[1],length(fams_scfea_fli[[1]])),
                            rep(names(fams_scfea_fli)[2],length(fams_scfea_fli[[2]])),
                            rep(names(fams_scfea_fli)[3],length(fams_scfea_fli[[3]])),
                            rep(names(fams_scfea_fli)[4],length(fams_scfea_fli[[4]])))
# fams_scfea_df <- fams_scfea_df[!duplicated(fams_scfea_df$Moudles),]
fams_scfea_df$TF <- !duplicated(fams_scfea_df$Moudles)
rownames(fams_scfea_df) <- 1:37

fams_scfea_df[!fams_scfea_df$TF,1]
fams_scfea_df[fams_scfea_df$Moudles %in% fams_scfea_df[!fams_scfea_df$TF,1],3] <- FALSE
fams_scfea_df[9,2] <- "ACKR1+ EC"
fams_scfea_df[23,2] <- "ITLN1+ EC"
fams_scfea_df[27,2] <- "ACKR1+ EC"
fams_scfea_df[29,2] <- "ACKR1+ EC"
fams_scfea_df <- fams_scfea_df[-c(1,37,7,10,12),] 

avg_fea <- AverageExpression(scw,assays = "scFEA",slot = "data",group.by = "celltype",
                             features = unique(unlist(fams_scfea_fli)))[[1]]

avg_fea <- t(scale(t(avg_fea)))
rownames(avg_fea)

hmcol <- scales::hue_pal()(26)[5:26]
names(hmcol) <- 1:22
row_su <- rowAnnotation(
  Super.Module = factor(scFEA_dict[match(rownames(avg_fea),scFEA_dict$X),6],levels = 1:22),
  col = list(
    Super.Module = hmcol)
)



fams_scfea_list_new <- split(fams_scfea_df$Moudles,fams_scfea_df$celltype)
UpSetR::upset(UpSetR::fromList(fams_scfea_fli))

avg_fea["M-169",]
avg_fea["M-61",]
fams_scfea_df$celltype <- factor(fams_scfea_df$celltype,
                                 levels = c("ITLN1+ EC", "TCIM+ EC","ACKR1+ EC", "Lymphatic EC"))

fams_scfea_df <- fams_scfea_df[order(fams_scfea_df$celltype),]
hmcol2 <- scales::hue_pal()(4) 
names(hmcol2) <- levels(factor(scw$celltype))
row_type <- rowAnnotation(
  celltype = factor(fams_scfea_df[match(rownames(avg_fea),fams_scfea_df$Moudles),2]),
  col = list(
    celltype = hmcol2)
)
avg_fea <- avg_fea[fams_scfea_df$Moudles,]
# pheatmap::pheatmap(avg_fea)
scFEA_hmplot <- ComplexHeatmap::Heatmap(avg_fea,cluster_columns = T,
                        cluster_rows = T,right_annotation = row_su,left_annotation = row_type)

# 1.创建
pdf( 
  file = "./Figures/Fig.5A.pdf", # 文件名称
  width = 6,           # 宽
  height = 6)              # 分辨率
# 2. 绘图
scFEA_hmplot
# 3. 关闭画布
dev.off()


scFEA_genes <- read.csv("scFEA/genes.csv",header = F,row.names = 1)

scFEA_count <- scFEA_genes[,1]
scFEA_genes <- scFEA_genes[,-1]
scFEA_genes <- as.data.frame(t(scFEA_genes))
scFEA_genes <- as.list(scFEA_genes)
for (i in 1:length(scFEA_genes)) {
  scFEA_genes[[i]] <- scFEA_genes[[i]][1:scFEA_count[i]]
}

scFEA_dict$old <- stringr::str_replace(rownames(scFEA_dict),"-","_")
setdiff(scFEA_genes$V1,scFEA_dict$old)
setdiff(scFEA_dict$old,scFEA_genes$V1)

fea_genes <- as.data.frame(unique(unlist(scFEA_genes)))
names(fea_genes) <- "Gene.Symbol"








DefaultAssay(scw) <- "RNA"
RidgePlot(scw,features = unique(unlist(scFEA_genes["M_58"])),stack = T)

load("E:/works/artery ECs/Track_gene_vasEC.RData")
Track_genes_NEW.fli <- cdst_track[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3) %>% arrange(.,desc(morans_I))
Track_genes_NEW.fli$clu <- FAMs_3d_scale_NEW[match(Track_genes_NEW.fli$gene_short_name,FAMs_3d_scale_NEW$gene),6]
Track_genes_NEW.fli[which(Track_genes_NEW.fli$clu %in% c(1,2,4,6,9,10)),8] <- "ACKR1+ EC"
Track_genes_NEW.fli[which(Track_genes_NEW.fli$clu %in% c(5,7)),8] <- "TCIM+ EC"
Track_genes_NEW.fli[which(Track_genes_NEW.fli$clu %in% c(3)),8] <- "ITLN1+ EC"
Track_genes_NEW.fli[which(Track_genes_NEW.fli$clu %in% c(8)),8] <- "Lymphatic EC"

Track_genes_NEW.fli_sig50 <- Track_genes_NEW.fli[1:50,]


Track_genes.fli$clu <- FAMs_3d_scale_NEW[match(Track_genes.fli$gene_short_name,FAMs_3d_scale_NEW$gene),6]
Track_genes.fli[which(Track_genes.fli$clu %in% c(1,2,4,6,9,10)),8] <- "ACKR1+ EC"
Track_genes.fli[which(Track_genes.fli$clu %in% c(5,7)),8] <- "TCIM+ EC"
Track_genes.fli[which(Track_genes.fli$clu %in% c(3)),8] <- "ITLN1+ EC"
Track_genes.fli[which(Track_genes.fli$clu %in% c(8)),8] <- "Lymphatic EC"







Track_genes.fli_sig50 <- Track_genes.fli[1:50,]




DefaultAssay(scw) <- "scFEA"
Nebulosa::plot_density(scw,features = "M-2",reduction = "scphere")
RidgePlot(scw,features = "M-58",assay = "FLUX",group.by = "NEW.clu")
scw@assays$RNA@key <- "rna_"
scw@assays$scFEA@key <- "scfea_"
scw@assays$SCT@key





