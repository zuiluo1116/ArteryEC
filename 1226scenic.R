library(SCENIC) 
library(Seurat) 
# devtools::install_github("aertslab/SCopeLoomR")
library(SCopeLoomR)

scenicLoomPath='E:/works/artery ECs/scenic/sample_SCENIC.loom' 
scenicLoomPath2='E:/works/artery ECs/scenic/sample.loom' 
# E:\works\artery ECs\scenic
# E:/works/GX/cal/1124/sample_SCENIC.loom
loom<-open_loom(scenicLoomPath)
loom2<-open_loom(scenicLoomPath2)
# Read i吁ormationfrom Loom file: 


regulons_incidMat <-get_regulons(loom, column.attr.name="Regulons") 
regulons <-regulonsToGeneLists(regulons_incidMat) 
regulonAUC <-get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
 
library(pheatmap) 
n=t(scale(t(getAUC(regulonAUC[,])))) #'scale'可以对Log-ratio数值进行归一化 n
n[n>2]=2 
n[n<-2]= -2 
n[1:4,1:4] 
dim(n) 
Idents(sct) <- "NEW.clu"
ac=data.frame(Cell.Type= as.character(Idents(sct))) 
n[1:4,1:4] 
n=n[,colnames(n) %in% colnames(sct)] 
rownames(ac)=colnames(n) 

ar <- data.frame(TFs = rownames(n))
rownames(ar) <- ar$TFs

ar$Anno <- stringr::str_sub(ar$TFs,1,nchar(ar$TFs)-3)

# cg=read.table('choose_tf.txt')[,1] 
# cg 
# cg_n=n[rownames(n) %in% cg,] 
# pheatmap(n,show_colnames =F,show_rownames = T, annotation_col=ac) 

table(ac$group) 
#尊重作者， 进行二分类！
p1 <- ComplexHeatmap::pheatmap(n,show_colnames =F,show_rownames = F, 
         annotation_col=ac,use_raster=F,annotation_colors = ann_colors)#, 
         # filename = 'heatmap_choose_regulon.pdf') 
dev.off() 
p2 <- p1 + rowAnnotation(link = anno_mark(at = which(rownames(n) %in% labs), 
                                   labels = labs, labels_gp = gpar(fontsize = 10)))
pdf( 
  file = "./Figures/sFig.6A.pdf", # 文件名称
  width = 6,           # 宽
  height = 6)
p2
dev.off() 

# 展示特定行名函数

labs <- ar[match(tmp$TFs,ar$Anno),1]
labs <- as.character(na.omit(labs))
add.flag(p1,kept.labels = labs,repel.degree = 0.2)


cellanno <- sct@meta.data$NEW.clu
RSS <- SCENIC::calcRSS(regulonAUC,cellAnnotation = cellanno)





regulons[unique(labs)]



RSS.sel <- RSS[unique(labs),]


col_anno2 <- melt(table(scw$celltype,scw$NEW.clu))
col_anno2 <- col_anno2[which(col_anno2$value!=0),]
rownames(col_anno2) <- col_anno2$Var2
col_anno2 <- as.data.frame(col_anno2[,1])
names(col_anno2) <- "Cell.Type"
ann_colors2 = list(
  Cell.Type= c(`ITLN1+ EC` = "#F8766D",
               `TCIM+ EC` = "#7CAE00",
               `ACKR1+ EC` = "#00BFC4",
               `Lymphatic EC` = "#C77CFF")
)

RSS.sel <- RSS.sel[,as.character(manual_order)]
pmtf <- ComplexHeatmap::pheatmap(RSS.sel, scale="column", clustering_method="ward.D2",
                         cellwidth = 20, cellheight = 10,border_color = "white",
                         annotation_col = col_anno2,
                         annotation_colors = ann_colors2,
                         cluster_cols = F,cluster_rows = F,show_rownames = T,use_raster =F)#+ 
  # rowAnnotation(link = anno_mark(at = which(rownames(n) %in% labs),
  #                                labels = labs, labels_gp = gpar(fontsize = 10)))

manual_order = c(3,5,7,6,9,4,1,10,2,8)
pmtf@column_order <- manual_order
pmea2 <- ggplotify::as.ggplot(pmea2)
save(RSS,file = "RSS.rdata")





RSS.df <- as.data.frame(RSS)

for (i in 1:nrow(RSS.df)) {
  # gene_module.df[i,7] <- max(gene_module.df[i,3:6])
  RSS.df[i,5] <- names(which.max(RSS.df[i,1:4]))
}

for (i in 1:nrow(RSS.df)) {
  l <- as.numeric(which.max(RSS.df[i,1:4]))
  k <- setdiff(1:4,l)
  RSS.df[i,6] <- (RSS.df[i,l]-RSS.df[i,k[1]])+
    (RSS.df[i,l]-RSS.df[i,k[2]])+
    (RSS.df[i,l]-RSS.df[i,k[3]])
}

RSS.df$V5 <- factor(RSS.df$V5,levels = levels(factor(scw$celltype)))

RSS.df$TFs <- rownames(RSS.df)
















#######################################
# ori
adj='E:/works/artery ECs/scenic/adj.sample.tsv' 

adj <- read.table(adj,header = T)

adj_fli <- adj[which(adj$importance>1),]
adj_split <- split(adj_fli$target,adj_fli$TF)
adj_split2 <- split(adj_fli$TF,adj_fli$target)
adj_split.sel <- adj_split[stringr::str_sub(names(regulons),1,nchar(names(regulons))-3)]


ITLN1_TF_Scenic <- adj_split2[FAMs_type$gene[which(FAMs_type$cluster=="ITLN1+ EC")]]
ITLN1_TF_Scenic.freq=as.data.frame(table(unlist(ITLN1_TF_Scenic)))
ITLN1_TF_Scenic.freq <- ITLN1_TF_Scenic.freq[order(ITLN1_TF_Scenic.freq$Freq,decreasing = T),]
ITLN1_TF_Scenic.freq$index <- ITLN1_TF_Scenic.freq$Freq/length(which(FAMs_type$cluster=="ITLN1+ EC"))

TCIM_TF_Scenic <- adj_split2[FAMs_type$gene[which(FAMs_type$cluster=="TCIM+ EC")]]
TCIM_TF_Scenic.freq=as.data.frame(table(unlist(TCIM_TF_Scenic)))
TCIM_TF_Scenic.freq <- TCIM_TF_Scenic.freq[order(TCIM_TF_Scenic.freq$Freq,decreasing = T),]
TCIM_TF_Scenic.freq$index <- TCIM_TF_Scenic.freq$Freq/length(which(FAMs_type$cluster=="TCIM+ EC"))


ACKR1_TF_Scenic <- adj_split2[FAMs_type$gene[which(FAMs_type$cluster=="ACKR1+ EC")]]
ACKR1_TF_Scenic.freq=as.data.frame(table(unlist(ACKR1_TF_Scenic)))
ACKR1_TF_Scenic.freq <- ACKR1_TF_Scenic.freq[order(ACKR1_TF_Scenic.freq$Freq,decreasing = T),]
ACKR1_TF_Scenic.freq$index <- ACKR1_TF_Scenic.freq$Freq/length(which(FAMs_type$cluster=="ACKR1+ EC"))


Lymphatic_TF_Scenic <- adj_split2[FAMs_type$gene[which(FAMs_type$cluster=="Lymphatic EC")]]
Lymphatic_TF_Scenic.freq=as.data.frame(table(unlist(Lymphatic_TF_Scenic)))
Lymphatic_TF_Scenic.freq <- Lymphatic_TF_Scenic.freq[order(Lymphatic_TF_Scenic.freq$Freq,decreasing = T),]
Lymphatic_TF_Scenic.freq$index <- Lymphatic_TF_Scenic.freq$Freq/length(which(FAMs_type$cluster=="Lymphatic EC"))



VlnPlot(sct,features = "ITLN1",group.by = "orig.ident")
FeaturePlot(sct,features = "ITLN1",split.by = "orig.ident")

DimPlot(sct,split.by = "orig.ident")
FeaturePlot(sct,features = "TCIM",split.by = "orig.ident")

VlnPlot(sct,features = "TCIM",split.by = "orig.ident",group.by = "orig.ident",assay = "RNA",slot = "count")




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
ggsave(p1,filename = "Figures/Fig.2B2.pdf",height = 8,width = 4)

save(scw,file = "scw20221227.RData")
rm(scw)
