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




Meta_out <- 
colnames(balance)
BAL_TCA <- BAL_tmp[Metabolites,]
rownames(BAL_tmp)
ComplexHeatmap::pheatmap(BAL_tmp, scale="none", clustering_method="ward.D2",
                         cellwidth = 20, cellheight = 10,border_color = "white",
                         #annotation_col = col_anno,annotation_colors = ann_colors,
                         cluster_cols = F,cluster_rows = F,show_rownames = T,use_raster =F)








scw$TCA2Gly <- FLU_all



VlnPlot(scw,"TCA2Gly")
RidgePlot(scw,"TCA2Gly")
FeaturePlot(scw,"TCA2Gly",reduction = "scphere")
Nebulosa::plot_density(scw,"TCA2Gly",reduction = "scphere")


























DimPlot(scw,group.by = "NEW.clu",reduction = "scphere",label = T)





flux <- t(matrix(rep(0,168)))
for (i in 1:length(GSEs)) {
  flux_tmp <- read.csv(paste0("E:/works/artery ECs//scFEA/",GSEs[i],"_flux.csv"), header = T, row.names = 1)
  flux_tmp <- data.matrix(flux_tmp)
  flux <- rbind(flux,flux_tmp)
}
colnames(flux) <- colnames(flux_tmp)
flux <- flux[scw@assays[["RNA"]]@counts@Dimnames[[2]],]

scw@assays$FLUX <- CreateAssayObject(counts = t(flux))
scw@assays[["FLUX"]]@scale.data <- scale(t(flux))





# pheatmap(exprTable, cluster_cols = col_cluster)
pmea <- ComplexHeatmap::pheatmap(BAL_TCA, scale="none", clustering_method="ward.D2",
                         cellwidth = 40, cellheight = 15,border_color = NA,
                         #annotation_row = row_anno,#annotation_col = ann_colors,
                         cluster_cols = F,cluster_rows = F,show_rownames = T,use_raster =F)










pmea@column_order <- manual_order
pmea@column_dend_param[["obj"]][["order"]] <- manual_order
pmea@column_dend_param[["reorder"]] <- T

Flow <- c(paste0("M-",1:14),paste0("M-",71:78),"M-34","M-35","M-105")
FLU_TCA <- FLU_tmp[Flow,]

setdiff(c("Glucose.in"),rownames(BAL_tmp))
ComplexHeatmap::pheatmap(FLU_TCA, scale="none", clustering_method="ward.D2",
                         cellwidth = 40, cellheight = 15,border_color = NA,
                         #annotation_row = row_anno,#annotation_col = ann_colors,
                         cluster_cols = F,cluster_rows = F,show_rownames = T,use_raster =F)











FLU_tmp <- AverageExpression(sct,assays = "FLUX",slot = "count")
FLU_tmp <- FLU_tmp[[1]]

FLU_Glycogen <- FLU_tmp["M-169",]


FLU_TCA <- FLU_tmp[1:14,]

celltype_df <- as.data.frame(sct$celltype)
celltype_df$barcode <- rownames(celltype_df)
celltype_df <- celltype_df[order(celltype_df[,1]),]
flux_TCA <- flux[rownames(celltype_df),1:14]


row_anno <- data.frame(celltype_df[,1])
rownames(row_anno) <- celltype_df$barcode
names(row_anno) <- "Cell.Type"
ann_colors = list(
  Cell.Type = c(`ITLN1+ EC` = "#F8766D",
                `TCIM+ EC` = "#7CAE00",
                `ACKR1+ EC` = "#00BFC4",
                `Lymphatic EC` = "#C77CFF")
)

ComplexHeatmap::pheatmap(flux_TCA, scale="none", clustering_method="ward.D2",
                   cellwidth = 20, cellheight = 0.05,border_color = NA,
                   annotation_row = row_anno,#annotation_col = ann_colors,
                   cluster_cols = F,cluster_rows = F,show_rownames = F,use_raster =F)


FeatureScatter(scw,feature1 = "M-5",feature2 = "M-6")
FeatureScatter(scw,feature1 = "M-35",feature2 = "M-34")
RidgePlot(scw,"M-6")











































BAL_sc <- t(scale(t(BAL_tmp)))


BALANCE_top5_list <- list(
  BALANCE_ITLN1 = BALANCE_ITLN1_rank$Name[1:5],
  BALANCE_TCIM = BALANCE_TCIM_rank$Name[1:5],
  BALANCE_ACKR1 = BALANCE_ACKR1_rank$Name[1:5],
  BALANCE_Lymphatic = BALANCE_Lymphatic_rank$Name[1:5]
)
un_FEAs <- unique(as.character(unlist(BALANCE_top5_list)))

BAL_sel <- BAL_sc[un_FEAs,]
pheatmap::pheatmap(BAL_sc, scale="none", clustering_method="ward.D2",
                   cellwidth = 10, cellheight = 10,border_color = "white",
                   cluster_cols = F)

VlnPlot(scw,features = "Glycogen")




DefaultAssay(scw) <- "BALANCE"
DimPlot(scw,reduction = "scphere",label=T,group.by = "NEW.clu")
Idents(scw) <- "celltype"
levels(Idents(scw))

BALANCE_ITLN1 <- FindConservedMarkers(scw,ident.1 = "ITLN1+ EC",grouping.var = "GSE",
                                      assay = "BALANCE",logfc.threshold = 0,min.pct = 0)
levels(factor(scw$GSE))

BALANCE_ITLN1_list <- list(
  GSE155468 = rownames(BALANCE_ITLN1[order(BALANCE_ITLN1$GSE155468_avg_log2FC,decreasing = T),]),
  GSE155514 = rownames(BALANCE_ITLN1[order(BALANCE_ITLN1$GSE155514_avg_log2FC,decreasing = T),]),
  GSE159677 = rownames(BALANCE_ITLN1[order(BALANCE_ITLN1$GSE159677_avg_log2FC,decreasing = T),]),
  GSE213740 = rownames(BALANCE_ITLN1[order(BALANCE_ITLN1$GSE213740_avg_log2FC,decreasing = T),]),
  GSE216860 = rownames(BALANCE_ITLN1[order(BALANCE_ITLN1$GSE216860_avg_log2FC,decreasing = T),])
)
BALANCE_ITLN1_rank=aggregateRanks(BALANCE_ITLN1_list)


BALANCE_TCIM <- FindConservedMarkers(scw,ident.1 = "TCIM+ EC",grouping.var = "GSE",
                                      assay = "BALANCE",logfc.threshold = 0,min.pct = 0)
BALANCE_TCIM_list <- list(
  GSE155468 = rownames(BALANCE_TCIM[order(BALANCE_TCIM$GSE155468_avg_log2FC,decreasing = T),]),
  GSE155514 = rownames(BALANCE_TCIM[order(BALANCE_TCIM$GSE155514_avg_log2FC,decreasing = T),]),
  GSE159677 = rownames(BALANCE_TCIM[order(BALANCE_TCIM$GSE159677_avg_log2FC,decreasing = T),]),
  GSE213740 = rownames(BALANCE_TCIM[order(BALANCE_TCIM$GSE213740_avg_log2FC,decreasing = T),]),
  GSE216860 = rownames(BALANCE_TCIM[order(BALANCE_TCIM$GSE216860_avg_log2FC,decreasing = T),])
)
BALANCE_TCIM_rank=aggregateRanks(BALANCE_TCIM_list)


BALANCE_ACKR1 <- FindConservedMarkers(scw,ident.1 = "ACKR1+ EC",grouping.var = "GSE",
                                      assay = "BALANCE",logfc.threshold = 0,min.pct = 0)
levels(factor(scw$GSE))

BALANCE_ACKR1_list <- list(
  GSE155468 = rownames(BALANCE_ACKR1[order(BALANCE_ACKR1$GSE155468_avg_log2FC,decreasing = T),]),
  GSE155514 = rownames(BALANCE_ACKR1[order(BALANCE_ACKR1$GSE155514_avg_log2FC,decreasing = T),]),
  GSE159677 = rownames(BALANCE_ACKR1[order(BALANCE_ACKR1$GSE159677_avg_log2FC,decreasing = T),]),
  GSE213740 = rownames(BALANCE_ACKR1[order(BALANCE_ACKR1$GSE213740_avg_log2FC,decreasing = T),]),
  GSE216860 = rownames(BALANCE_ACKR1[order(BALANCE_ACKR1$GSE216860_avg_log2FC,decreasing = T),])
)
BALANCE_ACKR1_rank=aggregateRanks(BALANCE_ACKR1_list)

BALANCE_Lymphatic <- FindConservedMarkers(scw,ident.1 = "Lymphatic EC",grouping.var = "GSE",
                                      assay = "BALANCE",logfc.threshold = 0,min.pct = 0)
levels(factor(scw$GSE))

BALANCE_Lymphatic_list <- list(
  GSE155468 = rownames(BALANCE_Lymphatic[order(BALANCE_Lymphatic$GSE155468_avg_log2FC,decreasing = T),]),
  GSE155514 = rownames(BALANCE_Lymphatic[order(BALANCE_Lymphatic$GSE155514_avg_log2FC,decreasing = T),]),
  GSE159677 = rownames(BALANCE_Lymphatic[order(BALANCE_Lymphatic$GSE159677_avg_log2FC,decreasing = T),]),
  GSE213740 = rownames(BALANCE_Lymphatic[order(BALANCE_Lymphatic$GSE213740_avg_log2FC,decreasing = T),]),
  GSE216860 = rownames(BALANCE_Lymphatic[order(BALANCE_Lymphatic$GSE216860_avg_log2FC,decreasing = T),])
)
BALANCE_Lymphatic_rank=aggregateRanks(BALANCE_Lymphatic_list)






ag$Freq=freq[match(ag$Name,freq$Var1),2]

BALANCE_TCIM <- FindConservedMarkers(scw,ident.1 = "TCIM+ EC",grouping.var = "GSE")
BALANCE_ACKR1 <- FindConservedMarkers(scw,ident.1 = "ACKR1+ EC",grouping.var = "GSE")
BALANCE_Lymphatic <- FindConservedMarkers(scw,ident.1 = "Lymphatic EC",grouping.var = "GSE")








fams_balance <- FindAllMarkers(scw,logfc.threshold = 0,min.pct = 0,only.pos = T)
fams_balance_fli <- fams_balance %>% dplyr::filter(avg_log2FC>log2(1) & p_val <0.05 & pct.1-pct.2>0.05)

DefaultAssay(scw) <- "FLUX"
fams_flux <- FindAllMarkers(scw,logfc.threshold = 0,min.pct = 0,only.pos = T)
fams_flux_fli <- fams_flux %>% dplyr::filter(avg_log2FC>log2(1) & p_val <0.05 & pct.1-pct.2>0.05)
log2(1.2)
