


FEA_Genes_tfs <- split(FEA_Genes_dict$Gene_Symbol,FEA_Genes_dict$Supermodule_id)
FEA_Genes_tfm <- split(FEA_Genes_dict$Gene_Symbol,FEA_Genes_dict$Module_ID)

GLY_TFs <- ChEA.TFs.RRA(FEA_Genes_tfm$M_6)

TCA_TFs <- ChEA.TFs.RRA(FEA_Genes_tfs$`1`)
FFA_TFs <- ChEA.TFs.RRA(FEA_Genes_tfs$`4`)
CHO_TFs <- ChEA.TFs.RRA(FEA_Genes_tfs$`22`)


meta_TFs <- list(
  "TCA_TFs" = TCA_TFs$Name[1:10],
  "FFA_TFs" <- FFA_TFs$Name[1:10],
  "CHO_TFs" <- CHO_TFs$Name[1:10],
  "ITLN1_EC_TFs" <- as.character(ITLN1_EC_TFs$TFs),
  "TCIM_EC_TFs"<- as.character(TCIM_EC_TFs$TFs),
  "ACKR1_EC_TFs" <- as.character(ACKR1_EC_TFs$TFs),
  "Lymphatic_EC_TFs" <- as.character(Lymphatic_EC_TFs$TFs)
)
names(meta_TFs) <- c("TCA_TFs","FFA_TFs","CHO_TFs",
                     "ITLN1_EC_TFs","TCIM_EC_TFs","ACKR1_EC_TFs","Lymphatic_EC_TFs")
UpSetR::upset(fromList(meta_TFs),keep.order = T,nsets = 7)

intersect(tmp$TFs,TCA_TFs$Name[1:50])
intersect(tmp$TFs,FFA_TFs$Name[1:50])
intersect(tmp$TFs,CHO_TFs$Name[1:50])
intersect(tmp$TFs,GLY_TFs$Name[1:50])

rm(tmp_probe)

Idents(scw)
FAMs_celltype <- (scw,only.pos = T)
library(tidyverse)
FAMs_celltype_sel <- FAMs_celltype %>% filter(pct.1 > 0.8)


gene_dict <- bitr(scw@assays[["RNA"]]@data@Dimnames[[1]],fromType = "SYMBOL",
                  toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)

FAMs_I <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = "ITLN1+ EC",only.pos = T,
                               logfc.threshold = 0.2,slot = "data")
FAMs_I$ENTREZID <- gene_dict[match(rownames(FAMs_I),gene_dict$SYMBOL),2]

FAMs_T <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = "TCIM+ EC",only.pos = T,
                               logfc.threshold = 0.2,slot = "data")
FAMs_T$ENTREZID <- gene_dict[match(rownames(FAMs_T),gene_dict$SYMBOL),2]

FAMs_A <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = "ACKR1+ EC",only.pos = T,
                               logfc.threshold = 0.2,slot = "data")
FAMs_A$ENTREZID <- gene_dict[match(rownames(FAMs_A),gene_dict$SYMBOL),2]

FAMs_L <- FindConservedMarkers(scw,grouping.var = "GSE",ident.1 = "Lymphatic EC",only.pos = T,
                               logfc.threshold = 0.2,slot = "data")
FAMs_L$ENTREZID <- gene_dict[match(rownames(FAMs_L),gene_dict$SYMBOL),2]

total_gene <- list(
  "FAMs_I" = rownames(FAMs_I),
  "FAMs_T" = rownames(FAMs_T),
  "FAMs_A" = rownames(FAMs_A),
  "FAMs_L" = rownames(FAMs_L)
)
comp_list <- list(
  "FAMs_I" = FAMs_I$ENTREZID,
  "FAMs_T" = FAMs_T$ENTREZID,
  "FAMs_A" = FAMs_A$ENTREZID,
  "FAMs_L" = FAMs_L$ENTREZID
)

library(org.Hs.eg.db)
library(enrichplot)
library(GOplot)
library(ggpubr)
library(ggrepel)
data(EC)


list_gene_I <- FindMarkers(scw,features = rownames(FAMs_I),ident.1 = "ITLN1+ EC",
                            min.pct = 0,min.cells.feature = 0,logfc.threshold = 0)
list_gene_T <- FindMarkers(scw,features = rownames(FAMs_T),ident.1 = "TCIM+ EC",
                           min.pct = 0,min.cells.feature = 0,logfc.threshold = 0)
list_gene_A <- FindMarkers(scw,features = rownames(FAMs_A),ident.1 = "ACKR1+ EC",
                           min.pct = 0,min.cells.feature = 0,logfc.threshold = 0)
list_gene_L <- FindMarkers(scw,features = rownames(FAMs_L),ident.1 = "Lymphatic EC",
                           min.pct = 0,min.cells.feature = 0,logfc.threshold = 0)

list_gene <- rbind(list_gene_I,list_gene_T,list_gene_A,list_gene_L)
list_gene$gene <- rownames(list_gene)
list_gene$diff <- list_gene$pct.1-list_gene$pct.2
list_gene <- list_gene[,c(6,2,3,4,1,5,7)]
names(list_gene) <- names(EC$genelist)

com_go.new <- compareCluster(geneClusters = comp_list,
                             fun = "enrichGO",OrgDb= org.Hs.eg.db,pvalueCutoff=0.5)
com_go.new <- simplify(com_go.new)
com_go.new <- pairwise_termsim(com_go.new)
com_go.new <- setReadable(com_go.new,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
com.GO_plot <- dotplot(com_go.new,label_format=40) +   theme_bw()

com_GO <- com_go.new@compareClusterResult
com_GO <- com_GO[,c(1,2,3,9,7)]
names(com_GO) <- names(EC$david)
com_GO$Genes <- stringr::str_replace_all(com_GO$Genes,pattern = "/",replacement = ",")

names(EC$genelist)

go_circ <- circle_dat(com_GO, list_gene)
scaleFUN <- function(x) sprintf("%.2f", x) 
go_circ$category <- factor(go_circ$category,levels = c("FAMs_A","FAMs_T","FAMs_I","FAMs_L"))
group <- levels(go_circ$category)

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
lmts.y <- c(0,range(sub$adj_pval)[2]+0.5)
terms <- sub[duplicated(sub$term),3]
sub[which(sub$term %in% terms),3]

g <- ggplot(sub, aes(zscore, adj_pval, fill = category, 
                     size = count)) + 
  labs(title = title, x = "z-score",
       y = "-log (adj p-value)") + 
  geom_point(shape = 21, col = "black",
             alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.5), col = "orange") + 
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


install.packages('presto')
library(presto)


scw <- scaleda
gsea_tmp <- AverageExpression(scw,group.by = "celltype",slot = "scale.data")
gsea_tmp <- gsea_tmp[["RNA"]]
gse_dict <- gene_dict[match(rownames(gsea_tmp),gene_dict$SYMBOL),1:2]
gse_dict <- na.omit(gse_dict)
gsea_sel <- gsea_tmp[gse_dict$SYMBOL,] 
rownames(gsea_sel) <- gse_dict$ENTREZID

gsea_A <- as.numeric(gsea_sel[order(gsea_sel[,3],decreasing = T),3])
names(gsea_A) <- rownames(gsea_sel)[order(gsea_sel[,3],decreasing = T)]

gsea_T <- as.numeric(gsea_sel[order(gsea_sel[,2],decreasing = T),2])
names(gsea_T) <- rownames(gsea_sel)[order(gsea_sel[,2],decreasing = T)]

gsea_I <- as.numeric(gsea_sel[order(gsea_sel[,1],decreasing = T),1])
names(gsea_I) <- rownames(gsea_sel)[order(gsea_sel[,1],decreasing = T)]

gsea_L <- as.numeric(gsea_sel[order(gsea_sel[,4],decreasing = T),4])
names(gsea_L) <- rownames(gsea_sel)[order(gsea_sel[,4],decreasing = T)]




gseGO_A <- gseGO(gsea_A,ont = "BP",OrgDb = org.Hs.eg.db,seed = 1116,keyType = "ENTREZID",pvalueCutoff = 0.05)
gseGO_T <- gseGO(gsea_T,ont = "BP",OrgDb = org.Hs.eg.db,seed = 1116,keyType = "ENTREZID",pvalueCutoff = 0.1)
gseGO_I <- gseGO(gsea_I,ont = "BP",OrgDb = org.Hs.eg.db,seed = 1116,keyType = "ENTREZID",pvalueCutoff = 0.1)
gseGO_L <- gseGO(gsea_L,ont = "BP",OrgDb = org.Hs.eg.db,seed = 1116,keyType = "ENTREZID",pvalueCutoff = 0.1)


gseGO_A.T <- simplify(gseGO_A)
gseGO_A.T <- pairwise_termsim(gseGO_A.T)
gseGO_A <- setReadable(gseGO_A,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

gseGO_A.sel <- gseGO_A@result

gseGO_A.sel$core_enrichment <- stringr::str_replace_all(gseGO_A.sel$core_enrichment,pattern = "/",replacement = ",")

gseGO_A.sel <- split(gseGO_A.sel$core_enrichment,gseGO_A.sel$ID)

gseGO_A.sel_fli <- lapply(gseGO_A.sel, function(x){return(ifelse(stringr::str_detect(x,pattern = "HLA"),NULL,x))})

gseGO_A.sel_fli <- gseGO_A.sel_fli[!is.na(gseGO_A.sel_fli)]

gseGO_A.sel <- gseGO_A.sel[which(gseGO_A.sel$ID %in% names(gseGO_A.sel_fli)),]

VlnPlot(scw,features = "NOP53")

gseKEGG_A <- gseKEGG(gsea_A,organism = "hsa",seed = 1116,keyType = "ENTREZID",pvalueCutoff = 0.05)








