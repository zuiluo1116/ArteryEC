# continue after 1130.R

# ECs analysis method
# librarys

library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(harmony)
library(COSG)
library(monocle3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

scw <- LoadH5Seurat("20221201.h5seurat")
# scw$scphere.clu
Idents(scw) <- "scphere.clu"
ident_markers <- COSG::cosg(scw,groups = "all",assay = "RNA",slot = "count",n_genes_user = 100)

# 构建monocle对象
data<-GetAssayData(scw,assay ='RNA',slot ='counts')
cell_metadata <-scw@meta.data
cell_metadata[,7:29] <- NULL

gene_annotation <-data.frame(gene_short_name =rownames(data))
rownames(gene_annotation)<-rownames(data)
cdsl <- new_cell_data_set(data,cell_metadata =cell_metadata,gene_metadata =gene_annotation)
cds <- new_cell_data_set(data,cell_metadata =cell_metadata,gene_metadata =gene_annotation)
# 主成分
cds <- preprocess_cds(cds, num_dim = 50,method = "PCA")
cdsl <- preprocess_cds(cdsl, num_dim = 50,method = "PCA")
plot_pc_variance_explained(cds)
# 降维 默认UMAP
cds <- reduce_dimension(cds,preprocess_method = "PCA",reduction_method = "UMAP") #preprocess_method默认是PCA
cdsl <- reduce_dimension(cdsl,preprocess_method = "PCA",reduction_method = "UMAP") #preprocess_method默认是PCA

# 替换UMAP坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scw, reduction = "scphere")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cdsl@int_colData$reducedDims$UMAP <- int.embed
# 识别分群（基于UMAP二维）
cds <- cluster_cells(cds,reduction_method = "UMAP")
cdsl <- cluster_cells(cdsl,reduction_method = "UMAP")
cds <- learn_graph(cds)
cdsl <- learn_graph(cdsl)

plot_cells(cds, reduction_method="UMAP")
plot_cells(cdsl, reduction_method="UMAP")

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "X3d.clu", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=T,group_label_size = 12)

Track_genesl <- graph_test(cdsl, neighbor_graph="principal_graph")

Track_genes <- graph_test(cds, neighbor_graph="principal_graph")
Track_genes.fli <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3) %>% arrange(desc(morans_I))
Track_genes_sig100 <- Track_genes.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Idents(scw)

DimPlot(scw,group.by = "scphere.clu",split.by = "scissor_LSS_split",ncol = 6,reduction = "scphere")

FAMs <- FindAllMarkers(scw,only.pos = T)
intersect(Track_genes_sig100,FAMs$gene)


save(list = c("Track_genes","Track_genes_convert","Track_genes.fli","Track_genes_convert.fli"),
     file = "track.RData")

genelist <- pull(Track_genes, gene_short_name) %>% as.character()
genelist_convert <- pull(Track_genes_convert, gene_short_name) %>% as.character()


gene_module <- find_gene_modules(cds[genelist,], resolution=5e-3, cores = 1)
gene_module1 <- find_gene_modules(cds[genelist,], resolution=5e-2, cores = 1)

gene_module_convert <- find_gene_modules(cds[genelist_convert,], resolution=5e-3, cores = 1)

save(list = c("genelist","genelist_convert","gene_module","gene_module_convert"),
     file = "gene.module.RData")

Track_match <- data.frame(
  "Track_genes_sig100" = Track_genes_sig100
)

Track_match[["cluster"]] = FAMs[match(Track_genes_sig100,FAMs$gene),"cluster"]
Track_match[["merge.cluster"]] = FAMs_merge[match(Track_genes_sig100,FAMs_merge$gene),"cluster"]
# Track_match[["module"]] = as.data.frame(gene_module[match(Track_genes_sig100,gene_module$id),"module"])[,1]
Track_match[["module"]] = as.data.frame(gene_module[match(Track_genes_sig100,gene_module$id),"module"])[,1]
write.csv(Track_match,"Track_match.csv")

FeaturePlot(scw,"RGCC",reduction = "scphere")


tmps <- cbind(table(Track_match$cluster,Track_match$module))
tmps <- tmps[,order(colSums(tmps),decreasing = T)]
tmps <- tmps[,which(colSums(tmps)!=0)]
levels(as.factor(scw$GSE))

colnames(tmps)
# [1] "23"  "3"   "8"   "82"  "26"  "19"  "48"  "51"  "81"  "13"  "25"  "34"  "52"  "69"  "72"  "77" 
# [17] "84"  "94"  "2"   "7"   "9"   "18"  "32"  "40"  "60"  "61"  "79"  "103"
plot_cells(cds,
           genes=gene_module %>% filter(module %in% c(36,5,24,8)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)
module.292t <- unlist(c(gene_module1[which(gene_module1$module == 292),1]))

module.23 <- unlist(c(gene_module[which(gene_module$module == 23),1]))
module.3 <- unlist(c(gene_module[which(gene_module$module == 3),1]))
module.8 <- unlist(c(gene_module[which(gene_module$module == 8),1]))
module.82 <- unlist(c(gene_module[which(gene_module$module == 82),1]))

module.3 <- bitr(module.3,fromType = "SYMBOL",
                 toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)
module.3 <- na.omit(module.3)
module.3 <- module.3[!duplicated(module.3$SYMBOL),]
module3.GO <- enrichGO(gene = module.3$ENSEMBL,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       pvalueCutoff = 0.05)
module3.GO <- simplify(module3.GO)
module3.GO <- pairwise_termsim(module3.GO)
module3.GO <- setReadable(module3.GO,OrgDb = org.Hs.eg.db)

module.23 <- bitr(module.23,fromType = "SYMBOL",
                 toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)
module.23 <- na.omit(module.23)
module.23 <- module.23[!duplicated(module.23$SYMBOL),]
module23.GO <- enrichGO(gene = module.23$ENSEMBL,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       pvalueCutoff = 0.5)
module23.GO <- simplify(module23.GO)
module23.GO <- pairwise_termsim(module23.GO)
module23.GO <- setReadable(module23.GO,OrgDb = org.Hs.eg.db)

module.8 <- bitr(module.8,fromType = "SYMBOL",
                  toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)
module.8 <- na.omit(module.8)
module.8 <- module.8[!duplicated(module.8$SYMBOL),]
module8.GO <- enrichGO(gene = module.8$ENSEMBL,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        pvalueCutoff = 0.5)
module8.GO <- simplify(module8.GO)
module8.GO <- pairwise_termsim(module8.GO)
module8.GO <- setReadable(module8.GO,OrgDb = org.Hs.eg.db)

module.82 <- bitr(module.82,fromType = "SYMBOL",
                 toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)
module.82 <- na.omit(module.82)
module.82 <- module.82[!duplicated(module.82$SYMBOL),]
module82.GO <- enrichGO(gene = module.82$ENSEMBL,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       pvalueCutoff = 0.5)
module82.GO <- simplify(module82.GO)
module82.GO <- pairwise_termsim(module82.GO)
module82.GO <- setReadable(module82.GO,OrgDb = org.Hs.eg.db)


################################################################################

library(GO.db)
library(org.Hs.eg.db)
ls("package:org.Hs.eg.db")
# [1] "GO"            "GO.db"         "GO_dbconn"     "GO_dbfile"     "GO_dbInfo"     "GO_dbschema"  
# "GOBPANCESTOR"  "GOBPCHILDREN"  "GOBPOFFSPRING" "GOBPPARENTS"   
# "GOCCANCESTOR"  "GOCCCHILDREN"  "GOCCOFFSPRING" "GOCCPARENTS"   
# "GOMAPCOUNTS"   
# "GOMFANCESTOR"  "GOMFCHILDREN"  "GOMFOFFSPRING" "GOMFPARENTS" 
# "GOOBSOLETE"    "GOSYNONYM"     "GOTERM"     

go.category <- as.list(GOTERM)
go.description <- as.list(GOTERM)
go.trans <- as.list(org.Hs.egGO2ALLEGS)
go.allterms <- as.list(org.Hs.egSYMBOL)
go.allterms <- as.data.frame(cbind(go.allterms))
go.allterms$eg <- rownames(go.allterms)

exact <- lapply(go.trans, find_out,a = go.allterms, b = 2, e = 1)
exact_un <- plyr::ldply(exact, data.frame)
names(exact_un) <- c("GO","Symbol")
exact_k <- split(exact_un$GO,exact_un$Symbol)

go.content <- exact_k


go.category <- lapply(go.category, function(x){
  tmp <- x@Term
  return(tmp)
})
go.category <- plyr::ldply(go.category, data.frame)
names(go.category) <- c("GO","Term")

exact.length <- as.data.frame(cbind(lapply(exact, function(x){return(length(x[]))})))
exact.length$GO <- rownames(exact.length)
exact_length <- exact_length %>% arrange(desc(as.numeric(V1)))

go.category$length <- find_out(exact.length,2,go.category$GO,1)

exact_simp <- exact_k
exact_simp <- lapply(exact_simp, find_out,a = go.category, b = 1, e = 2)


# 删除重复
exact_simp <- lapply(exact_simp, function(x){
  return(x <- x[!duplicated(x)])
})
# 删除目录条目
exact_simp <- lapply(exact_simp, function(x){
  return(x[-which(x %in% c("biological_process","cellular_component","molecular_function"))])
})
# 转变为","分割的一个字符串
exact_simp <- lapply(exact_simp, function(x){
  return(paste0(x,collapse = ","))
})

exact_sel <- exact_simp[Track_genes_sig100]
exact_sel <- plyr::ldply(exact_sel, data.frame)
names(exact_sel) <- c("Symbol","Term")
exact_sel <- split(exact_sel$Symbol,exact_sel$Term)
exact_sel_length <- as.data.frame(unlist(go.length[names(exact_sel),]))


exact_sel_length$Term <- names(exact_sel) 


module_p.GO <- enrichGO(gene = Track_genes_sig100,keyType = "SYMBOL",
                        )
sig100.con <- bitr(Track_genes_sig100,fromType = "SYMBOL",
                  toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)
sig100.con <- na.omit(sig100.con)
sig100.con <- sig100.con[!duplicated(sig100.con$SYMBOL),]
module_p.GO <- enrichGO(gene = sig100.con$ENSEMBL,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        pvalueCutoff = 1)
module_p.GO <- simplify(module_p.GO,cutoff=0.4)
module_p.GO <- pairwise_termsim(module_p.GO)
module_p.GO <- setReadable(module_p.GO,OrgDb = org.Hs.eg.db)

names(exact_sel)
exact_sel_tmp <- unlist(lapply(exact_sel, function(x){return(length(x))}))
names(exact_sel_length)[1] <- "total"

exact_sel_length$num <- unlist(lapply(exact_sel, function(x){return(length(x))}))

exact_sel_length$ratio <- exact_sel_length$num/exact_sel_length$`unlist(go.length[names(exact_sel), ])`
names(exact_sel_length)
exact_sel_length <- exact_sel_length %>% arrange(desc(ratio),total) %>% filter(num > 4)

exact_sel2 <- exact_sel[c("GO:0008009","GO:0005201","GO:0001525","GO:0045446")]
names(exact_sel2) <- find_out(go.category,1,names(exact_sel2),2)
intersect(exact_sel2$`endothelial cell differentiation`,exact_sel2$angiogenesis)

plot_genes_in_pseudotime(cds["HEY1",], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 5)
Nebulosa::plot_density(scw,features = "HEY1",reduction = "scphere")
exist <- read.csv("exist.csv",header = F)
nonexist <- setdiff(Track_genes_sig100,exist$V1)
write.csv(nonexist,"nonexist.csv")
query <- data.frame("genes" = Track_genes_sig100)
query <- Topp.gene(query,"genes","Topp")


exact_selk <- plyr::ldply(exact_sel2, data.frame)

write.csv(exact_selk,"exact_selk.csv")

geneterm <- read.csv("genelist.csv")
geneterm <- split(geneterm$genes,geneterm$description)

# str_split_fixed(x,pattern = ";",n=(str_count(x,pattern = ";")+1))

geneterms <- lapply(geneterm, function(x){return(as.character(str_split_fixed(x,pattern = ";",n=(str_count(x,pattern = ";")+1))))})
geneterms <- plyr::ldply(geneterms, data.frame)
class(geneterm)
geneterm <- as.data.frame(geneterms)
names(geneterms) <- c("Terms","Symbol")
geneterms$Term <- str_split_fixed(geneterms$Terms,pattern = "/",n=2)[,2]
geneterms$Database <- str_split_fixed(geneterms$Terms,pattern = "/",n=2)[,1]
geneterms$Terms <- NULL 
string_dict <- read.csv("string_dict.csv")
geneterms$ID <- find_out(string_dict,7,geneterms$Symbol,6)
write.csv(geneterms,"terms.csv")

pathway_sel1 <- read.csv("PATHWAYs.csv")
pathway_sel <- split(pathway_sel$matching.proteins.in.your.network..labels.,pathway_sel$term.description)
pathway_sel <- lapply(pathway_sel, function(x){return(as.character(str_split_fixed(x,pattern = ",",n=(str_count(x,pattern = ",")+1))))})
pathway_sel <- plyr::ldply(pathway_sel, data.frame)
pathway_sel$Term <- find_out(pathway_sel1,4,pathway_sel$.id,2)
names(pathway_sel) <- c("ID","Symbol","Term")
write.csv(pathway_sel,"pathway_sel.csv")


exist <- read.csv("exist.csv",header = F)
nonexist <- setdiff(Track_genes_sig100[1:50],exist$V1)
write.csv(nonexist,"nonexist.csv")

################################################################################
Idents(scw) <- "X3d.clu"
expr.clu <- AverageExpression(scw,slot = "count",assay = "RNA")[[1]]
expr.clu <- expr.clu[Track_genes_sig100[1:50],]
cor_X3d.clu <- cor(expr.clu, method = c("spearman"))
library(corrplot)
corrplot(cor_X3d.clu, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

ComplexHeatmap::Heatmap(cor_X3d.clu)




# bulkscore
library(GSVA)

DimPlot(scw,group.by = "X3d.clu",reduction = "scphere")
pattern_group <- data.frame("orig.clu"=scw$X3d.clu)
pattern_group$`merge.clu` <- ifelse(pattern_group$orig.clu == 1,"1",
                                    ifelse(pattern_group$orig.clu == 2,"2","3"))
scw$`merge.clu` <- factor(pattern_group$`merge.clu`,levels = 1:3)
class(scw$`merge.clu`)
Idents(scw) <- "merge.clu"
FAMs_merge <- FindAllMarkers(scw,only.pos = T,logfc.threshold = 0.4)


FAMs_fli <- FindAllMarkers(scw,only.pos = T)
FAMs_fli$cluster <- as.integer(FAMs_fli$cluster)
FAMs_fli <- FAMs_fli %>%  arrange("cluster")
FAMs_sel <- FAMs_fli %>% 
  group_by_("cluster") %>% 
  filter(pct.1-pct.2>0.2 & avg_log2FC>0.4) %>% 
  arrange("avg_log2FC") %>% 
  top_n(25)

setdiff(Track_genes_sig100,FAMs_sel$gene)






GSE167024 <- read.csv("./GSE167024/count.csv")
GSE167024 <- GSE167024[!duplicated(GSE167024$Symbol),]
rownames(GSE167024) <- GSE167024$Symbol
GSE167024 <- as.matrix(GSE167024[,3:6])


GSE104140 <- read.csv("./GSE104140/gene_count_matrix.csv")
GSE104140 <- GSE104140[!duplicated(GSE104140$gene_symbol),]
rownames(GSE104140) <- GSE104140$gene_symbol
GSE104140 <- as.matrix(GSE104140[,-c(1:2)])
GSE104140_dict <- read.csv("./GSE104140/check.csv")
# rm(GSE104140_dict)
colnames(GSE104140) <- find_out(GSE104140_dict,2,colnames(GSE104140),11)
GSE104140 <- as.matrix(GSE104140)

class(GSE104140)
class(GSE210522)
library(reshape2)
geneset <- split(FAMs_sel$gene,FAMs_sel$cluster)
es1 <- gsva(GSE120521, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es2 <- gsva(GSE199709, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es3 <- gsva(GSE210522_sel, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es4 <- gsva(GSE211662, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es5 <- gsva(GSE167024, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)
es6 <- gsva(GSE104140, geneset,verbose=TRUE,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE)



es1.data <- melt(t(as.data.frame(es1)))
es1.data$Var1 <- str_split_fixed(es1.data$Var1,pattern = "_",n=2)[,1]

es2.data <- melt(t(as.data.frame(es2)))
es2.data$Var1 <- str_split_fixed(es2.data$Var1,pattern = "_",n=2)[,1]

es3.data <- melt(t(as.data.frame(es3)))
es3.data$Var1 <- str_split_fixed(es3.data$Var1,pattern = "_",n=3)[,2]

es4.data <- melt(t(as.data.frame(es4)))
es4.data$Var1 <- str_split_fixed(es4.data$Var1,pattern = "_",n=2)[,1]

es5.data <- melt(t(as.data.frame(es5)))
es5.data$Var1 <- str_split_fixed(es5.data$Var1,pattern = "_",n=2)[,1]

es6.data <- melt(t(as.data.frame(es6)))
es6.data$Var1 <- str_split_fixed(es6.data$Var1,pattern = "_",n=2)[,1]

ggplot(es6.data, aes(x=Var1, y=value, fill = Var2)) + 
  geom_boxplot()+
  # geom_violin()+
  # geom_jitter(stroke = 0.1,width = 0.05)+
  # geom_point()+
  facet_wrap(~Var2, scale = "free_x") +
  scale_y_continuous(name = "Clusters Score") +
  scale_x_discrete(name = "Continent") +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))


ssgsea<- gsva(dat, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)













################################################################################




library(UpSetR)
upset(fromList(exact_sel2),order.by = c("freq"),decreasing = TRUE,mb.ratio = c(0.3,0.7))

stringr::str_count(geneterm$`Activation of immune response`,pattern = "'|'")
write.csv(Track_genes_sig100,"genes.csv",sep = ",")


exact_k_tmp <- exact_simp[1:10]


exact_k[["FAU"]]
Track_match



paste0(exact_simp)
c("go.category","exact_k")

tmptest <- module.3[1:10,]
tmptest <- split(tmptest$ENTREZID,tmptest$SYMBOL)
exact_simp <- exact[1:10]
# tmptest <- unsplit(tmptest,names(tmptest))
# data.table::rbindlist(tmptest)
df <- plyr::ldply (exact_simp, data.frame)





exact_un <- unsplit(exact,f = names(exact))
head(exact_un)[1:10]
find_out <- function(a,b,c,e){

  tmp <- a[match(c,a[,b]),e]
  tmp <- as.character(tmp)
  return(tmp)
}
find_out(go.allterms,2,"219736",1)


gene_module.list <- gene_module
gene_module.list$name <- paste("Module",gene_module.list$module,sep = " ")


gene_dict <- bitr(scw@assays[["RNA"]]@data@Dimnames[[1]],fromType = "SYMBOL",
                  toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db,drop = F)

gene_module.list$ENTREZID <- gene_dict[match(gene_module.list$id,gene_dict$SYMBOL),2]
gene_module.check <- na.omit(gene_module.list)
module.list <- split(gene_module.check$ENTREZID,gene_module.check$name)

module.sel <- module.list[c("Module 23","Module 3","Module 8")]
# module3.KEGG <- enrichKEGG(gene = module.3$ENTREZID,
#                            organism="hsa",
#                            pvalueCutoff = 0.05)

tmp.com <- compareCluster(module.sel,fun = "enrichKEGG",organism="hsa", pvalueCutoff=1)
tmp.com <- compareCluster(module.sel,fun = "enrichGO",OrgDb = org.Hs.eg.db, pvalueCutoff=1)
tmp.com <- simplify(tmp.com)
tmp.com <- pairwise_termsim(tmp.com)
tmp.com <- setReadable(tmp.com,OrgDb = org.Hs.eg.db)

dotplot(tmp.com,showCategory=3,includeAll=T)


DimPlot(scw,reduction = "umap",group.by = "seurat_clusters")
DimPlot(scw,reduction = "umap",group.by = "scissor",cols = c("gray","blue","red"))
DimPlot(scw,reduction = "scphere",group.by = "scissor",cols = c("gray","blue","red"))


scw1 <- scw
DefaultAssay(scw1) <- "SCT"
scw1 <- RunUMAP(scw1,min.dist = 2,n.neighbors = 200,reduction = "harmony",dims = 1:30,assay = "SCT")
DimPlot(scw1,reduction = "umap",split.by = "GSE")
DimPlot(scw1,reduction = "umap")

################################################################################
# loadanno
Anno <- read.delim("Anno.txt")
names(Anno)
Track_match$Gene.Ontology.terms <- Anno[match(Track_match$Track_genes_sig100,Anno$Ensembl.gene),"Gene.Ontology.terms"]



#坐标修改
max(int.embed[,1]) #180
min(int.embed[,1]) #-180

tmp <- as.data.frame(int.embed)
tmp[,1] <- ifelse(tmp[,1]>100,tmp[,1]-280,tmp[,1]+80)#拼凑clu.1
tmp[,1] <- ifelse(tmp[,1]>0,tmp[,1]-180,tmp[,1]+180)



###构建翻转坐标的轨迹
cds_convert <- cds
cds_convert@int_colData$reducedDims$UMAP <- as.matrix(tmp)
cds_convert <- cluster_cells(cds_convert,reduction_method = "UMAP")
cds_convert <- learn_graph(cds_convert)

plot_cells(cds_convert, color_cells_by = "X3d.clu", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=T,group_label_size = 12)

cds_convert <- order_cells(cds_convert)

Track_genes_convert <- graph_test(cds_convert, neighbor_graph="principal_graph")
Track_genes_convert.fli <- Track_genes_convert[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3) %>% arrange(desc(morans_I))
Track_genes_convert_sig100 <- Track_genes_convert.fli %>% top_n(n=100, morans_I) %>%
  pull(gene_short_name) %>% as.character()

comp.genes <- data.frame(
  "orig" = Track_genes_sig100,
  "conv" = Track_genes_convert_sig100
)


setdiff(comp.genes$orig,comp.genes$conv)



################################################################################

#Scissor

for (i in 1:length(GSEs)) {
  tmp.sc <- subset(scw,GSE==GSEs[i])
  cellnames <- tmp.sc@assays[["RNA"]]@data@Dimnames[[2]]
  matrixsel <-total[cellnames,cellnames] 
  tmp.sc@graphs[["SCT_snn"]] <- matrixsel
  infosx <- Scissor_SCT(bulk_dataset, tmp.sc, phenotype, alpha = 0.05, tag = tag,
                        family = "binomial", #二分类
                        Save_file = paste0("./Scissor/",GSEs[i],"_.Rdata"))
  infosy <- Scissor_SCT(bulk_dataset, tmp.sc, phenotype, alpha = NULL, cutoff = 0.1, 
                        family = "binomial", Load_file = paste0("./Scissor/",GSEs[i],"_.Rdata"))
  US_lists_A[[i]] <- infosx
  US_lists_B[[i]] <- infosy
  US_lists_pos <- c(US_lists_pos,infosy$Scissor_pos)
  US_lists_neg <- c(US_lists_neg,infosy$Scissor_neg)
}












