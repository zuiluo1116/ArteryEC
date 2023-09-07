library(scRNAtoolVis)
library(SeuratDisk)
library(Seurat)
library(SeuratObject)
library(tidyverse)

scw <- LoadH5Seurat("20221201.h5seurat")
load("E:/works/artery ECs/NEW.clu.RData")
scw@meta.data$`NEW.clu` <- df.3d$clusters
# DimPlot(scw,reduction = "scphere",group.by = "NEW.clu")

df.newtype <- data.frame(
  clu = scw$NEW.clu
)

df.newtype[which(df.newtype$clu %in% c(1,2,4,6,9,10)),2] <- "ACKR1+ EC"
df.newtype[which(df.newtype$clu %in% c(5,7)),2] <- "TCIM+ EC"
df.newtype[which(df.newtype$clu %in% c(3)),2] <- "ITLN1+ EC"
df.newtype[which(df.newtype$clu %in% c(8)),2] <- "Lymphatic EC"

scw$celltype <- factor(df.newtype$V2,levels = c("ITLN1+ EC","TCIM+ EC","ACKR1+ EC","Lymphatic EC"))

# DimPlot(scw,group.by = "celltype",reduction = "scphere")+NoLegend()


scw$state <- paste0(scw$GSE,"_",scw$orig.ident)
df.newstate <- data.frame(
  clu = scw$state
)
table(df.newstate$clu)
orig_states <- levels(factor(scw$state))

df.newstate[which(df.newstate$clu %in% orig_states[c(1:3,15,17,19,26:32)]),2] <- "Control"
df.newstate[is.na(df.newstate$V2),2] <- "Burden"
scw$state <- factor(df.newstate$V2,levels = c("Control","Burden"))



load("bulk.RData")

library(Scissor)
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

## GSE120521

phe_120521 <- c(rep(0,4),rep(1,4))
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

pa_1205211 <- Scissor_SCT(GSE120521, scw, phe_120521,cutoff = 0.5, tag = tag_120521,
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
        reduction = "scphere")+NoLegend()
ggsave(filename = "Figures/Fig.3B.pdf",height = 6,width = 8)

library(ggforce)
library(ggalluvial)
library(ggpubr)
df_132651 <- table(sct$GSE,sct$celltype,sct$state,sct$Sci_132651)
data <- reshape2::melt(df_132651)
# data$Var3 <- paste0("Clu_",data$Var3)
head(data)
# data <- data[which(data$Var4!="NR"),]
# data <- data[which(data$Var2!="Lymphatic EC"),]
# data <- data[which(data$Var1=="GSE159677"),]
data <- gather_set_data(data, 1:4)
data$Var4 <- factor(data$Var4,levels = c("Pos","Neg","NR"))


sankey_132651_R <- ggplot(data,
                          aes(axis1 = Var3, axis2 = Var2, axis3 = Var4,
                              y= value)) +
  scale_x_discrete(limits = c("Sample","Cell.type", "Scissor"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Var4)) +
  scale_fill_manual(values=c("#0000ff","#ff0000","grey"))+
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() +
  ggtitle("GSE132651",
          "Scissor for NL(-) or ABNL(+)")+NoLegend()
ggsave(sankey_132651_R,filename = "Figures/Fig.3E1.pdf",height = 6,width = 8)

sankey_132651_L <- ggplot(data,
                          aes(axis1 = Var3, axis2 = Var2, axis3 = Var4,
                              y= value)) +
  scale_x_discrete(limits = c("Sample","Cell.type", "Scissor"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Var2)) + 
  scale_fill_manual(values=scales::hue_pal()(4))+
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() +
  ggtitle("GSE132651",
          "Scissor for NL(-) or ABNL(+)")+NoLegend()
ggsave(sankey_132651_L,filename = "Figures/Fig.3E2.pdf",height = 6,width = 8)

DimPlot(sct,group.by = "Sci_132651",
        cols = c("#0000ff","#f0f0f0","#ff0000"),
        reduction = "scphere")+NoLegend()
ggsave(filename = "Figures/Fig.3F.pdf",height = 6,width = 8)

df_132651 <- table(sct$GSE,sct$celltype,sct$state,sct$Sci_132651)
data <- reshape2::melt(df_132651)
# data$Var3 <- paste0("Clu_",data$Var3)
head(data)
# data <- data[which(data$Var4!="NR"),]
# data <- data[which(data$Var2!="Lymphatic EC"),]
# data <- data[which(data$Var1=="GSE159677"),]
data <- gather_set_data(data, 1:4)
data$Var4 <- factor(data$Var4,levels = c("Pos","NR","Neg"))