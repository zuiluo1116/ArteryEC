# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE132651", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000001111111111111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- gset@assayData[["exprs"]]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("NL","ABNL"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)

colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
ex <- as.data.frame(ex)

tT <- topTable(fit2, adjust="fdr",number=nrow(ex))

tT <- dplyr::select(tT, c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.Symbol","Gene.Title"))

write.table(tT, file=stdout(), row.names=F, sep="\t")

ex <- ex[rownames(tT),]

GSE132651 <- cbind(ex,tT)

GSE132651$Gene.Symbol <- stringr::str_split_fixed(GSE132651$Gene.Symbol,pattern = "///",n=2)[,1]

GSE132651 <- GSE132651[!duplicated(GSE132651$Gene.Symbol),]

rownames(GSE132651) <- GSE132651$Gene.Symbol

save(GSE132651,file = "GSE132651_all.RData")
GSE132651 <- GSE132651[,1:19]
names(GSE132651) <- c(paste0(rep("NL",6),"_",1:6),paste0(rep("ABNL",13),"_",1:13))
save(GSE132651,file = "GSE132651.RData")





# 
# # Visualize and quality control test results.
# # Build histogram of P-values for all genes. Normal test
# # assumption is that most genes are not differentially expressed.
# tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
#      ylab = "Number of genes", main = "P-adj value distribution")
# 
# # summarize test results as "up", "down" or "not expressed"
# dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)
# 
# # Venn diagram of results
# vennDiagram(dT, circle.col=palette())
# 
# # create Q-Q plot for t-statistic
# t.good <- which(!is.na(fit2$F)) # filter out bad probes
# qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
# 
# # volcano plot (log P-value vs log fold change)
# colnames(fit2) # list contrast names
# ct <- 1        # choose contrast of interest
# volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
#             highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
# 
# # MD plot (log fold change vs mean log expression)
# # highlight statistically significant (p-adj < 0.05) probes
# plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
# abline(h=0)
# 
# ################################################################
# # General expression data analysis
# ex <- exprs(gset)
# 
# # box-and-whisker plot
# ord <- order(gs)  # order samples by group
# palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
#                    "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
#                    par(mar=c(7,4,2,1))
#                    title <- paste ("GSE132651", "/", annotation(gset), sep ="")
#                    boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
#                    legend("topleft", groups, fill=palette(), bty="n")
#                    
#                    # expression value distribution
#                    par(mar=c(4,4,2,1))
#                    title <- paste ("GSE132651", "/", annotation(gset), " value distribution", sep ="")
#                    plotDensities(ex, group=gs, main=title, legend ="topright")
#                    
#                    # UMAP plot (dimensionality reduction)
#                    ex <- na.omit(ex) # eliminate rows with NAs
#                    ex <- ex[!duplicated(ex), ]  # remove duplicates
#                    ump <- umap(t(ex), n_neighbors = 8, random_state = 123)
#                    par(mar=c(3,3,2,6), xpd=TRUE)
#                    plot(ump$layout, main="UMAP plot, nbrs=8", xlab="", ylab="", col=gs, pch=20, cex=1.5)
#                    legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
#                           col=1:nlevels(gs), title="Group", pt.cex=1.5)
#                    library("maptools")  # point labels without overlaps
#                    pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
#                    
#                    # mean-variance trend, helps to see if precision weights are needed
#                    plotSA(fit2, main="Mean variance trend, GSE132651")












##GSE199709
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
GSE199709 <- as.data.frame(GSE199709)

design=data.frame(case=c(0,0,0,1,1,1),
                  control=c(1,1,1,0,0,0))
ex <- scale(GSE199709)
colnames(design) <- c("LSS","Control")

fit <- lmFit(ex, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(colnames(design)[1], colnames(design)[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

ex <- as.data.frame(ex)

tT <- topTable(fit2, adjust="fdr",number=nrow(GSE199709))
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.5)


tT$change <- dT@.Data

tT <- tT[rownames(GSE199709),]



GSE118446 <- read.csv("GSE118446/GSE118446.txt",sep = "\t")
GSE118446 <- GSE118446[!duplicated(GSE118446$name),]
rownames(GSE118446) <- GSE118446$name
GSE118446 <- GSE118446[,8:25]
GSE118446 <- as.matrix(GSE118446)
GSE118446 <- GSE118446[,13:18]


design=data.frame(case=c(0,0,0,1,1,1),
                  control=c(1,1,1,0,0,0))
ex <- scale(GSE118446)
colnames(design) <- c("Control","MES")

fit <- lmFit(ex, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(colnames(design)[1], colnames(design)[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr",number=nrow(GSE118446))
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)


tT$change <- dT@.Data
tT <- tT[rownames(GSE118446),]

GSE118446_all <- cbind(GSE118446,tT)
GSE118446_all["CCL21",]




load("geneset_df.RData")
geneset_df$Cell.Type <- factor(geneset_df$Cell.Type,levels = c("ITLN1+ EC",
                                                               "TCIM+ EC",
                                                               "ACKR1+ EC",
                                                               "Lymphatic+ EC"))
plot_cells(cds = cds,genes = geneset_df)






GSE199709_all <- cbind(GSE199709,tT)
GSE199709_all$Gene.Symbol <- rownames(GSE199709_all)

TOPgenes_burden <- GSE199709_all %>% filter(P.Value<0.05 & abs(logFC)>0.58) %>% arrange(logFC)

TOPgenes_burden_high <- TOPgenes_burden[rownames(TOPgenes_burden) %in% rownames(Track_genes),] %>%
  filter(logFC>0) %>% pull(Gene.Symbol) %>% as.character()
TOPgenes_burden_low <- TOPgenes_burden[rownames(TOPgenes_burden) %in% rownames(Track_genes),] %>%
  filter(logFC<0) %>% pull(Gene.Symbol) %>% as.character()


GSE199709_geneset <- data.frame(
  "GeneSymbol" = c(TOPgenes_burden_high,TOPgenes_burden_low),
  "Group" = c(rep("High",length(TOPgenes_burden_high)),rep("Low",length(TOPgenes_burden_low)))
)

GSE199709_geneset <- data.frame(
  "GeneSymbol" = c(TOPgenes_burden_high[1:10],TOPgenes_burden_low[1:10]),
  "Group" = c(rep("High",10),rep("Low",10))
)


save(list = c("TOPgenes_burden_high","TOPgenes_burden_low"),file = "GSE199709_gene.RData")
load("GSE199709_gene.RData")

plot_cells(cds,
           genes=GSE199709_geneset,
           show_trajectory_graph=FALSE)

table(TOPgenes_burden$change)




load("E:/works/artery ECs/1213.Rdata")
load("E:/works/artery ECs/FAMS.Rdata")
FAMs_3d_scale_NEW$cluster <- as.character(FAMs_3d_scale_NEW$cluster)
FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$cluster %in% c(1,2,4,6,9,10)),6] <- "ACKR1+ EC"
FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$cluster %in% c(5,7)),6] <- "TCIM+ EC"
FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$cluster %in% c(3)),6] <- "ITLN1+ EC"
FAMs_3d_scale_NEW[which(FAMs_3d_scale_NEW$cluster %in% c(8)),6] <- "Lymphatic EC"

FAMs_sel <- FAMs_3d_scale_NEW%>% 
  group_by_("cluster") %>% 
  filter(pct.1-pct.2>0.2 & avg_log2FC>0.3) %>% 
  arrange("avg_log2FC") %>% 
  slice_max(order_by = avg_log2FC, n = 25)
geneset <- split(FAMs_sel$gene,FAMs_sel$cluster)
geneset <- geneset[c(2,4,1,3)]
names(geneset) <- paste0("Cluster ",names(geneset))

library(GSVA)
library(reshape2)
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
  geom_point()+
  facet_wrap(~Var2, scale = "free_x") +
  scale_y_continuous(name = "Clusters Score") +
  scale_x_discrete(name = "Continent") +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+theme_bw()


library(ggplot2)
library(GOplot)
data(EC)
head(EC$david)
head(EC$genelist)
circ <- circle_dat(EC$david, EC$genelist)
GOBar(subset(circ, category == 'BP'))
GOBar(circ, display = 'multiple')
# 加标题上色
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
# 分区
GOBubble(circ, labels = 3)
# 加标题, 改变环的颜色, f分区， 
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  
# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)  
# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)
#Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)
# Generate a circular visualization of selected terms
IDs <- c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886', 'GO:0008092', 'GO:0008047')
GOCircle(circ, nsub = IDs)
# Generate a circular visualization for 10 terms
GOCircle(circ, nsub = 10)
##??????##
head(EC$genes)
EC$process
chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)
# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = EC$genes)
# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = EC$process)
# Create the plot
chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(3, 0), gene.order = 'logFC')
##GOheatmap##
# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-8], nlfc = 0)
# Now we create the heatmap with logFC values and user-defined colour scale
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))
##GOclust##
GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)
GOCluster(circ, EC$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
##GOvn##
l1 <- subset(circ, term == 'heart development', c(genes,logFC))
l2 <- subset(circ, term == 'plasma membrane', c(genes,logFC))
l3 <- subset(circ, term == 'tissue morphogenesis', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('heart development', 'plasma membrane', 'tissue morphogenesis'))




