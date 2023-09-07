##construction
library(Seurat)
# GSE189179
GSE189179 <- list.files("./GSE189179_RAW")
GSE189179.list <- list()
for (i in 1:length(GSE189179)) {
  i=4
  GSE189179.list[[i]] <- CreateSeuratObject(counts = Read10X(paste0("./GSE189179_RAW/",GSE189179[i])))
}
run <- paste0("./GSE189179_RAW/",GSE189179[i])
barcode.loc <- file.path(run, "barcodes.tsv")
cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                            sep = "\t", row.names = NULL)
if (ncol(x = cell.barcodes) > 1)
cell.names <- readLines(con = barcode.loc)
features.loc <- file.path(run, "features.tsv")
feature.names <- read.delim(file = features.loc, header = FALSE, 
                            stringsAsFactors = FALSE)
counts <- readMM(paste0("./GSE189179_RAW/",GSE189179[i],"/matrix.mtx"))
counts@Dimnames[[1]] <- feature.names$V2
counts@Dimnames[[2]] <- cell.barcodes$V1
GSE189179.list[[4]] <- CreateSeuratObject(counts = counts) 
names(GSE189179.list) <- GSE189179
save(GSE189179.list,file = "GSE189179.RData")
rm(list = (ls()))

#GSE179159
GSE179159 <- as.data.frame(list.files("./GSE179159_RAW"))
GSE179159 <- as.data.frame(stringr::str_split_fixed(GSE179159$`list.files("./GSE179159_RAW")`,pattern = "[_.]",n=3))
GSE179159$V3 <- paste0(GSE179159$V2,"_",stringr::str_split_fixed(GSE179159$V3,pattern = ".tsv",n=3)[,1])
file.rename(paste0("./GSE179159_RAW/",list.files("./GSE179159_RAW")),paste0("./GSE179159_RAW/",GSE179159$V3,".tsv"))
GSE179159k <- list.files("./GSE179159_RAW")
GSE179159.list <- list()
for (i in 1:length(GSE179159)) {
  # i=1
  tmp <- fread(paste0("./GSE179159_RAW/",GSE179159k[i]))
  counts <- as(as.matrix(tmp[,2:ncol(tmp)]),"dgCMatrix")
  counts@Dimnames[[1]] <- tmp$Gene
  # head(counts)[1:4,1:4]
  GSE179159.list[[i]] <- CreateSeuratObject(counts = counts,project = )
}
names(GSE179159.list) <- stringr::str_split_fixed(GSE179159k,pattern = ".tsv",n=2)[,1]
save(GSE179159.list,file = "GSE179159.RData")

rm(list = (ls()))

#GSE213740
GSE213740 <- as.data.frame(list.files("./GSE213740_RAW"))
GSE213740 <- as.data.frame(stringr::str_split_fixed(GSE213740$`list.files("./GSE213740_RAW")`,pattern = "[_]",n=5))
for (i in 1:length(unique(GSE213740$V1))) {
  # i=1
  filename <- paste0(GSE213740[i*3,2],"_",GSE213740[i*3,4])
  dir.create(paste0("./GSE213740_RAW/",filename))
  file.rename(paste0("./GSE213740_RAW/",list.files("./GSE213740_RAW",pattern = paste0("^",GSE213740[i*3,1]))),
              paste0(paste0("./GSE213740_RAW/",filename,"/"),c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz")))
}
GSE213740.list <- list()
dirs <- list.files("./GSE213740_RAW/")
for (i in 1:length(unique(GSE213740$V1))) {
  GSE213740.list[[i]] <- CreateSeuratObject(counts = Read10X(paste0("./GSE213740_RAW/",dirs[i])))
}
names(GSE213740.list) <- dirs[1:9]

save(GSE213740.list,file = "GSE213740.RData")

rm(list = (ls()))

#GSE216860
GSE216860 <- as.data.frame(list.files("./GSE216860_RAW"))
GSE216860 <- as.data.frame(stringr::str_split_fixed(GSE216860$`list.files("./GSE216860_RAW")`,pattern = "[_]",n=5))

for (i in 1:length(unique(GSE216860$V1))) {
  # i=1
  filename <- paste0(GSE216860[i*3,2],"_",GSE216860[i*3,4])
  dir.create(paste0("./GSE216860_RAW/",filename))
  file.rename(paste0("./GSE216860_RAW/",list.files("./GSE216860_RAW",pattern = paste0("^",GSE216860[i*3,1]))),
              paste0(paste0("./GSE216860_RAW/",filename,"/"),c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz")))
}
GSE216860.list <- list()
dirs <- list.files("./GSE216860_RAW/")
for (i in 1:length(unique(GSE216860$V1))) {
  GSE216860.list[[i]] <- CreateSeuratObject(counts = Read10X(paste0("./GSE216860_RAW/",dirs[i])))
}
names(GSE216860.list) <- dirs[1:6]

save(GSE216860.list,file = "GSE216860.RData")

scRNAlist <- list()
scRNAlist <- GSE179159.list[1:8]
names(scRNAlist)[1:8] <- paste0("GSE179159_",names(GSE179159.list)[1:8])
scRNAlist[9:12] <- GSE189179.list[1:4]
names(scRNAlist)[9:12] <- paste0("GSE189179_",names(GSE189179.list)[1:4])
scRNAlist[13:21] <- GSE213740.list[1:9]
names(scRNAlist)[13:21] <- paste0("GSE213740_",names(GSE213740.list)[1:9])
scRNAlist[22:27] <- GSE216860.list[1:6]
names(scRNAlist)[22:27] <- paste0("GSE216860_",names(GSE216860.list)[1:6])
rm(list = c("GSE179159.list","GSE189179.list","GSE213740.list","GSE216860.list"))
scRNAlist <- parallel::mclapply(scRNAlist, FUN=function(x) SCTransform(x), mc.cores = 1) 

scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000)
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features) 
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)

save(scRNAlist)
