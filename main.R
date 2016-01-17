library(DESeq)
source("merge.R")


get_uuids <- function(filename){
    uuids <- read.table(filename, header=T,  colClasses="character")
}
uuids <- get_uuids("uuid.txt")

#put all datatsets in one file
source <- "sample_datasets"
all_files <- list.files(source, full.names=T)
data <- merge(all_files)

#get the same datasets
selected_datasets <- intersect(uuids$file_id, colnames(data))
selected <- data[,selected_datasets]

#perform differential gene expression analysis
condition <- rep("tumor", ncol(selected))
cds <- newCountDataSet(selected, condition)
cds <- estimateSizeFactors(cds)
datanorm <- counts(cds, normalized=TRUE)
tmp <- apply(datanorm, 1, var)
w <- which(tmp != 0)
datanormSel <- datanorm[w,]

#select different colors for each project
color <- vector(mode="character", length=length(colnames(datanormSel)))
projects <- unique(uuids[,1])
rownames(uuids) <- uuids$file_id
z <- uuids[selected_datasets,]
col_count=25
for(i in 1:length(projects)){
    w <- which(grepl(projects[i], z[,1]))
    color[w] = col_count
    col_count = col_count + 25
}


pc = prcomp(t(datanormSel), scale=T)
x <- pc$x
matrix <-x[, 1:3]
plot(matrix, xlab="Principal Component 1", ylab="Principal Component 2", pch=19, cex=0.4, col=color)
