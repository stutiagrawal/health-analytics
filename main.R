library(DESeq)
library(optparse)
source("merge.R")

get_uuids <- function(filename){
    uuids <- read.table(filename, header=T,  colClasses="character")
}

option_list <- list(
    make_option("--uuid", type="character", help="path to uuids"),
    make_option("--datasets", type="character", help="path to datasets"),
    make_option("--workdir", type="character", help="path to workdir")
)
parser <- OptionParser(usage="main.R [options] file", option_list=option_list)
args <- parse_args(parser)
print(args$uuid)
print(args$datasets)

#get association of project name and uuid
uuids <- get_uuids(args$uuid)

#select files associated with the selected datasets
desired_datasets <- vector(mode="character", length=0)
for(i in 1:length(uuids$file_id)){
    print(uuids$file_id[i])
    filename <- paste(args$datasets, uuids$file_id[i], ".htseq.counts", sep="")
    print(filename)
    if (file.exists(filename)){
        desired_datasets <- c(desired_datasets, filename)
    }
}


#put all datatsets in one file
selected <- merge(desired_datasets)

#get the same datasets
#selected_datasets <- intersect(uuids$file_id, colnames(data))
#selected <- data[,selected_datasets]

#normalize gene counts using DESeq
condition <- rep("tumor", ncol(selected))
cds <- newCountDataSet(selected, condition)
cds <- estimateSizeFactors(cds)
datanorm <- counts(cds, normalized=TRUE)
tmp <- apply(datanorm, 1, var) #get std-deviation of genes
w <- which(tmp != 0) #select genes with std-dev = 0
datanormSel <- datanorm[w,] #remove those genes as they cannot help with PCA

#select different colors for each project
color <- vector(mode="character", length=length(colnames(datanormSel)))
projects <- unique(uuids[,1])
rownames(uuids) <- uuids$file_id
z <- uuids[colnames(selected),]
col_count=25
for(i in 1:length(projects)){
    w <- which(grepl(projects[i], z[,1]))
    color[w] = col_count
    col_count = col_count + 25
}

#perform principal component analysis.
pc = prcomp(t(datanormSel), scale=T)
x <- pc$x
matrix <-x[, 1:3]

png(file=paste(args$workdir, "myplot.png", sep=""))
plot(matrix, xlab="Principal Component 1", ylab="Principal Component 2", pch=19, cex=0.4, col=color)
#dev.off()
