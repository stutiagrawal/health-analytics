merge <- function(all_files){

    for (i in 1:length(all_files)){
        id <- all_files[i]
        print(id)
        count <- read.table(id, colClasses="character")
        colnames(count) <- c("gene_id", "expression")
        if(i == 1){
            data <- matrix(nrow=nrow(count), ncol=length(all_files))
            columns <- vector(mode="character", length(all_files))
            rownames(data) <- count$gene_id
        }
        id <- unlist(strsplit(id, "/", fixed=T))
        id <- id[length(id)]
        id <- gsub(".htseq.counts", "", id)
        columns[i] <- id
        data[,i] <- as.numeric(count$expression)
    }   
    data <- data.frame(data)
    colnames(data) <- columns
    return(data)
}


