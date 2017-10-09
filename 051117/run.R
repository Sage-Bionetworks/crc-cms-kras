library(synapseClient)
library(data.table)
suppressPackageStartupMessages(library("hgu133plus2.db"))
synapseLogin()

## From https://github.com/chferte/KRAS_Analysis/blob/master/MEKi_prediction/MEK_framework/load_crc_tcga_data.R
combine_probes_2_gene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
    for(gene in ugenes){
    sub.expr <- as.matrix(expr[which(genes == gene),])
    if(dim(sub.expr)[2] == 1){
      M[gene,] <- sub.expr
    }else{
      tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
      tmp.c <- mean(cor(tmp, t(sub.expr)))
      #cat(gene," ", tmp.c, "\n")
      multiplier <- ifelse(tmp.c < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  M
}

doaffy <- function(synId){
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
    file <- getFileLocation(obj)
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    data.set <- fread(file, header=TRUE, fill=TRUE)
    df <- as.data.frame(data.set)
    rownames(df) <- df[,1]
    cols <- colnames(df)[1:(ncol(df)-1)]
    df <- df[,-1]
    colnames(df) <- cols
    data.set <- df
    rm(df)
    
    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    mask <- !is.na(symbol)
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    expr <- combine_probes_2_gene(data.set.m,symbol.m)
    colnames(expr) <- gsub("(.*?)_.*","\\1",colnames(expr))
    expr
}

french_expr <- doaffy("syn2363561")

synId <- "syn8533558"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
circ.genes <- read.table(file,as.is=TRUE)[,1]

genes <- circ.genes[circ.genes %in% rownames(french_expr)]
cat(paste0(length(genes), " of ", length(circ.genes), " CIRC genes in French data set:\n"))
cat(paste(genes, collapse=","))
cat("\n")

genes <- circ.genes[!(circ.genes %in% rownames(french_expr))]
cat(paste0(length(genes), " of ", length(circ.genes), " CIRC genes not in French data set:\n"))
cat(paste(genes, collapse=","))
cat("\n")

