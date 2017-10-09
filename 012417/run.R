suppressPackageStartupMessages(library("annotate"))
suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("hgu133plus2.db"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gtable"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("maxstat"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("BayesFactor"))
suppressPackageStartupMessages(library("ggbeeswarm"))

# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
num.processes <- 1
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  num.processes <- num.cores - 1
  cat(paste("Registering ", num.processes, " cores.\n", sep=""))
  registerDoMC(cores=num.processes)
}

suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("gplots"))

cms2mt.vs.wt.col <- "CMS2MT.vs.WT"
cms1mt.vs.cms2mt.col <- "CMS1MT.vs.CMS2MT"
cms3mt.vs.cms2mt.col <- "CMS3MT.vs.CMS2MT"
cms4mt.vs.cms2mt.col <- "CMS4MT.vs.CMS2MT"
nolblmt.vs.cms2mt.col <- "NOLBL.vs.CMS2MT"
kras.col <- "KRAS"

kras.col <- "krasmt"
cms2.col <- "cms2"
cms3.col <- "cms3"
cms4.col <- "cms4"

cms2.wt.vs.cms1.wt.col <- "cms2wt.vs.cms1wt"
cms3.wt.vs.cms1.wt.col <- "cms3wt.vs.cms1wt"
cms4.wt.vs.cms1.wt.col <- "cms4wt.vs.cms1wt"
cms2.mt.vs.cms2.wt.col <- "cms2mt.vs.cms2wt"
cms3.mt.vs.cms3.wt.col <- "cms3mt.vs.cms3wt"
cms4.mt.vs.cms4.wt.col <- "cms4mt.vs.cms4wt"

cms2.wt.vs.col.prefix <- "cms2wt.vs."
cms3.wt.vs.col.prefix <- "cms3wt.vs."
cms4.wt.vs.col.prefix <- "cms4wt.vs."
cms2.mt.vs.col.prefix <- "cms2mt.vs."
cms3.mt.vs.col.prefix <- "cms3mt.vs."
cms4.mt.vs.col.prefix <- "cms4mt.vs."

## Immune sets that Justin looked at.
immune.sets <- c("IMMUNE_ESTIMATE", "IMMUNE_RESP_GO_BP", "PD1_REACTOME", "IMMUNE_CD8MACRO_GALON", "IMMUNE_TH1_GALON", "IMMUNE_NKC_BREAST", "IMMUNE_THF_BREAST", "IMMUNE_TH17_GOUNARI", "IMMUNE_TREG_LUCAS", "IMMUNE_MDSC_ALBELDA")

system("mkdir output/")


## Optimize lambda and alpha
fit.elastic.net <- function(X, y) {
  N <- nrow(X)
  nfolds <- 5
  foldid <- sample(rep(seq(nfolds), length = N))
  cv.list <- list()
  alphas <- seq(from = 0, to = 1, by = 0.05)
  for(i in 1:length(alphas)) {
    cat(paste0("Fitting with alpha = ", alphas[i], "\n", sep=""))
    cv <- cv.glmnet(X, y, family = "binomial", type.measure = "auc", foldid = foldid, nfolds = nfolds, alpha = alphas[i], standardize = TRUE, parallel = TRUE)
    cv.list[[i]] <- cv
  }
  cv.list
}
# cvl <- fit.elastic.net(t(expr.m), as.factor(clin.m$kras))

##  fit elastic net for a given alpha.
##  NB: foldid argument guarantees the same folds are used across invocations of this function.
##      this allows us to use the same folds to evaluate different alphas.
##  Further: cv.glmnet randomly selects folds, so setting the seed is necessary to ensure comparability
##           across alpha values.
##  And, like, one more:  glmnet's standardize probably doesn't mean what you think it does.  You should standardize x
##                        explicitly if you need to (and be certain to do so before predicting as well).
##                        In particular, ?glmnet reports that 'The coefficients are always returned on the original scale.'
fit.elastic.net.alpha <- function(alpha, x, y, foldid, nfolds, seed = 1234) {
    set.seed(seed)
    N <- nrow(x)
    cv <- cv.glmnet(x, y, family = "binomial", type.measure = "auc", foldid = foldid, nfolds = nfolds, alpha = alpha, standardize = FALSE)    
#    cv <- cv.glmnet(x, y, family = "binomial", type.measure = "auc", foldid = foldid, nfolds = nfolds, alpha = alpha, standardize = TRUE)
    cv
}

cat("Logging in\n")
synapseLogin("brian.white")
cat("Logged in\n")

## Read in the clinical annotation data

cat("Reading clinical data\n")
obj <- synGet(id="syn2527101", downloadFile = TRUE, downloadLocation = ".")
clin.data.file <- getFileLocation(obj)
clin <- read.table(clin.data.file, sep=",", header=TRUE, as.is=TRUE)

## Read in the CMS clusters
cat("Reading CMS results\n")
obj <- synGet(id="syn4978511", downloadFile = TRUE, downloadLocation = ".")
cms.file <- getFileLocation(obj)
cms <- read.table(cms.file, header=TRUE, as.is=TRUE)

## CMS clusters that include private
obj <- synGet(id="syn2754865", downloadFile = TRUE, downloadLocation = ".")
cms.all.file <- getFileLocation(obj)
cms.all <- read.table(cms.all.file, header=TRUE, as.is=TRUE)
for(prefix in c("amc_ajccii.", "tcga_rnaseqAll.", "agendia_gse42284.", "agendia_vhb70.", "agendia_ico208.", "french.", "petacc.", "kfsyscc.", "mdanderson.")) {
    rownames(cms.all) <- gsub(pattern=prefix, replacement="", rownames(cms.all))
}
cms.all$sample <- rownames(cms.all)
cms.all$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples <- unlist(apply(cms.all[,c("CMS", "classifier")], 1, function(row) ifelse(!is.na(row[1]) & (row[1] != "UNK"), row[1], row[2])))
cms.all$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples <- unlist(lapply(cms.all$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples, function(x) ifelse(is.na(x), "NOLBL", x)))

## Check that the results from cms and cms.all are consistent for the overlapping data sets.
m <- merge(x=cms, y=cms.all, by="sample", suffixes=c(".public", ".all"))

fl <- m$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples.all != m$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples.public
any(fl)

## Read in the gene-set definitions used in Justin's Nat Med paper
obj <- synGet(id="syn2321865", downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
gene.sets <- GSA.read.gmt(file)

bindeaTbl <- read.table("input/bindea_immune.txt",sep='\t',header=TRUE)

## Restrict to a few cell types of interest
bindeaTbl <- bindeaTbl[bindeaTbl$CellType %in% c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells"),]

## Match the clinical annotations and the CMS clusters
idxs <- match(clin$sample, cms$sample)
clin <- clin[!is.na(idxs),]
clin$cms_label <- cms$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples[na.omit(idxs)]
not.cms2 <- clin$cms_label != "CMS2"

cat("Reading clinical data\n")
obj <- synGet(id="syn2754862", downloadFile = TRUE, downloadLocation = ".")
clin.all.file <- getFileLocation(obj)
clin.all <- read.table(clin.all.file, sep=",", header=TRUE, as.is=TRUE)
rownames(clin.all) <- clin.all$sample

idxs <- match(clin.all$sample, cms.all$sample)
clin.all <- clin.all[!is.na(idxs),]
clin.all$cms_label <- cms.all$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples[na.omit(idxs)]

do.coarse.grained.analysis <- FALSE
if(do.coarse.grained.analysis) {
    stop("not.cms2 is not defined if you are using clin.all (i.e., for agenda)")
    clin$cms_label[not.cms2] <- "CMSOther"
}

# From https://github.com/chferte/KRAS_Analysis/blob/master/MEKi_prediction/MEK_framework/load_crc_tcga_data.R
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

# wget https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-101/A-AFFY-101.adf.txt

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

entrez.to.symbol.mapping <- function(symbols) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), 
              filters = 'entrezgene',
              values = symbols, 
              mart = mart)
  names(bm) <- c("SYM", "INDEX")
  bm <- bm[!(bm$INDEX %in% c("")),]
  bm
}

do.custom.affy <- function(synId){
    ##    trans <- read.table("A-AFFY-101.adf.probe.gene.txt", sep="\t", header=TRUE)
    if(FALSE) {
        trans <- as.data.frame(read.table("A-AFFY-101.adf.probe.gene.txt", header=TRUE, fill=TRUE))
        trans <- trans[trans[,2] != "",]
        trans <- trans[!duplicated(trans[,1], fromLast=TRUE) & !duplicated(trans[,1], fromLast=FALSE),]
        rownames(trans) <- trans[,1]
    }
    
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
    file <- getFileLocation(obj)
    ## file <- "PETACC3_expression.tsv"
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    print(file)
    data.set <- fread(file, header=TRUE, fill=TRUE)
    df <- as.data.frame(data.set)
    rownames(df) <- df[,1]
    df <- df[,-1]
    data.set <- df
    rm(df)
    
    ## Array annotation file (this is from probe to entrez)
    annotationSynId <- "syn2199825"
    obj <- synGet(id=annotationSynId, downloadFile = TRUE, downloadLocation = ".")
    file <- getFileLocation(obj)
    trans <- read.table(file, sep="\t", header=TRUE, comment.char="")
    rownames(trans) <- trans[,1]

    trans <- trans[trans[,2] != "---",]
    trans <- trans[rownames(trans) %in% rownames(data.set),]
    
    ## Translate entrez to symbol
    sym.tbl <- entrez.to.symbol.mapping(trans[,2])
    trans <- merge(trans, sym.tbl, by.x="Entrez.GeneID", by.y="INDEX")
    dup.flag <- duplicated(trans$ProbesetID, fromLast=TRUE) | duplicated(trans$ProbesetID, fromLast=FALSE)
    trans <- trans[!dup.flag,]
    rownames(trans) <- trans$ProbesetID
    
    ##    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    symbol <- trans[rownames(data.set),"SYM"]
    mask <- !is.na(symbol) & (symbol != "")
    print(length(which(mask)))
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    print(head(data.set.m))
    expr <- combine_probes_2_gene(data.set.m,symbol.m)
    print(head(expr))
    colnames(expr) <- gsub("(.*?)_.*","\\1",colnames(expr))
    expr
}


do.discover.print.19742 <- function(synId) {
    ## Get the DiscoverPrint 19742 probe to gene translation table
    obj <- synGet(id="syn2192791", downloadFile = TRUE, downloadLocation = ".")
    zip.file <- getFileLocation(obj)    
    file <- gsub(pattern=".gz", replacement="", zip.file)
    if(!file.exists(file)) {
        system(paste0("gunzip ", zip.file))
    }
    trans <- read.table(file, sep="\t", header=TRUE)
    trans <- trans[!duplicated(trans$probe_id, fromLast=TRUE) & !duplicated(trans$probe_id, fromLast=FALSE),]
    rownames(trans) <- trans$probe_id
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
    zip.file <- getFileLocation(obj)
    file <- gsub(pattern=".zip", replacement=".txt", zip.file)
    if(!file.exists(file)) {
        system(paste0("unzip ", zip.file))
    }
    ## file <- "PETACC3_expression.tsv"
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    data.set <- fread(file, header=TRUE, fill=TRUE)
    df <- as.data.frame(data.set)
    df <- df[!duplicated(df[,1], fromLast=TRUE) & !duplicated(df[,1], fromLast=FALSE),]
    rownames(df) <- df[,1]
    df <- df[,-1]
    na.flag <- unlist(apply(df, 1, function(row) any(is.na(row))))
    data.set <- df[!na.flag,]
    rm(df)
    
    ##    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    symbol <- trans[rownames(data.set),"symbol"]
    mask <- !is.na(symbol)
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    expr <- combine_probes_2_gene(data.set.m,symbol.m)
    colnames(expr) <- gsub("(.*?)_.*","\\1",colnames(expr))
    expr
}

do.discover.print.32627 <- function(synId) {
    ## Get the DiscoverPrint 32627 probe to gene translation table
    obj <- synGet(id="syn2192801", downloadFile = TRUE, downloadLocation = ".")
    file <- getFileLocation(obj)    
    trans <- read.table(file, sep="\t", header=TRUE)
    trans <- trans[!duplicated(trans$ProbeName, fromLast=TRUE) & !duplicated(trans$ProbeName, fromLast=FALSE),]
    rownames(trans) <- trans$ProbeName
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
    zip.file <- getFileLocation(obj)
    file <- gsub(pattern=".zip", replacement=".txt", zip.file)
    if(!file.exists(file)) {
        system(paste0("unzip ", zip.file))
    }
    ## file <- "PETACC3_expression.tsv"
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    data.set <- fread(file, header=TRUE, fill=TRUE)
    df <- as.data.frame(data.set)
    df <- df[!duplicated(df[,1], fromLast=TRUE) & !duplicated(df[,1], fromLast=FALSE),]
    rownames(df) <- df[,1]
    df <- df[,-1]
    na.flag <- unlist(apply(df, 1, function(row) any(is.na(row))))
    data.set <- df[!na.flag,]
    rm(df)
    
    ##    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    symbol <- trans[rownames(data.set),"GeneName"]
    mask <- !is.na(symbol)
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    expr <- combine_probes_2_gene(data.set.m,symbol.m)
    colnames(expr) <- gsub("(.*?)_.*","\\1",colnames(expr))
    expr
}

do.discover.print.mdacc <- function(synId) {
    ## Get the DiscoverPrint 32627 probe to gene translation table
    ## This is the same 32627 file as above using this synId
    ## obj <- synGet(id="syn2192801", downloadFile = TRUE, downloadLocation = ".")    
    obj <- synGet(id="syn2233216", downloadFile = TRUE, downloadLocation = ".")
    file <- getFileLocation(obj)    
    trans <- read.table(file, sep="\t", header=TRUE)
    trans <- trans[!duplicated(trans$ProbeName, fromLast=TRUE) & !duplicated(trans$ProbeName, fromLast=FALSE),]
    rownames(trans) <- trans$ProbeName
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
    zip.file <- getFileLocation(obj)
    file <- zip.file
    if(grepl(x=zip.file, "zip")) {
        file <- gsub(pattern=".zip", replacement=".txt", zip.file)
        if(!file.exists(file)) {
            system(paste0("unzip ", zip.file))
        }
    } else if(grepl(x=zip.file, "gz")) {
        file <- gsub(pattern=".zip", replacement=".txt", zip.file)
        if(!file.exists(file)) {
            system(paste0("gunzip ", zip.file))
        }
    }
    ## file <- "PETACC3_expression.tsv"
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    data.set <- fread(file, header=TRUE, fill=TRUE)
    df <- as.data.frame(data.set)
    df <- df[!duplicated(df[,1], fromLast=TRUE) & !duplicated(df[,1], fromLast=FALSE),]
    rownames(df) <- df[,1]
    df <- df[,-1]
    na.flag <- unlist(apply(df, 1, function(row) any(is.na(row) | is.nan(row))))
    data.set <- df[!na.flag,]
    rm(df)
    
    ##    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    symbol <- trans[rownames(data.set),"GeneName"]
    mask <- !is.na(symbol)
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    expr <- combine_probes_2_gene(data.set.m,symbol.m)
    colnames(expr) <- gsub("(.*?)_.*","\\1",colnames(expr))
    expr
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

doaffy.prelim <- function(synId){
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
    file <- getFileLocation(obj)
    print(file)
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
    return(list(ds = data.set.m, sym = symbol.m))
}


## petacc3_expr <- doaffy("syn2175581")

## petacc3_expr <- do.custom.affy("syn2175581")

data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]

## data2npc <- function(x, range) 0

pval.to.text <- function(pval) {
    if(is.na(pval)) { return("n.s.") }
    if(pval > 0.05) {
        return("n.s.")
    } else if(pval <= 0.0001) {
        return("****")
    } else if(pval <= 0.001) {
        return("***")
    } else if(pval <= 0.01) {
        return("**")
    } else if(pval <= 0.05) {
        return("*")
    }
    die("Should not be here\n")
}

## Draw an error bar between the KRAS MT and the KRAS WT points within the
## same scatter plot facet.  There is one facet for each CMS cluster.
draw.mt.wt.err.bar <- function(st, df, tbl1b, ranges, panel.maxs, g, cmsLbl, cmsIndx, adj.pval.col, kras.states) {

    ## Select out the rows corresponding to this facet/CMS.
    mask <- df$cms == cmsLbl

    ## Get the adjusted pvalue (comparing MT vs WT samples) for expression
    ## of st (a gene or gene set) for this CMS cluster.
    ##    adj.pval.col <- paste(toString(cmsLbl), ".apval", sep="")
    if(!(adj.pval.col %in% colnames(tbl1b))) {
        cat(paste0("Could not find ", adj.pval.col, " in tbl1b\n"))
        print(colnames(tbl1b))
        stop("stop")
    }

    pval <- as.numeric(unique(tbl1b[st,adj.pval.col]))
    if(is.na(pval)) {
        cat(paste0("Adj pval in col ", adj.pval.col, " is NA for ", st, "\n"))
        print(tbl1b)
        print(st)
        stop("stop")
    }
    print(pval)

    ## Get the maximum and minimum values in this facet.
    mt.max <- panel.maxs[paste0(cmsLbl,"-",kras.states[1])]
    wt.max <- panel.maxs[paste0(cmsLbl,"-",kras.states[2])]
    m <- max(mt.max, wt.max)

    ## Get the length of the yaxis and define a step size, yoffset, with
    ## respect to it.  We will position the error bars in units of this
    ## step size.
    ymax <- max(ranges[[cmsIndx]][["y.range"]]) - min(ranges[[cmsIndx]][["y.range"]])
    yoffset <- 0.01 * ymax

    ## Use data2npc to translate from data space coordinates to
    ## coordinates used by ggplot within the facet.

    ## This is the vertical line from the top of the KRAS mutant expression
    ## data in this facet to what will be the horizontal line of the error bar
    start <- c(data2npc(1,ranges[[cmsIndx]][["x.range"]]),
               data2npc(mt.max + yoffset, ranges[[cmsIndx]][["y.range"]]))

    end <- c(data2npc(1,ranges[[cmsIndx]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))

    ## The first grob/panel is 4.
    l1 <- 2 + 2 * cmsIndx
    l2 <- 2 + 2 * cmsIndx
    ## This used to work
    t <- 4
    ## This is required in ggplot 2.2.0
    t <- 8
    
    ## Give the grob a random/unique name.  I don't know whether this is
    ## required.
    name=stringi::stri_rand_strings(1,5)
    ## I don't know why this delta is necessary--ggplot2 seems to be
    ## incorrectly confusing/reordering the grobs if I do not make them unique.
    delta <- runif(n=1, min=10^-5, max=10^-3)

    ## Set the current position to the start of the vertical line.
    g <- gtable_add_grob(g, grid.move.to(start[1],start[2],draw=FALSE,name=name), z = Inf, t = t + delta, l = l1, b = 4, r = l1)
    
    ## Draw line from the current position to the end of the vertical line
    ## (and set that end point as the current position)
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l1, r = l1, b = 4)

    ## Similarly, draw the horizontal line of the error bar.
    start <- c(data2npc(1,ranges[[cmsIndx]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    end <- c(data2npc(2,ranges[[cmsIndx]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)

    ## Finally, draw the vertical line on the "right side" of the error bar--
    ## this goes from the error bar to the maximum of the WT KRAS
    ## expression values.
    start <- c(data2npc(2,ranges[[cmsIndx]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    end <- c(data2npc(2,ranges[[cmsIndx]][["x.range"]]),
             data2npc(wt.max + 1 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)
    
    ## Update the maximum values used within the facet to reflect the
    ## newly added error bars
    panel.maxs[paste0(cmsLbl,"-",kras.states[1])] <- m + 4 * yoffset
    panel.maxs[paste0(cmsLbl,"-",kras.states[2])] <- m + 4 * yoffset     

    ## Add the asterisk designation of the pvalue.
    text <- pval.to.text(pval)
    
    ## Position the text in the middle of the error bar.
    pos <- c(data2npc(1.5, ranges[[cmsIndx]][["x.range"]]),
             data2npc(m + 4 * yoffset, ranges[[cmsIndx]][["y.range"]]))

    delta <- runif(n=1, min=10^-5, max=10^-3)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2]), t = t + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
}

## Draw an error bar between the KRAS MT expression values of two different facets
## (each corresponding to a different CMS cluster).
draw.mt.err.bar <- function(st, df, tbl2, ranges, panel.maxs, g, cmsLbl1, cmsLbl2, cmsIndx1, cmsIndx2, cmpLbl, kras.states) {

    ## Select out the rows corresponding to the two facets/CMS clusters to compare.
    mask1 <- df$cms == cmsLbl1
    mask2 <- df$cms == cmsLbl2

    ## Get the adjusted pvalue comparing expression of st (a gene or gene set)
    ## within KRAS MT samples across the two CMS clusters.
    adj.pval.col <- paste(toString(cmpLbl), ".apval", sep="")
    if(!(adj.pval.col %in% colnames(tbl2))) {
        cat(paste0("Could not find ", adj.pval.col, " in tbl2\n"))
    }
    pval <- as.numeric(unique(tbl2[st,adj.pval.col]))
    cat(paste0("Could not find ", adj.pval.col, " in tbl2\n"))
    print(pval)

    ## Get the maximum values across the two facets.
    m <- max(panel.maxs[paste0(cmsLbl1,"-",kras.states[1])], panel.maxs[paste0(cmsLbl2,"-",kras.states[1])], panel.maxs[paste0(cmsLbl1,"-",kras.states[2])], panel.maxs[paste0(cmsLbl2,"-",kras.states[2])])

    ## Get the length of the yaxis and define a step size, yoffset, with
    ## respect to it.  We will position the error bars in units of this
    ## step size.
    ymax <- max(ranges[[cmsIndx1]][["y.range"]]) - min(ranges[[cmsIndx1]][["y.range"]])
    yoffset <- 0.01 * ymax
    
    ## Use data2npc to translate from data space coordinates to
    ## coordinates used by ggplot within the facet.

    ## This is the vertical line from the top of the KRAS mutant expression
    ## data in the first CMS/facet to what will be the horizontal line of the error bar
    start <- c(data2npc(1,ranges[[cmsIndx1]][["x.range"]]),
               data2npc(m + yoffset, ranges[[cmsIndx1]][["y.range"]]))
    
    end <- c(data2npc(1,ranges[[cmsIndx1]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[cmsIndx1]][["y.range"]]))

    ## The first grob/panel is 4.
    l1 <- 2 + 2 * cmsIndx1
    l2 <- 2 + 2 * cmsIndx2
    ## This used to work
    t <- 4
    ## This is required in ggplot 2.2.0
    t <- 8
    
    ## Give the grob a random/unique name.  I don't know whether this is
    ## required.
    name=stringi::stri_rand_strings(1,5)
    ## I don't know why this delta is necessary--ggplot2 seems to be
    ## incorrectly confusing/reordering the grobs if I do not make them unique.
    delta <- runif(n=1, min=10^-5, max=10^-3)

    ## Set the current position to the start of the vertical line.
    g <- gtable_add_grob(g, grid.move.to(start[1],start[2],draw=FALSE,name=name), z = Inf, t = t + delta, l = l1, b = 4, r = l1)

    ## Draw line from the current position to the end of the vertical line
    ## (and set that end point as the current position)
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l1, r = l1, b = 4)

    ## Similarly, draw the horizontal line of the error bar--note that this spans from
    ## one facet (for CMS1) to another (for CMS2)
    start <- c(data2npc(1,ranges[[cmsIndx1]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    
    end <- c(data2npc(1,ranges[[cmsIndx2]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[cmsIndx2]][["y.range"]]))

    ## Finally, draw the vertical line on the "right side" of the error bar--
    ## this goes from the error bar to the maximum of the MT KRAS
    ## expression values in the second facet/CMS 2.
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)
    
    start <- c(data2npc(1,ranges[[cmsIndx2]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    end <- c(data2npc(1,ranges[[cmsIndx2]][["x.range"]]),
             data2npc(m + 1 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)

    ## Update the maximum values used within each facet to reflect the
    ## newly added error bars
    panel.maxs[paste0(cmsLbl1,"-",kras.states[1])] <- m + 4 * yoffset
    panel.maxs[paste0(cmsLbl2,"-",kras.states[1])] <- m + 4 * yoffset     
    panel.maxs[paste0(cmsLbl1,"-",kras.states[2])] <- m + 4 * yoffset
    panel.maxs[paste0(cmsLbl2,"-",kras.states[2])] <- m + 4 * yoffset     

    ## Add the asterisk designation of the pvalue.
    text <- pval.to.text(pval)

    ## Position the text in the middle of the error bar.
    xrange <- ranges[[cmsIndx1]][["x.range"]]
    xrange[2] <- xrange[2] + 3 * (cmsIndx2 - cmsIndx1)
    pos <- c(data2npc(1 + 0.5 * ( 1 + 3 * ( cmsIndx2 - cmsIndx1 ) - 1 ), xrange),
             data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    

    delta <- runif(n=1, min=10^-5, max=10^-3)
    print(text)
    print(pos)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2]), t = t + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
}

to.color <- function(x) {
  c <- "";
  if(x == "CMS1") { c <- "blue" }
  if(x == "CMS2") { c <- "yellow" }
  if(x == "CMS3") { c <- "red" }
  if(x == "CMS4") { c <- "green" }
  if(x == "NOLBL") { c <- "black" }
  c
}

to.number <- function(x) {
  n <- "";
  if(x == "CMS1") { n <- 1 }
  if(x == "CMS2") { n <- 2 }
  if(x == "CMS3") { n <- 3 }
  if(x == "CMS4") { n <- 4 }
  if(x == "NOLBL") { n <- 5 }
  n
}

## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
source("doanalysis.R")


do.harmonized.analysis <- FALSE
read.harmonized.data <- TRUE
harmonized_expr <- NULL
if(read.harmonized.data || do.harmonized.analysis) {
    cat("Reading harmonized data\n")
    obj <- synGet(id="syn2417854", downloadFile = TRUE, downloadLocation = ".")
    harmonized.file <- getFileLocation(obj)
    expr.mat.name <- load(harmonized.file)
    expr.mat <- get(expr.mat.name)

    colnames(expr.mat) <- getSYMBOL(colnames(expr.mat), data='org.Hs.eg')
    ## NB: transposing the matrix
    expr.mat <- t(expr.mat)
    
    # Sanitize the names of the samples so they match the annotation
    cnames <- colnames(expr.mat)
    
    flag <- grepl(pattern="GSM", x=cnames)
    cnames[flag] <- sapply(cnames[flag], function(x) {
      r <- regexpr(pattern="GSM\\d+", text=x)
      if(r==-1) {
        cat(paste("Found no match for GSMd+ in: ", x, "\n", sep=""))
        q(status=-1)
      }
      substr(x, r[1], r[1] + attr(r,"match.length")[1] - 1)
    })
    
    flag <- grepl(pattern=".CEL", x=cnames)
    cnames[flag] <- gsub(cnames[flag], replacement="", pattern=".CEL")
    
    flag <- grepl(pattern="kfsyscc-", x=cnames)
    cnames[flag] <- gsub(cnames[flag], replacement="", pattern="kfsyscc-")
    
    flag <- grepl(pattern="french", x=cnames)
    cnames[flag] <- gsub(cnames[flag], replacement="", pattern="french-")
    
    flag <- grepl(pattern="petacc", x=cnames)
    cnames[flag] <- gsub(cnames[flag], replacement="", pattern="petacc3-")
    
    flag <- grepl(pattern="tcga", x=cnames)
    cnames[flag] <- gsub(cnames[flag], replacement="", pattern="tcgacrc_merged-")
    
    if(any(duplicated(cnames))) {
      cat("Duplicate colnames!\n")
      flag <- duplicated(cnames, fromLast=TRUE) | duplicated(cnames, fromLast=FALSE)
      print(colnames(expr.mat)[flag])
      q(status=-1)
    }
    
    colnames(expr.mat) <- cnames
    harmonized_expr <- expr.mat
    rm(expr.mat)
}

do.tcga.analysis <- TRUE
do.kfs.analysis <- TRUE
do.rasness.analysis.based.on.harmonized.training <- FALSE

tcga_expr <- NULL
kfsyscc_expr <- NULL
tcga_cibersort_mat <- NULL
if(do.tcga.analysis || do.kfs.analysis) {
  cat("Reading TCGA data\n")
  ## TCGA data
  tcga_expr <- { 
    tcga_expr <- read.table(synGet("syn2325328")@filePath,sep="\t",header=TRUE,check.names=FALSE)
    ## tcga_expr <- read.table(opt$`tcga-expr-file`, sep="\t", header=TRUE, check.names=FALSE)
  
    expr <- as.matrix(tcga_expr[,-1])
    rownames(expr) <- tcga_expr[,1]
    expr
  }

  cat("Reading KFS data\n")
  kfsyscc_expr <- doaffy("syn2363564")
  ## kfsyscc_expr <- doaffy(read.table(opt$`kfs-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  colnames(kfsyscc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(kfsyscc_expr))

  tmp <- cbind(gene = rownames(kfsyscc_expr), kfsyscc_expr)
  write.table(file = "kfs-expr-for-cibersort.tsv", tmp, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  tcga_cibersort_mat <- read.table("input/CIBERSORT.Output_Job4_TCGA.tsv", sep="\t", header=TRUE, as.is=TRUE)
  plot.cibersort(tcga_cibersort_mat, clin, "tcga-cibersort")  
  
}


bindeaCellTypes <- unique(bindeaTbl$CellType)
bindeaGsets <- lapply(bindeaCellTypes, function(x){
    return(as.character(unique(bindeaTbl$Symbol[bindeaTbl$CellType==x])))
})
names(bindeaGsets) <- bindeaCellTypes
bindeaGsets <- bindeaGsets[sapply(bindeaGsets, length) > 10]
circ.genes <- read.table("input/circ_sig.txt",as.is=TRUE)[,1]
bindeaGsets[["circ"]] <- circ.genes
circ.genes.in.tcga <- circ.genes[circ.genes %in% rownames(tcga_expr)]
circ.genes.in.kfs <- circ.genes[circ.genes %in% rownames(kfsyscc_expr)]
circ.genes.in.harmonized <- circ.genes[circ.genes %in% rownames(harmonized_expr)]
circ.genes.in.tcga.or.kfs <- unique(c(circ.genes.in.tcga, circ.genes.in.kfs))
bindeaGsets[["harmcirc"]] <- circ.genes.in.harmonized
bindeaGsets[["mek"]] <- c("DUSP6","PHLDA1","SPRY2","DUSP4","ETV4")
## This TAK1 signature was read off of Fig 4A of Singh et al (2012) Cell
tak1.signature <- c("RBP1", "SEMA3A", "SYT1", "EMR2", "PROX1", "INHBB", "ABHD2", "C1orf116", "SNTB1", "TAF9B", "PRF1", "SLC2A1", "GAD1", "MSX2", "PELI2", "ITGB4", "C21orf96", "GPR56", "PDK3", "GLS", "ACSL1", "BIK", "RUNX1", "SYK", "RGL1", "NAV2", "FYN", "HSPA12A", "MBOAT2", "BAMBI", "BMP7", "GGH")
bindeaGsets[["tak1"]] <- tak1.signature

## Wnt target genes from PMID 17320548, as used by J. Guinney in his Nat Med pub.  
## wnt.targets <- c("ASCL2", "AXIN2", "BMP4", "C1orf33", "HIG2", "HSPC111", "KITLG", "LGR5", "MYC", "NOL1", "PPIF", "SOX4", "WRD71", "ZIC2", "ZNRF3")
wnt.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "WNT_FLIER")]]
bindeaGsets[["wnt"]] <- wnt.targets

for(immune.set in immune.sets) {
    bindeaGsets[[immune.set]] <- gene.sets$genesets[[which(gene.sets$geneset.names == immune.set)]]
}

immune.set.gene.symbols <- unique(as.vector(unlist(lapply(bindeaGsets[c(immune.sets,"circ")], function(x) as.vector(x)))))

## to.plot=c(immune.sets, "circ", "mek", "tak1", "myc.sig", "wnt", "BATF3", "ITGAE", "ATF3", "CCL4", "STAT1", "IRF1", "CXCL10", "CIITA", "CD3E", "CD4", "CD8A", "PML", "PIAS1", "MYC")

to.plot=c(immune.sets, "circ")

## Add any CIRC genes that are in TCGA or KFS but not harmonized
circ.genes.missing.in.harmonized <- circ.genes.in.tcga.or.kfs[!(circ.genes.in.tcga.or.kfs %in% circ.genes.in.harmonized)]
cat(paste("CIRC genes in TCGA:\n"))
cat(paste(circ.genes.in.tcga, collapse=" "), "\n")
cat(paste("CIRC genes in KFS:\n"))
cat(paste(circ.genes.in.kfs, collapse=" "), "\n")
cat(paste("CIRC genes in harmonized:\n"))
cat(paste(circ.genes.in.harmonized, collapse=" "), "\n")
cat(paste("CIRC genes in TCGA or KFS, but missing in harmonized: \n"))
cat(paste(circ.genes.missing.in.harmonized, collapse=" "), "\n")
to.plot <- c(to.plot, circ.genes.missing.in.harmonized)

## Also plot the SVGA of only those CIRC genes that are included in harmonized
## to.plot <- c("harmcirc", to.plot)

## to.plot=c("circ", "mek", "tak1", "myc", "wnt")
ylabels <- paste(to.plot, "Expression", sep=" ")
for(i in 1:(length(immune.sets) + 5)) {
    ylabels[i] <- paste0(toupper(to.plot[i]), " Enrichment Score")
}

to.plot <- c("IMMUNE_ESTIMATE", "circ")
num.enrichment <- 2

to.plot=c("circ", "myc.sig", "mek", "wnt", "STAT1", "CXCL10", "CIITA", "IRF1", "CD247")
num.enrichment <- which(to.plot == "wnt")

ylabels <- paste(to.plot, "Expression", sep=" ")
for(i in 1:num.enrichment) {
    ylabels[i] <- paste0(toupper(to.plot[i]), " Enrichment Score")    
}

pdl1 <- which(to.plot == "CD247")
ylabels[pdl1] <- "PD-L1 Expression"

## This is the MSI signature from Tian et al.
## There are 64 genes here--though one is "Unknown."
msi.signature <-  c("ACSL6", "AGR2", "ARID3A", "ASCL2", "ASXL1", "ATP9A", "AXIN", "BC000986", "C10orf47", "C13orf18", "C20orf11", "C20orf43", "CEACAM3", "CEACAM5", "CEP68", "DIDO1", "DUSP18", "DYNLRB1", "EP300", "EPDR1", "FBXO34", "GGA2", "GGT7", "GNG4", "GPR143", "GUCY2C", "HNRNPL", "KCNK5", "KHDRBS3", "KRT23", "LFNG", "LMO4", "LOC157860", "MDM2", "MLH1", "OIT3", "PLAGL2", "PPP1R3D", "PRR15", "QPRT", "RNF43", "ROCK2", "RPL22L1", "SHROOM2", "SHROOM4", "SLC25A22", "SMAD2", "SMCR7L", "SORBS1", "STRN3", "TCF7", "TFCP2L1", "TGFBR2", "TNFSF9", "TNNC2", "TNNT1", "TRIM7", "TSPAN6", "UNKL", "Unknown", "VAV3", "VNN2", "ZFP36L2", "ZSWIM3")


## plot heatmap of msi for (show braf, msi status, gene, # of mutations)
##    tcga
##    kfs
##    french
## correlation k-means things of msi for:
##    tcga
##    kfs
##    french
## svga of msi for:
##    tcga
##    kfs
##    french

stop("wait")

do.amc.analysis <- TRUE
amc_expr <- NULL
if(do.amc.analysis) {
    amc_expr <- doaffy("syn2363559")
    colnames(amc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(amc_expr))

    trans <- read.table("amc-translation-table.tsv", sep="\t", header=TRUE)
    amc_expr <- amc_expr[,colnames(amc_expr) %in% trans[,1]]
    amc_expr <- amc_expr[,trans[,1]]
    colnames(amc_expr) <- trans[,2]
    
    res <- doAnalysis(amc_expr, clin, "amc-all", to.plot=to.plot, ylabels=ylabels)

    ## Restrict to MSS cases
    clin.mss <- clin[!is.na(clin$msi) & (clin$msi == "mss"),]
    res <- doAnalysis(amc_expr, clin.mss, "amc-mss", to.plot=to.plot, ylabels=ylabels)    

    
}

do.agendia.analysis <- TRUE
agendia_expr <- NULL
if(do.agendia.analysis) {
    cat("Doing agendia\n")
    agendia_expr <- do.discover.print.19742("syn2192792")
    colnames(agendia_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(agendia_expr))

    res <- doAnalysis(agendia_expr, clin.all, "agendia-all", to.plot=to.plot, ylabels=ylabels)

    ## Restrict to MSS cases
    clin.mss <- clin.all[!is.na(clin.all$msi) & (clin.all$msi == "mss"),]
    res <- doAnalysis(agendia_expr, clin.mss, "agendia-mss", to.plot=to.plot, ylabels=ylabels)
    cat("Done with agendia\n")
}

do.agendia.ico208.analysis <- TRUE
ico208_expr <- NULL
if(do.agendia.ico208.analysis) {
    cat("Doing ico208\n")
    ico208_expr <- do.discover.print.19742("syn2192796")
    colnames(ico208_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(ico208_expr))

    res <- doAnalysis(ico208_expr, clin.all, "ico208-all", to.plot=to.plot, ylabels=ylabels)

    ## Restrict to MSS cases
    clin.mss <- clin.all[!is.na(clin.all$msi) & (clin.all$msi == "mss"),]
    res <- doAnalysis(ico208_expr, clin.mss, "ico208-mss", to.plot=to.plot, ylabels=ylabels)    
    cat("Done with ico208\n")
}

do.agendia.vh70.analysis <- TRUE
vh70_expr <- NULL
if(do.agendia.vh70.analysis) {
    cat("Doing vh70\n")    
    vh70_expr <- do.discover.print.19742("syn2192799")
    colnames(vh70_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(vh70_expr))

    res <- doAnalysis(vh70_expr, clin.all, "vh70-all", to.plot=to.plot, ylabels=ylabels)

    ## Restrict to MSS cases
    clin.mss <- clin.all[!is.na(clin.all$msi) & (clin.all$msi == "mss"),]
    ## This is crashing: after ~ cms.kras
    ## Error in wilcox.test.default(x = x[flag1], y = x[flag2]) : 
    ## not enough 'y' observations
    res <- doAnalysis(vh70_expr, clin.mss, "vh70-mss", to.plot=to.plot, ylabels=ylabels)    
    cat("Done with vh70\n")
}

do.mdacc.analysis <- TRUE
mdacc_expr <- NULL
if(do.mdacc.analysis) {
    cat("Doing mdacc\n")
    mdacc_expr <- do.discover.print.mdacc("syn2233387")
    colnames(mdacc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(mdacc_expr))

    cat("Doing mdacc clin.all\n")
    res <- doAnalysis(mdacc_expr, clin.all, "mdacc-all", to.plot=to.plot, ylabels=ylabels)

    ## Restrict to MSS cases
    clin.mss <- clin.all[!is.na(clin.all$msi) & (clin.all$msi == "mss"),]
    cat("Doing mdacc clin.mss\n")    
    res <- doAnalysis(mdacc_expr, clin.mss, "mdacc-mss", to.plot=to.plot, ylabels=ylabels)
    cat("Done with mdacc\n")
}

do.petacc.analysis <- TRUE

petacc_expr <- NULL
if(do.petacc.analysis) {
    cat("Reading petacc data\n")
    ## petacc_expr <- doaffy("syn2175581")

    petacc_expr <- do.custom.affy("syn2175581")    
}




## Merge vh70, ico208, and agendia 
common.genes <- intersect(rownames(agendia_expr), rownames(ico208_expr))
common.genes <- intersect(common.genes, rownames(vh70_expr))
common.genes <- intersect(common.genes, rownames(mdacc_expr))
## common.genes <- intersect(common.genes, rownames(petacc_expr))
common.genes <- intersect(common.genes, rownames(amc_expr))

circ.genes[circ.genes %in% common.genes]

merged_expr <-- cbind(agendia_expr[common.genes,], ico208_expr[common.genes,], vh70_expr[common.genes,], mdacc_expr[common.genes,], petacc_expr[common.genes,], amc_expr[common.genes,])


idxs <- match(clin.all$sample, colnames(merged_expr))
clin.mask <- clin.all[!is.na(idxs),]
expr.mask <- merged_expr[, na.omit(idxs)]

flag <- !is.na(clin.mask$msi) & (clin.mask$msi == "mss")

transformed_expr <- subtractOffDatasetEffectWithCombat(expr.mask[,flag], clin.mask[flag,], "subtract-dataset", to.plot=c("circ"), ylabels=c("circ"))

head(merged_expr[circ.genes[circ.genes %in% common.genes],c("6001220","6001222")])
head(transformed_expr[circ.genes[circ.genes %in% common.genes],c("6001220","6001222")])

res <- doAnalysis(transformed_expr, clin.all, "combat-expr", to.plot=c("circ"), ylabels=c("circ"))

include.datasets <- !(clin.all$dataset %in% c("tcga", "kfs", "french"))
res <- doAnalysis(harmonized_expr, clin.all[include.datasets,], "harm-expr", to.plot=c("circ"), ylabels=c("circ"))

include.mss <- !is.na(clin.all$msi) & (clin.all$msi == "mss") 
res <- doAnalysis(harmonized_expr, clin.all[include.datasets & include.mss,], "harm-expr-mss", to.plot=c("circ"), ylabels=c("circ"))

na.flag <- unlist(apply(clin.all[,c("msi", "kras", "cms_label", "site")], 1, function(row) any(is.na(row))))
res <- doAnalysis(harmonized_expr, clin.all[include.datasets & include.mss & !na.flag,], "harm-expr-mss-not-na", to.plot=c("circ"), ylabels=c("circ"))

agendia_es <- computeGSVA(agendia_expr, clin.all, "vh70-all", to.plot=to.plot, ylabels=ylabels)
ico208_es <- computeGSVA(ico208_expr, clin.all, "vh70-all", to.plot=to.plot, ylabels=ylabels)
vh70_es <- computeGSVA(vh70_expr, clin.all, "vh70-all", to.plot=to.plot, ylabels=ylabels)

merged_circ <- cbind(agendia_es["circ",,drop=F], ico208_es["circ",,drop=F], vh70_es["circ",,drop=F])

es.list <- list(agendia = agendia_es, ico208 = ico208_es, vh70 = vh70_es)
for(ds in c("agendia", "ico208", "vh70")) {
    png(paste0("circ-", ds, "-svga.png"))
    d <- dev.off()
}

common.es <- intersect(rownames(agendia_es), rownames(ico208_es))
common.es <- intersect(common.es, rownames(vh70_es))
merged_es <- cbind(agendia_es[common.es,,drop=F], ico208_es[common.es,,drop=F], vh70_es[common.es,,drop=F])

transformed.esn <- subtractOffDatasetEffect(merged_es, clin.all, "subtract-dataset", to.plot=c("circ"), ylabels=c("circ"))

doMergedDatasetAnalysisES(transformed.esn, clin.all, "merged-dataset", to.plot=c("circ"), ylabels=c("circ"), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT"))

doMergedDatasetAnalysisES(merged_es, clin.all, "merged-dataset", to.plot=c("circ"), ylabels=c("circ"), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) 

## Rather than apply combat, add a linear term for the dataset

## Fit the linear model and subtract off the effect of each dataset.
df <- data.frame(expr=merged_expr, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor, dataset=dataset.factor)

        formula <- "KRAS + CMS + site + dataset"
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            formula <- paste0(formula, " + msi")
            df$msi <- msi.factor
        }
        flag <- unlist(apply(df, 1, function(row) any(is.na(row))))
        df <- df[!flag,]

        formula.name <- "full-model-dataset"
        ## Now exclude NOLBL
        df.lm <- df[df$CMS != "NOLBL",]
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")


res <- doMergedDatasetAnalysis(merged_expr, clin.all, "agendia-merged", to.plot=to.plot, ylabels=ylabels)



## Just do linear model that includes dataset


if(FALSE) {
    ## For TCGA, the percent of correlation explained by PC1 is ~45-60%
    res <- doCorrelationAnalysis(tcga_expr, clin, "tcga-all", to.plot=to.plot, ylabels=ylabels)

    ## For KFS, 10-30%
    res <- doCorrelationAnalysis(kfsyscc_expr, clin, "kfs-all", to.plot=to.plot, ylabels=ylabels)    

    ## res <- doCorrelationAnalysis(french_expr, clin, "french-all", to.plot=to.plot, ylabels=ylabels)    
    ## res <- doCorrelationAnalysis(petacc_expr, clin, "petacc-all", to.plot=to.plot, ylabels=ylabels)

    ## For PETACC, 30-40%
    res <- doCorrelationAnalysis(petacc_harm_expr, clin, "petacc-harm-all", to.plot=to.plot, ylabels=ylabels)
}

if(do.tcga.analysis) {
    cat("Doing TCGA analysis\n")
    res <- doAnalysis(tcga_expr, clin, "tcga-all", to.plot=to.plot, ylabels=ylabels)

    tcga.mut.table <- read.table("data_mutations_extended.txt", sep="\t", header=TRUE, as.is=TRUE)
    clin.tcga.msi.inferred <- do.msi(tcga_expr, clin, "tcga-msi", mut.tbl = tcga.mut.table)
    table(clin.tcga.msi.inferred$msi, clin.tcga.msi.inferred$msi.inferred)
    clin.tcga.msi.inferred$msi <- clin.tcga.msi.inferred$msi.inferred
    res <- doAnalysis(tcga_expr, clin.tcga.msi.inferred, "tcga-msi-inferred-all", to.plot=to.plot, ylabels=ylabels)
    
    res <- doCorrelationAnalysis(tcga_expr, clin, "tcga-all", to.plot=to.plot, ylabels=ylabels)    

    ## Inferred MSI
    res <- doAnalysis(tcga_expr, clin.tcga.msi.inferred, "tcga-msi-inferred-all", to.plot=to.plot, ylabels=ylabels)    
    
    ## Restrict to (inferred) MSS cases
    clin.mss.inferred <- clin.tcga.msi.inferred[!is.na(clin.tcga.msi.inferred$msi) & (clin.tcga.msi.inferred$msi == "mss"),]
    res <- doAnalysis(tcga_expr, clin.mss.inferred, "tcga-mss-msi-inferred-all", to.plot=to.plot, ylabels=ylabels)    
    
    ## Restrict to MSS cases
    clin.mss <- clin[!is.na(clin$msi) & (clin$msi == "mss"),]
    res <- doAnalysis(tcga_expr, clin.mss, "tcga-mss", to.plot=to.plot, ylabels=ylabels)    

    ## Restrict to MSS cases and exclude CMS1
    clin.mss.and.no.cms1 <- clin[!is.na(clin$msi) & (clin$msi == "mss") & (is.na(clin$cms_label) | (clin$cms_label != "CMS1")), ]
    res <- doAnalysis(tcga_expr, clin.mss.and.no.cms1, "tcga-mss-no-cms1", to.plot=to.plot, ylabels=ylabels)        
    
##    dset <- "tcga"
##    for(tbl in c("kras.tbl", "cms.tbl", "kras.cms.tbl", "kras.cms.interaction.tbl", "kras.cms.no.interaction.tbl", "kras.mt.vs.wt.tbl")) {
##        tbl.underscore <- gsub(tbl, pattern="\\.", replacement="_")
##        write.table(res[[tbl]], file=paste0("Middleton_", dset, "_", tbl.underscore, ".xls"), sep="\t", quote=FALSE, row.names=FALSE)
##    }
}


if(do.petacc.analysis) {
    cat("Doing petacc analysis\n")
    res <- doAnalysis(petacc_expr, clin, "petacc-all", to.plot=to.plot, ylabels=ylabels)

    ## res <- doCorrelationAnalysis(petacc_expr, clin, "petacc-all", to.plot=to.plot, ylabels=ylabels)    
    
    ## Restrict to MSS cases
    clin.mss <- clin[!is.na(clin$msi) & (clin$msi == "mss"),]
    res <- doAnalysis(petacc_expr, clin.mss, "petacc-mss", to.plot=to.plot, ylabels=ylabels)    

    ## Restrict to MSS cases and exclude CMS1
    clin.mss.and.no.cms1 <- clin[!is.na(clin$msi) & (clin$msi == "mss") & (is.na(clin$cms_label) | (clin$cms_label != "CMS1")), ]
    res <- doAnalysis(petacc_expr, clin.mss.and.no.cms1, "petacc-mss-no-cms1", to.plot=to.plot, ylabels=ylabels)        
    
##    dset <- "tcga"
##    for(tbl in c("kras.tbl", "cms.tbl", "kras.cms.tbl", "kras.cms.interaction.tbl", "kras.cms.no.interaction.tbl", "kras.mt.vs.wt.tbl")) {
##        tbl.underscore <- gsub(tbl, pattern="\\.", replacement="_")
##        write.table(res[[tbl]], file=paste0("Middleton_", dset, "_", tbl.underscore, ".xls"), sep="\t", quote=FALSE, row.names=FALSE)
##    }
}

do.harmonized.petacc.analysis <- TRUE
if(do.harmonized.petacc.analysis) {

    petacc_harm_expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == "petacc"]]
    
    cat("Doing petacc harmonized analysis\n")
    res <- doAnalysis(petacc_harm_expr, clin, "petacc-harm-all", to.plot=to.plot, ylabels=ylabels)

    cat("Doing petacc harmonized analysis\n")
    res <- doCorrelationAnalysis(petacc_harm_expr, clin, "petacc-harm-all", to.plot=to.plot, ylabels=ylabels)
    
    ## Restrict to MSS cases
    clin.mss <- clin[!is.na(clin$msi) & (clin$msi == "mss"),]
    res <- doAnalysis(petacc_harm_expr, clin.mss, "petacc-harm-mss", to.plot=to.plot, ylabels=ylabels)    

    ## Restrict to MSS cases and exclude CMS1
    clin.mss.and.no.cms1 <- clin[!is.na(clin$msi) & (clin$msi == "mss") & (is.na(clin$cms_label) | (clin$cms_label != "CMS1")), ]
    res <- doAnalysis(petacc_harm_expr, clin.mss.and.no.cms1, "petacc-harm-mss-no-cms1", to.plot=to.plot, ylabels=ylabels)        

}

do.french.analysis <- TRUE
french_expr <- NULL
if(do.french.analysis) {
    cat("Reading French data\n")
    french_expr <- doaffy("syn2363561")

    tmp <- cbind(gene = rownames(french_expr), french_expr)
    write.table(file = "french-expr-for-cibersort.tsv", tmp, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    
    res <- doAnalysis(french_expr, clin, "french-all", to.plot=to.plot, ylabels=ylabels)

    clin.french.msi.inferred <- do.msi(french_expr, clin, "french-msi")
    clin.french.msi.inferred$msi <- clin.french.msi.inferred$msi.inferred
    stop("This isn't going to work--need to fix do.msi so that it handles outlier sample; i.e., when msi and mss don't segregate into 2 clusters")
    res <- doAnalysis(french_expr, clin.french.msi.inferred, "french-msi-inferred-all", to.plot=to.plot, ylabels=ylabels)    
    
    res <- doCorrelationAnalysis(french_expr, clin, "french-all", to.plot=to.plot, ylabels=ylabels)    
    
    ## Restrict to MSS cases
    clin.mss <- clin[!is.na(clin$msi) & (clin$msi == "mss"),]
    res <- doAnalysis(french_expr, clin.mss, "french-mss", to.plot=to.plot, ylabels=ylabels)    

    ## Restrict to MSS cases and exclude CMS1
    clin.mss.and.no.cms1 <- clin[!is.na(clin$msi) & (clin$msi == "mss") & (is.na(clin$cms_label) | (clin$cms_label != "CMS1")), ]
    res <- doAnalysis(french_expr, clin.mss.and.no.cms1, "french-mss-no-cms1", to.plot=to.plot, ylabels=ylabels)        

    
    ## french_expr <- doaffy(read.table(opt$`french-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
    
#    cat("Doing French analysis\n")
#    frenchR <- doAnalysis(french_expr, clin, "french-mt", to.plot=to.plot, ylabels=ylabels)
    
#    write.table(frenchR$tbl1, file="Middleton_french_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
                                        #    write.table(frenchR$tbl2, file="Middleton_french_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)

  french_cibersort_mat <- read.table("input/CIBERSORT.Output_Job3_French.tsv", sep="\t", header=TRUE, as.is=TRUE)
  plot.cibersort(french_cibersort_mat, clin, "french-cibersort")  
    
}

if(do.kfs.analysis) {
    cat("Doing KFS analysis\n")
    res <- doAnalysis(kfsyscc_expr, clin, "kfs-all", to.plot=to.plot, ylabels=ylabels)
    dset <- "kfs"

    clin.kfs.msi.inferred <- do.msi(kfsyscc_expr, clin, "kfs-msi")
    clin.kfs.msi.inferred$msi <- clin.kfs.msi.inferred$msi.inferred
    res <- doAnalysis(kfsyscc_expr, clin.kfs.msi.inferred, "kfs-msi-inferred-all", to.plot=to.plot, ylabels=ylabels)

    ## Restrict to (inferred) MSS cases
    clin.mss <- clin.kfs.msi.inferred[!is.na(clin.kfs.msi.inferred$msi) & (clin.kfs.msi.inferred$msi == "mss"),]
    res <- doAnalysis(kfsyscc_expr, clin.mss, "kfs-mss-msi-inferred-all", to.plot=to.plot, ylabels=ylabels)    
    
    res <- doCorrelationAnalysis(kfsyscc_expr, clin, "kfs-all", to.plot=to.plot, ylabels=ylabels)    
    
##    for(tbl in c("kras.tbl", "cms.tbl", "kras.cms.tbl", "kras.cms.interaction.tbl", "kras.cms.no.interaction.tbl", "kras.mt.vs.wt.tbl")) {
##        tbl.underscore <- gsub(tbl, pattern="\\.", replacement="_")
##        write.table(res[[tbl]], file=paste0("Middleton_", dset, "_", tbl.underscore, ".xls"), sep="\t", quote=FALSE, row.names=FALSE)
    ##    }

  kfs_cibersort_mat <- read.table("input/CIBERSORT.Output_Job2_KFS.tsv", sep="\t", header=TRUE, as.is=TRUE)
  plot.cibersort(kfs_cibersort_mat, clin, "kfs-cibersort")  
    
}

do.rasness.analysis.with.training <- FALSE
if(do.rasness.analysis.with.training) {
    cat("Doing rasness analysis\n")
    ## Fit elastic net to french, optimize the cutpoint on petacc, and
    ## test on tcga and kfs

    tcga_expr <- { 
        tcga_expr <- read.table(synGet("syn2325328")@filePath,sep="\t",header=TRUE,check.names=FALSE)
        ## tcga_expr <- read.table(opt$`tcga-expr-file`, sep="\t", header=TRUE, check.names=FALSE)
        
        expr <- as.matrix(tcga_expr[,-1])
        rownames(expr) <- tcga_expr[,1]
        expr
    }
    
    expr <- tcga_expr
    clin.kras <- clin[!is.na(clin$kras),]
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.tcga <- clin.kras[!is.na(idxs),]
    expr.tcga <- expr[, na.omit(idxs)]

    ## Perform TCGA analysis based on RAS mutant
    cat("Doing mutant-based analysis of TCGA\n")            
    tcga.mut.R <- doAnalysis(expr.tcga, clin.tcga, "tcga-mt", to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))
    
    cat("Reading KFS data\n")
    kfsyscc_expr <- doaffy("syn2363564")
    ## kfsyscc_expr <- doaffy(read.table(opt$`kfs-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
    colnames(kfsyscc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(kfsyscc_expr))

    expr <- kfsyscc_expr
    clin.kras <- clin[!is.na(clin$kras),]
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.kfs <- clin.kras[!is.na(idxs),]
    expr.kfs <- expr[, na.omit(idxs)]
    
    ## Perform TCGA analysis based on RAS mutant
    cat("Doing mutant-based analysis of KFS\n")    
    kfs.mut.R <- doAnalysis(expr.kfs, clin.kfs, "kfs-mt", to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))

    break
    
    ## We will train the elastic net on french, optimize the cut threshold on
    ## petacc, and apply to KFS and TCGA.  So read in French and Petacc
    ## and then take the intersection of all four data sets.
    cat("Reading PETACC data\n")
    petacc3_expr <- doaffy("syn2175581")    

    cat("Reading French data\n")
    french_expr <- doaffy("syn2363561")

    ## Since we will train on KFS, but apply to TCGA, restrict to genes that are in TCGA.
##    expr <- expr[kfs.tcga.inters,]
    inters <- rownames(tcga_expr)
    inters <- intersect(inters, rownames(kfsyscc_expr))
##    inters <- intersect(inters, rownames(petacc3_expr))
    inters <- intersect(inters, rownames(french_expr))

    expr.list <- list(kfs=kfsyscc_expr, tcga=tcga_expr, french=french_expr, petacc=petacc3_expr)
    expr.list <- list(kfs=kfsyscc_expr, tcga=tcga_expr, french=french_expr)    
    clin.list <- list()
    for(i in 1:length(expr.list)) {
        name <- names(expr.list)[i]
        expr <- expr.list[[name]]
        expr <- expr[inters, ]
        clin.kras <- clin[!is.na(clin$kras),]
        idxs <- match(clin.kras$sample, colnames(expr))
        clin.list[[name]] <- clin.kras[!is.na(idxs),]
        expr.list[[name]] <- expr[, na.omit(idxs)]
    }

    train.set <- "french"    
    
    alphas <- seq(from = 0, to = 1, by = 0.05)
    nfolds <- 5

    x.train <- t(expr.list[[train.set]])
    x.train <- scale(x.train)

    ## Split the training set in half; the rest will be used to
    ## optimize the cut point
    training.fraction <- 0.5
    set.seed(1234)
    train_ind <- as.vector(createDataPartition(clin.list[[train.set]]$kras, p=training.fraction, list = FALSE))
    x.train.model <- x.train[train_ind, ]
    clin.train.model <- clin.list[[train.set]][train_ind, ]
    x.train.cutpt <- x.train[-train_ind, ]
    clin.train.cutpt <- clin.list[[train.set]][-train_ind, ]

    x.train.model <- x.train
    x.train.cutpt <- x.train.model
    clin.train.model <- clin.list[[train.set]]
    clin.train.cutpt <- clin.train.model
    
    N <- length(clin.train.model$kras)
    foldid <- sample(rep(seq(nfolds), length = N))
    
    cat(paste0("Fitting elastic net to ", train.set, "\n"))
    ##    train.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.train, y = as.factor(clin.list[[train.set]]$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)
    train.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.train.model, y = as.factor(clin.train.model$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)    
    cat("Done training\n")
  
    ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
    pdf(paste0(train.set, "-aucs.pdf"))
    max.num.aucs <- max(unlist(lapply(train.models, function(x) length(x$cvm))))
    aucs <- as.data.frame(lapply(train.models, function(x) { if(length(x$cvm) == max.num.aucs) { x$cvm } else { c(x$cvm, rep(min(x$cvm), max.num.aucs-length(x$cvm))) }} ))
    colnames(aucs) <- 1:ncol(aucs)
    heatmap.2(as.matrix(aucs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
    
    ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
    best.aucs <- unlist(lapply(train.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
    plot(alphas, best.aucs)
    
    best.aucs <- unlist(lapply(train.models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
    plot(alphas, best.aucs)
    d <- dev.off()
    
    ## Best AUCs are achieved for lambdas in the range 0.05 to 0.15.
    ## best.alpha <- alphas[which(best.aucs == max(best.aucs))]
    ## Let's just use alpha = 0.1, which is what Justin published.
    best.alpha <- 0.1
    best.trained.model <- train.models[[which(alphas == best.alpha)]]
    
    ## Genes in our model
    train.genes.min <- rownames(x.train)[which(coef(best.trained.model, s="lambda.min") != 0)]
    train.genes.1se <- rownames(x.train)[which(coef(best.trained.model, s="lambda.1se") != 0)]  

##    cut.pt <- 0.5
    
    ## Optimize the cut point using petacc
##    optimize.set <- "petacc"
##    x.optimize <- t(expr.list[[optimize.set]])
##    x.optimize <- scale(x.optimize)
##    clin.optimize <- clin.list[[optimize.set]]
    
    optimize.set <- "french"
    x.optimize <- x.train.cutpt
    clin.optimize <- clin.train.cutpt
    cat(paste0("Optimizing elastic net cut point using ", optimize.set, "\n"))
    
    optimize.scores <- predict(best.trained.model, newx=x.optimize, s="lambda.1se", type="response")

    mxst <- maxstat(y = clin.optimize$kras, x = optimize.scores, smethod="Wilcoxon", pmethod="condMC")
    print(mxst)
    cut.pt <- as.numeric(mxst$maxstats[[1]]$estimate)
    cat(paste0("Cut point optimized on ", optimize.set, ": ", cut.pt, "\n"))

    ## Predict RIS on all data sets, including those we trained/optimized on
    for(i in 1:length(expr.list)) {
        name <- names(expr.list)[[i]]
        cat(paste0("Predicting RIS on ", name, " (trained on ", train.set, "; optimized on ", optimize.set, ")\n"))

        x.test <- t(expr.list[[name]])
        x.test <- scale(x.test)
        scores <- predict(best.trained.model, newx=x.test, s="lambda.1se", type="response")

        clin.tmp <- clin.list[[name]]
        clin.tmp$rasness <- rep(NA, nrow(clin.tmp))
        idxs <- match(clin.tmp$sample, rownames(scores))
        clin.tmp$rasness[!is.na(idxs)] <- scores[na.omit(idxs)] > cut.pt
        clin.tmp$score.rasness <- rep(NA, nrow(clin.tmp))
        clin.tmp$score.rasness[!is.na(idxs)] <- scores[na.omit(idxs)]
        clin.list[[name]] <- clin.tmp
        
        pdf(paste0(name, "-rasness-vs-mt.pdf"))
        df <- data.frame(score=clin.tmp$score.rasness, KRAS=as.factor(clin.tmp$kras))
        df <- df[!is.na(df$score) & !is.na(df$KRAS),]
        wilcox.test(score ~ KRAS, data = df)  
        p <- ggplot(data=df, aes(x=KRAS, y=score))
        p <- p + ggtitle(paste0(name, " RASness (trained on ", train.set, "; optimized on ", optimize.set, ")"))
        p <- p + ylab("RIS")
        
        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + geom_jitter()
        hline.data <- data.frame(y = cut.pt)
        p <- p + geom_hline(aes(yintercept = y), hline.data)
        print(p)
        d <- dev.off()
    }
    
    ## Perform analysis based on RASness (just for those we didn't
    ## train/optimize on)
    test.sets <- c("kfs", "tcga")
    for(i in 1:length(test.sets)) {
        name <- test.sets[[i]]
        cat(paste0("Doing rasness-based analysis of ", name, "\n"))
        rasness.R <- doAnalysis(expr.list[[name]], clin.list[[name]], paste0(name, "-rasness"), to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness", kras.score.field="score.rasness", kras.states=c("High", "Low"))        
    }
    
    ## Perform analysis based on rasness + RAS mutant (just for those we didn't
    ## train/optimize on)
    for(i in 1:length(test.sets)) {
        name <- test.sets[[i]]
        cat(paste0("Doing rasness+mutant-based analysis of ", name, "\n"))
        clin.tmp <- clin.list[[name]]
        clin.tmp$rasness.and.mut <- rep(NA, nrow(clin.tmp))
        flag <- unlist(apply(clin.tmp[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 1) && (row[2] == 1), FALSE)))
        clin.tmp$rasness.and.mut[flag] <- 1
        flag <- unlist(apply(clin.tmp[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 0) && (row[2] == 0), FALSE)))
        clin.tmp$rasness.and.mut[flag] <- 0
        flag <- !is.na(clin.tmp$braf) & (clin.tmp$braf == 1)
        clin.tmp$rasness.and.mut[flag] <- NA
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        clin.tmp$rasness.and.mut[flag] <- NA
        clin.list[[name]] <- clin.tmp
        rasness.and.mut.R <- doAnalysis(expr.list[[name]], clin.list[[name]], paste0(name, "-rasness-and-mut"), to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.and.mut", kras.states=c("High-MT", "Low-WT"))                
    }
}

do.harmonized.rasness.analysis <- FALSE
if(do.harmonized.rasness.analysis) {

    cat("Doing harmonized rasness analysis\n")
    ## unique(clin$dataset[!is.na(clin$kras)])
    ## [1] "kfs"    "nki"    "french" "petacc" "tcga"   "amc"  

    ## Fit elastic net to the _harmonized_ French data
    train.set <- "french"
    expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == train.set]]
    clin.kras <- clin[!is.na(clin$kras),]
    ## Exclude nras
    flag <- !is.na(clin.kras$nras) & (clin.kras$nras == 1)    
    clin.kras <- clin.kras[!flag, ]
    
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.train.harm <- clin.kras[!is.na(idxs),]
    expr.train.harm <- expr[, na.omit(idxs)]
    
    x.train <- t(expr.train.harm)
    x.train <- scale(x.train)
    
    alphas <- seq(from = 0, to = 1, by = 0.05)
    nfolds <- 5
    N <- length(clin.train.harm$kras)
    foldid <- sample(rep(seq(nfolds), length = N))
    cat(paste0("Fitting elastic net model to ", train.set, "\n"))
    train.harm.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.train, y = as.factor(clin.train.harm$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)
    
    ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
    pdf("train-harm-aucs.pdf")  
    aucs <- as.data.frame(lapply(train.harm.models, function(x) { return(x$cvm) }))
    colnames(aucs) <- 1:ncol(aucs)
    heatmap.2(as.matrix(aucs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
  
    ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
    best.aucs <- unlist(lapply(train.harm.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
    plot(alphas, best.aucs, main="lambda.1se")

    best.aucs <- unlist(lapply(train.harm.models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
    plot(alphas, best.aucs, main="lambda.min")
    d <- dev.off()
    
    ## Best AUCs are achieved for lambdas in the range 0.05 to 0.15.
    ## best.alpha <- alphas[which(best.aucs == max(best.aucs))]
    ## Let's just use alpha = 0.1, which is what Justin published.
    best.alpha <- 0.1
    best.train.harm.trained.model <- train.harm.models[[which(alphas == best.alpha)]]

    ## Optimize the cut point using nki
    optimize.set <- "petacc"
    expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == optimize.set]]
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.optimize.harm <- clin.kras[!is.na(idxs),]
    expr.optimize.harm <- expr[, na.omit(idxs)]
    
    x.optimize <- t(expr.optimize.harm)
    x.optimize <- scale(x.optimize)

    cat(paste0("Optimizing elastic net cutpoint for ", optimize.set, "\n"))
    optimize.harm.scores <- predict(best.train.harm.trained.model, newx=x.optimize, s="lambda.1se", type="response")

    mxst <- maxstat(y = clin.optimize.harm$kras, x = optimize.harm.scores, smethod="Wilcoxon", pmethod="condMC")
    print(mxst)
    cut.pt <- as.numeric(mxst$maxstats[[1]]$estimate)
    
    ## Let's look at the performance of each test set (include the training
    ## set--train--and everything (except train)
    test.sets <- c(unique(clin$dataset[!is.na(clin$kras)]), "all")
    test.sets <- test.sets[!(test.sets %in% c(train.set, optimize.set))]
    ## Also, drop nki, which is small
    ## test.sets <- test.sets[!(test.sets %in% c("nki"))]
    test.sets <- c("french", "petacc", "kfs", "tcga")

    clin$rasness <- rep(NA, nrow(clin))
    clin$score.rasness <- rep(NA, nrow(clin))    
    clin$rasness.and.mut <- rep(NA, nrow(clin))
    
    for(test.set in test.sets) {
        expr <- NULL
        if(test.set == "all") {
            expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[!(clin$dataset %in% c(train.set, optimize.set))]]
        } else {
            expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == test.set]]
        }
        if( ncol(expr) == 0 ) { next }
        idxs <- match(clin.kras$sample, colnames(expr))
        clin.test.harm <- clin.kras[!is.na(idxs),]
        expr.test.harm <- expr[, na.omit(idxs)]
    
        x.test <- t(expr.test.harm)
        x.test <- scale(x.test)

        cat(paste0("Predicting elastic net for ", test.set, "\n"))
        test.harm.scores <- predict(best.train.harm.trained.model, newx=x.test, s="lambda.1se", type="response")

        idxs <- match(clin$sample, rownames(test.harm.scores))
        clin$rasness[!is.na(idxs)] <- test.harm.scores[na.omit(idxs)] > cut.pt
        clin$score.rasness[!is.na(idxs)] <- test.harm.scores[na.omit(idxs)]

        pdf(paste0(test.set, "-harm-", train.set, "-trained-", optimize.set, "-optimized-rasness-vs-mt.pdf"))          
        df <- data.frame(score=clin$score.rasness, KRAS=as.factor(clin$kras))
        df <- df[!is.na(df$score) & !is.na(df$KRAS),]
        wilcox.test(score ~ KRAS, data = df)  
        hline.data <- data.frame(y = cut.pt)
        p <- p + geom_hline(aes(yintercept = y), hline.data)
        print(p)
        d <- dev.off()
    }

    flag <- unlist(apply(clin[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 1) && (row[2] == 1), FALSE)))
    clin$rasness.and.mut[flag] <- 1
    flag <- unlist(apply(clin[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 0) && (row[2] == 0), FALSE)))
    clin$rasness.and.mut[flag] <- 0
    flag <- !is.na(clin$nras) & (clin$nras == 1)
    clin$rasness.and.mut[flag] <- NA
    
}  ## do.harmonized.rasness.analysis

if(do.harmonized.analysis) {

    clin.kras <- clin[!is.na(clin$kras),]
    ## Exclude nras
    flag <- !is.na(clin.kras$nras) & (clin.kras$nras == 1)    
    clin.kras <- clin.kras[!flag, ]
    
    test.sets <- c("french", "petacc", "kfs", "tcga", "all")
    for(test.set in test.sets) {
        expr <- NULL
        if(test.set == "all") {
            expr <- harmonized_expr
        } else {
            expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == test.set]]
        }
        if( ncol(expr) == 0 ) { next }
        idxs <- match(clin.kras$sample, colnames(expr))
        clin.test.harm <- clin.kras[!is.na(idxs),]
        expr.test.harm <- expr[, na.omit(idxs)]

        ## Perform analysis based on RAS mutant
        cat(paste0("Performing harmonized mut-based analysis for ", test.set, "\n"))
        test.harm.mut.R <- doAnalysis(expr.test.harm, clin.test.harm, paste0(test.set, "-harm-mt"), to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))
        
    }
}

if(do.rasness.analysis.based.on.harmonized.training) {
    cat("Doing rasness analysis based on harmonized training\n")

    tcga_expr <- { 
        tcga_expr <- read.table(synGet("syn2325328")@filePath,sep="\t",header=TRUE,check.names=FALSE)
        ## tcga_expr <- read.table(opt$`tcga-expr-file`, sep="\t", header=TRUE, check.names=FALSE)
        
        expr <- as.matrix(tcga_expr[,-1])
        rownames(expr) <- tcga_expr[,1]
        expr
    }
    
    expr <- tcga_expr
    clin.kras <- clin[!is.na(clin$kras),]
    flag <- !is.na(clin.kras$nras) & (clin.kras$nras == 1)    
    clin.kras <- clin.kras[!flag, ]
    
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.tcga <- clin.kras[!is.na(idxs),]
    expr.tcga <- expr[, na.omit(idxs)]

    ## Perform TCGA analysis based on RAS mutant
    cat("Doing mutant-based analysis of TCGA\n")            
    tcga.mut.R <- doAnalysis(expr.tcga, clin.tcga, "tcga-mt", to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))
    
    cat("Reading KFS data\n")
    kfsyscc_expr <- doaffy("syn2363564")
    ## kfsyscc_expr <- doaffy(read.table(opt$`kfs-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
    colnames(kfsyscc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(kfsyscc_expr))

    expr <- kfsyscc_expr
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.kfs <- clin.kras[!is.na(idxs),]
    expr.kfs <- expr[, na.omit(idxs)]
    
    ## Perform TCGA analysis based on RAS mutant
    cat("Doing mutant-based analysis of KFS\n")    
    kfs.mut.R <- doAnalysis(expr.kfs, clin.kfs, "kfs-mt", to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))

    ## Perform analysis based on RASness (just for those we didn't
    ## train/optimize on)
    expr.list <- list(kfs=kfsyscc_expr, tcga=tcga_expr, french=french_expr)    
    test.sets <- c("kfs", "tcga")
    for(i in 1:length(test.sets)) {
        name <- test.sets[[i]]
        cat(paste0("Doing rasness-based analysis of ", name, "\n"))
        rasness.R <- doAnalysis(expr.list[[name]], clin.kras, paste0(name, "-rasness"), to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness", kras.score.field="score.rasness", kras.states=c("High", "Low"))        
    }
    
    ## Perform analysis based on rasness + RAS mutant (just for those we didn't
    ## train/optimize on)
    for(i in 1:length(test.sets)) {
        name <- test.sets[[i]]
        cat(paste0("Doing rasness+mutant-based analysis of ", name, "\n"))
        rasness.and.mut.R <- doAnalysis(expr.list[[name]], clin.kras, paste0(name, "-rasness-and-mut"), to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.and.mut", kras.states=c("High-MT", "Low-WT"))                
    }
}


## Look at the overlap in the immune signatures ...
immune.sigs <- list()
immune.sigs[["circ"]] <- read.table("input/circ_sig.txt",as.is=TRUE)[,1]
for(immune.set in immune.sets) {
    immune.sigs[[immune.set]] <- gene.sets$genesets[[which(gene.sets$geneset.names == immune.set)]]
}

m <- matrix(nrow=length(immune.sigs), ncol=length(immune.sigs), data=0)
for(i in 1:length(immune.sigs)) {
    for(j in 1:length(immune.sigs)) {
        m[i,j] <- length(intersect(immune.sigs[[i]], immune.sigs[[j]]))
    }
}
immune.sig.names <- names(immune.sigs)
immune.sig.names <- gsub(x=immune.sig.names, pattern="_", replacement="\\\n")
rownames(m) <- unlist(immune.sig.names)
colnames(m) <- unlist(immune.sig.names)

library(gridExtra)
library(grid)
pdf("immune-sig-overlap.pdf")
#grid.table(m, base_size = 5)
grid.newpage()
g3 <- tableGrob(m, theme = ttheme_default(base_size = 7, colhead=list(fg_params=list(rot=90))), rows=rownames(m), cols=colnames(m))
grid.draw(g3)
d <- dev.off()

## circ in kfs
## HLA-DQA1 CTLA4 PDCD1LG2 ICAM1 CD274 STAT1 IRF1 IFNG GNLY TBX21 CCL5 LAG3 CD247 ICOS IL18RAP CXCL9 CXCL10 HLA-DPB1 HLA-DPA1 HLA-DMB HLA-DRA HLA-DMA CD80 HLA-DOA CD4 HAVCR2

## circ in tcga
## HLA-DQA1 CTLA4 PDCD1LG2 ICAM1 CD274 STAT1 IRF1 IFNG GNLY TBX21 CCL5 LAG3 CD247 ICOS IL18RAP CXCL9 CXCL10 HLA-DPB1 HLA-DPA1 HLA-DMB HLA-DRA HLA-DMA CD80 HLA-DOA CD4 HAVCR2

## circ in harmonized
## ICAM1 CD274 STAT1 IRF1 GNLY CCL5 CD247 CXCL9 CXCL10 HLA-DPB1 HLA-DPA1 HLA-DMB HLA-DRA HLA-DMA CD80 HLA-DOA CD4 HAVCR2

## TODO
##    1b.  cms2 vs cmsoether 
## 3. cms x kras 
##    3b.  including coarse grained
## 4. cms x kras interaction
##    4b.  including coarse grained

## 5. DE analysis of cms2 mt vs highest
## 6. annotate with kegg--look for immune
## 7. use spia to map to immune
## 8. heatmap of DE genes of cms2 mt vs ... but showing all
##    select out immune, myc, etc.

## Make plot of cms2 mt vs all others
## Differential gene expression between cms2 mt and rest

cms.kras.to.color <- function(x) {
  c <- "";
  if(x == "CMS1-0") { c <- "blue" }
  if(x == "CMS1-1") { c <- "yellow" }
  if(x == "CMS2-0") { c <- "red" }
  if(x == "CMS2-1") { c <- "green" }
  if(x == "CMS3-0") { c <- "brown" }
  if(x == "CMS3-1") { c <- "orange" }
  if(x == "CMS4-0") { c <- "grey" }
  if(x == "CMS4-1") { c <- "yellowgreen" }
  if(x == "NOLBL-0") { c <- "black" }
  if(x == "NOLBL-1") { c <- "gray" }  
  c
}

cms.kras.to.number <- function(x) {
  n <- "";
  if(x == "CMS1-0") { n <- 1 }
  if(x == "CMS1-1") { n <- 2 }
  if(x == "CMS2-0") { n <- 3 }
  if(x == "CMS2-1") { n <- 4 }
  if(x == "CMS3-0") { n <- 5 }
  if(x == "CMS3-1") { n <- 6 }
  if(x == "CMS4-0") { n <- 7 }
  if(x == "CMS4-1") { n <- 8 }
  if(x == "NOLBL-0") { n <- 9 }
  if(x == "NOLBL-1") { n <- 10 }  
  n
}


## is <- immune.sigs[["circ"]][immune.sigs[["circ"]] %in% rownames(expr.m)]
## cms.kras <- apply(cbind(clin.m$cms_label, clin.m$kras), 1, function(row) paste(row[1], row[2], sep="-"))

## d <- order(cms.kras)
## cms.kras.color <- unlist(lapply(cms.kras, cms.kras.to.color))
## heatmap.2(expr.m[is,d],ColSideColors=cms.kras.color[d], scale="row", Colv=FALSE, density.info='none', trace='none', col=colorRampPalette(c('blue','yellow'))(12))

## Just to CMS2 MT vs everything else

## Do limma

## Run through SPIA

stop("stop")

## Check immune genes in (1) CIRC and (2) other
res=spia(de=ds, all=all, organism="hsa", nB=2000, plots=FALSE, beta=NULL)
ds <- rep(-3, length(is))
names(ds) <- is

## translate gene symbols to entrez

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

entrez.to.symbol.mapping <- function(symbols) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), 
              filters = 'entrezgene',
              values = symbols, 
              mart = mart)
  names(bm) <- c("SYM", "INDEX")
  bm <- bm[!(bm$INDEX %in% c("")),]
  bm
}


symbol.to.entrez.mapping <- function(symbols) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), 
              filters = 'hgnc_symbol',
              values = symbols, 
              mart = mart)
  names(bm) <- c("SYM", "INDEX")
  bm <- bm[!(bm$INDEX %in% c("")),]
  bm
}

symbol.to.uniprot.mapping <- function(symbols) {
    mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    bm <- getBM(attributes=c('hgnc_symbol', 'uniprot_swissprot'), 
                filters = 'hgnc_symbol',
                values = symbols, 
                mart = mart)
    names(bm) <- c("SYM", "INDEX")
    bm <- bm[!(bm$INDEX %in% c("")),]
    bm
}

res.ds <- symbol.to.entrez.mapping(names(ds))
res.all <- symbol.to.entrez.mapping(all)

ds.vec <- rep(-1, length(res.ds$INDEX))
names(ds.vec) <- res.ds$INDEX


res=spia(de=ds.vec, all=res.all$INDEX, organism="hsa", nB=2000, plots=FALSE, beta=NULL)


f <- grepl(pattern="immun", res$Name, ignore.case=TRUE)
res[f,]

suppressPackageStartupMessages(library("limma"))

cms.kras.limma <- cms.kras
flag <- cms.kras.limma == "CMS2-1"
cms.kras.limma[flag] <- "CMS2-MT"
cms.kras.limma[!flag] <- "Other"
design <- model.matrix(~ cms.kras.limma)

fit <- lmFit(expr.m, design)
efit <- eBayes(fit)
efit$p.value[is,2]

tt <- topTable(efit, number=Inf, coef="cms.kras.limmaOther", adjust.method="BH")

pc = prcomp( t ( expr.m ) )
plot( pc$x[ , 1:2 ] )

cms.kras <- apply(cbind(clin.m$cms_label, clin.m$kras), 1, function(row) paste(row[1], row[2], sep="-"))
flag <- cms.kras == "CMS2-1"
pval <- apply(expr.m, 1, function(row) wilcox.test(x=row[flag], y=row[!flag])$p.value)
df <- data.frame(gene=rownames(expr.m), pval=pval, padj=p.adjust(pval, method="BH"))

suppressPackageStartupMessages(library("limma"))

df.list <- list()
limma.list <- list()
test.sets <- c("french", "petacc", "kfs", "tcga", "all")
for(test.set in test.sets) {
    expr <- NULL
    if(test.set == "all") {
        expr <- harmonized_expr
    } else {
        expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == test.set]]
    }
    if( ncol(expr) == 0 ) { next }
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.test.harm <- clin.kras[!is.na(idxs),]
    expr.test.harm <- expr[, na.omit(idxs)]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")
    cat(paste0("Doing diff expr of ", test.set, "\n"))
    pval <- apply(expr.test.harm, 1, function(row) wilcox.test(x=row[flagx], y=row[flagy])$p.value)
    df <- data.frame(gene=rownames(expr.test.harm), pval=pval, padj=p.adjust(pval, method="BH"))
    df.list[[test.set]] <- df
    cat(paste0(test.set, ": ", length(which(df$padj < 0.05)), " of ", nrow(df), " genes are dysregulated (wilcox)\n"))

    expr.test.harm <- expr.test.harm[,flagx | flagy]
    clin.test.harm <- clin.test.harm[flagx | flagy,]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")
    
    cms.kras.limma <- rep("foo", nrow(clin.test.harm))
    cms.kras.limma[flagx] <- "CMS2MT"
    cms.kras.limma[flagy] <- "Other"
    cms.kras.limma <- factor(cms.kras.limma, levels=c("Other", "CMS2MT"))
    design <- model.matrix(~ cms.kras.limma)

    fit <- lmFit(expr.test.harm, design)
    efit <- eBayes(fit)
    ## efit$p.value[is,2]

    tt <- topTable(efit, number=Inf, coef="cms.kras.limmaCMS2MT", adjust.method="BH")
    limma.list[[test.set]] <- tt
    cat(paste0(test.set, ": ", length(which(tt$adj.P.Val < 0.05)), " of ", nrow(tt), " genes are dysregulated (limma)\n"))
}

for(test.set in test.sets) {
    df <- df.list[[test.set]]
    cat(paste0(test.set, ": ", length(which(df$padj < 0.05)), " of ", nrow(df), " genes are dysregulated\n"))
}

## 

suppressPackageStartupMessages(library("ToPASeq"))

species <- "hsapiens"
biocarta.pathways <- pathways(species, "biocarta")
kegg.pathways <- pathways(species, "kegg")
reactome.pathways <- pathways(species, "reactome")

## Let's collect all genes in the immune sets and pick out pathways
## that have more than 5 of them.

## biocarta.immune.pathway.flag <- unlist(lapply(biocarta.pathways, function(p) grepl(x=p@title, pattern="immun", ignore.case=TRUE) | grepl(x=p@title, pattern="PD-1", ignore.case=TRUE) | grepl(x=p@title, pattern="CTLA4", ignore.case=TRUE)))
## kegg.immune.pathway.flag <- unlist(lapply(kegg.pathways, function(p) grepl(x=p@title, pattern="immun", ignore.case=TRUE) | grepl(x=p@title, pattern="PD-1", ignore.case=TRUE) | grepl(x=p@title, pattern="CTLA4", ignore.case=TRUE)))
## reactome.immune.pathway.flag <- unlist(lapply(reactome.pathways, function(p) grepl(x=p@title, pattern="immun", ignore.case=TRUE) | grepl(x=p@title, pattern="PD-1", ignore.case=TRUE) | grepl(x=p@title, pattern="CTLA4", ignore.case=TRUE)))

## If there are 5 of more immune genes in a pathway, let's plot it.
min.num.immune.genes <- 5
biocarta.immune.pathway.flag <- unlist(lapply(biocarta.pathways, function(p) ifelse(length(which(nodes(convertIdentifiers(p, "SYMBOL")) %in% immune.set.gene.symbols)) > 5, TRUE, FALSE)))
kegg.immune.pathway.flag <- unlist(lapply(kegg.pathways, function(p) ifelse(length(which(nodes(convertIdentifiers(p, "SYMBOL")) %in% immune.set.gene.symbols)) > 5, TRUE, FALSE)))
reactome.immune.pathway.flag <- unlist(lapply(reactome.pathways, function(p) ifelse(length(which(nodes(convertIdentifiers(p, "SYMBOL")) %in% immune.set.gene.symbols)) > 5, TRUE, FALSE)))

## ls *spia*png | perl -ane '$l = $_; chomp($l); $l =~ s/-spia-all.png//g; $l =~ s/-spia-french.png//g; $l =~ s/-spia-kfs.png//g; $l =~ s/-spia-petacc.png//g; $l =~ s/-spia-tcga.png//g; print STDOUT $l . "\n";' | sort | uniq

## ls *spia-orig*png | perl -ane '$l = $_; chomp($l); $l =~ s/-spia-orig-all.png//g; $l =~ s/-spia-orig-french.png//g; $l =~ s/-spia-orig-kfs.png//g; $l =~ s/-spia-orig-petacc.png//g; $l =~ s/-spia-orig-tcga.png//g; print STDOUT $l . "\n";' | sort | uniq 

## Manually exclude some irrelevant pathways that would be included above",
irrelevant.pathways <- c("AGE-RAGE-signaling-pathway-in-diabetic-complications","Adherens-junction","Allograft-rejection","Chagas-disease-(American-trypanosomiasis)","ER-Phagosome-pathway","Epstein-Barr-virus-infection","Graft-versus-host-disease","Hepatitis-B","Hepatitis-C","Herpes-simplex-infection","Inflammatory-bowel-disease-(IBD)","Influenza-A","Legionellosis","Leishmaniasis","Measles","Pancreatic-cancer","Pertussis","Prostate-cancer","Rheumatoid-arthritis","Salmonella-infection","Semaphorin-interactions","Sphingolipid-signaling-pathway","Staphylococcus-aureus-infection","Systemic-lupus-erythematosus","Toxoplasmosis","Viral-myocarditis")

irrelevant.pathways <- unlist(lapply(irrelevant.pathways, function(x) gsub(x=x, pattern="-", replacement=" ")))

biocarta.immune.pathway.flag <- unlist(lapply(biocarta.pathways, function(p) ifelse((length(which(nodes(convertIdentifiers(p, "SYMBOL")) %in% immune.set.gene.symbols)) > 5) & !(p@title %in% irrelevant.pathways), TRUE, FALSE)))
kegg.immune.pathway.flag <- unlist(lapply(kegg.pathways, function(p) ifelse((length(which(nodes(convertIdentifiers(p, "SYMBOL")) %in% immune.set.gene.symbols)) > 5) & !(p@title %in% irrelevant.pathways), TRUE, FALSE)))
reactome.immune.pathway.flag <- unlist(lapply(reactome.pathways, function(p) ifelse((length(which(nodes(convertIdentifiers(p, "SYMBOL")) %in% immune.set.gene.symbols)) > 5) & !(p@title %in% irrelevant.pathways), TRUE, FALSE)))


names(biocarta.pathways[biocarta.immune.pathway.flag])
names(kegg.pathways[kegg.immune.pathway.flag])
names(reactome.pathways[reactome.immune.pathway.flag])

biocarta.results <- list()
kegg.results <- list()
reactome.results <- list()

cat("Translating gene names to uniprot ids\n")
uniprot.map <- symbol.to.uniprot.mapping(rownames(harmonized_expr))
uniprot.map <- uniprot.map[!is.na(uniprot.map$SYM),]
uniprot.map <- uniprot.map[!is.na(uniprot.map$INDEX),]
uniprot.map <- uniprot.map[!duplicated(uniprot.map$SYM, fromLast=TRUE) & !duplicated(uniprot.map$SYM, fromLast=FALSE),]

entrez.map <- symbol.to.entrez.mapping(rownames(harmonized_expr))
entrez.map <- entrez.map[!is.na(entrez.map$SYM),]
entrez.map <- entrez.map[!is.na(entrez.map$INDEX),]
entrez.map <- entrez.map[!duplicated(entrez.map$SYM, fromLast=TRUE) & !duplicated(entrez.map$SYM, fromLast=FALSE),]

save.image(".Rdata")

source("topaseq-plot.R")

for(test.set in test.sets) {
    tt <- limma.list[[test.set]]

    expr <- NULL
    if(test.set == "all") {
        expr <- harmonized_expr
    } else {
        expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == test.set]]
    }
    if( ncol(expr) == 0 ) { next }
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.test.harm <- clin.kras[!is.na(idxs),]
    expr.test.harm <- expr[, na.omit(idxs)]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")

    expr.test.harm <- expr.test.harm[,flagx | flagy]
    clin.test.harm <- clin.test.harm[flagx | flagy,]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")

    cms.kras.limma <- rep("foo", nrow(clin.test.harm))
    cms.kras.limma[flagx] <- "CMS2-MT"
    cms.kras.limma[flagy] <- "Other"

    group <- unique(cms.kras.limma)
    print(group)

    tt$ID <- rownames(tt)
    tt <- tt[,c("ID", "logFC", "t", "P.Value", "adj.P.Val")]
    names(tt) <- c("ID", "logFC", "t", "pval", "padj")

    tt.entrez.symbol <- merge(tt, entrez.map, by.x = "ID", by.y = "SYM")
    tt.entrez <- tt.entrez.symbol[,c("INDEX", "logFC", "t", "pval", "padj")]
    names(tt.entrez) <- c("ID", "logFC", "t", "pval", "padj")
    cat("Running SPIA agaist biocarta\n")
    biocarta.immune.pathways <- biocarta.pathways[biocarta.immune.pathway.flag]
    biocarta.res <- SPIA(tt.entrez, group=group, type="DEtable", norm.method="none", test.method="ignored", p.th = 0.05, logFC.th = -100, pathways = biocarta.immune.pathways, minEdges = 1, filterSPIA = FALSE)
    biocarta.res$degtable <- tt.entrez
    bad.indices <- unlist(lapply(biocarta.res$topo.sig, function(x) any(is.na(x))))
    biocarta.res$topo.sig <- biocarta.res$topo.sig[!bad.indices]
    biocarta.results[[test.set]] <- biocarta.res
    for(biocarta.immune.pathway in names(biocarta.immune.pathways)) {
##        if(biocarta.immune.pathway != "PD-1 signaling") { next }
        file <- paste0(gsub(biocarta.immune.pathway, pattern=" ", replacement="-"), "-spia-", test.set, ".png", sep="")
        file <- gsub(file, pattern="/", replacement="-")
        cat(paste0("Processing ", file, "\n"))
        if(!is.na(biocarta.res$res$results[biocarta.immune.pathway, "pPERT"]) && !is.na(biocarta.res$res$results[biocarta.immune.pathway, "pNDE"]) && !is.na(biocarta.res$res$results[biocarta.immune.pathway, "pG"])) {
            png(file)
            my.topaseq.plot(biocarta.results[[test.set]], biocarta.immune.pathway, biocarta.immune.pathways, convert=TRUE, graphIDs="SYMBOL", fontsize=40, remNodes = NULL, p.th = 0.001, logical=TRUE)
            d <- dev.off()

            file <- gsub(x=file, pattern=".png", replacement=".tsv")
            tsv <- tt.entrez.symbol[tt.entrez.symbol$INDEX %in% nodes(biocarta.immune.pathways[[biocarta.immune.pathway]]),]
            write.table(file=file, tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        }
    }
    
    cat("Running SPIA agaist kegg\n")
    kegg.immune.pathways <- kegg.pathways[kegg.immune.pathway.flag]
    kegg.res <- SPIA(tt.entrez, group=group, type="DEtable", norm.method="none", test.method="ignored", p.th = 0.05, logFC.th = -100, pathways = kegg.immune.pathways, minEdges = 1, filterSPIA = FALSE)
    kegg.res$degtable <- tt.entrez    
    bad.indices <- unlist(lapply(kegg.res$topo.sig, function(x) any(is.na(x))))
    kegg.res$topo.sig <- kegg.res$topo.sig[!bad.indices]
    kegg.results[[test.set]] <- kegg.res
    for(kegg.immune.pathway in names(kegg.immune.pathways)) {
##        if(kegg.immune.pathway != "PD-1 signaling") { next }        
        file <- paste0(gsub(kegg.immune.pathway, pattern=" ", replacement="-"), "-spia-", test.set, ".png", sep="")
        file <- gsub(file, pattern="/", replacement="-")        
        cat(paste0("Processing ", file, "\n"))
        if(!is.na(kegg.res$res$results[kegg.immune.pathway, "pPERT"]) && !is.na(kegg.res$res$results[kegg.immune.pathway, "pNDE"]) && !is.na(kegg.res$res$results[kegg.immune.pathway, "pG"])) {        
            png(file)
            my.topaseq.plot(kegg.results[[test.set]], kegg.immune.pathway, kegg.immune.pathways, convert=TRUE, graphIDs="SYMBOL", fontsize=40, remNodes = NULL, p.th = 0.001, logical=TRUE)
            d <- dev.off()
            
            file <- gsub(x=file, pattern=".png", replacement=".tsv")
            tsv <- tt.entrez.symbol[tt.entrez.symbol$INDEX %in% nodes(kegg.immune.pathways[[kegg.immune.pathway]]),]
            write.table(file=file, tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        }
    }

    tt.uniprot.symbol <- merge(tt, uniprot.map, by.x = "ID", by.y = "SYM")
    tt.uniprot <- tt.uniprot.symbol[,c("INDEX", "logFC", "t", "pval", "padj")]
    names(tt.uniprot) <- c("ID", "logFC", "t", "pval", "padj")
    
    cat("Running SPIA agaist reactome\n")
    reactome.immune.pathways <- reactome.pathways[reactome.immune.pathway.flag]
    reactome.res <- SPIA(tt.uniprot, group=group, type="DEtable", norm.method="none", test.method="ignored", p.th = 0.05, logFC.th = -100, pathways = reactome.immune.pathways, minEdges = 1, filterSPIA = FALSE)
    reactome.res$degtable <- tt.uniprot    
    bad.indices <- unlist(lapply(reactome.res$topo.sig, function(x) any(is.na(x))))
    reactome.res$topo.sig <- reactome.res$topo.sig[!bad.indices]
    reactome.results[[test.set]] <- reactome.res
    for(reactome.immune.pathway in names(reactome.immune.pathways)) {
##        if(reactome.immune.pathway != "PD-1 signaling") { next }        
        file <- paste0(gsub(reactome.immune.pathway, pattern=" ", replacement="-"), "-spia-", test.set, ".png", sep="")
        file <- gsub(file, pattern="/", replacement="-")
        cat(paste0("Processing ", file, "\n"))
        if(!is.na(reactome.res$res$results[reactome.immune.pathway, "pPERT"]) && !is.na(reactome.res$res$results[reactome.immune.pathway, "pNDE"]) && !is.na(reactome.res$res$results[reactome.immune.pathway, "pG"])) {
            png(file)        
            my.topaseq.plot(reactome.results[[test.set]], reactome.immune.pathway, reactome.immune.pathways, convert=TRUE, IDs="UNIPROT", graphIDs="SYMBOL", fontsize=40, remNodes = NULL, p.th = 0.001, logical=TRUE)
##            my.topaseq.plot(reactome.results[[test.set]], reactome.immune.pathway, reactome.immune.pathways, convert=TRUE, IDs="UNIPROT")
            d <- dev.off()

            file <- gsub(x=file, pattern=".png", replacement=".tsv")
            tsv <- tt.uniprot.symbol[tt.uniprot.symbol$INDEX %in% nodes(reactome.immune.pathways[[reactome.immune.pathway]]),]
            write.table(file=file, tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
            
        }
    }
}

save.image(".Rdata")

## Repeat the above analysis, but use the original tcga and kfs data sets, not the harmonized one.
## NB: harmonized does not have PD1 or CTLA4
df.orig.list <- list()
limma.orig.list <- list()
test.sets <- c("kfs", "tcga")
expr.list <- list(kfsyscc_expr, tcga_expr)
for(idx in 1:length(test.sets)) {
    test.set <- test.sets[idx]
    expr <- expr.list[[idx]]
    if( ncol(expr) == 0 ) { next }
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.test.harm <- clin.kras[!is.na(idxs),]
    expr.test.harm <- expr[, na.omit(idxs)]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")
    cat(paste0("Doing diff expr of ", test.set, "\n"))
    pval <- apply(expr.test.harm, 1, function(row) wilcox.test(x=row[flagx], y=row[flagy])$p.value)
    df <- data.frame(gene=rownames(expr.test.harm), pval=pval, padj=p.adjust(pval, method="BH"))
    df.orig.list[[test.set]] <- df
    cat(paste0(test.set, ": ", length(which(df$padj < 0.05)), " of ", nrow(df), " genes are dysregulated (wilcox)\n"))

    expr.test.harm <- expr.test.harm[,flagx | flagy]
    clin.test.harm <- clin.test.harm[flagx | flagy,]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")
    
    cms.kras.limma <- rep("foo", nrow(clin.test.harm))
    cms.kras.limma[flagx] <- "CMS2MT"
    cms.kras.limma[flagy] <- "Other"
    cms.kras.limma <- factor(cms.kras.limma, levels=c("Other", "CMS2MT"))
    design <- model.matrix(~ cms.kras.limma)

    fit <- lmFit(expr.test.harm, design)
    efit <- eBayes(fit)
    ## efit$p.value[is,2]

    tt <- topTable(efit, number=Inf, coef="cms.kras.limmaCMS2MT", adjust.method="BH")
    limma.orig.list[[test.set]] <- tt
    cat(paste0(test.set, ": ", length(which(tt$adj.P.Val < 0.05)), " of ", nrow(tt), " genes are dysregulated (limma)\n"))
}


source("topaseq-plot.R")

uniprot.orig.map <- symbol.to.uniprot.mapping(unique(c(rownames(kfsyscc_expr), rownames(tcga_expr))))
uniprot.orig.map <- uniprot.orig.map[!is.na(uniprot.orig.map$SYM),]
uniprot.orig.map <- uniprot.orig.map[!is.na(uniprot.orig.map$INDEX),]
uniprot.orig.map <- uniprot.orig.map[!duplicated(uniprot.orig.map$SYM, fromLast=TRUE) & !duplicated(uniprot.orig.map$SYM, fromLast=FALSE),]

entrez.orig.map <- symbol.to.entrez.mapping(unique(c(rownames(kfsyscc_expr), rownames(tcga_expr))))
entrez.orig.map <- entrez.orig.map[!is.na(entrez.orig.map$SYM),]
entrez.orig.map <- entrez.orig.map[!is.na(entrez.orig.map$INDEX),]
entrez.orig.map <- entrez.orig.map[!duplicated(entrez.orig.map$SYM, fromLast=TRUE) & !duplicated(entrez.orig.map$SYM, fromLast=FALSE),]

biocarta.orig.results <- list()
kegg.orig.results <- list()
reactome.orig.results <- list()


for(idx in 1:length(test.sets)) {
    test.set <- test.sets[idx]
    expr <- expr.list[[idx]]
    tt <- limma.orig.list[[test.set]]

    if( ncol(expr) == 0 ) { next }
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.test.harm <- clin.kras[!is.na(idxs),]
    expr.test.harm <- expr[, na.omit(idxs)]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")

    expr.test.harm <- expr.test.harm[,flagx | flagy]
    clin.test.harm <- clin.test.harm[flagx | flagy,]
    flagx <- (clin.test.harm$kras == 1) & (clin.test.harm$cms_label == "CMS2")
    flagy <- !flagx
    flagy <- (clin.test.harm$kras == 0) & (clin.test.harm$cms_label != "CMS2")

    cms.kras.limma <- rep("foo", nrow(clin.test.harm))
    cms.kras.limma[flagx] <- "CMS2-MT"
    cms.kras.limma[flagy] <- "Other"

    group <- unique(cms.kras.limma)
    print(group)

    tt$ID <- rownames(tt)
    tt <- tt[,c("ID", "logFC", "t", "P.Value", "adj.P.Val")]
    names(tt) <- c("ID", "logFC", "t", "pval", "padj")

    tt.entrez.symbol <- merge(tt, entrez.orig.map, by.x = "ID", by.y = "SYM")
    tt.entrez <- tt.entrez.symbol[,c("INDEX", "logFC", "t", "pval", "padj")]
    names(tt.entrez) <- c("ID", "logFC", "t", "pval", "padj")
    tt.entrez <- tt.entrez[!duplicated(tt.entrez$ID, fromLast=TRUE) & !duplicated(tt.entrez$ID, fromLast=FALSE),]
    cat("Running SPIA agaist biocarta\n")
    biocarta.immune.pathways <- biocarta.pathways[biocarta.immune.pathway.flag]
    biocarta.res <- SPIA(tt.entrez, group=group, type="DEtable", norm.method="none", test.method="ignored", p.th = 0.05, logFC.th = -100, pathways = biocarta.immune.pathways, minEdges = 1, filterSPIA = FALSE)
    biocarta.res$degtable <- tt.entrez
    bad.indices <- unlist(lapply(biocarta.res$topo.sig, function(x) any(is.na(x))))
    biocarta.res$topo.sig <- biocarta.res$topo.sig[!bad.indices]
    biocarta.orig.results[[test.set]] <- biocarta.res
    for(biocarta.immune.pathway in names(biocarta.immune.pathways)) {
##        if(biocarta.immune.pathway != "PD-1 signaling") { next }
        file <- paste0(gsub(biocarta.immune.pathway, pattern=" ", replacement="-"), "-spia-orig-", test.set, ".png", sep="")
        file <- gsub(file, pattern="/", replacement="-")
        cat(paste0("Processing ", file, "\n"))
        if(!is.na(biocarta.res$res$results[biocarta.immune.pathway, "pPERT"]) && !is.na(biocarta.res$res$results[biocarta.immune.pathway, "pNDE"]) && !is.na(biocarta.res$res$results[biocarta.immune.pathway, "pG"])) {
            png(file)
            my.topaseq.plot(biocarta.orig.results[[test.set]], biocarta.immune.pathway, biocarta.immune.pathways, convert=TRUE, graphIDs="SYMBOL", fontsize=40, remNodes = NULL, p.th = 0.001, logical=TRUE)
            d <- dev.off()

            file <- gsub(x=file, pattern=".png", replacement=".tsv")
            tsv <- tt.entrez.symbol[tt.entrez.symbol$INDEX %in% nodes(biocarta.immune.pathways[[biocarta.immune.pathway]]),]
            write.table(file=file, tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        }
    }
    
    cat("Running SPIA agaist kegg\n")
    kegg.immune.pathways <- kegg.pathways[kegg.immune.pathway.flag]
    kegg.res <- SPIA(tt.entrez, group=group, type="DEtable", norm.method="none", test.method="ignored", p.th = 0.05, logFC.th = -100, pathways = kegg.immune.pathways, minEdges = 1, filterSPIA = FALSE)
    kegg.res$degtable <- tt.entrez    
    bad.indices <- unlist(lapply(kegg.res$topo.sig, function(x) any(is.na(x))))
    kegg.res$topo.sig <- kegg.res$topo.sig[!bad.indices]
    kegg.orig.results[[test.set]] <- kegg.res
    for(kegg.immune.pathway in names(kegg.immune.pathways)) {
##        if(kegg.immune.pathway != "PD-1 signaling") { next }        
        file <- paste0(gsub(kegg.immune.pathway, pattern=" ", replacement="-"), "-spia-orig-", test.set, ".png", sep="")
        file <- gsub(file, pattern="/", replacement="-")        
        cat(paste0("Processing ", file, "\n"))
        if(!is.na(kegg.res$res$results[kegg.immune.pathway, "pPERT"]) && !is.na(kegg.res$res$results[kegg.immune.pathway, "pNDE"]) && !is.na(kegg.res$res$results[kegg.immune.pathway, "pG"])) {        
            png(file)
            my.topaseq.plot(kegg.orig.results[[test.set]], kegg.immune.pathway, kegg.immune.pathways, convert=TRUE, graphIDs="SYMBOL", fontsize=40, remNodes = NULL, p.th = 0.001, logical=TRUE)
            d <- dev.off()
            
            file <- gsub(x=file, pattern=".png", replacement=".tsv")
            tsv <- tt.entrez.symbol[tt.entrez.symbol$INDEX %in% nodes(kegg.immune.pathways[[kegg.immune.pathway]]),]
            write.table(file=file, tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        }
    }

    tt.uniprot.symbol <- merge(tt, uniprot.orig.map, by.x = "ID", by.y = "SYM")
    tt.uniprot <- tt.uniprot.symbol[,c("INDEX", "logFC", "t", "pval", "padj")]
    names(tt.uniprot) <- c("ID", "logFC", "t", "pval", "padj")
    tt.uniprot <- tt.uniprot[!duplicated(tt.uniprot$ID, fromLast=TRUE) & !duplicated(tt.uniprot$ID, fromLast=FALSE),]
    
    cat("Running SPIA agaist reactome\n")
    reactome.immune.pathways <- reactome.pathways[reactome.immune.pathway.flag]
    reactome.res <- SPIA(tt.uniprot, group=group, type="DEtable", norm.method="none", test.method="ignored", p.th = 0.05, logFC.th = -100, pathways = reactome.immune.pathways, minEdges = 1, filterSPIA = FALSE)
    reactome.res$degtable <- tt.uniprot    
    bad.indices <- unlist(lapply(reactome.res$topo.sig, function(x) any(is.na(x))))
    reactome.res$topo.sig <- reactome.res$topo.sig[!bad.indices]
    reactome.orig.results[[test.set]] <- reactome.res
    for(reactome.immune.pathway in names(reactome.immune.pathways)) {
##        if(reactome.immune.pathway != "PD-1 signaling") { next }        
        file <- paste0(gsub(reactome.immune.pathway, pattern=" ", replacement="-"), "-spia-orig-", test.set, ".png", sep="")
        file <- gsub(file, pattern="/", replacement="-")
        cat(paste0("Processing ", file, "\n"))
        if(!is.na(reactome.res$res$results[reactome.immune.pathway, "pPERT"]) && !is.na(reactome.res$res$results[reactome.immune.pathway, "pNDE"]) && !is.na(reactome.res$res$results[reactome.immune.pathway, "pG"])) {
            png(file)        
            my.topaseq.plot(reactome.orig.results[[test.set]], reactome.immune.pathway, reactome.immune.pathways, convert=TRUE, IDs="UNIPROT", graphIDs="SYMBOL", fontsize=40, remNodes = NULL, p.th = 0.001, logical=TRUE)
##            my.topaseq.plot(reactome.results[[test.set]], reactome.immune.pathway, reactome.immune.pathways, convert=TRUE, IDs="UNIPROT")
            d <- dev.off()

            file <- gsub(x=file, pattern=".png", replacement=".tsv")
            tsv <- tt.uniprot.symbol[tt.uniprot.symbol$INDEX %in% nodes(reactome.immune.pathways[[reactome.immune.pathway]]),]
            write.table(file=file, tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
            
        }
    }
}

save.image(".Rdata")

stop("stop")

my.topaseq.plot(reactome.results[["french"]], "PD-1 signaling", reactome.pathways, convert=TRUE)
path <- "PD-1 signaling"
my.topaseq.plot(reactome.results[["french"]], path, reactome.pathways, convert=TRUE, IDs="UNIPROT")
my.topaseq.plot(reactome.results[["french"]], path, reactome.pathways, convert=TRUE, IDs="UNIPROT", graphIDs="UNIPROT")
## plot(reactome.results[["french"]], "PD-1 signaling", reactome.pathways, convert=TRUE)

for(test.set in test.sets) {    
    if(test.set == "all") {
        expr <- harmonized_expr
    } else {
        expr <- harmonized_expr[, colnames(harmonized_expr) %in% clin$sample[clin$dataset == test.set]]
    }
    if( ncol(expr) == 0 ) { next }
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.test.harm <- clin.kras[!is.na(idxs),]
    expr.test.harm <- expr[, na.omit(idxs)]
    pc = prcomp( t ( expr.test.harm ) )
    labels <- clin.test.harm$cms_label
    color.labels <- unlist(lapply(labels, to.color))
    flag <- labels == "CMS2"
##    labels[!flag] <- ""
    pdf(paste0(test.set, "-pca.pdf"))
    u <- !duplicated(labels)
    plot( pc$x[ , 1:2 ], pch = 20, col = color.labels )
    legend("topright", legend = labels[u], fill = color.labels[u])
##    text(x = pc$x[,1], y = pc$x[,2], labels=labels)
    d <- dev.off()
}


  res <- NULL

  # Run DESeq2
  if(run.deseq) {
    cnt.mat <- tsv.entrez[,cnt.cols]
    rownames(cnt.mat) <- tsv.entrez$ID

    res <- SPIA(cnt.mat, group=group, type="RNASeq", norm.method="DESeq2", test.method="DESeq2", p.th = p.th, logFC.th = logFC.th, pathways=p, minEdges = 1, filterSPIA = FALSE)
  } else {
    rownames(tsv.entrez) <- tsv.entrez$ID  
    deg.table <- tsv.entrez[,c("ID", "logFC", "t", "pval", "padj")]
    res <- SPIA(tsv.entrez, group=group, type="DEtable", norm.method="DESeq2", test.method="DESeq2", p.th = p.th, logFC.th = logFC.th, pathways=p, minEdges = 1, filterSPIA = FALSE, deg.table.arg = deg.table)
  }


## TODO:
## Why is PD1 failing -- fixed?

## Fix plot
## Make sure we have the right entries
## Immune genes

## TODO:
## -- pathways with CD28/CTLA-4 (ligands: CD80 and CD86) and PD1
## -- immune response pathways
## -- are these missing from harmonized?
## -- what is gene atop the PD1 signaling cascade?
## -- role of toll-like receptors?
## -- TCR
## -- the-co-stimulatory-signal-during-t-cell-activation-spia-french.tsv ****
## -- t-cell-receptor-signaling-pathway 
## -- Costimulation-by-the-CD28-family
## -- why does CTLA4 not show up in the-co-stimulatory-signal-during-t-cell-activation-spia-french.tsv 
## -- PD1 does not show up anywhere

## -- add these pathways to the immune signatures that we look at.
## -- what does IMMUNE ESTIMATE do?
## -- is PD-L1 or PD-L2 expressed?
## -- is CD86 or CD80 expressed?
## -- is AKT or STAT3 upregulated in 'other,' which might lead to pd-1 upregulation?


