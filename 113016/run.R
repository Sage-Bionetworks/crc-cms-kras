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


## Immune sets that Justin looked at.
immune.sets <- c("IMMUNE_ESTIMATE", "IMMUNE_RESP_GO_BP", "PD1_REACTOME", "IMMUNE_CD8MACRO_GALON", "IMMUNE_TH1_GALON", "IMMUNE_NKC_BREAST", "IMMUNE_THF_BREAST", "IMMUNE_TH17_GOUNARI", "IMMUNE_TREG_LUCAS", "IMMUNE_MDSC_ALBELDA")

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

synapseLogin("brian.white")

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

## Read in the gene-set definitions used in Justin's Nat Med paper
obj <- synGet(id="syn2321865", downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
gene.sets <- GSA.read.gmt(file)

## Match the clinical annotations and the CMS clusters
idxs <- match(clin$sample, cms$sample)
clin <- clin[!is.na(idxs),]
clin$cms_label <- cms$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples[na.omit(idxs)]

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

## petacc3_expr <- doaffy("syn2175581")

data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]

## data2npc <- function(x, range) 0

pval.to.text <- function(pval) {
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
    pval <- as.numeric(unique(tbl1b[st,adj.pval.col]))

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
    pval <- as.numeric(unique(tbl2[st,adj.pval.col]))

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

## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
doAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.braf.and.nras.mutants <- FALSE
    if(exclude.braf.and.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
    }

    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin.tmp$sample, colnames(expr))
    clin.mask <- clin.tmp[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]
    
    ## Exclude any samples that do not have KRAS mutation status annotated
    mask <- !is.na(clin.mask[,kras.status.field])
    clin.m <- clin.mask[mask,]
    expr.m <- expr.mask[, mask]
    
    ## gene set prep
    env <- new.env()
    load("input/markersG.RData",envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
    bindeaTbl <- read.table("input/bindea_immune.txt",sep='\t',header=TRUE)
    
    ## Restrict to a few cell types of interest
    bindeaTbl <- bindeaTbl[bindeaTbl$CellType %in% c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells"),]
    
    bindeaCellTypes <- unique(bindeaTbl$CellType)
    bindeaGsets <- lapply(bindeaCellTypes, function(x){
        return(as.character(unique(bindeaTbl$Symbol[bindeaTbl$CellType==x])))
    })
    names(bindeaGsets) <- bindeaCellTypes
    bindeaGsets <- bindeaGsets[sapply(bindeaGsets, length) > 10]
    bindeaGsets[["circ"]] <- read.table("input/circ_sig.txt",as.is=TRUE)[,1]
    bindeaGsets[["mek"]] <- c("DUSP6","PHLDA1","SPRY2","DUSP4","ETV4")
    # This TAK1 signature was read off of Fig 4A of Singh et al (2012) Cell
    tak1.signature <- c("RBP1", "SEMA3A", "SYT1", "EMR2", "PROX1", "INHBB", "ABHD2", "C1orf116", "SNTB1", "TAF9B", "PRF1", "SLC2A1", "GAD1", "MSX2", "PELI2", "ITGB4", "C21orf96", "GPR56", "PDK3", "GLS", "ACSL1", "BIK", "RUNX1", "SYK", "RGL1", "NAV2", "FYN", "HSPA12A", "MBOAT2", "BAMBI", "BMP7", "GGH")
    bindeaGsets[["tak1"]] <- tak1.signature
    
    # Wnt target genes from PMID 17320548, as used by J. Guinney in his Nat Med pub.  
    # wnt.targets <- c("ASCL2", "AXIN2", "BMP4", "C1orf33", "HIG2", "HSPC111", "KITLG", "LGR5", "MYC", "NOL1", "PPIF", "SOX4", "WRD71", "ZIC2", "ZNRF3")
    wnt.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "WNT_FLIER")]]
    bindeaGsets[["wnt"]] <- wnt.targets

    for(immune.set in immune.sets) {
        bindeaGsets[[immune.set]] <- gene.sets$genesets[[which(gene.sets$geneset.names == immune.set)]]
    }

    # Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    # myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
    myc.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "MYC_TARGETS_ZELLER")]]
    bindeaGsets[["myc.sig"]] <- myc.targets
    
    cat(paste("Number of rows analyzed: ", nrow(expr.m), "\n", sep=""))
    cat(paste("Number of columns analyzed: ", ncol(expr.m), "\n", sep=""))

    ## save(expr.m, file="expr.m.Rdata")
    ## save(bindeaGsets, file="bind.Rdata")
    
    ## Output genes in each set that are in data set
    for(i in 1:length(bindeaGsets)) {
      genes <- bindeaGsets[[i]][ bindeaGsets[[i]] %in% rownames(expr.m) ]
      if(length(genes) == 0) { genes <- "None" }
      genes <- paste(genes, collapse=" ")
      cat(paste("Genes from data set ", analysis.name, " analyzed in set ", names(bindeaGsets)[i], ": ", genes, "\n", sep=""))
    }
    
    ## Perform GSVA analysis
    es <- gsva(expr.m, bindeaGsets,parallel.sz=num.processes,verbose=TRUE)$es.obs
    cat("Done with gsva.\n")
    
    ## Make sure that the genes/gene sets to plot are actually in the expression set.
    all.genes.and.sets <- c(rownames(expr.m), rownames(es))
    flag <- to.plot %in% all.genes.and.sets
    to.plot <- to.plot[flag]
    ylabels <- ylabels[flag]
    tmp <- rbind(es, expr.m[to.plot[to.plot %in% rownames(expr.m)],,drop=F])

    ## Exclude any samples that do not have KRAS mutation status annotated
    if(FALSE) {
        mask <- !is.na(clin.m[,kras.status.field])
        kras <- clin.m[mask, kras.status.field]
        kras.score <- clin.m[mask, kras.score.field]
        cms <- as.character(clin.m$cms_label[mask])
        esn <- tmp[, mask]
        expr.mask <- expr.m[,mask]
    } else {
        kras <- clin.m[, kras.status.field]
        kras.score <- clin.m[, kras.score.field]        
        cms <- as.character(clin.m$cms_label)
        msi <- clin.m$msi
        site <- clin.m$site
        esn <- tmp
        expr.mask <- expr.m
    }
    
    if(FALSE) {
      save(clin.m, file="clin.Rdata")
      save(kras, file="kras.Rdata")
      save(esn, file="esn.Rdata")
      save(bindeaGsets, file="bind.Rdata")
      save(expr.mask, file="expr.Rdata")
      save(cms, file="cms.Rdata")
    }

    ## Analysis 5:  look at correlation of
    ## a. wnt vs circ
    ## b. tak vs circ
    ## c. wnt vs myc
    sigs1 <- c("wnt", "tak1", "wnt")
    sigs2 <- c("circ", "circ", "myc.sig")
    for(i in 1:length(sigs1)) {
        sig1 <- sigs1[i]
        sig2 <- sigs2[i]
        out.pdf <- paste("scatter-", analysis.name, "-", sig1, "-vs-", sig2, ".pdf", sep="")
        pdf(out.pdf)
        ## par(mfrow=c(2,5))

        par(mfrow=c(2,1))
        cms.lbl <- "CMS1"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {        
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }

        cms.lbl <- "CMS2"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        cms.lbl <- "CMS3"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        cms.lbl <- "CMS4"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        if(any(flag)) {
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        cms.lbl <- "ALL"
        kras.status <- 1
        kras.indx <- 1
        flag <- kras == kras.status
        if(any(flag)) {        
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        
        kras.status <- 0
        kras.indx <- 2
        flag <- kras == kras.status
        if(any(flag)) {        
            plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
            lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
            lm.sum <- summary(lm.obj)
            coeffs <- coefficients(lm.sum)
            if((nrow(coeffs) == 2) && !any(is.na(coeffs[,1]))) {
                abline(lm.obj)
            }
        }
        d <- dev.off()
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]

    ## Analysis 3: linear modeling of biomarker ~ CMS + KRAS + CMS:KRAS
    tbl <- data.frame(variable=to.plot)
    rownames(tbl) <- to.plot
    cols <- c(cms2mt.vs.wt.col, cms1mt.vs.cms2mt.col, cms3mt.vs.cms2mt.col, cms4mt.vs.cms2mt.col, nolblmt.vs.cms2mt.col, kras.col)
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        tbl[,col] <- rep(NA, nrow(tbl))
    } 

    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        if(FALSE) {
        file <- paste0(st, "-", analysis.name, "-scatter.pdf")
        pdf(file)
        plot(kras.score, esn.st)
        d <- dev.off()

        df.tmp <- data.frame(esn=esn.st, cms=factor(cms), KRAS=(kras.score))
        for(cms.label in unique(df.tmp$cms)) {
            df.cms <- df.tmp[df.tmp$cms == cms.label,]
            file <- paste0(st, "-", analysis.name, "-", cms.label, "-scatter.pdf")
            pdf(file)
            plot(df.cms$KRAS, df.cms$esn, xlab="RASness", ylab=st, main=cms.label)
            
            d <- dev.off()            
        }
        }
        
        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        ## Create a data frame where each row is a sample and the columns are:
        ## (1) expression of the gene/gene set
        ## (2) CMS factor of the sample
        ## (3) KRAS mutation status of the sample
        ## Make CMS2 the first factor
        unique.cms <- unique(cms)
        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=c("CMS2", unique.cms[unique.cms != "CMS2"])), KRAS=factor(kras.label))

        formulas <- c("KRAS", "CMS", "KRAS + CMS", "KRAS + CMS + CMS:KRAS")
        formula.names <- c("kras", "cms", "kras-cms", "kras-cms-interaction")
        for(i in 1:length(formulas)) {
            df.lm <- df
            if(formula.names[i] != "kras") {
                ## Now exclude NOLBL
                ## df.lm <- df[df$CMS != "NOLBL",]                
            }
            lm.obj <- lm(as.formula(paste("expr", formulas[i], sep=" ~ ")), data=df.lm)
            lm.sum <- summary(lm.obj)
            coeffs <- as.data.frame(coefficients(lm.sum))
            lm.file <- paste(st, "-", analysis.name, "-", formula.names[i], "-lm.tsv", sep="")
            coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
            write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
            if(formula.names[i] == "kras") {
                kras.flag <- grepl(x=coeffs.2$coefficient, pattern="KRAS") & !grepl(x=coeffs.2$coefficient, pattern=":")
                col <- paste0(kras.col, ".pval")
                tbl[sti, col] <- coeffs.2[kras.flag, 5]
            }
            if(formula.names[i] == "kras-cms-interaction") {
                flag <- grepl(x=coeffs.2$coefficient, pattern="KRAS") & !grepl(x=coeffs.2$coefficient, pattern=":")
                col <- paste0(cms2mt.vs.wt.col, ".pval")
                tbl[sti, col] <- coeffs.2[flag, 5]
                flag <- grepl(x=coeffs.2$coefficient, pattern="CMS4") & !grepl(x=coeffs.2$coefficient, pattern=":")
                col <- paste0(cms4mt.vs.cms2mt.col, ".pval")
                tbl[sti, col] <- coeffs.2[flag, 5]
                flag <- grepl(x=coeffs.2$coefficient, pattern="CMS3") & !grepl(x=coeffs.2$coefficient, pattern=":")
                col <- paste0(cms3mt.vs.cms2mt.col, ".pval")            
                tbl[sti, col] <- coeffs.2[flag, 5]
                flag <- grepl(x=coeffs.2$coefficient, pattern="CMS1") & !grepl(x=coeffs.2$coefficient, pattern=":")
                col <- paste0(cms1mt.vs.cms2mt.col, ".pval")            
                tbl[sti, col] <- coeffs.2[flag, 5]
                flag <- grepl(x=coeffs.2$coefficient, pattern="NOLBL") & !grepl(x=coeffs.2$coefficient, pattern=":")
                col <- paste0(nolblmt.vs.cms2mt.col, ".pval")
                tbl[sti, col] <- coeffs.2[flag, 5]
            }
            
            ## Only output the anova for the interaction
            if(grepl(x=formula.names[i], pattern="interaction")) { 
                anova.df <- as.data.frame(anova(lm.obj))
                anova.file <- paste(st, "-", analysis.name, "-", formula.names[i], "-anova.tsv", sep="")
                write.table(file=anova.file, anova.df, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
            }
        }
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")        
        tbl[,apcol] <- p.adjust(tbl[,pcol], method="BH")
    }
    
    ## Creat the individual expression plots for each gene/gene set
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", kras.states[1], " vs ", kras.states[2], " for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        ## Create a data frame where each row is a sample and the columns are:
        ## (1) expression of the gene/gene set
        ## (2) CMS factor of the sample
        ## (3) KRAS mutation status of the sample
        df <- data.frame(esn=esn.st, cms=factor(cms), KRAS=factor(kras.label))

        ## Add a final column that has _all_ samples (independent of CMS label)
        df.all <- df
        df.all$cms <- rep("ALL", nrow(df.all))
        df <- rbind(df, df.all)

        ## Exclude the samples that do not have a CMS label (but not from their
        ## contribution to the all column)
        ## df <- df[df$cms != "NOLBL",]
        df$cms <- factor(df$cms)

        ## Set the random seed so that jittering is always the same.
        myseed <- 1234
        set.seed(myseed)

        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste(st, "-", analysis.name, "-", kras.states[1], "-vs-", kras.states[2], "-box-faceted.pdf", sep="")
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=esn))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + geom_jitter()

        ## ... faceted on CMS label.
        p <- p + facet_grid(~cms)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits) - min(ylimits) )
        ## And we will add as many as 3 error bars per facet.  So adjust the y axis
        ## to accommodate this shift.
        ylimits[2] <- ylimits[2] + 4 * 3 * yoffset
        p <- p + ylim(ylimits)

        gb <- ggplot_build(p)
        g <- ggplot_gtable(gb)

        ## ggplot2 doesn't use native units in data space
        ## instead, the data is rescaled to npc, i.e., from 0 to 1
        ## so we need to use the build info to convert from data to [0,1]
        ## ranges <- gb$panel$ranges
        ranges <- gb$layout$panel_ranges
        ## panel.maxs will track the maximum value currently plotted in each facet.
        ## This will allow us to know where to position the error bar in each facet.
        mt.panel.maxs <- sapply(c(paste0("CMS",1:4), "NOLBL", "ALL"), function(cmsLbl) {
            mask <- df$cms == cmsLbl & df$KRAS == kras.states[1]
            m <- max(df$esn[mask])
            m
        })
        names(mt.panel.maxs) <- sapply(paste0(c(paste0("CMS",1:4), "NOLBL", "ALL")), function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(c(paste0("CMS",1:4), "NOLBL", "ALL"), function(cmsLbl) {
            mask <- df$cms == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$esn[mask])
            m
        })
        names(wt.panel.maxs) <- sapply(paste0(c(paste0("CMS",1:4), "NOLBL", "ALL")), function(x) paste0(x, "-", kras.states[2]))
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)

        print(panel.maxs)
        ## Add the error bars between KRAS MT and KRAW WT expression values (within a facet)
        facets <- c(paste0("CMS",1:4), "NOLBL", "ALL")
        for(i in 1:length(facets)) {
            if( ( i == 2 ) || ( i == 6 ) ) {
                ap.col <- paste0(cols[i], ".apval")
                tmp <- draw.mt.wt.err.bar(st, df, tbl, ranges, panel.maxs, g, facets[i], i, ap.col, kras.states)
                g <- tmp[["g"]]
                panel.maxs <- tmp[["panel.maxs"]]
            }

        }
        
        ## Add the error bars between CMS2 MT and { CMS1-3 MT } -- i.e., across multiple facets.
        ## This follows from http://stackoverflow.com/questions/31690007/ggplot-drawing-line-between-points-across-facets
        lbl <- "CMS2mut_vs_CMS1mut"
        lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
        lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
        lbl <- cms1mt.vs.cms2mt.col
        tmp <- draw.mt.err.bar(st, df, tbl, ranges, panel.maxs, g, "CMS1", "CMS2", 1, 2, lbl, kras.states)
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]
        
        lbl <- "CMS2mut_vs_CMS3mut"
        lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
        lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
        lbl <- cms3mt.vs.cms2mt.col
        tmp <- draw.mt.err.bar(st, df, tbl, ranges, panel.maxs, g, "CMS2", "CMS3", 2, 3, lbl, kras.states)
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]

        lbl <- "CMS2mut_vs_CMS4mut"
        lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
        lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
        lbl <- cms4mt.vs.cms2mt.col
        tmp <- draw.mt.err.bar(st, df, tbl, ranges, panel.maxs, g, "CMS2", "CMS4", 2, 4, lbl, kras.states)        
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]

        lbl <- nolblmt.vs.cms2mt.col
        tmp <- draw.mt.err.bar(st, df, tbl, ranges, panel.maxs, g, "CMS2", "NOLBL", 2, 5, lbl, kras.states)        
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()
    })

    write.table(file=paste0(analysis.name, "-tbl.xls"), tbl, sep="\t", row.names=FALSE, col.names=TRUE) 
    return(list(tbl=tbl, kras=kras,cms=cms,es=esn))
}

to.plot=c(immune.sets, "circ", "mek", "tak1", "myc.sig", "wnt", "BATF3", "ITGAE", "ATF3", "CCL4", "STAT1", "IRF1", "CXCL10", "CIITA", "CD3E", "CD4", "CD8A", "PML", "PIAS1", "MYC")
## to.plot=c("circ", "mek", "tak1", "myc", "wnt")
ylabels <- paste(to.plot, "Expression", sep=" ")
for(i in 1:(length(immune.sets) + 5)) {
    ylabels[i] <- paste0(toupper(to.plot[i]), " Enrichment Score")
}

do.harmonized.analysis <- FALSE
if(do.harmonized.analysis) {
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

do.tcga.analysis <- FALSE
do.kfs.analysis <- FALSE
do.rasness.analysis <- TRUE

tcga_expr <- NULL
kfsyscc_expr <- NULL
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
}
  
if(do.tcga.analysis && FALSE) {
    cat("Doing TCGA analysis\n")
    tcgaR <- doAnalysis(tcga_expr, clin, "tcga-mt", to.plot=to.plot, ylabels=ylabels)
    
    ## rm(tcga_expr)
    write.table(tcgaR$tbl1, file="Middleton_tcga_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(tcgaR$tbl2, file="Middleton_tcga_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
}


do.french.analysis <- FALSE
french_expr <- NULL
if(do.french.analysis) {
    cat("Reading French data\n")
    french_expr <- doaffy("syn2363561")
    ## french_expr <- doaffy(read.table(opt$`french-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
    
#    cat("Doing French analysis\n")
#    frenchR <- doAnalysis(french_expr, clin, "french-mt", to.plot=to.plot, ylabels=ylabels)
    
#    write.table(frenchR$tbl1, file="Middleton_french_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
#    write.table(frenchR$tbl2, file="Middleton_french_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

if(do.kfs.analysis && FALSE) {
    cat("Doing KFS analysis\n")
    kfsR <- doAnalysis(kfsyscc_expr, clin, "kfs-mt", to.plot=to.plot, ylabels=ylabels)

    ## rm(kfsyscc_expr)
    write.table(kfsR$tbl2, file="Middleton_kfs_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(kfsR$tbl1, file="Middleton_kfs_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

if(do.rasness.analysis) {
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
    clin.kras <- clin[!is.na(clin$kras),]
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

        clin$rasness <- rep(NA, nrow(clin))
        clin$score.rasness <- rep(NA, nrow(clin))    
        idxs <- match(clin$sample, rownames(test.harm.scores))
        clin$rasness[!is.na(idxs)] <- test.harm.scores[na.omit(idxs)] > cut.pt
        clin$score.rasness[!is.na(idxs)] <- test.harm.scores[na.omit(idxs)]
        
        pdf(paste0(test.set, "-harm-", train.set, "-trained-", optimize.set, "-optimized-rasness-vs-mt.pdf"))          
        df <- data.frame(score=clin$score.rasness, KRAS=as.factor(clin$kras))
        df <- df[!is.na(df$score) & !is.na(df$KRAS),]
        wilcox.test(score ~ KRAS, data = df)  
        p <- ggplot(data=df, aes(x=KRAS, y=score))
        p <- p + ggtitle(paste0("Harmonized ", test.set, " RASness (trained on ", train.set, "; optimized on ", optimize.set, ") vs KRAS MT"))          
        p <- p + ylab("RIS")
    
        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + geom_jitter()
        hline.data <- data.frame(y = cut.pt)
        p <- p + geom_hline(aes(yintercept = y), hline.data)
        print(p)
        d <- dev.off()

        ## Perform analysis based on RAS mutant
        cat(paste0("Performing harmonized mut-based analysis for ", test.set, "\n"))
        test.harm.mut.R <- doAnalysis(expr.test.harm, paste0(test.set, "-harm-mt"), to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))
        
        ## Perform analysis based on rasness
        cat(paste0("Performing harmonized rasness-based analysis for ", test.set, "\n"))        
        test.harm.rasness.R <- doAnalysis(expr.test.harm, paste0(test.set, "-harm-rasness"), to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness", kras.states=c("High", "Low"))

        ## Perform analysis based on rasness + RAS mutant
        cat(paste0("Performing harmonized mut+rasness-based analysis for ", test.set, "\n"))        
        clin$rasness.and.mut <- rep(NA, nrow(clin))
        flag <- unlist(apply(clin[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 1) && (row[2] == 1), FALSE)))
        clin$rasness.and.mut[flag] <- 1
        flag <- unlist(apply(clin[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 0) && (row[2] == 0), FALSE)))
        clin$rasness.and.mut[flag] <- 0
        flag <- !is.na(clin$braf) & (clin$braf == 1)
        clin$rasness.and.mut[flag] <- NA
        flag <- !is.na(clin$nras) & (clin$nras == 1)
        clin$rasness.and.mut[flag] <- NA
        
        test.harm.rasness.and.mut.R <- doAnalysis(expr.test.harm, paste0(test.set, "-harm-rasness-and-mut"), to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.and.mut", kras.states=c("High-MT", "Low-WT"))
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
