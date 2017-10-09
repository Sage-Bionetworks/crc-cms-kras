suppressPackageStartupMessages(library("annotate"))
suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("hgu133plus2.db"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gtable"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("corrplot"))

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
    data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    mask <- !is.na(symbol)
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    expr <- combine_probes_2_gene(data.set.m,symbol.m)
    colnames(expr) <- gsub("(.*?)_.*","\\1",colnames(expr))
    expr
}

data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]

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
draw.mt.wt.err.bar <- function(st, df, tbl1b, ranges, panel.maxs, g, cmsLbl, cmsIndx, kras.states) {

    ## Select out the rows corresponding to this facet/CMS.
    mask <- df$cms == cmsLbl

    ## Get the adjusted pvalue (comparing MT vs WT samples) for expression
    ## of st (a gene or gene set) for this CMS cluster.
    adj.pval.col <- paste(toString(cmsLbl), ".apval", sep="")
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
    t1 <- 4
    l1 <- 2 + 2 * cmsIndx
    t2 <- 4
    l2 <- 2 + 2 * cmsIndx

    ## Give the grob a random/unique name.  I don't know whether this is
    ## required.
    name=stringi::stri_rand_strings(1,5)
    ## I don't know why this delta is necessary--ggplot2 seems to be
    ## incorrectly confusing/reordering the grobs if I do not make them unique.
    delta <- runif(n=1, min=10^-5, max=10^-3)

    ## Set the current position to the start of the vertical line.
    g <- gtable_add_grob(g, grid.move.to(start[1],start[2],draw=FALSE,name=name), z = Inf, t = 4 + delta, l = l1, b = 4, r = l1)
    
    ## Draw line from the current position to the end of the vertical line
    ## (and set that end point as the current position)
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = 4 + delta, l = l1, r = l1, b = 4)

    ## Similarly, draw the horizontal line of the error bar.
    start <- c(data2npc(1,ranges[[cmsIndx]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    end <- c(data2npc(2,ranges[[cmsIndx]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = 4 + delta, l = l2, r = l2, b = 4)

    ## Finally, draw the vertical line on the "right side" of the error bar--
    ## this goes from the error bar to the maximum of the WT KRAS
    ## expression values.
    start <- c(data2npc(2,ranges[[cmsIndx]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    end <- c(data2npc(2,ranges[[cmsIndx]][["x.range"]]),
             data2npc(wt.max + 1 * yoffset, ranges[[cmsIndx]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = 4 + delta, l = l2, r = l2, b = 4)
    
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
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2]), t = 4 + delta, l = l1, b = 4, r = l2)
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
    t1 <- 4
    l1 <- 2 + 2 * cmsIndx1
    t2 <- 4
    l2 <- 2 + 2 * cmsIndx2

    ## Give the grob a random/unique name.  I don't know whether this is
    ## required.
    name=stringi::stri_rand_strings(1,5)
    ## I don't know why this delta is necessary--ggplot2 seems to be
    ## incorrectly confusing/reordering the grobs if I do not make them unique.
    delta <- runif(n=1, min=10^-5, max=10^-3)

    ## Set the current position to the start of the vertical line.
    g <- gtable_add_grob(g, grid.move.to(start[1],start[2],draw=FALSE,name=name), z = Inf, t = 4 + delta, l = l1, b = 4, r = l1)

    ## Draw line from the current position to the end of the vertical line
    ## (and set that end point as the current position)
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = 4 + delta, l = l1, r = l1, b = 4)

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
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = 4 + delta, l = l2, r = l2, b = 4)
    
    start <- c(data2npc(1,ranges[[cmsIndx2]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    end <- c(data2npc(1,ranges[[cmsIndx2]][["x.range"]]),
             data2npc(m + 1 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = 4 + delta, l = l2, r = l2, b = 4)

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
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2]), t = 4 + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
}

## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
doAnalysis <- function(expr, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.states=c("MT","WT")) {

    ## Possibly exclude braf and nras mutants
    exclude.braf.and.nras.mutants <- TRUE
    if(exclude.braf.and.nras.mutants) {
        
    }

    ## Exclude any samples that do not have KRAS mutation status annotated
    mask <- !is.na(clin[,kras.status.field])
    clin.mask <- clin[mask,]
    expr.mask <- expr[, mask]
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin.mask$sample, colnames(expr.mask))
    clin.m <- clin.mask[!is.na(idxs),]
    expr.m <- expr.mask[, na.omit(idxs)]
    
    ## gene set prep
    env <- new.env()
    load("markersG.RData",envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
    bindeaTbl <- read.table("bindea_immune.txt",sep='\t',header=TRUE)
    
    ## Restrict to a few cell types of interest
    bindeaTbl <- bindeaTbl[bindeaTbl$CellType %in% c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells"),]
    
    bindeaCellTypes <- unique(bindeaTbl$CellType)
    bindeaGsets <- lapply(bindeaCellTypes, function(x){
        return(as.character(unique(bindeaTbl$Symbol[bindeaTbl$CellType==x])))
    })
    names(bindeaGsets) <- bindeaCellTypes
    bindeaGsets <- bindeaGsets[sapply(bindeaGsets, length) > 10]
    bindeaGsets[["circ"]] <- read.table("circ_sig.txt",as.is=TRUE)[,1]
    bindeaGsets[["mek"]] <- c("DUSP6","PHLDA1","SPRY2","DUSP4","ETV4")
    # This TAK1 signature was read off of Fig 4A of Singh et al (2012) Cell
    tak1.signature <- c("RBP1", "SEMA3A", "SYT1", "EMR2", "PROX1", "INHBB", "ABHD2", "C1orf116", "SNTB1", "TAF9B", "PRF1", "SLC2A1", "GAD1", "MSX2", "PELI2", "ITGB4", "C21orf96", "GPR56", "PDK3", "GLS", "ACSL1", "BIK", "RUNX1", "SYK", "RGL1", "NAV2", "FYN", "HSPA12A", "MBOAT2", "BAMBI", "BMP7", "GGH")
    bindeaGsets[["tak1"]] <- tak1.signature
    
    # Wnt target genes from PMID 17320548, as used by J. Guinney in his Nat Med pub.  
    # wnt.targets <- c("ASCL2", "AXIN2", "BMP4", "C1orf33", "HIG2", "HSPC111", "KITLG", "LGR5", "MYC", "NOL1", "PPIF", "SOX4", "WRD71", "ZIC2", "ZNRF3")
    wnt.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "WNT_FLIER")]]
    bindeaGsets[["wnt"]] <- wnt.targets

    # Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    # myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
    myc.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "MYC_TARGETS_ZELLER")]]
    bindeaGsets[["myc"]] <- myc.targets
    
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
    cms <- as.character(clin.m$cms_label[mask])
    esn <- tmp[, mask]
    expr.mask <- expr.m[,mask]
    } ekse {
    kras <- clin.m[, kras.status.field]
    cms <- as.character(clin.m$cms_label)
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
    sigs2 <- c("circ", "circ", "myc")
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
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)
        
        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        cms.lbl <- "CMS2"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        cms.lbl <- "CMS3"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        cms.lbl <- "CMS4"
        kras.status <- 1
        kras.indx <- 1
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        kras.status <- 0
        kras.indx <- 2
        flag <- cms == cms.lbl & kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        cms.lbl <- "ALL"
        kras.status <- 1
        kras.indx <- 1
        flag <- kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)

        kras.status <- 0
        kras.indx <- 2
        flag <- kras == kras.status
        plot(esn[sig1,flag], esn[sig2,flag], xlab=toupper(sig1), ylab=toupper(sig2), main=paste0(cms.lbl, " ", kras.states[kras.indx], sep=""))
        lm.obj <- lm(esn[sig1,flag] ~ esn[sig2,flag])
        abline(lm.obj)
        
        d <- dev.off()
    }
    
    ## Analysis 4:  look at the correlation of corelation of circ values across CMS
    circ.genes <- bindeaGsets[["circ"]][ bindeaGsets[["circ"]] %in% rownames(expr.mask) ]
    if(length(circ.genes) > 1) {
        cms.labels <- c(as.vector(sapply(paste0("CMS",1:4), function(x) paste0(x, "-", kras.states[1]))),as.vector(sapply(paste0("CMS",1:4), function(x) paste0(x, "-", kras.states[2]))))
        circ.cms.mat <- c()
        for(i in 1:length(cms.labels)) {
            cmslbl <- cms.labels[i]
            cmslbl <- gsub(x=cmslbl, pattern=paste0("-", kras.states[1]), replacement="")
            cmslbl <- gsub(x=cmslbl, pattern=paste0("-", kras.states[2]), replacement="")
            kras.status <- ifelse(grepl(pattern=kras.states[1], x=cms.labels[i]), 1, 0)
            # print(c(cmslbl,kras.status))
            mask <- (cmslbl == cms) & (kras == kras.status)
            expr.subset <- expr.mask[circ.genes,mask]
            # print(nrow(expr.subset))
            subset.cor <- cor(t(expr.subset))
            vec <- subset.cor[upper.tri(subset.cor)]
            # print(length(vec))
            circ.cms.mat <- rbind(circ.cms.mat, vec)
        }
        rownames(circ.cms.mat) <- cms.labels
        res <- cor(t(circ.cms.mat))
        # print(circ.cms.mat)
        # print(res)

        out.pdf <- paste("circ", "-", analysis.name, "-cor.pdf", sep="")
        pdf(out.pdf)
        corrplot(res, type="upper", order="hclust", tl.col="black", tl.srt=45)
        d <- dev.off()
        
    }

    
    if(FALSE){
    save(clin.m, file="clin.Rdata")
    save(kras, file="kras.Rdata")
    save(esn, file="esn.Rdata")
    save(bindeaGsets, file="bind.Rdata")
    save(expr.mask, file="expr.Rdata")
    save(cms, file="cms.Rdata")

    ## Analysis 4:  look at the correlation of mean circ values across CMS
    ## Create a matrix of CMS by mean gene expression
    circ.genes <- bindeaGsets[["circ"]][ bindeaGsets[["circ"]] %in% rownames(expr.mask) ]
    if(length(circ.genes) != 0) {
        cms.labels <- c(as.vector(sapply(paste0("CMS",1:4), function(x) paste0(x, "-", kras.states[1]))),as.vector(sapply(paste0("CMS",1:4), function(x) paste0(x, "-", kras.states[2]))))
        circ.cms.mat <- matrix(data = NA, nrow=length(cms.labels), ncol=length(circ.genes))
        rownames(circ.cms.mat) <- cms.labels
        colnames(circ.cms.mat) <- circ.genes
        for(i in 1:nrow(circ.cms.mat)) {
            cmslbl <- rownames(circ.cms.mat)[i]
            cmslbl <- gsub(x=cmslbl, pattern=paste0("-", kras.states[1]), replacement="")
            cmslbl <- gsub(x=cmslbl, pattern=paste0("-", kras.states[2]), replacement="")            
            kras.status <- ifelse(grepl(pattern=kras.states[1], x=rownames(circ.cms.mat)[i]), 1, 0)
            mask <- (cmslbl == cms) & (kras == kras.status)
            expr.subset <- expr.mask[,mask]
            print(nrow(expr.subset))
            for(j in 1:ncol(circ.cms.mat)) {
                gene <- colnames(circ.cms.mat)[j]
                print(gene)
                circ.cms.mat[i,j] <- mean(expr.subset[gene,])
            }
        }
        ## save(circ.cms.mat, file="circ.cms.mat.Rdata")
        res <- cor(t(circ.cms.mat))
        print(circ.cms.mat)
        print(res)

        out.pdf <- paste("circ", "-", analysis.name, "-cor.pdf", sep="")
        pdf(out.pdf)
        corrplot(res, type="upper", order="hclust", tl.col="black", tl.srt=45)
        d <- dev.off()
        
    }
    }
    ## analysis 1: compare mut vs wt within each CMS
    tbl1 <- do.call("cbind",lapply(paste0("CMS",1:4), function(cmsLbl){
        cat(paste("Computing ", kras.states[1], " vs ", kras.states[2], " for ", cmsLbl, "\n", sep=""))
        mask <- cmsLbl ==  cms
        pval <- apply(esn[,mask], 1, function(x) wilcox.test(x ~ kras[mask])$p.value)
        a <- p.adjust(pval,method="BH")
        b <- apply(esn[,mask], 1, function(x) (mean(x[kras[mask] == 1]) - mean(x[kras[mask] == 0])) > 0)
        cbind(pval=pval,apval=a, mutUp=b)
    }))

    colnames(tbl1) <- do.call(c, lapply(paste0("CMS",1:4), function(x) paste(x, c(".pval", ".apval", ".mutUp"), sep="")))
    
    ## analysis 1b: compare mut vs wt across ALL
    pval <- apply(esn, 1, function(x) wilcox.test(x ~ kras)$p.value)
    a <- p.adjust(pval,method="BH")
    b <- apply(esn, 1, function(x) (mean(x[kras == 1]) - mean(x[kras == 0])) > 0)
    tbl1b <- cbind(variable=rownames(tbl1), tbl1, cbind(ALL.pval=pval,ALL.apval=a,ALL.mutUp=b))
    
    ## analysis 2: compare CMS for ras mut
    
    doDiffCompare <- function(col.label.prefix, m1, m2){
        a <- t(apply(esn, 1, function(x) {
            pval <- wilcox.test(x[m1],x[m2])$p.value
            c(pval=pval,CMSup=(mean(x[m1]) - mean(x[m2]) > 0  ))
        }))
        ret <- cbind(pval=a[,1], apval=p.adjust(a[,1],method="BH"),CMSup=a[,2])
        colnames(ret) <- sapply(c(".pval",".apval",".CMSup"), function(x) paste(col.label.prefix, x, sep=""))
        return(ret)
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]
    
    lbl <- "CMS2mut_vs_CMS1mut"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t0 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, cms=="CMS1" & kras==1)    

    lbl <- "CMS2mut_vs_CMS3mut"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t1 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, cms=="CMS3" & kras==1)
    
    lbl <- "CMS2mut_vs_CMS4mut"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t2 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, cms=="CMS4" & kras==1)
    
    lbl <- "CMS2mut_vs_CMS34mut"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t3 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, (cms=="CMS3"|cms=="CMS4") & kras==1)
    
    lbl <- "CMS2wt_vs_CMS3wt"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t4 <- doDiffCompare(lbl, cms=="CMS2" & kras==0, cms=="CMS3" & kras==0)
    
    lbl <- "CMS2wt_vs_CMS3mut"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t5 <- doDiffCompare(lbl, cms=="CMS2" & kras==0, cms=="CMS3" & kras==1)
    
    lbl <- "CMS2wt_vs_CMS4wt"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t6 <- doDiffCompare(lbl, cms=="CMS2" & kras==0, cms=="CMS4" & kras==0)
    
    lbl <- "CMS2wt_vs_CMS4mut"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t7 <- doDiffCompare(lbl, cms=="CMS2" & kras==0, cms=="CMS4" & kras==1)
    
    lbl <- "CMS2mut_vs_CMS3"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t8 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, cms=="CMS3")
    
    lbl <- "CMS2mut_vs_CMS4"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t9 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, cms=="CMS4")
    
    lbl <- "CMS2mut_vs_CMS34"
    lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
    lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
    t10 <- doDiffCompare(lbl, cms=="CMS2" & kras==1, cms=="CMS3" | cms=="CMS4")
    
    tbl2 <- cbind(variable=rownames(t0), t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)

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
        df <- df[df$cms != "NOLBL",]
        df$cms <- factor(df$cms)

        ## Create a PDF for this gene/gene set (this is a non-faceted plot)
##        out.pdf <- paste(st, "-", analysis.name, "-mt-vs-wt-box.pdf", sep="")
##        pdf(out.pdf)

##        ## Create the plot
##        p <- ggplot(df, aes(x=cms, y=esn))
##        p <- p + ylab(ylab)
        
##        ## Hinges correspond to 25th and 75th percentile
##        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
##        print(p)
##        d <- dev.off()

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
        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range

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
        ranges <- gb$panel$ranges

        ## panel.maxs will track the maximum value currently plotted in each facet.
        ## This will allow us to know where to position the error bar in each facet.
        mt.panel.maxs <- sapply(c(paste0("CMS",1:4), "ALL"), function(cmsLbl) {
            mask <- df$cms == cmsLbl & df$KRAS == kras.states[1]
            m <- max(df$esn[mask])
            m
        })
        names(mt.panel.maxs) <- sapply(paste0(c(paste0("CMS",1:4), "ALL")), function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(c(paste0("CMS",1:4), "ALL"), function(cmsLbl) {
            mask <- df$cms == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$esn[mask])
            m
        })
        names(wt.panel.maxs) <- sapply(paste0(c(paste0("CMS",1:4), "ALL")), function(x) paste0(x, "-", kras.states[2]))
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)

        ## Add the error bars between KRAS MT and KRAW WT expression values (within a facet)
        facets <- c(paste0("CMS",1:4), "ALL")
        for(i in 1:length(facets)) {
            tmp <- draw.mt.wt.err.bar(st, df, tbl1b, ranges, panel.maxs, g, facets[i], i, kras.states)
            g <- tmp[["g"]]
            panel.maxs <- tmp[["panel.maxs"]]

        }
        
        ## Add the error bars between CMS2 MT and { CMS1-3 MT } -- i.e., across multiple facets.
        ## This follows from http://stackoverflow.com/questions/31690007/ggplot-drawing-line-between-points-across-facets
        lbl <- "CMS2mut_vs_CMS1mut"
        lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
        lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
        tmp <- draw.mt.err.bar(st, df, tbl2, ranges, panel.maxs, g, "CMS1", "CMS2", 1, 2, lbl, kras.states)
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]

        lbl <- "CMS2mut_vs_CMS3mut"
        lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
        lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
        tmp <- draw.mt.err.bar(st, df, tbl2, ranges, panel.maxs, g, "CMS2", "CMS3", 2, 3, lbl, kras.states)
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]

        lbl <- "CMS2mut_vs_CMS4mut"
        lbl <- gsub(x=lbl, pattern="mut", replacement=kras.pos)
        lbl <- gsub(x=lbl, pattern="wt", replacement=kras.neg)
        tmp <- draw.mt.err.bar(st, df, tbl2, ranges, panel.maxs, g, "CMS2", "CMS4", 2, 4, lbl, kras.states)        
        g <- tmp[["g"]]
        panel.maxs <- tmp[["panel.maxs"]]

        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()
        
    })

    ## Analysis 3: linear modeling of biomarker ~ CMS + KRAS + CMS:KRAS
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        ## Create a data frame where each row is a sample and the columns are:
        ## (1) expression of the gene/gene set
        ## (2) CMS factor of the sample
        ## (3) KRAS mutation status of the sample
        df <- data.frame(expr=esn.st, CMS=factor(cms), KRAS=factor(kras.label))

        df <- df[df$CMS != "NOLBL",]
        lm.obj <- lm(expr ~ CMS + KRAS + CMS:KRAS, data=df)
        lm.sum <- summary(lm.obj)
        print(lm.sum)

        coeffs <- coefficients(lm.sum)

        lm.file <- paste(st, "-", analysis.name, "-lm.tsv", sep="")
        coeffs.2 <- cbind(coefficient=rownames(coeffs), coeffs)
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        indices <- which(coeffs[,4] < 0.05)
        if(length(indices) > 0) {
            cat(paste(c(st, analysis.name, "covariate", colnames(coeffs), "SIGNIFICANT", "\n"), collapse="\t"))
            d <- unlist(lapply(indices, function(i) cat(paste(c(st, analysis.name, rownames(coeffs)[i], coeffs[i,], "SIGNIFICANT", "\n"), collapse="\t"))))
        }
    })

    return(list(tbl1=tbl1b, tbl2=tbl2, kras=kras,cms=cms,es=esn))
}

to.plot=c("circ", "mek", "tak1", "myc", "wnt", "BATF3", "ITGAE", "ATF3", "CCL4", "STAT1", "IRF1", "CXCL10", "CIITA", "CD3E", "CD4", "CD8A", "PML", "PIAS1", "MYC")
ylabels <- paste(to.plot, "Expression", sep=" ")
ylabels[1] <- "CIRC Enrichment Score"
ylabels[2] <- "MEK Enrichment Score"
ylabels[3] <- "TAK1 Enrichment Score"
ylabels[4] <- "MYC Enrichment Score"
ylabels[5] <- "WNT Enrichment Score"

do.harmonized.analysis <- TRUE
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
    
    ## harmonized_expr <- t(read.table(opt$`harmonized-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))

    cat("Doing harmonized analysis\n")
    harmonizedR <- doAnalysis(harmonized_expr, "harmonized-mt", to.plot=to.plot, ylabels=ylabels)

    ## rm(harmonized_expr)
    write.table(harmonizedR$tbl1, file="Middleton_harmonized_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(harmonizedR$tbl2, file="Middleton_harmonized_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

do.rppa.analysis <- FALSE
if(do.rppa.analysis) {
    cat("Reading RPPA data\n")
    rppa.file <- "data_rppa.txt"
    rppa_expr <- {
        rppa_expr <- read.table(opt$`rppa-expr-file`, sep="\t", header=TRUE, check.names=FALSE)
        ## Drop phospho probes and cleaved PARP--that should leave just one
        ## probe per gene
        flag <- grepl(pattern="_p", x=rppa_expr[,1]) | grepl(pattern="cleaved", x=rppa_expr[,1])
        rppa_expr <- rppa_expr[!flag,]
        proteins <- rppa_expr[,1]
        ## Parse out the protein name from the probe.  These are annotated as
        ## protein|probe
        proteins <- sapply(proteins, function(x) {
            r <- regexpr(pattern="^[^\\|]+\\|", text=x)
            if(r==-1) {
                cat(paste("Found no match for | in: ", x, "\n", sep=""))
                q(status=-1)
            }
            substr(x, r[1], r[1] + attr(r,"match.length")[1]-2)
        })
        if(any(duplicated(proteins))) {
            cat("Multiple RPPA expression values for the same protein!\n")
            flag <- duplicated(proteins, fromLast=TRUE) | duplicated(proteins, fromLast=FALSE)
            print(rppa_expr[flag,1])
            q(status=-1)
        }
        expr <- as.matrix(rppa_expr[,-1])
        rownames(expr) <- proteins
        
        ## Translate the column names from TCGA-G4-6317-01 -> TCGA-G4-6317
        cols <- colnames(expr)
        cols <- gsub(pattern="-\\d+1$", replacement="", x=cols)
        colnames(expr) <- cols
        expr
    }
    
    cat("Doing RPPA analysis\n")
    rppaR <- doAnalysis(rppa_expr, "rppa-mt", to.plot=to.plot, ylabels=ylabels)
    write.table(rppaR$tbl1, file="Middleton_rppa_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(rppaR$tbl2, file="Middleton_rppa_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

do.tcga.analysis <- TRUE
do.kfs.analysis <- TRUE

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
  
if(do.tcga.analysis) {
    cat("Doing TCGA analysis\n")
    tcgaR <- doAnalysis(tcga_expr, "tcga-mt", to.plot=to.plot, ylabels=ylabels)
    
    ## rm(tcga_expr)
    write.table(tcgaR$tbl1, file="Middleton_tcga_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(tcgaR$tbl2, file="Middleton_tcga_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

do.french.analysis <- FALSE
if(do.french.analysis) {
    cat("Reading French data\n")
    french_expr <- doaffy("syn2363561")
    ## french_expr <- doaffy(read.table(opt$`french-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
    
    cat("Doing French analysis\n")
    frenchR <- doAnalysis(french_expr, "french-mt", to.plot=to.plot, ylabels=ylabels)
    
    write.table(frenchR$tbl1, file="Middleton_french_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(frenchR$tbl2, file="Middleton_french_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

if(do.kfs.analysis) {
    cat("Doing KFS analysis\n")
    kfsR <- doAnalysis(kfsyscc_expr, "kfs-mt", to.plot=to.plot, ylabels=ylabels)

    ## rm(kfsyscc_expr)
    write.table(kfsR$tbl2, file="Middleton_kfs_tbl2.xls",sep="\t",quote=FALSE,row.names=FALSE)
    write.table(kfsR$tbl1, file="Middleton_kfs_tbl1.xls",sep="\t",quote=FALSE,row.names=FALSE)
}

do.harmonized.rasness.analysis <- TRUE
if(do.harmonized.rasness.analysis) {
  
    ## Fit elastic net to the _harmonized_ TCGA data
    ## expr <- harmonized_expr[,grepl(pattern="TCGA", x=colnames(harmonized_expr))]
    expr <- harmonized_expr[, colnames(harmonized_expr) %in% colnames(tcga_expr)]    
    clin.kras <- clin[!is.na(clin$kras),]
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.tcga.harm <- clin.kras[!is.na(idxs),]
    expr.tcga.harm <- expr[, na.omit(idxs)]
    
    x.tcga <- t(expr.tcga.harm)
    x.tcga <- scale(x.tcga)
    
    alphas <- seq(from = 0, to = 1, by = 0.05)
    nfolds <- 5
    N <- length(clin.tcga.harm$kras)
    foldid <- sample(rep(seq(nfolds), length = N))
    tcga.harm.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.tcga, y = as.factor(clin.tcga.harm$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)
    
    ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
    pdf("tcga-harm-aucs.pdf")  
    aucs <- as.data.frame(lapply(tcga.harm.models, function(x) { return(x$cvm) }))
    colnames(aucs) <- 1:ncol(aucs)
    heatmap.2(as.matrix(aucs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
  
    ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
    best.aucs <- unlist(lapply(tcga.harm.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
    plot(alphas, best.aucs, main="lambda.1se")

    best.aucs <- unlist(lapply(tcga.harm.models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
    plot(alphas, best.aucs, main="lambda.min")
    d <- dev.off()
    
    ## Best AUCs are achieved for lambdas in the range 0.05 to 0.15.
    ## best.alpha <- alphas[which(best.aucs == max(best.aucs))]
    ## Let's just use alpha = 0.1, which is what Justin published.
    best.alpha <- 0.1
    best.tcga.harm.trained.model <- tcga.harm.models[[which(alphas == best.alpha)]]

    ## Genes in our model
    tcga.harm.genes.min <- rownames(expr.tcga.harm)[which(coef(best.tcga.harm.trained.model, s="lambda.min") != 0)]
    tcga.harm.genes.1se <- rownames(expr.tcga.harm)[which(coef(best.tcga.harm.trained.model, s="lambda.1se") != 0)]  
  
    ## Fit elastic net to the _harmonized_ KFS data
    ##  expr <- harmonized_expr[,grepl(pattern="COLON", x=colnames(harmonized_expr))]
    expr <- harmonized_expr[, colnames(harmonized_expr) %in% colnames(kfsyscc_expr)]
    clin.kras <- clin[!is.na(clin$kras),]
    idxs <- match(clin.kras$sample, colnames(expr))
    clin.kfs.harm <- clin.kras[!is.na(idxs),]
    expr.kfs.harm <- expr[, na.omit(idxs)]
  
    alphas <- seq(from = 0, to = 1, by = 0.05)
    nfolds <- 5
    N <- length(clin.kfs.harm$kras)
    foldid <- sample(rep(seq(nfolds), length = N))
    x.kfs <- t(expr.kfs.harm)
    x.kfs <- scale(x.kfs)
    kfs.harm.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.kfs, y = as.factor(clin.kfs.harm$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)
  
    ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
    pdf("kfs-harm-aucs.pdf")  
    aucs <- as.data.frame(lapply(kfs.harm.models, function(x) { return(x$cvm) }))
    colnames(aucs) <- 1:ncol(aucs)
    heatmap.2(as.matrix(aucs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
  
    ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
    best.aucs <- unlist(lapply(kfs.harm.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
    plot(alphas, best.aucs)

    best.aucs <- unlist(lapply(kfs.harm.models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
    plot(alphas, best.aucs)
    d <- dev.off()
    
    ## Best AUCs are achieved for lambdas in the range 0.05 to 0.15.
    ## best.alpha <- alphas[which(best.aucs == max(best.aucs))]
    ## Let's just use alpha = 0.1, which is what Justin published.
    best.alpha <- 0.1
    best.kfs.harm.trained.model <- kfs.harm.models[[which(alphas == best.alpha)]]

    ## Genes in our model
    kfs.harm.genes.min <- rownames(expr.kfs.harm)[which(coef(best.kfs.harm.trained.model, s="lambda.min") != 0)]
    kfs.harm.genes.1se <- rownames(expr.kfs.harm)[which(coef(best.kfs.harm.trained.model, s="lambda.1se") != 0)]  

    ## Let's look at the intersection of KFS.harm and TCGA.harm genes
    print(intersect(tcga.harm.genes.min, kfs.harm.genes.min))
    print(intersect(tcga.harm.genes.1se, kfs.harm.genes.1se))

    ## Let's look at the performance of each test set
    kfs.harm.scores <- predict(best.tcga.harm.trained.model, newx=x.kfs, s="lambda.min", type="response")
    tcga.harm.scores <- predict(best.kfs.harm.trained.model, newx=x.tcga, s="lambda.min", type="response")

    clin$rasness.harm.kfs <- rep(NA, nrow(clin))
    clin$score.rasness.harm.kfs <- rep(NA, nrow(clin))    
    idxs <- match(clin$sample, rownames(kfs.harm.scores))
    clin$rasness.harm.kfs[!is.na(idxs)] <- kfs.harm.scores[na.omit(idxs)] > 0.5
    clin$score.rasness.harm.kfs[!is.na(idxs)] <- kfs.harm.scores[na.omit(idxs)]
    
    clin$rasness.harm.tcga <- rep(NA, nrow(clin))
    clin$score.rasness.harm.tcga <- rep(NA, nrow(clin))    
    idxs <- match(clin$sample, rownames(tcga.harm.scores))
    clin$rasness.harm.tcga[!is.na(idxs)] <- tcga.harm.scores[na.omit(idxs)] > 0.5
    clin$score.rasness.harm.tcga[!is.na(idxs)] <- tcga.harm.scores[na.omit(idxs)]
    
    table(clin$rasness.harm.tcga, clin$kras)
    
    tcga.harm.rasness.R <- doAnalysis(expr.tcga.harm, "tcga-harm-rasness", to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.harm.tcga", kras.states=c("High", "Low"))

    kfs.harm.rasness.R <- doAnalysis(expr.kfs.harm, "kfs-harm-rasness", to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.harm.kfs", kras.states=c("High", "Low"))
    
    pdf("kfs-harm-rasness-vs-mt.pdf")      
    df <- data.frame(score=clin$score.rasness.harm.kfs, KRAS=as.factor(clin$kras))
    df <- df[!is.na(df$score) & !is.na(df$KRAS),]
    wilcox.test(score ~ KRAS, data = df)
    p <- ggplot(data=df, aes(x=KRAS, y=score))
    p <- p + ylab("RIS")
    
    ## Create a box plot where the x axis is KRAS mutation status ...
    p <- p + geom_boxplot(aes(fill=KRAS))
    p <- p + geom_jitter()
    print(p)
    d <- dev.off()

    pdf("tcga-harm-rasness-vs-mt.pdf")          
    df <- data.frame(score=clin$score.rasness.harm.tcga, KRAS=as.factor(clin$kras))
    df <- df[!is.na(df$score) & !is.na(df$KRAS),]
    wilcox.test(score ~ KRAS, data = df)  
    p <- ggplot(data=df, aes(x=KRAS, y=score))
    p <- p + ylab("RIS")
    
    ## Create a box plot where the x axis is KRAS mutation status ...
    p <- p + geom_boxplot(aes(fill=KRAS))
    p <- p + geom_jitter()
    print(p)
    d <- dev.off()
    
    # kfs.harm.scores <- predict(best.tcga.harm.trained.model, newx=t(expr.kfs.harm), s="lambda.1se", type="class")  
    # tcga.harm.scores <- predict(best.kfs.harm.trained.model, newx=t(expr.tcga.harm), s="lambda.1se", type="class")    
}

do.rasness.analysis <- TRUE
if(do.rasness.analysis) {
  kfs.tcga.inters <- intersect(rownames(kfsyscc_expr), rownames(tcga_expr))
  
  ## Fit elastic net to TCGA data
  expr <- tcga_expr
  ## Since we will train on TCGA, but apply to KFS, restrict to genes that are in KFS.
  expr <- expr[kfs.tcga.inters,]
  clin.kras <- clin[!is.na(clin$kras),]
  idxs <- match(clin.kras$sample, colnames(expr))
  clin.tcga <- clin.kras[!is.na(idxs),]
  expr.tcga <- expr[, na.omit(idxs)]
  
  alphas <- seq(from = 0, to = 1, by = 0.05)
  nfolds <- 5
  N <- length(clin.tcga$kras)
  foldid <- sample(rep(seq(nfolds), length = N))
  x.tcga <- t(expr.tcga)
  x.tcga <- scale(x.tcga)
  tcga.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.tcga, y = as.factor(clin.tcga$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)
  
  ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
    pdf("tcga-aucs.pdf")    
  aucs <- as.data.frame(lapply(tcga.models, function(x) { return(x$cvm) }))
  colnames(aucs) <- 1:ncol(aucs)
  heatmap.2(as.matrix(aucs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
  
  ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
  best.aucs <- unlist(lapply(tcga.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
  plot(alphas, best.aucs)

  best.aucs <- unlist(lapply(tcga.models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
  plot(alphas, best.aucs)
  d <- dev.off()
  
  ## Best AUCs are achieved for lambdas in the range 0.05 to 0.15.
  ## best.alpha <- alphas[which(best.aucs == max(best.aucs))]
  ## Let's just use alpha = 0.1, which is what Justin published.
  best.alpha <- 0.1
  best.tcga.trained.model <- tcga.models[[which(alphas == best.alpha)]]

  ## Genes in our model
  tcga.genes.min <- rownames(expr.tcga)[which(coef(best.tcga.trained.model, s="lambda.min") != 0)]
  tcga.genes.1se <- rownames(expr.tcga)[which(coef(best.tcga.trained.model, s="lambda.1se") != 0)]  
  
  ## Fit elastic net to KFS data
  expr <- kfsyscc_expr
  ## Since we will train on KFS, but apply to TCGA, restrict to genes that are in TCGA.
  expr <- expr[kfs.tcga.inters,]
  clin.kras <- clin[!is.na(clin$kras),]
  idxs <- match(clin.kras$sample, colnames(expr))
  clin.kfs <- clin.kras[!is.na(idxs),]
  expr.kfs <- expr[, na.omit(idxs)]
  
  alphas <- seq(from = 0, to = 1, by = 0.05)
  nfolds <- 5
  N <- length(clin.kfs$kras)
  foldid <- sample(rep(seq(nfolds), length = N))
  x.kfs <- t(expr.kfs)
  x.kfs <- scale(x.kfs)
  kfs.models <- mclapply(alphas, fit.elastic.net.alpha, x = x.kfs, y = as.factor(clin.kfs$kras), foldid = foldid, nfolds = nfolds, mc.cores = num.processes)
  
  ## Generate a heatmap of mean cross-validation error (cvm) as a function of alpha and lambda _index_ (not lambda)
  pdf("kfs-aucs.pdf")      
  aucs <- as.data.frame(lapply(kfs.models, function(x) { return(x$cvm) }))
  colnames(aucs) <- 1:ncol(aucs)
  heatmap.2(as.matrix(aucs), Rowv = F, Colv = F, scale = "none", trace = "none", dendrogram = "none", labCol = alphas, xlab = "alpha", ylab = "lambda index")
  
  ## Find the max AUC based on lambda.1se (NB: this is returned in cvm)
  best.aucs <- unlist(lapply(kfs.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
  plot(alphas, best.aucs, main="lambda.1se")

  best.aucs <- unlist(lapply(kfs.models, function(x) x$cvm[which(x$lambda == x$lambda.min)]))
  plot(alphas, best.aucs, main="lambda.min")
  d <- dev.off()
  
  ## Best AUCs are achieved for lambdas in the range 0.05 to 0.15.
  ## best.alpha <- alphas[which(best.aucs == max(best.aucs))]
  ## Let's just use alpha = 0.1, which is what Justin published.
  best.alpha <- 0.1
  best.kfs.trained.model <- kfs.models[[which(alphas == best.alpha)]]

  ## Genes in our model
  kfs.genes.min <- rownames(expr.kfs)[which(coef(best.kfs.trained.model, s="lambda.min") != 0)]
  kfs.genes.1se <- rownames(expr.kfs)[which(coef(best.kfs.trained.model, s="lambda.1se") != 0)]  

  ## Let's look at the intersection of KFS and TCGA genes
  print(intersect(tcga.genes.min, kfs.genes.min))
  print(intersect(tcga.genes.1se, kfs.genes.1se))

  ## Let's look at the performance of each test set
  ## This is bimodal:
  ## kfs.scores <- predict(best.kfs.trained.model, newx=x.kfs, s="lambda.min", type="response")  
  kfs.scores <- predict(best.tcga.trained.model, newx=x.kfs, s="lambda.min", type="response")
  tcga.scores <- predict(best.kfs.trained.model, newx=x.tcga, s="lambda.min", type="response")

#  kfs.scores <- predict(best.tcga.trained.model, newx=x.kfs, s="lambda.1se", type="response")
#  tcga.scores <- predict(best.kfs.trained.model, newx=x.tcga, s="lambda.1se", type="response")

  clin$rasness.kfs <- rep(NA, nrow(clin))
  idxs <- match(clin$sample, rownames(kfs.harm.scores))
  clin$rasness.kfs[!is.na(idxs)] <- kfs.harm.scores[na.omit(idxs)] > 0.5
  clin$score.rasness.kfs <- rep(NA, nrow(clin))
  clin$score.rasness.kfs[!is.na(idxs)] <- kfs.scores[na.omit(idxs)]
  
  clin$rasness.tcga <- rep(NA, nrow(clin))
  idxs <- match(clin$sample, rownames(tcga.harm.scores))
  clin$rasness.tcga[!is.na(idxs)] <- tcga.harm.scores[na.omit(idxs)] > 0.5
  clin$score.rasness.tcga <- rep(NA, nrow(clin))
  clin$score.rasness.tcga[!is.na(idxs)] <- tcga.scores[na.omit(idxs)]

  tcga.rasness.R <- doAnalysis(expr.tcga, "tcga-rasness", to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.tcga", kras.states=c("High", "Low"))
  kfs.rasness.R <- doAnalysis(expr.kfs, "kfs-rasness", to.plot=to.plot, ylabels=ylabels, kras.status.field="rasness.kfs", kras.states=c("High", "Low"))

  pdf("tcga-rasness-vs-mt.pdf")        
  df <- data.frame(score=clin$score.rasness.tcga, KRAS=as.factor(clin$kras))
  df <- df[!is.na(df$score) & !is.na(df$KRAS),]
  wilcox.test(score ~ KRAS, data = df)
  p <- ggplot(data=df, aes(x=KRAS, y=score))
  p <- p + ylab("RIS")
  
  ## Create a box plot where the x axis is KRAS mutation status ...
  p <- p + geom_boxplot(aes(fill=KRAS))
  p <- p + geom_jitter()
  print(p)
  d <- dev.off()

  pdf("kfs-rasness-vs-mt.pdf")          
  df <- data.frame(score=clin$score.rasness.kfs, KRAS=as.factor(clin$kras))
  df <- df[!is.na(df$score) & !is.na(df$KRAS),]
  wilcox.test(score ~ KRAS, data = df)  
  p <- ggplot(data=df, aes(x=KRAS, y=score))
  p <- p + ylab("RIS")
  
  ## Create a box plot where the x axis is KRAS mutation status ...
  p <- p + geom_boxplot(aes(fill=KRAS))
  p <- p + geom_jitter()
  print(p)
  d <- dev.off()

#  kfs.scores <- predict(best.tcga.trained.model, newx=x.kfs, s="lambda.1se", type="class")  
#  tcga.scores <- predict(best.kfs.trained.model, newx=x.tcga, s="lambda.1se", type="class")    
}

