## devtools::install_github("NikNakk/forestmodel")                              
library("forestmodel")
suppressPackageStartupMessages(library("ggbeeswarm"))

inverse.norm.transform <- function(x) {
##    p <- 2*pnorm(abs(x), lower.tail=FALSE) 
##    x2 <- qnorm(p/2, lower.tail=FALSE)*sign(x)
    ##    x2
    qq <- qqnorm(x, plot.it = FALSE)
    trn <- ( mean(x) + ( sd(x) * qq$x ) )
    trn
}

doCorrelationAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {
    circ.genes.in.set <- circ.genes[circ.genes %in% rownames(expr)]
    cor.matrix <- c()
    for(cms.grp in c("CMS1", "CMS2", "CMS3", "CMS4")) {
        for(kras.status in c(0,1)) {
            grp.flag <- !is.na(clin$cms_label) & !is.na(clin$kras) & (clin$cms_label == cms.grp) & (clin$kras == kras.status)
            samples.grp <- clin$sample[grp.flag]
            samples.grp <- samples.grp[samples.grp %in% colnames(expr)]
            expr.subset <- expr[circ.genes.in.set,samples.grp]
            pca <- prcomp(expr.subset)
            cat(paste0("Percent of variance explained: ", pca$sdev[1]/sum(pca$sdev), "\n"))
            pc1 <- pca$rotation[,1]
            cors <- as.matrix(unlist(apply(expr.subset, 1, function(row) cor(row, pc1))))
            cors <- as.data.frame(cors)
            ##        colnames(cors) <- paste(cms.grp, kras.status, sep="-")
            colnames(cors) <- "expr"
            cors$KRAS <- ifelse(kras.status == 0, "WT", "MT")
            cors$CMS <- cms.grp
            cors$gene <- rownames(cors)
            cor.matrix <- rbind(cor.matrix, cors)
        }
    }
    
    kras.cor.matrix <- c()
    for(kras.status in c(0,1)) {
        grp.flag <- !is.na(clin$kras) & (clin$kras == kras.status)
        samples.grp <- clin$sample[grp.flag]
        samples.grp <- samples.grp[samples.grp %in% colnames(expr)]
        expr.subset <- expr[circ.genes.in.set,samples.grp]
        pca <- prcomp(expr.subset)
        cat(paste0("Percent of variance explained: ", pca$sdev[1]/sum(pca$sdev), "\n"))
        pc1 <- pca$rotation[,1]
        cors <- as.matrix(unlist(apply(expr.subset, 1, function(row) cor(row, pc1))))
        cors <- as.data.frame(cors)
        ##        colnames(cors) <- paste(cms.grp, kras.status, sep="-")
        colnames(cors) <- "expr"
        cors$KRAS <- ifelse(kras.status == 0, "WT", "MT")
        cors$gene <- rownames(cors)
        kras.cor.matrix <- rbind(kras.cor.matrix, cors)
    }

    out.png <- paste("output/", analysis.name, "-kras-cms-meta-gene", ".png", sep="")
    png(out.png)
    ## Create the plot
    ## p <- ggplot(data=cor.matrix, aes(x=KRAS, y=expr, colour=gene))
    p <- ggplot(data=cor.matrix, aes(x=KRAS, y=expr))
    
    ## Create a box plot where the x axis is KRAS mutation status ...
    p <- p + geom_boxplot(aes(fill=KRAS))
    ## p <- p + geom_jitter()
    p <- p + geom_beeswarm()
    p <- p + ylab("Correlation with metagene")
    
    ## ... faceted on CMS label.
    p <- p + facet_grid(~CMS)
    print(p)
    d <- dev.off()
    
    out.png <- paste("output/", analysis.name, "-kras-meta-gene", ".png", sep="")
    png(out.png)
    ## Create the plot
    ## p <- ggplot(data=kras.cor.matrix, aes(x=KRAS, y=expr, colour=gene))
    p <- ggplot(data=kras.cor.matrix, aes(x=KRAS, y=expr))
    
    ## Create a box plot where the x axis is KRAS mutation status ...
    p <- p + geom_boxplot(aes(fill=KRAS))
    ##    p <- p + geom_jitter()
    p <- p + geom_beeswarm()
    p <- p + ylab("Correlation with metagene")
    
    ## ... faceted on CMS label.
    ## p <- p + facet_grid(~CMS)
    print(p)
    d <- dev.off()
    
}

## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
doAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.nras.mutants <- TRUE
    if(exclude.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        ## flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
        clin.tmp <- clin.tmp[!flag, ]
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
    
    ## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    ## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
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
    kras <- clin.m[, kras.status.field]
    kras.score <- clin.m[, kras.score.field]        
    cms <- as.character(clin.m$cms_label)
    msi <- clin.m$msi
    site <- clin.m$site
    site.factor <- factor(site)
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)
    esn <- tmp
    expr.mask <- expr.m

    if(FALSE) {
      save(clin.m, file="clin.Rdata")
      save(kras, file="kras.Rdata")
      save(esn, file="esn.Rdata")
      save(bindeaGsets, file="bind.Rdata")
      save(expr.mask, file="expr.Rdata")
      save(cms, file="cms.Rdata")
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]

    ## Analysis 1: biomarker ~ kras
    num.biomarkers <- length(to.plot)

    kras.tbl <- data.frame(variable=to.plot, pval=rep(NA, num.biomarkers))
    rownames(kras.tbl) <- to.plot    
    for(sti in 1:num.biomarkers) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", st, " ~ ", "KRAS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=rep("ALL", length(esn.st)))        

        lm.obj <- lm(as.formula(paste("expr", "KRAS", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        ## Compute via wilcox below.
##        kras.flag <- grepl(x=coeffs.2$coefficient, pattern="KRAS") & !grepl(x=coeffs.2$coefficient, pattern=":")
##        kras.tbl$pval[sti] <- coeffs.2[kras.flag, 5]
    }

    ## Adjust the pvals for the kras analysis
##    kras.tbl$padj <- p.adjust(kras.tbl$pval, method="BH")

    pval <- apply(esn[to.plot,,drop=F], 1, function(x) wilcox.test(x ~ kras)$p.value)
    a <- p.adjust(pval,method="BH")
    b <- apply(esn[to.plot,,drop=F], 1, function(x) (mean(x[kras == 1]) - mean(x[kras == 0])) > 0)
    kras.tbl <- data.frame(variable=to.plot, pval=pval, apval=a, mutUp=b)
    rownames(kras.tbl) <- to.plot
    
    ## Create the individual expression plots for each gene/gene set for
    ## the kras analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", st, " ~ ", "KRAS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=rep("ALL", length(esn.st)))

        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], ".png", sep="")        
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## The rest of the ugliness below corresponds to the error bars.
        ## Use raw pvalue for univariate analysis
        ##        pval <- as.numeric(kras.tbl$apval[sti])
        pval <- as.numeric(kras.tbl$pval[sti])
        
        yoffset <- 0.01 * ( max(df$expr) - min(df$expr) )
        ## We only add one error bar.  Shift the ylimit to accomodate this shift.
        mask <- df$CMS == "ALL" & df$KRAS == kras.states[1]        
        mt.max <- max(df$expr[mask])

        mask <- df$CMS == "ALL" & df$KRAS == kras.states[2]        
        wt.max <- max(df$expr[mask])

        mt.wt.max <- max(mt.max, wt.max)
        
        ## Add the error bars between KRAS MT and KRAS WT expression values (within a facet)
        path.df <- data.frame(x = c(1, 1, 2, 2), y=c(mt.max + yoffset, mt.wt.max + 2 * yoffset, mt.wt.max + 2 * yoffset, wt.max + yoffset))
        p <- p + geom_path(data=path.df, aes(x, y))
        p <- p + annotate("text", x = 1.5, y = mt.wt.max + 4 * yoffset, label = pval.to.text(pval))
        print(p)
        d <- dev.off()
    })

    write.table(file=paste0("output/", analysis.name, "-kras-tbl.xls"), kras.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    ##    cms2.last.levels <- c(cmses[cmses != "CMS2"], "CMS2")
    cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    cms2.last.levels <- sort(cms.levels, decreasing=TRUE)
    if("NOLBL" %in% unique.cms) {
        cms2.last.levels <- c(cms2.last.levels, "NOLBL")
    }
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    ## Analysis 2: biomarker ~ cms
    cms.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {
        if(!(prefix %in% cms.levels)) { next }
        pval.col <- paste0(prefix, ".pval")
        cms.tbl[,pval.col] <- rep(NA, nrow(cms.tbl))
    }
    rownames(cms.tbl) <- to.plot
    ## Do the same analysis, but coarse grained.
    cms.coarse.tbl <- cms.tbl

    for(sti in 1:num.biomarkers) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", st, " ~ ", "CMS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=cms2.first.levels))
        df <- df[df$CMS != "NOLBL",]                        
        
        lm.obj <- lm(as.formula(paste("expr", "CMS", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-cms", "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cms2.first.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }            
            cms.str <- prefix
            pval.col <- paste0(cms.str, ".pval")
            pval.flag <- grepl(x=coeffs.2$coefficient, cms.str)
            ## Set this below using wilcox
##            cms.tbl[sti,pval.col] <- coeffs.2[pval.flag, 5]
        }

        ## Do a coarse grained analysis
        cms.coarse <- cms
        flag <- cms %in% c("CMS1", "CMS3", "CMS4")
        cms.coarse[flag] <- "CMSOther"
        flag <- cms %in% c("CMS2")
        cms.coarse[flag] <- "CMS2"
        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms.coarse, levels=c("CMS2", "CMSOther", "NOLBL")))
        df <- df[df$CMS != "NOLBL",]                        

        lm.obj <- lm(as.formula(paste("expr", "CMS", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-cms", "-sum-coarse.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum-coarse.tsv", replacement="-diag-coarse.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        cms.str <- "CMSOther"
        pval.col <- paste0(cms.str, ".pval")
        pval.flag <- grepl(x=coeffs.2$coefficient, cms.str)
            ## Set this below using wilcox        
##        cms.coarse.tbl[sti,pval.col] <- coeffs.2[pval.flag, 5]
        
    }

    ## Adjust the pvals for the cms analysis
##    for(cmsindx in c(1, 3, 4)) {
##        cms.str <- paste0("CMS", cmsindx)
##        pval.col <- paste0(cms.str, ".pval")
##        padj.col <- paste0(cms.str, ".apval")        
##        cms.tbl[,padj.col] <- p.adjust(cms.tbl[,pval.col], method="BH")
##    }
##    pval.col <- paste0("CMSOther", ".pval")
##    padj.col <- paste0("CMSOther", ".apval")        
##    cms.coarse.tbl[,padj.col] <- p.adjust(cms.coarse.tbl[,pval.col], method="BH")        

    cms.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {    
        if(!(prefix %in% cms.levels)) { next }        
        cms.str <- prefix
        flag <- ( cms == "CMS2") | ( cms == cms.str )
        cms.flag <- cms[flag]
        pval <- apply(esn[to.plot,flag,drop=F], 1, function(x) wilcox.test(x ~ cms.flag)$p.value)
        a <- p.adjust(pval,method="BH")
        b <- apply(esn[to.plot,flag,drop=F], 1, function(x) (mean(x[cms.flag == "CMS2"]) - mean(x[cms.flag == cms.str])) > 0)
        c.tmp <- cbind(pval=pval, apval=a, cms2Up=b)
        colnames(c.tmp) <- paste0(cms.str, ".", colnames(c.tmp))
        cms.tbl <- cbind(cms.tbl, c.tmp)
    }
    rownames(cms.tbl) <- to.plot
    
    ## Do the same analysis, but coarse grained.
    cms.coarse <- cms
    flag <- cms %in% c("CMS1", "CMS3", "CMS4")
    cms.coarse[flag] <- "CMSOther"
    flag <- cms %in% c("CMS2")
    cms.coarse[flag] <- "CMS2"
    flag <- ( cms.coarse == "CMS2" ) | ( cms.coarse == "CMSOther" )
    cms.coarse.flag <- cms.coarse[flag]
    pval <- apply(esn[to.plot,flag,drop=F], 1, function(x) wilcox.test(x ~ cms.coarse.flag)$p.value)
    a <- p.adjust(pval,method="BH")
    b <- apply(esn[to.plot,flag,drop=F], 1, function(x) (mean(x[cms.coarse.flag == "CMS2"]) - mean(x[cms.coarse.flag == "CMSOther"])) > 0)
    
    cms.coarse.tbl <- data.frame(variable=to.plot, CMSOther.pval=pval, CMSOther.apval=a, CMSOther.cms2Up=b)
    rownames(cms.coarse.tbl) <- to.plot

    ## Create the individual expression plots for each gene/gene set for
    ## the cms analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", st, " ~ ", "kras", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        df <- df[df$CMS != "NOLBL",]                        

        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-cms", ".png", sep="")
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## The rest of the ugliness below corresponds to the error bars.
        ## Let's include the unadjusted pvalues.
        yoffset <- 0.01 * ( max(df$expr) - min(df$expr) )

        ## panel.maxs will track the maximum value currently plotted in each facet.
        ## This will allow us to know where to position the error bar in each facet.
        panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl
            m <- max(df$expr[mask])
            m
        })
        names(panel.maxs) <- cms.levels
        
        ## Draw the error bars from CMS2 to CMS1, 3, and 4.
        for(prefix in cms2.first.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }            
            cms.str <- prefix
            ## Use raw pvalue for univariate analysis
            ## pval.col <- paste0(cms.str, ".apval")
            pval.col <- paste0(cms.str, ".pval")
            pval <- as.numeric(cms.tbl[sti, pval.col])
        
            cms2.max <- panel.maxs["CMS2"]
            cmsi.max <- panel.maxs[cms.str]
            
            cms2indx <- which(sort(cms2.first.levels) == "CMS2")[1]
            cmsindx <- which(sort(cms2.first.levels) == prefix)[1]
            
            both.max <- max(cms2.max, cmsi.max)
            xs <- c(cmsindx, cmsindx, cms2indx, cms2indx)
            ys <- c(cmsi.max + yoffset, both.max + 2 * yoffset, both.max + 2 * yoffset, cms2.max + yoffset)
            if(cmsindx > cms2indx) {
                xs <- c(cms2indx, cms2indx, cmsindx, cmsindx)
                ys <- c(cms2.max + yoffset, both.max + 2 * yoffset, both.max + 2 * yoffset, cmsi.max + yoffset)
            }

            panel.maxs["CMS2"] <- both.max + 4 * yoffset
            panel.maxs[cms.str] <- both.max + 4 * yoffset            
            
            path.df <- data.frame(x = xs, y = ys)
            p <- p + geom_path(data=path.df, aes(x, y))
            p <- p + annotate("text", x = 0.5 * (cms2indx + cmsindx), y = both.max + 4 * yoffset, label = pval.to.text(pval))
        }
        print(p)
        d <- dev.off()

        ## Plot the coarse grained tables
        cms.coarse <- cms
        flag <- cms %in% c("CMS1", "CMS3", "CMS4")
        cms.coarse[flag] <- "CMSOther"
        flag <- cms %in% c("CMS2")
        cms.coarse[flag] <- "CMS2"
        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms.coarse, levels=c("CMS2", "CMSOther", "NOLBL")))
        df <- df[df$CMS != "NOLBL",]                        

        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-cms-coarse", ".png", sep="")
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## The rest of the ugliness below corresponds to the error bars.
        ## Let's include the unadjusted pvalues.
        yoffset <- 0.01 * ( max(df$expr) - min(df$expr) )

        ## panel.maxs will track the maximum value currently plotted in each facet.
        ## This will allow us to know where to position the error bar in each facet.
        panel.maxs <- sapply(c("CMS2", "CMSOther"), function(cmsLbl) {
            mask <- df$CMS == cmsLbl
            m <- max(df$expr[mask])
            m
        })
        names(panel.maxs) <- c("CMS2", "CMSOther")
        
        ## Draw the error bars from CMS2 to CMS1, 3, and 4.
        cms.str <- "CMSOther"
        ## Use raw pvalue for univariate analysis        
        ## pval.col <- paste0(cms.str, ".apval")
        pval.col <- paste0(cms.str, ".pval")
        pval <- as.numeric(cms.coarse.tbl[sti, pval.col])
        
        cms2.max <- panel.maxs["CMS2"]
        cmsi.max <- panel.maxs[cms.str]

        both.max <- max(cms2.max, cmsi.max)
        xs <- c(1, 1, 2, 2)
        ys <- c(cms2.max + yoffset, both.max + 2 * yoffset, both.max + 2 * yoffset, cmsi.max + yoffset)

        path.df <- data.frame(x = xs, y = ys)
        p <- p + geom_path(data=path.df, aes(x, y))
        p <- p + annotate("text", x = 0.5 * (1 + 2), y = both.max + 4 * yoffset, label = pval.to.text(pval))
        print(p)
        d <- dev.off()
        
    })

    write.table(file=paste0("output/", analysis.name, "-cms-tbl.xls"), cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-cms-tbl-coarse.xls"), cms.coarse.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    ## End analysis 2

    ## Analysis 3: biomarker ~ cms1-kraswt + cms1-krasmt + cms2-kraswt + ... (with ref = cms2:krasMT)
    kras.cms.tbl <- data.frame(variable=to.plot)
    rownames(kras.cms.tbl) <- to.plot        
    base <- ".vs.CMS2MT"
    for(prefix in cms2.first.levels) {
        if(!(prefix %in% cms.levels)) { next }        
        for(genotype in kras.states) {
            if( ( prefix != "CMS2" ) || ( genotype != "MT" ) ) {
                pval.col <- paste0(prefix, genotype, base, ".pval")
                kras.cms.tbl[,pval.col] <- rep(NA, nrow(kras.cms.tbl))
            }
        }
    }

    for(sti in 1:num.biomarkers) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]
        
        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=cms2.first.levels))
        df <- df[df$CMS != "NOLBL",]                        

        df$cms.kras <- apply(cbind(as.character(df$CMS), as.character(df$KRAS)), 1, function(row) paste0(row[1], row[2]))
        df$cms.kras <- factor(df$cms.kras)
        df$cms.kras <- relevel(df$cms.kras, ref = "CMS2MT")
        
        cat(paste0("Computing ", st, " ~ ", "cms.kras", "\n"))
        
        lm.obj <- lm(as.formula(paste("expr", "cms.kras", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-cms-kras", "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()

        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        base <- ".vs.CMS2MT"
        for(prefix in cms2.first.levels) {        
            if(!(prefix %in% cms.levels)) { next }            
            for(genotype in kras.states) {
                if( ( prefix != "CMS2" ) || ( genotype != "MT" ) ) {
                    pval.col <- paste0(prefix, genotype, base, ".pval")
                    cms.kras.str <- paste0("cms.kras", prefix, genotype)
                    pval.flag <- grepl(x=coeffs.2$coefficient, cms.kras.str)
                    ## working
                    ## todo to do
                    ## replace this with wilcox
##                    kras.cms.tbl[sti,pval.col] <- coeffs.2[pval.flag, 5]
                }
            }
        }
    }

    ## Adjust the pvals
    base <- ".vs.CMS2MT"
    for(prefix in cms2.first.levels) {            
        if(!(prefix %in% cms.levels)) { next }        
        for(genotype in kras.states) {
            if( ( prefix != "CMS2" ) || ( genotype != "MT" ) ) {
                pval.col <- paste0(prefix, genotype, base, ".pval")
                flag1 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == "CMS2" ) & ( kras.label == "MT" ) ) )
                flag2 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == prefix ) & ( kras.label == genotype ) ) )
                pvals <- rep(1, length(to.plot))
                if((length(which(flag1)) >= 2) && (length(which(flag2)) >= 2)) {
                    pvals <- apply(esn[to.plot,,drop=F], 1, function(x) {
                        wilcox.test(x = x[flag1], y = x[flag2])$p.value
                    })
                }
                kras.cms.tbl[,pval.col] <- pvals
                padj.col <- paste0(prefix, genotype, base, ".apval")
                kras.cms.tbl[,padj.col] <- p.adjust(kras.cms.tbl[,pval.col], method="BH")
            }
        }
    }
    
    ## Create the individual expression plots for each gene/gene set for
    ## the cms analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Plotting ", st, " ~ ", "cms.kras", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        df <- df[df$CMS != "NOLBL",]                        

        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-kras-cms", ".png", sep="")
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits) - min(ylimits) )
        ## And we will add 7 error bars.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 7        
        ylimits[2] <- ylimits[2] + 6 * num.error.bars * yoffset
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
        mt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[1]
            m <- max(df$expr[mask])
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$expr[mask])
            m
        })
        names(wt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[2]))
        
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)

        ## Add the error bars between KRAS MT and KRAS WT expression values (within a facet)
        for(prefix in sort(cms2.last.levels)) {
            if(!(prefix %in% cms.levels)) { next }            
            for(genotype in kras.states) {
                if( ( prefix != "CMS2" ) || ( genotype != "MT" ) ) {
                    cmsLbl2 <- "CMS2"
                    kras.state2 <- "MT"
                    cmsIndx2 <- which(sort(cms2.last.levels) == cmsLbl2)[1]
                    cmsLbl1 <- prefix
                    kras.state1 <- genotype
                    cmsIndx1 <- which(sort(cms2.last.levels) == cmsLbl1)[1]                    
                    base <- ".vs.CMS2MT"
                    cmpLbl <- paste0(prefix, genotype, base)
                    tmp <- draw.err.bar(st, kras.cms.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states)
                    g <- tmp[["g"]]
                    panel.maxs <- tmp[["panel.maxs"]]
                }
            }
        }
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()
        
    })

    write.table(file=paste0("output/", analysis.name, "-kras-cms-tbl.xls"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    ## test.harm.mut.R <- doAnalysis(expr.test.harm, clin.test.harm, paste0(test.set, "-harm-mt"), to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))    
    ## End Analysis 3

    ## Analysis 4: biomarker ~ cms + kras + cms:kras (with ref = cms1, so cms2:krasmt shows up)
    interaction.tbl <- data.frame(variable=to.plot)
    rownames(interaction.tbl) <- to.plot
    cols <- c(kras.col)
    for(prefix in cms2.last.levels[-1]) {
        if(prefix == "NOLBL") { next }
        col <- paste0(tolower(prefix), ".wt.vs.", tolower(cms2.last.levels[1]), ".wt")
        cols <- c(cols, col)
        col <- paste0(tolower(prefix), ".mt.vs.", tolower(prefix), ".wt")
        cols <- c(cols, col)        
    }
    
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        interaction.tbl[,col] <- rep(NA, nrow(interaction.tbl))
    }

    transformed.interaction.tbl <- interaction.tbl

    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model with interaction for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        val <- esn.st
        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")))

        formula <- "KRAS + CMS + CMS:KRAS"
        formula.name<- "kras-cms-interaction"
        ## Now exclude NOLBL
        df.lm <- df[df$CMS != "NOLBL",]
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        ## library(lmPerm)
        ## lm.obj <- lmp(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        ## lm.obj <- glm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm, family=gaussian(link="log"))
        lm.sum <- summary(lm.obj)
        print(lm.sum)
        capture.output(lm.sum, file = sum.file)

        bf.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-bf.tsv")
        bf <- generalTestBF(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
##        save(bf, file="bf.Rd")
        capture.output(bf/bf["KRAS + CMS"], file=bf.file)
        
        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        qqnorm(df.lm$expr)
        qqline(df.lm$expr)
        plot(lm.obj)
        d <- dev.off()
        lm.sum <- summary(lm.obj)
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern="KRASMT") & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
        interaction.tbl[sti, col] <- coeffs.2[flag, 5]

        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }
            col <- paste0(tolower(prefix), ".wt.vs.", tolower(cms2.last.levels[1]), ".wt.pval")
            cmsstr <- prefix
            flag <- grepl(x=coeffs.2$coefficient, pattern=cmsstr) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(!any(flag)) {
                cat(paste0("Could not find ", cmsstr, "\n"))
                print(coeffs.2)
                next
            }
            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }
        
        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }
            col <- paste0(tolower(prefix), ".mt.vs.", tolower(prefix), ".wt.pval")
            cmsstr <- prefix
            flag <- grepl(x=coeffs.2$coefficient, pattern=cmsstr) & grepl(x=coeffs.2$coefficient, pattern=":")
            if(!any(flag)) {
                cat(paste0("Could not find ", cmsstr, "\n"))
                print(coeffs.2)
                next
            }
            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        ## Repeat the analysis for transformed data
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.sum <- summary(lm.obj)
        print(lm.sum)
        capture.output(lm.sum, file = sum.file)

        bf.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-bf.tsv")
        bf <- generalTestBF(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        capture.output(bf/bf["KRAS + CMS"], file=bf.file)        
        
        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        qqnorm(df.lm$expr)
        qqline(df.lm$expr)
        plot(lm.obj)
        d <- dev.off()
        lm.sum <- summary(lm.obj)
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern="KRASMT") & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
        transformed.interaction.tbl[sti, col] <- coeffs.2[flag, 5]

        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }
            col <- paste0(tolower(prefix), ".wt.vs.", tolower(cms2.last.levels[1]), ".wt.pval")
            cmsstr <- prefix
            flag <- grepl(x=coeffs.2$coefficient, pattern=cmsstr) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(!any(flag)) {
                cat(paste0("Could not find ", cmsstr, "\n"))
                print(coeffs.2)
                next
            }            
            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            transformed.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }
        
        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }
            col <- paste0(tolower(prefix), ".mt.vs.", tolower(prefix), ".wt.pval")
            cmsstr <- prefix
            flag <- grepl(x=coeffs.2$coefficient, pattern=cmsstr) & grepl(x=coeffs.2$coefficient, pattern=":")
            if(!any(flag)) {
                cat(paste0("Could not find ", cmsstr, "\n"))
                print(coeffs.2)
                next
            }            
            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            transformed.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
        cat(paste0("Attempting to padj to col ", apcol, "\n"))
        interaction.tbl[,apcol] <- p.adjust(interaction.tbl[,pcol], method="BH")
        transformed.interaction.tbl[,apcol] <- p.adjust(transformed.interaction.tbl[,pcol], method="BH")        
    }

    ## Create the individual expression plots for each gene/gene set for
    ## the cms analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing transformed ", st, " ~ ", "cms.kras", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        df <- df[df$CMS != "NOLBL",]
        df.transformed <- df
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-kras-cms-interaction", ".png", sep="")
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits) - min(ylimits) )
        ## And we will add 1 error bars.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 1       
        ylimits[2] <- ylimits[2] + 6 * num.error.bars * yoffset
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
        mt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[1]
            m <- max(df$expr[mask])
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$expr[mask])
            m
        })
        names(wt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[2]))
        
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)

        ## Add the error bars between MT and WT expression values (within a facet)
        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }
            cmsLbl <- prefix
            cmsIndx1 <- which(sort(cms2.last.levels) == prefix)
            cmsIndx2 <- which(sort(cms2.last.levels) == prefix)
            cmsLbl1 <- cmsLbl
            cmsLbl2 <- cmsLbl
            kras.state1 <- kras.states[1]
            kras.state2 <- kras.states[2]
            cmpLbl <- tolower(paste0(cmsLbl, ".mt.vs.", cmsLbl, ".wt"))
            tmp <- draw.err.bar(st, interaction.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states)
            g <- tmp[["g"]]
            panel.maxs <- tmp[["panel.maxs"]]
        }
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()

        ## Repeat for transformed data
        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-kras-cms-interaction-transformed", ".png", sep="")
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df.transformed, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits) - min(ylimits) )
        ## And we will add 1 error bars.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 1       
        ylimits[2] <- ylimits[2] + 6 * num.error.bars * yoffset
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
        mt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df.transformed$CMS == cmsLbl & df.transformed$KRAS == kras.states[1]
            m <- max(df.transformed$expr[mask])
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df.transformed$CMS == cmsLbl & df.transformed$KRAS == kras.states[2]
            m <- max(df.transformed$expr[mask])
            m
        })
        names(wt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[2]))
        
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)

        ## Add the error bars between MT and WT expression values (within a facet)
        for(prefix in sort(cms2.last.levels[-1])) {
            if(!(prefix %in% cms.levels)) { next }
            cmsLbl <- prefix
            cmsIndx1 <- which(sort(cms2.last.levels) == prefix)
            cmsIndx2 <- which(sort(cms2.last.levels) == prefix)
            cmsLbl1 <- cmsLbl
            cmsLbl2 <- cmsLbl
            kras.state1 <- kras.states[1]
            kras.state2 <- kras.states[2]
            cmpLbl <- tolower(paste0(cmsLbl, ".mt.vs.", cmsLbl, ".wt"))

            tmp <- draw.err.bar(st, transformed.interaction.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states)
            g <- tmp[["g"]]
            panel.maxs <- tmp[["panel.maxs"]]
        }
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()
        
    })

    write.table(file=paste0("output/", analysis.name, "-kras-cms-interaction-tbl.xls"), interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-kras-cms-transformed-interaction-tbl.xls"), transformed.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    
    ## End Analysis 4

    ## Analysis 5: biomarker ~ cms + kras (with ref = cms1, so cms2:krasmt shows up)
    no.interaction.tbl <- data.frame(variable=to.plot)
    rownames(no.interaction.tbl) <- to.plot
    cols <- c(kras.col, cms2.last.levels)
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        no.interaction.tbl[,col] <- rep(NA, nrow(no.interaction.tbl))
    } 

    transformed.no.interaction.tbl <- no.interaction.tbl
    
    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")))

        formula <- "KRAS + CMS"
        formula.name <- "kras-cms-no-interaction"
        ## Now exclude NOLBL
        df.lm <- df[df$CMS != "NOLBL",]
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern="KRASMT") & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
        no.interaction.tbl[sti, col] <- coeffs.2[flag, 5]

        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }            
            cmsstr <- prefix
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=cmsstr) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(!any(flag)) {
                cat(paste0("Could not find ", cmsstr, "\n"))
                print(coeffs.2)
                next
            }
            
            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            no.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        ## Repeat above for transformed data
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern="KRASMT") & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
        transformed.no.interaction.tbl[sti, col] <- coeffs.2[flag, 5]

        for(prefix in cms2.last.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }
            cmsstr <- prefix
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=cmsstr) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(!any(flag)) {
                cat(paste0("Could not find ", cmsstr, "\n"))
                print(coeffs.2)
                next
            }
            
            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            transformed.no.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
        cat(paste0("Attempting to padj to col ", apcol, "\n"))
        transformed.no.interaction.tbl[,apcol] <- p.adjust(transformed.no.interaction.tbl[,pcol], method="BH")
    }

    write.table(file=paste0("output/", analysis.name, "-kras-cms-no-interaction-tbl.xls"), no.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-kras-cms-transformed-no-interaction-tbl.xls"), transformed.no.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    
    ## End Analysis 5
    
    ## Analysis 6: biomarker ~ cms + kras (with error bars from kras-mt vs kras-wt within cms)

    mt.vs.wt.tbl <- do.call("cbind",lapply(cms.levels, function(cmsLbl){
        cat(paste("Computing ", kras.states[1], " vs ", kras.states[2], " for ", cmsLbl, "\n", sep=""))
        mask <- cmsLbl ==  cms
        print(length(which(mask)))
        pval <- rep(1, length(to.plot))
        if((length(which(mask)) >= 2) && (length(unique(kras[mask])) == 2)) {
            pval <- apply(esn[to.plot,mask,drop=F], 1, function(x) wilcox.test(x ~ kras[mask])$p.value)
        }
        a <- p.adjust(pval,method="BH")
        b <- apply(esn[to.plot,mask,drop=F], 1, function(x) (mean(x[kras[mask] == 1]) - mean(x[kras[mask] == 0])) > 0)
        cbind(pval=pval,apval=a, mutUp=b)
    }))

    colnames(mt.vs.wt.tbl) <- do.call(c, lapply(paste0(cms.levels, ".MT.vs.WT"), function(x) paste(x, c(".pval", ".apval", ".mutUp"), sep="")))
    
    ## Plot the mt-vs-wt wilcox error bars
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Plotting mt vs wt within CMS for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, unique.cms))
        df <- df[df$CMS != "NOLBL",]                        

        ## Create a PDF for this gene/gene set (with facets)
        out.png <- paste("output/", st, "-", analysis.name, "-kras-cms-faceted", ".png", sep="")
        png(out.png)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits) - min(ylimits) )
        ## And we will only add 1 error bar.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 1
        ylimits[2] <- ylimits[2] + 6 * num.error.bars * yoffset
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
        mt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[1]
            m <- max(df$expr[mask])
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$expr[mask])
            m
        })
        names(wt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[2]))
        
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)
        ## Add the error bars between KRAS MT and KRAS WT expression values (within a facet)
        for(prefix in cms.levels) {
            if(!(prefix %in% cms.levels)) { next }            
            cmsLbl <- prefix
            cmsIndx1 <- which(sort(cms.levels) == prefix)
            cmsIndx2 <- cmsIndx1
            cmsLbl1 <- cmsLbl
            cmsLbl2 <- cmsLbl1
            kras.state1 <- kras.states[1]
            kras.state2 <- kras.states[2]
            cmpLbl <- paste0(cmsLbl, ".MT.vs.WT")

            tmp <- draw.err.bar(st, mt.vs.wt.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states)
            g <- tmp[["g"]]
            panel.maxs <- tmp[["panel.maxs"]]
        }
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()
        
    })

    write.table(file=paste0("output/", analysis.name, "-kras-mt-vs-wt-tbl.xls"), mt.vs.wt.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    ## test.harm.mut.R <- doAnalysis(expr.test.harm, clin.test.harm, paste0(test.set, "-harm-mt"), to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT", "WT"))    
    ## End Analysis 6

    ## Analysis 7: biomarker ~ cms + krs + mss + site
    full.model.tbl <- data.frame(variable=to.plot)
    rownames(full.model.tbl) <- to.plot
    cols <- c(kras.col, cms2.last.levels.no.lbl[-1], site.levels[-1])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        cols <- c(cols, msi.levels[-1])
    }
    
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        full.model.tbl[,col] <- rep(NA, nrow(full.model.tbl))
    } 

    transformed.full.model.tbl <- full.model.tbl
    
    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor)
        formula <- "KRAS + CMS + site"
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            formula <- paste0(formula, " + msi")
            df$msi <- msi.factor
        }
        flag <- unlist(apply(df, 1, function(row) any(is.na(row))))
        df <- df[!flag,]

        formula.name <- "full-model"
        ## Now exclude NOLBL
        df.lm <- df[df$CMS != "NOLBL",]
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.png")
        png(forest.file)
        p <- forest_model(lm.obj)
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                print(col)
                print(coeffs.2)
                stop("Could not find coefficient\n")
            }
            full.model.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        ## Repeat above for transformed data
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.png")
        png(forest.file)
        p <- forest_model(lm.obj)
        print(p)
        d <- dev.off()

        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                print(col)
                print(coeffs.2)
                stop("Could not find coefficient\n")
            }
            transformed.full.model.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
        cat(paste0("Attempting to padj to col ", apcol, "\n"))
        transformed.full.model.tbl[,apcol] <- p.adjust(transformed.full.model.tbl[,pcol], method="BH")
        full.model.tbl[,apcol] <- p.adjust(full.model.tbl[,pcol], method="BH")        
    }

    write.table(file=paste0("output/", analysis.name, "-full-model-tbl.xls"), full.model.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-transformed-full-model-tbl.xls"), transformed.full.model.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    
    
    ## End Analysis 7

    ## Analysis 8 -- univariate msi and site
    formulas <- c("site")
    col.list <- list()
    base.list <- list()
    col.list[[1]] <- site.levels[-1]
    base.list[[1]] <- site.levels[1]
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {    
        formulas <- c(formulas, "msi")
        col.list[[2]] <- msi.levels[-1]
        base.list[[2]] <- msi.levels[1]
    }
    formula.tbls <- list()
    for(i in 1:length(formulas)) {
        formula <- formulas[i]
        formula.name <- formula
        
        formula.tbls[[i]] <- data.frame(variable=to.plot)
        rownames(formula.tbls[[i]]) <- to.plot
        cols <- col.list[[i]]
        base <- base.list[[i]]
        
        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            formula.tbls[[i]][,col] <- rep(NA, nrow(formula.tbls[[i]]))
        } 

        df <- data.frame(CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor)
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$msi <- msi.factor
        }

        for(lvl in cols) {
            cat(paste("Computing ", lvl, " vs ", base, "\n", sep=""))
            vec <- as.vector(df[,formula,drop=T])
            flag1 <- !is.na(vec) & ( vec == base)
            flag2 <- !is.na(vec) & ( vec == lvl)            
            pvals <- apply(esn[to.plot,,drop=F], 1, function(x) wilcox.test(x = x[flag1], y = x[flag2])$p.value)
            
            pval.col <- paste0(lvl, ".pval")
            padj.col <- paste0(lvl, ".apval")
            
            formula.tbls[[i]][,pval.col] <- pvals
            formula.tbls[[i]][,padj.col] <- p.adjust(formula.tbls[[i]][,pval.col], method="BH")
        }
            
        write.table(file=paste0("output/", analysis.name, "-", formula.name, "-tbl.xls"), formula.tbls[[i]], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    

        ## Create the individual expression plots for each gene/gene set for
        ## the cms analysis
        lapply(1:length(to.plot), function(sti){

            ## st will be the gene/gene set to plot
            st <- to.plot[sti]

            ## The y axis label for this gene/gene set
            ylab <- ylabels[sti]
            
            cat(paste("Plotting ", st, " ~ ", formula, "\n", sep=""))

            ## Subset the data to the gene/gene set of interest
            esn.st <- as.vector(esn[st,])

            ## Get the KRAS mutation status of each sample
            kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
            
            df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor)
            if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {            
                df$msi <- msi.factor
            }
            flag <- unlist(apply(df[,formula,drop=F], 1, function(row) any(is.na(row))))
            df <- df[!flag,]

            ## Create a png for this gene/gene set (with facets)
            out.png <- paste("output/", st, "-", analysis.name, "-", formula, ".png", sep="")
            png(out.png)

            ## Create the plot
            formula.indx <- which(colnames(df) == formula)[1]
            expr.indx <- which(colnames(df) == "expr")[1]
            p <- ggplot(data=df, aes_string(x=colnames(df)[formula.indx], y=colnames(df)[expr.indx]))
            p <- p + ylab(ylab)

            ## Create a box plot where the x axis is CMS ...
            p <- p + geom_boxplot(aes_string(fill=colnames(df)[formula.indx]))
##            p <- p + geom_jitter()
            p <- p + geom_beeswarm()
            
            ## The rest of the ugliness below corresponds to the error bars.
            ## Let's include the unadjusted pvalues.
            yoffset <- 0.01 * ( max(df$expr) - min(df$expr) )

            ## panel.maxs will track the maximum value currently plotted in each facet.
            ## This will allow us to know where to position the error bar in each facet.
            lvls <- levels(df[,formula,drop=T])
            panel.maxs <- sapply(lvls, function(lvl) {
                mask <- df[,formula,drop=T] == lvl
                m <- max(df$expr[mask])
                m
            })
            names(panel.maxs) <- lvls
        
            ## Draw the error bars from the first level to each of the others
            for(prefix in cols) {
                str <- prefix
                ## Use raw pvalue for univariate analysis                        
                ## pval.col <- paste0(str, ".apval")
                pval.col <- paste0(str, ".pval")
                pval <- as.numeric(formula.tbls[[i]][sti, pval.col])
        
                base.max <- panel.maxs[base]
                col.max <- panel.maxs[str]
            
                base.indx <- which(lvls == base)[1]
                col.indx <- which(lvls == prefix)[1]
            
                both.max <- max(base.max, col.max)
                xs <- c(col.indx, col.indx, base.indx, base.indx)
                ys <- c(col.max + yoffset, both.max + 2 * yoffset, both.max + 2 * yoffset, base.max + yoffset)
                if(col.indx > base.indx) {
                    xs <- c(base.indx, base.indx, col.indx, col.indx)
                    ys <- c(base.max + yoffset, both.max + 2 * yoffset, both.max + 2 * yoffset, col.max + yoffset)
                }

                panel.maxs[base] <- both.max + 4 * yoffset
                panel.maxs[str] <- both.max + 4 * yoffset            
                
                path.df <- data.frame(x = xs, y = ys)
                p <- p + geom_path(data=path.df, aes(x, y))
                p <- p + annotate("text", x = 0.5 * (base.indx + col.indx), y = both.max + 4 * yoffset, label = pval.to.text(pval))
            }
            print(p)
            d <- dev.off()
            
        })
        
    }

    ## End Analysis 8

    ret.list <- list(kras.tbl=kras.tbl, cms.tbl=cms.tbl, kras.cms.tbl=kras.cms.tbl, kras.cms.interaction.tbl=interaction.tbl, kras.cms.no.interaction.tbl=no.interaction.tbl, kras.mt.vs.wt.tbl=mt.vs.wt.tbl, site.tbl=formula.tbls[[1]])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        ret.list[["msi.tbl"]] <- formula.tbls[[2]]
    }
    return(ret.list)

    
}

## Draw an error bar between two facets (that may or may not be the same).
draw.err.bar <- function(st, tbl2, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states) {

    ## Get the corresponding adjusted pvalue comparing expression of st (a gene or gene set)
    ## across the two conditions specified by cmsLbl1 and cmsLbl2
    ## Plot raw pvalues!!!
    ## adj.pval.col <- paste(toString(cmpLbl), ".apval", sep="")
    adj.pval.col <- paste(toString(cmpLbl), ".pval", sep="")
    if(!(adj.pval.col %in% colnames(tbl2))) {
        cat(paste0("Could not find ", adj.pval.col, " in tbl2\n"))
        print(colnames(tbl2))
        stop("stop")
    }
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

    ## This is the vertical line from the top of the KRAS mutant or WT expression
    ## data, which are x offset 1 or 2 within the CMS/facet, to what will be the horizontal line of the error bar.
    x1.index <- which(kras.state1 == kras.states)
    x2.index <- which(kras.state2 == kras.states)    
    start <- c(data2npc(x1.index,ranges[[cmsIndx1]][["x.range"]]),
               data2npc(m + yoffset, ranges[[cmsIndx1]][["y.range"]]))
    
    end <- c(data2npc(x1.index,ranges[[cmsIndx1]][["x.range"]]),
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
    start <- c(data2npc(x1.index,ranges[[cmsIndx1]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    
    end <- c(data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[cmsIndx2]][["y.range"]]))

    ## Finally, draw the vertical line on the "right side" of the error bar--
    ## this goes from the error bar to the maximum of the MT KRAS
    ## expression values in the second facet/CMS 2.
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)
    
    start <- c(data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    end <- c(data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]),
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
    ## I believe this 3 comes from the fact that the high range is 2.6 and I was
    ## padding for the space between the facets
    xrange[2] <- xrange[2] + 3 * abs(cmsIndx2 - cmsIndx1)
    ## pos <- c(0.5 * ( data2npc(x1.index,ranges[[cmsIndx1]][["x.range"]]) + data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]) ), data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    pos <- c(data2npc(0.5 * ( x1.index + 3 * abs(cmsIndx2 - cmsIndx1) + x2.index), xrange),
             data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    ## pos <- c(data2npc(0.5, xrange), data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    

    delta <- runif(n=1, min=10^-5, max=10^-3)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2]), t = t + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
}

## TODO
## -- univariate: site and msi
## -- read
## -- what does gary actually want to show in terms of comparisons??
## -- what does gary actually want to show in terms of genes/signatures??
## -- write


## -- how many CMS1 are MSS and non-KRAS NA: (tcga: MT=4, WT=4; french: MT=4, WT=16)


doMergedDatasetAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.nras.mutants <- TRUE
    if(exclude.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        ## flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
        clin.tmp <- clin.tmp[!flag, ]
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
    
    ## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    ## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
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
    kras <- clin.m[, kras.status.field]
    kras.score <- clin.m[, kras.score.field]        
    cms <- as.character(clin.m$cms_label)
    msi <- clin.m$msi
    site <- clin.m$site
    datasets <- clin.m$dataset
    site.factor <- factor(site)
    dataset.factor <- factor(datasets)
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)
    esn <- tmp
    expr.mask <- expr.m

    if(FALSE) {
      save(clin.m, file="clin.Rdata")
      save(kras, file="kras.Rdata")
      save(esn, file="esn.Rdata")
      save(bindeaGsets, file="bind.Rdata")
      save(expr.mask, file="expr.Rdata")
      save(cms, file="cms.Rdata")
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]

    ## Analysis 1: biomarker ~ kras
    num.biomarkers <- length(to.plot)

    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    ##    cms2.last.levels <- c(cmses[cmses != "CMS2"], "CMS2")
    cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    cms2.last.levels <- sort(cms.levels, decreasing=TRUE)
    if("NOLBL" %in% unique.cms) {
        cms2.last.levels <- c(cms2.last.levels, "NOLBL")
    }
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    ## Analysis 7: biomarker ~ cms + krs + mss + site + dataset
    full.model.tbl <- data.frame(variable=to.plot)
    rownames(full.model.tbl) <- to.plot
    cols <- c(kras.col, cms2.last.levels.no.lbl[-1], site.levels[-1])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        cols <- c(cols, msi.levels[-1])
    }
    
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        full.model.tbl[,col] <- rep(NA, nrow(full.model.tbl))
    } 

    transformed.full.model.tbl <- full.model.tbl
    
    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor, dataset=dataset.factor)
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
        capture.output(lm.sum, file = sum.file)

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.png")
        png(forest.file)
        p <- forest_model(lm.obj)
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        png.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-facet.tsv", sep="")        
        png(png.file)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        d <- dev.off()
        
        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                print(col)
                print(coeffs.2)
                stop("Could not find coefficient\n")
            }
            full.model.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        ## Repeat above for transformed data
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.png")
        png(forest.file)
        p <- forest_model(lm.obj)
        print(p)
        d <- dev.off()

        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                print(col)
                print(coeffs.2)
                stop("Could not find coefficient\n")
            }
            transformed.full.model.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
        cat(paste0("Attempting to padj to col ", apcol, "\n"))
        transformed.full.model.tbl[,apcol] <- p.adjust(transformed.full.model.tbl[,pcol], method="BH")
        full.model.tbl[,apcol] <- p.adjust(full.model.tbl[,pcol], method="BH")        
    }

    write.table(file=paste0("output/", analysis.name, "-full-model-tbl.xls"), full.model.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-transformed-full-model-tbl.xls"), transformed.full.model.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    
    
    ## End Analysis 7

    
}

doMergedDatasetAnalysisES <- function(es, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.nras.mutants <- TRUE
    if(exclude.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        ## flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
        clin.tmp <- clin.tmp[!flag, ]
    }
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin.tmp$sample, colnames(es))
    clin.mask <- clin.tmp[!is.na(idxs),]
    expr.mask <- es[, na.omit(idxs)]
    
    ## Exclude any samples that do not have KRAS mutation status annotated
    mask <- !is.na(clin.mask[,kras.status.field])
    clin.m <- clin.mask[mask,]
    expr.m <- expr.mask[, mask]
    
    ## gene set prep
    env <- new.env()
    load("input/markersG.RData",envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
    ## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    ## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
    myc.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "MYC_TARGETS_ZELLER")]]
    bindeaGsets[["myc.sig"]] <- myc.targets
    
    cat(paste("Number of rows analyzed: ", nrow(expr.m), "\n", sep=""))
    cat(paste("Number of columns analyzed: ", ncol(expr.m), "\n", sep=""))

    all.genes.and.sets <- c(rownames(es))
    flag <- to.plot %in% all.genes.and.sets
    to.plot <- to.plot[flag]
    ylabels <- ylabels[flag]
    tmp <- es

    ## Exclude any samples that do not have KRAS mutation status annotated
    kras <- clin.m[, kras.status.field]
    kras.score <- clin.m[, kras.score.field]        
    cms <- as.character(clin.m$cms_label)
    msi <- clin.m$msi
    site <- clin.m$site
    datasets <- clin.m$dataset
    site.factor <- factor(site)
    dataset.factor <- factor(datasets)
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)
    esn <- tmp
    expr.mask <- expr.m

    if(FALSE) {
      save(clin.m, file="clin.Rdata")
      save(kras, file="kras.Rdata")
      save(esn, file="esn.Rdata")
      save(bindeaGsets, file="bind.Rdata")
      save(expr.mask, file="expr.Rdata")
      save(cms, file="cms.Rdata")
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]

    ## Analysis 1: biomarker ~ kras
    num.biomarkers <- length(to.plot)

    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    ##    cms2.last.levels <- c(cmses[cmses != "CMS2"], "CMS2")
    cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    cms2.last.levels <- sort(cms.levels, decreasing=TRUE)
    if("NOLBL" %in% unique.cms) {
        cms2.last.levels <- c(cms2.last.levels, "NOLBL")
    }
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    ## Analysis 7: biomarker ~ cms + krs + mss + site + dataset
    full.model.tbl <- data.frame(variable=to.plot)
    rownames(full.model.tbl) <- to.plot
    cols <- c(kras.col, cms2.last.levels.no.lbl[-1], site.levels[-1])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        cols <- c(cols, msi.levels[-1])
    }
    
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        full.model.tbl[,col] <- rep(NA, nrow(full.model.tbl))
    } 

    transformed.full.model.tbl <- full.model.tbl
    
    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor, dataset=dataset.factor)
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
        capture.output(lm.sum, file = sum.file)

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.png")
        png(forest.file)
        p <- forest_model(lm.obj)
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                print(col)
                print(coeffs.2)
                stop("Could not find coefficient\n")
            }
            full.model.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        png.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-facet.png", sep="")        
        cat(paste0("Creating png ", png.file, "\n"))
        png(png.file)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)
        print(p)
        d <- dev.off()

        
        ## Repeat above for transformed data
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        print(plot(lm.obj))
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.png")
        png(forest.file)
        p <- forest_model(lm.obj)
        print(p)
        d <- dev.off()

        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                print(col)
                print(coeffs.2)
                stop("Could not find coefficient\n")
            }
            transformed.full.model.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
        cat(paste0("Attempting to padj to col ", apcol, "\n"))
        transformed.full.model.tbl[,apcol] <- p.adjust(transformed.full.model.tbl[,pcol], method="BH")
        full.model.tbl[,apcol] <- p.adjust(full.model.tbl[,pcol], method="BH")        
    }

    write.table(file=paste0("output/", analysis.name, "-full-model-tbl.xls"), full.model.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-transformed-full-model-tbl.xls"), transformed.full.model.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    
    
    ## End Analysis 7

    
}

subtractOffDatasetEffect <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.nras.mutants <- TRUE
    if(exclude.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        ## flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
        clin.tmp <- clin.tmp[!flag, ]
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
    
    ## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    ## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
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
    kras <- clin.m[, kras.status.field]
    kras.score <- clin.m[, kras.score.field]        
    cms <- as.character(clin.m$cms_label)
    msi <- clin.m$msi
    site <- clin.m$site
    datasets <- clin.m$dataset
    site.factor <- factor(site)
    dataset.factor <- factor(datasets)
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)
    esn <- tmp
    expr.mask <- expr.m

    if(FALSE) {
      save(clin.m, file="clin.Rdata")
      save(kras, file="kras.Rdata")
      save(esn, file="esn.Rdata")
      save(bindeaGsets, file="bind.Rdata")
      save(expr.mask, file="expr.Rdata")
      save(cms, file="cms.Rdata")
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]

    ## Analysis 1: biomarker ~ kras
    num.biomarkers <- length(to.plot)

    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    ##    cms2.last.levels <- c(cmses[cmses != "CMS2"], "CMS2")
    cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    cms2.last.levels <- sort(cms.levels, decreasing=TRUE)
    if("NOLBL" %in% unique.cms) {
        cms2.last.levels <- c(cms2.last.levels, "NOLBL")
    }
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    ## Analysis 7: biomarker ~ cms + krs + mss + site + dataset
    full.model.tbl <- data.frame(variable=to.plot)
    rownames(full.model.tbl) <- to.plot
    cols <- c(kras.col, cms2.last.levels.no.lbl[-1], site.levels[-1])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        cols <- c(cols, msi.levels[-1])
    }
    
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        full.model.tbl[,col] <- rep(NA, nrow(full.model.tbl))
    } 

    transformed.full.model.tbl <- full.model.tbl

    transformed.esn <- esn
    
    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor, dataset=dataset.factor)
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
        print(lm.sum)

        coeffs <- coefficients(lm.sum)
        for(dataset.coeff in rownames(coeffs)[grepl(x=rownames(coeffs), pattern="dataset")]) {
            dataset.name <- gsub(dataset.coeff, pattern="dataset", replacement="")
            flag <- dataset.factor == dataset.name
            effect <- coeffs[dataset.coeff,1]
            cat(paste0("Transforming data for ", dataset.name, "by subtracting off effect", effect, "\n"))
            print(head(transformed.esn[st,flag]))
            transformed.esn[st,flag] <- transformed.esn[st,flag] - effect
            print(head(transformed.esn[st,flag]))
        }
    }
    transformed.esn
}

subtractOffDatasetEffectWithCombat <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    suppressPackageStartupMessages(library(sva))
    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.nras.mutants <- TRUE
    if(exclude.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        ## flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
        clin.tmp <- clin.tmp[!flag, ]
    }

    na.rows <- unlist(apply(clin.tmp[,c("site","msi","cms_label",kras.status.field)], 1, function(row) any(is.na(row))))
    clin.tmp <- clin.tmp[!na.rows,]
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin.tmp$sample, colnames(expr))
    clin.mask <- clin.tmp[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]
    
    ## Exclude any samples that do not have KRAS mutation status annotated
    mask <- !is.na(clin.mask[,kras.status.field])
    clin.m <- clin.mask[mask,]
    expr.m <- expr.mask[, mask]

    formula <- paste0("~ site + msi + cms_label + ", kras.status.field)
    formula <- paste0("~ 1")
##    formula <- paste0("~ msi")
    combat.model <- model.matrix(as.formula(formula), data=clin.m)
    transformed.expr <- ComBat(dat=expr.m, batch=clin.m$dataset, par.prior=TRUE, mod=combat.model)
    return(transformed.expr)
    
    ## gene set prep
    env <- new.env()
    load("input/markersG.RData",envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
    ## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    ## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
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
    kras <- clin.m[, kras.status.field]
    kras.score <- clin.m[, kras.score.field]        
    cms <- as.character(clin.m$cms_label)
    msi <- clin.m$msi
    site <- clin.m$site
    datasets <- clin.m$dataset
    site.factor <- factor(site)
    dataset.factor <- factor(datasets)
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)
    esn <- tmp
    expr.mask <- expr.m

    if(FALSE) {
      save(clin.m, file="clin.Rdata")
      save(kras, file="kras.Rdata")
      save(esn, file="esn.Rdata")
      save(bindeaGsets, file="bind.Rdata")
      save(expr.mask, file="expr.Rdata")
      save(cms, file="cms.Rdata")
    }
    
    kras.pos <- kras.states[1]
    kras.neg <- kras.states[2]

    ## Analysis 1: biomarker ~ kras
    num.biomarkers <- length(to.plot)

    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    ##    cms2.last.levels <- c(cmses[cmses != "CMS2"], "CMS2")
    cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    cms2.last.levels <- sort(cms.levels, decreasing=TRUE)
    if("NOLBL" %in% unique.cms) {
        cms2.last.levels <- c(cms2.last.levels, "NOLBL")
    }
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    ## Analysis 7: biomarker ~ cms + krs + mss + site + dataset
    full.model.tbl <- data.frame(variable=to.plot)
    rownames(full.model.tbl) <- to.plot
    cols <- c(kras.col, cms2.last.levels.no.lbl[-1], site.levels[-1])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        cols <- c(cols, msi.levels[-1])
    }
    
    for(prefix in cols) {
        col <- paste0(prefix, ".pval")
        full.model.tbl[,col] <- rep(NA, nrow(full.model.tbl))
    } 

    transformed.full.model.tbl <- full.model.tbl

    transformed.esn <- esn
    
    for(sti in 1:length(to.plot)) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing linear model for ", st, "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=c("WT", "MT")), site=site.factor, dataset=dataset.factor)
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
        print(lm.sum)

        coeffs <- coefficients(lm.sum)
        for(dataset.coeff in rownames(coeffs)[grepl(x=rownames(coeffs), pattern="dataset")]) {
            dataset.name <- gsub(dataset.coeff, pattern="dataset", replacement="")
            flag <- dataset.factor == dataset.name
            effect <- coeffs[dataset.coeff,1]
            cat(paste0("Transforming data for ", dataset.name, "by subtracting off effect", effect, "\n"))
            print(head(transformed.esn[st,flag]))
            transformed.esn[st,flag] <- transformed.esn[st,flag] - effect
            print(head(transformed.esn[st,flag]))
        }
    }
    transformed.esn
}


computeGSVA <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.score.field="kras", kras.states=c("MT","WT")) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
    exclude.nras.mutants <- TRUE
    if(exclude.nras.mutants) {
        flag <- !is.na(clin.tmp$nras) & (clin.tmp$nras == 1)
        ## flag <- flag | ( !is.na(clin.tmp$braf) & ( clin.tmp$braf == 1 ) )
        clin.tmp <- clin.tmp[!flag, ]
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
    
    ## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
    ## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
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
    tmp
}

do.msi.old <- function(expr, clin, analysis.name) {

    ## TODO:
    ## matrix_2 label
    ## msi label
    ## braf label
    ## braf values
    ## row/column scale -- what did paper do?
    ## expression should be continuous
    ## add mutation count
    ## cut dendr
    
    ##    suppressPackageStartupMessages(library("NMF"))
    suppressPackageStartupMessages(library("ComplexHeatmap"))
    suppressPackageStartupMessages(library("circlize"))
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(expr))
    clin.mask <- clin[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]
    
    has.msi <- any(!is.na(clin.mask$msi))
    if(has.msi) {
        clin.mask <- clin.mask[!is.na(clin.mask$msi),]
    }

    idxs <- match(clin.mask$sample, colnames(expr.mask))
    clin.mask <- clin.mask[!is.na(idxs),]
    expr.mask <- expr.mask[, na.omit(idxs)]
    
    msi.common <- intersect(msi.signature, rownames(expr.mask))
    print(msi.common)
    print(unique(clin.mask$msi))
    covariates <- c()
    covariate.colors <- list()
    covariate.legend.params <- list()
    if(has.msi) {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row", ColSideColors=unlist(lapply(clin.mask$msi, function(x) ifelse(x=="mss", "black", "red"))))
        ##        covariates <- cbind(covariates, msi=clin.mask$msi)
        covariates <- cbind(msi=clin.mask$msi)
        covariates <- as.data.frame(covariates)        
        covariate.colors[[length(covariate.colors)+1]] <- c("msi" = "black", "mss" = "white")
        names(covariate.colors)[length(covariate.colors)] <- "msi"
        covariate.legend.params[[length(covariate.legend.params)+1]] <- list()
        names(covariate.legend.params[length(covariate.legend.params)]) <- "msi"
    } else {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row")
        ## aheatmap(expr.mask[msi.common,])
    }
    h <- Heatmap(expr.mask[msi.common,])
    if(any(!is.na(clin.mask$braf))) {
        covariates <- cbind(covariates, braf=clin.mask$braf)
        ## covariates <- cbind(braf=as.character(clin.mask$braf))
        covariate.colors[[length(covariate.colors)+1]] <- c("1" = "black", "0" = "white")
        names(covariate.colors)[length(covariate.colors)] <- "braf"
        covariate.legend.params[[length(covariate.legend.params)+1]] <- list()
        names(covariate.legend.params[length(covariate.legend.params)]) <- "braf"
    }
    if("TYMS" %in% rownames(expr.mask)) {
        tyms.expr <- expr.mask["TYMS",]
        covariates <- cbind(covariates, TYMS=tyms.expr)
        ## covariates <- as.data.frame(covariates)
        ## covariates$TYMS <- as.numeric(as.character(covariates$TYMS))
        covariate.colors[[length(covariate.colors)+1]] <- colorRamp2(c(min(tyms.expr), max(tyms.expr)), c("blue", "red"))
        ## covariate.colors[[length(covariate.colors)+1]] <- colorRamp2(c(min(tyms.expr), max(tyms.expr)), c("white", "black"))
        names(covariate.colors)[length(covariate.colors)] <- "TYMS"
        rng <- max(tyms.expr) - min(tyms.expr)
        print(seq(from=min(tyms.expr), to=max(tyms.expr), by=rng/10))
        covariate.legend.params[[length(covariate.legend.params)+1]] <- list(at = seq(from=min(tyms.expr), to=max(tyms.expr), by=rng/10), labels = round(seq(from=min(tyms.expr), to=max(tyms.expr), by=rng/10), digits=1))
        names(covariate.legend.params[length(covariate.legend.params)]) <- "TYMS"
    }
    if(!is.null(ncol(covariates))) {
        covariates <- as.data.frame(covariates)
        if("TYMS" %in% colnames(covariates)) {
            covariates$TYMS <- as.numeric(as.character(covariates$TYMS))
        }
        print(head(covariates))
        ha <- HeatmapAnnotation(df = covariates, col = covariate.colors, annotation_legend_param = covariate.legend.params)
        print(Heatmap(expr.mask[msi.common,], top_annotation = ha, show_column_names = FALSE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 6)))
##        aheatmap(expr.mask[msi.common,], annCol=covariates)
    } else {
##        aheatmap(expr.mask[msi.common,])
##        print(Heatmap(expr.mask[msi.common,]))
    }
}


do.msi.aheatmap <- function(expr, clin, analysis.name, mut.tbl = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification") {

    ## TODO:
    ## msi label
    ## braf label
    ## braf values
    ## rotate
    ## add mutation count
    ## cut dendr

    suppressPackageStartupMessages(library("NMF"))
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(expr))
    clin.mask <- clin[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]
    
    has.msi <- any(!is.na(clin.mask$msi))
    if(has.msi) {
        clin.mask <- clin.mask[!is.na(clin.mask$msi),]
    }

    idxs <- match(clin.mask$sample, colnames(expr.mask))
    clin.mask <- clin.mask[!is.na(idxs),]
    expr.mask <- expr.mask[, na.omit(idxs)]

    mutation.cnt.by.sample <- NULL
    if(!is.null(mut.tbl)) {
        ## 3'Flank
        ## 3'UTR
        ## 5'Flank
        ## 5'UTR
        ## Frame_Shift_Del
        ## Frame_Shift_Ins
        ## IGR
        ## In_Frame_Del
        ## In_Frame_Ins
        ## Intron
        ## Missense_Mutation
        ## Nonsense_Mutation
        ## Nonstop_Mutation
        ## RNA
        ## Silent
        ## Splice_Region
        ## Splice_Site
        ## Translation_Start_Site
        ## Variant_Classification
        flag <- !is.na(mut.tbl[,variant.type.col]) & (mut.tbl[,variant.type.col] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation"))
        mut.tbl <- mut.tbl[flag,]
        print(unique(mut.tbl[,variant.type.col]))

        mut.tbl <- mut.tbl[!duplicated(mut.tbl[,c(sample.col, gene.col)]),]

        suppressPackageStartupMessages(library(reshape2))
##        gene.by.sample.mutations <- acast(data=mut.tbl, formula=gene.col ~ sample.col, value.var=gene.col, fill = "NA")
##        gene.by.sample.mutations <- acast(data=mut.tbl, formula=SYMBOL ~ Tumor_Sample_Barcode, fill = "NA")        

        gene.by.sample.mutations <- acast(data=mut.tbl, formula=as.formula(paste0(gene.col, " ~ ", sample.col)), fill = "NOMUTATION")
        gene.by.sample.mutations.binarized <- apply(gene.by.sample.mutations, c(1,2), function(x) ifelse(x == "NOMUTATION", 0, 1))
        mutation.cnt.by.sample <- colSums(gene.by.sample.mutations.binarized)
        names(mutation.cnt.by.sample) <- gsub(names(mutation.cnt.by.sample), pattern="-01$", replacement="")

##        common.samples <- intersect(names(mutation.cnt.by.sample), colnames(expr.mask))
##        mutation.cnt.by.sample <- mutation.cnt.by.sample[common.samples]
    }
    


    
    msi.common <- intersect(msi.signature, rownames(expr.mask))
    print(msi.common)
    print(unique(clin.mask$msi))
    covariates <- c()
    covariate.colors <- list()
    covariate.legend.params <- list()
    if(has.msi) {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row", ColSideColors=unlist(lapply(clin.mask$msi, function(x) ifelse(x=="mss", "black", "red"))))
        ##        covariates <- cbind(covariates, msi=clin.mask$msi)
        covariates <- cbind(msi=clin.mask$msi)
    } else {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row")
        ## aheatmap(expr.mask[msi.common,])
    }
    if(any(!is.na(clin.mask$braf))) {
        covariates <- cbind(covariates, braf=clin.mask$braf)
    }
    if("TYMS" %in% rownames(expr.mask)) {
        tyms.expr <- expr.mask["TYMS",]
        covariates <- cbind(covariates, TYMS=tyms.expr)
    }
    if(!is.null(ncol(covariates))) {
        covariates <- as.data.frame(covariates)
        if("TYMS" %in% colnames(covariates)) {
            covariates$TYMS <- as.numeric(as.character(covariates$TYMS))
        }
        print(head(covariates))
        ##        aheatmap(expr.mask[msi.common,], annCol=covariates)
        if(!is.null(mutation.cnt.by.sample)) {
            if(TRUE) {
                library('gridExtra')
                layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
                print(aheatmap(t(expr.mask[msi.common,]), annRow=covariates, Colv=FALSE))
                ##            print(plot(x=mutation.cnt.by.sample, y=names(mutation.cnt.by.sample)))
                mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask)]
                names(mutation.cnt.by.sample) <- colnames(expr.mask)
                mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
                print(barplot(mutation.cnt.by.sample, horiz=TRUE))
            } else {
                mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask)]
                names(mutation.cnt.by.sample) <- colnames(expr.mask)
                mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
                mutation.cnt.by.sample[(mutation.cnt.by.sample) == 0] <- 1
                mutation.cnt.by.sample <- log(mutation.cnt.by.sample)
                covariates$mnt.cnt <- unname(mutation.cnt.by.sample)
                print(aheatmap(t(expr.mask[msi.common,]), annRow=covariates, Colv=FALSE))
            }
            
        } else {
            aheatmap(t(expr.mask[msi.common,]), annRow=covariates)
        }
    } else {
##        aheatmap(expr.mask[msi.common,])
##        print(Heatmap(expr.mask[msi.common,]))
    }
}

do.msi.aheatmap2 <- function(expr, clin, analysis.name, mut.tbl = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification") {

    ## TODO:
    ## use ggplot2 from here and here
    ##     http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    ##     http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically/17371177#17371177
    ## msi label
    ## braf label
    ## braf values
    ## rotate
    ## add mutation count
    ## cut dendr

    suppressPackageStartupMessages(library("NMF"))
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(expr))
    clin.mask <- clin[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]
    
    has.msi <- any(!is.na(clin.mask$msi))
    if(has.msi) {
        clin.mask <- clin.mask[!is.na(clin.mask$msi),]
    }

    idxs <- match(clin.mask$sample, colnames(expr.mask))
    clin.mask <- clin.mask[!is.na(idxs),]
    expr.mask <- expr.mask[, na.omit(idxs)]

    mutation.cnt.by.sample <- NULL
    if(!is.null(mut.tbl)) {
        ## 3'Flank
        ## 3'UTR
        ## 5'Flank
        ## 5'UTR
        ## Frame_Shift_Del
        ## Frame_Shift_Ins
        ## IGR
        ## In_Frame_Del
        ## In_Frame_Ins
        ## Intron
        ## Missense_Mutation
        ## Nonsense_Mutation
        ## Nonstop_Mutation
        ## RNA
        ## Silent
        ## Splice_Region
        ## Splice_Site
        ## Translation_Start_Site
        ## Variant_Classification
        flag <- !is.na(mut.tbl[,variant.type.col]) & (mut.tbl[,variant.type.col] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation"))
        mut.tbl <- mut.tbl[flag,]
        print(unique(mut.tbl[,variant.type.col]))

        mut.tbl <- mut.tbl[!duplicated(mut.tbl[,c(sample.col, gene.col)]),]

        suppressPackageStartupMessages(library(reshape2))
##        gene.by.sample.mutations <- acast(data=mut.tbl, formula=gene.col ~ sample.col, value.var=gene.col, fill = "NA")
##        gene.by.sample.mutations <- acast(data=mut.tbl, formula=SYMBOL ~ Tumor_Sample_Barcode, fill = "NA")        

        gene.by.sample.mutations <- acast(data=mut.tbl, formula=as.formula(paste0(gene.col, " ~ ", sample.col)), fill = "NOMUTATION")
        gene.by.sample.mutations.binarized <- apply(gene.by.sample.mutations, c(1,2), function(x) ifelse(x == "NOMUTATION", 0, 1))
        mutation.cnt.by.sample <- colSums(gene.by.sample.mutations.binarized)
        names(mutation.cnt.by.sample) <- gsub(names(mutation.cnt.by.sample), pattern="-01$", replacement="")

##        common.samples <- intersect(names(mutation.cnt.by.sample), colnames(expr.mask))
##        mutation.cnt.by.sample <- mutation.cnt.by.sample[common.samples]
    }
    


    
    msi.common <- intersect(msi.signature, rownames(expr.mask))
    print(msi.common)
    print(unique(clin.mask$msi))
    covariates <- c()
    covariate.colors <- list()
    covariate.legend.params <- list()
    if(has.msi) {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row", ColSideColors=unlist(lapply(clin.mask$msi, function(x) ifelse(x=="mss", "black", "red"))))
        ##        covariates <- cbind(covariates, msi=clin.mask$msi)
        covariates <- cbind(msi=clin.mask$msi)
    } else {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row")
        ## aheatmap(expr.mask[msi.common,])
    }
    if(any(!is.na(clin.mask$braf))) {
        covariates <- cbind(covariates, braf=clin.mask$braf)
    }
    if("TYMS" %in% rownames(expr.mask)) {
        tyms.expr <- expr.mask["TYMS",]
        covariates <- cbind(covariates, TYMS=tyms.expr)
    }
    if(!is.null(ncol(covariates))) {
        covariates <- as.data.frame(covariates)
        if("TYMS" %in% colnames(covariates)) {
            covariates$TYMS <- as.numeric(as.character(covariates$TYMS))
        }
        print(head(covariates))
        ##        aheatmap(expr.mask[msi.common,], annCol=covariates)
        if(!is.null(mutation.cnt.by.sample)) {
            if(TRUE) {
                library('gridExtra')
                layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
                ah <- aheatmap(t(expr.mask[msi.common,]), annRow=covariates, Colv=FALSE)
                clusters <- cutree(as.hclust(ah$Rowv), k=2)
                ## save(ah, file="ah.Rd")
                print(ah)
                ##            print(plot(x=mutation.cnt.by.sample, y=names(mutation.cnt.by.sample)))
                mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask)]
                names(mutation.cnt.by.sample) <- colnames(expr.mask)
##                mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
                print(barplot(mutation.cnt.by.sample, horiz=TRUE))
                df <- data.frame(cnt=unname(mutation.cnt.by.sample[!is.na(mutation.cnt.by.sample)]))
                rownames(df) <- names(mutation.cnt.by.sample[!is.na(mutation.cnt.by.sample)])
                cls.df <- data.frame(cluster=clusters)
                rownames(cls.df) <- names(clusters)
                clusters <- merge(cls.df, df, by="row.names", both=TRUE)
                df <- clusters
                df$cluster <- factor(df$cluster)
                df <- df[df$cnt > 0,]
                df$cnt <- log(df$cnt)
                save(df, file="df.Rd")
                print(df)
                p <- ggplot(data=df, aes(x=cluster, y=cnt))
    
                ## Create a box plot where the x axis is KRAS mutation status ...
                p <- p + geom_boxplot(aes(fill=cluster))
                ##                p <- p + geom_jitter()
                p <- p + geom_beeswarm()
                print(p)
            } else {
                mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask)]
                names(mutation.cnt.by.sample) <- colnames(expr.mask)
                mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
                mutation.cnt.by.sample[(mutation.cnt.by.sample) == 0] <- 1
                mutation.cnt.by.sample <- log(mutation.cnt.by.sample)
                covariates$mnt.cnt <- unname(mutation.cnt.by.sample)
                print(aheatmap(t(expr.mask[msi.common,]), annRow=covariates, Colv=FALSE))
            }
            
        } else {
            aheatmap(t(expr.mask[msi.common,]), annRow=covariates)
        }
    } else {
##        aheatmap(expr.mask[msi.common,])
##        print(Heatmap(expr.mask[msi.common,]))
    }
}

do.msi.grob <- function(expr, clin, analysis.name, mut.tbl = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification") {

    ## TODO:
    ## use ggplot2 from here and here
    ##     http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    ##     http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically/17371177#17371177
    ## msi label
    ## braf label
    ## braf values
    ## rotate
    ## add mutation count
    ## cut dendr

    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("gtable"))    
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(expr))
    clin.mask <- clin[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]
    
    has.msi <- any(!is.na(clin.mask$msi))
    if(has.msi) {
        clin.mask <- clin.mask[!is.na(clin.mask$msi),]
    }

    idxs <- match(clin.mask$sample, colnames(expr.mask))
    clin.mask <- clin.mask[!is.na(idxs),]
    expr.mask <- expr.mask[, na.omit(idxs)]

    msi.common <- intersect(msi.signature, rownames(expr.mask))
    print(msi.common)

    ## Z-score the genes/rows
    expr.mask <- expr.mask[msi.common,]
    expr.mask <- t(scale(t(expr.mask), center=TRUE, scale=TRUE))

    dd.row <- as.dendrogram(hclust(dist(expr.mask)))
    row.ord <- order.dendrogram(dd.row)
    
    dd.col <- as.dendrogram(hclust(dist(t(expr.mask))))
    col.ord <- order.dendrogram(dd.col)

##    ddata_x <- dendro_data(dd.row)
##    ddata_y <- dendro_data(dd.col)
    
    expr.melt <- melt(expr.mask[row.ord, col.ord])
    p <- ggplot(data = expr.melt, aes(x = Var1, y=Var2))
    p <- p + geom_tile(aes(fill = value), show.legend = FALSE)
##    p <- p + scale_fill_gradient2(high = "green", mid = "black", low = "red") 
##    print(g)
    
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())    
    
    g <- ggplotGrob(p)

    mutation.cnt.by.sample <- NULL
    if(!is.null(mut.tbl)) {
        ## 3'Flank
        ## 3'UTR
        ## 5'Flank
        ## 5'UTR
        ## Frame_Shift_Del
        ## Frame_Shift_Ins
        ## IGR
        ## In_Frame_Del
        ## In_Frame_Ins
        ## Intron
        ## Missense_Mutation
        ## Nonsense_Mutation
        ## Nonstop_Mutation
        ## RNA
        ## Silent
        ## Splice_Region
        ## Splice_Site
        ## Translation_Start_Site
        ## Variant_Classification
        flag <- !is.na(mut.tbl[,variant.type.col]) & (mut.tbl[,variant.type.col] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation"))
        mut.tbl <- mut.tbl[flag,]
        print(unique(mut.tbl[,variant.type.col]))

        mut.tbl <- mut.tbl[!duplicated(mut.tbl[,c(sample.col, gene.col)]),]

        suppressPackageStartupMessages(library(reshape2))
        gene.by.sample.mutations <- acast(data=mut.tbl, formula=as.formula(paste0(gene.col, " ~ ", sample.col)), fill = "NOMUTATION")
        gene.by.sample.mutations.binarized <- apply(gene.by.sample.mutations, c(1,2), function(x) ifelse(x == "NOMUTATION", 0, 1))
        mutation.cnt.by.sample <- colSums(gene.by.sample.mutations.binarized)
        names(mutation.cnt.by.sample) <- gsub(names(mutation.cnt.by.sample), pattern="-01$", replacement="")
        mutation.cnt.by.sample <- data.frame(x = names(mutation.cnt.by.sample), y = unname(mutation.cnt.by.sample))
    }
    
    print(g)    
    g <- gtable_add_cols(g, unit(1,"cm"))
##    g <- gtable_add_grob(g, rectGrob(gp=gpar(fill="red")), t = 6, l=ncol(g), b=6, r=ncol(g))
    p.mut <- ggplot(data = mutation.cnt.by.sample, aes(x = x, y = y))
    p.mut <- p.mut + geom_bar(stat = "identity")
    g.mut <- ggplotGrob(p.mut)
    
    
    g <- gtable_add_cols(g, unit(1,"cm"))
    ##    g <- gtable_add_grob(g, rectGrob(gp=gpar(fill="black")), t = 6, l=ncol(g), b=6, r=ncol(g))
    g <- gtable_add_grob(g, g.mut, t = 6, l=ncol(g), b=6, r=ncol(g))

    g <- gtable_add_rows(g, unit(1,"in"), 0)
    g <- gtable_add_grob(g, rectGrob(gp=gpar(fill="blue")),
                         t = 1, l=4, b=1, r=4)
    print(g)
    print(ncol(g))
    print(nrow(g))    
    grid.newpage()
    grid.draw(g)
    
    return()
    
    print(unique(clin.mask$msi))
    covariates <- c()
    covariate.colors <- list()
    covariate.legend.params <- list()
    if(has.msi) {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row", ColSideColors=unlist(lapply(clin.mask$msi, function(x) ifelse(x=="mss", "black", "red"))))
        ##        covariates <- cbind(covariates, msi=clin.mask$msi)
        covariates <- cbind(msi=clin.mask$msi)
    } else {
        ##        heatmap.2(expr.mask[msi.common,], trace="none", scale="row")
        ## aheatmap(expr.mask[msi.common,])
    }
    if(any(!is.na(clin.mask$braf))) {
        covariates <- cbind(covariates, braf=clin.mask$braf)
    }
    if("TYMS" %in% rownames(expr.mask)) {
        tyms.expr <- expr.mask["TYMS",]
        covariates <- cbind(covariates, TYMS=tyms.expr)
    }
    if(!is.null(ncol(covariates))) {
        covariates <- as.data.frame(covariates)
        if("TYMS" %in% colnames(covariates)) {
            covariates$TYMS <- as.numeric(as.character(covariates$TYMS))
        }
        print(head(covariates))
        ##        aheatmap(expr.mask[msi.common,], annCol=covariates)
        if(!is.null(mutation.cnt.by.sample)) {
            if(TRUE) {
                library('gridExtra')
                layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
                ah <- aheatmap(t(expr.mask[msi.common,]), annRow=covariates, Colv=FALSE)
                clusters <- cutree(as.hclust(ah$Rowv), k=2)
                ## save(ah, file="ah.Rd")
                print(ah)
                ##            print(plot(x=mutation.cnt.by.sample, y=names(mutation.cnt.by.sample)))
                mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask)]
                names(mutation.cnt.by.sample) <- colnames(expr.mask)
##                mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
                print(barplot(mutation.cnt.by.sample, horiz=TRUE))
                df <- data.frame(cnt=unname(mutation.cnt.by.sample[!is.na(mutation.cnt.by.sample)]))
                rownames(df) <- names(mutation.cnt.by.sample[!is.na(mutation.cnt.by.sample)])
                cls.df <- data.frame(cluster=clusters)
                rownames(cls.df) <- names(clusters)
                clusters <- merge(cls.df, df, by="row.names", both=TRUE)
                df <- clusters
                df$cluster <- factor(df$cluster)
                df <- df[df$cnt > 0,]
                df$cnt <- log(df$cnt)
                save(df, file="df.Rd")
                print(df)
                p <- ggplot(data=df, aes(x=cluster, y=cnt))
    
                ## Create a box plot where the x axis is KRAS mutation status ...
                p <- p + geom_boxplot(aes(fill=cluster))
                ##                p <- p + geom_jitter()
                p <- p + geom_beeswarm()
                print(p)
            } else {
                mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask)]
                names(mutation.cnt.by.sample) <- colnames(expr.mask)
                mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
                mutation.cnt.by.sample[(mutation.cnt.by.sample) == 0] <- 1
                mutation.cnt.by.sample <- log(mutation.cnt.by.sample)
                covariates$mnt.cnt <- unname(mutation.cnt.by.sample)
                print(aheatmap(t(expr.mask[msi.common,]), annRow=covariates, Colv=FALSE))
            }
            
        } else {
            aheatmap(t(expr.mask[msi.common,]), annRow=covariates)
        }
    } else {
##        aheatmap(expr.mask[msi.common,])
##        print(Heatmap(expr.mask[msi.common,]))
    }
}



do.msi <- function(expr, clin, file, mut.tbl = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification") {

    ## Based on:
    ##     http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    ##     http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically/17371177#17371177

    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("gtable"))
    suppressPackageStartupMessages(library("ggdendro"))    

    clin.orig <- clin
    
    mut.tbl[,sample.col] <- gsub(mut.tbl[,sample.col], pattern="-01$", replacement="")
    limit.size <- FALSE
##    limit.size <- TRUE
    if(limit.size) {
        expr <- expr[,1:200]
    }
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(expr))
    clin.mask <- clin[!is.na(idxs),]
    expr.mask <- expr[, na.omit(idxs)]

    has.msi <- any(!is.na(clin.mask$msi))

    if(has.msi) {
        clin.mask <- clin.mask[!is.na(clin.mask$msi),]
    }

    has.braf <- any(!is.na(clin.mask$braf))
    has.tyms <- "TYMS" %in% rownames(expr.mask)
    has.genomic <- !is.null(mut.tbl)
    num.annotation.legends <- 0
    if(has.msi) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    if(has.braf) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    if(has.tyms) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    num.cols <- num.annotation.legends
    if(has.genomic) {
        num.cols <- num.cols + 1
    }
    ## Add to column count for expression heatmap and dendrogram
    num.cols <- num.cols + 2
    ## And one for the legends
    if(num.annotation.legends > 0) {
        num.cols <- num.cols + 1
    }
    
    idxs <- match(clin.mask$sample, colnames(expr.mask))
    clin.mask <- clin.mask[!is.na(idxs),]
    expr.mask <- expr.mask[, na.omit(idxs)]

    ## Scale/z-score the rows/genes of the matrix
    ## Don't cluster scaled values.
    ##    expr.mask <- t(scale(t(expr.mask), center=TRUE, scale=TRUE))
    ##    print(sum(expr.mask[1,]))
    ##    print(sum(expr.mask[2,]))
    ##    expr.mask <- (scale((expr.mask), center=TRUE, scale=TRUE))

    expr.mask.with.tyms <- expr.mask
    msi.common <- intersect(msi.signature, rownames(expr.mask))

    expr.mask <- expr.mask[msi.common,]

    method <- "ward.D"
    method <- "complete"
    dd.row <- as.dendrogram(hclust(dist(expr.mask), method=method))
    row.ord <- order.dendrogram(dd.row)

    hc.col <- hclust(dist(t(expr.mask)), method=method)
    dd.col <- as.dendrogram(hc.col)
    col.ord <- order.dendrogram(dd.col)

    ## Find the separation between the MSI and MSS cases
    clusters <- cutree(hc.col, k=2)
    clusters <- clusters[col.ord]
    cluster.freq <- as.data.frame(table(clusters), stringsAsFactors=FALSE)
    print(cluster.freq)
    msi.cluster <- 1
    mss.cluster <- 2
    ## Assume (for now) that cluster with more members is the MSS cluster.
    ## In the future, look at the expression of the genes
    min.freq.cluster <- cluster.freq$clusters[which(cluster.freq$Freq == min(cluster.freq$Freq))]
    msi.cluster <- as.numeric(min.freq.cluster)
    cat(paste0("MSI cluster: ", msi.cluster, "\n"))
    mss.cluster <- 2 - msi.cluster + 1
    cat(paste0("MSS cluster: ", mss.cluster, "\n"))    

    msi.assign.tbl <- data.frame(sample=names(clusters), msi.inferred=unname(clusters), stringsAsFactors=FALSE)
    msi.assign.tbl$msi.inferred[msi.assign.tbl$msi.inferred==msi.cluster] <- "msi"
    msi.assign.tbl$msi.inferred[msi.assign.tbl$msi.inferred==mss.cluster] <- "mss"    

    clin.orig <- merge(clin.orig, msi.assign.tbl, by="sample", all=TRUE)
    
    cross.over.sample <- NULL
    for(i in 2:length(clusters)) {
        if(clusters[i-1] != clusters[i]) {
            cross.over.sample <- names(clusters)[i]
        }
    }

    mutation.cnt.by.sample <- NULL
    if(!is.null(mut.tbl)) {
        ## 3'Flank
        ## 3'UTR
        ## 5'Flank
        ## 5'UTR
        ## Frame_Shift_Del
        ## Frame_Shift_Ins
        ## IGR
        ## In_Frame_Del
        ## In_Frame_Ins
        ## Intron
        ## Missense_Mutation
        ## Nonsense_Mutation
        ## Nonstop_Mutation
        ## RNA
        ## Silent
        ## Splice_Region
        ## Splice_Site
        ## Translation_Start_Site
        ## Variant_Classification
        mut.tbl <- mut.tbl[mut.tbl[,sample.col] %in% colnames(expr),]
        flag <- !is.na(mut.tbl[,variant.type.col]) & (mut.tbl[,variant.type.col] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation"))
        mut.tbl <- mut.tbl[flag,]

        mut.tbl <- mut.tbl[!duplicated(mut.tbl[,c(sample.col, gene.col)]),]

        suppressPackageStartupMessages(library(reshape2))
        gene.by.sample.mutations <- acast(data=mut.tbl, formula=as.formula(paste0(gene.col, " ~ ", sample.col)), fill = "NOMUTATION")
        gene.by.sample.mutations.binarized <- apply(gene.by.sample.mutations, c(1,2), function(x) ifelse(x == "NOMUTATION", 0, 1))
        mutation.cnt.by.sample <- colSums(gene.by.sample.mutations.binarized)
        names(mutation.cnt.by.sample) <- gsub(names(mutation.cnt.by.sample), pattern="-01$", replacement="")
        mutation.cnt.by.sample <- mutation.cnt.by.sample[colnames(expr.mask[, col.ord])]
        mutation.cnt.by.sample[is.na(mutation.cnt.by.sample)] <- 0
        names(mutation.cnt.by.sample) <- colnames(expr.mask[, col.ord])
    }

    ## Re-order the expression matrix based on the clustering
    expr.mask <- (expr.mask[row.ord, col.ord])
    expr.mask.with.tyms <- expr.mask.with.tyms[, col.ord]    
    clin.mask <- clin.mask[col.ord, ]
    
    remove.msi.labels <- TRUE
    remove.labels <- TRUE
    remove.braf.labels <- TRUE
    remove.tyms.labels <- TRUE
    remove.mut.labels <- TRUE

    if(!remove.msi.labels || !remove.braf.labels || !remove.tyms.labels || !remove.mut.labels) {
        ## Set all put every 20th to just a number so we can read
        flag <- rep_len(c(FALSE,rep(TRUE,9)), ncol(expr.mask))
        print(any(is.na(colnames(expr.mask))))
        colnames(expr.mask)[flag] <- as.character((1:ncol(expr.mask))[flag])
        colnames(expr.mask.with.tyms)[flag] <- as.character((1:ncol(expr.mask.with.tyms))[flag])        
        print(head(colnames(expr.mask),20))
        print(any(is.na(colnames(expr.mask))))        
        clin.mask$sample[flag] <- as.character((1:length(clin.mask$sample))[flag])
        print(any(is.na(clin.mask$sample[flag])))
    }
    
    ## Scale the rows/genes of the matrix
    ## NB: we are doing it here--after clustering!  So it will only effect
    ## the visuals.
    ## Scale/z-score the rows/genes of the matrix
    expr.mask <- t(scale(t(expr.mask), center=TRUE, scale=TRUE))
    print(sum(expr.mask[1,]))
    print(sum(expr.mask[2,]))

    ## It is important that the samples be factors to ensure consistency
    ## across the plots
    ## MATRIX TRANSPOSE!
    xx <- t(expr.mask)
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$sample <- xx_names[[1]]
    
    levels <- df$sample
    df$sample <- with(df, factor(sample, levels=levels, ordered=TRUE))

    ## Find the index of the sample at which the dendrogram transitions
    ## between MSI and MSS
    cross.over.sample.indx <- which(levels(df$sample) == cross.over.sample)

    if(has.genomic) {
        ## Make the mutation samples factors
        mutation.cnt.by.sample <- mutation.cnt.by.sample[levels(df$sample)]
        mutation.cnt.by.sample <- data.frame(sample = names(mutation.cnt.by.sample), cnt = unname(mutation.cnt.by.sample))
        mutation.cnt.by.sample$sample <- factor(mutation.cnt.by.sample$sample, levels=levels, ordered=TRUE)
    }
    
    ddata_y <- dendro_data(dd.col)

    ## Create the dendrogram plot--this will be in the first column

    p.sample.dendro <- ggplot(segment(ddata_y)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
    p.sample.dendro <- p.sample.dendro + coord_flip() + scale_y_reverse()
    p.sample.dendro <- p.sample.dendro + geom_vline(xintercept = cross.over.sample.indx, linetype = 'dashed')    
    p.sample.dendro <- p.sample.dendro + theme(plot.margin=unit(c(0,0,0,0), "cm"))

    segs <- ddata_y$segments
    ## Remove the extra space so we can align with the heatmap
    p.sample.dendro <- p.sample.dendro + scale_x_continuous(expand = c(0,0), limits=range(segs[,1]))
    p.sample.dendro <- p.sample.dendro + theme(panel.grid=element_blank(), panel.background=element_rect(fill = "transparent",colour = NA), panel.border=element_blank())

    dendro.grob <- ggplotGrob(p.sample.dendro)
    dendro.grob <- gtable::gtable_filter(dendro.grob, "panel")    

    ## Build up the rows for the gtable.
    indx <- 1
    title.row <- c(NA)    ## No title for the dendrogram
    plot.row <- c(indx)
    grobs <- list()
    grobs[[indx]] <- dendro.grob
    indx <- indx + 1

    fs <- 10
    heat.just <- c(0.5, 2)
    just <- c(0.5, 0.5)
    just <- c(1, 1)
    just <- c(1, 0)
    just <- c(0, 0)
    just <- c(0.5, 0)
    just <- c(0.5, 0.5)                
    
    ## Create the heatmap and its title (second column)
    heat.label <- textGrob("Genes in MSI Signature", gp=gpar(fontsize=fs), just=heat.just)

    expr.melt <- melt(df, id.vars="sample")
    p.heat <- ggplot(data = expr.melt, aes(x = variable, y=sample))
    p.heat <- p.heat + geom_tile(aes(fill = value), show.legend = TRUE)
    p.heat <- p.heat + geom_hline(yintercept = cross.over.sample.indx, linetype = 'dashed')
    p.heat <- p.heat + guides(fill = guide_colourbar(title = "Expression", direction = "horizontal", title.position = "top", title.hjust = 0.5))
    p.heat <- p.heat + theme(legend.text = element_text(size = fs))
    
    if(remove.labels) {
        p.heat <- p.heat + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        p.heat <- p.heat + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())    
    }
    
    g.heat <- ggplotGrob(p.heat)
    g.legend.heat <- gtable::gtable_filter(g.heat, "guide-box")
    if(remove.labels) {
        g.heat <- gtable::gtable_filter(g.heat, "panel")    
    }
    
    title.row <- c(title.row, indx)
    grobs[[indx]] <- heat.label
    indx <- indx + 1
    
    plot.row <- c(plot.row, indx)
    grobs[[indx]] <- g.heat
    indx <- indx + 1

    g.legend.msi <- NULL
    if(has.msi) {
        cat("Has MSI\n")
        msi.tbl <- data.frame(sample=clin.mask$sample, msi=clin.mask$msi)
        msi.tbl$msi <- as.character(msi.tbl$msi)
        msi.tbl$msi[msi.tbl$msi == "msi"] <- "MSI"
        msi.tbl$msi[msi.tbl$msi == "mss"] <- "MSS"

        msi.tbl$sample <- factor(msi.tbl$sample, levels=levels, ordered=TRUE)
        
        msi.melt <- melt(msi.tbl, id.vars="sample")
        p.msi <- ggplot(data = msi.melt, aes(x = variable, y=sample))
        p.msi <- p.msi + geom_tile(aes(fill = value), show.legend = TRUE)
        p.msi <- p.msi + scale_fill_manual(values = c("black", "gray"), na.value = "white")
        p.msi <- p.msi + guides(fill = guide_legend(title="MSI\nStatus", title.hjust=0.5))
        p.msi <- p.msi + theme(legend.text = element_text(size = fs))
        p.msi <- p.msi + geom_hline(yintercept = cross.over.sample.indx, linetype = 'dashed')    
        g.msi <- ggplotGrob(p.msi)
        g.legend.msi <- gtable::gtable_filter(g.msi, "guide-box")    
        g.msi <- gtable::gtable_filter(g.msi, "panel")

        if(!remove.msi.labels) {
            g.msi <- ggplotGrob(p.msi)
        } 
        
        msi.label <- textGrob("MSI", gp=gpar(fontsize=fs), just=just, rot=90)

        title.row <- c(title.row, indx)
        grobs[[indx]] <- msi.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.msi
        indx <- indx + 1
    }

    g.legend.braf <- NULL
    if(has.braf) {
        cat("Has BRAF\n")        
        braf.tbl <- data.frame(sample=clin.mask$sample, braf=clin.mask$braf)
        braf.tbl$sample <- factor(braf.tbl$sample, levels=levels, ordered=TRUE)
        braf.tbl$braf <- as.character(braf.tbl$braf)
        braf.tbl$braf[braf.tbl$braf == "1"] <- "MT"
        braf.tbl$braf[braf.tbl$braf == "0"] <- "WT"
        
        braf.melt <- melt(braf.tbl, id.vars="sample")
        p.braf <- ggplot(data = braf.melt, aes(x = variable, y=sample))
        p.braf <- p.braf + geom_tile(aes(fill = value), show.legend = TRUE)
        p.braf <- p.braf + scale_fill_manual(values = c("black", "gray"), na.value = "white")    
        p.braf <- p.braf + guides(fill = guide_legend(title="BRAF", title.hjust = 0.5, title.theme=element_text(face = 'italic', angle = 0)))
        p.braf <- p.braf + theme(legend.text = element_text(size = fs))        
        p.braf <- p.braf + geom_hline(yintercept = cross.over.sample.indx, linetype = 'dashed')    
        g.braf <- ggplotGrob(p.braf)
        g.legend.braf <- gtable::gtable_filter(g.braf, "guide-box")
        g.braf <- gtable::gtable_filter(g.braf, "panel")

        if(!remove.braf.labels) {
            g.braf <- ggplotGrob(p.braf)
        }
        
        braf.label <- textGrob("BRAF", gp=gpar(fontsize=fs, fontface='italic'), just=just, rot=90)
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- braf.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.braf
        indx <- indx + 1
    }

    g.legend.tyms <- NULL
    if(has.tyms) {
        cat("Has TYMS\n")        
        tyms.expr <- expr.mask.with.tyms["TYMS",]
        tyms.expr.sd <- sd(tyms.expr)
        tyms.expr <- ( tyms.expr - mean(tyms.expr) ) / tyms.expr.sd
        tyms.tbl <- data.frame(sample=names(tyms.expr), tyms=tyms.expr)
        tyms.tbl$sample <- factor(tyms.tbl$sample, levels=levels, ordered=TRUE)
        
        tyms.melt <- melt(tyms.tbl, id.vars="sample")
        p.tyms <- ggplot(data = tyms.melt, aes(x = variable, y=sample))
        p.tyms <- p.tyms + geom_tile(aes(fill = value), show.legend = TRUE)
        p.tyms <- p.tyms + scale_fill_gradient2(high = "green", mid = "black", low = "red")         
        ##        p.tyms <- p.tyms + guides(fill = guide_colourbar(title = Expression(italic(TYMS)~Expression), direction = "vertical", title.position = "top", title.hjust = 0.5))
        p.tyms <- p.tyms + guides(fill = guide_colourbar(title = "TYMS", direction = "vertical", title.position = "top", title.hjust = 0.5, title.theme=element_text(face = 'italic', angle = 0)))        
        p.tyms <- p.tyms + theme(legend.text = element_text(size = fs))
        p.tyms <- p.tyms + geom_hline(yintercept = cross.over.sample.indx, linetype = 'dashed')
        g.tyms <- ggplotGrob(p.tyms)
        g.legend.tyms <- gtable::gtable_filter(g.tyms, "guide-box")        
        g.tyms <- gtable::gtable_filter(g.tyms, "panel")

        if(!remove.tyms.labels) {
            g.tyms <- ggplotGrob(p.tyms)
        }
        
        tyms.label <- textGrob("TYMS", gp=gpar(fontsize=fs, fontface='italic'), just=just, rot=90)
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- tyms.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.tyms
        indx <- indx + 1
    }
    
    if(has.genomic) {
        cat("Has Genomic\n")        
        p.mut <- ggplot(data = mutation.cnt.by.sample, aes(x = sample, y = cnt))
        p.mut <- p.mut + geom_bar(stat = "identity")
        p.mut <- p.mut + coord_flip()
        if(remove.mut.labels) {
            p.mut <- p.mut + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
            p.mut <- p.mut + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
        }
        p.mut <- p.mut + geom_vline(xintercept = cross.over.sample.indx, linetype = 'dashed')    
        
        g.mut <- ggplotGrob(p.mut)
        if(remove.mut.labels) {        
            g.mut <- gtable::gtable_filter(g.mut, "panel")
        }
        str <- paste(strwrap("Num Mutations", width=12), collapse="\n")
        mut.label <- textGrob(str, gp=gpar(fontsize=fs))
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- mut.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.mut
        indx <- indx + 1
    }

    ## Introduce a spacer row between plots and legends and a second
    ## below the legends
    spacer.row <- rep(NA, num.cols)
    bottom.legend.row <- spacer.row
    bottom.legend.row[2] <- indx
    grobs[[indx]] <- g.legend.heat
    indx <- indx + 1

    if(num.annotation.legends > 0) {
        title.row <- c(title.row, NA)
    }
    
    layout <- title.row

    widths <- c(1, 5)
    annotation.width <- 0.25
    
    ## Add the annotation legends
    if(has.msi) {
        row <- c(plot.row, indx)
        print(layout)
        layout <- rbind(layout, row)
        print(layout)
        grobs[[indx]] <- g.legend.msi
        indx <- indx + 1
        widths <- c(widths, annotation.width)
    }
    if(has.braf) {
        row <- c(plot.row, indx)
##        row <- spacer.row
##        row[num.cols] <- indx
        layout <- rbind(layout, row)
        grobs[[indx]] <- g.legend.braf
        indx <- indx + 1
        widths <- c(widths, annotation.width)        
    }

    if(has.tyms) {
        row <- c(plot.row, indx)
##        row <- spacer.row
##        row[num.cols] <- indx
        layout <- rbind(layout, row)
        grobs[[indx]] <- g.legend.tyms
        indx <- indx + 1
        widths <- c(widths, annotation.width)        
    }
    if(has.genomic) {
        widths <- c(widths, 1)
    }
    if(num.annotation.legends > 0) {
        widths <- c(widths, 2)
    }
    layout <- rbind(layout, spacer.row, bottom.legend.row, spacer.row)
    rownames(layout) <- NULL
    heights <- c(1,rep(10/max(1,num.annotation.legends), max(1,num.annotation.legends)),0.5,1,0.5)
    
##    lay <- cbind(c(NA,7,NA,NA,NA),c(2,8,NA,13,NA),c(3,9,NA,14,NA),c(4,10,NA,15,NA),c(5,11,NA,NA,NA),c(6,12,NA,NA,NA))
##    lay <- t(lay)
##    lay <- rbind(c(NA,2,3,4,5,6),c(7,8,9,10,11,12), c())

    print(widths)
    print(heights)
    print(layout)
    print(length(grobs))
##    grobs <- list(heat.label, msi.label, braf.label, tyms.label, mut.label, dendro.grob, g.heat, g.msi, g.braf, g.tyms, g.mut, g.legend.heat, g.legend.msi, g.legend.braf)
    ##    grob.table <- grid.arrange(grobs=grobs, layout_matrix=lay, heights=c(1,10,0.5,1,0.5), widths=c(1, 5, 0.5, 0.5, 0.5, 1))
    png(paste0(file, ".png"))
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    d <- dev.off()
    
    ## tyms different color
    return(clin.orig)
}




plot.cibersort <- function(cibersort.mat, clin, file) {

    ## Based on:
    ##     http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    ##     http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically/17371177#17371177

    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("gtable"))
    suppressPackageStartupMessages(library("ggdendro"))    

    clin.orig <- clin
    tmp <- clin
    tmp <- tmp[!is.na(tmp$cms_label) & (tmp$cms_label %in% c("CMS1", "CMS2", "CMS3", "CMS4")),]
    tmp$sample <- gsub(tmp$sample, pattern="-", replacement=".")

    cibersort.tmp <- cibersort.mat
##    cibersort.tmp <- cibersort.tmp[cibersort.tmp$P.value < 0.1,]
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(tmp$sample, cibersort.tmp$Input.Sample)
    clin.mask <- tmp[!is.na(idxs),]
    cibersort.mask <- cibersort.tmp[na.omit(idxs), ]
    rownames(cibersort.mask) <- cibersort.mask$Input.Sample
    cibersort.data <- cibersort.mask[,!(colnames(cibersort.mask) %in% c("Input.Sample", "P.value", "Pearson.Correlation", "RMSE"))]

    ## We will annotate with RAS status and CMS
    num.annotation.legends <- 2
    num.cols <- num.annotation.legends

    ## Add to column count for cibersort heatmap and dendrogram
    num.cols <- num.cols + 2
    ## And one for the legends
    if(num.annotation.legends > 0) {
        num.cols <- num.cols + 1
    }
    
    method <- "ward.D"
    method <- "complete"
    dd.row <- as.dendrogram(hclust(dist(cibersort.data), method=method))
    row.ord <- order.dendrogram(dd.row)

    hc.col <- hclust(dist(t(cibersort.data)), method=method)
    dd.col <- as.dendrogram(hc.col)
    col.ord <- order.dendrogram(dd.col)

    ## Re-order the samples in the cibersort matrix and clin based on the clustering
    cibersort.mask <- (cibersort.mask[row.ord, ])
    cibersort.data <- (cibersort.data[row.ord, ])    
    clin.mask <- clin.mask[row.ord, ]

    remove.labels <- TRUE
    remove.kras.labels <- TRUE
    remove.cms.labels <- TRUE

    if(!remove.labels || !remove.kras.labels || !remove.cms.labels) {
        ## Set all put every 20th to just a number so we can read
        flag <- rep_len(c(FALSE,rep(TRUE,9)), nrow(cibersort.mask))
        rownames(cibersort.mask)[flag] <- as.character((1:nrow(cibersort.mask))[flag])
        rownames(cibersort.data)[flag] <- as.character((1:nrow(cibersort.data))[flag])        
        clin.mask$sample[flag] <- as.character((1:length(clin.mask$sample))[flag])
    }

    ## NB: we are not scaling the columns/rows of the cibersort data

    ## Scale the rows/genes of the matrix
    ## NB: we are doing it here--after clustering!  So it will only effect
    ## the visuals.
    ## Scale/z-score the rows/genes of the matrix
##    cibersort.data <- t(scale(t(cibersort.data), center=TRUE, scale=TRUE))
##    cibersort.data <- (scale((cibersort.data), center=TRUE, scale=TRUE))    

    
    ## It is important that the samples be factors to ensure consistency
    ## across the plots
    xx <- cibersort.data
    xx_names <- dimnames(xx)
##    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$sample <- xx_names[[1]]
    
    levels <- df$sample

    df$sample <- with(df, factor(sample, levels=levels, ordered=TRUE))

    ddata_y <- dendro_data(dd.row)

    
    ## Create the dendrogram plot--this will be in the first column

    p.sample.dendro <- ggplot(segment(ddata_y)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
    p.sample.dendro <- p.sample.dendro + coord_flip() + scale_y_reverse()
    p.sample.dendro <- p.sample.dendro + theme(plot.margin=unit(c(0,0,0,0), "cm"))

    segs <- ddata_y$segments
    ## Remove the extra space so we can align with the heatmap
    p.sample.dendro <- p.sample.dendro + scale_x_continuous(expand = c(0,0), limits=range(segs[,1]))
    p.sample.dendro <- p.sample.dendro + theme(panel.grid=element_blank(), panel.background=element_rect(fill = "transparent",colour = NA), panel.border=element_blank())

    dendro.grob <- ggplotGrob(p.sample.dendro)
    dendro.grob <- gtable::gtable_filter(dendro.grob, "panel")    

    ## Build up the rows for the gtable.
    indx <- 1
    title.row <- c(NA)    ## No title for the dendrogram
    plot.row <- c(indx)
    grobs <- list()
    grobs[[indx]] <- dendro.grob
    indx <- indx + 1

    fs <- 10
    heat.just <- c(0.5, 2)
    just <- c(0.5, 0.5)
    just <- c(1, 1)
    just <- c(1, 0)
    just <- c(0, 0)
    just <- c(0.5, 0)
    just <- c(0.5, 0.5)                
    just <- c(1, 1)
    
    ## Create the heatmap and its title (second column)
    heat.label <- textGrob("Genes in MSI Signature", gp=gpar(fontsize=fs), just=heat.just)

    expr.melt <- melt(df, id.vars="sample")
    p.heat <- ggplot(data = expr.melt, aes(x = variable, y=sample))
    p.heat <- p.heat + geom_tile(aes(fill = value), show.legend = TRUE)
    p.heat <- p.heat + guides(fill = guide_colourbar(title = "Cellular Fraction", direction = "horizontal", title.position = "top", title.hjust = 0.5))
    p.heat <- p.heat + theme(legend.text = element_text(size = fs))
    
    if(remove.labels) {
        ##        p.heat <- p.heat + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        p.heat <- p.heat + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
##        p.heat <- p.heat + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())    
    }

    g.heat <- ggplotGrob(p.heat)
    heat.label <- gtable::gtable_filter(g.heat, "axis-b")
    g.legend.heat <- gtable::gtable_filter(g.heat, "guide-box")
    print(g.heat)
    if(remove.labels) {
        g.heat <- gtable::gtable_filter(g.heat, "panel")    
    }
    
    title.row <- c(title.row, indx)
    grobs[[indx]] <- heat.label
    indx <- indx + 1
    
    plot.row <- c(plot.row, indx)
    grobs[[indx]] <- g.heat
    indx <- indx + 1

    g.legend.ras <- NULL
    has.kras <- TRUE
    if(has.kras) {
        cat("Has KRAS\n")
        kras.tbl <- data.frame(sample=clin.mask$sample, kras=clin.mask$kras)
        kras.tbl$kras <- as.character(kras.tbl$kras)
        kras.tbl$kras[kras.tbl$kras == "0"] <- "WT"
        kras.tbl$kras[kras.tbl$kras == "1"] <- "MT"

        kras.tbl$sample <- factor(kras.tbl$sample, levels=levels, ordered=TRUE)
        
        kras.melt <- melt(kras.tbl, id.vars="sample")
        p.kras <- ggplot(data = kras.melt, aes(x = variable, y=sample))
        p.kras <- p.kras + geom_tile(aes(fill = value), show.legend = TRUE)
        p.kras <- p.kras + scale_fill_manual(values = c("black", "gray"), na.value = "white")
        p.kras <- p.kras + guides(fill = guide_legend(title="KRAS\nStatus", title.hjust=0.5))
        p.kras <- p.kras + theme(legend.text = element_text(size = fs))
        p.kras <- p.kras + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())        
        g.kras <- ggplotGrob(p.kras)
        g.legend.kras <- gtable::gtable_filter(g.kras, "guide-box")
        kras.label <- textGrob("KRAS", gp=gpar(fontsize=fs), just=just, rot=90)

##        kras.label <- gtable::gtable_filter(g.kras, "axis-b")
        
        g.kras <- gtable::gtable_filter(g.kras, "panel")

        if(!remove.kras.labels) {
            g.kras <- ggplotGrob(p.kras)
        } 
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- kras.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.kras
        indx <- indx + 1
    }

    g.legend.cms <- NULL
    has.cms <- TRUE
    if(has.cms) {
        cat("Has CMS\n")        
        cms.tbl <- data.frame(sample=clin.mask$sample, cms=clin.mask$cms_label)

        cms.tbl$cms <- as.character(cms.tbl$cms)

        cms.tbl$sample <- factor(cms.tbl$sample, levels=levels, ordered=TRUE)
        
        cms.melt <- melt(cms.tbl, id.vars="sample")
        p.cms <- ggplot(data = cms.melt, aes(x = variable, y=sample))
        p.cms <- p.cms + geom_tile(aes(fill = value), show.legend = TRUE)
        ##        p.cms <- p.cms + scale_fill_manual(values = c("black", "gray"), na.value = "white")
        p.cms <- p.cms + scale_fill_manual(values = c("black", "gray", "blue", "green"), na.value = "white")        
        p.cms <- p.cms + guides(fill = guide_legend(title="CMS\nLabel", title.hjust=0.5))
        p.cms <- p.cms + theme(legend.text = element_text(size = fs))
        p.cms <- p.cms + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
        g.cms <- ggplotGrob(p.cms)
        g.legend.cms <- gtable::gtable_filter(g.cms, "guide-box")

        cms.label <- textGrob("CMS", gp=gpar(fontsize=fs), just=just, rot=90)
##        cms.label <- gtable::gtable_filter(g.cms, "axis-b")        

        g.cms <- gtable::gtable_filter(g.cms, "panel")

        if(!remove.cms.labels) {
            g.cms <- ggplotGrob(p.cms)
        } 
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- cms.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.cms
        indx <- indx + 1
        
    }

    ## Introduce a spacer row between plots and legends and a second
    ## below the legends
    spacer.row <- rep(NA, num.cols)
    bottom.legend.row <- spacer.row
    bottom.legend.row[2] <- indx
    grobs[[indx]] <- g.legend.heat
    indx <- indx + 1

    if(num.annotation.legends > 0) {
        title.row <- c(title.row, NA)
    }
    
    layout <- title.row

    widths <- c(1, 5)
    annotation.width <- 0.25
    
    ## Add the annotation legends
    if(has.kras) {
        row <- c(plot.row, indx)
        layout <- rbind(layout, row)
        grobs[[indx]] <- g.legend.kras
        indx <- indx + 1
        widths <- c(widths, annotation.width)        
    }

    if(has.cms) {
        row <- c(plot.row, indx)
        layout <- rbind(layout, row)
        grobs[[indx]] <- g.legend.cms
        indx <- indx + 1
        widths <- c(widths, annotation.width)        
    }

    if(num.annotation.legends > 0) {
        widths <- c(widths, 2)
    }
    layout <- rbind(layout, spacer.row, bottom.legend.row, spacer.row)
    rownames(layout) <- NULL
    heights <- c(5,rep(10/max(1,num.annotation.legends), max(1,num.annotation.legends)),0.5,1,0.5)
    
##    lay <- cbind(c(NA,7,NA,NA,NA),c(2,8,NA,13,NA),c(3,9,NA,14,NA),c(4,10,NA,15,NA),c(5,11,NA,NA,NA),c(6,12,NA,NA,NA))
##    lay <- t(lay)
##    lay <- rbind(c(NA,2,3,4,5,6),c(7,8,9,10,11,12), c())

    print(widths)
    print(heights)
    print(layout)
    print(length(grobs))
##    grobs <- list(heat.label, msi.label, braf.label, tyms.label, mut.label, dendro.grob, g.heat, g.msi, g.braf, g.tyms, g.mut, g.legend.heat, g.legend.msi, g.legend.braf)
    ##    grob.table <- grid.arrange(grobs=grobs, layout_matrix=lay, heights=c(1,10,0.5,1,0.5), widths=c(1, 5, 0.5, 0.5, 0.5, 1))
    png(paste0(file, ".png"))
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    d <- dev.off()
    
}

