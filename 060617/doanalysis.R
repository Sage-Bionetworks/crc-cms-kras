

bee.dodge.width <- 0.75
bee.sz <- 3

## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
doAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.states=c("MT","WT"), only.do.kras.analysis = FALSE, plot.adjusted.pvalues = FALSE) {

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
    
    cat(paste("Number of rows analyzed: ", nrow(expr.m), "\n", sep=""))
    cat(paste("Number of columns analyzed: ", ncol(expr.m), "\n", sep=""))

    ## save(expr.m, file="expr.m.Rdata")
    ## save(bindeaGsets, file="bind.Rdata")
    
    ## Output genes in each set that are in data set
    for(i in 1:length(bindeaGsets)) {
        if(!(names(bindeaGsets)[i] %in% to.plot)) { next }
        flag <- bindeaGsets[[i]] %in% rownames(expr.m) 
        genes <- bindeaGsets[[i]][ flag ]
        missing.genes <- bindeaGsets[[i]][ !flag ]
        num.genes <- length(genes)
        num.missing.genes <- length(missing.genes)
        if(length(genes) == 0) { genes <- "None" }
        if(length(missing.genes) == 0) { missing.genes <- "None" }
        tot.genes <- length(bindeaGsets[[i]])
        genes <- paste(genes, collapse=" ")
        missing.genes <- paste(missing.genes, collapse=" ")        
        cat(paste(analysis.name, ": ", num.genes, " of ", tot.genes, " genes from data set analyzed in set ", names(bindeaGsets)[i], ": ", genes, "\n", sep=""))
        if(num.missing.genes > 0) {
            cat(paste(analysis.name, ": ", num.missing.genes, " of ", tot.genes, " genes missing from data set analyzed in set ", names(bindeaGsets)[i], ": ", missing.genes, "\n", sep=""))
        }
    }

    ## Make a heatmap showing the overlap (in terms of number of genes)
    ## in various gene sets -- bindea immune and CIRC
    populations <- c("B cells", "T cells", "Th1 cells", "Th2 cells", "Cytotoxic cells", "iDC", "Neutrophils")    
    gene.sets <- list()
    for(population in populations) {
        genes <- bindeaGsets[[population]][ bindeaGsets[[population]] %in% rownames(expr.m) ]
        gene.sets[[population]] <- genes
    }
    genes <- bindeaGsets[["CIRC"]][ bindeaGsets[["CIRC"]] %in% rownames(expr.m) ]
    gene.sets[["CIRC"]] <- genes

    m <- matrix(nrow=length(gene.sets), ncol=length(gene.sets), data=NA)
    for(i in 1:length(gene.sets)) {
        for(j in i:length(gene.sets)) {
            m[i,j] <- length(intersect(gene.sets[[i]], gene.sets[[j]]))
        }
    }
    rownames(m) <- names(gene.sets)
    colnames(m) <- names(gene.sets)

    m <- na.omit(melt(m))
    colnames(m) <- c("x", "y", "value")
    m$y <- factor(m$y, levels=rev(levels(m$y)))

    out.pdf <- paste("output/", "gene-set-heatmap-", analysis.name, ".pdf", sep="")
    pdf(out.pdf, useDingbats = FALSE)
    g <- ggplot(m, aes(x, y))
    g <- g + theme_bw() + xlab("") + ylab("")
    g <- g + geom_tile(aes(fill = value), color='white')
    nz.data <- m[m$value > 0,]
    g <- g + geom_text(data = nz.data, aes(x = x, y = y, label = value))
    g <- g + scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab', guide = guide_colorbar(title = "# Genes"))
    g <-g + theme(axis.text.x=element_text(angle=90),
                  axis.ticks=element_blank(),
                  axis.line=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_line(color='#eeeeee'))
    print(g)
    d <- dev.off()
    
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
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    kras.codon <- clin.m$kras_codon
    kras.code <- clin.m$kras_code   
    cms <- as.character(clin.m$cms_label)
    msi <- clin.m$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    site <- clin.m$site
    site.factor <- factor(site)
    neoantigens <- clin.m$neoantigens
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

    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    ##    cms2.last.levels <- c(cmses[cmses != "CMS2"], "CMS2")
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }
    ## Whatever is first will be used as the baseline;
    ## we don't need cms2 last, just make sure it's not first.
    ## Actually--used CMS4 as the baseline
    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]

    ## Analysis 5:  look at correlation of
    ## a. wnt vs circ
    ## b. tak vs circ
    ## c. wnt vs myc
    ## d. circ vs neoantigens
    sigs1 <- c("wnt", "tak1", "wnt", "CIRC")
    sigs2 <- c("CIRC", "CIRC", "myc.sig", "neoantigens")
    sig1.names <- c("wnt", "tak1", "wnt", "CIRC Enrichment Score")
    sig2.names <- c("CIRC", "CIRC", "myc.sig", "Log10 ( # Neoantigens + 1 ) ")
    for(i in 1:length(sigs1)) {
        sig1 <- sigs1[i]
        sig2 <- sigs2[i]
        ## par(mfrow=c(2,5))

        if(!(sig1 %in% rownames(esn)) || !(sig2 %in% rownames(esn))) { next }

        ## Here, don't facet by KRAS.
        cat(paste0(analysis.name, " ", sig1, " vs ", sig2, ": "))

        out.pdf <- paste("output/", "scatter-", analysis.name, "-", sig1, "-vs-", sig2, ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)
        df <- data.frame(x = esn[sig2,], y = esn[sig1,])

        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        
        g <- ggplot(data = df, aes(x = x, y = y, label = y), size = bee.sz)
        df <- df[!is.na(df$x),]
        df <- df[!is.na(df$y),]        
        if("status" %in% colnames(df)) {
            g <- g + geom_point(aes(colour = status))
            g <- g + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))            
        } else {
            g <- g + geom_point()
        }
        g <- g + ylab(sig1.names[i])
        g <- g + xlab(sig2.names[i])
        g <- g + theme(text = element_text(size = text.size))        
        if((sig1 == "CIRC") && (sig2 == "neoantigens")) {
            ## When comparing CIRC to neoantigens, don't plot the linear
            ## fit.  Instead, break into quadrants based on densities.
            ## (But y axis should be separated at CIRC enrichment = 0)
            ## Also, add MSI status.
            ## Calculate fisher's.
            d <- density(na.omit(esn[sig2,]))
            x.upper <- 2.75
            x.lower <- 2.0
            flag <- (d$x < x.upper) & (d$x > x.lower)
            x.flag <- d$x[flag]
            y.flag <- d$y[flag]
            sep.line <- x.flag[which(min(y.flag)==y.flag)[1]]
            res <- fisher.test(table(df$x < sep.line, df$y < 0))
            cat(paste0(analysis.name, ": Analyzing ", sig1, " vs ", sig2, " using ", res$method, ": p = ", res$p.value, " n = ", nrow(df), "\n"))
            
            g <- g + geom_hline(yintercept = 0, linetype = 'dashed')
            g <- g + geom_vline(xintercept = sep.line, linetype = 'dashed')
        } else {
            g <- g + stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse=TRUE)
            g <- g + geom_smooth(method = "lm", se = FALSE)
        }
        print(g)
        d <- dev.off()
    }
    cat("\n")
    
    ## Analysis 1: biomarker ~ kras
    num.biomarkers <- length(to.plot)

    kras.tbl <- data.frame(variable=to.plot, pval=rep(NA, num.biomarkers))
    rownames(kras.tbl) <- to.plot    
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

        lm.obj <- lm(as.formula(paste("expr", "KRAS", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        plot(lm.obj)
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

##    pval <- apply(esn[to.plot,,drop=F], 1, function(x) wilcox.test(x ~ kras)$p.value)
    pval <- unlist(lapply(to.plot, function(var) {
        x <- esn[var,,drop=T]
        non.na <- !is.na(x) & !is.na(kras)
        n.mt <- length(which(kras[non.na] == 1))
        n.wt <- length(which(kras[non.na] == 0))        
        res <- wilcox.test(x ~ kras)
        p <- res$p.value
        cat(paste0(analysis.name, ": Analyzing ", var, " vs KRAS using ", res$method, ": W = ", res$statistic, ", n_WT = ", n.wt, ", n_MT = ", n.mt, ", p = ", p, "\n"))
        p
    }))
    cat("\n")

    a <- p.adjust(pval,method="BH")
    b <- apply(esn[to.plot,,drop=F], 1, function(x) (mean(x[kras == 1]) - mean(x[kras == 0])) > 0)
    kras.tbl <- data.frame(variable=to.plot, pval=pval, apval=a, mutUp=b)
    rownames(kras.tbl) <- to.plot

    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    ## Create the individual expression plots for each gene/gene set for
    ## the kras analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
##        cat(paste("Computing ", st, " ~ ", "KRAS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=cms2.first.levels))
        
##        kras.label <- kras.code
        ##    df <- data.frame(expr=esn.st, KRAS=factor(kras.label), CMS=rep("ALL", length(esn.st)))
##        df <- data.frame(expr=esn.st, KRAS=factor(kras.codon), CMS=rep("ALL", length(esn.st)))
        
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], ".pdf", sep="")        
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))                
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))        
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        
        
        ## The rest of the ugliness below corresponds to the error bars.
        ## Use raw pvalue for univariate analysis
        pval <- as.numeric(kras.tbl[sti, pval.suffix])
        yoffset <- 0.01 * ( max(df$expr, na.rm=TRUE) - min(df$expr, na.rm=TRUE) )
        ## We only add one error bar.  Shift the ylimit to accomodate this shift.
        mask <- !is.na(df$KRAS) & (df$KRAS == kras.states[1])        
        mt.max <- max(df$expr[mask], na.rm=TRUE)

        mask <- !is.na(df$KRAS) & (df$KRAS == kras.states[2])                
        wt.max <- max(df$expr[mask], na.rm=TRUE)

        mt.wt.max <- max(mt.max, wt.max)
        
        ## Add the error bars between KRAS MT and KRAS WT expression values (within a facet)
        path.df <- data.frame(x = c(1, 1, 2, 2), y=c(mt.max + yoffset, mt.wt.max + 2 * yoffset, mt.wt.max + 2 * yoffset, wt.max + yoffset))
        p <- p + geom_path(data=path.df, aes(x, y))
        num.dy <- 1
        sz <- star.size
        text <- pval.to.text(pval)
        if(text == "n.s.") {
            num.dy <- 3
            sz <- ns.size
        }
        p <- p + annotate("text", x = 1.5, y = mt.wt.max + num.dy * yoffset, label = text, size = sz, vjust = 0)
        print(p)
        d <- dev.off()
    })

    write.table(file=paste0("output/", analysis.name, "-kras-tbl.xls"), kras.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    if(only.do.kras.analysis) {
        return(list(esn=esn, kras.tbl=kras.tbl))
    }
    
    ## Analysis 2: biomarker ~ cms
    cms.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {
        if(!(prefix %in% cms.levels)) { next }
        pval.col <- paste0(prefix, ".pval")
        cms.tbl[,pval.col] <- rep(NA, nrow(cms.tbl))
    }
    rownames(cms.tbl) <- to.plot
    
    cms.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {    
        if(!(prefix %in% cms.levels)) { next }        
        cms.str <- prefix
        flag <- !is.na(cms) & ( ( cms == "CMS2") | ( cms == cms.str ) )
        cms.flag <- cms[flag]
        ##        pval <- apply(esn[to.plot,flag,drop=F], 1, function(x) wilcox.test(x ~ cms.flag)$p.value)
        pval <- unlist(lapply(to.plot, function(var) {
            x <- esn[var,flag,drop=T]
            n.cms2 <- length(which(cms.flag == "CMS2"))
            n.cmso <- length(which(cms.flag == cms.str))
            res <- wilcox.test(x ~ cms.flag)
            p <- res$p.value
            cat(paste0(analysis.name, ": Analyzing ", var, " CMS2 vs ", cms.str, " using ", res$method, ": W = ", res$statistic, ", n_CMS2 = ", n.cms2, ", n_", cms.str, " = ", n.cmso, ", p = ", p, "\n"))
            p
        }))
        
        a <- p.adjust(pval,method="BH")
        b <- apply(esn[to.plot,flag,drop=F], 1, function(x) (mean(x[cms.flag == "CMS2"]) - mean(x[cms.flag == cms.str])) > 0)
        c.tmp <- cbind(pval=pval, apval=a, cms2Up=b)
        colnames(c.tmp) <- paste0(cms.str, ".", colnames(c.tmp))
        cms.tbl <- cbind(cms.tbl, c.tmp)
    }
    rownames(cms.tbl) <- to.plot
    cat("\n")
    
    ## Create the individual expression plots for each gene/gene set for
    ## the cms analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        if(exclude.no.lbl) {        
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-cms", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))
        ##        p <- p + geom_jitter()

        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))        
        if("status" %in% colnames(df)) {            
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        

        ## The rest of the ugliness below corresponds to the error bars.
        ## Let's include the unadjusted pvalues.
        yoffset <- 0.01 * ( max(df$expr, na.rm=TRUE) - min(df$expr, na.rm=TRUE) )
        ## panel.maxs will track the maximum value currently plotted in each facet.
        ## This will allow us to know where to position the error bar in each facet.
        panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(panel.maxs) <- cms.levels
        
        ## Draw the error bars from CMS2 to CMS1, 3, and 4.
        for(prefix in cms2.first.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }            
            cms.str <- prefix
            ## Use raw pvalue for univariate analysis
            ## pval.col <- paste0(cms.str, ".apval")
            pval.col <- paste0(cms.str, ".", pval.suffix)
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

            panel.maxs["CMS2"] <- both.max + 2 * 6 * yoffset
            panel.maxs[cms.str] <- both.max + 2 * 6 * yoffset            
            
            path.df <- data.frame(x = xs, y = ys)
            p <- p + geom_path(data=path.df, aes(x, y))
            text <- pval.to.text(pval)
            num.dy <- 1
            ## "-cms"
            sz <- star.size
            if(text == "n.s.") {
                num.dy <- 3
                sz <- ns.size
            }
            p <- p + annotate("text", x = 0.5 * (cms2indx + cmsindx), y = both.max + num.dy * yoffset, label = text, size = sz, vjust = 0)
        }
        print(p)
        d <- dev.off()

    })

    write.table(file=paste0("output/", analysis.name, "-cms-tbl.xls"), cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    ## End analysis 2

    ## Analysis 2b: biomarker ~ cms mt (only)
    cms.mt.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {
        if(!(prefix %in% cms.levels)) { next }
        pval.col <- paste0(prefix, ".pval")
        cms.mt.tbl[,pval.col] <- rep(NA, nrow(cms.mt.tbl))
    }
    rownames(cms.mt.tbl) <- to.plot
    
    for(sti in 1:num.biomarkers) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(paste0(cms, " MT"), levels=paste0(cms2.first.levels, " MT")))
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL MT",]                        
        }
        ## Restrict to KRAS MT only
        df <- df[!is.na(df$KRAS) & (kras.label == "MT"),]
        
        do.fit.and.qc(df, st, "CMS", analysis.name, "CMS-mt")
        for(prefix in cms2.first.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }            
            cms.str <- prefix
            pval.col <- paste0(cms.str, ".pval")
            pval.flag <- grepl(x=coeffs.2$coefficient, cms.str)
            ## Set this below using wilcox
            ##            cms.mt.tbl[sti,pval.col] <- coeffs.2[pval.flag, 5]
        }

        
    }
    cat("\n")

    cms.mt.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {    
        if(!(prefix %in% cms.levels)) { next }        
        cms.str <- prefix
        flag <- !is.na(kras.label) & !is.na(cms) & ( kras.label == "MT" ) & ( ( cms == "CMS2") | ( cms == cms.str ) )
        cms.flag <- cms[flag]
        ##        pval <- apply(esn[to.plot,flag,drop=F], 1, function(x) wilcox.test(x ~ cms.flag)$p.value)
        pval <- unlist(lapply(to.plot, function(var) {
            x <- esn[var,flag,drop=T]
            n.cms2 <- length(which(cms.flag == "CMS2"))
            n.cmso <- length(which(cms.flag == cms.str))
            res <- wilcox.test(x ~ cms.flag)
            p <- res$p.value
             cat(paste0(analysis.name, ": Analyzing ", var, " CMS2 vs ", cms.str, " (KRAS MT only) using ", res$method, ": W = ", res$statistic, ", n_CMS2 = ", n.cms2, ", n_", cms.str, " = ", n.cmso, ", p = ", p, "\n"))
            p
        }))
        
        a <- p.adjust(pval,method="BH")
        b <- apply(esn[to.plot,flag,drop=F], 1, function(x) (mean(x[cms.flag == "CMS2"], na.rm=TRUE) - mean(x[cms.flag == cms.str], na.rm=TRUE)) > 0)
        c.tmp <- cbind(pval=pval, apval=a, cms2Up=b)
        colnames(c.tmp) <- paste0(cms.str, ".", colnames(c.tmp))
        cms.mt.tbl <- cbind(cms.mt.tbl, c.tmp)
    }
    rownames(cms.mt.tbl) <- to.plot
    cat("\n")
    
    ## Create the individual expression plots for each gene/gene set for
    ## the cms analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(paste0(cms, " MT"), levels=paste0(unique.cms, " MT")))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        if(exclude.no.lbl) {        
            df <- df[df$CMS != "NOLBL MT",]                        
        }
        ## Restrict to KRAS MT
        df <- df[kras.label == "MT", ]
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-cms-mt", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))        
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        
        
        ## The rest of the ugliness below corresponds to the error bars.
        ## Let's include the unadjusted pvalues.
        yoffset <- 0.01 * ( max(df$expr, na.rm=TRUE) - min(df$expr, na.rm=TRUE) )

        ## panel.maxs will track the maximum value currently plotted in each facet.
        ## This will allow us to know where to position the error bar in each facet.
        panel.maxs <- sapply(paste0(cms.levels, " MT"), function(cmsLbl) {
            mask <- df$CMS == cmsLbl
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(panel.maxs) <- cms.levels
        
        ## Draw the error bars from CMS2 to CMS1, 3, and 4.
        for(prefix in cms2.first.levels[-1]) {
            if(!(prefix %in% cms.levels)) { next }            
            cms.str <- prefix
            ## Use raw pvalue for univariate analysis
            ## pval.col <- paste0(cms.str, ".apval")
            pval.col <- paste0(cms.str, ".", pval.suffix)
            pval <- as.numeric(cms.mt.tbl[sti, pval.col])
        
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

            panel.maxs["CMS2"] <- both.max + 2 * 6 * yoffset
            panel.maxs[cms.str] <- both.max + 2 * 6 * yoffset            
            
            path.df <- data.frame(x = xs, y = ys)
            p <- p + geom_path(data=path.df, aes(x, y))
            num.dy <- 1
            sz <- star.size
            text <- pval.to.text(pval)
            if(text == "n.s.") {
                num.dy <- 3
                sz <- ns.size
            }            
            p <- p + annotate("text", x = 0.5 * (cms2indx + cmsindx), y = both.max + num.dy * yoffset, label = text, size = sz, vjust = 0)
        }
        print(p)
        d <- dev.off()

    })

    write.table(file=paste0("output/", analysis.name, "-cms-mt-tbl.xls"), cms.mt.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    ## End analysis 2b

    ## Analysis 3: biomarker ~ cms1-kraswt + cms1-krasmt + cms2-kraswt + ... (with ref = cms2:krasMT)
    kras.cms.tbl <- data.frame(variable=to.plot)
    rownames(kras.cms.tbl) <- to.plot        
    base <- paste0(".vs.CMS2", kras.states[1])
    for(prefix in cms2.first.levels) {
        if(!(prefix %in% cms.levels)) { next }        
        for(genotype in kras.states) {
            if( ( prefix != "CMS2" ) || ( genotype != kras.states[1] ) ) {
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
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        df$cms.kras <- apply(cbind(as.character(df$CMS), as.character(df$KRAS)), 1, function(row) paste0(row[1], row[2]))
        df$cms.kras <- factor(df$cms.kras)
        df$cms.kras <- relevel(df$cms.kras, ref = paste0("CMS2", kras.states[1]))
        
        lm.obj <- lm(as.formula(paste("expr", "cms.kras", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-cms-kras", "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        plot(lm.obj)
        d <- dev.off()

        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        base <- paste0(".vs.CMS2", kras.states[1])
        for(prefix in cms2.first.levels) {        
            if(!(prefix %in% cms.levels)) { next }            
            for(genotype in kras.states) {
                if( ( prefix != "CMS2" ) || ( genotype != kras.states[1] ) ) {
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
    base <- paste0(".vs.CMS2", kras.states[1])
    for(prefix in cms2.first.levels) {            
        if(!(prefix %in% cms.levels)) { next }        
        for(genotype in kras.states) {
            if( ( prefix != "CMS2" ) || ( genotype != kras.states[1] ) ) {
                pval.col <- paste0(prefix, genotype, base, ".pval")
                flag1 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == "CMS2" ) & ( kras.label == kras.states[1] ) ) )
                flag2 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == prefix ) & ( kras.label == genotype ) ) )
                pvals <- rep(1, length(to.plot))
                
                if((length(which(flag1)) >= 2) && (length(which(flag2)) >= 2)) {
                    ## pvals <- apply(esn[to.plot,,drop=F], 1, function(x) { wilcox.test(x = x[flag1], y = x[flag2])$p.value })
                    pvals <- unlist(lapply(to.plot, function(var) {
                        x <- esn[var,,drop=T]
                        n.cms2mt <- length(which(flag1))
                        n.other <- length(which(flag2))
                        res <- wilcox.test(x = x[flag1], y = x[flag2])
                        p <- res$p.value
                        cat(paste0(analysis.name, ": Analyzing ", var, " CMS2-", kras.states[1], " vs ", prefix, "-", genotype, " using ", res$method, ": W = ", res$statistic, ", n_CMS2-", kras.states[1], " = ", n.cms2mt, ", n_", prefix, "-", genotype, " = ", n.other, ", p = ", p, "\n"))
                        p
                    }))
                    
                }
                kras.cms.tbl[,pval.col] <- pvals
                padj.col <- paste0(prefix, genotype, base, ".apval")
                kras.cms.tbl[,padj.col] <- p.adjust(kras.cms.tbl[,pval.col], method="BH")
            }
        }
    }
    cat("\n")
    
    ## Create the individual expression plots for each gene/gene set for
    ## the cms analysis
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))        
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits, na.rm=TRUE) - min(ylimits, na.rm=TRUE) )
        ## We will add 2 error bars for each cms label excluding cms2
        ## and one more between cms2 mt and wt.
        num.error.bars <- 2 * ( length(unique(df$CMS)) - 1 ) + 1
        ylimits[2] <- ylimits[2] + 2 * 7 * num.error.bars * yoffset
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
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(wt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[2]))
        
        panel.maxs <- c(mt.panel.maxs, wt.panel.maxs)

        ## Add the error bars between KRAS MT and KRAS WT expression values (within a facet)
        for(prefix in sort(cms2.last.levels)) {
            if(!(prefix %in% cms.levels)) { next }            
            for(genotype in kras.states) {
                if( ( prefix != "CMS2" ) || ( genotype != kras.states[1] ) ) {
                    cmsLbl2 <- "CMS2"
                    kras.state2 <- kras.states[1]
                    cmsIndx2 <- which(sort(cms2.last.levels) == cmsLbl2)[1]
                    cmsLbl1 <- prefix
                    kras.state1 <- genotype
                    cmsIndx1 <- which(sort(cms2.last.levels) == cmsLbl1)[1]                    
                    base <- paste0(".vs.CMS2", kras.states[1])
                    cmpLbl <- paste0(prefix, genotype, base)
                    ## "-kras-cms"
                    tmp <- draw.err.bar(st, kras.cms.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix, yoffset)
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
    kras.col <- paste0("kras", kras.states[1])
    cols <- c(kras.col)
    for(prefix in cms2.last.levels[-1]) {
        if(exclude.no.lbl) {
            if(prefix == "NOLBL") { next }
        }
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
        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)))

        
        formula <- "KRAS + CMS + CMS:KRAS"
        formula.no.interaction <- "KRAS + CMS"        
        formula.name<- "kras-cms-interaction"
        ## Now exclude NOLBL
        df.lm <- df
        if(exclude.no.lbl) {
            df.lm <- df[df$CMS != "NOLBL",]
        }
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        lm.obj.no.interaction <- lm(as.formula(paste("expr", formula.no.interaction, sep=" ~ ")), data=df.lm)        
        ## library(lmPerm)
        ## lm.obj <- lmp(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        ## lm.obj <- glm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm, family=gaussian(link="log"))
        lm.sum <- summary(lm.obj)
##        print(lm.sum)
        capture.output(lm.sum, file = sum.file)

        anova.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-anova.tsv")        
        anova.out <- anova(lm.obj.no.interaction, lm.obj)              
##        anova.df <- as.data.frame(anova.out)
##        write.table(file=anova.file, anova.df, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
        capture.output(anova.out, file = anova.file)
        
        flag <- unlist(apply(df.lm, 1, function(row) any(is.na(row))))
        bf.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-bf.tsv")
        bf <- generalTestBF(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm[!flag, ])
##        save(bf, file="bf.Rd")
        capture.output(bf/bf["KRAS + CMS"], file=bf.file)
        
        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        qqnorm(df.lm$expr)
        qqline(df.lm$expr)
        plot(lm.obj)
        d <- dev.off()
        lm.sum <- summary(lm.obj)
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern=paste0("KRAS", kras.states[1])) & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
##        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", st, "\n"))                
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
##            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", st, "\n"))                
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
##            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", st, "\n"))                
            interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        ## Repeat the analysis for transformed data
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.obj.no.interaction <- lm(as.formula(paste("expr", formula.no.interaction, sep=" ~ ")), data=df.transformed)        
        lm.sum <- summary(lm.obj)
##        print(lm.sum)
        capture.output(lm.sum, file = sum.file)

        anova.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-anova.tsv")        
        anova.out <- anova(lm.obj.no.interaction, lm.obj)              
##        anova.df <- as.data.frame(anova.out)
##        write.table(file=anova.file, anova.df, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
        capture.output(anova.out, file = anova.file)
        
        flag <- unlist(apply(df.lm, 1, function(row) any(is.na(row))))
        bf.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-bf.tsv")
        bf <- generalTestBF(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm[!flag,])
        capture.output(bf/bf["KRAS + CMS"], file=bf.file)        
        
        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        qqnorm(df.lm$expr)
        qqline(df.lm$expr)
        plot(lm.obj)
        d <- dev.off()
        lm.sum <- summary(lm.obj)
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern=paste0("KRAS", kras.states[1])) & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
##        cat(paste0("Attempting to write transformed ", coeffs.2[flag, 5], " to col ", col, " for ", st, "\n"))                
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
##            cat(paste0("Attempting to write transformed ", coeffs.2[flag, 5], " to col ", col, " for ", st, "\n"))                
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
##            cat(paste0("Attempting to write transformed ", coeffs.2[flag, 5], " to col ", col, " for ", st, "\n"))                
            transformed.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
##        cat(paste0("Attempting to padj to col ", apcol, "\n"))
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
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        ##        df <- data.frame(expr=esn.st, KRAS=factor(kras.code), CMS=factor(cms, levels=unique.cms))
##        df <- data.frame(expr=esn.st, KRAS=factor(kras.codon), CMS=factor(cms, levels=unique.cms))        


        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]
        }
        df.transformed <- df
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms-interaction", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))        
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits, na.rm=TRUE) - min(ylimits, na.rm=TRUE) )
        ## And we will add 1 error bars.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 1       
        ylimits[2] <- ylimits[2] + 2 * 6 * num.error.bars * yoffset
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
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$expr[mask], na.rm=TRUE)
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
##            print(cmpLbl)
            tmp <- draw.err.bar(st, interaction.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix, yoffset)
            g <- tmp[["g"]]
            panel.maxs <- tmp[["panel.maxs"]]
        }
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()

        ## Repeat for transformed data
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms-interaction-transformed", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df.transformed, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))        
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {            
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        
        
        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits, na.rm=TRUE) - min(ylimits, na.rm=TRUE) )
        ## And we will add 1 error bars.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 1       
        ylimits[2] <- ylimits[2] + 2 * 6 * num.error.bars * yoffset
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
            m <- max(df.transformed$expr[mask], na.rm=TRUE)
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df.transformed$CMS == cmsLbl & df.transformed$KRAS == kras.states[2]
            m <- max(df.transformed$expr[mask], na.rm=TRUE)
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

            tmp <- draw.err.bar(st, transformed.interaction.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix, yoffset)
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

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)))

        formula <- "KRAS + CMS"
        formula.name <- "kras-cms-no-interaction"
        ## Now exclude NOLBL
        df.lm <- df
        if(exclude.no.lbl) {
            df.lm <- df[df$CMS != "NOLBL",]
        }
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        anova.df <- as.data.frame(anova(lm.obj))
        anova.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-anova.tsv")        
        write.table(file=anova.file, anova.df, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

        
        
        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        plot(lm.obj)
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern=paste0("KRAS", kras.states[1])) & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
##        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
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
            
##            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            no.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        ## Repeat above for transformed data
        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.transformed)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-transformed-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        anova.df <- as.data.frame(anova(lm.obj))
        anova.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-anova.tsv")        
        write.table(file=anova.file, anova.df, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
        
        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        plot(lm.obj)
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        flag <- grepl(x=coeffs.2$coefficient, pattern=paste0("KRAS", kras.states[1])) & !grepl(x=coeffs.2$coefficient, pattern=":")
        col <- paste0(kras.col, ".pval")
##        cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
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
            
##            cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))                
            transformed.no.interaction.tbl[sti, col] <- coeffs.2[flag, 5]
        }

        
    }

    for(prefix in cols) {
        pcol <- paste0(prefix, ".pval")
        apcol <- paste0(prefix, ".apval")
##        cat(paste0("Attempting to padj to col ", apcol, "\n"))
        transformed.no.interaction.tbl[,apcol] <- p.adjust(transformed.no.interaction.tbl[,pcol], method="BH")
    }

    write.table(file=paste0("output/", analysis.name, "-kras-cms-no-interaction-tbl.xls"), no.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-kras-cms-transformed-no-interaction-tbl.xls"), transformed.no.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)    
    ## End Analysis 5

    ## Analysis 6: biomarker ~ cms + kras (with error bars from kras-mt vs kras-wt within cms)

    mt.vs.wt.tbl <- do.call("cbind",lapply(cms.levels, function(cmsLbl){
        mask <- cmsLbl ==  cms
        pval <- rep(1, length(to.plot))
        
        if((length(which(mask)) >= 2) && (length(unique(kras[mask])) == 2)) {
            ## pval <- apply(esn[to.plot,mask,drop=F], 1, function(x) wilcox.test(x ~ kras[mask])$p.value)
            pval <- unlist(lapply(to.plot, function(var) {
                x <- esn[var,mask,drop=T]
                n.mt <- length(which(kras[mask] == 1))
                n.wt <- length(which(kras[mask] == 0))                
                res <- wilcox.test(x ~ kras[mask])
                p <- res$p.value
                cat(paste0(analysis.name, ": Analyzing ", var, " ", cmsLbl, " MT vs WT using ", res$method, ": W = ", res$statistic, ", n_MT = ", n.mt, ", n_WT = ", n.wt, ", p = ", p, "\n"))
                p
            }))
        }
        a <- p.adjust(pval,method="BH")
        b <- apply(esn[to.plot,mask,drop=F], 1, function(x) (mean(x[kras[mask] == 1], na.rm=TRUE) - mean(x[kras[mask] == 0], na.rm=TRUE)) > 0)
        ret <- cbind(pval=pval,apval=a, mutUp=b)
        ret
    }))
    cat("\n")
    
    mt.vs.wt.tbl <- as.data.frame(mt.vs.wt.tbl)
    mt.vs.wt.tbl$variable <- to.plot
    colnames(mt.vs.wt.tbl) <- do.call(c, lapply(paste0(cms.levels, ".", kras.states[1], ".vs.", kras.states[2]), function(x) paste(x, c(".pval", ".apval", ".mutUp"), sep="")))

    ## Plot the mt-vs-wt wilcox error bars
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms-faceted", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))        
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        

        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        ## The rest of the ugliness below corresponds to the error bars.

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(ylimits, na.rm=TRUE) - min(ylimits, na.rm=TRUE) )
        ## And we will only add 1 error bar.  So adjust the y axis
        ## to accommodate this shift.
        num.error.bars <- 1
        ylimits[2] <- ylimits[2] + 2 * 6 * num.error.bars * yoffset
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
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(mt.panel.maxs) <- sapply(cms.levels, function(x) paste0(x, "-", kras.states[1]))

        wt.panel.maxs <- sapply(cms.levels, function(cmsLbl) {
            mask <- df$CMS == cmsLbl & df$KRAS == kras.states[2]
            m <- max(df$expr[mask], na.rm=TRUE)
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
            cmpLbl <- paste0(cmsLbl, ".", kras.states[1], ".vs.", kras.states[2])

            ## "-kras-cms-faceted"
            tmp <- draw.err.bar(st, mt.vs.wt.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix, yoffset)
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

        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)), site=site.factor)
        formula <- "KRAS + CMS + site"
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            formula <- paste0(formula, " + status")
            df$status <- msi.factor
        }
        if(!all(is.na(neoantigens))) {
            formula <- paste0(formula, " + neoantigens")
            df$neoantigens <- neoantigens
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        flag <- unlist(apply(df, 1, function(row) any(is.na(row))))
        df <- df[!flag,]

        formula.name <- "full-model"
        ## Now exclude NOLBL
        df.lm <- df
        if(exclude.no.lbl) {
            df.lm <- df[df$CMS != "NOLBL",]
        }
        df.transformed <- df.lm
        df.transformed$expr <- inverse.norm.transform(df.transformed$expr)

        lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file, useDingbats = FALSE)
        p <- suppressWarnings(forest_model(lm.obj))
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file, useDingbats = FALSE)
        plot(lm.obj)
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        for(prefix in cols) {
            col <- paste0(prefix, ".pval")
            flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
            if(length(which(flag)) != 1) {
##                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
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
        pdf(diag.file, useDingbats = FALSE)
        plot(lm.obj)
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file, useDingbats = FALSE)
        p <- suppressWarnings(forest_model(lm.obj))
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
##                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
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
##        cat(paste0("Attempting to padj to col ", apcol, "\n"))
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
        formulas <- c(formulas, "status")
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

        df <- data.frame(CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)), site=site.factor)
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
                df$status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }

        for(lvl in cols) {
            cat(paste("Computing ", lvl, " vs ", base, "\n", sep=""))
            vec <- as.vector(df[,formula,drop=T])
            flag1 <- !is.na(vec) & ( vec == base)
            flag2 <- !is.na(vec) & ( vec == lvl)
            
            ## pvals <- apply(esn[to.plot,,drop=F], 1, function(x) wilcox.test(x = x[flag1], y = x[flag2])$p.value)
            pvals <- unlist(lapply(to.plot, function(var) {
                x <- esn[var,,drop=T]
                n.lvl <- length(which(flag2))
                n.base <- length(which(flag1))
                res <- wilcox.test(x = x[flag1], y = x[flag2])
                p <- res$p.value
                cat(paste0(analysis.name, ": Analyzing ", var, " ", formula, ": ", lvl, " vs ", base, res$method, ": W = ", res$statistic, ", n_", lvl, " = ", n.lvl, ", n_", base, " = ", n.base, ", p = ", p, "\n"))
                p
            }))
            
            pval.col <- paste0(lvl, ".pval")
            padj.col <- paste0(lvl, ".apval")
            
            formula.tbls[[i]][,pval.col] <- pvals
            formula.tbls[[i]][,padj.col] <- p.adjust(formula.tbls[[i]][,pval.col], method="BH")
        }
        cat("\n")
            
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
            
            df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)), site=site.factor)
            if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {            
                df$status <- msi.factor
            }
            if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
                df$site <- site.factor
            }
            
            flag <- unlist(apply(df[,formula,drop=F], 1, function(row) any(is.na(row))))
            df <- df[!flag,]

            ## Create a pdf for this gene/gene set (with facets)
            out.pdf <- paste("output/", st, "-", analysis.name, "-", formula, ".pdf", sep="")
            pdf(out.pdf, useDingbats = FALSE)

            ## Create the plot
            formula.indx <- which(colnames(df) == formula)[1]
            expr.indx <- which(colnames(df) == "expr")[1]
            formula.col <- colnames(df)[formula.indx]
            expr.col <- "expr"
            p <- ggplot(data=df, aes_string(x = formula.col, y = expr.col))
            p <- p + theme(text = element_text(size = text.size))            
            p <- p + ylab(ylab)

            ## Create a box plot where the x axis is CMS ...
            p <- p + geom_boxplot(aes_string(fill = formula.col), outlier.shape = NA)
            ##            p <- p + geom_jitter()
            ## Do not annotate MSI/MSS when we are comparing expr to MSI status
            p <- p + guides(fill = guide_legend(order = 1))
            nxt.order <- 2
            if(("status" %in% colnames(df)) && ("site" %in% colnames(df)) && (formula != "site") && (formula != "status")) {
                p <- p + geom_beeswarm(aes(colour = status, shape = site), dodge.width = bee.dodge.width, size = bee.sz)
                p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
                p <- p + guides(colour = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                
                p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
                p <- p + guides(shape = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1
            } else if(("status" %in% colnames(df)) && (formula != "status")) {
                p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
                p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
                p <- p + guides(colour = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                                
            } else if(("site" %in% colnames(df)) && (formula != "site")) {
                p <- p + geom_beeswarm(aes(shape = site), dodge.width = bee.dodge.width, size = bee.sz)
                p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
                p <- p + guides(shape = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                                
            } else {
                p <- p + geom_beeswarm(size = bee.sz)
            }
            
            ## The rest of the ugliness below corresponds to the error bars.
            ## Let's include the unadjusted pvalues.
            yoffset <- 0.01 * ( max(df$expr, na.rm=TRUE) - min(df$expr, na.rm=TRUE) )
            
            ## panel.maxs will track the maximum value currently plotted in each facet.
            ## This will allow us to know where to position the error bar in each facet.
            lvls <- levels(df[,formula,drop=T])
            panel.maxs <- sapply(lvls, function(lvl) {
                mask <- df[,formula,drop=T] == lvl
                m <- max(df$expr[mask], na.rm=TRUE)
                m
            })
            names(panel.maxs) <- lvls

            ## BAD here
            
            ## Draw the error bars from the first level to each of the others
            for(prefix in cols) {
                str <- prefix
                ## Use raw pvalue for univariate analysis                        
                ## pval.col <- paste0(str, ".apval")
                pval.col <- paste0(str, ".", pval.suffix)
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

                panel.maxs[base] <- both.max + 2 * 6 * yoffset
                panel.maxs[str] <- both.max + 2 * 6 * yoffset            
                
                path.df <- data.frame(x = xs, y = ys)
                p <- p + geom_path(data=path.df, aes(x, y))
                ## "-site"
                num.dy <- 1
                sz <- star.size
                text <- pval.to.text(pval)
                if(text == "n.s.") {
                    num.dy <- 3
                    sz <- ns.size
                }                
                p <- p + annotate("text", x = 0.5 * (base.indx + col.indx), y = both.max + num.dy * yoffset, label = text, size = sz, vjust = 0)
            }
            print(p)
            d <- dev.off()
            
        })
        
    }

    ## End Analysis 8

    ret.list <- list(esn=esn, kras.tbl=kras.tbl, cms.tbl=cms.tbl, kras.cms.tbl=kras.cms.tbl, kras.cms.interaction.tbl=interaction.tbl, kras.cms.no.interaction.tbl=no.interaction.tbl, kras.mt.vs.wt.tbl=mt.vs.wt.tbl, site.tbl=formula.tbls[[1]])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        ret.list[["msi.tbl"]] <- formula.tbls[[2]]
    }
    return(ret.list)

    
}

doUnivariateKrasAnalysis <- function(expr, clin, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT")) {

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
    
    cat(paste("Number of rows analyzed: ", nrow(expr.m), "\n", sep=""))
    cat(paste("Number of columns analyzed: ", ncol(expr.m), "\n", sep=""))


    if(length(to.plot) == 0) { to.plot <- rownames(expr.m) }
    expr.m <- expr.m[rownames(expr.m) %in% to.plot,]
    to.plot <- rownames(expr.m)

    ## Exclude any samples that do not have KRAS mutation status annotated
    kras <- clin.m[, kras.status.field]

    ## Switch this to use limma
    ## Filtered to have cpm > 1
##    ./121816/run.R:suppressPackageStartupMessages(library("limma"))
    
##cms.kras.limma <- cms.kras
##flag <- cms.kras.limma == "CMS2-1"
##cms.kras.limma[flag] <- "CMS2-MT"
##cms.kras.limma[!flag] <- "Other"
##design <- model.matrix(~ cms.kras.limma)

    ## Should these be log transformed?  grcma is used in limma user's guide.

    ##fit <- lmFit(expr.m, design)
##efit <- eBayes(fit, trend=TRUE)
    ##efit$p.value[is,2]

    ## Filtering described here: See for example Case studies 15.3 or 15.4 in the limma User's Guide.

##tt <- topTable(efit, number=Inf, coef="cms.kras.limmaOther", adjust.method="BH")

    pval <- ldply(to.plot, .parallel = TRUE, .fun = function(var) {
        x <- expr.m[var,,drop=T]
        non.na <- !is.na(x) & !is.na(kras)
        mt.flag <- non.na & (kras == 1)
        effect <- mean(x[mt.flag], na.rm=TRUE) - mean(x[!mt.flag], na.rm=TRUE)
        n.mt <- length(which(kras[non.na] == 1))
        n.wt <- length(which(kras[non.na] == 0))        
        res <- wilcox.test(x ~ kras)
        p <- res$p.value
        vec <- c(var, p, effect)
        names(vec) <- c("gene", "p", "effect")
        return(vec)
    })
    pval

}

## pvals <- doUnivariateKrasAnalysis(tcga_expr, clin, to.plot=intersect(rownames(tcga_expr), rownames(kfsyscc_expr)))
