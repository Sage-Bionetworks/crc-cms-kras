

bee.dodge.width <- 0.75
bee.sz <- 3

## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
doAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.states=c("MT","WT"), only.do.kras.analysis = FALSE, plot.adjusted.pvalues = FALSE, main = NULL, do.gene.set.heatmap = TRUE, do.gsva = TRUE) {

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
    if(do.gene.set.heatmap) {
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
        if(!is.null(main)) { g <- g + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }    
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
    }
    
    ## Perform GSVA analysis
    es <- expr.m
    if(do.gsva) {
        es <- gsva(expr.m, bindeaGsets,parallel.sz=num.processes,verbose=TRUE)$es.obs
        cat("Done with gsva.\n")
    }
    
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
        if(!is.null(main)) { g <- g + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
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
        
        cat(paste("Computing ", st, " ~ ", "KRAS", "\n", sep=""))

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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
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
    col <- "f.pval"
    full.model.tbl[,col] <- rep(NA, nrow(full.model.tbl))    

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

        ## Also output the overall pvalue
        f <- lm.sum$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        full.model.tbl[sti, "f.pval"] <- p
        
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

        ## Also output the overall pvalue
        f <- lm.sum$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        transformed.full.model.tbl[sti, "f.pval"] <- p
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
            if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                    
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
    num.annotation.legends <- 0
    has.kras <- TRUE
    if(has.kras) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    has.cms <- TRUE
    has.cms <- any(!is.na(clin.mask$cms_label)) && (length(unique(clin.mask$cms_label)) > 1)
    if(has.cms) {
        num.annotation.legends <- num.annotation.legends + 1
    }
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
    heat.label <- textGrob("Genes in RAS Dependency Signature", gp=gpar(fontsize=fs), just=heat.just)

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
    pdf(paste0(file, ".pdf"))
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    d <- dev.off()
    
}

cluster.annotated.matrix <- function(mat, clin, file) {

    ## Based on:
    ##     http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    ##     http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically/17371177#17371177

    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("gtable"))
    suppressPackageStartupMessages(library("ggdendro"))    

    clin.orig <- clin
    tmp <- clin
    tmp <- tmp[!is.na(tmp$cms_label) & (tmp$cms_label %in% c("CMS1", "CMS2", "CMS3", "CMS4")),]
##    tmp$sample <- gsub(tmp$sample, pattern="-", replacement=".")

    mat.tmp <- mat
##    mat.tmp <- mat.tmp[mat.tmp$P.value < 0.1,]
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(tmp$sample, colnames(mat.tmp))
    clin.mask <- tmp[!is.na(idxs),]
    mat.data <- t(mat.tmp[,na.omit(idxs)])
    print(head(idxs))
    print(head(mat.tmp))
    print(head(mat.data))

    ## We will annotate with RAS status and CMS
    num.annotation.legends <- 0
    has.kras <- TRUE
    if(has.kras) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    has.cms <- TRUE
    has.cms <- any(!is.na(clin.mask$cms_label)) && (length(unique(clin.mask$cms_label)) > 1)
    if(has.cms) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    num.cols <- num.annotation.legends

    ## Add to column count for data heatmap and dendrogram
    num.cols <- num.cols + 2
    ## And one for the legends
    if(num.annotation.legends > 0) {
        num.cols <- num.cols + 1
    }
    
    method <- "ward.D"
    method <- "complete"
    dd.row <- as.dendrogram(hclust(dist(mat.data), method=method))
    row.ord <- order.dendrogram(dd.row)

    hc.col <- hclust(dist(t(mat.data)), method=method)
    dd.col <- as.dendrogram(hc.col)
    col.ord <- order.dendrogram(dd.col)

    ## Re-order the samples in the data matrix and clin based on the clustering
    mat.data <- (mat.data[row.ord, ])    
    clin.mask <- clin.mask[row.ord, ]

    remove.labels <- TRUE
    remove.kras.labels <- TRUE
    remove.cms.labels <- TRUE

    if(!remove.labels || !remove.kras.labels || !remove.cms.labels) {
        ## Set all put every 20th to just a number so we can read
        flag <- rep_len(c(FALSE,rep(TRUE,9)), nrow(mat.mask))
        rownames(mat.data)[flag] <- as.character((1:nrow(mat.data))[flag])        
        clin.mask$sample[flag] <- as.character((1:length(clin.mask$sample))[flag])
    }

    ## NB: we are not scaling the columns/rows of the data

    ## Scale the rows/genes of the matrix
    ## NB: we are doing it here--after clustering!  So it will only effect
    ## the visuals.
    ## Scale/z-score the rows/genes of the matrix
    ## mat.data <- t(scale(t(mat.data)))
    
    ## It is important that the samples be factors to ensure consistency
    ## across the plots
    xx <- mat.data
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
    heat.label <- textGrob("Genes in RAS Dependency Signature", gp=gpar(fontsize=fs), just=heat.just)

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
    pdf(paste0(file, ".pdf"))
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    d <- dev.off()
    
}

## Generate permutations
## For each permutation
##    For each gene set
##        Calculate all pvalues and sign
## For each gene set
##    Perform test
##    Calculate FDR
##    Output table
## Generate figure by plotting table

doGSVA <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", exclude.nras.mutants = TRUE) {

    clin.tmp <- clin
     
    ## Possibly exclude braf and nras mutants
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

    ## Perform GSVA analysis
    es <- gsva(expr.m, bindeaGsets,parallel.sz=num.processes,verbose=TRUE)$es.obs
    cat("Done with gsva.\n")
    
    ## Make sure that the genes/gene sets to plot are actually in the expression set.
    all.genes.and.sets <- c(rownames(expr.m), rownames(es))
    flag <- to.plot %in% all.genes.and.sets
    to.plot <- to.plot[flag]
    tmp <- rbind(es, expr.m[to.plot[to.plot %in% rownames(expr.m)],,drop=F])
    tmp
    list(es = tmp, clin = clin.m)
}


doKrasCMSAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 1000, calc.fdr.with.bh = FALSE, calc.fdr = TRUE, seed = 1234) {

    set.seed(seed)

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    cms <- as.character(clin$cms_label)    
    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }
    
    permuted.tests <- NULL
    if( (num.permutations > 0) && (calc.fdr.with.bh == FALSE) ) {
        ## Run permutated tests of expr ~ kras + cms in which we compare everything to CMS2 MT
        permuted.tests <- ldply(1:num.permutations, .parallel = TRUE,
                                .fun = function(i) {
                                    ## Permute the kras and cms labels
                                    clin.perm <- clin
                                    permutation <- sample.int(n = nrow(clin.perm), replace = FALSE)
                                    clin.perm[, kras.status.field] <- clin.perm[permutation, kras.status.field]
                                    clin.perm$cms_label <- clin.perm$cms_label[permutation]
                                    
                                    kras <- clin.perm[, kras.status.field]
                                    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
                                    cms <- as.character(clin.perm$cms_label)
                                    
                                    ## Run the analysis over the permuted labels
                                    res.i <- ldply(genes, .parallel = FALSE,
                                                   .fun = function(gene) {
                                                       res.gene <- unlist(llply(cms.levels, .parallel = FALSE,
                                                                                .fun = function(prefix) {
                                                                                    llply(kras.states, .parallel = FALSE,
                                                                                          .fun = function(genotype) {
                                                                                              if( ( prefix == "CMS2" ) && ( genotype == kras.states[1] ) ) { return(NULL) }
                                                                                              flag1 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == "CMS2" ) & ( kras.label == kras.states[1] ) ) )
                                                                                              flag2 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == prefix ) & ( kras.label == genotype ) ) )
                                                                                              gene.expr <- es[gene,,drop=T]
                                                                                              res <- wilcox.test(x = gene.expr[flag1], y = gene.expr[flag2])
                                                                                              dir <- sign(mean(gene.expr[flag1], na.rm=TRUE) - mean(gene.expr[flag2], na.rm=TRUE))
                                                                                              p <- res$p.value
                                                                                              w <- res$statistic
                                                                                              ret <- c(p, w, dir)
                                                                                              names(ret) <- c("pval", "w", "dir")
                                                                                              names(ret) <- paste0(prefix, genotype, ".vs.CMS2MT.", names(ret))
                                                                                              ret
                                                                                          })
                                                                                }))
                                                       ret.gene <- c(gene, res.gene)
                                                       names(ret.gene) <- c("gene.set", names(res.gene))
                                                       ret.gene
                                                   })
                                    ret.i <- cbind(rep(i, nrow(res.i)), res.i)
                                    names(ret.i) <- c("perm.i", names(res.i))
                                    ret.i
                                })
    }

    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    
    tests <- ldply(genes, .parallel = FALSE,
                   .fun = function(gene) {
                       res.gene <- unlist(llply(cms.levels, .parallel = FALSE,
                                                .fun = function(prefix) {
                                                    llply(kras.states, .parallel = FALSE,
                                                          .fun = function(genotype) {
                                                              if( ( prefix == "CMS2" ) && ( genotype == kras.states[1] ) ) { return(NULL) }
                                                              flag1 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == "CMS2" ) & ( kras.label == kras.states[1] ) ) )
                                                              flag2 <- !is.na(cms) & !is.na(kras) & ( ( ( cms == prefix ) & ( kras.label == genotype ) ) )
                                                              gene.expr <- es[gene,,drop=T]
                                                              res <- wilcox.test(x = gene.expr[flag1], y = gene.expr[flag2])
                                                              dir <- sign(mean(gene.expr[flag1], na.rm=TRUE) - mean(gene.expr[flag2], na.rm=TRUE))
                                                              p <- res$p.value
                                                              w <- res$statistic
                                                              ret <- c(p, w, dir)
                                                              names(ret) <- c("pval", "w", "dir")
                                                              names(ret) <- paste0(prefix, genotype, ".vs.CMS2MT.", names(ret))
                                                              ret
                                                          })
                                                }))
                       ret.gene <- c(gene, res.gene)
                       names(ret.gene) <- c("gene.set", names(res.gene))
                       ret.gene
                   })

    rownames(tests) <- tests$gene.set
    if(!calc.fdr) { return(list(permuted.tests = NULL, tests = tests)) }
    
    ## Calculate the FDR based on the permuted results
    for(pval.col in colnames(tests)[grepl(pattern="pval", colnames(tests))]) {
        padj.col <- gsub(x=pval.col, pattern="pval", replacement="apval")
        tests[, padj.col] <- NA
    }

    if( (num.permutations > 0) && (calc.fdr.with.bh == FALSE) ) {
        for(pval.col in colnames(tests)[grepl(pattern="pval", colnames(tests)) & !grepl(pattern="apval", colnames(tests))]) {
            padj.col <- gsub(x=pval.col, pattern="pval", replacement="apval")
            dir.col <- gsub(pval.col, pattern="pval", replacement="dir")
            
            tests[, pval.col] <- as.numeric(tests[, pval.col])
            permuted.tests[, pval.col] <- as.numeric(permuted.tests[, pval.col])
            tests[, dir.col] <- as.numeric(tests[, dir.col])
            permuted.tests[, dir.col] <- as.numeric(permuted.tests[, dir.col])
            
            ## Calculate the number of genes that had a positive (negative) direction of change
            pos.flag <- tests[, dir.col] > 0
            neg.flag <- tests[, dir.col] < 0        
            n.pos.genes <- length(which(pos.flag)) + 1
            n.neg.genes <- length(which(neg.flag)) + 1
            
            ## Calculate the fraction of positive (negative) direction of _permuted_ tests with pvalue less than the observed pvalue
            pos.permuted.flag <- permuted.tests[, dir.col] > 0
            neg.permuted.flag <- permuted.tests[, dir.col] < 0        
            n.pos.permuted.genes <- length(which(pos.permuted.flag)) + 1
            n.neg.permuted.genes <- length(which(neg.permuted.flag)) + 1
            
            frac.permuted.more.extreme <- rep(0, nrow(tests))
            frac.permuted.more.extreme[pos.flag] <- unlist(lapply(tests[pos.flag, pval.col], function(p) { ( length(which(permuted.tests[pos.permuted.flag, pval.col] <= p)) + 1 ) / n.pos.permuted.genes } ))
            frac.permuted.more.extreme[neg.flag] <- unlist(lapply(tests[neg.flag, pval.col], function(p) { ( length(which(permuted.tests[neg.permuted.flag, pval.col] <= p)) + 1 ) / n.neg.permuted.genes } ))        
            
            frac.more.extreme <- rep(0, nrow(tests))
            frac.more.extreme[pos.flag] <- unlist(lapply(tests[pos.flag, pval.col], function(p) { ( length(which(tests[pos.flag, pval.col] <= p)) + 1 ) / n.pos.genes } ))
            frac.more.extreme[neg.flag] <- unlist(lapply(tests[neg.flag, pval.col], function(p) { ( length(which(tests[neg.flag, pval.col] <= p)) + 1 ) / n.neg.genes } ))        
            
            ## Clamp max FDR to 1
            tests[, padj.col] <- unlist(lapply(frac.permuted.more.extreme / frac.more.extreme, function(x) ifelse(x < 1, x, 1)))
            
        }
    } else if(calc.fdr.with.bh == TRUE) {
        for(pval.col in colnames(tests)[grepl(pattern="pval", colnames(tests)) & !grepl(pattern="apval", colnames(tests))]) {
            padj.col <- gsub(x=pval.col, pattern="pval", replacement="apval")
            dir.col <- gsub(pval.col, pattern="pval", replacement="dir")

            tests[, pval.col] <- as.numeric(tests[, pval.col])
            tests[, dir.col] <- as.numeric(tests[, dir.col])

            tests[, padj.col] <- p.adjust(tests[, pval.col], method = "BH")
        }
    }
    list(permuted.tests = permuted.tests, tests = tests)
}


doKrasCMSInteractionAnalysis <- function(es.arg, clin.arg, to.plot=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 1000, calc.fdr.with.bh = FALSE, calc.fdr = TRUE, seed = 1234) {

    set.seed(seed)

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    permuted.tests <- NULL

    cms <- as.character(clin$cms_label)    
    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }
    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    
    tests <- ldply(genes, .parallel = FALSE,
                   .fun = function(gene) {

                       ## Subset the data to the gene/gene set of interest
                       esn.st <- as.vector(es[gene,])

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

                       sum.file <- paste("output/", gene, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
                       lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
                       lm.obj.no.interaction <- lm(as.formula(paste("expr", formula.no.interaction, sep=" ~ ")), data=df.lm)        
                       ## library(lmPerm)
                       ## lm.obj <- lmp(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
                       ## lm.obj <- glm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm, family=gaussian(link="log"))
                       lm.sum <- summary(lm.obj)
                       ##        print(lm.sum)
                       capture.output(lm.sum, file = sum.file)

                       ## Also output the overall pvalue
                       f <- lm.sum$fstatistic
                       f.pval <- pf(f[1],f[2],f[3],lower.tail=F)
                       
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
        
                       diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diagnostic.pdf")
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

                       pval.lst <- list(f.pval = as.numeric(unname(f.pval)))
                       for(row.i in 1:nrow(coeffs.2)) {
                           coef <- as.character(coeffs.2$coefficient[row.i])
                           pval <- coeffs.2[row.i,5]
                           pval.lst[[coef]] <- pval
                       }
                       ret.gene <- c(gene, unlist(pval.lst))
                       names(ret.gene) <- c("gene.set", names(pval.lst))
                       ret.gene
                       
                   })

    rownames(tests) <- tests$gene.set
    if(!calc.fdr) { return(list(permuted.tests = NULL, tests = tests)) }
}

plotKrasCMSAnalysis <- function(es.arg, clin.arg, kras.cms.tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL) {

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

    cms <- as.character(clin$cms_label)    
    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }
    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    msi <- clin$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    site <- clin$site
    site.factor <- factor(site)
    neoantigens <- clin$neoantigens
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    exclude.no.lbl <- TRUE
    
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]
        
        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(es[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$Status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$Site <- site.factor
        }
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
        p <- p + ylab(ylab)

        tx.size <- text.size
        if(nchar(ylab) > 20) { tx.size <- 20 }
        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = tx.size))        
##        p <- p + geom_jitter()
        if(("Status" %in% colnames(df)) && ("Site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = Status, shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("Status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = Status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("Site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("Status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("Site" %in% colnames(df)) {        
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
        for(prefix in cms.levels) {
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
                    tmp <- draw.err.bar(st, kras.cms.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix, yoffset, plot.pvals.as.stars = plot.pvals.as.stars, stat.name = ifelse(plot.adjusted.pvalues, "q", "p"))
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

}

doKrasCMSWithinCMSAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, seed = 1234) {

    if( (calc.fdr.with.bh == TRUE) || (calc.fdr == TRUE) || (num.permutations > 0) ) {
        cat("Multiple-testing correction not implemented in doForestAnalysis\n")
        q(status=-1)
    }
    
    set.seed(seed)

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]

    cms <- as.character(clin$cms_label)    
    unique.cms <- sort(unique(cms))
    cmses <- c("CMS1", "CMS4", "CMS2", "CMS3", "NOLBL")
    cmses <- cmses[cmses %in% unique.cms]
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }
    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    permuted.tests <- NULL

    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    
    tests <- ldply(genes, .parallel = FALSE,
                   .fun = function(gene) {
                       res.gene <- unlist(llply(cms.levels, .parallel = FALSE,
                                                .fun = function(cmsLbl) {

                                                    mask <- cmsLbl ==  cms
                                                    pval <- 1
                                                    b <- 0
        
                                                    if((length(which(mask)) >= 2) && (length(unique(kras[mask])) == 2)) {
                                                        ## pval <- apply(esn[to.plot,mask,drop=F], 1, function(x) wilcox.test(x ~ kras[mask])$p.value)
                                                        x <- es[gene,mask,drop=T]
                                                        res <- wilcox.test(x ~ kras[mask])
                                                        pval <- res$p.value
                                                        b <- (mean(x[kras[mask] == 1], na.rm=TRUE) - mean(x[kras[mask] == 0], na.rm=TRUE)) > 0
                                                    }
                                                    ret <- c(pval, b)
                                                    names(ret) <- c("pval", "mutUp")
                                                    names(ret) <- paste0(cmsLbl, ".", kras.states[1], ".vs.", kras.states[2], ".", names(ret))                                                    
                                                    ret
                                                }))
                       ret.gene <- c(gene, res.gene)
                       names(ret.gene) <- c("gene.set", names(res.gene))
                       ret.gene
                   })
    rownames(tests) <- tests$gene.set
    list(permuted.tests = permuted.tests, tests = tests)
}

plotKrasCMSWithinCMSAnalysis <- function(es.arg, clin.arg, kras.cms.within.cms.tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL) {

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    unique.cms <- sort(unique(cms))

    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }

    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    msi <- clin$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    site <- clin$site
    site.factor <- factor(site)
    neoantigens <- clin$neoantigens
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    exclude.no.lbl <- TRUE
    
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]
        
        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(es[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$Status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$Site <- site.factor
        }
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms-within-cms", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))        
##        p <- p + geom_jitter()
        if(("Status" %in% colnames(df)) && ("Site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = Status, shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("Status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = Status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("Site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("Status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("Site" %in% colnames(df)) {        
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
        ## And we will only add 1 error bar (per facet).  So adjust the y axis
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
            tmp <- draw.err.bar(st, kras.cms.within.cms.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix, yoffset)
            g <- tmp[["g"]]
            panel.maxs <- tmp[["panel.maxs"]]
        }
        
        ## Turn clipping off to see the line across the panels
        g$layout$clip <- "off"

        grid.draw(g)
        d <- dev.off()
        
    })
}


doCMSAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, seed = 1234) {

    if( (calc.fdr.with.bh == TRUE) || (calc.fdr == TRUE) || (num.permutations > 0) ) {
        cat("Multiple-testing correction not implemented in doForestAnalysis\n")
        q(status=-1)
    }
    
    set.seed(seed)

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]

    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    permuted.tests <- NULL

    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    cms <- as.character(clin$cms_label)
    
    cms.levels <- unique(cms)
    unique.cms <- sort(unique(cms))    
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    
    tests <- ldply(genes, .parallel = FALSE,
                   .fun = function(gene) {
                       res.gene <- unlist(llply(cms2.first.levels[-1], .parallel = FALSE,
                                                .fun = function(prefix) {
                                                    if(!(prefix %in% cms.levels)) { return(NULL) }        
                                                    cms.str <- prefix
                                                    flag <- !is.na(cms) & ( ( cms == "CMS2") | ( cms == cms.str ) )
                                                    cms.flag <- cms[flag]
                                                    x <- es[gene,flag,drop=T]
                                                    res <- wilcox.test(x ~ cms.flag)
                                                    pval <- res$p.value
                                                    b <- (mean(x[cms.flag == "CMS2"]) - mean(x[cms.flag == cms.str])) > 0
                                                    ret <- c(pval, b)
                                                    names(ret) <- c("pval", "mutUp")
                                                    names(ret) <- paste0(cms.str, ".vs.CMS2.", names(ret))
                                                    ret
                                                }))
                       ret.gene <- c(gene, res.gene)
                       names(ret.gene) <- c("gene.set", names(res.gene))
                       ret.gene
                   })
    rownames(tests) <- tests$gene.set
    list(permuted.tests = permuted.tests, tests = tests)
}

plotCMSAnalysis <- function(es.arg, clin.arg, cms.tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL) {

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    unique.cms <- sort(unique(cms))    
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }

    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    msi <- clin$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    site <- clin$site
    site.factor <- factor(site)
    neoantigens <- clin$neoantigens
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    exclude.no.lbl <- TRUE
    
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]
        
        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(es[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=unique.cms))
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$Status <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$Site <- site.factor
        }
        if(exclude.no.lbl) {        
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-cms", ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = text.size))
        ##        p <- p + geom_jitter()

        if(("Status" %in% colnames(df)) && ("Site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = Status, shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("Status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = Status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("Site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else {
            p <- p + geom_beeswarm(size = bee.sz)
        }
        p <- p + guides(fill = guide_legend(order = 1))        
        if("Status" %in% colnames(df)) {            
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("Site" %in% colnames(df)) {        
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
            pval.col <- paste0(cms.str, ".vs.CMS2.", pval.suffix)
            pval <- as.numeric(cms.tbl[st, pval.col])
        
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
}


doForestAnalysis <- function(es.arg, clin.arg, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, seed = 1234, plot.adjusted.pvalues = FALSE, main = NULL) {

    if( (calc.fdr.with.bh == TRUE) || (calc.fdr == TRUE) || (num.permutations > 0) ) {
        cat("Multiple-testing correction not implemented in doForestAnalysis\n")
        q(status=-1)
    }
    
    set.seed(seed)
    
    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    unique.cms <- sort(unique(cms))

    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }

    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    msi <- clin$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    site <- clin$site
    site.factor <- factor(site)
    neoantigens <- clin$neoantigens
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    exclude.no.lbl <- TRUE

    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    
    kras.col <- paste0("kras", kras.states[1])
    
    cols <- c(kras.col, cms2.last.levels.no.lbl[-1], site.levels[-1])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        cols <- c(cols, msi.levels[-1])
    }
    
    tests <- ldply(1:length(genes), .parallel = FALSE,
                   .fun = function(sti) {

                       ## st will be the gene/gene set to plot
                       st <- genes[sti]

                       ## The y axis label for this gene/gene set
                       ylab <- ylabels[sti]
        
                       cat(paste("Creating forest plot for ", st, "\n", sep=""))

                       ## Subset the data to the gene/gene set of interest
                       esn.st <- as.vector(es[st,])

                       ## Get the KRAS mutation status of each sample
                       kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

                       df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)), site=site.factor)
                       formula <- "KRAS + CMS + site"
                       if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
                           formula <- paste0(formula, " + status")
                           df$Status <- msi.factor
                       }
                       if(!all(is.na(neoantigens))) {
                           formula <- paste0(formula, " + neoantigens")
                           df$neoantigens <- neoantigens
                       }
                       if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
                           df$Site <- site.factor
                       }
                       flag <- unlist(apply(df, 1, function(row) any(is.na(row))))
                       df <- df[!flag,]

                       formula.name <- "full-model"
                       ## Now exclude NOLBL
                       df.lm <- df
                       if(exclude.no.lbl) {
                           df.lm <- df[df$CMS != "NOLBL",]
                       }

                       lm.obj <- lm(as.formula(paste("expr", formula, sep=" ~ ")), data=df.lm)
                       lm.sum <- summary(lm.obj)
                       sum.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-sum.tsv", sep="")
                       capture.output(lm.sum, file = sum.file)
                       
                       forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
                       pdf(forest.file, useDingbats = FALSE)
                       p <- suppressWarnings(forest_model(lm.obj))
                       if(!is.null(main)) {
                           p <- p + ggtitle(main)  + theme(plot.title = element_text(hjust = 0.5))
                       }
                       print(p)
                       d <- dev.off()

                       diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diagnostic.pdf")
                       pdf(diag.file, useDingbats = FALSE)
                       plot(lm.obj)
                       d <- dev.off()
        
                       lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
                       coeffs <- as.data.frame(coefficients(lm.sum))
                       coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
                       write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

                       pval.lst <- list()
                       for(prefix in cols) {
                           col <- paste0(prefix, ".pval")
                           flag <- grepl(x=coeffs.2$coefficient, pattern=prefix, ignore.case=TRUE) & !grepl(x=coeffs.2$coefficient, pattern=":")
                           if(length(which(flag)) != 1) {
                               ##                cat(paste0("Attempting to write ", coeffs.2[flag, 5], " to col ", col, " for ", sti, "\n"))
                               print(col)
                               print(coeffs.2)
                               stop("Could not find coefficient\n")
                           }
                           pval.lst[length(pval.lst)+1] <- coeffs.2[flag, 5]
                       }
                       
                       ## Also output the overall pvalue
                       f <- lm.sum$fstatistic
                       p <- pf(f[1],f[2],f[3],lower.tail=F)
                       vec <- c("gene.set", unlist(pval.lst), p)
                       names(vec) <- c(st, cols, "f.pval")
                       vec
                   })
    rownames(tests) <- tests$gene.set
    tests
}

doKrasAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 1000, calc.fdr.with.bh = FALSE, calc.fdr = TRUE, seed = 1234) {

    if( (calc.fdr.with.bh == TRUE) || (calc.fdr == TRUE) || (num.permutations > 0) ) {
        cat("Multiple-testing correction not implemented in doForestAnalysis\n")
        q(status=-1)
    }

    set.seed(seed)

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)

    tests <- ldply(genes, .parallel = FALSE,
                   .fun = function(gene) {
                       x <- es[gene,,drop=T]
                       res <- wilcox.test(x ~ kras)
                       p <- res$p.value
                       b <- (mean(x[kras == 1]) - mean(x[kras == 0])) > 0
                       vec <- c(gene, p, b)
                       names(vec) <- c("gene.set", "pval", "mutUp")
                       vec
                   })
    rownames(tests) <- tests$gene.set
    list(tests = tests, permuted.tests = NULL)
}


plotKrasAnalysis <- function(es.arg, clin.arg, kras.tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL) {

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    unique.cms <- sort(unique(cms))
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }

    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    msi <- clin$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    site <- clin$site
    site.factor <- factor(site)
    neoantigens <- clin$neoantigens
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    exclude.no.lbl <- TRUE
    
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", st, " ~ ", "KRAS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(es[st,])

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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
        p <- p + ylab(ylab)

        tx.size <- text.size
        if(nchar(ylab) > 20) { tx.size <- 20 }

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS), outlier.shape = NA)
        p <- p + theme(text = element_text(size = tx.size))                
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status), dodge.width = bee.dodge.width, size = bee.sz)
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
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

        ## Get the length of the y axis
        ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
        ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range

        ## The rest of the ugliness below corresponds to the error bars.
        ## Use raw pvalue for univariate analysis
        pval <- as.numeric(kras.tbl[st, pval.suffix])
        ## When we add error bars, the spacing between them will be in units of yoffset
        yoffset <- 0.01 * ( max(df$expr, na.rm=TRUE) - min(df$expr, na.rm=TRUE) )
        ## We only add one error bar.  Shift the ylimit to accomodate this shift.
        mask <- !is.na(df$KRAS) & (df$KRAS == kras.states[1])        
        mt.max <- max(df$expr[mask], na.rm=TRUE)

        ## We will add 1 error bar
        num.error.bars <- 1
        ylimits[2] <- ylimits[2] + 6 * num.error.bars * yoffset
        p <- p + ylim(ylimits)
        
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
        ##        p <- p + annotate("text", x = 1.5, y = mt.wt.max + num.dy * yoffset, label = text, size = sz, vjust = 0)
        p <- p + geom_text( x = 1.5, y = mt.wt.max + num.dy * yoffset, label = text, size = sz, vjust = 0)
        print(p)
        d <- dev.off()
    })

}

doUnivariateAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, seed = 1234, var.name = NULL) {

    if( (calc.fdr.with.bh == TRUE) || (calc.fdr == TRUE) || (num.permutations > 0) ) {
        cat("Multiple-testing correction not implemented in doForestAnalysis\n")
        q(status=-1)
    }

    set.seed(seed)

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    genes <- to.plot
    if(length(genes) == 0) { genes <- rownames(es) }
    genes <- intersect(genes, rownames(es))

    ## Perform the test on the unpermuted data
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    unique.cms <- unique(cms)
    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")    
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    
    msi <- clin$msi
    clin$msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    clin$msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    msi <- clin$msi
    
    site <- clin$site
    site.factor <- factor(site)
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    var <- clin[,var.name]
    var.levels <- levels(factor(var))
    cols <- var.levels[-1]
    base <- var.levels[1]

    df <- data.frame(CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)), site=site.factor)
    
    tests <- ldply(genes, .parallel = FALSE,
                   .fun = function(gene) {
                       vec <- c(unlist(llply(cols, .parallel = FALSE,
                                             .fun = function(lvl) {
                                                 flag1 <- !is.na(var) & (var == base)
                                                 flag2 <- !is.na(var) & (var == lvl)
                                                 x <- es[gene,,drop=F]
                                                 res <- wilcox.test(x = x[flag1], y = x[flag2])
                                                 p <- res$p.value
                                                 names(p) <- paste0(lvl, ".vs.", base, ".pval")
                                                 p
                                             })))
                       ret <- c(gene, vec)
                       names(ret) <- c("gene.set", names(vec))
                       ret
                   })
    rownames(tests) <- tests$gene.set
    list(tests = tests, permuted.tests = NULL)
}

plotUnivariateAnalysis <- function(es.arg, clin.arg, tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL, var.name = NULL) {

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    kras <- clin[, kras.status.field]
    kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
    cms <- as.character(clin$cms_label)
    unique.cms <- sort(unique(cms))
    cms2.first.levels <- c("CMS2", unique.cms[unique.cms != "CMS2"])
    
    cms.levels <- unique(cms)
    exclude.no.lbl <- TRUE
    if(exclude.no.lbl) {
        cms.levels <- unique(cms)[unique(cms) != "NOLBL"]
    }

    ordered.cms <- c("CMS4", "CMS3", "CMS2", "CMS1", "NOLBL")
    cms2.last.levels <- ordered.cms[ordered.cms %in% unique.cms]
    cms2.last.levels.no.lbl <- cms2.last.levels[cms2.last.levels != "NOLBL"]
    
    msi <- clin$msi
    clin$msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    clin$msi[!is.na(msi) & (msi == "mss")] <- "MSS"
    msi <- clin$msi
    site <- clin$site
    site.factor <- factor(site)
    neoantigens <- clin$neoantigens
    msi.factor <- factor(msi)
    msi.levels <- levels(msi.factor)
    site.levels <- levels(site.factor)

    var <- clin[,var.name]
    var.levels <- levels(factor(na.omit(var)))
    cols <- var.levels[-1]
    base <- var.levels[1]
    
    pval.suffix <- "pval"
    if(plot.adjusted.pvalues) {
        pval.suffix <- "apval"
    }
    
    exclude.no.lbl <- TRUE
    
    lapply(1:length(to.plot), function(sti){

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
        cat(paste("Computing ", st, " ~ ", "KRAS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(es[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))
        
        df <- data.frame(expr=esn.st, CMS=factor(cms, levels=cms2.last.levels), KRAS=factor(kras.label, levels=rev(kras.states)), site=site.factor)
        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {            
            df$msi <- msi.factor
        }
        if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
            df$site <- site.factor
        }
        
        flag <- unlist(apply(df[,var.name,drop=F], 1, function(row) any(is.na(row))))
        df <- df[!flag,]

        ## Create the plot
        expr.col <- "expr"
        p <- ggplot(data=df, aes_string(x = var.name, y = expr.col))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }                    
        p <- p + theme(text = element_text(size = text.size))            
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes_string(fill = var.name), outlier.shape = NA)
        ##            p <- p + geom_jitter()
        ## Do not annotate MSI/MSS when we are comparing expr to MSI status
        if(var.name == "msi") {
            p <- p + xlab("status")
            p <- p + guides(fill = guide_legend(title = "status", order = 1))
        } else {
            p <- p + guides(fill = guide_legend(order = 1))
        }
        nxt.order <- 2
        if(("msi" %in% colnames(df)) && ("site" %in% colnames(df)) && (var.name != "site") && (var.name != "msi")) {
            p <- p + geom_beeswarm(aes(colour = msi, shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(title = "status", order = nxt.order))
            nxt.order <- nxt.order + 1                
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = nxt.order))
            nxt.order <- nxt.order + 1
        } else if(("msi" %in% colnames(df)) && (var.name != "msi")) {
            p <- p + geom_beeswarm(aes(colour = msi), dodge.width = bee.dodge.width, size = bee.sz)
            p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(title = "status", order = nxt.order))
            nxt.order <- nxt.order + 1                                
        } else if(("site" %in% colnames(df)) && (var.name != "site")) {
            p <- p + geom_beeswarm(aes(shape = Site), dodge.width = bee.dodge.width, size = bee.sz)
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
        lvls <- levels(df[,var.name,drop=T])
        panel.maxs <- sapply(lvls, function(lvl) {
            mask <- df[,var.name,drop=T] == lvl
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(panel.maxs) <- lvls

        ## Draw the error bars from the first level to each of the others
        for(prefix in cols) {
            str <- prefix
            ## Use raw pvalue for univariate analysis                        
            ## pval.col <- paste0(str, ".apval")
            pval.col <- paste0(str, ".vs.", base, ".", pval.suffix)
            pval <- as.numeric(tbl[st, pval.col])
            
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
        
        ## Create a pdf for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-", var.name, ".pdf", sep="")
        pdf(out.pdf, useDingbats = FALSE)
        
        print(p)
        d <- dev.off()
        
    })
}

doSiteAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, seed = 1234) {
    doUnivariateAnalysis(es.arg, clin.arg, to.plot=to.plot, kras.status.field=kras.status.field, kras.states=kras.states, num.permutations = num.permutations, calc.fdr.with.bh = calc.fdr.with.bh, calc.fdr = calc.fdr, seed = seed, var.name = "site")
}

plotSiteAnalysis <- function(es.arg, clin.arg, tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL) {
    plotUnivariateAnalysis(es.arg, clin.arg, tbl, to.plot=to.plot, ylabels=ylabels, analysis.name=analysis.name, kras.status.field=kras.status.field, kras.states=kras.states, plot.adjusted.pvalues = plot.adjusted.pvalues, plot.pvals.as.stars = plot.pvals.as.stars, main = main, var.name = "site")
}

doMSIAnalysis <- function(es.arg, clin.arg, to.plot=c(), kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, seed = 1234) {
    doUnivariateAnalysis(es.arg, clin.arg, to.plot=to.plot, kras.status.field=kras.status.field, kras.states=kras.states, num.permutations = num.permutations, calc.fdr.with.bh = calc.fdr.with.bh, calc.fdr = calc.fdr, seed = seed, var.name = "msi")
}

plotMSIAnalysis <- function(es.arg, clin.arg, tbl, to.plot=c(), ylabels=c(), analysis.name, kras.status.field="kras", kras.states=c("MT","WT"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = FALSE, main = NULL) {
    plotUnivariateAnalysis(es.arg, clin.arg, tbl, to.plot=to.plot, ylabels=ylabels, analysis.name=analysis.name, kras.status.field=kras.status.field, kras.states=kras.states, plot.adjusted.pvalues = plot.adjusted.pvalues, plot.pvals.as.stars = plot.pvals.as.stars, main = main, var.name = "msi")
}


plotOverlapOfExpressedGeneSets <- function(expr.m, all.gene.sets, gene.set.names, main = NULL) { 
    gene.sets <- list()
    for(population in gene.set.names) {
        genes <- all.gene.sets[[population]][ all.gene.sets[[population]] %in% rownames(expr.m) ]
        gene.sets[[population]] <- genes
    }
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
    
    g <- ggplot(m, aes(x, y))
    g <- g + theme_bw() + xlab("") + ylab("")
    if(!is.null(main)) { g <- g + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }    
    g <- g + geom_tile(aes(fill = value), color='white')
    nz.data <- m[m$value > 0,]
    g <- g + geom_text(data = nz.data, aes(x = x, y = y, label = value))
    g <- g + scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab', guide = guide_colorbar(title = "# Genes"))
    g <- g + theme(axis.text.x=element_text(angle=90),
                  axis.ticks=element_blank(),
                  axis.line=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_line(color='#eeeeee'))
    g
}

## The following comes from Deducer (https://github.com/ifellows/Deducer/blob/master/Deducer/R/util.R)
## which I could not install.
## Renamed this function from d
ddf<-function(..., row.names = NULL, check.rows = FALSE,
            check.names = FALSE,
            stringsAsFactors = FALSE){
  data.frame(...,row.names=row.names,check.rows=check.rows,check.names=check.names,stringsAsFactors=stringsAsFactors)
}

cor.matrix<-function(variables,with.variables,data=NULL,test=cor.test,...){
  arguments <- as.list(match.call()[-1])
  variables<-eval(substitute(variables),data,parent.frame())
  if(length(dim(variables))<1.5){
    variables<-ddf(variables)
    fn <- arguments$variables
    names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
  }
  if(missing(with.variables))
    with.variables <-variables
  else{
    with.variables<-eval(substitute(with.variables),data,parent.frame())
    if(length(dim(with.variables))<1.5){
      with.variables<-ddf(with.variables)
      fn <- arguments$with.variables
      names(with.variables)<-if(is.call(fn)) format(fn) else as.character(fn)
    }		
  }
  cors<-list()
  for(var1 in colnames(variables)){
    cors[[var1]]<-list()
    for(var2  in colnames(with.variables)){
      tmp<-na.omit(data.frame(as.numeric(variables[[var1]]),as.numeric(with.variables[[var2]])))
      names(tmp)<-c(var1,var2)
      cors[[var1]][[var2]]<-test(tmp[[1]],tmp[[2]],...)
      attr(cors[[var1]][[var2]],"N")<-nrow(tmp)
    }
  }
  class(cors)<-"cor.matrix"
  cors
}

plot.smooth.scatter <- function(data, x.col, y.col, x.label, y.label) {
  g <- ggplot(data, aes_string(x = x.col, y = y.col))
  g <- g + stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200)
  colramp = colorRampPalette(c("white", blues9))
  print(colramp)
  ##  g <- g + scale_fill_continuous(low = "white", high = "dodgerblue4", guide=FALSE)
    g <- g + scale_fill_continuous(low = "white", high = blues9[length(blues9)], guide=FALSE)
##  g <- g + geom_point(alpha = 0.1, shape = 20)
  g <- g + xlab(x.label)
  g <- g + ylab(y.label)
  g
}

my_fn <- function(data, mapping, ...){
##    geom_point() + 
##    scale_fill_continuous(low = "white", high = "dodgerblue4", guide=FALSE) +
p <- ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 128) +
    scale_fill_continuous(low = "white", high = blues9[length(blues9)], guide=FALSE) +    
    geom_smooth(method=lm, fill="blue", color="blue", ...)
p
}

my_fn2 <- function(data, mapping, ...){
##    geom_point() + 
    ##    scale_fill_continuous(low = "white", high = "dodgerblue4", guide=FALSE) +
    orig.data <- data
    tmp <- ggplot2:::ggplot.data.frame(data, mapping)
    x.lab <- tmp$labels$x
    y.lab <- tmp$labels$y
    data <- data[,c(x.lab,y.lab)]
    colnames(data) <- c("x", "y")
    nbin <- 128
    bandwidth <- NULL
    map <- grDevices:::.smoothScatterCalcDensity(data, nbin)
    xm <- map$x1
    ym <- map$x2
    dens <- map$fhat
    transformation = function(x) x^0.25
    transformation = function(x) x^0.5  
    dens[] <- transformation(dens)
    dens.df <- as.data.frame(dens)
    colnames(dens.df) <- ym
    ##    dens.df$x <- xm
    x.lab <- "x"
    y.lab <- "y"
    dens.df[,x.lab] <- xm
    ##    dens <- gather_(dens.df, y.lab, "dens", c(y.lab,"dens"))
    dens <- gather(dens.df, y, dens, -x)
    dens[,x.lab] <- as.numeric(dens[,x.lab])
    dens[,y.lab] <- as.numeric(dens[,y.lab])    
    dens$dens <- as.numeric(dens$dens)        
##    dens <- dens[, !(colnames(dens) == "id")]
    ##    data <- data.frame(x = xm, y = ym, dens = dens)
    colramp <- colorRampPalette(c("white", blues9))
    cr <- colramp(256)
    ## theme(legend.position = "none", panel.grid.major = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank() ,axis.text.x = element_blank(), panel.border = element_rect(linetype = "dashed", colour = "black", fill=NA))

    p <- ggplot(data = dens, aes_string(x = x.lab, y = y.lab)) + geom_raster(aes(fill = dens)) +
        scale_fill_continuous(low = cr[1], high = cr[length(cr)], guide=FALSE) +
        geom_smooth(data = orig.data, mapping = mapping, method=lm, color="black", ...) +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank())

p
}

ggally_mycor <- function(data, mapping, ...){
    p <- ggally_cor(data, mapping, ...)
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank())
    p
}

ggally_mysmooth <- function(data, mapping, ...){
    ##               theme(legend.position = "none", panel.grid.major = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank() ,axis.text.x = element_blank(), panel.border = element_rect(colour = "black", size = 5, fill=NA))
    g <- ggally_myDiagAxis(data, mapping, labelSize = 3, ...)
    g <- g +       theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank())

    return(g)
  ggplot(data = data, mapping=mapping) +
      geom_density() +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank())
}

ggally_myDiagAxis <-
function (data, mapping, label = mapping$x, labelSize = 3, labelXPercent = 0.5, 
    labelYPercent = 0.55, labelHJust = 0.5, labelVJust = 0.5, 
    gridLabelSize = 4, ...) 
{
    if (is.null(mapping$x)) {
        stop("mapping$x is null.  There must be a column value in this location.")
    }
    mapping$y <- NULL
    numer <- !GGally:::is_horizontal(data, mapping, "x")
    if (!is.character(label)) {
        label <- deparse(mapping$x)
    }
    label <- gsub(label, pattern="\\.", replacement="\n")
    xData <- GGally:::eval_data_col(data, mapping$x)
    if (numer) {
        xmin <- min(xData, na.rm = TRUE)
        xmax <- max(xData, na.rm = TRUE)
        xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
            (xmax - xmin))
        p <- ggally_text(label = label, mapping = aes(col = "grey50"), 
            xrange = xrange, yrange = xrange, size = labelSize, 
            xP = labelXPercent, yP = labelYPercent, hjust = labelHJust, 
            vjust = labelVJust)
    }
    else {
        breakLabels <- levels(as.factor(xData))
        numLvls <- length(breakLabels)
        p <- ggally_text(label = label, mapping = aes(col = "grey50"), 
            xrange = c(0, 1), yrange = c(0, 1), size = labelSize, 
            yP = labelYPercent, xP = labelXPercent, hjust = labelHJust, 
            vjust = labelVJust)
    }
    p
}


ggcorplot <- function(cor.mat,data=NULL,lines=TRUE,line.method=c("lm","loess"),type="points",
                      alpha=.25,main="auto",var_text_size=5,
                      cor_text_limits=c(5,25),level=.05){
  x_var <- y_var <- trans <- rsq <- p <- x_label <- NULL
  #define a helper function (borrowed from the "ez" package)
  ezLev<-function(x,new_order){
    for(i in rev(new_order)){
      x<-relevel(x,ref=i)
    }
    return(x)
  }							
  
  if(all(line.method==c("lm","loess")))
    line.method<-"lm"	
  
  nm <- names(cor.mat)
  for(i in 1:length(nm))
    dat <- if(i==1) ddf(eval(parse(text=nm[i]),data,parent.frame())) else ddf(dat, eval(parse(text=nm[i]),data,parent.frame()))
  data <- dat
  names(data) <- nm
  # normalize data
  for(i in 1:length(data)){
    data[,i]<-as.numeric(data[,i])
    data[,i]<-(data[,i]-mean(data[,i],na.rm=TRUE))/sd(data[,i],na.rm=TRUE)
  }
  # obtain new data frame
  z<-data.frame()
  i <- 1
  j <- i
  while(i<=length(data)){
    if(j>length(data)){
      i<-i+1
      j<-i
    }else{
      x <- data[,i]
      y <- data[,j]
      temp<-as.data.frame((cbind(x,y)))
      temp<-cbind(temp,names(data)[i],names(data)[j])
      z<-rbind(z,temp)
      
      j<-j+1
    }
  }
  z<-cbind(z,alpha)
  names(z)=c('x_var','y_var','x_label','y_label','trans')
  z$x_label <- ezLev(factor(z$x_label),names(data))
  z$y_label <- ezLev(factor(z$y_label),names(data))
  z=z[z$x_label!=z$y_label,]
  #obtain correlation values
  z_cor <- data.frame()
  i <- 1
  j <- i
  while(i<=length(data)){
    if(j>length(data)){
      i<-i+1
      j<-i
    }else{
      x <- na.omit(data[,i])
      y <- na.omit(data[,j])
      x_mid <- min(x)+diff(range(x))/2
      y_mid <- min(y)+diff(range(y))/2
      this_cor <- cor.mat[[i]][[j]]$estimate
      this_cor.test <- cor.mat[[i]][[j]]
      this_col <- ifelse(this_cor.test$p.value<level,"red"
                         ,"blue")
      this_size <- (this_cor)^2
      cor_text <- ifelse(
        this_cor>0
        ,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
        ,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
      )
      b<-as.data.frame(cor_text)
      b<-cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor<-rbind(z_cor,b)
      j<-j+1
    }
  }
  names(z_cor)<-c('cor','x_mid','y_mid','p','rsq','x_label','y_label')
  z_cor$x_label <- ezLev(factor(z_cor$x_label),names(data))
  z_cor$y_label <- ezLev(factor(z_cor$y_label),names(data))
  diag <- z_cor[z_cor$x_label==z_cor$y_label,]
  z_cor<-z_cor[z_cor$x_label!=z_cor$y_label,]
  
  #start creating layers	
  points_layer <- geom_point(aes(x = x_var, y = y_var, alpha=trans),data = z)
  ## points_layer <- stat_density2d(data = z, aes(x = x_var, y = y_var, fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200)
  points_layer <- stat_density2d(data = z, aes(x = x_var, y = y_var, fill = ..density..^0.25), geom = "tile", contour = FALSE)
  points_layer2 <- scale_fill_continuous(low = "white", high = "dodgerblue4", guide=FALSE)
  print(head(z))
  bin_layer<-geom_hex(data = z, mapping = aes(x = x_var, y = y_var,alpha=trans),bins=10)
  lm_line_layer <- stat_smooth(aes(x=x_var, y=y_var), method=line.method)
  cor_text <- geom_text(aes(x=y_mid, y=x_mid, label=cor, size = rsq, colour = p), data=z_cor)
  var_text <- geom_text(aes(x=y_mid, y=x_mid, label=x_label), data=diag, size=var_text_size)
  
  f <- facet_grid(y_label~x_label,scales='free')
  o <- theme(
    panel.grid.minor = element_blank()
    ,panel.grid.major = element_blank()
    ,axis.ticks = element_blank()
    ,axis.text.y = element_blank()
    ,axis.text.x = element_blank()
    ,axis.title.y = element_blank()
    ,axis.title.x = element_blank()
    ,legend.position='none'
  )
  
  size_scale <- scale_size(limits = c(0,1),range=cor_text_limits)
  the.plot<-ggplot(data=z)
  if(type=="bins")
    the.plot<-the.plot+bin_layer
  else if(type=="points") {
      ##      the.plot<-the.plot+points_layer + points_layer2 + scale_alpha_identity()
      the.plot<-the.plot+points_layer + points_layer2
  }
  the.plot<-the.plot+var_text+
    cor_text+
    f+
    o+
    size_scale
  if(type=="bins")
    the.plot<-the.plot+scale_fill_gradient(low="grey", high="black")
  if(lines)
    the.plot<-the.plot+lm_line_layer
  if(main=="auto")
    main<-cor.mat[[1]][[1]]$method
  the.plot<-the.plot+ggtitle(main)
  return(the.plot)
}

plot.pairwise.correlation <- function(data, main = NULL) {

    df <- data
    colnames(df) <- make.names(colnames(data))
    corr.mat1 <- tryCatch({cor.matrix(variables=do.call("ddf", as.list(df)),,
                                      dat=df,
                                      test=cor.test,
                                      method='spearman',
                                      alternative="two.sided",exact=FALSE)},
                          error = function(e) { 
                              cat("Caught error in cor.matrix")
                              return(NULL) 
                          })
    if(is.null(corr.mat1)) { return(NULL) }
    
    p <- ggcorplot(corr.mat1, data = df, var_text_size = 3)
    p
}

plot.fimm.panel.of.pairs <- function(data, main = NULL) {
  gs <- list()
  for(metric in c("ic50", "auc", "dss1", "dss2")) {
    tmp <- data[,c("SCREEN_ID", "DRUG_ID", metric)]
    flag <- is.finite(tmp[,metric])
    tmp <- tmp[flag,]
    dat <- spread_(data=tmp, key="DRUG_ID", value=metric)
    cols <- unique(tmp$DRUG_ID)
    df <- dat[,cols]
    for(col in 1:ncol(df)) {
      df[,col] <- remove.outliers(df[,col])
    }
##    print(head(df))
##    print(length(which(complete.cases(df))))
##    if(length(which(complete.cases(df))) < 5) { return(NULL) }
    colnames(df) <- make.names(colnames(df))
    corr.mat1<-tryCatch({cor.matrix(variables=do.call("ddf", as.list(df)),,
                          dat=df,
                          test=cor.test,
                          method='spearman',
                          alternative="two.sided",exact=FALSE)},
                        error = function(e) {
                          cat("Caught error in cor.matrix\n")
                          return(NULL)
                        })
    if(is.null(corr.mat1)) { next }
    p <- ggcorplot(corr.mat1,data = df, var_text_size = 3)
    if(is.null(main)) {
      p <- p + ggtitle(toupper(metric))
    } else {
      p <- p + ggtitle(paste0(main, ":\n", toupper(metric)))
    }
    gs[[metric]] <- p
  }
  ## do.call("grid.arrange", gs)
  return(do.call("arrangeGrob", list(grobs = gs, ncol = 2)))
}

my_ggdendro <- function(mat, method = "average", main = NULL) {

    dd <- as.dist((1 - cor(mat))/2)
    hc <- hclust(dd, method=method)
    dd <- as.dendrogram(hc)
    
    suppressPackageStartupMessages(library("dendextend"))
    dend <- hang.dendrogram(dd)
    y.ticks <- seq(from = 0, to = max(get_branches_heights(dend)), by = 0.1)
    g <- ggplot(as.ggdend(dend), offset_labels = - 0.01, theme = theme_minimal())
    g <- g + ylim(-.15, max(get_branches_heights(dend)) )
    g <- g + scale_y_continuous(breaks=y.ticks, limits=c(-0.15, max(get_branches_heights(dend))))
    g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("")
    g <- g + ylab("(1 - Correlation) / 2")
    if(!is.null(main)) {
        g <- g + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5))
    }
    return(g)
    
    ddata <- dendro_data(dd)
    suppressPackageStartupMessages(library("ggdendro"))
    
    g <- ggdendrogram(ddata)
    return(g)
    
##    g <- ggplot(segment(ddata)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
##    g <- g + geom_text(aes(x = x, y = y, label = label, angle = -90, hjust = 0), data = label(ddata)) + scale_y_continuous(expand = c(0.3, 0))
##    g
}

plot.pval.heatmap <- function(tbl, x, y, value, xlab, ylab, use.italics = TRUE, show.values = TRUE, log.transform = TRUE) {
    g <- ggplot(data = tbl, aes_string(x = x, y = y, fill = value))
    g <- g + geom_tile(color="white")
    
    ##    g <- g + scale_fill_gradient(low = "blue", high = "white", labels = c(0, 0.05), space = "Lab", name = "p-value", breaks = c(0, 0.05))

    g <- g + theme_minimal()
    g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), legend.title.align=0.5) 
    if(use.italics) {
        g <- g + theme(axis.text.y = element_text(face = 'italic', angle = 0))
    }
    g <- g + xlab(xlab)
    g <- g + ylab(ylab)
    g <- g + coord_fixed()
    max.value <- ceiling(max(tbl[,value], na.rm=TRUE))
    min.value <- floor(min(tbl[,value], na.rm=TRUE))    
    ##    s <- seq(from = -log10(0.1), to=max.value, by=1)
    s <- seq(from = min.value, to=max.value, by=1)
    labels <- s
##    labels[1] <- "< 1"
    s <- seq(from = -max(abs(min.value), abs(max.value)), to=max(abs(min.value), abs(max.value)), by=1)
    if(log.transform == FALSE) {
        s <- seq(from = -1, to = 1, by = 1)
    }
    print(s)
    labels <- abs(s)
    ct <- 0.2

    darkneg <- "darkblue"
    lightneg <- "lightblue"
    darkpos <- "darkred"
    lightpos <- "pink"
    na.value <- "black"

    if(FALSE) {
        if(log.transform) {
            g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,"white", darkpos), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
        } else {
            g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg,"white", lightpos, darkpos), values = scales::rescale(c(min(s), min(s) + ct, 0, max(s) - ct, max(s)), from = c(min(s), max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
        }
    }
    title1 <- NULL
    if(log.transform) {
        ##        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,"white"), values = scales::rescale(c(0, min(s))), limit = c(min(s),0), space="Lab", name="p-value", labels = rev(labels[s <= 0]), breaks = s[s <= 0], guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Upregulated, -Log[10]~italic('p')))))
        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,"white"), values = scales::rescale(c(0, min(s))), limit = c(min(s),0), space="Lab", name="p-value", labels = rev(labels[s <= 0]), breaks = s[s <= 0], guide=guide_colorbar(title=""))
    title1 <- textGrob(bquote(atop(italic('KRAS')~Upregulated, -Log[10]~italic('p'))))
    } else {
        ##        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg,"white"), values = scales::rescale(c(min(s), min(s) + ct, 0), from = c(min(s), 0)), limit = c(min(s),0), space="Lab", labels = c(0, ct, 1), breaks = c(min(s), min(s) + ct, 0), name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Upregulated, italic('q')-value))))
                g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg,"white"), values = scales::rescale(c(min(s), min(s) + ct, 0), from = c(min(s), 0)), limit = c(min(s),0), space="Lab", labels = c(0, ct, 1), breaks = c(min(s), min(s) + ct, 0), name="upq-value", guide=guide_colorbar(title=""))
    title1 <- textGrob(bquote(atop(italic('KRAS')~Upregulated, italic('q')-value)))
    }
    non.na.tbl <- na.omit(tbl)
    if(log.transform) {
        non.na.tbl$value <- round(10^-abs(non.na.tbl$value), digits=2)
    } else {
        non.na.tbl$value <- round(-(abs(non.na.tbl$value)-1), digits=2)
    }
    if(show.values) {
        g <- g + geom_text(data = non.na.tbl, aes_string(x = x, y = y, label = value))
    }

    gt <- ggplotGrob(g)
    g.legend1 <- gtable::gtable_filter(gt, "guide-box")

    if(log.transform) {
        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,"white", darkpos), values = scales::rescale(c(min(s), 0, max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
    } else {
        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg, "white", lightpos, darkpos), values = scales::rescale(c(min(s), min(s) + ct, 0, max(s) - ct , max(s)), from = c(min(s), max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
    }
    gt <- ggplotGrob(g)
    g2 <- ggplot(data = tbl, aes_string(x = x, y = y, fill = value))
    g2 <- g2 + geom_tile(color="white")
    g2 <- g2 + theme_minimal()
    g2 <- g2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), legend.title.align=0.5)
  
    
    title2 <- NULL
    if(log.transform) {
        ##        g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, "white"), values = scales::rescale(c(max(s), 0)), limit = c(0,max(s)), space="Lab", name="p-value2", labels = (labels[s >= 0]), breaks = (s[s >= 0]), guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Downregulated, -Log[10]~italic('p')))))
                g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, "white"), values = scales::rescale(c(max(s), 0)), limit = c(0,max(s)), space="Lab", name="p-value2", labels = (labels[s >= 0]), breaks = (s[s >= 0]), guide=guide_colorbar(title=""))
    title2 <- textGrob(bquote(atop(italic('KRAS')~Downregulated, -Log[10]~italic('p'))))    
## + guides(colour = guide_legend(title.hjust=-0.2)) 
    } else {
        ##        g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, lightpos,"white"), values = scales::rescale(c(0, ct, 1), from = c(0, 1)), labels = c(0, ct, 1), breaks = c(0, ct, 1), limit = c(0,1), space="Lab", name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Downregulated, italic('q')-value))))
                g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, lightpos,"white"), values = scales::rescale(c(0, ct, 1), from = c(0, 1)), labels = c(0, ct, 1), breaks = c(0, ct, 1), limit = c(0,1), space="Lab", name="upq-value", guide=guide_colorbar(title=""))
    title2 <- textGrob(bquote(atop(italic('KRAS')~Downregulated, italic('q')-value)))
    }

    gt2 <- ggplotGrob(g2)
    
    g.legend2 <- gtable::gtable_filter(gt2, "guide-box")    
    
##    g <- g + scale_colour_gradient(limits=c(0,5))
    ##    g
    ##    grid.arrange(arrangeGrob(gt, g.legend1, g.legend1, nrow=1))
    g.legend.title1 <- textGrob("hi")
    h <- grobHeight(g.legend1)
    w <- grobWidth(g.legend2)
    title <- textGrob("Title", y=unit(0.5,"npc") + 0.5*h, 
                      vjust=0, gp=gpar(fontsize=20))
    ##    g.legend1 <- gTree(children=gList(title, g.legend1))
    ##    g.legend1 <- gtable_add_rows(g.legend1, heights = grobHeight(title), pos = 0)
##    g.legend1 <- gtable_add_rows(g.legend1, heights = grobHeight(title), pos = 0)    
##    g.legend1 <- gtable_add_grob(g.legend1, list(title), t=c(1), l=c(1), r=ncol(g.legend1))
    grobs <- list(gt, title1, g.legend1, title2, g.legend2)
    matrix <- rbind(c(1,2),c(1,3),c(1,NA),c(1,4),c(1,5),c(1,NA))
    ## grid.arrange(grobs = grobs, layout_matrix = matrix)
    arrangeGrob(grobs = grobs, layout_matrix = matrix)
}

plot.circ.vs.neoantigens <- function(es.arg, clin.arg, main = NULL) {
    sig1 <- "CIRC"
    sig2 <- "neoantigens"

    sig1.name <- "CIRC Enrichment Score"
    sig2.name <- "Log10 ( # Neoantigens + 1 ) "

    clin <- clin.arg
    es <- es.arg
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(clin$sample, colnames(es))
    clin <- clin[!is.na(idxs),]
    es <- es[, na.omit(idxs)]
    
    msi <- clin$msi
    msi[!is.na(msi) & (msi == "msi")] <- "MSI"
    msi[!is.na(msi) & (msi == "mss")] <- "MSS"    
    msi.factor <- factor(msi)
    
    df <- data.frame(x = es[sig2,], y = es[sig1,])

    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        df$Status <- msi.factor
    }
    
    g <- ggplot(data = df, aes(x = x, y = y, label = y), size = bee.sz)
    if(!is.null(main)) { g <- g + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
    df <- df[!is.na(df$x),]
    df <- df[!is.na(df$y),]        
    if("status" %in% colnames(df)) {
        g <- g + geom_point(aes(colour = status))
        g <- g + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))            
    } else {
        g <- g + geom_point()
    }
    g <- g + ylab(sig1.name)
    g <- g + xlab(sig2.name)
    g <- g + theme(text = element_text(size = text.size))        

    ## When comparing CIRC to neoantigens, don't plot the linear
    ## fit.  Instead, break into quadrants based on densities.
    ## (But y axis should be separated at CIRC enrichment = 0)
    ## Also, add MSI status.
    ## Calculate fisher's.
    d <- density(na.omit(es[sig2,]))
    x.upper <- 2.75
    x.lower <- 2.0
    flag <- (d$x < x.upper) & (d$x > x.lower)
    x.flag <- d$x[flag]
    y.flag <- d$y[flag]
    sep.line <- x.flag[which(min(y.flag)==y.flag)[1]]
    res <- fisher.test(table(df$x < sep.line, df$y < 0))
    cat(paste0("Analyzing ", sig1, " vs ", sig2, " using ", res$method, ": p = ", res$p.value, " n = ", nrow(df), "\n"))
            
    g <- g + geom_hline(yintercept = 0, linetype = 'dashed')
    g <- g + geom_vline(xintercept = sep.line, linetype = 'dashed')
    g
    
}


## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.

do.msi <- function(expr, clin, file, neoantigen.cnts = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification", num.clusters = 2) {

    ## This is the MSI signature from Tian et al.
    ## There are 64 genes here--though one is "Unknown."
    msi.signature <-  c("ACSL6", "AGR2", "ARID3A", "ASCL2", "ASXL1", "ATP9A", "AXIN", "BC000986", "C10orf47", "C13orf18", "C20orf11", "C20orf43", "CEACAM3", "CEACAM5", "CEP68", "DIDO1", "DUSP18", "DYNLRB1", "EP300", "EPDR1", "FBXO34", "GGA2", "GGT7", "GNG4", "GPR143", "GUCY2C", "HNRNPL", "KCNK5", "KHDRBS3", "KRT23", "LFNG", "LMO4", "LOC157860", "MDM2", "MLH1", "OIT3", "PLAGL2", "PPP1R3D", "PRR15", "QPRT", "RNF43", "ROCK2", "RPL22L1", "SHROOM2", "SHROOM4", "SLC25A22", "SMAD2", "SMCR7L", "SORBS1", "STRN3", "TCF7", "TFCP2L1", "TGFBR2", "TNFSF9", "TNNC2", "TNNT1", "TRIM7", "TSPAN6", "UNKL", "Unknown", "VAV3", "VNN2", "ZFP36L2", "ZSWIM3")
    
    ## Based on:
    ##     http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
    ##     http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically/17371177#17371177

    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("gtable"))
    suppressPackageStartupMessages(library("ggdendro"))    

    clin.orig <- clin
    
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

    cat("MSI nrow = ", nrow(clin.mask), "\n")
    
    has.cms <- any(!is.na(clin.mask$cms_label))    
    has.braf <- any(!is.na(clin.mask$braf))
    has.tyms <- "TYMS" %in% rownames(expr.mask)
    has.neoantigens <- !is.null(neoantigen.cnts)
    num.annotation.legends <- 0
    if(has.msi) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    if(has.braf) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    if(has.cms) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    if(has.tyms) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    num.cols <- num.annotation.legends
    if(has.neoantigens) {
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
    clusters <- cutree(hc.col, k = num.clusters)
    clusters <- clusters[col.ord]
    cluster.freq <- as.data.frame(table(clusters), stringsAsFactors=FALSE)
    print(cluster.freq)

    ## Assume (for now) that cluster with more members is the MSS cluster.
    ## In the future, look at the expression of the genes
    min.freq.cluster <- NULL
    msi.clusters <- c()
    mss.clusters <- c()
    if(num.clusters == 2) {
        min.freq.cluster <- cluster.freq$clusters[which(cluster.freq$Freq == min(cluster.freq$Freq))]
        msi.clusters <- c(as.numeric(min.freq.cluster))
        mss.clusters <- c(2 - msi.clusters[1] + 1)
    } else {
        ## Assign the largest cluster to MSS; the second largest to MSI.
        ## And the remain to which is closest to their mean
        cluster.freq <- cluster.freq[order(cluster.freq$Freq, decreasing=TRUE),]
        mss.clusters <- c(cluster.freq$clusters[1])
        msi.clusters <- c(cluster.freq$clusters[2])        

        mss.samples <- names(clusters)[clusters == mss.clusters[1]]
        msi.samples <- names(clusters)[clusters == msi.clusters[1]]
        mss.mean.expr <- unname(rowMeans(expr[,mss.samples]))
        msi.mean.expr <- unname(rowMeans(expr[,msi.samples]))
        
        for(i in 3:nrow(cluster.freq)) {
            cluster.samples <- names(clusters)[clusters == cluster.freq$clusters[i]]
            mean.expr <- unname(rowMeans(expr[,cluster.samples,drop=F]))
            diff.mss <- mss.mean.expr - mean.expr
            diff.msi <- msi.mean.expr - mean.expr
            print(length(mss.mean.expr))
            if(sum(diff.mss*diff.mss) < sum(diff.msi*diff.msi)) {
                mss.clusters <- c(mss.clusters, cluster.freq$clusters[i])
            } else {
                msi.clusters <- c(msi.clusters, cluster.freq$clusters[i])
            }
        }
    }
    cat(paste0("MSI clusters: ", paste(msi.clusters, collapse=","), "\n"))
    cat(paste0("MSS clusters: ", paste(mss.clusters, collapse=","), "\n"))
    
    msi.assign.tbl <- data.frame(sample=names(clusters), msi.inferred=unname(clusters), stringsAsFactors=FALSE)
    msi.assign.tbl$msi.inferred[msi.assign.tbl$msi.inferred %in% msi.clusters] <- "msi"
    msi.assign.tbl$msi.inferred[msi.assign.tbl$msi.inferred %in% mss.clusters] <- "mss"    

    clin.orig <- merge(clin.orig, msi.assign.tbl, by="sample", all=TRUE)
    
    cross.over.sample <- NULL
    for(i in 2:length(clusters)) {
        if(clusters[i-1] != clusters[i]) {
            cross.over.sample <- names(clusters)[i]
        }
    }

    ## Re-order the expression matrix based on the clustering
    expr.mask <- (expr.mask[row.ord, col.ord])
    expr.mask.with.tyms <- expr.mask.with.tyms[, col.ord]    
    clin.mask <- clin.mask[col.ord, ]
    
    remove.msi.labels <- TRUE
    remove.labels <- TRUE
    remove.braf.labels <- TRUE
    remove.cms.labels <- TRUE
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

    if(has.neoantigens) {
        ## Make the mutation samples factors
        mutation.cnt.by.sample <- neoantigen.cnts[levels(df$sample)]
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
    msi.tbl <- NULL
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

    g.legend.cms <- NULL
    if(has.cms) {
        cat("Has CMS\n")        
        cms.tbl <- data.frame(sample=clin.mask$sample, cms=clin.mask$cms_label)
        cms.tbl$sample <- factor(cms.tbl$sample, levels=levels, ordered=TRUE)
        cms.tbl$cms <- as.character(cms.tbl$cms)
        
        cms.melt <- melt(cms.tbl, id.vars="sample")
        p.cms <- ggplot(data = cms.melt, aes(x = variable, y=sample))
        p.cms <- p.cms + geom_tile(aes(fill = value), show.legend = TRUE)
        p.cms <- p.cms + scale_fill_manual(values = c("black", "gray", "red", "blue", "yellow"), na.value = "white")    
        p.cms <- p.cms + guides(fill = guide_legend(title="CMS", title.hjust = 0.5))
        p.cms <- p.cms + theme(legend.text = element_text(size = fs))        
        p.cms <- p.cms + geom_hline(yintercept = cross.over.sample.indx, linetype = 'dashed')    
        g.cms <- ggplotGrob(p.cms)
        g.legend.cms <- gtable::gtable_filter(g.cms, "guide-box")
        g.cms <- gtable::gtable_filter(g.cms, "panel")

        if(!remove.cms.labels) {
            g.cms <- ggplotGrob(p.cms)
        }
        
        cms.label <- textGrob("CMS", gp=gpar(fontsize=fs), just=just, rot=90)
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- cms.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.cms
        indx <- indx + 1
    }

    g.legend.tyms <- NULL
    tyms.tbl <- NULL
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
        ##        p.tyms <- p.tyms + scale_fill_gradient2(high = "green", mid = "black", low = "red")
        p.tyms <- p.tyms + scale_fill_gradient2(high = "green", mid = "black", low = "blue")                 
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
    
    if(has.neoantigens) {
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
        str <- paste(strwrap("Num Neoantigens", width=12), collapse="\n")
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
    if(has.cms) {
        row <- c(plot.row, indx)
##        row <- spacer.row
##        row[num.cols] <- indx
        layout <- rbind(layout, row)
        grobs[[indx]] <- g.legend.cms
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
    if(has.neoantigens) {
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
    ## pdf(paste0("output/", file, ".pdf"), useDingbats = FALSE)
    ## grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    msi.g <- arrangeGrob(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    ret.lst <- list("msi.g" = msi.g, "tyms.msi.inferred.g" = NULL, "tyms.msi.annotated.g" = NULL)

    ## d <- dev.off()

    if(has.tyms) {
        ## Plot TYMS vs inferred MSI (and annotated MSI, if available)
        if(has.msi) {
            df <- merge(tyms.tbl, clin.orig, by="sample", all=FALSE)
            df.imputed <- data.frame(Status = df$msi.inferred, TYMS = df$tyms, Annotation = rep("Imputed", nrow(df)))
            df.annotated <- data.frame(Status = df$msi, TYMS = df$tyms, Annotation = rep("Annotated", nrow(df)))
            df <- rbind(df.imputed, df.annotated)
            df$Status <- factor(toupper(df$Status))
            df$Annotation <- factor(df$Annotation)
            
            p <- ggplot(data=df, aes(x=Status, y=TYMS))
            p <- p + theme(text = element_text(size = text.size))            
            p <- p + ylab(expression(italic(TYMS)~"Expression"))

            ## Create a box plot where the x axis is status ...
            p <- p + geom_boxplot(aes(fill=Status))
            p <- p + geom_beeswarm()
            p <- p + facet_grid(~ Annotation)

            ## The rest of the ugliness below corresponds to the error bars.
            
            ## Get the length of the y axis
            ##        ylimits <- ggplot_build(p)$panel$ranges[[1]]$y.range
            ylimits <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
            
            ## When we add error bars, the spacing between them will be in units of yoffset
            yoffset <- 0.01 * ( max(ylimits, na.rm=TRUE) - min(ylimits, na.rm=TRUE) )
            ## And we will only add 1 error bar.  So adjust the y axis
            ## to accommodate this shift.
            num.error.bars <- 1
            ylimits[2] <- ylimits[2] + 6 * num.error.bars * yoffset
            p <- p + ylim(ylimits)
            
            gb <- ggplot_build(p)
            g <- ggplot_gtable(gb)

            ranges <- gb$layout$panel_ranges
            
            facets <- levels(df$Annotation)
            statuses <- levels(df$Status)
            
            msi.panel.maxs <- sapply(facets, function(lbl) {
                mask <- df$Annotation == lbl & df$Status == "MSI"
                m <- max(df$TYMS[mask], na.rm=TRUE)
                m
            })

            names(msi.panel.maxs) <- sapply(facets, function(x) paste0(x, "-", "MSI"))

            mss.panel.maxs <- sapply(facets, function(lbl) {
                mask <- df$Annotation == lbl & df$Status == "MSS"
                m <- max(df$TYMS[mask], na.rm=TRUE)
                m
            })

            names(mss.panel.maxs) <- sapply(facets, function(x) paste0(x, "-", "MSS"))
            
            panel.maxs <- c(msi.panel.maxs, mss.panel.maxs)

            ## Add the error bars between MSI and MSS expression values (within a facet)

            pval.tbl <- data.frame(variable=c("tyms"))
            rownames(pval.tbl) <- pval.tbl$variable
            for(facet in facets) {
                df.facet <- df[!is.na(df$Annotation) & (df$Annotation == facet),]
                res <- wilcox.test(as.formula(paste("TYMS", "Status", sep=" ~ ")), data = df.facet)
                pval <- as.numeric(res$p.value)
                cat(paste0(file, ": Analyzing TYMS vs ", facet, " Status using ", res$method, ": W = ", res$statistic, " p = ", pval, " ", paste(unlist(lapply(unique(df.facet[,"Status"]), function(var) { paste0("n_", var, " = ", length(which(df.facet[,"Status"] == var))) })), collapse=", "), "\n"))
                pval.col <- paste0(facet, ".pval")
                pval.tbl["tyms", pval.col] <- pval
            }

            for(indx in 1:length(facets)) {
                lbl1 <- facets[indx]
                lbl2 <- facets[indx]
                cmpLbl <- facets[indx]
##                if(indx == 2) { indx <- 3 }
                indx1 <- indx
                indx2 <- indx
                states <- statuses
                tmp <- draw.err.bar("tyms", pval.tbl, ranges, panel.maxs, g, lbl1, statuses[1], lbl2, statuses[2], indx1, indx2, cmpLbl, states, pval.suffix = "pval", yoffset)
                g <- tmp[["g"]]
                panel.maxs <- tmp[["panel.maxs"]]
            }
            
            ## Turn clipping off to see the line across the panels
            g$layout$clip <- "off"

##            pdf(paste0("output/", file, "-tyms-vs-msi.pdf"), useDingbats = FALSE)
##            grid.draw(g)
##            d <- dev.off()

            ret.lst[["tyms.msi.annotated.g"]] <- g
            
        } else {
            df <- merge(tyms.tbl, clin.orig, by="sample", all=FALSE)
            df <- df[,c("msi.inferred", "tyms")]
            colnames(df) <- c("Inferred.Status", "TYMS")
            df$Inferred.Status <- factor(toupper(df$Inferred.Status))
            ## pdf(paste0("output/", file, "-tyms-vs-annotated-msi.pdf"), useDingbats = FALSE)
            plt <- plot.beeswarm.with.err(df, x.name = "Inferred.Status", y.name = "TYMS", x.lab = "Inferred Status", y.lab = expression(italic(TYMS)~"Expression"))
            ## print(plt)
            ## d <- dev.off()
            ret.lst[["tyms.msi.inferred.g"]] <- plt
            
        }

        ## Plot TYMS vs annotated MSI
    }

    ret.lst[["clin"]] <- clin.orig
    return(ret.lst)
}


do.all.analyses <- function(es, clin.es, dataset, gene.sets.to.analyze, gene.sets.to.plot, gene.sets.to.plot.labels) {

    res <- doKrasCMSInteractionAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, analysis.name = dataset, num.permutations = 0, calc.fdr = FALSE)
    ## NB: <<- modifies a global variable
    kras.cms.interaction.res[[dataset]] <<- res$tests
    kras.cms.interaction.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms-interaction.tsv"), kras.cms.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    res <- doKrasCMSAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.res[[dataset]] <<- res$tests
    kras.cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms.tsv"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doKrasCMSWithinCMSAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.within.cms.res[[dataset]] <<- res$tests
    kras.cms.within.cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms-within-cms.tsv"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    res <- doKrasAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, num.permutations = 0, calc.fdr = FALSE)
    kras.res[[dataset]] <<- res$tests
    kras.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras.tsv"), kras.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doCMSAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, num.permutations = 0, calc.fdr = FALSE)
    cms.res[[dataset]] <<- res$tests
    cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-cms.tsv"), cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doSiteAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, num.permutations = 0, calc.fdr = FALSE)
    site.res[[dataset]] <<- res$tests
    site.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-site.tsv"), site.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


    if( ("msi" %in% colnames(clin.es)) & (!all(is.na(clin.es$msi))) ) {
        res <- doMSIAnalysis(es, clin.es, to.plot=gene.sets.to.analyze, num.permutations = 0, calc.fdr = FALSE)
        msi.res[[dataset]] <<- res$tests
        msi.tbl <- res$tests
        write.table(file = paste0(dataset, "-immune-hallmark-msi.tsv"), msi.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    
    plotKrasCMSAnalysis(es, clin.es, kras.cms.tbl, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

    plotKrasCMSWithinCMSAnalysis(es, clin.es, kras.cms.within.cms.tbl, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))    

    plotKrasAnalysis(es, clin.es, kras.tbl, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

    plotCMSAnalysis(es, clin.es, cms.tbl, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

    plotSiteAnalysis(es, clin.es, site.tbl, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))        

    if( ("msi" %in% colnames(clin.es)) & (!all(is.na(clin.es$msi))) ) {
        plotMSIAnalysis(es, clin.es, msi.tbl, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))    
    }
    
    forest.tbl <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=gene.sets.to.plot, ylabels=gene.sets.to.plot.labels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
    
}

## TODO:

## univariate kras (scatter of fdr-corrected)
##   - bindea
##   - mcpcounter
##   - hallmarks

## how best to show multivariate analysis (fdr correction?)
##   - bindea
##   - mcpcounter
##   - hallmarks
## revise plot to count # and to correct

## copy to new directory
## remove: doAnalyses
## remove any unnecessary code
## find proper way to incorporate mcpcounter

## label figures and save them to appropriate fig/supp file name

## check to see how figures look:
##   adjust height of error bars
##   star size consistency

## check all analysis

## check everything into github
