## devtools::install_github("NikNakk/forestmodel")                              
library("forestmodel")
suppressPackageStartupMessages(library("ggbeeswarm"))

inverse.norm.transform <- function(x) {
##    p <- 2*pnorm(abs(x), lower.tail=FALSE) 
##    x2 <- qnorm(p/2, lower.tail=FALSE)*sign(x)
    ##    x2
    qq <- qqnorm(x, plot.it = FALSE)
    trn <- ( mean(x, na.rm=TRUE) + ( sd(x, na.rm=TRUE) * qq$x ) )
    trn
}

## train.set will be used to optimize alpha, unless use.published.alpha
## is TRUE.
fit.rasness.harmonized <- function(clin, harmonized_expr, train.set, optimize.set, prefix, use.published.alpha = TRUE) {
        
    data.sets <- unique(clin$dataset[!is.na(clin$kras)])
    ## [1] "kfs"    "nki"    "french" "petacc" "tcga"   "amc"  

    if(!(train.set %in% data.sets)) {
        stop(paste0("Training set ", train.set, " not in harmonized data sets: ", paste(data.sets, collapse=","), "\n"))
    }

    if(!(optimize.set %in% data.sets)) {
        stop(paste0("Training set ", train.set, " not in harmonized data sets: ", paste(data.sets, collapse=","), "\n"))
    }

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
    pdf(paste0(prefix, "-aucs.pdf"))
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
    best.alpha <- alphas[which(best.aucs == max(best.aucs))]
    ## Let's just use alpha = 0.1, which is what Justin published.
    if(use.published.alpha) {    
        best.alpha <- 0.1
    }
    best.train.harm.trained.model <- train.harm.models[[which(alphas == best.alpha)]]

    ## Optimize the cut point using the optimize.set
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
    
    pdf(paste0(prefix, "-rasness-vs-mt.pdf"))
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

        df <- data.frame(score=clin$score.rasness, kras.not.factor=clin$kras, KRAS=as.factor(clin$kras))
        df <- df[!is.na(df$score) & !is.na(df$KRAS),]
        wilcox.test(score ~ KRAS, data = df)
        p <- ggplot(data=df, aes(x=KRAS, y=score))
        p <- p + ggtitle(paste0("Harmonized ", test.set, " RASness (trained on ", train.set, "; optimized on ", optimize.set, ") vs KRAS MT"))
        p <- p + ylab("RIS")

        ## Create a box plot where the x axis is KRAS mutation status ...                                                                                                   
        p <- p + geom_boxplot(aes(fill=KRAS))
        ##        p <- p + geom_jitter()
        p <- p + geom_beeswarm()        
        hline.data <- data.frame(y = cut.pt)
        p <- p + geom_hline(aes(yintercept = y), hline.data)

        ## The rest of the ugliness below corresponds to the error bars.
        ## Use raw pvalue for univariate analysis
        ##        pval <- as.numeric(kras.tbl$apval[sti])
        pval <- as.numeric(wilcox.test(score ~ KRAS, data = df)$p.value)

        
        yoffset <- 0.01 * ( max(df$score) - min(df$score) )
        ## We only add one error bar.  Shift the ylimit to accomodate this shift.
        mask <- df$kras.not.factor == 1
        mt.max <- max(df$score[mask])

        mask <- df$kras.not.factor == 0
        wt.max <- max(df$score[mask])

        mt.wt.max <- max(mt.max, wt.max)
        
        ## Add the error bars between KRAS MT and KRAS WT expression values (within a facet)
        path.df <- data.frame(x = c(1, 1, 2, 2), y=c(mt.max + yoffset, mt.wt.max + 2 * yoffset, mt.wt.max + 2 * yoffset, wt.max + yoffset))
        p <- p + geom_path(data=path.df, aes(x, y))
        p <- p + annotate("text", x = 1.5, y = mt.wt.max + 4 * yoffset, label = pval.to.text(pval))
        
        print(p)
    }
    d <- dev.off()

    flag <- unlist(apply(clin[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 1) && (row[2] == 1), FALSE)))
    clin$rasness.and.mut[flag] <- 1
    flag <- unlist(apply(clin[,c("rasness", "kras")], 1, function(row) ifelse(!any(is.na(row)), (row[1] == 0) && (row[2] == 0), FALSE)))
    clin$rasness.and.mut[flag] <- 0
    flag <- !is.na(clin$nras) & (clin$nras == 1)
    clin$rasness.and.mut[flag] <- NA
    clin
}

plot.beeswarm.with.err <- function(df, x.name, y.name, x.lab = NULL, y.lab = NULL) {
    if(is.null(x.lab)) { x.lab <- x.name }
    if(is.null(y.lab)) { y.lab <- y.name }
    df <- df[!is.na(df[,x.name]) & !is.na(df[,y.name]),]
    p <- ggplot(data=df, aes_string(x=x.name, y=y.name))
##    p <- p + ggtitle(paste0(y.name, " vs ", x.name))
    p <- p + ylab(y.lab)
    p <- p + xlab(x.lab)
    
    p <- p + geom_boxplot(aes_string(fill=x.name))
    p <- p + geom_beeswarm()        
    p <- p + guides(fill = guide_legend(title = x.lab))
    
    ## The rest of the ugliness below corresponds to the error bars.
    ## Use raw pvalue for univariate analysis
    ##        pval <- as.numeric(kras.tbl$apval[sti])
    res <- wilcox.test(as.formula(paste(y.name, x.name, sep=" ~ ")), data = df)
    pval <- as.numeric(res$p.value)
    cat(paste0("Beeswarm: ", y.name, " vs ", x.name, " W = ", res$statistic, " p = ", pval, " ", paste(unlist(lapply(unique(df[,x.name]), function(var) { paste0("n_", var, " = ", length(which(df[,x.name] == var))) })), collapse=", "), "\n"))

    yoffset <- 0.01 * ( max(df[,y.name], na.rm=TRUE) - min(df[,y.name], na.rm=TRUE) )
    mask <- !is.na(df[, x.name]) & (df[, x.name] == levels(df[, x.name])[1])
    max1 <- max(df[mask, y.name], na.rm=TRUE)

    mask <- !is.na(df[, x.name]) & (df[, x.name] == levels(df[, x.name])[2])
    max2 <- max(df[mask, y.name], na.rm=TRUE)

    max12 <- max(max1, max2)
        
    ## Add the error bars
    path.df <- data.frame(x = c(1, 1, 2, 2), y=c(max1 + yoffset, max12 + 2 * yoffset, max12 + 2 * yoffset, max2 + yoffset))
    p <- p + geom_path(data=path.df, aes(x, y))
    p <- p + annotate("text", x = 1.5, y = max12 + 4 * yoffset, label = pval.to.text(pval))
    p
}

plot.cibersort.beeswarm <- function(cibersort.mat, clin, file) {

    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("ggdendro"))    

    clin.orig <- clin
    tmp <- clin
    tmp <- tmp[!is.na(tmp$cms_label) & (tmp$cms_label %in% c("CMS1", "CMS2", "CMS3", "CMS4")),]
    tmp$sample <- gsub(tmp$sample, pattern="-", replacement=".")

    cibersort.tmp <- cibersort.mat
    cibersort.tmp <- cibersort.tmp[cibersort.tmp$P.value < 0.01,]
    
    ## Match up the clinical annotations and the expression data.
    idxs <- match(tmp$sample, cibersort.tmp$Input.Sample)
    clin.mask <- tmp[!is.na(idxs),]
    cibersort.mask <- cibersort.tmp[na.omit(idxs), ]
    rownames(cibersort.mask) <- cibersort.mask$Input.Sample
    cibersort.data <- cibersort.mask[,!(colnames(cibersort.mask) %in% c("Input.Sample", "P.value", "Pearson.Correlation", "RMSE"))]

    df <- clin.mask[, c("kras"), drop=F]
    rownames(df) <- clin.mask$sample
    df <- df[!is.na(df$kras),,drop=F]
    flag <- df$kras == 0
    df$kras[flag] <- "WT"
    df$kras[!flag] <- "MT"
    df$kras <- factor(df$kras)
    
    df <- merge(cibersort.data, df, by.x = "row.names", by.y = "row.names")
    pdf(paste0(file, ".pdf"))
    for(col in colnames(cibersort.data)) {
        cat(paste0("Making scatter plot for cibersort data: ", col, "\n"))
        print(unique(df$kras))
        p <- plot.beeswarm.with.err(df, "kras", col)
        print(p)
    }
    d <- dev.off()
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

    out.pdf <- paste("output/", analysis.name, "-kras-cms-meta-gene", ".pdf", sep="")
    pdf(out.pdf)
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
    
    out.pdf <- paste("output/", analysis.name, "-kras-meta-gene", ".pdf", sep="")
    pdf(out.pdf)
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

## e.g., variable = CMS
do.fit.and.qc <- function(df, st, variable, analysis.name, file.variable) {
  lm.obj <- lm(as.formula(paste("expr", variable, sep=" ~ ")), data=df)
  lm.sum <- summary(lm.obj)
  sum.file <- paste("output/", st, "-", analysis.name, "-", file.variable, "-sum.tsv", sep="")
  capture.output(lm.sum, file = sum.file)

  diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
  pdf(diag.file)
  ##  print(plot(lm.obj))
  plot(lm.obj)
  d <- dev.off()
        
  lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
  coeffs <- as.data.frame(coefficients(lm.sum))
  coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
  write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

}

## This modified from
##        library(devtools)
##        source_gist("524eade46135f6348140", filename="ggplot_smooth_func.R")
## to inclue pvalue
stat_smooth_func <- function(mapping = NULL, data = NULL,
                        geom = "smooth", position = "identity",
                        ...,
                        method = "auto",
                        formula = y ~ x,
                        se = TRUE,
                        n = 80,
                        span = 0.75,
                        fullrange = FALSE,
                        level = 0.95,
                        method.args = list(),
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE,
                        xpos = NULL,
                        ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,
                      
                      setup_params = function(data, params) {
                        # Figure out what type of smoothing to do: loess for small datasets,
                        # gam with a cubic regression basis for large data
                        # This is based on the size of the _largest_ group.
                        if (identical(params$method, "auto")) {
                          max_group <- max(table(data$group))
                          
                          if (max_group < 1000) {
                            params$method <- "loess"
                          } else {
                            params$method <- "gam"
                            params$formula <- y ~ s(x, bs = "cs")
                          }
                        }
                        if (identical(params$method, "gam")) {
                          params$method <- mgcv::gam
                        }
                        
                        params
                      },
                      
                      compute_group = function(data, scales, method = "auto", formula = y~x,
                                               se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                               xseq = NULL, level = 0.95, method.args = list(),
                                               na.rm = FALSE, xpos=NULL, ypos=NULL) {
                        if (length(unique(data$x)) < 2) {
                          # Not enough data to perform fit
                          return(data.frame())
                        }
                        
                        if (is.null(data$weight)) data$weight <- 1
                        
                        if (is.null(xseq)) {
                          if (is.integer(data$x)) {
                            if (fullrange) {
                              xseq <- scales$x$dimension()
                            } else {
                              xseq <- sort(unique(data$x))
                            }
                          } else {
                            if (fullrange) {
                              range <- scales$x$dimension()
                            } else {
                              range <- range(data$x, na.rm = TRUE)
                            }
                            xseq <- seq(range[1], range[2], length.out = n)
                          }
                        }
                        # Special case span because it's the most commonly used model argument
                        if (identical(method, "loess")) {
                          method.args$span <- span
                        }
                        
                        if (is.character(method)) method <- match.fun(method)
                        
                        base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                        model <- do.call(method, c(base.args, method.args))
                        
                        m = model
                        f <- summary(m)$fstatistic

                        ##                            ct <- cor.test(formula = y ~ x, data = data)
                        ct <- cor.test(x = data$x, y = data$y)
##                        print(ct)
                        corr <- ct$estimate
                        pval <- ct$p.value
                        stat <- ct$statistic
                        cat(paste0("t = ", stat, " df = ", ct$parameter, " nrow = ", nrow(na.omit(data[,c("x","y")])), " cor = ", corr, " pval = ", pval, "\n"))

                        ## This includes eqn of line and r2
                        ####                        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pval, list(a = format(coef(m)[1], digits = 1), b = format(coef(m)[2], digits = 1), pval = format(pf(f[1],f[2],f[3],lower.tail=F), digits = 1), r2 = format(summary(m)$r.squared, digits = 1)))
                        ## eq <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~pval, list(a = format(coef(m)[1], digits = 1), b = format(coef(m)[2], digits = 1), pval = format(pf(f[1],f[2],f[3],lower.tail=F), digits = 1), r2 = format(summary(m)$r.squared, digits = 1)))
                        eq <- substitute(italic(r)~"="~rval*","~~italic(p)~"="~pval, list(pval = format(pval, digits = 1), rval = format(corr, digits = 1)))
                        func_string = as.character(as.expression(eq))
                        if(is.null(xpos)) xpos = min(data$x)*0.9
                        if(is.null(ypos)) ypos = max(data$y)*0.9
                        data.frame(x=xpos, y=ypos, label=func_string)
                        
                      },
                      
                      required_aes = c("x", "y")
)


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

    
    ## gene set prep
    env <- new.env()
    load(paste0(synapse.repo.dir, "/markersG.RData"),envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
    cat(paste("Number of rows analyzed: ", nrow(expr.m), "\n", sep=""))
    cat(paste("Number of columns analyzed: ", ncol(expr.m), "\n", sep=""))

##    save(clin.m, file="clin.m.Rd")
##    save(expr.m, file="expr.m.Rd")    
    
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
    pdf(out.pdf)
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
    exclude.no.lbl <- FALSE
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
    sigs1 <- c("wnt", "tak1", "wnt", "pos.ras.dep", "pos.ras.dep", "neg.ras.dep", "pos.ras.dep", "neg.ras.dep", "CIRC")
    sigs2 <- c("CIRC", "CIRC", "myc.sig", "neg.ras.dep", "CIRC", "CIRC", "tak1", "tak1", "neoantigens")
    sig1.names <- c("wnt", "tak1", "wnt", "pos.ras.dep", "pos.ras.dep", "neg.ras.dep", "pos.ras.dep", "neg.ras.dep", "CIRC Enrichment Score")
    sig2.names <- c("CIRC", "CIRC", "myc.sig", "neg.ras.dep", "CIRC", "CIRC", "tak1", "tak1", "Log10 ( # Neoantigens + 1 ) ")
    for(i in 1:length(sigs1)) {
        sig1 <- sigs1[i]
        sig2 <- sigs2[i]
        ## par(mfrow=c(2,5))

        if(!(sig1 %in% rownames(esn)) || !(sig2 %in% rownames(esn))) { next }

        if(FALSE) {
        for(cms.lbl in c("CMS1", "CMS2", "CMS3", "CMS4")) {
            out.pdf <- paste("output/", "scatter-", analysis.name, "-", sig1, "-vs-", sig2, "-", cms.lbl, ".pdf", sep="")
            pdf(out.pdf)
            flag <- cms == cms.lbl & !is.na(cms) & !is.na(kras.label)
            df <- data.frame(x = esn[sig2, flag], y = esn[sig1, flag], KRAS = paste0("KRAS ", kras.label[flag]))
            
            g <- ggplot(data = df, aes(x = x, y = y, label = y))
            g <- g + stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse=TRUE)
            g <- g + geom_smooth(method = "lm", se = FALSE)
            g <- g + geom_point()
            g <- g + facet_wrap(~ KRAS)
            g <- g + ylab(sig1.names[i])
            g <- g + xlab(sig2.names[i])
            print(g)
            d <- dev.off()
        }

        ## Here, facet'ing by KRAS.  
        cms.lbl <- "ALL"
        out.pdf <- paste("output/", "scatter-", analysis.name, "-", sig1, "-vs-", sig2, "-", cms.lbl, ".pdf", sep="")
        pdf(out.pdf)
        flag <- !is.na(kras.label)
        df <- data.frame(x = esn[sig2, flag], y = esn[sig1, flag], KRAS = paste0("KRAS ", kras.label[flag]))

        g <- ggplot(data = df, aes(x = x, y = y, label = y))
        g <- g + stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse=TRUE)
        g <- g + geom_smooth(method = "lm", se = FALSE)
        g <- g + geom_point()
        g <- g + facet_wrap(~ KRAS)
        g <- g + ylab(sig1.names[i])
        g <- g + xlab(sig2.names[i])
        print(g)
        d <- dev.off()
        }
        
        ## Here, don't facet by KRAS.
        cat(paste0(analysis.name, " ", sig1, " vs ", sig2, ": "))

        out.pdf <- paste("output/", "scatter-", analysis.name, "-", sig1, "-vs-", sig2, ".pdf", sep="")
        pdf(out.pdf)
        df <- data.frame(x = esn[sig2,], y = esn[sig1,])

        g <- ggplot(data = df, aes(x = x, y = y, label = y))
        g <- g + stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse=TRUE)
        g <- g + geom_smooth(method = "lm", se = FALSE)
        g <- g + geom_point()
        g <- g + ylab(sig1.names[i])
        g <- g + xlab(sig2.names[i])
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
        if(FALSE) {        
            if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
                df$status <- msi.factor
            }
            if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
                df$site <- site.factor
            }
            if(exclude.no.lbl) {        
                df <- df[df$CMS != "NOLBL",]                        
            }
        }

        lm.obj <- lm(as.formula(paste("expr", "KRAS", sep=" ~ ")), data=df)
        lm.sum <- summary(lm.obj)
        sum.file <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], "-sum.tsv", sep="")
        capture.output(lm.sum, file = sum.file)

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
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
        if(FALSE) {
            if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
                df$status <- msi.factor
            }
            if(!all(is.na(site)) && !(length(unique(site[!is.na(site)])) == 1)) {
                df$site <- site.factor
            }
            if(exclude.no.lbl) {        
                df <- df[df$CMS != "NOLBL",]                        
            }
        }

        
##        kras.label <- kras.code
        ##    df <- data.frame(expr=esn.st, KRAS=factor(kras.label), CMS=rep("ALL", length(esn.st)))
##        df <- data.frame(expr=esn.st, KRAS=factor(kras.codon), CMS=rep("ALL", length(esn.st)))
        
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], ".pdf", sep="")        
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))        
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
        p <- p + annotate("text", x = 1.5, y = mt.wt.max + 4 * yoffset, label = pval.to.text(pval))
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
    ## Do the same analysis, but coarse grained.
    cms.coarse.tbl <- cms.tbl

    for(sti in 1:num.biomarkers) {

        ## st will be the gene/gene set to plot
        st <- to.plot[sti]

        ## The y axis label for this gene/gene set
        ylab <- ylabels[sti]
        
##        cat(paste("Computing ", st, " ~ ", "CMS", "\n", sep=""))

        ## Subset the data to the gene/gene set of interest
        esn.st <- as.vector(esn[st,])

        ## Get the KRAS mutation status of each sample
        kras.label <- sapply(kras, function(x) ifelse(x==1, kras.states[1], kras.states[2]))

        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms, levels=cms2.first.levels))
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        
        do.fit.and.qc(df, st, "CMS", analysis.name, "CMS")
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
        flag <- cms %in% c("CMS1", "CMS3", "CMS4", "NOLBL")
        cms.coarse[flag] <- "CMSOther"
        flag <- cms %in% c("CMS2")
        cms.coarse[flag] <- "CMS2"
        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms.coarse, levels=c("CMS2", "CMSOther", "NOLBL")))
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }

        do.fit.and.qc(df, st, "CMS", analysis.name, "CMS")
        
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
    
    ## Do the same analysis, but coarse grained.
    cms.coarse <- cms
    flag <- cms %in% c("CMS1", "CMS3", "CMS4", "NOLBL")
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
        ##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))        
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
        flag <- cms %in% c("CMS1", "CMS3", "CMS4", "NOLBL")
        cms.coarse[flag] <- "CMSOther"
        flag <- cms %in% c("CMS2")
        cms.coarse[flag] <- "CMS2"
        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(cms.coarse, levels=c("CMS2", "CMSOther", "NOLBL")))
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
        out.pdf <- paste("output/", st, "-", analysis.name, "-cms-coarse", ".pdf", sep="")
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))        
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
        panel.maxs <- sapply(c("CMS2", "CMSOther"), function(cmsLbl) {
            mask <- df$CMS == cmsLbl
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(panel.maxs) <- c("CMS2", "CMSOther")
        
        ## Draw the error bars from CMS2 to CMS1, 3, and 4.
        cms.str <- "CMSOther"
        ## Use raw pvalue for univariate analysis        
        ## pval.col <- paste0(cms.str, ".apval")
        pval.col <- paste0(cms.str, ".", pval.suffix)
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

    ## Analysis 2b: biomarker ~ cms mt (only)
    cms.mt.tbl <- data.frame(variable=to.plot)
    for(prefix in cms2.first.levels[-1]) {
        if(!(prefix %in% cms.levels)) { next }
        pval.col <- paste0(prefix, ".pval")
        cms.mt.tbl[,pval.col] <- rep(NA, nrow(cms.mt.tbl))
    }
    rownames(cms.mt.tbl) <- to.plot
    ## Do the same analysis, but coarse grained.
    cms.mt.coarse.tbl <- cms.mt.tbl

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

        ## Do a coarse grained analysis
        cms.coarse <- cms
        flag <- cms %in% c("CMS1", "CMS3", "CMS4", "NOLBL")
        cms.coarse[flag] <- "CMSOther"
        flag <- cms %in% c("CMS2")
        cms.coarse[flag] <- "CMS2"
        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(paste0(cms.coarse, " MT"), levels=paste0(c("CMS2", "CMSOther", "NOLBL"), " MT")))
        if(exclude.no.lbl) {
            df <- df[df$CMS != "NOLBL",]                        
        }
        df <- df[!is.na(df$KRAS) & (kras.label == "MT"),]
        
        do.fit.and.qc(df, st, "CMS", analysis.name, "CMS-mt")
        
        cms.str <- "CMSOther"
        pval.col <- paste0(cms.str, ".pval")
        pval.flag <- grepl(x=coeffs.2$coefficient, cms.str)
            ## Set this below using wilcox        
##        cms.coarse.tbl[sti,pval.col] <- coeffs.2[pval.flag, 5]
        
    }
    cat("\n")
    
    ## Adjust the pvals for the cms analysis
##    for(cmsindx in c(1, 3, 4)) {
##        cms.str <- paste0("CMS", cmsindx)
##        pval.col <- paste0(cms.str, ".pval")
##        padj.col <- paste0(cms.str, ".apval")        
##        cms.mt.tbl[,padj.col] <- p.adjust(cms.mt.tbl[,pval.col], method="BH")
##    }
##    pval.col <- paste0("CMSOther", ".pval")
##    padj.col <- paste0("CMSOther", ".apval")        
##    cms.mt.coarse.tbl[,padj.col] <- p.adjust(cms.mt.coarse.tbl[,pval.col], method="BH")        

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
    
    ## Do the same analysis, but coarse grained.
    cms.coarse <- cms
    flag <- cms %in% c("CMS1", "CMS3", "CMS4", "NOLBL")
    cms.coarse[flag] <- "CMSOther"
    flag <- cms %in% c("CMS2")
    cms.coarse[flag] <- "CMS2"
    flag <- ( kras.label == "MT" ) & ( ( cms.coarse == "CMS2" ) | ( cms.coarse == "CMSOther" ) )
    cms.coarse.flag <- cms.coarse[flag]
    pval <- apply(esn[to.plot,flag,drop=F], 1, function(x) wilcox.test(x ~ cms.coarse.flag)$p.value)
    a <- p.adjust(pval,method="BH")
    b <- apply(esn[to.plot,flag,drop=F], 1, function(x) (mean(x[cms.coarse.flag == "CMS2"], na.rm=TRUE) - mean(x[cms.coarse.flag == "CMSOther"], na.rm=TRUE)) > 0)
    
    cms.mt.coarse.tbl <- data.frame(variable=to.plot, CMSOther.pval=pval, CMSOther.apval=a, CMSOther.cms2Up=b)
    rownames(cms.mt.coarse.tbl) <- to.plot

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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
        flag <- cms %in% c("CMS1", "CMS3", "CMS4", "NOLBL")
        cms.coarse[flag] <- "CMSOther"
        flag <- cms %in% c("CMS2")
        cms.coarse[flag] <- "CMS2"
        df <- data.frame(expr=esn.st, KRAS=factor(kras.label, levels=c(kras.states[1], kras.states[2])), CMS=factor(paste0(cms.coarse, " MT"), levels=paste0(c("CMS2", "CMSOther", "NOLBL"), " MT")))
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
        out.pdf <- paste("output/", st, "-", analysis.name, "-cms-mt-coarse", ".pdf", sep="")
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
        panel.maxs <- sapply(c("CMS2 MT", "CMSOther MT"), function(cmsLbl) {
            mask <- df$CMS == cmsLbl
            m <- max(df$expr[mask], na.rm=TRUE)
            m
        })
        names(panel.maxs) <- c("CMS2", "CMSOther")

        ## Draw the error bars from CMS2 to CMS1, 3, and 4.
        cms.str <- "CMSOther"
        ## Use raw pvalue for univariate analysis        
        ## pval.col <- paste0(cms.str, ".apval")
        pval.col <- paste0(cms.str, ".", pval.suffix)
        pval <- as.numeric(cms.mt.coarse.tbl[sti, pval.col])
        
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

    write.table(file=paste0("output/", analysis.name, "-cms-mt-tbl.xls"), cms.mt.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(file=paste0("output/", analysis.name, "-cms-mt-tbl-coarse.xls"), cms.mt.coarse.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
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
        pdf(diag.file)
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
        ## And we will add 7 error bars.  So adjust the y axis
        ## to accommodate this shift.
        ylimits[2] <- ylimits[2] + 7 * num.error.bars * yoffset
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
                    tmp <- draw.err.bar(st, kras.cms.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix)
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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
            tmp <- draw.err.bar(st, interaction.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix)
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df.transformed, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {            
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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

            tmp <- draw.err.bar(st, transformed.interaction.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix)
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
        pdf(diag.file)
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
        pdf(diag.file)
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
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

            tmp <- draw.err.bar(st, mt.vs.wt.tbl, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix = pval.suffix)
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
        pdf(forest.file)
        p <- suppressWarnings(forest_model(lm.obj))
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
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
        pdf(diag.file)
        plot(lm.obj)
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file)
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

            ## BAD here
            
            ## Create a pdf for this gene/gene set (with facets)
            out.pdf <- paste("output/", st, "-", analysis.name, "-", formula, ".pdf", sep="")
            pdf(out.pdf)

            ## Create the plot
            formula.indx <- which(colnames(df) == formula)[1]
            expr.indx <- which(colnames(df) == "expr")[1]
            formula.col <- colnames(df)[formula.indx]
            expr.col <- "expr"
            p <- ggplot(data=df, aes_string(x = formula.col, y = expr.col))
            p <- p + ylab(ylab)

            ## Create a box plot where the x axis is CMS ...
            p <- p + geom_boxplot(aes_string(fill = formula.col))
            ##            p <- p + geom_jitter()
            ## Do not annotate MSI/MSS when we are comparing expr to MSI status
            p <- p + guides(fill = guide_legend(order = 1))
            nxt.order <- 2
            if(("status" %in% colnames(df)) && ("site" %in% colnames(df)) && (formula != "site") && (formula != "status")) {
                p <- p + geom_beeswarm(aes(colour = status, shape = site))
                p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
                p <- p + guides(colour = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                
                p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
                p <- p + guides(shape = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1
            } else if(("status" %in% colnames(df)) && (formula != "status")) {
                p <- p + geom_beeswarm(aes(colour = status))
                p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
                p <- p + guides(colour = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                                
            } else if(("site" %in% colnames(df)) && (formula != "site")) {
                p <- p + geom_beeswarm(aes(shape = site))
                p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
                p <- p + guides(shape = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                                
            } else {
                p <- p + geom_beeswarm()
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

    ret.list <- list(esn=esn, kras.tbl=kras.tbl, cms.tbl=cms.tbl, kras.cms.tbl=kras.cms.tbl, kras.cms.interaction.tbl=interaction.tbl, kras.cms.no.interaction.tbl=no.interaction.tbl, kras.mt.vs.wt.tbl=mt.vs.wt.tbl, site.tbl=formula.tbls[[1]])
    if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
        ret.list[["msi.tbl"]] <- formula.tbls[[2]]
    }
    return(ret.list)

    
}

## Draw an error bar between two facets (that may or may not be the same).
draw.err.bar <- function(st, tbl2, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix) {

    ## Get the corresponding adjusted pvalue comparing expression of st (a gene or gene set)
    ## across the two conditions specified by cmsLbl1 and cmsLbl2
    ## Plot raw pvalues!!!
    ## adj.pval.col <- paste(toString(cmpLbl), ".apval", sep="")
    adj.pval.col <- paste(toString(cmpLbl), ".", pval.suffix, sep="")
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
    load(paste0(synapse.repo.dir, "/markersG.RData"),envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
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
    neoantigens <- clin.m$neoantigens
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
##    cms.levels <- unique(cms)
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

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file)
        p <- suppressWarnings(forest_model(lm.obj))
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
        plot(lm.obj)
        d <- dev.off()
        
        lm.file <- gsub(x = sum.file, pattern="-sum", replacement="-lm")
        coeffs <- as.data.frame(coefficients(lm.sum))
        coeffs.2 <- as.data.frame(cbind(coefficient=rownames(coeffs), coeffs))
        write.table(file=lm.file, coeffs.2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        pdf.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-facet.tsv", sep="")        
        pdf(pdf.file)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
##        p <- p + geom_jitter()
        if(("status" %in% colnames(df)) && ("site" %in% colnames(df))) {
            p <- p + geom_beeswarm(aes(colour = status, shape = site))
        } else if("status" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(colour = status))
        } else if("site" %in% colnames(df)) {
            p <- p + geom_beeswarm(aes(shape = site))
        } else {
            p <- p + geom_beeswarm()
        }
        p <- p + guides(fill = guide_legend(order = 1))                
        if("status" %in% colnames(df)) {        
            p <- p + scale_colour_manual(values = c("MSI" = "blue", "MSS" = "black"))
            p <- p + guides(colour = guide_legend(order = 2))
        }        
        if("site" %in% colnames(df)) {        
            p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
            p <- p + guides(shape = guide_legend(order = 3))            
        }        

        ## ... faceted on CMS label.
        p <- p + facet_grid(~CMS)

        d <- dev.off()
        
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
        pdf(diag.file)
        plot(lm.obj)
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file)
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
    load(paste0(synapse.repo.dir, "/markersG.RData"),envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
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

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file)
        p <- suppressWarnings(forest_model(lm.obj))
        print(p)
        d <- dev.off()

        diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
        pdf(diag.file)
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

        pdf.file <- paste("output/", st, "-", analysis.name, "-", formula.name, "-facet.pdf", sep="")        
        cat(paste0("Creating pdf ", pdf.file, "\n"))
        pdf(pdf.file)

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
        plot(lm.obj)
        d <- dev.off()

        forest.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-forest.pdf")
        pdf(forest.file)
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
    load(paste0(synapse.repo.dir, "/markersG.RData"),envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
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
    load(paste0(synapse.repo.dir, "/markersG.RData"),envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
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
    load(paste0(synapse.repo.dir, "/markersG.RData"),envir=env)
    myeloid <- na.omit(unlist(mget(x=env$markersG$Myeloid, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    monocytic <- na.omit(unlist(mget(x=env$markersG$Monocyte_derived, envir=org.Hs.egSYMBOL2EG,ifnotfound=NA)))
    
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

do.msi.mut.tbl <- function(expr, clin, file, mut.tbl = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification", num.clusters = 2) {

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
    pdf(paste0(file, ".pdf"))
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    d <- dev.off()
    
    ## tyms different color
    return(clin.orig)
}


do.msi <- function(expr, clin, file, neoantigen.cnts = NULL, sample.col = "Tumor_Sample_Barcode", gene.col = "SYMBOL", variant.type.col = "Variant_Classification", num.clusters = 2) {

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
    pdf(paste0("output/", file, ".pdf"))
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
    d <- dev.off()

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
            p <- p + ylab("TYMS Expression")

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
                cat(paste0("Beeswarm ", facet, ": ", "TYMS", " vs ", "Status", " W = ", res$statistic, " p = ", pval, " ", paste(unlist(lapply(unique(df.facet[,"Status"]), function(var) { paste0("n_", var, " = ", length(which(df.facet[,"Status"] == var))) })), collapse=", "), "\n"))
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
                tmp <- draw.err.bar("tyms", pval.tbl, ranges, panel.maxs, g, lbl1, statuses[1], lbl2, statuses[2], indx1, indx2, cmpLbl, states, pval.suffix = "pval")
                g <- tmp[["g"]]
                panel.maxs <- tmp[["panel.maxs"]]
            }
            
            ## Turn clipping off to see the line across the panels
            g$layout$clip <- "off"

            pdf(paste0("output/", file, "-tyms-vs-msi.pdf"))
            grid.draw(g)
            d <- dev.off()
            
        } else {
            df <- merge(tyms.tbl, clin.orig, by="sample", all=FALSE)
            df <- df[,c("msi.inferred", "tyms")]
            colnames(df) <- c("Inferred.Status", "TYMS")
            df$Inferred.Status <- factor(toupper(df$Inferred.Status))
            pdf(paste0("output/", file, "-tyms-vs-annotated-msi.pdf"))
            plt <- plot.beeswarm.with.err(df, x.name = "Inferred.Status", y.name = "TYMS", x.lab = "Inferred Status", y.lab = "TYMS Expression")
            print(plt)
            d <- dev.off()
        }

        ## Plot TYMS vs annotated MSI
    }

    return(clin.orig)
}


do.ras.dependency <- function(expr, clin, ras.dependent.signature, file) {

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

    num.annotation.legends <- 0
    has.kras <- any(!is.na(clin.mask$kras))
    if(has.kras) {
        num.annotation.legends <- num.annotation.legends + 1
    }
    num.cols <- num.annotation.legends

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
    ras.dependent.common <- intersect(ras.dependent.signature, rownames(expr.mask))

    expr.mask <- expr.mask[ras.dependent.common,]

    method <- "ward.D"
    method <- "complete"
    dd.row <- as.dendrogram(hclust(dist(expr.mask), method=method))
    row.ord <- order.dendrogram(dd.row)

    hc.col <- hclust(dist(t(expr.mask)), method=method)
    dd.col <- as.dendrogram(hc.col)
    col.ord <- order.dendrogram(dd.col)

    ## Find the separation between the MSI and MSS cases
    num.clusters <- 2
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

    mutation.cnt.by.sample <- NULL

    ## Re-order the expression matrix based on the clustering
    expr.mask <- (expr.mask[row.ord, col.ord])
    expr.mask.with.tyms <- expr.mask.with.tyms[, col.ord]    
    clin.mask <- clin.mask[col.ord, ]
    
    remove.labels <- TRUE
    remove.kras.labels <- TRUE

    if(!remove.kras.labels) {
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
    heat.label <- textGrob("Genes in RAS Dependency Signature", gp=gpar(fontsize=fs), just=heat.just)

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

    g.legend.kras <- NULL
    if(has.kras) {
        cat("Has kras\n")        
        kras.tbl <- data.frame(sample=clin.mask$sample, kras=clin.mask$kras)
        kras.tbl$sample <- factor(kras.tbl$sample, levels=levels, ordered=TRUE)
        kras.tbl$kras <- as.character(kras.tbl$kras)
        kras.tbl$kras[kras.tbl$kras == "1"] <- "MT"
        kras.tbl$kras[kras.tbl$kras == "0"] <- "WT"
        
        kras.melt <- melt(kras.tbl, id.vars="sample")
        p.kras <- ggplot(data = kras.melt, aes(x = variable, y=sample))
        p.kras <- p.kras + geom_tile(aes(fill = value), show.legend = TRUE)
        p.kras <- p.kras + scale_fill_manual(values = c("black", "gray"), na.value = "white")    
        p.kras <- p.kras + guides(fill = guide_legend(title="kras", title.hjust = 0.5, title.theme=element_text(face = 'italic', angle = 0)))
        p.kras <- p.kras + theme(legend.text = element_text(size = fs))        
        p.kras <- p.kras + geom_hline(yintercept = cross.over.sample.indx, linetype = 'dashed')    
        g.kras <- ggplotGrob(p.kras)
        g.legend.kras <- gtable::gtable_filter(g.kras, "guide-box")
        g.kras <- gtable::gtable_filter(g.kras, "panel")

        if(!remove.kras.labels) {
            g.kras <- ggplotGrob(p.kras)
        }
        
        kras.label <- textGrob("KRAS", gp=gpar(fontsize=fs, fontface='italic'), just=just, rot=90)
        
        title.row <- c(title.row, indx)
        grobs[[indx]] <- kras.label
        indx <- indx + 1
        
        plot.row <- c(plot.row, indx)
        grobs[[indx]] <- g.kras
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
##        row <- spacer.row
##        row[num.cols] <- indx
        layout <- rbind(layout, row)
        grobs[[indx]] <- g.legend.kras
        indx <- indx + 1
        widths <- c(widths, annotation.width)        
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
    pdf(paste0(file, ".pdf"))
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

## TODO/TO SHOW (all-comers and mss)
## out.pdf <- paste("output/", st, "-", analysis.name, "-kras-cms-faceted", ".pdf", sep="")
## univariate -cms.pdf
## univeriate -kras.pdf  <- problem
## forest
## -no-interaction
## look at fraction of CMS2 (and CMS3) in tcga and kfs

## cp circ-tcga-*kras-cms-faceted.pdf afigs/
## cp circ-kfs-*kras-cms-faceted.pdf afigs/
## cp circ-tcga-all-cms.pdf afigs/
## cp circ-tcga-mss-cms.pdf afigs/
## cp circ-tcga-mss-msi-inferred-cms.pdf afigs/
## cp circ-kfs-all-cms.pdf afigs/
## cp circ-kfs-mss-msi-inferred-cms.pdf afigs/

## cp circ-tcga-all-kras-MT-vs-WT.pdf afigs/
## cp circ-tcga-mss-kras-MT-vs-WT.pdf afigs/
## cp circ-tcga-mss-msi-inferred-kras-MT-vs-WT.pdf afigs/

## cp circ-kfs-all-kras-MT-vs-WT.pdf afigs/
## cp circ-kfs-mss-msi-inferred-kras-MT-vs-WT.pdf afigs/


plot.pval.heatmap <- function(tbl, x, y, value, xlab, ylab, use.italics = TRUE, show.values = TRUE, log.transform = TRUE) {
    g <- ggplot(data = tbl, aes_string(x = x, y = y, fill = value))
    g <- g + geom_tile(color="white")
    
    ##    g <- g + scale_fill_gradient(low = "blue", high = "white", labels = c(0, 0.05), space = "Lab", name = "p-value", breaks = c(0, 0.05))

    g <- g + theme_minimal()
    g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
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
    if(log.transform) {
        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,"white"), values = scales::rescale(c(0, min(s))), limit = c(min(s),0), space="Lab", name="p-value", labels = rev(labels[s <= 0]), breaks = s[s <= 0], guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Upregulated, -Log[10]~italic('p')))))
    } else {
        g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg,"white"), values = scales::rescale(c(min(s), min(s) + ct, 0), from = c(min(s), 0)), limit = c(min(s),0), space="Lab", labels = c(0, ct, 1), breaks = c(min(s), min(s) + ct, 0), name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Upregulated, italic('q')-value))))
    }
    non.na.tbl <- na.omit(tbl)
    if(log.transform) {
        non.na.tbl$value <- round(10^-abs(non.na.tbl$value), digits=2)
    } else {
        non.na.tbl$value <- round(-(abs(non.na.tbl$value)-1), digits=2)
    }
    print(non.na.tbl)
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

    if(log.transform) {
        print(rev(s[s>=0]))
        print(rev(labels[s>=0]))
        
        g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, "white"), values = scales::rescale(c(max(s), 0)), limit = c(0,max(s)), space="Lab", name="p-value2", labels = (labels[s >= 0]), breaks = (s[s >= 0]), guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Downregulated, -Log[10]~italic('p')))))
    } else {
        g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, lightpos,"white"), values = scales::rescale(c(0, ct, 1), from = c(0, 1)), labels = c(0, ct, 1), breaks = c(0, ct, 1), limit = c(0,1), space="Lab", name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Downregulated, italic('q')-value))))
    }

    gt2 <- ggplotGrob(g2)
    
    g.legend2 <- gtable::gtable_filter(gt2, "guide-box")    
    
    
##    g <- g + scale_colour_gradient(limits=c(0,5))
    ##    g
    ##    grid.arrange(arrangeGrob(gt, g.legend1, g.legend1, nrow=1))
    grobs <- list(gt, g.legend1, g.legend2)
    matrix <- rbind(c(1,2),c(1,3))
    grid.arrange(grobs = grobs, layout_matrix = matrix)
}

if(FALSE) {
        ##    g <- g + scale_fill_gradientn(na.value="black", colours = c(lightneg,"white", "white", "white", lightpos), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide=guide_colorbar(title=expression(-Log[10]~italic('p'))), na.value="black")
    ##    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(lightneg,"white", "white", "white", lightpos), values = scales::rescale(c(min(s), -1, 0, 1, max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide=guide_colorbar(title=expression(-Log[10]~italic('p'))), na.value="black")
    ##    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(lightneg,"white", lightpos), values = scales::rescale(c(min(s), 0, max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide=guide_colorbar(title=expression(-Log[10]~italic('p'))), na.value="black")
###    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(lightneg,"white", lightpos), values = scales::rescale(c(min(s), 0, max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg,"white", lightpos, darkpos), values = scales::rescale(c(min(s), min(s) + ct, 0, max(s) - ct, max(s)), from = c(min(s), max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
    ##    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(lightneg,"white"), values = scales::rescale(c(min(s), 0)), limit = c(min(s),0), space="Lab", name="p-value", labels = labels[s <= 0], breaks = s[s <= 0], guide=guide_colorbar(title=expression(-Log[10]~italic('p'))), na.value="black")
    lines <- list(bquote("first"),bquote("second"))
    ##    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(lightneg,"white"), values = scales::rescale(c(min(s), 0)), limit = c(min(s),0), space="Lab", name="p-value", labels = labels[s <= 0], breaks = s[s <= 0], guide=guide_colorbar(title=expression(italic('KRAS')~upregulated)), na.value="black")
##    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg, lightneg,"white"), limit = c(min(s), min(s) + ct, 0), space="Lab", name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Upregulated, italic('q')-value))))
    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg,"white"), values = scales::rescale(c(min(s), min(s) + ct, 0), from = c(min(s), 0)), limit = c(min(s),0), space="Lab", labels = c(0, ct, 1), breaks = c(min(s), min(s) + ct, 0), name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Upregulated, italic('q')-value))))
    non.na.tbl <- na.omit(tbl)
    if(log.transform) {
        non.na.tbl$value <- round(10^-abs(non.na.tbl$value), digits=2)
    } else {
        non.na.tbl$value <- round(-(abs(non.na.tbl$value)-1), digits=2)
    }
    print(non.na.tbl)
    if(show.values) {
        g <- g + geom_text(data = non.na.tbl, aes_string(x = x, y = y, label = value))
    }

    gt <- ggplotGrob(g)
    g.legend1 <- gtable::gtable_filter(gt, "guide-box")

###    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(lightneg,"white", lightpos), values = scales::rescale(c(min(s), 0, max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
    g <- g + scale_fill_gradientn(na.value = na.value, colours = c(darkneg,lightneg, "white", lightpos, darkpos), values = scales::rescale(c(min(s), min(s) + ct, 0, max(s) - ct , max(s)), from = c(min(s), max(s))), limit = c(min(s),max(s)), space="Lab", name="p-value", labels = labels, breaks = s, guide = "none")
    gt <- ggplotGrob(g)

    g2 <- ggplot(data = tbl, aes_string(x = x, y = y, fill = value))
    g2 <- g2 + geom_tile(color="white")

    ##    g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(lightpos,"white"), limit = c(0,1), space="Lab", name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Downregulated, italic('q')-value))))
    g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c(darkpos, lightpos,"white"), values = scales::rescale(c(0, ct, 1), from = c(0, 1)), labels = c(0, ct, 1), breaks = c(0, ct, 1), limit = c(0,1), space="Lab", name="upq-value", guide=guide_colorbar(title=bquote(atop(italic('KRAS')~Downregulated, italic('q')-value))))
    
##     g2 <- g2 + scale_fill_gradientn(na.value = na.value, colours = c("white",lightpos), values = scales::rescale(c(0, max(s))), limit = c(0,max(s)), space="Lab", name="p-value2", labels = labels[s >= 0], breaks = s[s >= 0], guide=guide_colorbar(title=expression(-Log[10]~italic('p2'))), na.value="black")
    gt2 <- ggplotGrob(g2)
    
    g.legend2 <- gtable::gtable_filter(gt2, "guide-box")    
    
    
##    g <- g + scale_colour_gradient(limits=c(0,5))
    ##    g
    ##    grid.arrange(arrangeGrob(gt, g.legend1, g.legend1, nrow=1))
    grobs <- list(gt, g.legend1, g.legend2)
    matrix <- rbind(c(1,2),c(1,3))
    grid.arrange(grobs = grobs, layout_matrix = matrix)

}
