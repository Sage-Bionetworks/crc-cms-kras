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
suppressPackageStartupMessages(library("gridExtra"))
library("forestmodel")

cat("Examine metabolic pathways/genes\n")

# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

num.processes <- 1

synapse.repo.dir <- "../input/"

set.seed(1234)

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

text.size <- 25
star.size <- 12
ns.size <- 6
gp.star.size <- 25
gp.ns.size <- 12

## Immune sets that Justin looked at.
immune.sets <- c("IMMUNE_ESTIMATE", "IMMUNE_RESP_GO_BP", "PD1_REACTOME", "IMMUNE_CD8MACRO_GALON", "IMMUNE_TH1_GALON", "IMMUNE_NKC_BREAST", "IMMUNE_THF_BREAST", "IMMUNE_TH17_GOUNARI", "IMMUNE_TREG_LUCAS", "IMMUNE_MDSC_ALBELDA")

system("mkdir output/")



## Draw an error bar between two facets (that may or may not be the same).
draw.err.bar <- function(st, tbl2, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix, yoffset) {

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
##    ymax <- max(ranges[[cmsIndx1]][["y.range"]]) - min(ranges[[cmsIndx1]][["y.range"]])
##    yoffset <- 0.01 * ymax
    
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

    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)
    
    ## Finally, draw the vertical line on the "right side" of the error bar--
    ## this goes from the error bar to the maximum of the MT KRAS
    ## expression values in the second facet/CMS 2.
    start <- c(data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    end <- c(data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]),
             data2npc(m + 1 * yoffset, ranges[[cmsIndx2]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)

    ## Update the maximum values used within each facet to reflect the
    ## newly added error bars
    panel.maxs[paste0(cmsLbl1,"-",kras.states[1])] <- m + 2 * 6 * yoffset
    panel.maxs[paste0(cmsLbl2,"-",kras.states[1])] <- m + 2 * 6 * yoffset     
    panel.maxs[paste0(cmsLbl1,"-",kras.states[2])] <- m + 2 * 6 * yoffset
    panel.maxs[paste0(cmsLbl2,"-",kras.states[2])] <- m + 2 * 6 * yoffset     

    ## Add the asterisk designation of the pvalue.
    text <- pval.to.text(pval)

    ## Position the text in the middle of the error bar.
    xrange <- ranges[[cmsIndx1]][["x.range"]]
    ## I believe this 3 comes from the fact that the high range is 2.6 and I was
    ## padding for the space between the facets
    xrange[2] <- xrange[2] + 3 * abs(cmsIndx2 - cmsIndx1)
    ## pos <- c(0.5 * ( data2npc(x1.index,ranges[[cmsIndx1]][["x.range"]]) + data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]) ), data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    num.dy <- 2
    sz <- gp.star.size
    if(text == "n.s.") {
        num.dy <- 4
        sz <- gp.ns.size
    }
    pos <- c(data2npc(0.5 * ( x1.index + 3 * abs(cmsIndx2 - cmsIndx1) + x2.index), xrange),
             data2npc(m + num.dy * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    ## pos <- c(data2npc(0.5, xrange), data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    

    delta <- runif(n=1, min=10^-5, max=10^-3)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2], gp=gpar(fontsize = sz), vjust = 0), t = t + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
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
    num.dy <- 2
    sz <- gp.star.size
    if(text == "n.s.") {
        num.dy <- 4
        sz <- gp.ns.size
    }
    pos <- c(data2npc(1.5, ranges[[cmsIndx]][["x.range"]]),
             data2npc(m + num.dy * yoffset, ranges[[cmsIndx]][["y.range"]]))

    delta <- runif(n=1, min=10^-5, max=10^-3)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2], gp=gpar(fontsize = sz), vjust = 0), t = t + delta, l = l1, b = 4, r = l2)
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
    num.dy <- 2
    sz <- gp.star.size
    if(text == "n.s.") {
        num.dy <- 4
        sz <- gp.ns.size
    }
    pos <- c(data2npc(1 + 0.5 * ( 1 + 3 * ( cmsIndx2 - cmsIndx1 ) - 1 ), xrange),
             data2npc(m + num.dy * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    

    delta <- runif(n=1, min=10^-5, max=10^-3)
    print(text)
    print(pos)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2], gp=gpar(fontsize = sz), vjust = 0), t = t + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
}


plot.beeswarm.with.err <- function(df, x.name, y.name, x.lab = NULL, y.lab = NULL) {
    if(is.null(x.lab)) { x.lab <- x.name }
    if(is.null(y.lab)) { y.lab <- y.name }
    df <- df[!is.na(df[,x.name]) & !is.na(df[,y.name]),]
    p <- ggplot(data=df, aes_string(x=x.name, y=y.name))
##    p <- p + ggtitle(paste0(y.name, " vs ", x.name))
    p <- p + ylab(y.lab)
    p <- p + xlab(x.lab)
    p <- p + theme(text = element_text(size = text.size))
    
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
    text <- pval.to.text(pval)
    num.dy <- 1
    sz <- star.size
    if(text == "n.s.") {
        num.dy <- 3
        sz <- ns.size
    }
    p <- p + annotate("text", x = 1.5, y = max12 + num.dy * yoffset, size = sz, label = text, vjust = 0)
    p
}


## For an expression data set, expr, plot the expression of a gene/gene set as a function
## of CMS cluster and KRAS mutation status.  Do such independently for each gene/gene set
## specified in to.plot and ylabels.
doAnalysis <- function(expr, clin, analysis.name, to.plot=c(), ylabels=c(), kras.status.field="kras", kras.states=c("MT","WT"), only.do.kras.analysis = FALSE, plot.adjusted.pvalues = FALSE, main = NULL) {

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
    pdf(out.pdf)
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
        pdf(out.pdf)
        df <- data.frame(x = esn[sig2,], y = esn[sig1,])

        if(!all(is.na(msi)) && !(length(unique(msi[!is.na(msi)])) == 1)) {
            df$status <- msi.factor
        }
        
        g <- ggplot(data = df, aes(x = x, y = y, label = y))
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
        
##        kras.label <- kras.code
        ##    df <- data.frame(expr=esn.st, KRAS=factor(kras.label), CMS=rep("ALL", length(esn.st)))
##        df <- data.frame(expr=esn.st, KRAS=factor(kras.codon), CMS=rep("ALL", length(esn.st)))
        
        
        ## Create a PDF for this gene/gene set (with facets)
        out.pdf <- paste("output/", st, "-", analysis.name, "-kras-", kras.states[1], "-vs-", kras.states[2], ".pdf", sep="")        
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=KRAS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + theme(text = element_text(size = text.size))                
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
        p <- p + theme(text = element_text(size = text.size))
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df, aes(x=CMS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is CMS ...
        p <- p + geom_boxplot(aes(fill=CMS))
        p <- p + theme(text = element_text(size = text.size))        
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + theme(text = element_text(size = text.size))        
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + theme(text = element_text(size = text.size))        
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
        pdf(out.pdf)

        ## Create the plot
        p <- ggplot(data=df.transformed, aes(x=KRAS, y=expr))
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + theme(text = element_text(size = text.size))        
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
        if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }        
        p <- p + ylab(ylab)

        ## Create a box plot where the x axis is KRAS mutation status ...
        p <- p + geom_boxplot(aes(fill=KRAS))
        p <- p + theme(text = element_text(size = text.size))        
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
            if(!is.null(main)) { p <- p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) }            
            p <- p + theme(text = element_text(size = text.size))            
            p <- p + ylab(ylab)

            ## Create a box plot where the x axis is CMS ...
            p <- p + geom_boxplot(aes_string(fill = formula.col))
            ##            p <- p + geom_jitter()
            ## Do not annotate MSI/MSS when we are comparing expr to MSI status
            p <- p + guides(fill = guide_legend(order = 1))
            nxt.order <- 2
            if(("status" %in% colnames(df)) && ("site" %in% colnames(df)) && (formula != "site") && (formula != "status")) {
                p <- p + geom_beeswarm(aes(colour = status, shape = site))
                p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
                p <- p + guides(colour = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1                
                p <- p + scale_shape_manual(values = c("left" = 15, "right" = 16, "rectum" = 17), na.value = 8)
                p <- p + guides(shape = guide_legend(order = nxt.order))
                nxt.order <- nxt.order + 1
            } else if(("status" %in% colnames(df)) && (formula != "status")) {
                p <- p + geom_beeswarm(aes(colour = status))
                p <- p + scale_colour_manual(na.value = "gray", values = c("MSI" = "blue", "MSS" = "black"))
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


inverse.norm.transform <- function(x) {
##    p <- 2*pnorm(abs(x), lower.tail=FALSE) 
##    x2 <- qnorm(p/2, lower.tail=FALSE)*sign(x)
    ##    x2
    qq <- qqnorm(x, plot.it = FALSE)
    trn <- ( mean(x, na.rm=TRUE) + ( sd(x, na.rm=TRUE) * qq$x ) )
    trn
}

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
            p <- p + theme(text = element_text(size = text.size))            
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

cat("Logging in\n")
synapseLogin("brian.white")
cat("Logged in\n")

## Read in the clinical annotation data

cat("Reading clinical data\n")
obj <- synGet(id="syn8533554", downloadFile = TRUE, downloadLocation = synapse.repo.dir)
clin.data.file <- getFileLocation(obj)
clin <- read.table(clin.data.file, sep=",", header=TRUE, as.is=TRUE)

## Read in the CMS clusters
cat("Reading CMS results\n")
obj <- synGet(id="syn4978511", downloadFile = TRUE, downloadLocation = synapse.repo.dir)
cms.file <- getFileLocation(obj)
cms <- read.table(cms.file, header=TRUE, as.is=TRUE)

## Read in the gene-set definitions used in Justin's Nat Med paper
obj <- synGet(id="syn2321865", downloadFile = TRUE, downloadLocation = synapse.repo.dir)
file <- getFileLocation(obj)
gene.sets <- GSA.read.gmt(file)

cat("Reading bindea\n")
synId <- "syn8533557"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = synapse.repo.dir)
file <- getFileLocation(obj)
bindeaTbl <- read.table(file,sep='\t',header=TRUE)

## Restrict to a few cell types of interest
bindeaTbl <- bindeaTbl[bindeaTbl$CellType %in% c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells"),]

bindeaCellTypes <- unique(bindeaTbl$CellType)
bindeaGsets <- lapply(bindeaCellTypes, function(x){
    return(as.character(unique(bindeaTbl$Symbol[bindeaTbl$CellType==x])))
})
names(bindeaGsets) <- bindeaCellTypes
bindeaGsets <- bindeaGsets[sapply(bindeaGsets, length) > 10]
cat("Reading CIRC signature\n")
synId <- "syn8533558"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = synapse.repo.dir)
file <- getFileLocation(obj)
circ.genes <- read.table(file,as.is=TRUE)[,1]
bindeaGsets[["CIRC"]] <- circ.genes

bindeaGsets[["mek"]] <- c("DUSP6","PHLDA1","SPRY2","DUSP4","ETV4")
## This TAK1 signature was read off of Fig 4A of Singh et al (2012) Cell
tak1.signature <- c("RBP1", "SEMA3A", "SYT1", "EMR2", "PROX1", "INHBB", "ABHD2", "C1orf116", "SNTB1", "TAF9B", "PRF1", "SLC2A1", "GAD1", "MSX2", "PELI2", "ITGB4", "C21orf96", "GPR56", "PDK3", "GLS", "ACSL1", "BIK", "RUNX1", "SYK", "RGL1", "NAV2", "FYN", "HSPA12A", "MBOAT2", "BAMBI", "BMP7", "GGH")
bindeaGsets[["tak1"]] <- tak1.signature

## MDSC signature from Angelova
mdsc.sig <- c("ADORA3", "AG2", "ARG1", "BIN2", "C1orf162", "CAPS", "CD117", "CD11B", "CD11C", "CD124", "CD14", "CD15", "CD163L1", "CD1D1", "CD21", "CD23", "CD274", "CD31", "CD35", "CD40", "CD43", "CD44", "CD66B", "COX2", "CR2", "CTR9", "EBP", "FAM48A", "FAM70B", "FCER2", "FCGRT", "FERMT3", "FLOT1", "FLT1", "GIMAP7", "GLI4", "GNA15", "GPR34", "GPSM3", "HLA-DR", "IDO", "IKZF1", "IL12", "IL13", "IL18BP", "IL1R", "IL4RA", "INPP5D", "ITGA3", "KDR", "KRIT1", "LGALS3", "MGAT4A", "NAIP", "NEK3", "NFSF13", "NOG", "PARVG", "PDRG1", "PECAM1", "PIK3R5", "PPP1R2P4", "PSAP", "PTGES2", "PTPRE", "RNASE1", "RP11", "S100A8", "S100A9", "SELPLG", "SLA", "SLC36A1", "SLC44A1", "ST8SIA4", "STAT3", "STAT6", "TBXAS1", "TFGFB1", "TFGFB2", "TFGFB3", "TFGFB5", "TFRC", "TGFB2", "TPP1", "VTCN1")
bindeaGsets[["MDSC"]] <- mdsc.sig

nfkb.pathway <- gene.sets$genesets[[which(gene.sets$geneset.names == "NFKB_BIOCARTA")]]
bindeaGsets[["NFKB"]] <- nfkb.pathway

## Add Glutathione and Glutamine pathways examined in Nat Medicine paper--
## these are from KEGG and GO BP, respectively (see online methods)
bindeaGsets[["Glutathione"]] <- gene.sets$genesets[[which(gene.sets$geneset.names == "GLUTATHIONE_KEGG")]]

bindeaGsets[["Glutamine"]] <- gene.sets$genesets[[which(gene.sets$geneset.names == "GLUTAMINE_GO_BP")]]


## Myc targets from PMID 14519204, as used by J. Guinney in his Nat Med pub. 
## myc.targets <- c("APEX", "CAD", "CCNA2", "CCND2", "CCNE1", "CDK4", "CDKN1A", "CDKN2B", "CHC1", "DDX18", "DUSP1", "EIF4E", "ENO1", "FASN", "FKBP4", "FN1", "GADD45A", "HSPA4", "HSPCAL3", "HSPD1", "HSPE1", "LDHA", "MGST1", "MYC", "NCL", "NME1", "NME2", "NPM1", "ODC1", "PPAT", "PTMA", "RPL23", "RPL3", "RPL6", "RPS15A", "SRM", "TERT", "TFRC", "THBS1", "TNFSF6", "TP53", "TPM1")
myc.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "MYC_TARGETS_ZELLER")]]
bindeaGsets[["myc.sig"]] <- myc.targets


## Wnt target genes from PMID 17320548, as used by J. Guinney in his Nat Med pub.  
## wnt.targets <- c("ASCL2", "AXIN2", "BMP4", "C1orf33", "HIG2", "HSPC111", "KITLG", "LGR5", "MYC", "NOL1", "PPIF", "SOX4", "WRD71", "ZIC2", "ZNRF3")
wnt.targets <- gene.sets$genesets[[which(gene.sets$geneset.names == "WNT_FLIER")]]
bindeaGsets[["wnt"]] <- wnt.targets

for(immune.set in immune.sets) {
    bindeaGsets[[immune.set]] <- gene.sets$genesets[[which(gene.sets$geneset.names == immune.set)]]
}

immune.set.gene.symbols <- unique(as.vector(unlist(lapply(bindeaGsets[c(immune.sets,"CIRC")], function(x) as.vector(x)))))


## Match the clinical annotations and the CMS clusters
idxs <- match(clin$sample, cms$sample)
clin <- clin[!is.na(idxs),]
clin$cms_label <- cms$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples[na.omit(idxs)]

cat("Reading TCIA Neoantigens\n")
synId <- "syn8533559"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = synapse.repo.dir)
file <- getFileLocation(obj)
neoantigen.tbl <- read.table(file, sep="\t", header=TRUE)

## Count the number of neoantigens
## NB: this allows multiple per gene, which happens frequently.
neoantigen.cnts <- as.data.frame(table(neoantigen.tbl$patientBarcode))
clin$neoantigens <- rep(NA, nrow(clin))
neoantigen.df <- data.frame(sample = neoantigen.cnts$Var1, neoantigens=log10(unname(neoantigen.cnts$Freq) + 1))

ns <- neoantigen.cnts$Var1
neoantigen.cnts <- neoantigen.cnts$Freq
names(neoantigen.cnts) <- ns

idxs <- match(clin$sample, neoantigen.df$sample)
clin$neoantigens[!is.na(idxs)] <- neoantigen.df$neoantigens[na.omit(idxs)]

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


doaffy <- function(synId){
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = synapse.repo.dir)
    file <- getFileLocation(obj)
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    cat(paste0("doaffy reading file: ", file, "\n"))
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

data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]

pval.to.text <- function(pval) {
    if(is.na(pval)) { return("n.s.") }
    if(pval >= 0.05) {
        return("n.s.")
    } else if(pval < 0.0001) {
        return("****")
    } else if(pval < 0.001) {
        return("***")
    } else if(pval < 0.01) {
        return("**")
    } else if(pval < 0.05) {
        return("*")
    }
    die("Should not be here\n")
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

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  ##  r <- cor(x, y)
  ct <- cor.test(x, y)
  r <- ct$estimate

  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r = ", txt, sep = "")
  text(0.5, 0.6, txt)

  # p-value calculation
  ## p <- cor.test(x, y)$p.value
  p <- ct$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  if(p < 0.01) {
      txt2 <- format(p, digits = digits, scientific = TRUE)
  }
  txt2 <- paste("p = ", txt2, sep = "")
##  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

do.tcga.analysis <- TRUE
do.kfs.analysis <- TRUE

tcga_expr <- NULL
kfsyscc_expr <- NULL
if(do.tcga.analysis || do.kfs.analysis) {
  cat("Reading TCGA data\n")
  ## TCGA data
  tcga_expr <- {
      obj <- synGet(id="syn2325328", downloadFile = TRUE, downloadLocation = synapse.repo.dir)
      tcga.file <- getFileLocation(obj)
      
      tcga_expr <- read.table(tcga.file, sep="\t", header=TRUE, check.names=FALSE)
  
    expr <- as.matrix(tcga_expr[,-1])
    rownames(expr) <- tcga_expr[,1]
    expr
  }

  median.expr <- unlist(apply(tcga_expr, 1, median))
  names(median.expr) <- rownames(tcga_expr)
  ## CIITA median expr is in the 46% percentile.
  length(which(median.expr <= median.expr["CIITA"]))/length(median.expr)

  
  ## Append the counts to the expression matrix
  tcga_expr <- rbind(tcga_expr, neoantigens=rep(NA, ncol(tcga_expr)))
  neoantigen.cnts <- neoantigen.cnts[names(neoantigen.cnts) %in% colnames(tcga_expr)]
  tcga_expr["neoantigens",names(neoantigen.cnts)] <- unname(neoantigen.cnts)
  tcga_expr["neoantigens",names(neoantigen.cnts)] <- log10(tcga_expr["neoantigens",names(neoantigen.cnts)] + 1)
  
  cat("Reading KFS data\n")
  kfsyscc_expr <- doaffy("syn2363564")
  ## kfsyscc_expr <- doaffy(read.table(opt$`kfs-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  colnames(kfsyscc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(kfsyscc_expr))

  ## Look at the relative expression of CIITA
  median.expr <- unlist(apply(kfsyscc_expr, 1, median))
  names(median.expr) <- rownames(kfsyscc_expr)
  ## CIITA median expr is in the 5% percentile.
  length(which(median.expr <= median.expr["CIITA"]))/length(median.expr)
  
}

bindeaGsets[["Glutathione"]] <- gene.sets$genesets[[which(gene.sets$geneset.names == "GLUTATHIONE_KEGG")]]

bindeaGsets[["Glutamine"]] <- gene.sets$genesets[[which(gene.sets$geneset.names == "GLUTAMINE_GO_BP")]]

to.plot <- c("Glutathione", "Glutamine", "HK2", "LDHA", "MGLL", "ANGPTL4", "VEGFA", "ACSL5", "PCK2", "AGPAT7", "SLCO1B3")

to.plot <- c("Glutathione", "Glutamine", "HK2", "LDHA", "MGLL", "ANGPTL4", "VEGFA", "ACSL5", "PCK2", "LPCAT4", "SLCO1B3")



ylabels <- unlist(lapply(to.plot, function(x) paste0(x, " Expression")))

ylabels[1:2] <- c("Glutathione Enrichment Score", "Glutamine Enrichment Score")

if(do.tcga.analysis) {
    
    cat("Doing TCGA analysis\n")

    res <- doAnalysis(tcga_expr, clin, "tcga-all", to.plot=to.plot, ylabels=ylabels, main = "TCGA")

    clin.tcga.msi.inferred <- do.msi(tcga_expr, clin, "tcga-msi", neoantigen.cnts = neoantigen.cnts)

    tcga.samples <- clin.tcga.msi.inferred$dataset == "tcga"

    cat("Overlap of MSI inferred and annotated in TCGA:\n")
    print(table(clin.tcga.msi.inferred$msi[tcga.samples], clin.tcga.msi.inferred$msi.inferred[tcga.samples]))
    
    clin.tcga.msi.inferred$msi <- clin.tcga.msi.inferred$msi.inferred
    ## Inferred MSI
##    res <- doAnalysis(tcga_expr, clin.tcga.msi.inferred, "tcga-all-msi-inferred", to.plot=to.plot, ylabels=ylabels, main = "TCGA")
    
    ## Restrict to (inferred) MSS cases
    clin.mss.inferred <- clin.tcga.msi.inferred[!is.na(clin.tcga.msi.inferred$msi) & (clin.tcga.msi.inferred$msi == "mss"),]
##    res <- doAnalysis(tcga_expr, clin.mss.inferred, "tcga-mss-msi-inferred", to.plot=to.plot, ylabels=ylabels, main = "TCGA")    

    ## How many MSI cases are CMS1?
    cat("Overlap of MSI and CMS1 in TCGA inferred:\n")
    print(table(clin.tcga.msi.inferred$msi, clin.tcga.msi.inferred$cms_label))
    
    ## Restrict to (inferred) MSS cases and exclude CMS1
    clin.mss.inferred <- clin.tcga.msi.inferred[!is.na(clin.tcga.msi.inferred$msi) & (clin.tcga.msi.inferred$msi == "mss") & !is.na(clin.tcga.msi.inferred$cms_label) & (clin.tcga.msi.inferred$cms_label != "CMS1"),]
##    res <- doAnalysis(tcga_expr, clin.mss.inferred, "tcga-mss-msi-inferred-no-cms", to.plot=to.plot, ylabels=ylabels, main = "TCGA")    
    
    
    ## Restrict to MSS cases
    clin.mss <- clin[!is.na(clin$msi) & (clin$msi == "mss"),]
    res <- doAnalysis(tcga_expr, clin.mss, "tcga-mss-msi-annotated", to.plot=to.plot, ylabels=ylabels, main = "TCGA")    

    
    ## How many MSI cases are CMS1?
    cat("Overlap of MSI and CMS1 in TCGA annotated:\n")
    print(table(clin$msi, clin$cms_label))
    
    ## Restrict to MSS annotated cases and exclude CMS1
    clin.mss.and.no.cms1 <- clin[!is.na(clin$msi) & (clin$msi == "mss") & !is.na(clin$cms_label) & (clin$cms_label != "CMS1"), ]
    res <- doAnalysis(tcga_expr, clin.mss.and.no.cms1, "tcga-mss-msi-annotated-no-cms1", to.plot=to.plot, ylabels=ylabels, main = "TCGA")        
    
##    dset <- "tcga"
##    for(tbl in c("kras.tbl", "cms.tbl", "kras.cms.tbl", "kras.cms.interaction.tbl", "kras.cms.no.interaction.tbl", "kras.mt.vs.wt.tbl")) {
##        tbl.underscore <- gsub(tbl, pattern="\\.", replacement="_")
##        write.table(res[[tbl]], file=paste0("Middleton_", dset, "_", tbl.underscore, ".xls"), sep="\t", quote=FALSE, row.names=FALSE)
##    }
}


if(do.kfs.analysis) {

    cat("Doing KFS analysis\n")
    res <- doAnalysis(kfsyscc_expr, clin, "kfs-all", to.plot=to.plot, ylabels=ylabels, main = "KFSYSCC")
    dset <- "kfs"

    clin.kfs.msi.inferred <- do.msi(kfsyscc_expr, clin, "kfs-msi")
    clin.kfs.msi.inferred$msi <- clin.kfs.msi.inferred$msi.inferred
    res <- doAnalysis(kfsyscc_expr, clin.kfs.msi.inferred, "kfs-all-msi-inferred", to.plot=to.plot, ylabels=ylabels, main = "KFSYSCC")

    ## Restrict to (inferred) MSS cases
    clin.mss <- clin.kfs.msi.inferred[!is.na(clin.kfs.msi.inferred$msi) & (clin.kfs.msi.inferred$msi == "mss"),]
    ## res <- doAnalysis(kfsyscc_expr, clin.mss, "kfs-mss-msi-inferred", to.plot=to.plot, ylabels=ylabels, main = "KFSYSCC")    

    ## How many MSI cases are CMS1?
    cat("Overlap of MSI and CMS1 in KFS inferred:\n")
    kfs.samples <- clin.kfs.msi.inferred$dataset == "kfs"
    print(table(clin.kfs.msi.inferred$msi[kfs.samples], clin.kfs.msi.inferred$cms_label[kfs.samples]))
    
    ## Restrict to (inferred) MSS cases and exclude CMS1
    clin.mss <- clin.kfs.msi.inferred[!is.na(clin.kfs.msi.inferred$msi) & (clin.kfs.msi.inferred$msi == "mss") & !is.na(clin.kfs.msi.inferred$cms_label) & (clin.kfs.msi.inferred$cms_label != "CMS1"),]
    res <- doAnalysis(kfsyscc_expr, clin.mss, "kfs-mss-msi-inferred-no-cms1", to.plot=to.plot, ylabels=ylabels, main = "KFSYSCC")
    

}

## Look at the overlap in the immune signatures ...
immune.sigs <- list()
cat("Reading CIRC signature again\n")
immune.sigs[["CIRC"]] <- circ.genes
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

print(sessionInfo())

save.image(".Rdata")


