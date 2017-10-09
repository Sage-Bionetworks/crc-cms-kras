
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
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("boot"))
suppressPackageStartupMessages(library("cowplot"))
library("forestmodel")

# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("devtools"))
if(!("MCPcounter" %in% .packages(all.available=TRUE))) {
    devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")
}
suppressPackageStartupMessages(library("MCPcounter"))
              
num.processes <- 1

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapse.repo.dir <- "../input/"

set.seed(1234)

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
draw.err.bar <- function(st, tbl2, ranges, panel.maxs, g, cmsLbl1, kras.state1, cmsLbl2, kras.state2, cmsIndx1, cmsIndx2, cmpLbl, kras.states, pval.suffix, yoffset, plot.pvals.as.stars = TRUE, stat.name = "p") {

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
    text <- pval.to.text(pval, plot.pvals.as.stars = plot.pvals.as.stars, stat.name = stat.name)

    ## Position the text in the middle of the error bar.
    xrange <- ranges[[cmsIndx1]][["x.range"]]
    ## I believe this 3 comes from the fact that the high range is 2.6 and I was
    ## padding for the space between the facets
    xrange[2] <- xrange[2] + 3 * abs(cmsIndx2 - cmsIndx1)
    ## pos <- c(0.5 * ( data2npc(x1.index,ranges[[cmsIndx1]][["x.range"]]) + data2npc(x2.index,ranges[[cmsIndx2]][["x.range"]]) ), data2npc(m + 4 * yoffset, ranges[[cmsIndx1]][["y.range"]]))
    num.dy <- 2
    sz <- gp.star.size
    if((text == "n.s.") || (plot.pvals.as.stars)) {
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


plot.correlation <- function(x, y, x.name, y.name) {
    g <- ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y))
    g <- g + geom_point()
    g <- g + stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse=TRUE)
    g <- g + geom_smooth(method = "lm", se = FALSE)
    g <- g + xlab(x.name)
    g <- g + ylab(y.name)
    g
}

           
## End doanalysis
source("doanalysis.R")

inverse.norm.transform <- function(x) {
##    p <- 2*pnorm(abs(x), lower.tail=FALSE) 
##    x2 <- qnorm(p/2, lower.tail=FALSE)*sign(x)
    ##    x2
    qq <- qqnorm(x, plot.it = FALSE)
    trn <- ( mean(x, na.rm=TRUE) + ( sd(x, na.rm=TRUE) * qq$x ) )
    trn
}


## e.g., variable = CMS
do.fit.and.qc <- function(df, st, variable, analysis.name, file.variable) {
  lm.obj <- lm(as.formula(paste("expr", variable, sep=" ~ ")), data=df)
  lm.sum <- summary(lm.obj)
  sum.file <- paste("output/", st, "-", analysis.name, "-", file.variable, "-sum.tsv", sep="")
  capture.output(lm.sum, file = sum.file)

  diag.file <- gsub(x = sum.file, pattern="-sum.tsv", replacement="-diag.pdf")
  pdf(diag.file, useDingbats = FALSE)
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
immune.populations <- c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells")

immune.population.labels <- unlist(lapply(immune.populations, function(str) paste0(gsub(str, pattern="cells", replacement="Cell"), " Enrichment Score")))
bindeaTbl <- bindeaTbl[bindeaTbl$CellType %in% immune.populations,]

circ.signature <- "CIRC"
circ.signature.label <- "CIRC Enrichment Score"

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

file <- "../input/c2.cp.biocarta.v6.0.symbols.gmt"
gsea.gene.sets <- GSA.read.gmt(file)
for(set.name in gsea.gene.sets$geneset.names) {
    bindeaGsets[[set.name]] <- gsea.gene.sets$genesets[[which(gsea.gene.sets$geneset.names == set.name)]]
}

file <- "../input/c2.cp.kegg.v6.0.symbols.gmt"
gsea.gene.sets <- GSA.read.gmt(file)
for(set.name in gsea.gene.sets$geneset.names) {
    bindeaGsets[[set.name]] <- gsea.gene.sets$genesets[[which(gsea.gene.sets$geneset.names == set.name)]]
}

file <- "../input/c2.cp.reactome.v6.0.symbols.gmt"
gsea.gene.sets <- GSA.read.gmt(file)
for(set.name in gsea.gene.sets$geneset.names) {
    bindeaGsets[[set.name]] <- gsea.gene.sets$genesets[[which(gsea.gene.sets$geneset.names == set.name)]]
}

file <- "../input/h.all.v6.0.symbols.gmt"
gsea.gene.sets <- GSA.read.gmt(file)
for(set.name in gsea.gene.sets$geneset.names) {
    bindeaGsets[[set.name]] <- gsea.gene.sets$genesets[[which(gsea.gene.sets$geneset.names == set.name)]]
}




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
combine_probes_2_gene_svd <- function(expr, genes){
  
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

combine_rows <- function(matrx, map, fun = "mean") {
    map <- subset(map, from %in% rownames(matrx))
    rownames(map) <- map$from
    matrx <- matrx[rownames(matrx) %in% map$from, ]
    map <- map[rownames(matrx),]
    suppressPackageStartupMessages(library("Matrix.utils"))    
    matrx <- as.matrix((aggregate.Matrix((matrx), groupings=list(map$to), fun = fun)))
    matrx
}


combine_probes_2_gene <- function(expr, genes, method="svd"){
    if(method == "svd") {
        cat("Combining probes using SVD\n")        
        return(combine_probes_2_gene_svd(expr, genes))
    } else if(method == "mean") {
        cat("Combining probes using mean\n")
        map <- data.frame(from = rownames(expr), to = genes)
##        save(list = ls(all.names = TRUE, envir=environment()), file = ".RData.combine", envir = environment())
        return(combine_rows(expr, map, fun = method))
    } else {
        cat(paste0("combine_probes_2_gene does implement: ", method, "\n"))
        q(status=-1)
    }
}

# wget https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-101/A-AFFY-101.adf.txt

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))


synapse.fread <- function(synId) {
    obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = synapse.repo.dir)
    file <- getFileLocation(obj)
  # data.set <- synId
    # data.set <- as.matrix(data.set)
    ## data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    ## A little nonsense because fread does not seem to deal well with
    ## rows
    cat(paste0("doaffy reading file: ", file, "\n"))
    use.fread <- FALSE
    data.set <- NULL
    if(use.fread) {
      ## I don't believe this works
      data.set <- fread(file, header=TRUE, fill=TRUE)
      df <- as.data.frame(data.set)
      rownames(df) <- df[,1]
      cols <- colnames(df)[1:(ncol(df)-1)]
      df <- df[,-1]
      colnames(df) <- cols
      data.set <- df
      rm(df)
    } else {
        data.set <- read.table(file, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
        cat("Head of data set\n")
        print(head(data.set[, c(1:5, ((ncol(data.set)-5):ncol(data.set)))]))
    }
    colnames(data.set) <- gsub("(.*?)_.*","\\1",colnames(data.set))
    data.set
}

convert_hg133_probes_to_genes <- function(data.set, method = "svd") {
    symbol <- unlist(mget(rownames(data.set), envir=hgu133plus2SYMBOL, ifnotfound=NA))
    mask <- !is.na(symbol)
    data.set.m <- data.set[mask,]
    symbol.m <- symbol[mask]
    expr <- combine_probes_2_gene(data.set.m, symbol.m, method = method)
    expr
}

doaffy <- function(synId, method = "svd"){
    data.set <- synapse.fread(synId)
    expr <- convert_hg133_probes_to_genes(data.set, method = method)
    expr
}

data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]

pval.to.text <- function(pval, plot.pvals.as.stars = TRUE, stat.name = "p") {
    if(is.na(pval)) { return("n.s.") }
    if(!plot.pvals.as.stars) {
        eq <- substitute(italic(var)~"="~pval, list(var = stat.name, pval = ifelse(pval < 0.001, format(pval, digits = 1, scientific=TRUE), format(pval, digits = 1))))
        return(eq)
##        return(as.character(signif(pval, digits=2)))
    }
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

## Append the neoantigen counts to the expression matrix
tcga_expr <- rbind(tcga_expr, neoantigens=rep(NA, ncol(tcga_expr)))
neoantigen.cnts <- neoantigen.cnts[names(neoantigen.cnts) %in% colnames(tcga_expr)]
tcga_expr["neoantigens",names(neoantigen.cnts)] <- unname(neoantigen.cnts)
tcga_expr["neoantigens",names(neoantigen.cnts)] <- log10(tcga_expr["neoantigens",names(neoantigen.cnts)] + 1)

cat("Reading KFS data\n")
kfsyscc_microarray_expr <- synapse.fread("syn2363564")
colnames(kfsyscc_microarray_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(kfsyscc_microarray_expr))  
kfsyscc_expr <- doaffy("syn2363564", method = "mean")
## kfsyscc_expr <- doaffy(read.table(opt$`kfs-expr-file`, sep="\t", header=TRUE, row.names=1, as.is=TRUE))
colnames(kfsyscc_expr) <- gsub("(.*?)\\.CEL","\\1", colnames(kfsyscc_expr))

to.plot <- c("CIRC", "HIF1A", "BHLHE40", "CIITA", "wnt", "myc.sig", "neoantigens")
ylabels <- c("CIRC Enrichment Score", "HIF1A Expression", "BHLHE40 Expression", "CIITA Expression", "Wnt Signature Enrichment Score", "Myc Signature Enrichment Score", "Log10 ( # Neoantigens + 1 ) ")

to.plot <- c("CIRC", "CIITA", "wnt", "myc.sig", "neoantigens", "BATF3")
ylabels <- c("CIRC Enrichment Score", "CIITA Expression", "Wnt Signature Enrichment Score", "Myc Signature Enrichment Score", "Log10 ( # Neoantigens + 1 ) ", "BATF3 Expression")

## to.plot <- c("CIRC", "CIITA", "wnt", "myc.sig", "neoantigens", "BATF3", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION")
## ylabels <- c("CIRC Enrichment Score", "CIITA Expression", "Wnt Signature Enrichment Score", "Myc Signature Enrichment Score", "Log10 ( # Neoantigens + 1 ) ", "BATF3 Expression", "HALLMARK_IL6_JAK_STAT3_SIGNALING\nEnrichment Score", "HALLMARK_INTERFERON_GAMMA_RESPONSE\nEnrichment Score", "HALLMARK_COMPLEMENT\nEnrichment Score", "HALLMARK_INFLAMMATORY_RESPONSE\nEnrichment Score", "HALLMARK_IL2_STAT5_SIGNALING\nEnrichment Score", "HALLMARK_ALLOGRAFT_REJECTION\nEnrichment Score")

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
        intersection <- intersect(immune.sigs[[i]], immune.sigs[[j]])
        m[i,j] <- length(intersection)
        if(length(intersection) > 0) {
        }
    }
}
immune.sig.names <- names(immune.sigs)
immune.sig.names <- gsub(x=immune.sig.names, pattern="_", replacement="\\\n")
rownames(m) <- unlist(immune.sig.names)
colnames(m) <- unlist(immune.sig.names)

library(gridExtra)
library(grid)
pdf("immune-sig-overlap.pdf", useDingbats = FALSE)
#grid.table(m, base_size = 5)
grid.newpage()
g3 <- tableGrob(m, theme = ttheme_default(base_size = 7, colhead=list(fg_params=list(rot=90))), rows=rownames(m), cols=colnames(m))
grid.draw(g3)
d <- dev.off()


## BEGIN

## IFN signatures in Figure 2
## Yes B2 is those signatures but the reviewer seems to want to see the HALLMARK GSEA IFN signature..

## HALLMARK
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/h.all.v6.0.symbols.gmt
## HALLMARK_INTERFERON_GAMMA_RESPONSE

## BIOCARTA
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.biocarta.v6.0.symbols.gmt

## KEGG
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.kegg.v6.0.symbols.gmt

## REACTOME
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.reactome.v6.0.symbols.gmt
## REACTOME_INTERFERON_GAMMA_SIGNALING
## REACTOME_INTERFERON_SIGNALING


## Print the overlap in the CIRC, and IFNG data sets
set1 <- "CIRC"
set2 <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
set1.genes <- bindeaGsets[[set1]]
set2.genes <- bindeaGsets[[set2]]
cat(paste0("Overlap between ", set1, " (n = ", length(unique(set1.genes)), ") and ", set2, " (n = ", length(unique(set2.genes)), ") = ", length(unique(intersect(set1.genes, set2.genes))), "\n"))

set1 <- "CIRC"
set2 <- "REACTOME_INTERFERON_GAMMA_SIGNALING"
set1.genes <- bindeaGsets[[set1]]
set2.genes <- bindeaGsets[[set2]]
cat(paste0("Overlap between ", set1, " (n = ", length(unique(set1.genes)), ") and ", set2, " (n = ", length(unique(set2.genes)), ") = ", length(unique(intersect(set1.genes, set2.genes))), "\n"))


to.plot <- c("CIRC", "neoantigens", "BATF3", "THBD", "CCL4", "ATF3", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_INTERFERON_GAMMA_SIGNALING", "IFNG")
ylabels <- c("CIRC Enrichment Score", "Log10 ( # Neoantigens + 1 ) ", "BATF3 Expression", "THBD Expression", "CCL4 Expression", "ATF3 Expression", "Hallmark IFNG Enrichment Score", "Reactome IFNG Enrichment Score", "IFNG Expression")

hallmark.kras.cms.res <- list()
immune.kras.cms.res <- list()
circ.kras.cms.res <- list()
circ.kras.res <- list()
hallmark.kras.res <- list()
bindea.kras.res <- list()
gsva.es.lst <- list()
clin.es.lst <- list()
mcp.lst <- list()

## BEGIN ANALYSIS

hallmark.gene.sets <- names(bindeaGsets)[grepl(pattern="HALLMARK", names(bindeaGsets))]
## hallmark.gene.sets <- c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_NOTCH_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_PANCREAS_BETA_CELLS")
immune.gene.sets <- c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells")

gsva.es.lst <- list()
clin.es.lst <- list()
mcp.populations <- c()
mcp.population.labels <- c()

other.genes <- c("STAT1", "CIITA", "CXCL10", "IRF1", "IL6")
other.labels <- unlist(lapply(other.genes, function(str) paste0(str, " Expression")))
hallmark.ylabels <- unlist(lapply(hallmark.gene.sets, function(str) paste0(str, "\nEnrichment Score")))

for(dataset in c("tcga", "kfs")) {
    ds <- tcga_expr
    if(dataset == "kfs") {
        ds <- kfsyscc_expr
    }
    tmp <- doGSVA(ds, clin, dataset, to.plot=unique(c(circ.signature, hallmark.gene.sets, immune.gene.sets, other.genes)))
    es <- tmp[["es"]]
    clin.es <- tmp[["clin"]]

    mcp <- NULL
    if(dataset == "kfs") {
        mcp <- MCPcounter.estimate(kfsyscc_microarray_expr, featuresType="affy133P2_probesets")
    } else {
        mcp <- MCPcounter.estimate(tcga_expr, featuresType="HUGO_symbols")        
    }
    mcp.populations <- rownames(mcp)
    mcp.population.labels <- unlist(lapply(mcp.populations, function(x) paste0(x, " Population Percentage")))
    ## rownames(mcp) <- unlist(lapply(rownames(mcp), function(x) paste0("MCP ", x)))

    common.cols <- intersect(colnames(mcp), colnames(es))
    es <- as.matrix(rbind(mcp[, common.cols], es[, common.cols]))

    gsva.es.lst[[dataset]] <- es
    clin.es.lst[[dataset]] <- clin.es
}

## Add neoantigens to tcga
clin <- clin.es.lst[["tcga"]]
rownames(clin) <- clin$sample
es <- gsva.es.lst[["tcga"]]
cols <- intersect(colnames(es), rownames(clin))
es <- rbind(es, neoantigens=rep(NA, ncol(es)))
es["neoantigens",cols] <- clin[cols,"neoantigens"]
gsva.es.lst[["tcga"]] <- es

## Infer MSI status of KFS (where it is not annotated) and TCGA (where it is)
clin.inferred.msi.lst <- list()
tmp <- do.msi(tcga_expr, clin.es.lst[["tcga"]], "tcga-msi", neoantigen.cnts = neoantigen.cnts)
clin.inferred.msi.lst[["tcga"]] <- tmp[["clin"]]
tcga.msi.figa <- tmp[["msi.g"]]
tcga.msi.figb <- tmp[["tyms.msi.annotated.g"]]

tmp <- do.msi(kfsyscc_expr, clin.es.lst[["kfs"]], "kfs-msi")
clin.inferred.msi.lst[["kfs"]] <- tmp[["clin"]]
kfs.msi.figa <- tmp[["msi.g"]]
kfs.msi.figb <- tmp[["tyms.msi.inferred.g"]]


neoantigen.signature <- "neoantigens"
neoantigen.label <- "Log10 ( # Neoantigens + 1 )"

kras.cms.res <- list()
kras.cms.interaction.res <- list()
kras.cms.within.cms.res <- list()
kras.res <- list()
cms.res <- list()
site.res <- list()
msi.res <- list()

for(dataset in c("tcga", "kfs")) {
    ds <- tcga_expr
    if(dataset == "kfs") {
        ds <- kfsyscc_expr
    }

    es <- gsva.es.lst[[dataset]]
    clin.es <- clin.es.lst[[dataset]]

    to.plot <- c(circ.signature, immune.populations, hallmark.gene.sets, other.genes, neoantigen.signature, mcp.populations)
    ylabels <- c(circ.signature.label, immune.population.labels, hallmark.ylabels, other.labels, neoantigen.label, mcp.population.labels)

    to.plot <- circ.signature
    ylabels <- circ.signature.label
    
    res <- doKrasCMSInteractionAnalysis(es, clin.es, to.plot=to.plot, analysis.name = dataset, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.interaction.res[[dataset]] <- res$tests
    kras.cms.interaction.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms-interaction.tsv"), kras.cms.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    res <- doKrasCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.res[[dataset]] <- res$tests
    kras.cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms.tsv"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doKrasCMSWithinCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.within.cms.res[[dataset]] <- res$tests
    kras.cms.within.cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms-within-cms.tsv"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.res[[dataset]] <- res$tests
    kras.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras.tsv"), kras.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    cms.res[[dataset]] <- res$tests
    cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-cms.tsv"), cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doSiteAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    site.res[[dataset]] <- res$tests
    site.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-site.tsv"), site.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    if(dataset == "kfs") {
        clin.es <- clin.inferred.msi.lst[["kfs"]]
        clin.es$msi <- clin.es$msi.inferred
    }
    res <- doMSIAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    msi.res[[dataset]] <- res$tests
    msi.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-msi.tsv"), msi.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
}


num.bootstraps <- 10000
num.bootstraps <- 500
num.bootstraps <- 100000
to.plot <- immune.gene.sets
to.plot <- c("CIRC")
to.plot <- c(hallmark.gene.sets, immune.gene.sets)

run.kras.effect.size <- function(es1, clin.es1, es1.name, es2, clin.es2, es2.name, to.plot, num.bootstraps = 10, seed = 1234) {
    res1 <- doKrasEffectSizeAnalysis(es1, clin.es1, to.plot = to.plot, num.bootstraps = num.bootstraps, seed = seed)

    res2 <- doKrasEffectSizeAnalysis(es2, clin.es2, to.plot = to.plot, num.bootstraps = num.bootstraps, seed = seed)

    both.res <- merge(res1, res2, by = "gene.set", suffixes = c(paste0(".", es1.name), paste0(".", es2.name)))
    both.res
}

kras.results <- list()

plot.lst <- list("hallmark" = hallmark.gene.sets, "immune" = immune.gene.sets, "mcp" = rownames(tcga.mcp), "circ" = circ.genes)
plot.seeds <- c(1234,4321)

ns <- names(plot.lst)
ns <- "circ"
for(name in ns) {
    es1 <- gsva.es.lst[["tcga"]]
    clin1 <- clin.es.lst[["tcga"]]
    es.name1 <- "tcga"

    es2 <- gsva.es.lst[["kfs"]]
    clin2 <- clin.es.lst[["kfs"]]
    es.name2 <- "kfs"

    if(name == "mcp") {
        es1 <- tcga.mcp
        es2 <- kfs.mcp
    }
    
    kras.results[[name]] <- list()

    for(i in 1:length(plot.seeds)) {
        cat(paste0("Running KRAS on ", name, " with seed: ", plot.seeds[i], "\n"))
        res <- run.kras.effect.size(es1, clin1, es.name1, es2, clin2, es.name2, to.plot = plot.lst[[name]], num.bootstraps = num.bootstraps, seed = plot.seeds[i])
        kras.results[[name]][[i]] <- res
    }
}

## End analysis

save.image(".Rdata")

##
##    2a. make plot function -- like krascms -- at least align
##    2b. hallmark
##    2c. bindea
##    2e. mcp counter


stop("stop")

run.kras.cms.effect.size <- function(es1, clin.es1, es1.name, es2, clin.es2, es2.name, to.plot, num.bootstraps = 10, seed = 1234) {
    res1 <- doKrasCMSEffectSizeAnalysis(es1, clin.es1, to.plot = to.plot, num.bootstraps = num.bootstraps, seed = seed)

    res2 <- doKrasCMSEffectSizeAnalysis(es2, clin.es2, to.plot = to.plot, num.bootstraps = num.bootstraps, seed = seed)

    both.res <- merge(res1, res2, by = "gene.set", suffixes = c(paste0(".", es1.name), paste0(".", es2.name)))
    both.res
}

num.bootstraps <- 500
num.bootstraps <- 100000

kras.cms.results <- list()
for(name in names(plot.lst)) {
    es1 <- gsva.es.lst[["tcga"]]
    clin1 <- clin.es.lst[["tcga"]]
    es.name1 <- "tcga"

    es2 <- gsva.es.lst[["kfs"]]
    clin2 <- clin.es.lst[["kfs"]]
    es.name2 <- "kfs"

    if(name == "mcp") {
        es1 <- tcga.mcp
        es2 <- kfs.mcp
    }
    
    kras.cms.results[[name]] <- list()

    for(i in 1:length(plot.seeds)) {
        cat(paste0("Running KRAS on ", name, " with seed: ", plot.seeds[i], "\n"))
        res <- run.kras.cms.effect.size(es1, clin1, es.name1, es2, clin2, es.name2, to.plot = plot.lst[[name]], num.bootstraps = num.bootstraps, seed = plot.seeds[i])
        kras.cms.results[[name]][[i]] <- res
    }
}

## Compare the two bootstraps to see whether results are stable
## Report the maximum relative error--it is 1.2% for KRAS CMS mcp
for(name in names(plot.lst)) {
    cat(paste0("KRAS CMS ", name, "\n"))
    kc1 <- kras.cms.results[[name]][[1]]
    kc2 <- kras.cms.results[[name]][[2]]
    kc1 <- kc1[, !(colnames(kc1) %in% "gene.set")]
    kc2 <- kc2[, !(colnames(kc2) %in% "gene.set")]
    kc1 <- as.data.frame(sapply(kc1, as.numeric))
    kc2 <- as.data.frame(sapply(kc2, as.numeric))
    ## print(max(abs(kc1-kc2)))
    out <- outer(1:nrow(kc1), 1:ncol(kc1), FUN = function(i, j) abs(kc1[i,j] - kc2[i,j])/abs(min(kc1[i,j], kc2[i,j])))
    print(max(out))
    
}

for(name in names(plot.lst)) {
    cat(paste0("KRAS ", name, "\n"))
    kc1 <- kras.results[[name]][[1]]
    kc2 <- kras.results[[name]][[2]]
    kc1 <- kc1[, !(colnames(kc1) %in% "gene.set")]
    kc2 <- kc2[, !(colnames(kc2) %in% "gene.set")]
    kc1 <- as.data.frame(sapply(kc1, as.numeric))
    kc2 <- as.data.frame(sapply(kc2, as.numeric))
    ## print(max(abs(kc1-kc2)))
    out <- outer(1:nrow(kc1), 1:ncol(kc1), FUN = function(i, j) abs(kc1[i,j] - kc2[i,j])/abs(min(kc1[i,j], kc2[i,j])))
    print(max(out))
}


to.plot <- circ.signature
ylabels <- circ.signature.label

to.plot <- c(circ.signature, "IL6", "IRF1")
ylabels <- c(circ.signature.label, "IL6 Expression", "IRF1 Expression")

## TCGA CIRC vs KRAS + CMS
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
circ.kras.cms.res[[dataset]] <- res$tests

## KFS CIRC vs KRAS + CMS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
circ.kras.cms.res[[dataset]] <- res$tests

## TCGA CIRC vs KRAS
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
circ.kras.res[[dataset]] <- res$tests

## KFS CIRC vs KRAS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
circ.kras.res[[dataset]] <- res$tests

## Volcano plot of Bindea vs KRAS

## TCGA Bindea vs KRAS
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
to.plot <- immune.gene.sets
res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
bindea.kras.res[[dataset]] <- res$tests
tbl <- kras.results[["immune"]][[1]]
tmp <- bindea.kras.res[[dataset]]
tbl <- merge(tbl, tmp, by = "gene.set")

## KFS Binda vs KRAS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
bindea.kras.res[[dataset]] <- res$tests
tmp <- bindea.kras.res[[dataset]]
tbl <- merge(tbl, tmp, by = "gene.set", suffixes = c(".tcga", ".kfs"))

tbl.tcga <- tbl[,c(colnames(tbl)[grepl(colnames(tbl), pattern=".tcga")], "gene.set")]
colnames(tbl.tcga) <- gsub(colnames(tbl.tcga), pattern=".tcga", replacement="")
tbl.tcga$Dataset <- "TCGA"
tbl.kfs <- tbl[,c(colnames(tbl)[grepl(colnames(tbl), pattern=".kfs")], "gene.set")]
colnames(tbl.kfs) <- gsub(colnames(tbl.kfs), pattern=".kfs", replacement="")
tbl.kfs$Dataset <- "KFSYSCC"
tbl.kfs <- tbl.kfs[, colnames(tbl.tcga)]
tbl <- rbind(tbl.tcga, tbl.kfs)

tbl$pval <- as.numeric(tbl$pval)
tbl$mle <- as.numeric(tbl$mle)
tbl$low <- as.numeric(tbl$low)
tbl$high <- as.numeric(tbl$high)
tbl$Dataset <- factor(tbl$Dataset, levels = c("TCGA", "KFSYSCC"))
g <- ggplot(data = tbl, aes(x = mle, y = - log10(pval), colour = Dataset))
g <- g + geom_point()
## g <- g + geom_errorbarh(aes(xmin=low, xmax=high), height = 0.05)
tbl.sub <- subset(tbl, gene.set %in% c("Neutrophils", "Cytotoxic cells", "Th1 cells"))
tbl.sub <- tbl
## Just jitter Th1 KFS and iDC TCGA
flag <- ( ( tbl.sub$Dataset == "KFSYSCC") & ( tbl.sub$gene.set == "Th1 cells" ) ) | ( ( tbl.sub$Dataset == "TCGA") & ( tbl.sub$gene.set == "iDC" ) )
flag <- ( ( tbl.sub$Dataset == "KFSYSCC") & ( tbl.sub$gene.set == "Th1 cells" ) )
tbl.no.jit <- tbl[!flag, ]
g <- g + geom_text(data = tbl.no.jit, aes(x = mle, y = - log10(pval), label = gene.set), hjust = -0.1, show.legend = FALSE)
tbl.jit <- tbl[flag, ]
## set.seed(1234)
## Manually move Th1 out of the way
g <- g + geom_text(data = tbl.jit, aes(x = mle, y = - log10(pval), label = gene.set), hjust = 1.1, show.legend = FALSE)
## g <- g + geom_text(data = tbl.jit, aes(x = mle, y = - log10(pval), label = gene.set), hjust = -0.1, position = position_jitter(width = 0.02, height = 0.05))
g <- g + xlab("KRAS WT - KRAS MT")
g <- g + ylab(bquote(-Log[10]~italic('p')-value))
## lim <- max(max(abs(tbl$high)), max(abs(tbl$low)))
lim <- max(abs(tbl$mle))
g <- g + xlim(c(-lim, lim))
print(g)

## End volcano plot of Bindea vs KRAS


## Volcano plot of Hallmark vs KRAS

## TCGA Bindea vs KRAS
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
to.plot <- hallmark.gene.sets
res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
hallmark.kras.res[[dataset]] <- res$tests
tbl <- kras.results[["hallmark"]][[1]]
tmp <- hallmark.kras.res[[dataset]]
tbl <- merge(tbl, tmp, by = "gene.set")

## KFS Binda vs KRAS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
hallmark.kras.res[[dataset]] <- res$tests
tmp <- hallmark.kras.res[[dataset]]
tbl <- merge(tbl, tmp, by = "gene.set", suffixes = c(".tcga", ".kfs"))

tbl.tcga <- tbl[,c(colnames(tbl)[grepl(colnames(tbl), pattern=".tcga")], "gene.set")]
colnames(tbl.tcga) <- gsub(colnames(tbl.tcga), pattern=".tcga", replacement="")
tbl.tcga$Dataset <- "TCGA"
tbl.kfs <- tbl[,c(colnames(tbl)[grepl(colnames(tbl), pattern=".kfs")], "gene.set")]
colnames(tbl.kfs) <- gsub(colnames(tbl.kfs), pattern=".kfs", replacement="")
tbl.kfs$Dataset <- "KFSYSCC"
tbl.kfs <- tbl.kfs[, colnames(tbl.tcga)]
tbl <- rbind(tbl.tcga, tbl.kfs)

tbl$pval <- as.numeric(tbl$pval)
tbl$mle <- as.numeric(tbl$mle)
tbl$low <- as.numeric(tbl$low)
tbl$high <- as.numeric(tbl$high)
tbl$Dataset <- factor(tbl$Dataset, levels = c("TCGA", "KFSYSCC"))
tbl$gene.set <- gsub(tbl$gene.set, pattern="HALLMARK_", replacement="")
g.hall <- ggplot(data = tbl, aes(x = mle, y = - log10(pval), colour = Dataset))
g.hall <- g.hall + geom_point()

## Manually move two labels:
## TCGA UV_RESPONSE_UP
## KFS ESTRORGEN_RESPONSE_LATE
tbl.neg <- subset(tbl, pval < 0.1 & mle < 0)
tbl.pos <- subset(tbl, pval < 0.1 & mle > 0)
g.hall <- g.hall + geom_text(data = tbl.neg, aes(x = mle, y = - log10(pval), label = gene.set), hjust = 1.05, show.legend = FALSE, size = 2)
g.hall <- g.hall + geom_text(data = tbl.pos, aes(x = mle, y = - log10(pval), label = gene.set), hjust = -0.05, show.legend = FALSE, size = 2)
g.hall <- g.hall + xlab("KRAS WT - KRAS MT")
g.hall <- g.hall + ylab(bquote(-Log[10]~italic('p')-value))
lim <- max(abs(tbl$mle))
lim <- 0.4
g.hall <- g.hall + xlim(c(-lim, lim))
print(g.hall)


## End volcano plot of Hallmark vs KRAS


## Make univariate plots of kras mt - kras wt effect size for the following data sets:
## hallmark, bindea immune, and mcp (Figures 1, 2, and Supp Figure 5)

add.figure.labels <- TRUE

## Plot Supp Figure 1: CIRC vs KRAS
to.plot <- circ.signature
ylabels <- circ.signature.label

dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.tbl <- circ.kras.res[[dataset]]
p.tcga <- plotKrasAnalysis(es, clin.es, kras.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.tbl <- circ.kras.res[[dataset]]
p.kfs <- plotKrasAnalysis(es, clin.es, kras.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

pg <- plot_grid(p.tcga$CIRC, p.kfs$CIRC, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 1", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig1.tiff", width=14)

## Plot the pairwise overlap of the immune signatures and CIRC, restricting to expressed genes.
gene.set.names <- c("B cells", "T cells", "Th1 cells", "Th2 cells", "Cytotoxic cells", "iDC", "Neutrophils", "CIRC")
tcga.immune.overlap <- plotOverlapOfExpressedGeneSets(tcga_expr, bindeaGsets, gene.set.names = gene.set.names, main = "TCGA")
## TCGA: overlap between B cells and CIRC = HLA-DQA1
## TCGA: overlap between Th1 cells and CIRC = IFNG,CTLA4
## TCGA: overlap between Cytotoxic cells and CIRC = GNLY

kfs.immune.overlap <- plotOverlapOfExpressedGeneSets(kfsyscc_expr, bindeaGsets, gene.set.names = gene.set.names, main = "KFSYSCC")
## KFSYSCC: overlap between B cells and CIRC = HLA-DQA1
## KFSYSCC: overlap between Th1 cells and CIRC = IFNG,CTLA4
## KFSYSCC: overlap between Cytotoxic cells and CIRC = GNLY

pg <- plot_grid(tcga.immune.overlap, kfs.immune.overlap, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 2", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig2.tiff", width=14)


## Plot the pairwise correlation of the immune signatures and CIRC (Supp Figure 3)

tiff("supp-fig3a.tiff")
pairs(t(gsva.es.lst[["tcga"]][gene.set.names,]), upper.panel = panel.cor, panel = function(x, y, ...) { smoothScatter(x, y, ..., nrpoints = 0, add = TRUE); abline(lm(y~x), lwd=3); }, main="TCGA")
mtext("A", line=2.5, at = 0, cex = 1.5)
d <- dev.off()

tiff("supp-fig3b.tiff")
pairs(t(gsva.es.lst[["kfs"]][gene.set.names,]), upper.panel = panel.cor, panel = function(x, y, ...) { smoothScatter(x, y, ..., nrpoints = 0, add = TRUE); abline(lm(y~x), lwd=3); }, main="KFSYSCC")
mtext("B", line=2.5, at = 0, cex = 1.5)
if(add.figure.labels) {
    mtext("Supp Fig 3", line=2.5, at = 1, cex = 1.5)
}
d <- dev.off()

if(FALSE) {
    ## ggpairs does not create a gtable or a ggplot that we can use with grid.arrange.
    ## So, I think we're stuck with base graphics above
    df <- as.data.frame(t(gsva.es.lst[["kfs"]][gene.set.names,]))
    cols <- gsub(colnames(df), pattern=" ", replacement="\n")
    colnames(df) <- make.names(colnames(df))
    pm <- ggpairs(df, lower = list(continuous = my_fn2), diag = list(continuous = ggally_mysmooth), upper = list(continuous = ggally_mycor), axisLabels="show", columnLabels=cols)
    ## pm <- ggpairs(df, lower = list(continuous = my_fn2), diag = list(continuous = wrap('diagAxis', labelSize = 3)), upper = list(continuous = ggally_mycor), axisLabels="show", columnLabels=cols)
    pm$showXAxisPlotLabels <- FALSE
    pm$showYAxisPlotLabels <- FALSE
    pm$xAxisLabels <- NULL
    pm$yAxisLabels <- NULL
    pm
}


## Plot the dendrograms of the Bindea immune pathways
## Nah -- this confuses the issue; the scatter plots are more informative.

pdf("immune-dendros.pdf")
methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
for(method in methods) {
    par(mfrow=c(1,2))    
    tcga.dendro.g <- my_ggdendro(t(gsva.es.lst[["tcga"]][c(circ.signature, immune.populations),]), method = method, main = paste0("TCGA: ", method))
    kfs.dendro.g <- my_ggdendro(t(gsva.es.lst[["kfs"]][c(circ.signature, immune.populations),]), method = method, main = paste0("KFSYSCC: ", method))
    grid.arrange(tcga.dendro.g, kfs.dendro.g)
##    plot_grid(tcga.dendro.g, kfs.dendro.g, labels = c("A", "B"), nrow=1)
}
d <- dev.off()

tcga.dendro.g <- my_ggdendro(t(gsva.es.lst[["tcga"]][c(circ.signature, immune.populations),]), method = "average", main = "TCGA")
kfs.dendro.g <- my_ggdendro(t(gsva.es.lst[["kfs"]][c(circ.signature, immune.populations),]), method = "average", main = "KFSYSCC")
if(FALSE) {
pg <- plot_grid(tcga.dendro.g, kfs.dendro.g, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 4", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig4.tiff", width=14)
}


ns <- c("immune", "hallmark", "circ")
files <- c("fig1", "fig2", "supp-fig4")
labels <- c("Figure 1", "Figure 2", "Supp Figure 5")
if(!add.figure.labels) {
    labels <- c(NULL, NULL, NULL)
}

for(i in 1:length(ns)) {
    name <- ns[i]
    tbl <- kras.results[[name]][[1]]
    tbl$gene.set <- gsub(tbl$gene.set, pattern="HALLMARK_", replacement="")
    tbl1.mle.col <- "mle.tcga"
    tbl1.lb.col <- "low.tcga"
    tbl1.ub.col <- "high.tcga"
    tbl2.mle.col <- "mle.kfs"
    tbl2.lb.col <- "low.kfs"
    tbl2.ub.col <- "high.kfs"
    main1 <- "TCGA"
    main2 <- "KFSYSCC"
    ## Add extra space for hallmark
    if(name == "hallmark") { main2 <- "\nKFSYSCC" }
    text.size <- NULL
    if(name == "hallmark") { text.size <- 5 }
    switch.sign <- FALSE
    ## We re-ran circ with the fixed code, so no need to invert
    if(name != "circ") {
        tbl[, tbl1.mle.col] <- - as.numeric(tbl[, tbl1.mle.col])
        tbl[, tbl1.lb.col] <- - as.numeric(tbl[, tbl1.lb.col])
        tbl[, tbl1.ub.col] <- - as.numeric(tbl[, tbl1.ub.col])                
        tbl[, tbl2.mle.col] <- - as.numeric(tbl[, tbl2.mle.col])
        tbl[, tbl2.lb.col] <- - as.numeric(tbl[, tbl2.lb.col])
        tbl[, tbl2.ub.col] <- - as.numeric(tbl[, tbl2.ub.col])
        tbl <- tbl[order(tbl$high.tcga,decreasing=TRUE),]
    } else {
        tbl <- tbl[order(tbl$low.tcga,decreasing=TRUE),]
    }    
    tbl$gene.set <- factor(tbl$gene.set, levels = tbl$gene.set)
    l <- plot.stacked.effect.size.figures(tbl, tbl, tbl1.mle.col, tbl1.lb.col, tbl1.up.col, tbl2.mle.col, tbl2.lb.col, tbl2.ub.col, main1, main2, text.size = text.size, switch.sign = switch.sign)

    pg <- plot_grid(l[[1]], l[[2]], labels = c("A", "B"), nrow=2)
    label <- labels[i]
    if(!is.null(label)) {
        pg <- pg + draw_label(label, x=1, y=1, vjust=1, hjust=1)
    }
    pg
    filename <- paste0(name, "-kras-effect.tiff")
    filename <- paste0(files[i], ".tiff")
    ggsave(filename)
}

## Plot univariate analysis for TCGA (Supp Figure 5) and for KFS (Supp Figure 8)

## TCGA
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
cms.tbl <- cms.res[[dataset]]
site.tbl <- site.res[[dataset]]
msi.tbl <- msi.res[[dataset]]

to.plot <- circ.signature
ylabels <- circ.signature.label

main <- ifelse(dataset == "tcga", "TCGA", "KFSYSCC")
main <- NULL

p.cms.tcga <- plotCMSAnalysis(es, clin.es, cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

p.site.tcga <- plotSiteAnalysis(es, clin.es, site.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

p.msi.tcga <- plotMSIAnalysis(es, clin.es, msi.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

## Plot CIRC vs neoantigens for TCGA
p.neo.tcga <- plot.circ.vs.neoantigens(es, clin, main = main)

grid.newpage()
pg <- plot_grid(p.site.tcga$CIRC, p.cms.tcga$CIRC, p.neo.tcga, p.msi.tcga$CIRC, labels = c("A", "B", "C", "D"), nrow=2, ncol=2)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 5", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig5.tiff", width=14, height=14)


## KFS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]

es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
cms.tbl <- cms.res[[dataset]]
site.tbl <- site.res[[dataset]]
msi.tbl <- msi.res[[dataset]]

to.plot <- circ.signature
ylabels <- circ.signature.label

main <- ifelse(dataset == "tcga", "TCGA", "KFSYSCC")
main <- NULL

p.cms.kfs <- plotCMSAnalysis(es, clin.es, cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

p.site.kfs <- plotSiteAnalysis(es, clin.es, site.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

clin.es <- clin.inferred.msi.lst[["kfs"]]
clin.es$msi <- clin.es$msi.inferred

p.msi.kfs <- plotMSIAnalysis(es, clin.es, msi.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

grid.newpage()
pg <- plot_grid(p.site.kfs$CIRC, p.cms.kfs$CIRC, p.msi.kfs$CIRC, labels = c("A", "B", "C"), nrow=2, ncol=2)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 8", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig8.tiff", width=14, height=14)

## Plot the imputed MSI status for TCGA (Supp Figure 6)
grid.newpage()
pg <- plot_grid(tcga.msi.figa, tcga.msi.figb, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 6", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig6.tiff", width=14)

## Plot the imputed MSI status for KFS (Supp Figure 7)
grid.newpage()
pg <- plot_grid(kfs.msi.figa, kfs.msi.figb, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 7", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig7.tiff", width=14)

## Plot univariate analysis for TCGA (Supp Figure 8)

## Plot Figure 3 CIRC vs KRAS + CMS and Figure 4 CIRC forest
to.plot <- circ.signature
ylabels <- circ.signature.label

to.plot <- circ.signature
ylabels <- circ.signature.label

## TCGA
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.cms.tbl <- circ.kras.cms.res[[dataset]]
p.tcga <- plotKrasCMSAnalysis(es, clin.es, kras.cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
p.forest.tcga <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

## KFS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.cms.tbl <- circ.kras.cms.res[[dataset]]
p.kfs <- plotKrasCMSAnalysis(es, clin.es, kras.cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
p.kfs.tcga <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
p.forest.kfs <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

gene <- "CIRC"
pg <- plot_grid(p.tcga[[gene]], p.kfs[[gene]], labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Figure 3", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("fig3.tiff", width=14)

pg <- plot_grid(p.forest.tcga[[gene]]$plot, p.forest.kfs[[gene]]$plot, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Figure 4", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("fig4.tiff", width=14)


## Make multivariate cms-kras plots of hallmark, bindea immune, and mcp (Figures 5 and 6)
ns <- c("immune", "hallmark")
files <- c("fig5", "fig6")
labels <- c("Figure 5", "Figure 6")
if(!add.figure.labels) {
    labels <- c(NULL, NULL)
}
for(i in 1:length(ns)) {
    name <- ns[i]
    tbl <- kras.cms.results[[name]][[1]]

    tbl$gene.set <- gsub(tbl$gene.set, pattern="HALLMARK_", replacement="")
    
    tbl.kfs <- tbl[, grepl(colnames(tbl), pattern="kfs") | (colnames(tbl) == "gene.set")]
    colnames(tbl.kfs) <- gsub(colnames(tbl.kfs), pattern=".kfs", replacement="")
    tbl.kfs <- tbl.kfs  %>% gather(key, value, -gene.set) %>% extract(key, c("cms", "kras", "interval"), "(CMS.)(.+)\\.(.+)") %>% spread(interval, value) 
    
    tbl.tcga <- tbl[, grepl(colnames(tbl), pattern="tcga") | (colnames(tbl) == "gene.set")]
    colnames(tbl.tcga) <- gsub(colnames(tbl.tcga), pattern=".tcga", replacement="")
    tbl.tcga <- tbl.tcga  %>% gather(key, value, -gene.set) %>% extract(key, c("cms", "kras", "interval"), "(CMS.)(.+)\\.(.+)") %>% spread(interval, value) 
    
    stbl <- ddply(tbl.tcga[, c("gene.set", "low")], .variables = "gene.set", .fun = function(df) min(as.numeric((df$low))))
    colnames(stbl) <- c("gene.set", "min.effect")
    stbl <- stbl[order(as.numeric(stbl$min.effect), decreasing=TRUE),]
    
    gene.set.order <- unique(stbl$gene.set)
    
    tbl.kfs$gene.set <- factor(tbl.kfs$gene.set, levels = gene.set.order)
    tbl.tcga$gene.set <- factor(tbl.tcga$gene.set, levels = gene.set.order)
    tbl.kfs$cms <- factor(tbl.kfs$cms, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
    tbl.tcga$cms <- factor(tbl.tcga$cms, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
    tbl.kfs$kras <- factor(tbl.kfs$kras, levels = c("MT", "WT"))
    tbl.tcga$kras <- factor(tbl.tcga$kras, levels = c("MT", "WT"))
    
    tbl1.mle.col <- "mle"
    tbl1.lb.col <- "low"
    tbl1.ub.col <- "high"
    tbl2.mle.col <- "mle"
    tbl2.lb.col <- "low"
    tbl2.ub.col <- "high"
    main1 <- "TCGA"
    main2 <- "KFSYSCC"

    text.size <- NULL
    if(name == "hallmark") { text.size <- 5 }

    l <- plot.stacked.faceted.effect.size.figures(tbl.tcga, tbl.kfs, tbl1.mle.col, tbl1.lb.col, tbl1.up.col, tbl2.mle.col, tbl2.lb.col, tbl2.ub.col, main1, main2, text.size = text.size)
##    do.call("grid.arrange", l)
    pg <- plot_grid(l[[1]], l[[3]], l[[2]], labels = c("A", "", "B"), nrow=2, rel_widths = c(10,1.5), align="none")
    label <- labels[i]
    if(!is.null(label)) {
        pg <- pg + draw_label(label, x=1, y=1, vjust=1, hjust=1)
    }
    pg
    filename <- paste0(name, "-kras-effect.tiff")
    filename <- paste0(files[i], ".tiff")
    ggsave(filename)
    ## pg <- plot_grid(l[[1]], NULL, l[[2]], labels = c("A", "", "B"), nrow=2, rel_widths = c(10,1.5), align="none")
    ## pg + draw_grob(l[[3]], x = 0, y = 0.5)
    ## ggsave(paste0(name, "-kras-cms-effect.tiff"))
    if(FALSE) {
        pg <- plot_grid(l[[1]], l[[2]], labels = c("A", "B"), nrow=2, rel_widths = c(10,1.5), align="none")
        tiff(paste0(name, "-kras-cms-effect.tiff"))
        grid.arrange(pg, l[[3]], nrow=1, widths = c(10,1.5))
        d <- dev.off()
    }
}

hall <- kras.cms.results[["hallmark"]][[1]]
hall <- hall[order(hall$high.tcga),]
flag <- unlist(apply(hall[, grepl(colnames(hall), pattern="tcga")], 1, function(row) all(row > 0)))
hall[flag, "gene.set"]
## [1] "HALLMARK_ALLOGRAFT_REJECTION"       "HALLMARK_COMPLEMENT"               
## [3] "HALLMARK_IL2_STAT5_SIGNALING"       "HALLMARK_IL6_JAK_STAT3_SIGNALING"  
## [5] "HALLMARK_INFLAMMATORY_RESPONSE"     "HALLMARK_INTERFERON_ALPHA_RESPONSE"
## [7] "HALLMARK_INTERFERON_GAMMA_RESPONSE"
sig <- hall[flag,]
tbl <- table(unlist(lapply(sig[,"gene.set"], function(x) bindeaGsets[[as.character(x)]])))
sort(tbl[tbl > 2])
## Genes enriched in hallmark pathways
##
##     B2M     BST2      C1S    CASP1    CASP3     CCL2     CD40     CD74 
##       3        3        3        3        3        3        3        3 
##    CD86  EIF2AK2      FAS     GBP4     GZMA    HIF1A    ICAM1   IFITM3 
##       3        3        3        3        3        3        3        3 
##  IFNGR1   IFNGR2     IL10   IL10RA   IL15RA   IL18R1     IL1B    IL2RA 
##       3        3        3        3        3        3        3        3 
##    IRF2     IRF4     IRF8     IRF9     JAK2     LAP3      LCK      LIF 
##       3        3        3        3        3        3        3        3 
##     LTB     LY6E      LYN      NMI   PIK3R5    PSMB9     RTP4     SELL 
##       3        3        3        3        3        3        3        3 
##   STAT1    STAT2     TAP1    TAPBP    TIMP1     TLR2 TNFRSF1B  TNFSF10 
##       3        3        3        3        3        3        3        3 
##    WARS     CCL5     CCL7   CXCL11    CXCL9     IL15    IL2RB      IL7 
##       3        4        4        4        4        4        4        4 
##    LCP2     PIM1   PLSCR1    RIPK2    SOCS1     CSF1   CXCL10      IL6 
##       4        4        4        4        4        5        5        5 
##    IRF1     IRF7     IL4R 
##       5        5        6 
## as.data.frame(sort(tbl[tbl > 2]))
stop("stop")


## Convert TIFFs to PDFs
## for x in `ls *.tiff`; do echo $x; sips -s format pdf $x --out `basename $x .tiff`.pdf ; done

## "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" --output crc-cms-kras-figures.pdf supp-fig1.pdf supp-fig2.pdf supp-fig3a.pdf supp-fig3b.pdf fig1.pdf fig2.pdf supp-fig4.pdf supp-fig5.pdf supp-fig6.pdf supp-fig7.pdf supp-fig8.pdf fig3.pdf fig4.pdf fig5.pdf fig6.pdf

hall <- kras.results[["hallmark"]][[1]]
hall <- hall[order(hall$high.tcga),]

## Hallmark with KRAS 95% CI < 0 (almost)
## hall[1:5, "gene.set"]
## [1] HALLMARK_ALLOGRAFT_REJECTION       HALLMARK_INTERFERON_GAMMA_RESPONSE
## [3] HALLMARK_COMPLEMENT                HALLMARK_UV_RESPONSE_UP           
## [5] HALLMARK_IL6_JAK_STAT3_SIGNALING

## tbl <- table(unlist(lapply(hall[(hall$high.tcga < 10^-2) & (hall$high.kfs < 10^-2),"gene.set"], function(x) bindeaGsets[[as.character(x)]])))
## sort(tbl[tbl > 2])
## CCL5  CCL7 CXCL9   FAS  GZMA  IL4R   IL7  IRF1  IRF7  JAK2  LCP2  PIM1 SOCS1 
##    3     3     3     3     3     3     3     3     3     3     3     3     3 
##STAT1   IL6 
##    3     4 

## Genes enriched in hallmark pathways
tbl <- table(unlist(lapply(hall[1:5,"gene.set"], function(x) bindeaGsets[[as.character(x)]])))
sort(tbl[tbl > 2])

## CASP3  CCL5  CCL7 CXCL9   FAS  GZMA ICAM1  IL4R   IL7  IRF7  JAK2  LCP2   LYN 
##     3     3     3     3     3     3     3     3     3     3     3     3     3 
##  PIM1 SOCS1 STAT1  TAP1  IRF1   IL6 
## 3     3     3     3     4     5

save.image(".Rdata")

stop("stop with cms kras")

## both.kras1 <- run.kras.effect.size(gsva.es.lst[["tcga"]], clin.es.lst[["tcga"]], "tcga", gsva.es.lst[["kfs"]], clin.es.lst[["kfs"]], "kfs", to.plot = to.plot, num.bootstraps = 500, seed = 1234)

dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
tcga.kras1 <- doKrasEffectSizeAnalysis(es, clin.es, to.plot=to.plot, num.bootstraps = num.bootstraps, seed = 1234)

tcga.kras2 <- doKrasEffectSizeAnalysis(es, clin.es, to.plot=to.plot, num.bootstraps = num.bootstraps, seed = 4321)

dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kfs.kras <- doKrasEffectSizeAnalysis(es, clin.es, to.plot=to.plot, num.bootstraps = num.bootstraps)

both.kras <- merge(tcga.kras, kfs.kras, by = "gene.set", suffixes = c(".tcga", ".kfs"))

if(FALSE) {
g <- ggplot(both.kras, aes(x = mle.tcga, y = mle.kfs))
g <- g + geom_point()
g <- g + geom_errorbar(aes(ymin=low.kfs, ymax=high.kfs), width=0.01)
g <- g + geom_errorbarh(aes(xmin=low.tcga, xmax=high.tcga), height=0.01) 
g <- g + geom_text(aes(label = gene.set), hjust = -0.1, vjust = 0)
g <- g + xlab("TCGA")
g <- g + ylab("KFSYSCC")

print(g)
}

stop("stop")

## todo to do
## 3. dokrascmseffect
##    3a. make plot that stacks
##        3a1.  bindea
##        3a3.  mcp counter
##    3b. make plot for hallmark
##        3b1.  scatter by hand on x axis
##        3b2.  order from CMS1 MT, WT ... CMS WT
##        3b3.  x label is "" or hallmark
##        3b4.  plot hallmark

if(FALSE) {
plotme <- c(hallmark.gene.sets, immune.gene.sets)
plotme <- immune.gene.sets
res <- doKrasCMSEffectSizeAnalysis(es, clin.es, to.plot=plotme, num.bootstraps = 10000)

tbl <- res$tests[!grepl(res$tests$gene.set, pattern="HALLMARK"),]  %>% gather(key, value, -gene.set) %>% extract(key, c("comparison", "direction"), "(.+.CMS2MT)\\.(.+)") %>% spread(direction, value)

tbl$mn <- as.numeric(tbl$mn)
tbl$hi <- as.numeric(tbl$hi)
tbl$low <- as.numeric(tbl$low)

## Order by maximum CMS2WT.vs.CMS2MT.pval in TCGA
stbl <- ddply(tbl[, c("gene.set", "hi")], .variables = "gene.set", .fun = function(df) max(as.numeric((df$hi))))
colnames(stbl) <- c("gene.set", "min.effect")
stbl <- stbl[order(as.numeric(stbl$min.effect), decreasing=FALSE),]

sorted.levels <- unique(stbl$gene.set)
## tbl$gene.set <- factor(tbl$gene.set, levels = sorted.levels)


tbl$gene.set <- factor(tbl$gene.set, levels = sorted.levels)

ggplot(tbl, aes(x=gene.set, y=mn, colour=comparison)) + 
    geom_errorbar(aes(ymin=low, ymax=hi), width=1) +
    geom_line() +
    geom_point() + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, hjust = 1), axis.ticks.x=element_blank())
head(sorted.levels)

ggplot(tbl, aes(x=comparison, y=mn, colour=comparison)) + 
    geom_errorbar(aes(ymin=low, ymax=hi), width=1) +
    geom_line() +
    geom_point() + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, hjust = 1), axis.ticks.x=element_blank()) + facet_wrap(~gene.set, nrow=1)
head(sorted.levels) 
}



res <- doKrasCMSEffectSizeAnalysis(es, clin.es, to.plot=hallmark.gene.sets, num.bootstraps = 10000)

res <- doKrasCMSAnalysis(es, clin.es, paste0(dataset, "-hallmark-kras-cms-fdr"), to.plot=hallmark.gene.sets, num.permutations = 10000)
hallmark.kras.cms.res[[dataset]] <- res$tests
write.table(file = paste0(dataset, "-hallmark-kras-cms-fdr.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## Call FDR q < 25% significant
sig.flag <- unlist(apply(res$tests[, grepl(pattern="apval", colnames(res$tests))], 1, function(row) all(row < 0.25)))
tcga.sig.hallmark.gene.sets <- res$tests[sig.flag,"gene"]
cat("Hallmark gene sets that are significant in TCGA:\n")
print(tcga.sig.hallmark.gene.sets)

.simpleCap <- function(x) {
         s <- strsplit(x, " ")[[1]]
         paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
}

to.plot <- tcga.sig.hallmark.gene.sets
## ylabels <- unlist(lapply(to.plot, function(str) paste0(.simpleCap(tolower(gsub(gsub(str, pattern="_", replacement=" "), pattern="HALLfooMARK ", replacement=""))), "\nEnrichment Score")))
ylabels <- unlist(lapply(to.plot, function(str) paste0(str, "\nEnrichment Score")))
plotKrasCMSAnalysis(es, clin.es, hallmark.kras.cms.res[[dataset]], to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-hallmark-kras-cms-fdr"), plot.adjusted.pvalues = TRUE, plot.pvals.as.stars = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

for(gs in tcga.sig.hallmark.gene.sets) {
    cat(paste0(gs, ": ", paste(sort(bindeaGsets[[gs]]), collapse=", "), "\n\n"))
}
## Now test in KFS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasCMSAnalysis(es, clin.es, paste0(dataset, "-hallmark-kras-cms-fdr"), to.plot=tcga.sig.hallmark.gene.sets, num.permutations = 10000)
hallmark.kras.cms.res[[dataset]] <- res$tests
write.table(file = paste0(dataset, "-hallmark-kras-cms-fdr.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.flag <- unlist(apply(res$tests[, grepl(pattern="apval", colnames(res$tests))], 1, function(row) all(row < 0.25)))
kfs.sig.hallmark.gene.sets <- res$tests[sig.flag,"gene"]
cat("Hallmark gene sets that are significant in KFS:\n")
print(kfs.sig.hallmark.gene.sets)

to.plot <- kfs.sig.hallmark.gene.sets
## ylabels <- unlist(lapply(to.plot, function(str) paste0(.simpleCap(tolower(gsub(gsub(str, pattern="_", replacement=" "), pattern="HALLfooMARK ", replacement=""))), "\nEnrichment Score")))
ylabels <- unlist(lapply(to.plot, function(str) paste0(str, "\nEnrichment Score")))
plotKrasCMSAnalysis(es, clin.es, hallmark.kras.cms.res[[dataset]], to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-hallmark-kras-cms-fdr"), plot.adjusted.pvalues = TRUE, plot.pvals.as.stars = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

## Do the same for immune signatures--discover in TCGA and test in KFS.
## But these sets are independent, so we can calculate FDR using BH.
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasCMSAnalysis(es, clin.es, paste0(dataset, "-immune-kras-cms-fdr"), to.plot=immune.gene.sets, num.permutations = 0, calc.fdr.with.bh = TRUE)
immune.kras.cms.res[[dataset]] <- res$tests
write.table(file = paste0(dataset, "-immune-kras-cms-fdr.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.flag <- unlist(apply(res$tests[, grepl(pattern="apval", colnames(res$tests))], 1, function(row) all(row < 0.25)))
tcga.sig.immune.gene.sets <- res$tests[sig.flag,"gene"]
cat("Immune gene sets that are significant in TCGA:\n")
print(tcga.sig.immune.gene.sets)

to.plot <- tcga.sig.immune.gene.sets
## ylabels <- unlist(lapply(to.plot, function(str) paste0(.simpleCap(tolower(gsub(gsub(str, pattern="_", replacement=" "), pattern="HALLfooMARK ", replacement=""))), "\nEnrichment Score")))
ylabels <- unlist(lapply(to.plot, function(str) paste0(str, "\nEnrichment Score")))
plotKrasCMSAnalysis(es, clin.es, immune.kras.cms.res[[dataset]], to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-immune-kras-cms-fdr"), plot.adjusted.pvalues = TRUE, plot.pvals.as.stars = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))


dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
res <- doKrasCMSAnalysis(es, clin.es, paste0(dataset, "-immune-kras-cms-fdr"), to.plot=tcga.sig.immune.gene.sets, num.permutations = 0, calc.fdr.with.bh = TRUE)
immune.kras.cms.res[[dataset]] <- res$tests
write.table(file = paste0(dataset, "-immune-kras-cms-fdr.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.flag <- unlist(apply(res$tests[, grepl(pattern="apval", colnames(res$tests))], 1, function(row) all(row < 0.25)))
kfs.sig.immune.gene.sets <- res$tests[sig.flag,"gene"]
cat("Immune gene sets that are significant in KFS:\n")
print(kfs.sig.immune.gene.sets)

to.plot <- tcga.sig.immune.gene.sets
## ylabels <- unlist(lapply(to.plot, function(str) paste0(.simpleCap(tolower(gsub(gsub(str, pattern="_", replacement=" "), pattern="HALLfooMARK ", replacement=""))), "\nEnrichment Score")))
ylabels <- unlist(lapply(to.plot, function(str) paste0(str, "\nEnrichment Score")))
plotKrasCMSAnalysis(es, clin.es, immune.kras.cms.res[[dataset]], to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-immune-kras-cms-fdr"), plot.adjusted.pvalues = TRUE, plot.pvals.as.stars = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

for(dataset in c("tcga", "kfs")) {
    es <- NULL
    ds <- tcga_expr
    if(dataset == "kfs") {
        ds <- kfsyscc_expr
    }

    if(FALSE) {
    es <- gsva.es.lst[[dataset]]
    res <- doKrasCMSAnalysis(es, clin.es, paste0(dataset, "-hallmark-kras-cms-fdr"), to.plot=hallmark.gene.sets, num.permutations = 1000)
    hallmark.kras.cms.res[[dataset]] <- res$tests
    write.table(file = paste0(dataset, "-hallmark-kras-cms-fdr.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doKrasCMSAnalysis(es, clin.es, paste0(dataset, "-immune-kras-cms-fdr"), to.plot=immune.gene.sets, num.permutations = 1000)
    immune.kras.cms.res[[dataset]] <- res$tests
    write.table(file = paste0(dataset, "-immune-kras-cms-fdr.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    next
    }
    
    main = toupper(dataset)
    if(main == "KFS") { main <- "KFSYSCC" }
    res <- doAnalysis(ds, clin, paste0(dataset, "-chemo-infg"), to.plot=to.plot, ylabels=ylabels, main = main)
    
    g1 <- plot.correlation(es["CIRC", ], es["IFNG",], "circ", "IFNG")
    g2 <- plot.correlation(es["CIRC", ], es["HALLMARK_INTERFERON_GAMMA_RESPONSE",], "circ", "hallmark IFNG")
    g3 <- plot.correlation(es["CIRC", ], es["REACTOME_INTERFERON_GAMMA_SIGNALING",], "circ", "reactome IFNG")
    pdf(paste0(dataset, "-", "circ", "-vs-ifng.pdf"))
    l <- list(g1, g2, g3)
    do.call("grid.arrange", l)
    d <- dev.off()

    
    ## Correlation of CIRC with iDC and chemokines
    genes <- c("CIRC", "iDC", "BATF3", "THBD", "CCL4", "ATF3")
    for(gene in genes) {
        l <- llply(genes[!(genes == gene)], .parallel = TRUE,
               .fun = function(gene2) {
                   g <- plot.correlation(es[gene, ], es[gene2, ], x.name = gene, y.name = gene2)
                   g
               })
        pdf(paste0(dataset, "-", gene, "-vs-chemokines.pdf"))
        do.call("grid.arrange", l)
        d <- dev.off()
    }
}

## Plot figure
## Plot cytotoxic
## Plot hallmark infg

save.image(".Rdata")


mcp.kras.cms.res <- list()

dataset <- "tcga"
es <- gsva.es.lst[[dataset]]

tcga.mcp <- MCPcounter.estimate(tcga_expr, featuresType="HUGO_symbols")
to.plot <- rownames(tcga.mcp)
ylabels <- unlist(lapply(to.plot, function(x) paste0(x, " Population Percentage")))
ylabels <- unlist(lapply(to.plot, function(x) gsub(paste0(x, " Population Percentage"), pattern="MCP ", replacement="")))
rownames(tcga.mcp) <- unlist(lapply(rownames(tcga.mcp), function(x) paste0("MCP ", x)))
to.plot <- rownames(tcga.mcp)
## common.cols <- intersect(colnames(tcga.mcp), colnames(tcga_expr))
common.cols <- intersect(colnames(tcga.mcp), colnames(es))
tcga.mcp <- tcga.mcp[, common.cols]
## tcga_expr <- tcga_expr[, common.cols]
## expr.and.mcp <- rbind(tcga.mcp, tcga_expr)
expr.and.mcp <- rbind(tcga.mcp, es[, common.cols])

hallmark.ylabels <- unlist(lapply(hallmark.gene.sets, function(str) paste0(str, "\nEnrichment Score")))

to.plot <- c(circ.signature, immune.populations, hallmark.gene.sets)
ylabels <- c(circ.signature.label, immune.population.labels, hallmark.ylabels)

## res <- doAnalysis(expr.and.mcp, clin, "tcga-mcp", to.plot=to.plot, ylabels=ylabels, main="TCGA", do.gene.set.heatmap = FALSE)
## cluster.annotated.matrix(tcga.mcp, clin, "tcga-mcp")

dataset <- "tcga"
res <- doKrasCMSAnalysis(expr.and.mcp, clin, paste0(dataset, "-mcp-kras-cms-fdr"), to.plot=to.plot, calc.fdr = FALSE)
mcp.kras.cms.res[[dataset]] <- res$tests
tbl <- res$tests
write.table(file = paste0(dataset, "-immune-hallmark-kras-cms.tsv"), tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.flag <- unlist(apply(res$tests[, grepl(pattern="apval", colnames(res$tests))], 1, function(row) all(row < 0.25)))
tcga.sig.mcp.gene.sets <- res$tests[sig.flag,"gene"]

to.plot <- tcga.sig.mcp.gene.sets
ylabels <- unlist(lapply(to.plot, function(x) paste0(x, " Population Percentage")))
ylabels <- unlist(lapply(to.plot, function(x) gsub(paste0(x, " Population Percentage"), pattern="MCP ", replacement="")))

plotKrasCMSAnalysis(expr.and.mcp, clin, res$tests, to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-mcp-kras-cms-fdr"), plot.adjusted.pvalues = TRUE, plot.pvals.as.stars = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))



dataset <- "kfs"
es <- gsva.es.lst[[dataset]]

## kfs.mcp <- MCPcounter.estimate(kfsyscc_expr, featuresType="HUGO_symbols")
kfs.mcp <- MCPcounter.estimate(kfsyscc_microarray_expr, featuresType="affy133P2_probesets")
to.plot <- rownames(kfs.mcp)
kfs.expr.mat <- kfsyscc_microarray_expr
ylabels <- unlist(lapply(to.plot, function(x) paste0(x, " Population Percentage")))
rownames(kfs.mcp) <- unlist(lapply(rownames(kfs.mcp), function(x) paste0("MCP ", x)))
to.plot <- rownames(kfs.mcp)
## Only plot/test those sig in TCGA
## to.plot <- tcga.sig.mcp.gene.sets
## common.cols <- intersect(colnames(kfs.mcp), colnames(kfs.expr.mat))
common.cols <- intersect(colnames(kfs.mcp), colnames(es))
kfs.mcp <- kfs.mcp[, common.cols]
kfs.expr.mat <- kfs.expr.mat[, common.cols]
## expr.and.mcp <- as.matrix(rbind(kfs.mcp, kfs.expr.mat))
expr.and.mcp <- as.matrix(rbind(kfs.mcp, es[, common.cols]))

to.plot <- c(to.plot, hallmark.gene.sets)
ylabels <- c(ylabels, hallmark.ylabels)

res <- doAnalysis(expr.and.mcp, clin, "kfs-mcp", to.plot=to.plot, ylabels=ylabels, main="KFSYSCC", do.gene.set.heatmap = FALSE, do.gsva = FALSE)
## cluster.annotated.matrix(kfs.mcp, clin, "kfs-mcp")

dataset <- "kfs"

to.plot <- c(circ.signature, immune.populations, hallmark.gene.sets)
ylabels <- c(circ.signature.label, immune.population.labels, hallmark.ylabels)

res <- doKrasCMSAnalysis(expr.and.mcp, clin, paste0(dataset, "-mcp-kras-cms-fdr"), to.plot=to.plot, num.permutations = 0, calc.fdr.with.bh = TRUE)
mcp.kras.cms.res[[dataset]] <- res$tests
tbl <- res$tests
write.table(file = paste0(dataset, "-immune-hallmark-kras-cms.tsv"), res$tests, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.flag <- unlist(apply(res$tests[, grepl(pattern="apval", colnames(res$tests))], 1, function(row) all(row < 0.25)))
kfs.sig.mcp.gene.sets <- res$tests[sig.flag,"gene"]
to.plot <- kfs.sig.mcp.gene.sets
ylabels <- unlist(lapply(to.plot, function(x) paste0(x, " Population Percentage")))

## None are significant--just plot those from tcg
to.plot <- tcga.sig.mcp.gene.sets
ylabels <- unlist(lapply(to.plot, function(x) gsub(paste0(x, " Population Percentage"), pattern="MCP ", replacement="")))

plotKrasCMSAnalysis(expr.and.mcp, clin, res$tests, to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-mcp-kras-cms-fdr"), plot.adjusted.pvalues = TRUE, plot.pvals.as.stars = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))


print(sessionInfo())

save.image(".Rdata")

## Find the largest independent set of gene sets
find.largest.independent.sets <- function(sets) {
    suppressPackageStartupMessages(library("igraph"))
    set.names <- names(sets)
    set.connectivity.tbl <- ldply(set.names,
                                  .fun = function(set.name1) {
                                      ldply(set.names, .fun = function(set.name2) {
                                          connected <- ifelse(length(intersect(sets[[set.name1]], sets[[set.name2]])) > 0, 1, 0)
                                          vec <- c(set.name1, set.name2, connected)
                                          names(vec) <- c("set1", "set2", "connected")
                                          vec
                                      })
                                  })
    set.connectivity.tbl <- subset(set.connectivity.tbl, connected == 1)

    ## Create an igraph, where nodes are sets and edges indicated that the
    ## two sets shared at least one gene.
    ig <- graph_from_data_frame(set.connectivity.tbl[,c("set1", "set2")], directed = FALSE)

##    cat(paste0("Largest independent vertex set(s) size: ", independence.number(ig), "\n"))
    
    ## Calculate the independent vertex set
    ##    largest_ivs(ig)
    ivs(ig)
}

## independent.hallmark.sets <- find.largest.independent.sets(bindeaGsets[grepl(pattern="HALLMARK", names(bindeaGsets)) | (names(bindeaGsets) == "CIRC")])

independent.hallmark.sets <- find.largest.independent.sets(bindeaGsets[grepl(pattern="HALLMARK", names(bindeaGsets))])


pval.col <- "CMS2.pval"
pval.col <- "krasMT.pval"
tcga.tbl <- read.table("tcga-mcp-full-model-tbl.xls", sep="\t", header=TRUE)
kfs.tbl <- read.table("kfs-mcp-full-model-tbl.xls", sep="\t", header=TRUE)
tcga.tbl <- tcga.tbl[order(as.numeric(tcga.tbl[,pval.col])),]
tcga.tbl$dataset <- "tcga"
kfs.tbl$dataset <- "kfs"
common.cols <- intersect(colnames(tcga.tbl), colnames(kfs.tbl))

tbl <- rbind(tcga.tbl[,common.cols],kfs.tbl[,common.cols])
tbl <- subset(tbl, grepl(variable, pattern="HALLMARK"))
tbl$variable <- factor(tbl$variable, levels = unique(tbl$variable))
pval.col <- "CMS2.pval"
pval.col <- "krasMT.pval"
g <- ggplot(data = tbl)
g <- g + geom_point(aes(x = variable, y = -log10(krasMT.pval), col = dataset))
g <- g + geom_point(aes(x = variable, y = log(CMS2.pval, base=100000000), col = dataset))
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
print(g)

tbl <- rbind(tcga.tbl[,common.cols],kfs.tbl[,common.cols])
tbl <- subset(tbl, grepl(variable, pattern="HALLMARK"))
tbl$variable <- factor(tbl$variable, levels = unique(tbl$variable))
g1 <- ggplot(data = tbl)
g1 <- g1 + geom_point(aes(x = variable, y = -log10(krasMT.pval), col = dataset))
## g1 <- g1 + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
g1 <- g1 + theme(axis.text.x=element_blank())

g2 <- ggplot(data = tbl)
g2 <- g2 + geom_point(aes(x = variable, y = -log10(CMS2.pval), col = dataset))
## g2 <- g2 + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
g2 <- g2 + theme(axis.text.x=element_blank())

grid.arrange(g1,g2)

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

tbl <- rbind(tcga.tbl[,common.cols],kfs.tbl[,common.cols])
tbl <- subset(tbl, !grepl(variable, pattern="HALLMARK"))
tbl$variable <- factor(tbl$variable, levels = unique(tbl$variable))
g1 <- ggplot(data = tbl)
## g1 <- g1 + geom_point(aes(x = variable, y = -log10(krasMT.pval), col = dataset))
g1 <- g1 + geom_point(aes(x = variable, y = krasMT.pval, col = dataset))
## g1 <- g1 + scale_y_reverse() + scale_y_log10()
g1 <- g1 + scale_y_continuous(trans=reverselog_trans(10))
g1 <- g1 + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
## g1 <- g1 + theme(axis.text.x=element_blank())
g1 <- g1 + geom_hline(yintercept = 0.01)

g2 <- ggplot(data = tbl)
## g2 <- g2 + geom_point(aes(x = variable, y = -log10(CMS2.pval), col = dataset))
g2 <- g2 + geom_point(aes(x = variable, y = CMS2.pval, col = dataset))
g2 <- g2 + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
## g2 <- g2 + theme(axis.text.x=element_blank())
g2 <- g2 + scale_y_continuous(trans=reverselog_trans(10))
g2 <- g2 + geom_hline(yintercept = 0.01)

grid.arrange(g1,g2)


tbl <- rbind(tcga.tbl[,common.cols],kfs.tbl[,common.cols])
tbl <- subset(tbl, grepl(variable, pattern="HALLMARK"))
tbl$variable <- factor(tbl$variable, levels = unique(tbl$variable))
g1 <- ggplot(data = tbl)
## g1 <- g1 + geom_point(aes(x = variable, y = -log10(krasMT.pval), col = dataset))
g1 <- g1 + geom_point(aes(x = variable, y = krasMT.pval, col = dataset))
## g1 <- g1 + scale_y_reverse() + scale_y_log10()
g1 <- g1 + scale_y_continuous(trans=reverselog_trans(10))
g1 <- g1 + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
g1 <- g1 + theme(axis.text.x=element_blank())
g1 <- g1 + geom_hline(yintercept = 0.01)

g2 <- ggplot(data = tbl)
## g2 <- g2 + geom_point(aes(x = variable, y = -log10(CMS2.pval), col = dataset))
g2 <- g2 + geom_point(aes(x = variable, y = CMS2.pval, col = dataset))
g2 <- g2 + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
## g2 <- g2 + theme(axis.text.x=element_blank())
g2 <- g2 + scale_y_continuous(trans=reverselog_trans(10))
g2 <- g2 + geom_hline(yintercept = 0.01)

grid.arrange(g1,g2)


tbl <- rbind(tcga.tbl[,common.cols],kfs.tbl[,common.cols])
tbl <- subset(tbl, !grepl(variable, pattern="HALLMARK"))
tbl$variable <- factor(tbl$variable, levels = unique(tbl$variable))
pval.col <- "CMS2.pval"
pval.col <- "krasMT.pval"
g <- ggplot(data = tbl)
g <- g + geom_point(aes(x = variable, y = -log10(krasMT.pval), col = dataset))
g <- g + geom_point(aes(x = variable, y = log(CMS2.pval, base=10000000), col = dataset))
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 0), axis.ticks.x=element_blank())
print(g)


plot.stacked.pval.figures <- function(tbl1, tbl2, main1, main2, text.size = NULL) {
    g1 <- ggplot(data = tbl1)
    ## g1 <- g1 + geom_point(aes(x = gene.set, y = -log10(pval), col = Direction))
    g1 <- g1 + geom_point(aes(x = gene.set, y = pval, col = Direction))
    g1 <- g1 + scale_y_continuous(trans=reverselog_trans(10))
    g1 <- g1 + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, hjust = 1), axis.ticks.x=element_blank())
    if(!is.null(text.size)) {
        g1 <- g1 + theme(axis.text.x=element_text(angle = 45, hjust = 1, size = text.size))
    }
    g1 <- g1 + geom_hline(yintercept = 0.01)
    g1 <- g1 + ylab(bquote(italic('p')-value))
    g1 <- g1 + ggtitle(main1) + theme(plot.title = element_text(hjust = 0.5))
    ## Grab the axis labels
    gxlabels <- gtable::gtable_filter(ggplotGrob(g1), "axis-b")
    ## Then remove it
    g1 <- g1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
    ## gr1 <- ggplotGrob(g1)
    
    g2 <- ggplot(data = tbl2)
    ## g2 <- g2 + geom_point(aes(x = gene.set, y = -log10(pval), col = Direction))
    g2 <- g2 + geom_point(aes(x = gene.set, y = pval, col = Direction))
    g2 <- g2 + scale_y_continuous(trans=reverselog_trans(10))
    g2 <- g2 + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, hjust = 1), axis.ticks.x=element_blank())
    if(!is.null(text.size)) {
        g2 <- g2 + theme(axis.text.x=element_text(angle = 45, hjust = 1, size = text.size))
    }
    g2 <- g2 + geom_hline(yintercept = 0.01)
    g2 <- g2 + ylab(bquote(italic('p')-value))
    g2 <- g2 + ggtitle(main2) + theme(plot.title = element_text(hjust = 0.5))

    gr1 <- ggplot_gtable(ggplot_build(g1))
    gr2 <- ggplot_gtable(ggplot_build(g2))
    
    ## Set the heights to be the same for both plots
    gr1$heights[which(gr1$layout$name == "panel")] <- unit(5, units="cm")
    gr2$heights[which(gr2$layout$name == "panel")] <- unit(5, units="cm")    
    ##    grid.arrange(gr1, gr2)
    return(list(gr1, gr2))    
}
    

gsva.es.lst <- list()
clin.es.lst <- list()

for(dataset in c("tcga", "kfs")) {
    ds <- tcga_expr
    if(dataset == "kfs") {
        ds <- kfsyscc_expr
    }
    tmp <- doGSVA(ds, clin, dataset, to.plot=unique(c(circ.signature, hallmark.gene.sets, immune.gene.sets, other.genes)))
    es <- tmp[["es"]]
    clin.es <- tmp[["clin"]]

    mcp <- NULL
    if(dataset == "kfs") {
        mcp <- MCPcounter.estimate(kfsyscc_microarray_expr, featuresType="affy133P2_probesets")
    } else {
        mcp <- MCPcounter.estimate(tcga_expr, featuresType="HUGO_symbols")        
    }
    mcp.populations <- rownames(mcp)
    mcp.population.labels <- unlist(lapply(mcp.populations, function(x) paste0(x, " Population Percentage")))
    ## rownames(mcp) <- unlist(lapply(rownames(mcp), function(x) paste0("MCP ", x)))

    common.cols <- intersect(colnames(mcp), colnames(es))
    es <- as.matrix(rbind(mcp[, common.cols], es[, common.cols]))

    gsva.es.lst[[dataset]] <- es
    clin.es.lst[[dataset]] <- clin.es
}

## Add neoantigens to tcga
clin <- clin.es.lst[["tcga"]]
rownames(clin) <- clin$sample
es <- gsva.es.lst[["tcga"]]
cols <- intersect(colnames(es), rownames(clin))
es <- rbind(es, neoantigens=rep(NA, ncol(es)))
es["neoantigens",cols] <- clin[cols,"neoantigens"]
clin.es.lst[["tcga"]] <- es

neoantigen.signature <- "neoantigens"
neoantigen.label <- "Log10 ( # Neoantigens + 1 )"

kras.cms.res <- list()
kras.cms.interaction.res <- list()
kras.cms.within.cms.res <- list()
kras.res <- list()
cms.res <- list()
site.res <- list()
msi.res <- list()

for(dataset in c("tcga", "kfs")) {
    ds <- tcga_expr
    if(dataset == "kfs") {
        ds <- kfsyscc_expr
    }

    es <- gsva.es.lst[[dataset]]
    clin.es <- clin.es.lst[[dataset]]

    to.plot <- c(circ.signature, immune.populations, hallmark.gene.sets, other.genes, neoantigen.signature)
    ylabels <- c(circ.signature.label, immune.population.labels, hallmark.ylabels, other.labels, neoantigen.label)

    res <- doKrasCMSInteractionAnalysis(es, clin.es, to.plot=to.plot, analysis.name = dataset, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.interaction.res[[dataset]] <- res$tests
    kras.cms.interaction.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms-interaction.tsv"), kras.cms.interaction.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    res <- doKrasCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.res[[dataset]] <- res$tests
    kras.cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms.tsv"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doKrasCMSWithinCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.cms.within.cms.res[[dataset]] <- res$tests
    kras.cms.within.cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras-cms-within-cms.tsv"), kras.cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.res[[dataset]] <- res$tests
    kras.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras.tsv"), kras.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    cms.res[[dataset]] <- res$tests
    cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-cms.tsv"), cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doSiteAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    site.res[[dataset]] <- res$tests
    site.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-site.tsv"), site.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    if(dataset != "kfs") {
        res <- doMSIAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
        msi.res[[dataset]] <- res$tests
        msi.tbl <- res$tests
        write.table(file = paste0(dataset, "-immune-hallmark-msi.tsv"), msi.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    
    flag <- to.plot %in% c("CIRC", "STAT1", "CIITA", "CXCL10", "Cytotoxic cells", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", neoantigen.signature)
    to.plot <- to.plot[flag]
    ylabels <- ylabels[flag]
    plotKrasCMSAnalysis(es, clin.es, kras.cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

    plotKrasCMSWithinCMSAnalysis(es, clin.es, kras.cms.within.cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))    

    plotKrasAnalysis(es, clin.es, kras.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

    plotCMSAnalysis(es, clin.es, cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

    plotSiteAnalysis(es, clin.es, site.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))        

    if(dataset != "kfs") {
        plotMSIAnalysis(es, clin.es, msi.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))    
    }
    
    forest.tbl <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
    
}

gene.sets.to.analyze <- c(circ.signature, immune.populations, hallmark.gene.sets, other.genes, neoantigen.signature)
gene.sets.to.analyze.labels <- c(circ.signature.label, immune.population.labels, hallmark.ylabels, other.labels, neoantigen.label)

flag <- gene.sets.to.analyze %in% c("CIRC", "STAT1", "CIITA", "CXCL10", "Cytotoxic cells", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", neoantigen.signature)
gene.sets.to.plot <- gene.sets.to.analyze[flag]
gene.sets.to.plot.labels <- gene.sets.to.plot.labels[flag]

dataset <- "tcga"
do.all.analyses(gsva.es.lst[[dataset]], clin.es.lst[[dataset]], dataset, gene.sets.to.analyze = gene.sets.to.analyze, gene.sets.to.plot = gene.sets.to.plot, gene.sets.to.plot.labels = gene.sets.to.plot.labels)

dataset <- "kfs"
do.all.analyses(gsva.es.lst[[dataset]], clin.es.lst[[dataset]], dataset, gene.sets.to.analyze = gene.sets.to.analyze, gene.sets.to.plot = gene.sets.to.plot, gene.sets.to.plot.labels = gene.sets.to.plot.labels)


## Redo the above analyses for MSI-annotated KFS

dataset <- "kfs"
clin.tmp <- clin.inferred.msi.lst[[dataset]]
clin.tmp$msi <- clin.tmp$msi.inferred
do.all.analyses(gsva.es.lst[[dataset]], clin.tmp, paste0(dataset, "-msi"), gene.sets.to.analyze = gene.sets.to.analyze, gene.sets.to.plot = circ.signature, gene.sets.to.plot.labels = circ.signature.label)

## Redo the above analysis for MMS-only

dataset <- "tcga"
clin.tmp <- clin.inferred.msi.lst[[dataset]]
clin.tmp$msi <- clin.tmp$msi.inferred
clin.tmp <- clin.tmp[!is.na(clin.tmp$msi) & (clin.tmp$msi == "mss") & !is.na(clin.tmp$cms_label) & (clin.tmp$cms_label != "CMS1"), ]
do.all.analyses(gsva.es.lst[[dataset]], clin.tmp, paste0(dataset, "-mss-only"), gene.sets.to.analyze = gene.sets.to.analyze, gene.sets.to.plot = circ.signature, gene.sets.to.plot.labels = circ.signature.label)

dataset <- "kfs"
clin.tmp <- clin.inferred.msi.lst[[dataset]]
clin.tmp$msi <- clin.tmp$msi.inferred
clin.tmp <- clin.tmp[!is.na(clin.tmp$msi) & (clin.tmp$msi == "mss") & !is.na(clin.tmp$cms_label) & (clin.tmp$cms_label != "CMS1"), ]
do.all.analyses(gsva.es.lst[[dataset]], clin.tmp, paste0(dataset, "-mss-only"), gene.sets.to.analyze = gene.sets.to.analyze, gene.sets.to.plot = circ.signature, gene.sets.to.plot.labels = circ.signature.label)


## Plot the p-values for the immune and hallmark signatures

##tcga.tbl <- read.table("tcga-immune-hallmark-kras-cms.tsv", sep="\t", header=TRUE)
##kfs.tbl <- read.table("kfs-immune-hallmark-kras-cms.tsv", sep="\t", header=TRUE)

tcga.tbl <- kras.cms.res[["tcga"]]
kfs.tbl <- kras.cms.res[["kfs"]]

tcga.tbl$dataset <- "TCGA"
kfs.tbl$dataset <- "KFSYSCC"

common.cols <- intersect(colnames(tcga.tbl), colnames(kfs.tbl))
tbl <- rbind(tcga.tbl[, common.cols], kfs.tbl[, common.cols])

tbl$gene.set <- gsub(tbl$gene.set, pattern="HALLMARK_", replacement="")

pval.tbl <- tbl[, c("dataset", "gene.set", colnames(tbl)[grepl(colnames(tbl), pattern=".pval")])]
colnames(pval.tbl) <- gsub(colnames(pval.tbl), pattern=".pval", replacement="")
pval.tbl <- gather_(data = pval.tbl, key_col = "comparison", value_col = "pval", colnames(pval.tbl)[!(colnames(pval.tbl) %in% c("dataset", "gene.set"))])

dir.tbl <- tbl[, c("dataset", "gene.set", colnames(tbl)[grepl(colnames(tbl), pattern=".dir")])]
colnames(dir.tbl) <- gsub(colnames(dir.tbl), pattern=".dir", replacement="")
dir.tbl <- gather_(data = dir.tbl, key_col = "comparison", value_col = "dir", colnames(dir.tbl)[!(colnames(dir.tbl) %in% c("dataset", "gene.set"))])

tbl <- merge(pval.tbl, dir.tbl, by = c("dataset", "gene.set", "comparison"))

immune.flag <- tbl$gene.set %in% immune.populations
hallmark.flag <- tbl$gene.set %in% gsub(hallmark.gene.sets, pattern="HALLMARK_", replacement="")

tbl$Direction <- unlist(lapply(tbl$dir, function(x) ifelse(x > 0, "Increased", "Decreased")))

## Order by maximum CMS2WT.vs.CMS2MT.pval in TCGA
stbl <- ddply(tbl[tbl$dataset == "TCGA", c("gene.set", "pval")], .variables = "gene.set", .fun = function(df) max(as.numeric(df$pval)))
colnames(stbl) <- c("gene.set", "max.pval")
stbl <- stbl[order(as.numeric(stbl$max.pval), decreasing=FALSE),]

sorted.levels <- unique(stbl$gene.set)
tbl$gene.set <- factor(tbl$gene.set, levels = sorted.levels)

tbl$pval <- as.numeric(tbl$pval)
tcga.flag <- tbl$dataset == "TCGA"
tcga.tbl <- tbl[hallmark.flag & tcga.flag, ]
kfs.tbl <- tbl[hallmark.flag & !tcga.flag, ]
## pdf("hallmark-pvalues.pdf", onefile=FALSE)
l <- plot.stacked.pval.figures(tcga.tbl, kfs.tbl, "TCGA", "KFSYSCC", text.size = 5)
plot_grid(l[[1]], l[[2]], labels = c("A", "B"), nrow=2)
ggsave("hallmark-pvalues.tiff")
## d <- dev.off()

tcga.flag <- tbl$dataset == "TCGA"
tcga.tbl <- tbl[immune.flag & tcga.flag, ]
kfs.tbl <- tbl[immune.flag & !tcga.flag, ]
## pdf("immune-pvalues.pdf", onefile=FALSE)
l <- plot.stacked.pval.figures(tcga.tbl, kfs.tbl, "TCGA", "KFSYSCC")
plot_grid(l[[1]], l[[2]], labels = c("A", "B"), nrow=2)
ggsave("immune-pvalues.tiff")
## d <- dev.off()



## WORKING

## Plot the heatmap of the Bindea gene sets and CIRC
plot.pvals <- TRUE
log.transform <- TRUE

tcga.gene.set.kras.tbl <- kras.res[["tcga"]][c(immune.populations),]
kfs.gene.set.kras.tbl <- kras.res[["kfs"]][c(immune.populations),]
merged.kras.tbl <- merge(tcga.gene.set.kras.tbl, kfs.gene.set.kras.tbl, by = "gene.set", all = TRUE, suffixes = c(".tcga", ".kfs"))
rownames(merged.kras.tbl) <- merged.kras.tbl$gene.set

merged.kras.tbl$val.tcga <- merged.kras.tbl$pval.tcga
merged.kras.tbl$val.kfs <- merged.kras.tbl$pval.kfs  
if(!plot.pvals) {
    merged.kras.tbl$val.tcga <- merged.kras.tbl$apval.tcga
    merged.kras.tbl$val.kfs <- merged.kras.tbl$apval.kfs  
}

merged.kras.tbl$val.tcga <- as.numeric(merged.kras.tbl$val.tcga)
merged.kras.tbl$val.kfs <- as.numeric(merged.kras.tbl$val.kfs)

## Make sure everything is in the direction of circ if it significant
for(suffix in c("tcga", "kfs")) {
    pval.col <- paste0("val.", suffix)
    dir.col <- paste0("mutUp.", suffix)
    flag <- !is.na(merged.kras.tbl[,pval.col]) & (merged.kras.tbl[,pval.col] < 1)
    ##      inconsistent <- merged.kras.tbl[flag, dir.col] != merged.kras.tbl["CIRC", dir.col]
    inconsistent <- merged.kras.tbl[flag, dir.col] != FALSE
    if(any(inconsistent)) {
        cat("The following circ components are in the direction opposite the circ metagene\n")
        write.table(merged.kras.tbl[inconsistent,], row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    }
}
## Make mutUp == FALSE be -log10 pvalues and mutUp == TRUE by log10 pvalues
if(log.transform) {
    merged.kras.tbl$val.tcga <- -log10(merged.kras.tbl$val.tcga)
} else {
    merged.kras.tbl$val.tcga <- 1 - merged.kras.tbl$val.tcga
}
flag <- !is.na(merged.kras.tbl$mutUp.tcga) & (merged.kras.tbl$mutUp.tcga == TRUE)
merged.kras.tbl$val.tcga[flag] <- - merged.kras.tbl$val.tcga[flag]

if(log.transform) {  
    merged.kras.tbl$val.kfs <- -log10(merged.kras.tbl$val.kfs)
} else {
    merged.kras.tbl$val.kfs <- 1 - merged.kras.tbl$val.kfs
}
flag <- !is.na(merged.kras.tbl$mutUp.kfs) & (merged.kras.tbl$mutUp.kfs == TRUE)  
merged.kras.tbl$val.kfs[flag] <- - merged.kras.tbl$val.kfs[flag]

merged.kras.tbl <- merged.kras.tbl[,c("gene.set", "val.tcga", "val.kfs")]  
colnames(merged.kras.tbl) <- c("gene.set", "TCGA", "KFSYSCC") 
melted.kras.tbl <- melt(merged.kras.tbl)
colnames(melted.kras.tbl) <- c("GeneSet", "Dataset", "value")

tiff("immune-gene-set-heatmap.tiff")
g <- plot.pval.heatmap(melted.kras.tbl, y = "GeneSet", x = "Dataset", value = "value", xlab = "Data Set", ylab="Gene Set", use.italics = FALSE, log.transform = log.transform)
grid.draw(g)
dev.off()


q(status=0)


library(Biobase)
library(GEOquery)

## Load the metastatic dataset
gset <- getGEO("GSE5851", GSEMatrix =TRUE, getGPL=FALSE)
gset <- gset[[1]]
exprs(gset)
pData(gset)
map <- pData(gset)[,c("title"),drop=FALSE]
map$sample_id <- rownames(map)

## Code taken/modified from Mike Mason
## https://github.com/Sage-Bionetworks/CRC_DownStreamAnalysis/blob/master/2017ASCOanalysis.R
processFromGEO <- function(geo, GSMS=NULL, method="rma", url, download.path) ## GSMS is for subsetting to specific samples if necessary
{
    library("pd.hg.u133a.2")
    library("hgu133a2frmavecs")
    cat(paste0("Processing ", url, "\n"))
    file <- unlist(strsplit(url, split="/"))
    gse <- file[length(file)-1]
    file <- file[length(file)]
    full.path <- paste0(download.path, "/", file)
    if(!file.exists(full.path)) {
        download.file(url, full.path)
    }
    gse.dir <- paste0(download.path, "/", gse)
    system(paste0("mkdir ", gse.dir))
    system(paste0("tar -xvf ", full.path, " -C ", gse.dir))
    celfiles      <- list.files(gse.dir, full.names = T)
    celfilesShort <- list.files(gse.dir, full.names = F)
    celInds       <- grep("\\.CEL", celfiles)  
    celfiles      <- celfiles[celInds]
    celfilesShort <- celfilesShort[celInds]
    locGSMs = gsub("\\..*$","",celfilesShort)
    if(!is.null(GSMS)){ celfiles = celfiles[locGSMs %in% GSMS] }
    if(method == "rma")
    {
        library(oligo)                
        temp     = read.celfiles(celfiles)
        exprData = oligo::rma(temp)
    }else if(method ==  "frma")
    {
        library(frma)
        library(oligo)        
        temp = read.celfiles(filenames = celfiles); print(class(temp))
        if(class(temp) == "ExpressionFeatureSet"){temp <- read.affybatch(filenames = celfiles)}
        exprData = frma(temp)
    }
    rm(temp);gc()
    return(exprData)
}


url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/GSE5851/GSE5851_RAW.tar"
kf.expr.geo <- processFromGEO(geo = "GSE5851", GSMS=NULL, method="frma", url, synapse.repo.dir)

## Load the table with KRAS status
## Supp Table 1 from online table 1

kf.tbl <- read.table(paste0(synapse.repo.dir, "/Khambata-Ford-Table1.tsv"), sep="\t", header=TRUE, as.is=TRUE, comment.char="")
kf.tbl$id <- gsub(kf.tbl$Affymetrix.ID, pattern="FL", replacement="")

kf.clin <- merge(map, kf.tbl, by.x="title", by.y="id", all=FALSE)

## The majority of the mets are liver mets (61 of 80)--to reduce heterogeneity, let's just keep those.
kf.clin$Biopsy.tissue <- gsub(kf.clin$Biopsy.tissue, pattern="Liver ", replacement="Liver")
table(kf.clin$Biopsy.tissue)
kf.clin <- kf.clin[kf.clin$Biopsy.tissue == "Liver",]

rownames(kf.clin) <- kf.clin$sample_id
kf.clin$kras <- kf.clin$KRAS.Mutation
kras.mut.flag <- grep(kf.clin$kras, pattern="c.")
kf.clin$kras[kras.mut.flag] <- "MT"

## Excludde NAs:
kf.clin <- kf.clin[(kf.clin$kras == "MT") | (kf.clin$kras == "WT"),]

kf.expr <- exprs(gset)
kf.expr <- kf.expr[, colnames(kf.expr) %in% kf.clin$sample_id]

kf.clin <- kf.clin[colnames(kf.expr),]

kf.expr <- convert_hg133_probes_to_genes(kf.expr, method = "mean")

## Set up kf.clin like clin
kras.mut <- kf.clin$kras == "MT"
kf.clin$kras[kras.mut] <- 1
kf.clin$kras[!kras.mut] <- 0
kf.clin$kras <- as.numeric(kf.clin$kras)
kf.clin$sample <- kf.clin$sample_id
## Fake CMS for the moment
kf.clin$cms_label <- "CMS1"
kf.clin$msi <- NA
kf.clin$site <- NA
kf.clin$neoantigens <- NA

ds <- kf.expr
dataset <- "kf"
## We don't have NRAS mutation data.
tmp <- doGSVA(ds, kf.clin, dataset, to.plot=unique(c(circ.signature, hallmark.gene.sets, immune.gene.sets)), exclude.nras.mutants = FALSE)
es <- tmp[["es"]]
clin.es <- tmp[["clin"]]

res <- doKrasAnalysis(es, clin.es, paste0(dataset, "-kras"), to.plot=c("CIRC"), num.permutations = 0, calc.fdr = FALSE)
kras.res[[dataset]] <- res$tests
tbl <- res$tests
write.table(file = paste0(dataset, "-immune-hallmark-kras.tsv"), tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

to.plot <- "CIRC"
ylabels <- "CIRC Enrichment Score"
plotKrasAnalysis(es, clin.es, tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = paste0(dataset, "-kras-cms"), plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = "KF")


## install CMS classifier
## run classifier
## look at distribution of kf relative to TCGA and KFS
## doKRASCMS
## plotKRASCMS
## 

## Clean up run.R
## Make effect size plot
##    clear up; stack; add kfs
## Make plots for bindea
## Make plots sfor mcp counter



## 1. Add bootstrap to kras + cms (bootstrap within group to calculate effect size) -- i.e., bootstrap under the alternative hypothesis
## 2. Use 2-step in kras + cms
## 3. Run _separately_ on
##    3a. STAT1, and pathways
##    3b. bindea
##    3c. mcp counter
## 4. Create plot showing
##    4a. overall sig
##    4b. boxplot of effect size (which effect) in color
