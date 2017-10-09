my.convert.vec <- function(x, species, from, to) {

    db <- graphite:::loadDb(species)
    to <- graphite:::destIdentifier(to, db)

    tbl <- graphite:::lookupKeys(db, x, from, to)

    if (is.null(tbl)) {
      return(x)
    }

##    print(tbl)
    r <- as.vector(unlist(tbl))
    names(r) <- names(tbl)
##    print(r)
    r <- r[x]
    if(length(x) != length(r)) {
        save(r, file="r.Rd")
        save(x, file="x.Rd")
        cat(paste0("length(x) = ", length(x), "\n"))
        cat(paste0("length(r) = ", length(r), "\n"))        
        stop("Trouble in my.convert.vec\n")
    }
    ##    r <- as.vector(r[x])
    d <- duplicated(r)
    if(any(d)) {

    }
        
    nxt <- 1
    for(i in 1:length(r)) {
        if(is.na(r[i])) {
            nai <- paste0("NA", nxt)
            r[i] <- nai
            nxt <- nxt + 1
        }
    }
##    print(length(r))
##    print(length(x))
##    print(r)
    return(r)
##   return(tbl[,2])

#df <- data.frame(src=as.character(x), stringsAsFactors=FALSE)
    #r <- my.convert(db, df, "src", from, to)
    #return(r$src)
}

## ToPASeq:::makeNodeTable
## Sometimes we have non-unique mapping from one gene id to a second type of gene id.
## These non-unique names/symbols/ids would be collapsed by nodes(g)
## Hence, just take the list of node names instead of g, gc, and gp--
## these may be duplicated
my.makeNodeTable <-
function (g.nodes, gc.nodes, gp.nodes, breaks, deg.table, sigpal, tsig.whole, tsig, 
    mis.col = "white", p.th = 0.05, col.lim = NULL) 
{
    rownames(deg.table) <- as.character(deg.table[, 1])
    namesOrig <- g.nodes
    namesConv <- gc.nodes
    namesPlot <- gp.nodes
    isPresent <- gc.nodes %in% as.character(deg.table[, 1])
    nodeStat <- deg.table[gc.nodes, 2]
    nodePval <- deg.table[gc.nodes, "pval"]
    isSig <- nodePval < p.th
    border <- ifelse(isSig | is.na(isSig), "white", "black")
    if (is.null(col.lim)) 
        b <- ToPASeq:::makeBreaks(nodeStat[!is.na(nodeStat)], breaks[1])
    else b <- ToPASeq:::makeBreaks(col.lim, breaks[1])
    nodeCol <- as.character(cut(nodeStat, b, labels = sigpal(length(b) - 
        1), include.lowest = TRUE))
    nodeCol[!isPresent] <- mis.col
    if (all(tsig.whole == 1)) {
        tsigCat <- factor(tsig)
    }
    else {
        b2 <- ToPASeq:::makeBreaks(tsig.whole, breaks[2])
        tsigCat <- cut(tsig, b2)
    }
    nodeTable <- data.frame(namesOrig = namesOrig, namesConv = namesConv, 
        namesPlot = namesPlot, isPresent = isPresent, nodeStat = nodeStat, 
        nodePval = nodePval, isSig = isSig, border = border, 
        nodeCol = nodeCol, tsigCat = tsigCat, stringsAsFactors = FALSE)
    return(nodeTable)
}

### ToPASeq:::adjustAttr
## Changed this function not to turn off edges between insignificant nodes
my.adjustAttr <-
function (g, NodeTable, EdgeList, stats, cols, remNodes = NULL) 
{
    if (is.character(NodeTable$isSig)) {
        sig <- NodeTable$isSig == "TRUE"
        names(sig) <- NodeTable$namesPlot
    }
    else {
        sig <- NodeTable$isSig
        names(sig) <- nodes(g)
    }
    E <- edges(g)
    selectValid <- selectInvalid <- rep(FALSE, nrow(E))
    th <- 0
    selectValid[stats[E[, 1]] > th & stats[E[, 2]] > th & (regexpr("activation", 
        as.character(E[, 4])) != -1 | regexpr("expression", as.character(E[, 
        4])) != -1)] <- TRUE
    selectValid[stats[E[, 1]] > th & stats[E[, 2]] < -th & (regexpr("inhibition", 
        as.character(E[, 4])) != -1 | regexpr("repression", as.character(E[, 
        4])) != -1)] <- TRUE
    selectValid[stats[E[, 1]] < -th & stats[E[, 2]] > th & (regexpr("inhibition", 
        as.character(E[, 4])) != -1 | regexpr("repression", as.character(E[, 
        4])) != -1)] <- TRUE
    selectValid[stats[E[, 1]] < -th & stats[E[, 2]] < -th & (regexpr("activation", 
        as.character(E[, 4])) != -1 | regexpr("expression", as.character(E[, 
        4])) != -1)] <- TRUE
    ## Commented this out so that we still color edges between insignificant nodes
    ## selectValid[!sig[E[, 1]] | !sig[E[, 2]]] <- FALSE

    selectInvalid[stats[E[, 1]] > th & stats[E[, 2]] > th & (regexpr("inhibition", 
        as.character(E[, 4])) != -1 | regexpr("repression", as.character(E[, 
        4])) != -1)] <- TRUE
    selectInvalid[stats[E[, 1]] > th & stats[E[, 2]] < -th & 
        (regexpr("activation", as.character(E[, 4])) != -1 | 
            regexpr("expression", as.character(E[, 4])) != -1)] <- TRUE
    selectInvalid[stats[E[, 1]] < -th & stats[E[, 2]] > th & 
        (regexpr("activation", as.character(E[, 4])) != -1 | 
            regexpr("expression", as.character(E[, 4])) != -1)] <- TRUE
    selectInvalid[stats[E[, 1]] < -th & stats[E[, 2]] < -th & 
        (regexpr("inhibition", as.character(E[, 4])) != -1 | 
            regexpr("repression", as.character(E[, 4])) != -1)] <- TRUE
    ## Commented this out so that we still color edges between insignificant nodes
    ## selectInvalid[!sig[E[, 1]] | !sig[E[, 2]]] <- FALSE
    
    valid.edges <- paste(E[, 1], E[, 2], sep = "~")[selectValid]
    not.valid.edges <- paste(E[, 1], E[, 2], sep = "~")[selectInvalid]
    valid.nodes <- unique(c(E[, 1][selectValid | selectInvalid], 
        E[, 2][selectValid | selectInvalid]))

    indecisive.edges <- paste(E[, 1], E[, 2], sep = "~")[!selectValid & 
        !selectInvalid]
    edg.col <- c(setNames(rep(cols[1], length(valid.edges)), 
        valid.edges), setNames(rep(cols[2], length(not.valid.edges)), 
        not.valid.edges), setNames(rep(cols[3], length(indecisive.edges)), 
        indecisive.edges))
    EdgeList$Ecol <- edg.col
    if (!is.null(remNodes)) {
        if (!is.character(remNodes)) 
            stop("Invalid color specification for nodes with neutral interacions or with small expression change.")
        valid.nodes <- valid.nodes[valid.nodes %in% NodeTable$namesPlot]
        cat("namesPlot\n")
        invalid.nodes <- as.character(NodeTable$namesPlot[!NodeTable$namesPlot %in% 
            valid.nodes])
        nodecol <- as.character(NodeTable$nodeCol)
        nodecol[NodeTable$namesPlot %in% invalid.nodes] <- remNodes
        NodeTable$nodeCol <- nodecol
    }
    out <- list(NodeTable, EdgeList)
    return(out)
}

if(FALSE){
selectValid[ stats[E[,1]] >  th & stats[E[,2]] >  th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1 )]<-TRUE
selectValid[ stats[E[,1]] >  th & stats[E[,2]] < -th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1)]<-TRUE
selectValid[ stats[E[,1]] < -th & stats[E[,2]] >  th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1)] <-TRUE
selectValid[ stats[E[,1]] < -th & stats[E[,2]] < -th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1)]<-TRUE
selectValid[!sig[E[,1] ] | !sig[E[,2] ] ]<-FALSE


selectInvalid[ stats[E[,1]] >  th & stats[E[,2]] >  th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1 )]<-TRUE
selectInvalid[ stats[E[,1]] >  th & stats[E[,2]] < -th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1)]<-TRUE
selectInvalid[ stats[E[,1]] < -th & stats[E[,2]] >  th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1)] <-TRUE
selectInvalid[ stats[E[,1]] < -th & stats[E[,2]] < -th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1)]<-TRUE
selectInvalid[!sig[E[,1] ] | !sig[E[,2] ] ]<-FALSE

valid.edges<-paste(E[,1], E[,2], sep="~")[selectValid]
not.valid.edges<-paste(E[,1], E[,2], sep="~")[selectInvalid]

valid.nodes<-unique(c(E[,1][selectValid | selectInvalid], E[,2][selectValid | selectInvalid]))
indecisive.edges<-paste(E[,1], E[,2], sep="~")[!selectValid & !selectInvalid] 

edg.col<-c(setNames(rep(cols[1],length(valid.edges)), valid.edges), 
           setNames(rep(cols[2],length(not.valid.edges)), not.valid.edges),
           setNames(rep(cols[3],length(indecisive.edges)), indecisive.edges) ) 
EdgeList$Ecol<-edg.col

if (!is.null(remNodes)) {
if (!is.character(remNodes)) stop("Invalid color specification for nodes with neutral interacions or with small expression change.")
valid.nodes<-valid.nodes[valid.nodes %in% NodeTable$namesPlot]
invalid.nodes<-as.character(NodeTable$namesPlot[ ! NodeTable$namesPlot %in% valid.nodes])
nodecol<-as.character(NodeTable$nodeCol)
nodecol[NodeTable$namesPlot %in% invalid.nodes]<-remNodes
NodeTable$nodeCol<-nodecol
}
}

### ToPASeq:::plot.topResult
## Added remNodes and passed through to setAttr
## Added p.th and passed through makeNodeTable
## Evidently (by playing with p.th) dashed lines mean insignificant
my.topaseq.plot <-
function (x, which, graphs, stats = "logFC", convert = TRUE, 
    IDs = "ENTREZID", graphIDs = "SYMBOL", col.lim = NULL, reduction = list(), 
    agg.fun = function(x) mean(x, na.rm = TRUE), logical = TRUE, 
    sig.th = 0.1, title = TRUE, cex.main = 1, breaks = c(100, 
        5), pallete.colors = c("blue", "white", "red"), na.col = "grey87", 
    cli.color = "red", layout = "dot", nodesize = 1, fontsize = 14, 
    alpha = 0.05, add.legend = TRUE, statName = "Log fold change", 
    cex.legend = 1, remNodes = "white", p.th = 0.05, ...) 
{

    res <- x
    g <- graphs[[which]]
    ## OLD code
    if (convert) 
        gc <- convertIdentifiers(graphs[[which]], IDs)
    else gc <- g
    gp <- convertIdentifiers(graphs[[which]], graphIDs)
    ## Begin bug fix
    g.nodes <- nodes(graphs[[which]])
    gc.nodes <- g.nodes
    gp.nodes <- g.nodes
    if(TRUE) {
        orig.g <- graphs[[which]]
        if(convert) {
            gc <- orig.g
            gc.nodes <- my.convert.vec(nodes(gc), orig.g@species, from=orig.g@identifier, to=IDs)
            ##            nodes(gc) <- new.nodes
            gc <- convertIdentifiersByVector(graphs[[which]], gc.nodes, "unknown")
        }
        gp <- graphs[[which]]
###    new.nodes <- my.convert.vec(g.nodes, orig.g@species, from=orig.g@identifier, to=gp@identifier)
        gp.nodes <- my.convert.vec(nodes(gp), orig.g@species, from=orig.g@identifier, to=graphIDs)    
        ##    save(new.nodes, file="nn.Rd")
        ##    save(g.nodes, file="gn.Rd")
        ## nodes(gp) <- new.nodes
        gp <- convertIdentifiersByVector(graphs[[which]], gp.nodes, "unknown")
    }
    ## End bug fix
    deg.table <- x$degtable
    sigpal <- colorRampPalette(pallete.colors)
    na.col <- colorRampPalette(c(na.col))(1)
    defaultEdgeAttrs <- ToPASeq:::makeDefaultEdgeData()
    if ("topResultC" %in% class(res)) {
        cliq <- res$topo.sig[[which]]
        NodeTable <- ToPASeq:::makeNodeTable(g, gc, gp, breaks, deg.table, 
            sigpal, tsig.whole = 1, tsig = 1, mis.col = na.col, 
            p.th = alpha, col.lim = col.lim)
        EdgeList <- ToPASeq:::makeEdgeList(gp, defaultEdgeAttrs)
        cols <- cli.color
        att <- ToPASeq:::adjustAttrCli(gc, NodeTable, EdgeList, cliq[[1]], 
            cliq[[2]], cols, alpha, remNodes = NULL)
        xxg <- ToPASeq:::renderOrig(gp, NodeTable, EdgeList, nodesize, 
            fontsize)
        xxred <- ToPASeq:::renderReduced(gp, reduction, att[[1]], att[[2]], 
            xxg, nodesize, fontsize, agg.fun)
        ToPASeq:::drawGraph(xxred, res, which, NodeTable, nodesize, fontsize, 
            statName = statName, cex.main = cex.main, col.lim = col.lim, 
            breaks, legend = add.legend)
    }
    if ("topResultW" %in% class(res)) {
        if (is.null(dim(res$topo.sig[[which]]))) 
            tsig <- res$topo.sig[[which]]
        else tsig <- res$topo.sig[[which]][1, ]
        if (length(tsig) == 0) 
            stop("This pathway was not analysed")
        if (is.null(dim(res$topo.sig[[which]]))) 
            tsig.whole <- unlist(sapply(res$topo.sig, function(x) x))
        else tsig.whole <- unlist(sapply(res$topo.sig, function(x) x[1, 
                                                                     ]))
    }
    if ("topResultE" %in% class(res)) {
        nod <- nodes(graphs[[which]])
        tsig <- setNames(rep(1, length(nod)), nod)
        tsig.whole <- rep(1, 100)
    }
    if ("topResultW" %in% class(res) | "topResultE" %in% class(res)) {
        ## NodeTable <- ToPASeq:::makeNodeTable(g, gc, gp, breaks, deg.table, sigpal, tsig.whole, tsig, col.lim = col.lim, p.th = p.th)
        NodeTable <- my.makeNodeTable(as.vector(g.nodes), as.vector(gc.nodes), as.vector(gp.nodes), breaks, deg.table, sigpal, tsig.whole, tsig, col.lim = col.lim)
        EdgeList <- ToPASeq:::makeEdgeList(gp, defaultEdgeAttrs)
        if (length(reduction) > 0) {
            gpr <- reduceGraph(gp, reduction)
        }
        else gpr <- gp
        NodeTable.red <- ToPASeq:::applyReduction(reduction, NodeTable, 
            agg.fun)
        EdgeList.red <- ToPASeq:::makeEdgeList(gpr, defaultEdgeAttrs)
        if (logical) {
            stats <- setNames(as.numeric(NodeTable.red$nodeStat), 
                              NodeTable.red$namesPlot)
            ## Changed this--added remNodes to this function
            ## att <- ToPASeq:::adjustAttr(gpr, NodeTable.red, EdgeList.red, stats, cols = c("black", "red", "grey87"), remNodes = "white")
            ## att <- ToPASeq:::adjustAttr(gpr, NodeTable.red, EdgeList.red, stats, cols = c("black", "red", "grey87"), remNodes = NULL)
            att <- my.adjustAttr(gpr, NodeTable.red, EdgeList.red, stats, cols = c("black", "red", "grey87"), remNodes = NULL)                        
        }
        else att <- list(NodeTable, EdgeList)
        xxg <- ToPASeq:::renderOrig(gp, NodeTable, EdgeList, nodesize, 
            fontsize)
        xxred <- ToPASeq:::renderReduced(gp, reduction, att[[1]], att[[2]], 
            xxg, nodesize = nodesize, fontsize = fontsize, agg.fun = agg.fun)
        ToPASeq:::drawGraph(xxred, res, which, NodeTable, nodesize, fontsize, 
            statName = statName, cex.main = cex.main, col.lim = col.lim, 
            breaks = breaks, sigpal = sigpal, legend = add.legend, 
            cex.legend = cex.legend)
    }
}
