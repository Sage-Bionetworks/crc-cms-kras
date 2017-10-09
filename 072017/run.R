
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
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))
              
num.processes <- 1

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapse.repo.dir <- "../input/"

add.figure.labels <- TRUE

set.seed(1234)

text.size <- 25
star.size <- 12
ns.size <- 6
gp.star.size <- 25
gp.ns.size <- 12

gp.star.size <- 35
gp.ns.size <- 22

system("mkdir output/")

           
source("doanalysis.R")

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

cat("Reading Bindea population definitions (i.e., genes in each immune population)\n")
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
print(bindeaCellTypes)
bindeaGsets <- lapply(bindeaCellTypes, function(x){
    return(as.character(unique(bindeaTbl$Symbol[bindeaTbl$CellType==x])))
})
names(bindeaGsets) <- bindeaCellTypes

## bindeaGsets <- bindeaGsets[sapply(bindeaGsets, length) > 10]
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

## Immune sets used in CMS Nat Med paper
immune.sets <- c("IMMUNE_ESTIMATE", "IMMUNE_RESP_GO_BP", "PD1_REACTOME", "IMMUNE_CD8MACRO_GALON", "IMMUNE_TH1_GALON", "IMMUNE_NKC_BREAST", "IMMUNE_THF_BREAST", "IMMUNE_TH17_GOUNARI", "IMMUNE_TREG_LUCAS", "IMMUNE_MDSC_ALBELDA")

for(immune.set in immune.sets) {
    bindeaGsets[[immune.set]] <- gene.sets$genesets[[which(gene.sets$geneset.names == immune.set)]]
}

immune.set.gene.symbols <- unique(as.vector(unlist(lapply(bindeaGsets[c(immune.sets,"CIRC")], function(x) as.vector(x)))))

## Download biocarta, kegg, reactome, and hallmark data sets
## HALLMARK
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/h.all.v6.0.symbols.gmt
## HALLMARK_INTERFERON_GAMMA_RESPONSE

## BIOCARTA
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.biocarta.v6.0.symbols.gmt

## KEGG
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.kegg.v6.0.symbols.gmt

## REACTOME
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.reactome.v6.0.symbols.gmt

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

to.plot <- c("CIRC", "neoantigens", "BATF3", "THBD", "CCL4", "ATF3", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_INTERFERON_GAMMA_SIGNALING", "IFNG")
ylabels <- c("CIRC Enrichment Score", "Log10 ( # Neoantigens + 1 ) ", "BATF3 Expression", "THBD Expression", "CCL4 Expression", "ATF3 Expression", "Hallmark IFNG Enrichment Score", "Reactome IFNG Enrichment Score", "IFNG Expression")

circ.kras.cms.res <- list()
circ.kras.res <- list()
gsva.es.lst <- list()
clin.es.lst <- list()

## BEGIN ANALYSIS

hallmark.gene.sets <- names(bindeaGsets)[grepl(pattern="HALLMARK", names(bindeaGsets))]
immune.gene.sets <- c("B cells", "T cells", "iDC", "Neutrophils", "Th1 cells", "Th2 cells", "Cytotoxic cells")

gsva.es.lst <- list()
clin.es.lst <- list()

other.genes <- c("STAT1", "CIITA", "CXCL10", "IRF1", "IL6", "IFNG", "CD247")
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

    gsva.es.lst[[dataset]] <- es
    clin.es.lst[[dataset]] <- clin.es
}

neoantigen.signature <- "neoantigens"
neoantigen.label <- "Log10 ( # Neoantigens + 1 )"

kras.cms.res <- list()
kras.cms.interaction.res <- list()
kras.res <- list()
cms.res <- list()

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

    res <- doKrasAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    kras.res[[dataset]] <- res$tests
    kras.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-kras.tsv"), kras.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    res <- doCMSAnalysis(es, clin.es, to.plot=to.plot, num.permutations = 0, calc.fdr = FALSE)
    cms.res[[dataset]] <- res$tests
    cms.tbl <- res$tests
    write.table(file = paste0(dataset, "-immune-hallmark-cms.tsv"), cms.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

}

to.plot <- circ.signature
ylabels <- circ.signature.label

## Begin volcano plot of Bindea vs KRAS (Fig 1a)

tbl.tcga <- kras.res[["tcga"]][immune.gene.sets,]
tbl.tcga$Dataset <- "TCGA"
rownames(tbl.tcga) <- NULL
tbl.kfs <- kras.res[["kfs"]][immune.gene.sets,]
tbl.kfs$Dataset <- "KFSYSCC"
rownames(tbl.kfs) <- NULL
tbl <- rbind(tbl.tcga, tbl.kfs)

tbl$pval <- as.numeric(tbl$pval)
tbl$estimate <- as.numeric(tbl$estimate)
tbl$estimate.lb <- as.numeric(tbl$estimate.lb)
tbl$estimate.ub <- as.numeric(tbl$estimate.ub)
tbl$Dataset <- factor(tbl$Dataset, levels = c("TCGA", "KFSYSCC"))
g.bindea <- ggplot(data = tbl, aes(x = estimate, y = - log10(pval), colour = Dataset))
g.bindea <- g.bindea + geom_point()
## g.bindea <- g.bindea + geom_errorbarh(aes(xmin=estimate.lb, xmax=estimate.ub), height = 0.05)
## Just jitter Th1 KFS and iDC TCGA
j1 <- ( tbl$Dataset == "TCGA" ) & ( tbl$gene.set == "iDC" )
j2 <- ( tbl$Dataset == "KFSYSCC" ) & ( tbl$gene.set == "Th1 cells" )
tbl.neg <- subset(tbl, ( estimate < 0 & !j1 ) | j2)
tbl.pos <- subset(tbl, ( estimate > 0 & !j2 ) | j1)
g.bindea <- g.bindea + geom_text(data = tbl.pos, aes(x = estimate, y = - log10(pval), label = gene.set), hjust = -0.1, show.legend = FALSE)
g.bindea <- g.bindea + geom_text(data = tbl.neg, aes(x = estimate, y = - log10(pval), label = gene.set), hjust = 1.1, show.legend = FALSE)
g.bindea <- g.bindea + xlab("KRAS WT - KRAS MT")
g.bindea <- g.bindea + ylab(bquote(-Log[10]~italic('p')-value))
## Don't show legend, since it will be shown in panel B
## g.bindea <- g.bindea + theme(legend.position = "none")
## lim <- max(max(abs(tbl$estimate.ub)), max(abs(tbl$estimate.lb)))
## lim <- max(abs(tbl$estimate)) + 0.05
lim <- 0.2
g.bindea <- g.bindea + xlim(c(-lim, lim))

## End volcano plot of Bindea vs KRAS

## Begin volcano plot of Hallmark vs KRAS

tbl.tcga <- kras.res[["tcga"]][hallmark.gene.sets,]
tbl.tcga$Dataset <- "TCGA"
rownames(tbl.tcga) <- NULL
tbl.kfs <- kras.res[["kfs"]][hallmark.gene.sets,]
tbl.kfs$Dataset <- "KFSYSCC"
rownames(tbl.kfs) <- NULL
tbl <- rbind(tbl.tcga, tbl.kfs)

tbl$pval <- as.numeric(tbl$pval)
tbl$estimate <- as.numeric(tbl$estimate)
tbl$estimate.lb <- as.numeric(tbl$estimate.lb)
tbl$estimate.ub <- as.numeric(tbl$estimate.ub)
tbl$Dataset <- factor(tbl$Dataset, levels = c("TCGA", "KFSYSCC"))
tbl$gene.set <- gsub(tbl$gene.set, pattern="HALLMARK_", replacement="")

g.hall <- ggplot(data = tbl, aes(x = estimate, y = - log10(pval), colour = Dataset))
g.hall <- g.hall + geom_point()
## Manually move two labels:
## TCGA UV_RESPONSE_UP
## KFS ESTRORGEN_RESPONSE_LATE
j1 <- ( tbl$Dataset == "TCGA" ) & ( tbl$gene.set == "UV_RESPONSE_UP" )
j1 <- j1 | ( tbl$Dataset == "KFSYSCC" ) & ( tbl$gene.set == "INTERFERON_GAMMA_RESPONSE" )
j1 <- j1 | ( tbl$Dataset == "KFSYSCC" ) & ( tbl$gene.set == "INTERFERON_ALPHA_RESPONSE" )
j2 <- ( tbl$Dataset == "KFSYSCC" ) & ( tbl$gene.set == "ESTROGEN_RESPONSE_LATE" )
tbl.neg <- subset(tbl, ( pval < 0.1 & estimate < 0 & !j2 ) | j1)
tbl.pos <- subset(tbl, ( pval < 0.1 & estimate > 0 & !j1 ) | j2)
g.hall <- g.hall + geom_text(data = tbl.neg, aes(x = estimate, y = - log10(pval), label = gene.set), hjust = 1.05, show.legend = FALSE, size = 2.5)
g.hall <- g.hall + geom_text(data = tbl.pos, aes(x = estimate, y = - log10(pval), label = gene.set), hjust = -0.05, show.legend = FALSE, size = 2.5)
g.hall <- g.hall + xlab("KRAS WT - KRAS MT")
g.hall <- g.hall + ylab(bquote(-Log[10]~italic('p')-value))
## lim <- max(abs(tbl$estimate)) + 0.05
lim <- 0.25
g.hall <- g.hall + xlim(c(-lim, lim))

## End volcano plot of Hallmark vs KRAS

pg <- plot_grid(g.bindea, g.hall, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Figure 1", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("fig1.tiff", width=14)

## Begin volcano plot of IFNG, STAT1, CIITA vs KRAS

ifng.genes <- c("STAT1", "CIITA", "CXCL10")

tbl.tcga <- kras.res[["tcga"]][ifng.genes,]
tbl.tcga$Dataset <- "TCGA"
rownames(tbl.tcga) <- NULL
tbl.kfs <- kras.res[["kfs"]][ifng.genes,]
tbl.kfs$Dataset <- "KFSYSCC"
rownames(tbl.kfs) <- NULL
tbl <- rbind(tbl.tcga, tbl.kfs)

tbl$pval <- as.numeric(tbl$pval)
tbl$estimate <- as.numeric(tbl$estimate)
tbl$estimate.lb <- as.numeric(tbl$estimate.lb)
tbl$estimate.ub <- as.numeric(tbl$estimate.ub)
tbl$Dataset <- factor(tbl$Dataset, levels = c("TCGA", "KFSYSCC"))
tbl$gene.set <- gsub(tbl$gene.set, pattern="HALLMARK_", replacement="")

g.infg <- ggplot(data = tbl, aes(x = estimate, y = - log10(pval), colour = Dataset))
g.infg <- g.infg + geom_point()
## Manually move two labels:
## TCGA UV_RESPONSE_UP
## KFS ESTRORGEN_RESPONSE_LATE
g.infg <- g.infg + geom_text(data = tbl, aes(x = estimate, y = - log10(pval), label = gene.set), hjust = -0.1, show.legend = FALSE, size = 5)
g.infg <- g.infg + xlab("KRAS WT - KRAS MT")
g.infg <- g.infg + ylab(bquote(-Log[10]~italic('p')-value))
## lim <- max(abs(tbl$estimate)) + 0.05

g.infg <- g.infg + xlim(c(-0.1,0.75))
if(add.figure.labels) {
g.infg <- g.infg + ggtitle("Supp Figure 6") + theme(plot.title = element_text(hjust = 1))
}
print(g.infg)
ggsave("supp-fig6.tiff")

## End volcano plot of IFNG, STAT1, CIITA vs KRAS



## Plot Supp Figure 1: CIRC vs KRAS
to.plot <- circ.signature
ylabels <- circ.signature.label

dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.tbl <- kras.res[[dataset]]
p.tcga <- plotKrasAnalysis(es, clin.es, kras.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.tbl <- kras.res[[dataset]]
p.kfs <- plotKrasAnalysis(es, clin.es, kras.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

pg <- plot_grid(p.tcga$CIRC, p.kfs$CIRC, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 1", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig1.tiff", width=14)

## Plot the pairwise overlap of the immune signatures and CIRC, restricting to expressed genes.
## Supp Figure 2
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

## Plot CMS analysis for TCGA and KFS (Supp Figure 5)

## TCGA
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
cms.tbl <- cms.res[[dataset]]

to.plot <- circ.signature
ylabels <- circ.signature.label

main <- ifelse(dataset == "tcga", "TCGA", "KFSYSCC")

p.cms.tcga <- plotCMSAnalysis(es, clin.es, cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

## KFS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]

es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
cms.tbl <- cms.res[[dataset]]

main <- ifelse(dataset == "tcga", "TCGA", "KFSYSCC")

p.cms.kfs <- plotCMSAnalysis(es, clin.es, cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = main)

grid.newpage()
pg <- plot_grid(p.cms.tcga$CIRC, p.cms.kfs$CIRC, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Supp Figure 10", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("supp-fig10.tiff", width=14)


## Plot Figure 3 CIRC vs KRAS + CMS and Figure 4 CIRC forest
to.plot <- circ.signature
ylabels <- circ.signature.label

to.plot <- c(to.plot, "STAT1", "CXCL10", "CIITA", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
ylabels <- c(ylabels, "STAT1 Expression", "CXCL10 Expression", "CIITA Expression", "INTERFERON_GAMMA_RESPONSE Enrichment")

## TCGA
dataset <- "tcga"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.cms.tbl <- kras.cms.res[[dataset]]
p.tcga <- plotKrasCMSAnalysis(es, clin.es, kras.cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
p.forest.tcga <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

## KFS
dataset <- "kfs"
es <- gsva.es.lst[[dataset]]
clin.es <- clin.es.lst[[dataset]]
kras.cms.tbl <- kras.cms.res[[dataset]]
p.kfs <- plotKrasCMSAnalysis(es, clin.es, kras.cms.tbl, to.plot=to.plot, ylabels=ylabels, analysis.name = dataset, plot.adjusted.pvalues = FALSE, plot.pvals.as.stars = TRUE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))
p.forest.kfs <- doForestAnalysis(es, clin.es, analysis.name = dataset, to.plot=to.plot, ylabels=ylabels, kras.status.field="kras", kras.states=c("MT","WT"), num.permutations = 0, calc.fdr.with.bh = FALSE, calc.fdr = FALSE, main = ifelse(dataset == "tcga", "TCGA", "KFSYSCC"))

gene <- "CIRC"
pg <- plot_grid(p.tcga[[gene]], p.kfs[[gene]], labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Figure 2", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("fig2.tiff", width=14)

pg <- plot_grid(p.forest.tcga[[gene]]$plot, p.forest.kfs[[gene]]$plot, labels = c("A", "B"), nrow=1)
if(add.figure.labels) {
    pg <- pg + draw_label("Figure 3", x=1, y=1, vjust=1, hjust=1)
}
pg
ggsave("fig3.tiff", width=14)

## Make multivariate cms-kras plots of hallmark and bindea immune (Figures 5 and 6)
ns <- c("immune", "hallmark", "ifng")
files <- c("fig4", "fig5", "supp-fig11")
labels <- c("Figure 4", "Figure 5", "Supp Figure 11")
labels <- c(NA, NA, "Supp Figure 11")
if(!add.figure.labels) {
    labels <- c(NA, NA, NA)
}
plots <- list()

hallmark.levels.to.plot <- c("INTERFERON_GAMMA_RESPONSE", "INFLAMMATORY_RESPONSE", "IL6_JAK_STAT3_SIGNALING", "COMPLEMENT", "INTERFERON_ALPHA_RESPONSE")
immune.levels.to.plot <- c("Cytotoxic cells", "Neutrophils", "Th1 cells")



for(i in 1:length(ns)) {
    name <- ns[i]

    tbl.tcga <- kras.cms.res[["tcga"]]
    tbl.kfs <- kras.cms.res[["kfs"]]    
    rownames(tbl.tcga) <- NULL
    rownames(tbl.kfs) <- NULL    
    
    tbl.tcga$gene.set <- gsub(tbl.tcga$gene.set, pattern="HALLMARK_", replacement="")
    tbl.kfs$gene.set <- gsub(tbl.kfs$gene.set, pattern="HALLMARK_", replacement="")    
    
    ## Restrict to subset of pathways/populations that came out of our analysis in Fig 1 (KRAS)
    gene.set.order <- NULL
    if(name == "immune") {
        gene.set.order <- immune.levels.to.plot
    }
    if(name == "hallmark") {
        gene.set.order <- hallmark.levels.to.plot
    }
    if(name == "ifng") {
        gene.set.order <- ifng.genes
    }

    tbl.tcga <- subset(tbl.tcga, gene.set %in% gene.set.order)
    tbl.kfs <- subset(tbl.kfs, gene.set %in% gene.set.order)    
    
    tbl.kfs <- tbl.kfs  %>% gather(key, value, -gene.set) %>% extract(key, c("cms", "kras", "interval"), "(CMS.)(.+).vs.CMS2MT\\.(.+)") %>% spread(interval, value)
    tbl.tcga <- tbl.tcga  %>% gather(key, value, -gene.set) %>% extract(key, c("cms", "kras", "interval"), "(CMS.)(.+).vs.CMS2MT\\.(.+)") %>% spread(interval, value)     
    
    tbl.kfs$gene.set <- factor(tbl.kfs$gene.set, levels = gene.set.order)
    tbl.tcga$gene.set <- factor(tbl.tcga$gene.set, levels = gene.set.order)
    tbl.kfs$cms <- factor(tbl.kfs$cms, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
    tbl.tcga$cms <- factor(tbl.tcga$cms, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
    tbl.kfs$kras <- factor(tbl.kfs$kras, levels = c("MT", "WT"))
    tbl.tcga$kras <- factor(tbl.tcga$kras, levels = c("MT", "WT"))
    
    tbl1.mle.col <- "estimate"
    tbl1.lb.col <- "estimate.lb"
    tbl1.ub.col <- "estimate.ub"
    tbl2.mle.col <- "estimate"
    tbl2.lb.col <- "estimate.lb"
    tbl2.ub.col <- "estimate.ub"
    main1 <- "TCGA"
    main2 <- "KFSYSCC"

    text.size <- NULL
    if(name == "hallmark") { text.size <- 7 }

    l <- plot.stacked.faceted.effect.size.figures(tbl.tcga, tbl.kfs, tbl1.mle.col, tbl1.lb.col, tbl1.up.col, tbl2.mle.col, tbl2.lb.col, tbl2.ub.col, main1, main2, text.size = text.size)
    ##    do.call("grid.arrange", l)
    plots[[name]] <- l
    if(name == "ifng") {
        pg <- plot_grid(l[[1]], l[[3]], l[[2]], labels = c("A", "", "B"), nrow=2, rel_widths = c(10,1.5), align="none")
        label <- labels[i]
        if(!is.na(label)) {
            pg <- pg + draw_label(label, x=1, y=1, vjust=1, hjust=1)
        }
        pg
        filename <- paste0(name, "-kras-effect.tiff")
        filename <- paste0(files[i], ".tiff")
        ggsave(filename)
    }
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

pg <- plot_grid(plots[["immune"]][[1]], plots[["hallmark"]][[1]], plots[["hallmark"]][[3]], plots[["immune"]][[2]], plots[["hallmark"]][[2]], labels = c("A", "C", "", "B", "D"), nrow=2, rel_widths = c(10,10,2.5), align="none")
label <- "Figure 4"
if(!is.null(label)) {
    pg <- pg + draw_label(label, x=1, y=1, vjust=1, hjust=1)
}
pg
filename <- paste0("immune-hallmark", "-kras-effect.tiff")
filename <- paste0("fig4.tiff")
ggsave(filename)

tbl <- table(unlist(lapply(paste0("HALLMARK_", hallmark.levels.to.plot), function(x) bindeaGsets[[as.character(x)]])))
tbl <- tbl[tbl > 1]
tbl <- as.data.frame(tbl)
colnames(tbl) <- c("gene", "freq")
tbl <- tbl[order(tbl$freq, decreasing=TRUE),]
write.table(file="genes-enriched-in-hallmark-immune-pathways.tsv", tbl, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
cat(paste0("Genes enriched in Hallmark immune pathways: ", paste(hallmark.levels.to.plot, collapse=", "), "\n"))
write.table(tbl, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## Convert TIFFs to PDFs
## for x in `ls *.tiff`; do echo $x; sips -s format pdf $x --out `basename $x .tiff`.pdf ; done
## sips -s format pdf --out Supp\ fig\ 4.pdf Supp\ fig\ 4.tif
## sips -s format pdf --out Supp\ fig\ 5.pdf Supp\ fig\ 5.tif
## sips -s format pdf --out Supp\ fig\ 7.pdf Supp\ fig\ 7.tif
## sips -s format pdf --out Supp\ fig\ 8.pdf Supp\ fig\ 8.tif
## sips -s format pdf --out Supp\ fig\ 9.pdf Supp\ fig\ 9.tif 

## cp Supp\ fig\ 4.pdf supp-fig4.pdf
## cp Supp\ fig\ 5.pdf supp-fig5.pdf
## cp Supp\ fig\ 7.pdf supp-fig7.pdf
## cp Supp\ fig\ 8.pdf supp-fig8.pdf
## cp Supp\ fig\ 9.pdf supp-fig9.pdf

## "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" --output crc-cms-kras-figures.pdf supp-fig1.pdf supp-fig2.pdf supp-fig3a.pdf supp-fig3b.pdf fig1.pdf Supp\ fig\ 4.pdf Supp\ fig\ 5.pdf supp-fig6.pdf Supp\ fig\ 7.pdf Supp\ fig\ 8.pdf Supp\ fig\ 9.pdf supp-fig10.pdf fig2.pdf fig3.pdf fig4.pdf supp-fig11.pdf

save.image(".Rdata")
cat("Successfully completed.\n")
