#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("annotate"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

option_list <- list(
    make_option(c("--data-input"), action="store",
                default=NULL,
                help=".rdata input file holding expression"),
    make_option(c("--data-output"), action="store",
                default=NULL,
                help="TSV output file of expression data with sample names harmonized to clinical annotations")
)


required.args <- c("data-input", "data-output")

descr <- "\
   Translate the entrez ID colnames in the .rdata file to gene symbols."

parser <- OptionParser(usage = "%prog [options]", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if ( length(arguments$args) != 0 ) {
  print_help(parser)
  q(status=1)
}

missing.args <- required.args[!(required.args %in% names(opt))]
if( length(missing.args) > 0 ) {
    cat("Need to specify arguments:\n")
    sapply(missing.args, function(x) cat(paste("--", x, "\n", sep="")))
    q(status=-1)
}

expr.tsv <- opt$`data-input`
expr.data.output <- opt$`data-output`

## Load in the expression matrix.  We expect this to be the only object
# ni the .rdata file
expr.mat <- read.table(expr.tsv, sep="\t", header=TRUE, as.is=TRUE)

rnames <- rownames(expr.mat)

flag <- grepl(pattern="GSM", x=rnames)
rnames[flag] <- sapply(rnames[flag], function(x) {
    r <- regexpr(pattern="GSM\\d+", text=x)
    if(r==-1) {
        cat(paste("Found no match for GSMd+ in: ", x, "\n", sep=""))
        q(status=-1)
    }
    substr(x, r[1], r[1] + attr(r,"match.length")[1] - 1)
})

flag <- grepl(pattern=".CEL", x=rnames)
rnames[flag] <- gsub(rnames[flag], replacement="", pattern=".CEL")

flag <- grepl(pattern="kfsyscc-", x=rnames)
rnames[flag] <- gsub(rnames[flag], replacement="", pattern="kfsyscc-")

flag <- grepl(pattern="french", x=rnames)
rnames[flag] <- gsub(rnames[flag], replacement="", pattern="french-")

flag <- grepl(pattern="petacc", x=rnames)
rnames[flag] <- gsub(rnames[flag], replacement="", pattern="petacc3-")

flag <- grepl(pattern="tcga", x=rnames)
rnames[flag] <- gsub(rnames[flag], replacement="", pattern="tcgacrc_merged-")

if(any(duplicated(rnames))) {
    cat("Duplicate rownames!\n")
    flag <- duplicated(rnames, fromLast=TRUE) | duplicated(rnames, fromLast=FALSE)
    print(rownames(expr.mat)[flag])
    q(status=-1)
}

rownames(expr.mat) <- rnames

write.table(expr.mat, file=expr.data.output, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##more cms_labels_public_all.txt  | cut -f2 | sort | uniq
##dataset
##gse13067
##gse13294
##gse14333
##gse17536
##gse20916
##gse2109
##gse23878
##gse33113
##gse35896
##gse37892
##gse39582
