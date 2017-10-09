#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("annotate"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

option_list <- list(
    make_option(c("--expression-data-input"), action="store",
                default=NULL,
                help=".rdata input file holding expression"),
    make_option(c("--expression-data-output"), action="store",
                default=NULL,
                help="TSV output file of expression data with sample names harmonized to clinical annotations")
)


required.args <- c("expression-data-input", "expression-data-output")

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

expr.rdata <- opt$`expression-data-input`
expr.data.output <- opt$`expression-data-output`

# file <- arguments$args[1]
# expr.rdata <- arguments$args[2]

## Load in the expression matrix.  We expect this to be the only object
# ni the .rdata file
expr.mat.name <- load(expr.rdata)
expr.mat <- get(expr.mat.name)
rm(expr.mat.name)

colnames(expr.mat) <- getSYMBOL(colnames(expr.mat), data='org.Hs.eg')

write.table(expr.mat, file=expr.data.output, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

