# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

logcounts <- function(object) {
    # Return log-transformed counts
    UseMethod("logcounts")
}

logcounts.DESeqDataSet <- function(object) {
    # DESeq2 method
    rld <- DESeq2::rlog(object, blind = FALSE)
    return(SummarizedExperiment::assay(rld))
}

logcounts.DGEList <- function(object) {
    # edgeR method
    return(edgeR::cpm(object, log = TRUE))
}

logcounts.EList <- function(object) {
    # limma method
    return(object$E)
}

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    obj <- readRDS(input$rds)

    mat <- logcounts(object = obj)

    write.csv(mat, file = output$csv)

}

main(snakemake@input, snakemake@output, snakemake@log)
