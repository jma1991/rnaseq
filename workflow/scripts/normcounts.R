# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

normcounts <- function(object) {
    # Return normalized counts
    UseMethod("normcounts")
}

normcounts.DESeqDataSet <- function(object) {
    # DESeq2 method
    return(DESeq2::counts(object, normalized = TRUE))
}

normcounts.DGEList <- function(object) {
    # edgeR method
    return(edgeR::cpm(object))
}

normcounts.EList <- function(object) {
    # limma method
    return(2 ^ object$E)
}

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    obj <- readRDS(input$rds)

    mat <- normcounts(object = obj)

    write.csv(mat, file = output$csv)

}

main(snakemake@input, snakemake@output, snakemake@log)
