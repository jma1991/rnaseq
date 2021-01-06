# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

DESeqDataSet <- function(object, ...) {
    # Dispatch method
    UseMethod("DESeqDataSet")
}

DESeqDataSet.tximport <- function(object, colData, design) {
    # Create DESeqDataSet object from tximport
    DESeqDataSetFromTximport(txi = object, colData = colData, design = design)
}

DESeqDataSet.rsubread <- function(object, colData, design) {
    # Create DESeqDataSet object from rsubread
    DESeqDataSetFromMatrix(countData = object$counts, colData = colData, design = design)   
}

main <- function(input, output, log, config) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(DESeq2)

    obj <- readRDS(input$rds)

    dat <- read.csv("config/samples.csv", row.names = "sample", stringsAsFactors = FALSE)

    dds <- DESeqDataSet(object = obj, colData = dat, design = ~ condition)

    dds <- DESeq(dds)

    saveRDS(dds, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@config)
