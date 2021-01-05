# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

counts <- function(object) {
    # Return raw counts
    UseMethod("counts")
}

counts.DESeqDataSet <- function(object) {
    # DESeq2 method
    return(DESeq2::counts(object))
}

counts.DGEList <- function(object) {
    # edgeR method
    return(object$counts)

}

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    obj <- readRDS(input$rds)

    mat <- counts(object = obj)

    write.csv(mat, file = output$csv)

}

main(snakemake@input, snakemake@output, snakemake@log)
