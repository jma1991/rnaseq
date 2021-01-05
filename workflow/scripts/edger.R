# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

DGEList <- function(object, ...) {
    UseMethod("DGEList")
}

DGEList.rsubread <- function(object, samples, group) {
    edgeR::DGEList(counts = object$counts, samples = samples, group = group)
}

DGEList.tximport <- function(object, samples, group) {

    cts <- object$counts
    
    len <- object$length
    
    len <- len / exp(rowMeans(log(len)))
    
    lib <- log(edgeR::calcNormFactors(cts/len)) + log(colSums(cts/len))
    
    dge <- edgeR::DGEList(cts, samples = samples, group = group)
    
    edgeR::scaleOffset(dge, t(t(log(len)) + lib))

}

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(edgeR)

    obj <- readRDS(input$rds)

    dat <- read.csv("config/sample_table.csv", row.names = "sample_name", stringsAsFactors = FALSE)

    dge <- DGEList(object = obj, samples = dat, group = dat$condition)

    idx <- filterByExpr(dge, group = dat$condition)

    dge <- dge[idx, , keep.lib.sizes = FALSE]

    dge <- calcNormFactors(dge)

    dge <- estimateDisp(dge)

    saveRDS(dge, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
