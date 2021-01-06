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

    normMat <- object$length

    normMat <- normMat/exp(rowMeans(log(normMat)))

    normCts <- cts/normMat

    eff.lib <- edgeR::calcNormFactors(normCts) * colSums(normCts)

    normMat <- sweep(normMat, 2, eff.lib, "*")

    normMat <- log(normMat)

    y <- edgeR::DGEList(cts, samples = samples, group = group)

    edgeR::scaleOffset(y, normMat)

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

    dat <- read.csv("config/samples.csv", row.names = "sample")

    dge <- DGEList(object = obj, samples = dat, group = dat$condition)

    ind <- filterByExpr(dge, group = dat$condition)

    dge <- dge[ind, , keep.lib.sizes = FALSE]

    dge <- calcNormFactors(dge)

    dge <- estimateDisp(dge)

    saveRDS(dge, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
