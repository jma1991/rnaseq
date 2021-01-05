#!/usr/bin/env Rscript

mean.baseMean <- function(x) {
    matrixStats::rowMeans(sapply(x, "[[", "baseMean"))
}

mean.baseMeanA <- function(x) {
    matrixStats::rowMeans(sapply(x, "[[", "baseMeanA"))
}

mean.baseMeanB <- function(x) {
    matrixStats::rowMeans(sapply(x, "[[", "baseMeanB"))
}

mean.foldChange <- function(x) {
    matrixStats::rowMeans(sapply(x, "[[", "foldChange"))
}

mean.log2FoldChange <- function(x) {
    matrixStats::rowMeans(sapply(x, "[[", "log2FoldChange"))
}

sd.log2FoldChange <- function(x) {
    matrixStats::rowSds(sapply(x, "[[", "log2FoldChange"))
}

consensusDE.edger_adj_p <- function(x) {
    x[["edger"]]["PAdj"]
}

consensusDE.deseq2_adj_p <- function(x) {
    x[["deseq2"]]["PAdj"]
}

consensusDE.voom_adj_p <- function(x) {
    x[["limma"]]["PAdj"]
}

consensusDE.edger_rank <- function(x) {
    rank(x[["edger"]]["PValue"])
}

consensusDE.deseq2_rank <- function(x) {
    rank(x[["deseq2"]]["PValue"])
}

consensusDE.limma_rank <- function(x) {
    rank(x[["limma"]]["PValue"])
}

consensusDE.rank_sum <- function(x) {
    matrixStats::rowSums(apply(sapply(x, "[[", "PValue"), 2, rank))
}

consensusDE.p_intersect <- function(x) {
    matrixStats::rowMaxs(sapply(x, "[[", "PValue"))
}

consensusDE.p_union <- function(x) {
    matrixStats::rowMins(sapply(x, "[[", "PValue"))
}


main <- function(input, output) {

    res <- lapply(input$csv, read.csv)

    names(res) <- c("deseq2", "edger", "limma")

    ids <- Reduce(intersect, lapply(res, rownames))

    res <- lapply(res, "[[", i = ids, )
    
    dat <- data.frame(
        geneId           = res[[1]][["geneId"]],
        geneName         = res[[1]][["geneName"]],
        baseMean         = mean.baseMean(res),
        baseMeanA        = mean.baseMeanA(res),
        baseMeanB        = mean.baseMeanB(res),
        foldChange       = mean.foldChange(res),
        log2FoldChange   = mean.log2FoldChange(res),
        log2FoldChangeSD = sd.log2FoldChange(res),
    )

    write.csv(dat, file = output$csv)

}

main(snakemake@input, snakemake@output)
