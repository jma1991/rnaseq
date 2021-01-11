# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

pheatmap.mat <- function(x) {

    # Return distance matrix

    as.matrix(dist(t(x), method = "euclidean"))

}

pheatmap.color <- function(x) {

    # Return color vector

    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = x)))(100)

}

pheatmap.breaks <- function(x) {

    # Return breaks vector

    seq(0, max(x), length.out = 101)

}

pheatmap.cluster_rows <- function(x) {

    # Return hclust object for rows

    hclust(as.dist(x), method = "complete")

}

pheatmap.cluster_cols <- function(x) {

    # Return hclust object for columns

    hclust(as.dist(x), method = "complete")

}

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(pheatmap)

    library(RColorBrewer)

    dat <- read.csv(input$csv, row.names = 1)

    mat <- pheatmap.mat(dat)

    pheatmap(
        mat = mat,
        color = pheatmap.color("Blues"),
        breaks = pheatmap.breaks(mat),
        cluster_rows = pheatmap.cluster_rows(mat),
        cluster_cols = pheatmap.cluster_cols(mat),
        filename = output$pdf,
        width = 7, 
        height = 7
    )

}

main(snakemake@input, snakemake@output, snakemake@log)
