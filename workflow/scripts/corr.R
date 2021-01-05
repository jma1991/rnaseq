# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

pheatmap.mat <- function(x) {

    # Compute spearman correlation

    as.matrix(cor(x, method = "spearman"))

}

pheatmap.color <- function(x) {

    # Return color vector

    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = x)))(100)

}

pheatmap.breaks <- function(x) {

    # Return breaks vector

    seq(-1, +1, length.out = 101)

}

pheatmap.cluster_rows <- function(x) {

    # Return hclust object for rows

    hclust(dist(x, method = "euclidean"), method = "complete")

}

pheatmap.cluster_cols <- function(x) {

    # Return hclust object for columns

    hclust(dist(t(x), method = "euclidean"), method = "complete")

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

    rho <- pheatmap.mat(dat)

    pheatmap(
        mat = rho, 
        color = pheatmap.color("RdBu"), 
        breaks = pheatmap.breaks(rho), 
        cluster_rows = pheatmap.cluster_rows(rho), 
        cluster_cols = pheatmap.cluster_cols(rho), 
        filename = output$pdf, 
        width = 7, 
        height = 7
    )

    # Image function

    library(magick)

    pdf <- image_read_pdf(output$pdf)
    
    pdf <- image_trim(pdf)

    pdf <- image_border(pdf, color = "#FFFFFF", geometry = "50x50")
    
    pdf <- image_write(pdf, path = output$pdf, format = "pdf")

}

main(snakemake@input, snakemake@output, snakemake@log)
