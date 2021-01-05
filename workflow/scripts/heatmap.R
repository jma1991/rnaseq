# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

selectCPM <- function(x) {

    id1 <- which(colnames(x) == "falsePos") + 1

    id2 <- ncol(x)

    idx <- seq(from = id1, to = id2)

    cpm <- x[, idx]

    cpm <- as.matrix(cpm)

    return(cpm)

}

pheatmap.mat <- function(x) {

    # Scale rows by 'variance-aware' Z-transformation

    M <- rowMeans(x, na.rm = TRUE)

    DF <- ncol(x) - 1

    isNA <- is.na(x)

    anyNA <- any(isNA)

    if (anyNA) {

        mode(isNA) <- "integer"

        DF <-  DF - rowSums(isNA)

        DF[DF == 0] <- 1

    }

    x <- x - M

    V <- rowSums(x^2, na.rm = TRUE) / DF

    x <- x / sqrt(V + 0.01)

}

pheatmap.color <- function(x) {

    # Return color vector

    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = x)))(100)

}

pheatmap.breaks <- function(x) {

    # Return breaks vector

    abs <- max(abs(x))

    abs <- min(abs, 5)

    seq(-abs, +abs, length.out = 101)

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

    res <- read.csv(input$csv, row.names = 1)

    ind <- order(res$PValue)[seq_len(50)]

    res <- res[ind, ]

    cpm <- selectCPM(res)

    mat <- pheatmap.mat(cpm)
    
    pheatmap(
        mat = mat, 
        color = pheatmap.color("RdBu"), 
        breaks = pheatmap.breaks(mat), 
        cluster_rows = pheatmap.cluster_rows(mat), 
        cluster_cols = pheatmap.cluster_cols(cpm),
        labels_row = res$geneName,
        filename = output$pdf, 
        width = 6, 
        height = 8
    )

    # Image function

    library(magick)

    pdf <- image_read_pdf(output$pdf)
    
    pdf <- image_trim(pdf)

    pdf <- image_border(pdf, color = "#FFFFFF", geometry = "50x50")
    
    pdf <- image_write(pdf, path = output$pdf, format = "pdf")

}

main(snakemake@input, snakemake@output, snakemake@log)
