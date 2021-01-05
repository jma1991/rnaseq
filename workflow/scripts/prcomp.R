# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

calculatePCA <- function(x, n = 500) {

    var <- apply(x, 1, var)

    num <- min(length(var), n)
    
    sel <- order(var, decreasing = TRUE)[seq_len(num)]

    sub <- x[sel, ]
    
    pca <- prcomp(t(sub))
    
    dev <- pca$sdev^2 / sum(pca$sdev^2)

    dat <- data.frame("PCA.1" = pca$x[, 1], "PCA.2" = pca$x[, 2], row.names = colnames(sub))

    attr(dat, "percentVar") <- dev[1:2]

    return(dat)

}

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(ggplot2)

    dat <- read.csv("config/sample_table.csv", stringsAsFactors = FALSE)

    mat <- read.csv(input$csv, row.names = 1)

    pca <- calculatePCA(mat, n = 500)

    dat <- merge(dat, pca, by.x = "sample_name", by.y = "row.names")

    var <- attr(pca, "percentVar")

    plt <- ggplot(dat, aes(PCA.1, PCA.2, color = condition)) + 
        geom_point() + 
        xlab(paste0("PC1: ", round(var[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(var[2] * 100), "% variance")) +
        coord_fixed() + 
        theme_bw() + 
        theme(
            axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "lines")),
            axis.title.y = element_text(margin = unit(c(0, 1, 0, 0), "lines"))
        )

    ggsave(output$pdf, plot = plt, width = 7, height = 7)

    # Image function

    library(magick)

    pdf <- image_read_pdf(output$pdf)
    
    pdf <- image_trim(pdf)

    pdf <- image_border(pdf, color = "#FFFFFF", geometry = "50x50")
    
    pdf <- image_write(pdf, path = output$pdf, format = "pdf")

}

main(snakemake@input, snakemake@output, snakemake@log)
