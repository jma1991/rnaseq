# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

add.status <- function(x, fdr.threshold = 0.05) {

    # Determine the status of each gene in the results table

    x$status <- factor("NS", levels = c("Up", "NS", "Down"))

    x$status[x$log2FoldChange > 0 & x$FDR < fdr.threshold] <- "Up"

    x$status[x$log2FoldChange < 0 & x$FDR < fdr.threshold] <- "Down"

    return(x)

}

add.symbol <- function(x, n = 50) {

    # Show symbol for the top N genes in the results table

    i <- sort(x$FDR, decreasing = FALSE, index.return = TRUE)$ix[seq_len(n)]

    x$symbol <- ""
    
    x$symbol[i] <- x$geneName[i]

    return(x)

}

main <- function(input, output, log ) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function
    
    library(ggplot2)

    library(ggrepel)

    library(scales)
    
    res <- read.csv(input$csv, row.names = 1)

    res <- add.status(res, fdr.threshold = 0.05)

    res <- add.symbol(res, n = 50)

    res <- res[order(res$FDR, decreasing = TRUE), ]
    
    col <- c("Up" = "#d8b365", "NS" = "#f5f5f5", "Down" = "#5ab4ac")

    lab <- c(
        "Up" = sprintf("Up (%s)", comma(sum(res$status == "Up"))),
        "NS" = sprintf("NS (%s)", comma(sum(res$status == "NS"))),
        "Down" = sprintf("Down (%s)", comma(sum(res$status == "Down")))
    )
    
    plt <- ggplot(res, aes(log10(baseMean), log2FoldChange, colour = status, label = symbol)) + 
        geom_point() +
        scale_colour_manual(values = col, labels = lab) +
        geom_text_repel(size = 1.88, colour = "#000000", show.legend = FALSE, max.overlaps = Inf) +  
        labs(x = "Normalized mean (log10)", y = "Fold change (log2)", colour = "Status") + 
        theme_bw() + 
        theme(
            aspect.ratio = 1,
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
