# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

calculateMDS <- function(x) {

    # Perform multi-dimensional scaling

    d <- dist(t(x))
    
    d <- cmdscale(d)

    d <- as.data.frame(d)

    colnames(d) <- c("MDS.1", "MDS.2")

    return(d)

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

    mds <- calculateMDS(mat)

    dat <- merge(dat, mds, by.x = "sample_name", by.y = "row.names")

    plt <- ggplot(dat, aes(MDS.1, MDS.2, colour = condition)) + 
        geom_point() + 
        labs(x = "MDS 1", y = "MDS 2") + 
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
