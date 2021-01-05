# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, params) {

    library(ggplot2)

    dat <- lapply(input$txt, read.delim)

    dat <- do.call(rbind, dat)

    dat <- dat / 1e6

    dat$id <- params$ids
    
    plt <- ggplot(dat, aes(id, median_estimated_unobs)) + 
        geom_errorbar(width = 0.1, aes(ymin = lower_ci, ymax = upper_ci), colour = "#E15759") + 
        geom_point(shape = 21, size = 3, fill = "#BAB0AC") + 
        labs(x = "Sample", y = "Species Richness (M)") + 
        coord_flip() +
        theme_classic() + 
        theme(aspect.ratio = 1,
              axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "lines")),
              axis.title.y = element_text(margin = unit(c(0, 1, 0, 0), "lines")))
    
    ggsave(output$pdf, plot = plt, width = 4, height = 3)

}

main(snakemake@input, snakemake@output, snakemake@params)