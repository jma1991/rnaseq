# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output) {

    library(ggplot2)

    c_curve <- read.delim(input$txt[1])

    lc_extrap <- read.delim(input$txt[2])

    # Scale counts per million

    c_curve <- c_curve / 1e6

    lc_extrap <- lc_extrap / 1e6

    # Plot library complexity
    
    plt <- ggplot(lc_extrap, aes(x = TOTAL_READS, y = EXPECTED_DISTINCT)) + 
        geom_ribbon(aes(ymin = LOWER_0.95CI, ymax = UPPER_0.95CI), colour = "#4E79A7", fill = "#A0CBE8") + 
        geom_line(colour = "#E15759", linetype = "dashed") + 
        geom_line(data = c_curve, aes(total_reads, distinct_reads), colour = "#000000") + 
        labs(x = "Sequenced reads (M)", y = "Distinct reads (M)") + 
        theme_classic() + 
        theme(aspect.ratio = 1,
              axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "lines")),
              axis.title.y = element_text(margin = unit(c(0, 1, 0, 0), "lines")))
    
    ggsave(output$pdf, plot = plt, width = 3, height = 3)

}

main(snakemake@input, snakemake@output)