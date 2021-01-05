# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output) {

    library(dupRadar)

    mat <- read.csv(input$csv, stringsAsFactors = FALSE)

    pdf(output$pdf)

    expressionHist(DupMat = mat)

    dev.off()

}

main(snakemake@input, snakemake@output)