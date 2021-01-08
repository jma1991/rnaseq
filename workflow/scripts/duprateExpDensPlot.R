# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(dupRadar)

    mat <- read.csv(input$csv)

    pdf(output$pdf)

    duprateExpDensPlot(DupMat = mat)

    dev.off()

}

main(snakemake@input, snakemake@output, snakemake@log)
