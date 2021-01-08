# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, params, log, threads) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(dupRadar)

    mat <- analyzeDuprates(
        bam = input$bam,
        gtf = input$gtf,
        stranded = params$stranded,
        paired = params$paired,
        threads = threads
    )

    write.csv(mat, file = output$csv)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log, snakemake@threads)
