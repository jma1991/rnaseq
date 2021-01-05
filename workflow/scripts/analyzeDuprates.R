# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, params, threads) {

    library(dupRadar)

    mat <- analyzeDuprates(input$bam, input$gtf, stranded = params$stranded, paired = TRUE, threads = threads)

    write.csv(mat, file = output$csv, quote = FALSE, row.names = FALSE)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@threads)