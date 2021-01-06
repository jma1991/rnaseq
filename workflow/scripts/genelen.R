#!/usr/bin/env Rscript

main <- function(input, output, log) {

	library(GenomicFeatures)

	txi <- makeTxDbFromGFF(input$gtf, format = "gtf")

	rng <- transcriptsBy(x, "gene")

	len <- median(width(rng))

	save(len, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
