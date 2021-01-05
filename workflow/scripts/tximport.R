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

    library(tximport)
    
    dat <- read.csv("config/sample_table.csv", colClasses = "character", stringsAsFactors = FALSE)

    hdf <- file.path("results/kallisto/quant", dat$sample_name, "abundance.h5")

    names(hdf) <- dat$sample_name

    tbl <- read.table(input$tsv, col.names = c("transcript_id", "gene_id", "gene_name"))

    txi <- tximport(files = hdf, type = "kallisto", tx2gene = tbl)

    class(txi) <- c("list", "tximport")

    saveRDS(txi, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
