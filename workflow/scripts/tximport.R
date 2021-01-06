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
    
    dat <- read.csv("config/samples.csv")

    hdf <- file.path("results/kallisto/quant", dat$sample, "abundance.h5")

    names(hdf) <- dat$sample

    ann <- read.table(input$tsv)

    txi <- tximport(files = hdf, type = "kallisto", tx2gene = ann)

    class(txi) <- c("list", "tximport")

    saveRDS(txi, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
