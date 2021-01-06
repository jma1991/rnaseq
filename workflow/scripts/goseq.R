#!/usr/bin/env Rscript

main <- function(input, output, params, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function    
    
    library(goseq)

    res <- read.csv(input$csv, stringsAsFactors = FALSE)

    len <- readRDS(input$rds)

    vec <- as.integer(res$FDR < 0.05)

    names(vec) <- res$geneId

    map <- mapGenomeBuilds(params$genome, style = "UCSC")
    
    pwf <- nullp(vec, map$ucscID, "ensGene")

    out <- goseq(pwf, map$ucscID, "ensGene")

    write.csv(out, file = output$csv)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)
