#!/usr/bin/env Rscript

goseq.genome <- function(x) {

    # Return genome identifier

    x <- GenomeInfoDb::mapGenomeBuilds(x, style = "UCSC")

    x <- unique(x$ucscID)

    return(x)

}

goseq.id <- function(x) {

    # Return gene identifier

    d <- c("UCSC" = "knownGene", "Ensembl" = "ensGene", "Vega" = "vegaGene", "geneSymbol" = "refGene")

    v <- d[x]

    return(v)

}

main <- function(input, output, log, config) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function    
    
    library(goseq)

    res <- read.csv(input$csv, stringsAsFactors = FALSE)

    vec <- as.integer(res$FDR < 0.05)

    names(vec) <- res$geneId

    genome <- goseq.genome(config$txome$genome)

    id <- goseq.id(config$txome$source)

    pwf <- nullp(vec, genome, id)

    out <- goseq(pwf, genome, id)

    write.csv(out, file = output$csv)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@config)
