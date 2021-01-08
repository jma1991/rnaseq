# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

featureCounts.strandSpecific <- function(x) {

    # Return strand 

    d <- list("U" = 0, "F" = 1, "R" = 2)
    
    v <- d[x]

    return(v)

}

main <- function(input, output, log, threads) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(Rsubread)

    dat <- read.csv("config/samples.csv", colClasses = "character", stringsAsFactors = FALSE)

    out <- featureCounts(
        files = file.path("results/star/align", dat$sample, "Aligned.sortedByCoord.out.bam"),
        annot.ext = input$gtf,
        isGTFAnnotationFile = TRUE,
        strandSpecific = featureCounts.strandSpecific(dat$stranded),
        isPairedEnd = TRUE,
        nthreads = threads
    )

    colnames(out$counts) <- dat$sample

    class(out) <- c("list", "rsubread")

    saveRDS(out, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@threads)
