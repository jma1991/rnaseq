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

    library(edgeR)

    library(limma)

    dge <- readRDS(input$rds)
    
    mod <- model.matrix(~ 0 + condition, data = dge$samples)
    
    colnames(mod) <- unique(dge$samples$condition)
    
    els <- voom(dge, design = mod)

    saveRDS(els, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)
