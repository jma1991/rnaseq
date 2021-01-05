# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, log, wildcards) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function
    
    library(limma)
    
    obj <- readRDS(input$rds)

    ann <- read.delim(input$tsv, header = FALSE, col.names = c("transcript_id", "gene_id", "gene_name"))

    fit <- lmFit(obj)

    str <- paste(wildcards$A, wildcards$B, sep = "-")

    lvl <- colnames(fit$design)

    con <- makeContrasts(contrasts = str, levels = lvl)
    
    fit <- contrasts.fit(fit, contrasts = con)
    
    fit <- eBayes(fit)
    
    res <- topTable(fit, number = Inf, sort.by = "P")
    
    cts <- 2 ^ obj$E[rownames(res), ]
    
    idx <- list(
        A = which(fit$design[, wildcards$A] == 1),
        B = which(fit$design[, wildcards$B] == 1)
    )
    
    dat <- data.frame(
        geneId         = ann$gene_id[match(rownames(res), ann$gene_id)],
        geneName       = ann$gene_name[match(rownames(res), ann$gene_id)],
        baseMean       = rowMeans(cts[, c(idx$A, idx$B)]),
        baseMeanA      = rowMeans(cts[, idx$A]),
        baseMeanB      = rowMeans(cts[, idx$B]),
        foldChange     = 2 ^ res$logFC,
        log2FoldChange = res$logFC,
        PValue         = res$P.Value,
        PAdj           = p.adjust(res$P.Value, method = "hochberg"),
        FDR            = res$adj.P.Val,
        falsePos       = round(seq_along(res$adj.P.Val) * res$adj.P.Val)
    )
    
    out <- merge(dat, cts[, c(idx$A, idx$B)], by.x = "geneId", by.y = "row.names", sort = FALSE)
    
    write.csv(out, file = output$csv, quote = FALSE, row.names = FALSE)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@wildcards)