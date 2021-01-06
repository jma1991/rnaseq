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
    
    library(DESeq2)
    
    dds <- readRDS(input$rds)

    ann <- read.delim(input$tsv, header = FALSE, col.names = c("gene_id", "gene_name"))
 
    con <- c("condition", wildcards$A, wildcards$B)
    
    res <- results(dds, contrast = con)
    
    res <- res[order(res$pvalue), ]
    
    cpm <- counts(dds, normalized = TRUE)[rownames(res), ]
    
    idx <- list(
        A = which(dds$condition %in% wildcards$A),
        B = which(dds$condition %in% wildcards$B)
    )

    dat <- data.frame(
        geneId         = ann$gene_id[match(rownames(res), ann$gene_id)],
        geneName       = ann$gene_name[match(rownames(res), ann$gene_id)],
        baseMean       = rowMeans(cpm[, c(idx$A, idx$B)]),
        baseMeanA      = rowMeans(cpm[, idx$A]),
        baseMeanB      = rowMeans(cpm[, idx$B]),
        foldChange     = 2 ^ res$log2FoldChange,
        log2FoldChange = res$log2FoldChange,
        PValue         = res$pvalue,
        PAdj           = p.adjust(res$pvalue, method = "hochberg"),
        FDR            = res$padj,
        falsePos       = round(seq_along(res$padj) * res$padj)
    )
    
    out <- merge(dat, cpm[, c(idx$A, idx$B)], by.x = "geneId", by.y = "row.names", sort = FALSE)
    
    write.csv(out, file = output$csv, quote = FALSE, row.names = FALSE)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@wildcards)
