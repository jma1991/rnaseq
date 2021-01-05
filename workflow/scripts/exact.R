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

    library(edgeR)

    dge <- readRDS(input$rds)

    ann <- read.delim(input$tsv, header = FALSE, col.names = c("transcript_id", "gene_id", "gene_name"))

    con <- c(wildcards$B, wildcards$A)

    res <- exactTest(dge, pair = con)

    res <- topTags(res, n = Inf, sort.by = "PValue")
    
    cts <- cpm(dge)[rownames(res), ]
    
    idx <- list(
        "A" = which(dge$samples$group %in% wildcards$A),
        "B" = which(dge$samples$group %in% wildcards$B)
    )

    dat <- data.frame(
        geneId         = ann$gene_id[match(rownames(res), ann$gene_id)],
        geneName       = ann$gene_name[match(rownames(res), ann$gene_id)],
        baseMean       = rowMeans(cts[, c(idx$A, idx$B)]),
        baseMeanA      = rowMeans(cts[, idx$A]),
        baseMeanB      = rowMeans(cts[, idx$B]),
        foldChange     = 2 ^ res$table$logFC,
        log2FoldChange = res$table$logFC,
        PValue         = res$table$PValue,
        PAdj           = p.adjust(res$table$PValue, method = "hochberg"),
        FDR            = res$table$FDR,
        falsePos       = round(seq_along(res$table$FDR) * res$table$FDR)
    )
    
    tbl <- merge(dat, cts[, c(idx$A, idx$B)], by.x = "geneId", by.y = "row.names", sort = FALSE)

    write.csv(tbl, file = output$csv, quote = FALSE, row.names = FALSE)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@wildcards)