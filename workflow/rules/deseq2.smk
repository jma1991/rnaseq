# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule deseq2_init:
    input:
        rds = "results/{result}/counts.rds"
    output:
        rds = "results/deseq2/dds.{result}.rds"
    log:
        out = "results/deseq2/dds.{result}.out",
        err = "results/deseq2/dds.{result}.err"
    message:
        "[DESeq2] Create a DESeqDataSet object from {wildcards.result} output: {input.rds}"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"

rule deseq2_counts:
    input:
        rds = "results/deseq2/dds.{result}.rds"
    output:
        csv = "results/deseq2/counts.{result}.csv"
    log:
        out = "results/deseq2/counts.{result}.out",
        err = "results/deseq2/counts.{result}.err"
    message:
        "[DESeq2] Write counts matrix to disk"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/counts.R"

rule deseq2_normcounts:
    input:
        rds = "results/deseq2/dds.{result}.rds"
    output:
        csv = "results/deseq2/normcounts.{result}.csv"
    log:
        out = "results/deseq2/normcounts.{result}.out",
        err = "results/deseq2/normcounts.{result}.err"
    message:
        "[DESeq2] Write normalized counts matrix to disk"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/normcounts.R"

rule deseq2_logcounts:
    input:
        rds = "results/deseq2/dds.{result}.rds"
    output:
        csv = "results/deseq2/logcounts.{result}.csv"
    log:
        out = "results/deseq2/logcounts.{result}.out",
        err = "results/deseq2/logcounts.{result}.err"
    message:
        "[DESeq2] Write log normalized counts matrix to disk"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/logcounts.R"

rule deseq2_results:
    input:
        rds = "results/deseq2/dds.{result}.rds",
        tsv = expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = pep.sample_table["genome"].unique())
    output:
        csv = "results/deseq2/condition_{A}_vs_{B}.{result}.csv"
    log:
        out = "results/deseq2/condition_{A}_vs_{B}.{result}.out",
        err = "results/deseq2/condition_{A}_vs_{B}.{result}.err"
    message:
        "[DESeq2] Extract results from a DESeq analysis: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/results.R"
