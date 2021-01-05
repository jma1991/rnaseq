# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule edger_init:
    input:
        rds = "results/{result}/counts.rds"
    output:
        rds = "results/edger/dge.{result}.rds"
    log:
        out = "results/edger/dge.{result}.out",
        err = "results/edger/dge.{result}.err"
    message:
        "[edgeR] Construct a DGEList object from {wildcards.result} output: {input.rds}"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/edger.R"

rule edger_counts:
    input:
        rds = "results/edger/dge.{result}.rds"
    output:
        csv = "results/edger/counts.{result}.csv"
    log:
        out = "results/edger/counts.{result}.out",
        err = "results/edger/counts.{result}.err"
    message:
        "[edgeR] Write counts matrix to disk"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/counts.R"

rule edger_normcounts:
    input:
        rds = "results/edger/dge.{result}.rds"
    output:
        csv = "results/edger/normcounts.{result}.csv"
    log:
        out = "results/edger/normcounts.{result}.out",
        err = "results/edger/normcounts.{result}.err"
    message:
        "[edgeR] Write normalized counts matrix to disk"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/normcounts.R"

rule edger_logcounts:
    input:
        rds = "results/edger/dge.{result}.rds"
    output:
        csv = "results/edger/logcounts.{result}.csv"
    log:
        out = "results/edger/logcounts.{result}.out",
        err = "results/edger/logcounts.{result}.err"
    message:
        "[edgeR] Write log normalized counts matrix to disk"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/logcounts.R"

rule edger_results:
    input:
        rds = "results/edger/dge.{result}.rds",
        tsv = expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = pep.sample_table["genome"].unique())
    output:
        csv = "results/edger/condition_{A}_vs_{B}.{result}.csv"
    log:
        out = "results/edger/condition_{A}_vs_{B}.{result}.out",
        err = "results/edger/condition_{A}_vs_{B}.{result}.err"
    message:
        "[edgeR] Compute gene-wise exact test: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/exact.R"
