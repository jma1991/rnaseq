# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule deseq2_object:
    input:
        rds = "results/{result}/object.rds"
    output:
        rds = "results/deseq2/object.{result}.rds"
    log:
        out = "results/deseq2/object.{result}.out",
        err = "results/deseq2/object.{result}.err"
    message:
        "[DESeq2] Create a DESeqDataSet object from {wildcards.result} object"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"

rule deseq2_counts:
    input:
        rds = "results/deseq2/object.{result}.rds"
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
        rds = "results/deseq2/object.{result}.rds"
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
        rds = "results/deseq2/object.{result}.rds"
    output:
        csv = "results/deseq2/logcounts.{result}.csv"
    log:
        out = "results/deseq2/logcounts.{result}.out",
        err = "results/deseq2/logcounts.{result}.err"
    message:
        "[DESeq2] Write log-normalized counts matrix to disk"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/logcounts.R"

rule deseq2_results:
    input:
        rds = "results/deseq2/object.{result}.rds",
        tsv = expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = config["genome"])
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
