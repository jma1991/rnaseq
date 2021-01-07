# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule edger_object:
    input:
        rds = "results/{result}/object.rds"
    output:
        rds = "results/edger/object.{result}.rds"
    log:
        out = "results/edger/object.{result}.out",
        err = "results/edger/object.{result}.err"
    message:
        "[edgeR] Construct a DGEList object from {wildcards.result} output: {input.rds}"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/edger.R"

rule edger_counts:
    input:
        rds = "results/edger/object.{result}.rds"
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
        rds = "results/edger/object.{result}.rds"
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
        rds = "results/edger/object.{result}.rds"
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
        rds = "results/edger/object.{result}.rds",
        tsv = expand("results/gffread/{genome}/{genome}.id2name.tsv", genome = config["genome"])
    output:
        csv = "results/edger/results_{A}_vs_{B}.{result}.csv"
    log:
        out = "results/edger/results_{A}_vs_{B}.{result}.out",
        err = "results/edger/results_{A}_vs_{B}.{result}.err"
    message:
        "[edgeR] Compute gene-wise exact test: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/toptags.R"

