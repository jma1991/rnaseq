# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule limma_object:
    input:
        rds = "results/{result}/object.rds"
    output:
        rds = "results/limma/object.{result}.rds"
    log:
        out = "results/limma/object.{result}.out",
        err = "results/limma/object.{result}.err"
    message:
        "[limma] Create a DGEList object from {wildcards.result} output: {input.rds}"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/limma.R"

rule limma_voom:
    input:
        rds = "results/limma/object.{result}.rds"
    output:
        rds = "results/limma/voom.{result}.rds"
    log:
        out = "results/limma/voom.{result}.out",
        err = "results/limma/voom.{result}.err"
    message:
        "[limma] Create a EList object from {wildcards.result} output: {input.rds}"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/voom.R"

rule limma_counts:
    input:
        rds = "results/limma/object.{result}.rds"
    output:
        csv = "results/limma/counts.{result}.csv"
    log:
        out = "results/limma/counts.{result}.out",
        err = "results/limma/counts.{result}.err"
    message:
        "[limma] Write counts matrix to disk"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/counts.R"

rule limma_normcounts:
    input:
        rds = "results/limma/voom.{result}.rds"
    output:
        csv = "results/limma/normcounts.{result}.csv"
    log:
        out = "results/limma/normcounts.{result}.out",
        err = "results/limma/normcounts.{result}.err"
    message:
        "[limma] Write normalized counts matrix to disk"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/normcounts.R"

rule limma_logcounts:
    input:
        rds = "results/limma/voom.{result}.rds"
    output:
        csv = "results/limma/logcounts.{result}.csv"
    log:
        out = "results/limma/logcounts.{result}.out",
        err = "results/limma/logcounts.{result}.err"
    message:
        "[limma] Write log normalized counts matrix to disk"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/logcounts.R"

rule limma_results:
    input:
        rds = "results/limma/voom.{result}.rds",
        tsv = expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = config["genome"])
    output:
        csv = "results/limma/results_{A}_vs_{B}.{result}.csv"
    log:
        out = "results/limma/results_{A}_vs_{B}.{result}.out",
        err = "results/limma/results_{A}_vs_{B}.{result}.err"
    message:
        "[limma] Extract a table of the top-ranked genes: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/toptable.R"
