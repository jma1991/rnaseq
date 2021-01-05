# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule tximport:
    input:
        hdf = expand("results/kallisto/quant/{sample}/abundance.h5", sample = pep.sample_table.index),
        tsv = expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = pep.sample_table["genome"].unique())
    output:
        rds = "results/tximport/counts.rds"
    log:
        out = "results/tximport/log.out",
        err = "results/tximport/log.err"
    message:
        "[tximport] Import and summarize transcript-level estimates from Kallisto"
    conda:
        "../envs/tximport.yaml"
    script:
        "../scripts/tximport.R"
