# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule tximport:
    input:
        hdf = expand("results/kallisto/quant/{sample}/abundance.h5", sample = project.samples["sample"]),
        tsv = expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = config["genome"])
    output:
        rds = "results/tximport/object.rds"
    log:
        out = "results/tximport/object.out",
        err = "results/tximport/object.err"
    message:
        "[tximport] Import and summarize transcript-level estimates from Kallisto"
    conda:
        "../envs/tximport.yaml"
    script:
        "../scripts/tximport.R"
