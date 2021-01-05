# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule gffcompare:
    input:
        ref = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = GENOME),
        qry = "results/stringtie/merged.gtf"
    output:
        ext = multiext("results/gffcompare/merged", ".annotated.gtf", ".loci", ".stats", ".tracking")
    params:
        out = "results/gffcompare/merged"
    message:
        "Examine how the transcripts compare with the reference annotation:"
    conda:
        "envs/gffcompare.yaml"
    shell:
        "gffcompare -r {input.ref} -o {params.out} {input.qry}"