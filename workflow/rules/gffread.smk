# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule gffread_transcripts:
    input:
        "results/genomepy/{genome}/{genome}.fa",
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.transcripts.fa"
    log:
        "results/gffread/{genome}/{genome}.transcripts.log"
    message:
        "[gffread] Extract transcript sequences from {wildcards.genome} reference genome"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -w {output} -g {input} 2> {log}"

rule gffread_tx2gene:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.tx2gene.tsv"
    log:
        "results/gffread/{genome}/{genome}.tx2gene.log"
    message:
        "[gffread] Output a tx2gene annotation table from {wildcards.genome} reference genome"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table transcript_id,gene_id,gene_name | sort -u 1> {output} 2> {log}"

rule gffread_id2name:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.id2name.tsv"
    log:
        "results/gffread/{genome}/{genome}.id2name.log"
    message:
        "[gffread] Output a id2name annotation table from {wildcards.genome} reference genome"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table gene_id,gene_name | sort -u 1> {output} 2> {log}"
