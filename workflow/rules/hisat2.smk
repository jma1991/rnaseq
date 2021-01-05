# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule hisat2_extract_splice_sites:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/genomepy/{genome}/index/hisat2/{genome}.ss"
    message:
        "[HISAT2] Extract splice-site information from the gene annotation file:"
    conda:
        "envs/hisat2.yaml"
    shell:
        "extract_splice_sites.py {input} > {output}"

rule hisat2_extract_exons:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/genomepy/{genome}/index/hisat2/{genome}.exon"
    message:
        "[HISAT2] Extract exon information from the gene annotation file:"
    conda:
        "envs/hisat2.yaml"
    shell:
        "extract_exons.py {input} > {output}"

rule hisat2_build:
    input:
        fas = "results/genomepy/{genome}/{genome}.fa",
        ss = "results/genomepy/{genome}/index/hisat2/{genome}.ss",
        exon = "results/genomepy/{genome}/index/hisat2/{genome}.exon"
    output:
        dir = directory("results/genomepy/{genome}/index/hisat2/{genome}")
    params:
        idx = "results/genomepy/{genome}/index/hisat2/{genome}"
    message:
        "[HISAT2] Build a HISAT2 index:"
    threads:
        16
    conda:
        "envs/hisat2.yaml"
    shell:
        "hisat2-build -p {threads} --ss {input.ss} --exon {input.exon} {input.fas} {params.idx}"

rule hisat2_align:
    input:
        idx = expand("results/genomepy/{genome}/index/hisat2/{genome}", genome = GENOME),
        fq1 = lambda wildcards: expand("results/cutadapt/{sample}/{run}_1.fastq.gz", sample = wildcards.sample, run = runs.loc[runs["sample"] == wildcards.sample, "run"]),
        fq2 = lambda wildcards: expand("results/cutadapt/{sample}/{run}_2.fastq.gz", sample = wildcards.sample, run = runs.loc[runs["sample"] == wildcards.sample, "run"])
    output:
        sam = "results/hisat2/{sample}.sam"
    message:
        "[HISAT2] Map the reads for each sample to the reference genome:"
    threads:
        4
    conda:
        "envs/hisat2.yaml"
    shell:
        "hisat2 -p {threads} -x {input.idx} -1 {input.fq1} -2 {input.fq2} -S {output.sam}"

rule samtools_sort:
    input:
        "results/hisat2/{sample}.sam"
    output:
        "results/hisat2/{sample}.bam"
    message:
        "[Samtools] Sort and convert SAM to BAM"
    threads:
        4
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule samtools_index:
    input:
        "results/hisat2/{sample}.bam"
    output:
        "results/hisat2/{sample}.bam.bai"
    message:
        "[Samtools] Index coordinate-sorted BAM"
    threads:
        4
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input}"
