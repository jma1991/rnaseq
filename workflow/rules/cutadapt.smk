# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule cutadapt_se:
    input:
        fq1 = "results/fastq/{sample}/{unit}.fastq.gz"
    output:
        fq1 = "results/cutadapt/{sample}/{unit}.fastq.gz"
    log:
        out = "results/cutadapt/{sample}/{unit}_cutadapt.out",
        err = "results/cutadapt/{sample}/{unit}_cutadapt.err"
    message:
        "[cutadapt] Remove adapter sequences from single-end library: {input.fq1}"
    threads:
        1
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 33 -e 0.005 -O 7 -o {output.fq1} {input.fq1} 1> {log.out} 2> {log.err}"

rule cutadapt_pe:
    input:
        fq1 = "results/fastq/{sample}/{unit}_1.fastq.gz",
        fq2 = "results/fastq/{sample}/{unit}_2.fastq.gz"
    output:
        fq1 = "results/cutadapt/{sample}/{unit}_1.fastq.gz",
        fq2 = "results/cutadapt/{sample}/{unit}_2.fastq.gz"
    log:
        out = "results/cutadapt/{sample}/{unit}_cutadapt.out",
        err = "results/cutadapt/{sample}/{unit}_cutadapt.err"
    message:
        "[cutadapt] Remove adapter sequences from paired-end library: {input.fq1} & {input.fq2}"
    threads:
        1
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 33 -e 0.005 -O 7 -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} 1> {log.out} 2> {log.err}"
