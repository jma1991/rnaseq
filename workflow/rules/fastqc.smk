# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule fastqc_fqz:
    input:
        fqz = "results/fastq/{sample}/{run}.fastq.gz"
    output:
        ext = multiext("results/fastqc/{sample}/{run}_fastqc", ".html", ".zip")
    params:
        dir = "results/fastqc/{sample}"
    log:
        out = "results/fastqc/{sample}/{run}_fastqc.out",
        err = "results/fastqc/{sample}/{run}_fastqc.err"
    message:
        "[fastqc] Quality control raw sequence data: {input.fqz}"
    threads:
        4
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc -o {params.dir} -t {threads} {input.fqz} 1> {log.out} 2> {log.err}"

rule fastqc_fq1:
    input:
        fq1 = "results/fastq/{sample}/{run}_1.fastq.gz"
    output:
        ext = multiext("results/fastqc/{sample}/{run}_1_fastqc", ".html", ".zip")
    params:
        dir = "results/fastqc/{sample}"
    log:
        out = "results/fastqc/{sample}/{run}_1_fastqc.out",
        err = "results/fastqc/{sample}/{run}_1_fastqc.err"
    message:
        "[fastqc] Quality control raw sequence data: {input.fq1}"
    threads:
        4
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc -o {params.dir} -t {threads} {input.fq1} 1> {log.out} 2> {log.err}"

rule fastqc_fq2:
    input:
        fq2 = "results/fastq/{sample}/{run}_2.fastq.gz"
    output:
        ext = multiext("results/fastqc/{sample}/{run}_2_fastqc", ".html", ".zip")
    params:
        dir = "results/fastqc/{sample}"
    log:
        out = "results/fastqc/{sample}/{run}_2_fastqc.out",
        err = "results/fastqc/{sample}/{run}_2_fastqc.err"
    message:
        "[fastqc] Quality control raw sequence data: {input.fq2}"
    threads:
        4
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc -o {params.dir} -t {threads} {input.fq2} 1> {log.out} 2> {log.err}"
