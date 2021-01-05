# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule cutadapt_single:
    input:
        unpack(get_fastqs)
    output:
        fq1 = "results/cutadapt/{sample_name}/{subsample_name}.fastq.gz"
    log:
        out = "results/cutadapt/{sample_name}/{subsample_name}_cutadapt.out",
        err = "results/cutadapt/{sample_name}/{subsample_name}_cutadapt.err"
    message:
        "[cutadapt] Remove adapter sequences from single-end library: {input.fq1}"
    threads:
        1
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 33 -e 0.005 -O 7 -o {output.fq1} {input.fq1} 1> {log.out} 2> {log.err}"

rule cutadapt_paired:
    input:
        unpack(get_fastqs)
    output:
        fq1 = "results/cutadapt/{sample_name}/{subsample_name}_1.fastq.gz",
        fq2 = "results/cutadapt/{sample_name}/{subsample_name}_2.fastq.gz"
    log:
        out = "results/cutadapt/{sample_name}/{subsample_name}_cutadapt.out",
        err = "results/cutadapt/{sample_name}/{subsample_name}_cutadapt.err"
    message:
        "[cutadapt] Remove adapter sequences from paired-end library: {input.fq1} & {input.fq2}"
    threads:
        1
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -j {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 33 -e 0.005 -O 7 -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} 1> {log.out} 2> {log.err}"
