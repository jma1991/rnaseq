# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule kallisto_index:
    input:
        fas = "results/gffread/{genome}/{genome}.transcripts.fa"
    output:
        idx = "results/kallisto/index/{genome}/{genome}.idx"
    log:
        out = "results/kallisto/index/{genome}/{genome}.out",
        err = "results/kallisto/index/{genome}/{genome}.err"
    message:
        "[Kallisto] Build a Kallisto index: {wildcards.genome}"
    threads:
        16
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output.idx} {input.fas} 1> {log.out} 2> {log.err}"

rule kallisto_quant:
    input:
        idx = expand("results/kallisto/index/{genome}/{genome}.idx", genome = config["genome"]),
        fqz = kallisto_quant_fastq
    output:
        ext = multiext("results/kallisto/quant/{sample}/", "abundance.h5", "abundance.tsv", "run_info.json")
    log:
        out = "results/kallisto/quant/{sample}/kallisto_quant.out",
        err = "results/kallisto/quant/{sample}/kallisto_quant.err"
    params:
        out = "results/kallisto/quant/{sample}",
        arg = kallisto_quant_arguments
    message:
        "[kallisto]"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {params.out} -t {threads} {params.arg} {input.fqz} 1> {log.out} 2> {log.err}"
