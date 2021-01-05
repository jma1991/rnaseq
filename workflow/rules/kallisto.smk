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
    shell:
        "kallisto index -i {output.idx} {input.fas} 1> {log.out} 2> {log.err}"

def kallisto_quant_fastq(wildcards):

    df = project.units

    ix = df["sample"] == wildcards.sample

    df = df[ix]

    se = df["read2"].isnull().all()

    pe = df["read2"].notnull().all()

    if se:
        return expand("results/cutadapt/{sample}/{unit}.fastq.gz", sample = wildcards.sample, unit = df["unit"])

    if pe:
        return expand("results/cutadapt/{sample}/{unit}{read}.fastq.gz", sample = wildcards.sample, unit = df["unit"], read = ["_1", "_2"])
    
    raise ValueError("You weren't supposed to be able to get here you know.")


def kallisto_quant_stranded(sample):

    d = {"U": "", "F": "--fr-stranded", "R": "--rf-stranded"}

    i = project.samples["sample"] == sample

    k = project.samples.loc[i, "stranded"].item()

    return d.get(k)

def kallisto_quant_arguments(wildcards):

    df = project.units

    ix = df["sample"] == wildcards.sample

    df = df[ix]

    se = df["read2"].isnull().all()

    pe = df["read2"].notnull().all()

    st = kallisto_quant_stranded(wildcards.sample)

    if se:
        return f"--single {st} -l 200 -s 20"

    if pe:
        return f"{st}"

    raise ValueError("You weren't supposed to be able to get here you know.")


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
    shell:
        "kallisto quant -i {input.idx} -o {params.out} -t {threads} {params.arg} {input.fqz} 1> {log.out} 2> {log.err}"
