# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule star_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.fa",
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        dir = directory("results/star/index/{genome}")
    message:
        "[STAR] Creating a STAR index: {wildcards.genome}"
    threads:
        16
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fas} --sjdbGTFfile {input.gtf} --sjdbOverhang 74"

def star_align_fastq(wildcards):

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


def star_param_fastq(wildcards):

    df = project.units

    ix = df["sample"] == wildcards.sample

    df = df[ix]

    se = df["read2"].isnull().all()

    pe = df["read2"].notnull().all()

    if se:
        r1 = expand("results/cutadapt/{sample}/{unit}.fastq.gz", sample = wildcards.sample, unit = df["unit"])
        return r1

    if pe:
        r1 = expand("results/cutadapt/{sample}/{unit}_1.fastq.gz", sample = wildcards.sample, unit = df["unit"])
        r2 = expand("results/cutadapt/{sample}/{unit}_2.fastq.gz", sample = wildcards.sample, unit = df["unit"])
        return ",".join(r1) + " " + ",".join(r2)

    raise ValueError("You weren't supposed to be able to get here you know.")


rule star_align:
    input:
        idx = expand("results/star/index/{genome}", genome = config["genome"]),
        fqz = star_align_fastq
    output:
        bam = "results/star/align/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        fqz = star_param_fastq,
        out = "results/star/align/{sample}/"
    message:
        "[STAR] Align single-end library to genome: {wildcards.sample}"
    threads:
        16
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {params.fqz} --readFilesCommand gunzip -c --outFileNamePrefix {params.out} --outSAMtype BAM SortedByCoordinate --outTmpDir /tmp/TMPDIR/{wildcards.sample}"
