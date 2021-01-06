# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

import pandas as pd

wildcard_constraints:
    sample = "|".join(project.samples["sample"]),
    unit = "|".join(project.units["unit"])

def fastq_input(wildcards):

    df = project.units

    ix = (df["sample"] == wildcards.sample) & (df["unit"] == wildcards.unit)

    df = df.loc[ix, ]

    if not wildcards.read:

        return df["read1"]

    elif wildcards.read == "_1":

        return df["read1"]

    elif wildcards.read == "_2":

        return df["read2"]

    else:

        raise ValueError("You weren't supposed to be able to get here you know.")

def all_single_end(sample):
    u = project.units.loc[sample, ]
    m = [is_single_end(sample, x) for x in u["unit"]]
    return all(m)

def all_paired_end(sample):
    u = project.units.loc[sample, ]
    m = [is_paired_end(sample, x) for x in u["unit"]]
    return all(m)

def is_single_end(sample, unit):
    obj = project.units.loc[(sample, unit), "read2"]
    if isinstance(obj, pd.core.series.Series):
        raise ValueError()
    return pd.isnull(obj)

def is_paired_end(sample, unit):
    obj = project.units.loc[(sample, unit), "read2"]
    if isinstance(obj, pd.core.series.Series):
        raise ValueError()
    return pd.notnull(obj)

def get_fastqs(wildcards):
    if is_single_end(wildcards.sample, wildcards.unit):
        u = project.units.loc[(wildcards.sample, wildcards.unit), "read1"]
        return {"fq1": u}
    else:
        u = project.units.loc[(wildcards.sample, wildcards.unit), ["read1", "read2"]].dropna()
        return {"fq1": u.read1, "fq2": u.read2}

def kallisto_quant_index(wildcards):
    x = project.samples.loc[wildcards.sample, "genome"]
    return f"results/kallisto/index/{x}/{x}.idx"

def kallisto_quant_fastq(wildcards):
    if all_single_end(wildcards.sample):
        return expand("results/cutadapt/{sample}/{unit}.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"])
    elif all_paired_end(wildcards.sample):
        return expand("results/cutadapt/{sample}/{unit}_{read_name}.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"], read_name = ["1", "2"])
    else:
        raise ValueError()

def kallisto_quant_stranded(sample):
    d = {"U": "", "F": "--fr-stranded", "R": "--rf-stranded"}
    k = project.samples.loc[sample, "stranded"]
    return d.get(k)

def kallisto_quant_arguments(wildcards):
    str = kallisto_quant_stranded(wildcards.sample)
    if all_single_end(wildcards.sample):
        return f"--single {str} --fragment-length=200 --sd=20"
    elif all_paired_end(wildcards.sample):
        return f"{str}"
    else:
        raise ValueError()

def star_align_index(wildcards):
    x = project.samples.loc[wildcards.sample, "genome"]
    return f"results/star/index/{x}"

def star_align_fastq(wildcards):
    if all_single_end(wildcards.sample):
        return expand("results/cutadapt/{sample}/{unit}.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"])
    elif all_paired_end(wildcards.sample):
        return expand("results/cutadapt/{sample}/{unit}_{read_name}.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"], read_name = ["1", "2"])
    else:
        raise ValueError()

def star_align_reads(wildcards):
    if all_single_end(wildcards.sample):
        fqz = expand("results/cutadapt/{sample}/{unit}.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"])
        return ",".join(fqz)
    elif all_paired_end(wildcards.sample):
        fq1 = expand("results/cutadapt/{sample}/{unit}_1.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"])
        fq2 = expand("results/cutadapt/{sample}/{unit}_2.fastq.gz", sample = wildcards.sample, unit = project.units.loc[wildcards.sample, "unit"])
        return ",".join(fq1) + " " + ",".join(fq2)
    else:
        raise ValueError()
