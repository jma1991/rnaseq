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
    u = project.units.loc[project.units["sample"] == sample, ]
    m = [is_single_end(sample, x) for x in u["unit"]]
    return all(m)

def all_paired_end(sample):
    u = project.units.loc[project.units["sample"] == sample, ]
    m = [is_paired_end(sample, x) for x in u["unit"]]
    return all(m)

def is_single_end(sample, unit):
    ind = (project.units["sample"] == sample) & (project.units["unit"] == unit)
    obj = project.units.loc[ind, "read2"].item()
    return pd.isnull(obj)

def is_paired_end(sample, unit):
    ind = (project.units["sample"] == sample) & (project.units["unit"] == unit)
    obj = project.units.loc[ind, "read2"].item()
    return pd.notnull(obj)

def get_fastqs(wildcards):
    if is_single_end(wildcards.sample, wildcards.unit):
        u = project.units.loc[(wildcards.sample, wildcards.unit), "read1"]
        return {"fq1": u}
    else:
        u = project.units.loc[(wildcards.sample, wildcards.unit), ["read1", "read2"]].dropna()
        return {"fq1": u.read1, "fq2": u.read2}

def qualimap_rnaseq_params(wildcards):
    
    arg = []

    d = {"U": "non-strand-specific", "F": "strand-specific-forward", "R": "strand-specific-reverse"}
    k = project.samples.loc[project.samples["sample"] == wildcards.sample, "stranded"].item()
    v = "-p" + " " + d.get(k)
    arg.append(v)
    
    if all_single_end(wildcards.sample):
        arg.append("")
    elif all_paired_end(wildcards.sample):
        arg.append("-pe")
    else:
        raise ValueError("You weren't supposed to be able to get here you know.")

    return " ".join(arg)

def preseq_params(wildcards):
    if all_single_end(wildcards.sample):
        return ""
    if all_paired_end(wildcards.sample):
        return "-P"
    raise ValueError("You weren't supposed to be able to get here you know.")

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
