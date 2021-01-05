# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

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

rule fastq:
    input:
        fastq_input
    output:
        "results/fastq/{sample}/{unit}{read,.*}.fastq.gz"
    shell:
        "ln -s {input} {output}"
