# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule fastq:
    input:
        fqz = fastq_input
    output:
        fqz = "results/fastq/{sample}/{unit}{read,.*}.fastq.gz"
    message:
        "[coreutils] Create a symbolic link"
    threads:
        1
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -s {input.fqz} {output.fqz}"
