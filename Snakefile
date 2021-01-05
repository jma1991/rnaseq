# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

configfile: "config/config.yaml"

include: "rules/project.smk"
include: "rules/fastq.smk"
include: "rules/fastqc.smk"
include: "rules/cutadapt.smk"
include: "rules/genomepy.smk"
include: "rules/gffread.smk"

project = Project(config)

rule all:
    input:
        project.fastq_output(),
        project.fastqc_output(),
        project.cutadapt_output(),
        project.genomepy_output(),
        project.gffread_output()
