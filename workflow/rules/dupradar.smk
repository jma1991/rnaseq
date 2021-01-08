# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

def dupradar_params_stranded(wildcards):

    k = project.samples.loc[project.samples["sample"] == wildcards.sample, "stranded"].item()

    d = {"U": 0, "F": 1, "R": 2}

    v = d.get(k)

    return v


def dupradar_params_paired(wildcards):

    if all_single_end(wildcards.sample):

        return False

    if all_paired_end(wildcards.sample):

        return True

    return ValueError("ValueError")


rule dupradar_analyze:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = config["genome"])
    output:
        csv = "results/dupradar/{sample}.analyzeDuprates.csv"
    log:
        out = "results/dupradar/{sample}.analyzeDuprates.out",
        err = "results/dupradar/{sample}.analyzeDuprates.err"
    params:
        stranded = dupradar_params_stranded,
        paired = dupradar_params_paired
    conda:
        "../envs/dupradar.yaml"
    script:
        "../scripts/analyzeDuprates.R"

rule dupradar_densplot:
    input:
        csv = "results/dupradar/{sample}.analyzeDuprates.csv"
    output:
        pdf = "results/dupradar/{sample}.duprateExpDensPlot.pdf"
    log:
        out = "results/dupradar/{sample}.duprateExpDensPlot.out",
        err = "results/dupradar/{sample}.duprateExpDensPlot.err"
    conda:
        "../envs/dupradar.yaml"
    script:
        "../scripts/duprateExpDensPlot.R"
