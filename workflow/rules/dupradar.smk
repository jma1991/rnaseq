# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule dupradar_analyze_duprates:
    input:
        bam = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam",
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = GENOME)
    output:
        csv = "results/dupradar/{sample}.analyzeDuprates.csv"
    log:
        "results/dupradar/{sample}.analyzeDuprates.log"
    params:
        stranded = dupradar_stranded
    conda:
        "../envs/dupradar.yaml"
    script:
        "../scripts/analyzeDuprates.R"

rule dupradar_duprate_densplot:
    input:
        csv = "results/dupradar/{sample}.analyzeDuprates.csv"
    output:
        pdf = "results/dupradar/{sample}.duprateExpDensPlot.pdf"
    log:
        "results/dupradar/{sample}.duprateExpDensPlot.log"
    conda:
        "../envs/dupradar.yaml"
    script:
        "../scripts/duprateExpDensPlot.R"

rule dupradar_expression_hist:
    input:
        csv = "results/dupradar/{sample}.analyzeDuprates.csv"
    output:
        pdf = "results/dupradar/{sample}.expressionHist.pdf"
    log:
        "results/dupradar/{sample}.expressionHist.log"
    conda:
        "../envs/dupradar.yaml"
    script:
        "../scripts/expressionHist.R"

rule dupradar_duprate_boxplot:
    input:
        csv = "results/dupradar/{sample}.analyzeDuprates.csv"
    output:
        pdf = "results/dupradar/{sample}.duprateExpBoxplot.pdf"
    log:
        "results/dupradar/{sample}.duprateExpBoxplot.log"
    conda:
        "../envs/dupradar.yaml"
    script:
        "../scripts/duprateExpBoxplot.R"
