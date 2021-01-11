# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule dist:
    input:
        csv = "results/{result}/logcounts.{quant}.csv"
    output:
        pdf = "results/{result}/dist.{quant}.pdf"
    log:
        out = "results/{result}/dist.{quant}.out",
        err = "results/{result}/dist.{quant}.err"
    message:
        "[Plotting] Sample distance"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/dist.R"

rule prcomp:
    input:
        csv = "results/{result}/logcounts.{quant}.csv"
    output:
        pdf = "results/{result}/prcomp.{quant}.pdf"
    log:
        out = "results/{result}/prcomp.{quant}.out",
        err = "results/{result}/prcomp.{quant}.err"
    message:
        "[Plotting] Principal components analysis"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/prcomp.R"

rule cmdscale:
    input:
        csv = "results/{result}/logcounts.{quant}.csv"
    output:
        pdf = "results/{result}/cmdscale.{quant}.pdf"
    log:
        out = "results/{result}/cmdscale.{quant}.out",
        err = "results/{result}/cmdscale.{quant}.err"
    message:
        "[Plotting] Multi-dimensional scaling"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/cmdscale.R"

rule diff:
    input:
        csv = "results/{result}/results_{A}_vs_{B}.{quant}.csv"
    output:
        pdf = "results/{result}/results_{A}_vs_{B}.{quant}.diff.pdf"
    log:
        out = "results/{result}/results_{A}_vs_{B}.{quant}.diff.out",
        err = "results/{result}/results_{A}_vs_{B}.{quant}.diff.err"
    message:
        "[Plotting] Sample difference"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/diff.R"

rule pvalue:
    input:
        csv = "results/{result}/results_{A}_vs_{B}.{quant}.csv"
    output:
        pdf = "results/{result}/results_{A}_vs_{B}.{quant}.pvalue.pdf"
    log:
        out = "results/{result}/results_{A}_vs_{B}.{quant}.pvalue.out",
        err = "results/{result}/results_{A}_vs_{B}.{quant}.pvalue.err"
    message:
        "[Plotting] Plot p-values of differentially expressed genes: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/pvalue.R"

rule volcano:
    input:
        csv = "results/{result}/results_{A}_vs_{B}.{quant}.csv"
    output:
        pdf = "results/{result}/results_{A}_vs_{B}.{quant}.volcano.pdf"
    log:
        out = "results/{result}/results_{A}_vs_{B}.{quant}.volcano.out",
        err = "results/{result}/results_{A}_vs_{B}.{quant}.volcano.err"
    message:
        "[Plotting] Plot volcano of differentially expressed genes: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/volcano.R"

rule heatmap:
    input:
        csv = "results/{result}/results_{A}_vs_{B}.{quant}.csv"
    output:
        pdf = "results/{result}/results_{A}_vs_{B}.{quant}.heatmap.pdf"
    log:
        out = "results/{result}/results_{A}_vs_{B}.{quant}.heatmap.out",
        err = "results/{result}/results_{A}_vs_{B}.{quant}.heatmap.err"
    message:
        "[Plotting] Plot heatmap of differentially expressed genes: {wildcards.A} vs {wildcards.B}"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/heatmap.R"
