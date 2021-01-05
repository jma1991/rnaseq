# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

def preseq_params(wildcards):
    return None

rule preseq_c_curve:
    input:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        txt = "results/preseq/{sample}_complexity_curve.txt"
    log:
        "results/preseq/{sample}_complexity_curve.log"
    params:
        len = 50818468
    message:
        "[preseq] Compute complexity curve"
    conda:
        "../envs/preseq.yaml"
    shell:
        "preseq c_curve -o {output.txt} -l {params.len} -s 100 -P -B {input.bam} > {log}"

rule preseq_lc_extrap:
    input:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        txt = "results/preseq/{sample}_future_yield.txt"
    log:
        "results/preseq/{sample}_future_yield.log"
    params:
        ext = 1e+08,
        len = 50818468
    message:
        "[preseq] Estimate future yield"
    conda:
        "../envs/preseq.yaml"
    shell:
        "preseq lc_extrap -o {output.txt} -e {params.ext} -l {params.len} -B -P {input.bam} > {log}"

rule library_complexity:
    input:
        txt = ["results/preseq/{sample}_complexity_curve.txt", "results/preseq/{sample}_future_yield.txt"]
    output:
        pdf = "results/preseq/{sample}_library_complexity.pdf"
    script:
        "../scripts/library_complexity.R"

rule preseq_bound_pop:
    input:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        txt = "results/preseq/{sample}_species_richness.txt"
    log:
        "results/preseq/{sample}_species_richness.log"
    params:
        len = 50818468
    message:
        "[preseq] Estimate species richness"
    conda:
        "../envs/preseq.yaml"
    shell:
        "preseq bound_pop -o {output.txt} -l {params.len} -P -B {input.bam} > {log}"

rule species_richness:
    input:
        txt = expand("results/preseq/{sample}_species_richness.txt", sample = SAMPLES)
    output:
        pdf = "results/preseq/species_richness.pdf"
    params:
        ids = SAMPLES
    script:
        "../scripts/species_richness.R"