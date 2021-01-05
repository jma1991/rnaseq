# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

ruleorder: qualimap_rnaseq_pe > qualimap_rnaseq_se

rule qualimap_rnaseq_se:
    input:
        bam = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam",
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = GENOME)
    output:
        dir = directory("results/qualimap/{sample}")
    params:
        arg = qualimap_params
    log:
        out = "results/qualimap/{sample}/qualimap_rnaseq.out",
        err = "results/qualimap/{sample}/qualimap_rnaseq.err"
    message:
        "[Qualimap] Perform RNA-seq QC analysis: {input.bam}"
    wildcard_constraints:
        sample = "|".join(runs.loc[SINGLE, "sample"])
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap rnaseq -bam {input.bam} -gtf {input.gtf} -outdir {output.dir} {params.arg} 1> {log.out} 2> {log.err}"

rule qualimap_rnaseq_pe:
    input:
        bam = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam",
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = GENOME)
    output:
        dir = directory("results/qualimap/{sample}")
    params:
        arg = qualimap_params
    log:
        out = "results/qualimap/{sample}/qualimap_rnaseq.out",
        err = "results/qualimap/{sample}/qualimap_rnaseq.err"
    message:
        "[Qualimap] Perform RNA-seq QC analysis: {input.bam}"
    wildcard_constraints:
        sample = "|".join(runs.loc[PAIRED, "sample"])
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap rnaseq -bam {input.bam} -gtf {input.gtf} -outdir {output.dir} {params.arg} -pe 1> {log.out} 2> {log.err}"
