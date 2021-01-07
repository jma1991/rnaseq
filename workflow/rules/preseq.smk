# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule preseq_bam2mr:
    input:
        "results/star/align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/preseq/bam2mr/{sample}.mr"
    shell:
        "to-mr {input} | sort -k 1,1, -k 2,2n -k 3,3n > {output}"

rule preseq_c_curve:
    input:
        "results/preseq/bam2mr/{sample}.mr"
    output:
        "results/preseq/c_curve/{sample}.c_curve.txt"
    log:
        "results/preseq/c_curve/{sample}.c_curve.log"
    params:
        preseq_params
    message:
        "[preseq] Compute complexity curve"
    conda:
        "../envs/preseq.yaml"
    shell:
        "preseq c_curve -o {output} {params} {input} 2> {log}"
        
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
