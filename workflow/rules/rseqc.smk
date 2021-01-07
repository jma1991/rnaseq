# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule rseqc_bam_stat:
    input:
        "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/rseqc/bam_stat/{sample}.bam_stat.txt"
    log:
        "results/rseqc/bam_stat/{sample}.bam_stat.log"
    message:
        "[RSeQC] Summarize mapping statistics of BAM file"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} 1> {output} 2> {log}"

rule rseqc_inner_distance:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bed = expand("results/genomepy/{genome}/{genome}.annotation.bed", genome = config["txome"]["genome"])
    output:
        ext = multiext("results/rseqc/inner_distance/{sample}", ".inner_distance.txt", ".inner_distance_freq.txt", ".inner_distance_plot.pdf", ".inner_distance_plot.r")
    log:
        out = "results/rseqc/inner_distance/{sample}.inner_distance.out",
        err = "results/rseqc/inner_distance/{sample}.inner_distance.err"
    params:
        out = "results/rseqc/inner_distance/{sample}"
    message:
        "[RSeQC] Calculate inner distance between read pairs"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -i {input.bam} -r {input.bed} -o {params.out} 1> {log.out} 2> {log.err}"

rule rseqc_infer_experiment:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bed = expand("results/genomepy/{genome}/{genome}.annotation.bed", genome = config["txome"]["genome"])
    output:
        txt = "results/rseqc/infer_experiment/{sample}.infer_experiment.txt"
    log:
        err = "results/rseqc/infer_experiment/{sample}.infer_experiment.log"
    message:
        "[RSeQC] Infer experiment"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -i {input.bam} -r {input.bed} 1> {output.txt} 2> {log.err}"

rule rseqc_junction_annotation:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bed = expand("results/genomepy/{genome}/{genome}.annotation.bed", genome = config["txome"]["genome"])
    output:
        ext = multiext("results/rseqc/junction_annotation/{sample}", ".junction.Interact.bed", ".junction.bed", ".junction.xls", ".junction_plot.r", ".splice_events.pdf", ".splice_junction.pdf")
    log:
        err = "results/rseqc/junction_annotation/{sample}.junction_annotation.log"
    params:
        out = "results/rseqc/junction_annotation/{sample}"
    message:
        "[RSeQC] Compare detected splice junctions to reference gene model"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py -i {input.bam} -r {input.bed} -o {params.out} 2> {log.err}"

rule rseqc_junction_saturation:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bed = expand("results/genomepy/{genome}/{genome}.annotation.bed", genome = config["txome"]["genome"])
    output:
        ext = multiext("results/rseqc/junction_saturation/{sample}", ".junctionSaturation_plot.pdf", ".junctionSaturation_plot.r")
    log:
        err = "results/rseqc/junction_saturation/{sample}.junction_saturation.log"
    params:
        out = "results/rseqc/junction_saturation/{sample}"
    message:
        "[RSeQC] Check sequencing depth for alternative splicing"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py -i {input.bam} -r {input.bed} -o {params.out} 2> {log.err}"
        
rule rseqc_read_distribution:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bed = expand("results/genomepy/{genome}/{genome}.annotation.bed", genome = config["txome"]["genome"])
    output:
        txt = "results/rseqc/read_distribution/{sample}.read_distribution.txt"
    log:
        err = "results/rseqc/read_distribution/{sample}.read_distribution.log"
    message:
        "[RSeQC] Calculate read distribution over genome feature"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} 1> {output.txt} 2> {log.err}"

rule rseqc_read_duplication:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        ext = multiext("results/rseqc/read_duplication/{sample}", ".DupRate_plot.pdf", ".DupRate_plot.r", ".pos.DupRate.xls", ".seq.DupRate.xls")
    log:
        err = "results/rseqc/read_duplication/{sample}.read_duplication.log"
    params:
        out = "results/rseqc/read_duplication/{sample}"
    message:
        "[RSeQC] Calculate read duplication rate"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input.bam} -o {params.out} 2> {log.err}"


        
rule rseqc_gene_body_coverage:
    input:
        bam = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam",
        bed = expand("results/genomepy/{genome}/{genome}.annotation.bed", genome = config["txome"]["genome"])
    output:
        ext = multiext("results/rseqc/{sample}", ".geneBodyCoverage.curves.pdf", ".geneBodyCoverage.r", ".geneBodyCoverage.txt")
    log:
        err = "results/rseqc/{sample}.geneBodyCoverage.log"
    params:
        out = "results/rseqc/{sample}"
    message:
        "[RSeQC] Calculate read coverage over gene body"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.out} 2> {log.err}"

rule rseqc_read_gc:
    input:
        bam = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam"
    output:
        ext = multiext("results/rseqc/{sample}", ".GC.xls", ".GC_plot.pdf", ".GC_plot.r")
    log:
        err = "results/rseqc/{sample}.GC.log"
    params:
        out = "results/rseqc/{sample}"
    message:
        "[RSeQC] Calculate GC content distribution of reads"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input.bam} -o {params.out} 2> {log.err}"
