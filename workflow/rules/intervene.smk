# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule csvtk_filter:
    input:
        csv = "results/{result}/condition_{A}_vs_{B}.{quant}.csv"
    output:
        csv = temp("results/{result}/condition_{A}_vs_{B}.{quant}.filtered.csv")
    params:
        fdr = 0.05
    message:
        "[csvtk] Filter results by FDR threshold: {params.fdr}"
    conda:
        "../envs/csvtk.yaml"
    shell:
        "cat {input.csv} | csvtk filter -f 'FDR<{params.fdr}' | csvtk cut -f 'geneId,geneName' | csvtk del-header > {output.csv}"

rule intervene_result:
    input:
        expand("results/{result}/condition_{{A}}_vs_{{B}}.{{quant}}.filtered.csv", result = ["deseq2", "edger", "limma"])
    output:
        directory("results/intervene/condition_{A}_vs_{B}.{quant}")
    message:
        "[intervene] Compare results between DESeq2, edgeR, and limma for {wildcards.quant}"
    conda:
        "../envs/intervene.yaml"
    shell:
        "intervene upset --input {input} --type list --names deseq2,edger,limma --output {output} --save-overlaps --figsize 4 3"

rule intervene_quant:
    input:
        expand("results/{{result}}/condition_{{A}}_vs_{{B}}.{quant}.filtered.csv", quant = ["rsubread", "tximport"])
    output:
        directory("results/intervene/condition_{A}_vs_{B}.{result}")
    message:
        "[intervene] Compare results between Rsubread and tximport for {wildcards.result}"
    conda:
        "../envs/intervene.yaml"
    shell:
        "intervene upset --input {input} --type list --names rsubread,tximport --output {output} --save-overlaps --figsize 4 3"