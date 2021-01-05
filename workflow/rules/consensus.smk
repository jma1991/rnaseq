# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule consensusDE:
    input:
        csv = expand("results/{result}/condition_{{A}}_vs_{{B}}.{{quant}}.csv", result = ["deseq2", "edger", "limma"])
    output:
        csv = "results/consensus/condition_{A}_vs_{B}.{quant}.consensus.csv"
    log:
        out = "results/consensus/condition_{A}_vs_{B}.{quant}.consensus.out",
        err = "results/consensus/condition_{A}_vs_{B}.{quant}.consensus.err"
    message:
        "[intervene] Generate consensus results between DESeq2, edgeR, and limma for {wildcards.quant} quantification"
    script:
        "../scripts/consensus.R"
