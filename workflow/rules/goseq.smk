# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule goseq:
    input:
        csv = "results/{result}/results_{A}_vs_{B}.{type}.csv"
    output:
        csv = "results/{result}/goseq_{A}_vs_{B}.{type}.csv"
    log:
        out = "results/{result}/goseq_{A}_vs_{B}.{type}.out",
        err = "results/{result}/goseq_{A}_vs_{B}.{type}.err"
    message:
        "[bioconductor-goseq]"
    conda:
        "../envs/bioconductor-goseq.yaml"
    script:
        "../scripts/goseq.R"
