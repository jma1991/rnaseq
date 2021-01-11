# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule multiqc:
    input:
        dir = [
            "results/cutadapt",
            "results/fastqc",
            "results/kallisto",
            "results/star",
            "results/deeptools",
            "results/qualimap",
            "results/rseqc"
        ]
    output:
        ext = multiext("results/multiqc/multiqc_", "report.html", "data.zip")
    log:
        out = "results/multiqc/multiqc.out",
        err = "results/multiqc/multiqc.err"
    params:
        out = "results/multiqc"
    message:
        "[multiqc] Aggregate results into a single report"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc -o {params.out} -z {input.dir} 1> {log.out} 2> {log.err}"
