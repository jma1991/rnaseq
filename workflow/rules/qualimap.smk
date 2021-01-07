# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule qualimap_rnaseq:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = config["txome"]["genome"])
    output:
        dir = directory("results/qualimap/rnaseq/{sample}")
    params:
        arg = qualimap_rnaseq_params
    log:
        out = "results/qualimap/rnaseq/{sample}/qualimap_rnaseq.out",
        err = "results/qualimap/rnaseq/{sample}/qualimap_rnaseq.err"
    message:
        "[Qualimap] Perform RNA-seq QC analysis: {input.bam}"
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap rnaseq -bam {input.bam} -gtf {input.gtf} -outdir {output.dir} {params.arg} 1> {log.out} 2> {log.err}"
