# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule sambamba_markdup:
    input:
        bam = "results/star/align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam.bai"
    log:
        out = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.out",
        err = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.err"
    message:
        "[sambamba] Mark duplicate reads in BAM file: {input}"
    threads:
        4
    conda:
        "../envs/sambamba.yaml"
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.bam} 1> {log.out} 2> {log.err}"
