# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule sambamba_markdup:
    input:
        bam = "results/star/align/{sample_name}/Aligned.sortedByCoord.out.bam"
    output:
        bam = "results/sambamba/{sample_name}/Aligned.sortedByCoord.out.markdup.bam",
        bai = "results/sambamba/{sample_name}/Aligned.sortedByCoord.out.markdup.bam.bai"
    log:
        out = "results/sambamba/{sample_name}/Aligned.sortedByCoord.out.markdup.out",
        err = "results/sambamba/{sample_name}/Aligned.sortedByCoord.out.markdup.err"
    message:
        "[sambamba] Mark duplicate reads in BAM file: {input}"
    threads:
        4
    conda:
        "../envs/sambamba.yaml"
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.bam} 1> {log.out} 2> {log.err}"
