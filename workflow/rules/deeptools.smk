# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule bam_coverage:
    input:
        bam = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        wig = "results/deeptools/coverage/{sample}.RPKM.bigWig"
    log:
        out = "results/deeptools/coverage/{sample}.RPKM.out",
        err = "results/deeptools/coverage/{sample}.RPKM.err"
    message:
        "[deepTools] Generate a genome coverage track: {wildcards.sample}"
    conda:
        "../envs/deeptools.yaml"
    threads:
        4
    shell:
        "bamCoverage --bam {input.bam} --outFileName {output.wig} --numberOfProcessors {threads} --normalizeUsing RPKM --skipNonCoveredRegions 1> {log.out} 2> {log.err}"
