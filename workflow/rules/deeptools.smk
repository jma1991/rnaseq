# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule bam_coverage:
    input:
        bam = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam",
        bai = "results/sambamba/{sample}/Aligned.sortedByCoord.out.markdup.bam.bai"
    output:
        wig = "results/deeptools/{sample}.RPKM.bigWig"
    log:
        "results/deeptools/{sample}.RPKM.log"
    message:
        "[deepTools] Generate a genome coverage track: {wildcards.sample}"
    conda:
        "../envs/deeptools.yaml"
    threads:
        4
    shell:
        "bamCoverage --bam {input.bam} --outFileName {output.wig} --numberOfProcessors {threads} --normalizeUsing RPKM --skipNonCoveredRegions"
