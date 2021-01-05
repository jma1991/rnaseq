# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule rsubread:
    input:
        bam = expand("results/star/align/{sample_name}/Aligned.sortedByCoord.out.bam", sample_name = pep.sample_table.index),
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = pep.sample_table["genome"].unique())
    output:
        rds = "results/rsubread/counts.rds"
    log:
        out = "results/rsubread/counts.out",
        err = "results/rsubread/counts.err"
    message:
        "[Rsubread] Counting reads to genomic features"
    threads:
        16
    conda:
        "../envs/rsubread.yaml"
    script:
        "../scripts/rsubread.R"
