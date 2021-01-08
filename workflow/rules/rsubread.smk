# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule rsubread:
    input:
        bam = expand("results/star/align/{sample}/Aligned.sortedByCoord.out.bam", sample = project.samples["sample"]),
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = config["genome"])
    output:
        rds = "results/rsubread/object.rds"
    log:
        out = "results/rsubread/object.out",
        err = "results/rsubread/object.err"
    message:
        "[Rsubread] Counting reads to genomic features"
    threads:
        16
    conda:
        "../envs/rsubread.yaml"
    script:
        "../scripts/rsubread.R"
