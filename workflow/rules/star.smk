# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule star_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.fa",
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        dir = directory("results/star/index/{genome}")
    log:
        out = "results/star/index/{genome}/STAR.out",
        err = "results/star/index/{genome}/STAR.err"
    message:
        "[STAR] Creating a STAR index: {wildcards.genome}"
    threads:
        16
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fas} --sjdbGTFfile {input.gtf} 1> {log.out} 2> {log.err}"

rule star_align:
    input:
        idx = expand("results/star/index/{genome}", genome = config["genome"]),
        fqz = star_align_fastq
    output:
        bam = "results/star/align/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        out = "results/star/align/{sample}/STAR.out",
        err = "results/star/align/{sample}/STAR.err"
    params:
        fqz = star_param_fastq,
        out = "results/star/align/{sample}/"
    message:
        "[STAR] Align single-end library to genome: {wildcards.sample}"
    threads:
        16
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {params.fqz} --readFilesCommand gunzip -c --outFileNamePrefix {params.out} --outSAMtype BAM SortedByCoordinate --outTmpDir /tmp/TMPDIR/{wildcards.sample} 1> {log.out} 2> {log.err}"
