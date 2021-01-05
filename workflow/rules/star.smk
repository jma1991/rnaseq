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
    message:
        "[STAR] Creating a STAR index: {wildcards.genome}"
    threads:
        16
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fas} --sjdbGTFfile {input.gtf} --sjdbOverhang 74"

rule star_align:
    input:
        idx = star_align_index,
        fqz = star_align_fastq
    output:
        bam = "results/star/align/{sample_name}/Aligned.sortedByCoord.out.bam"
    params:
        fqz = star_align_reads,
        out = "results/star/align/{sample_name}/"
    message:
        "[STAR] Align single-end library to genome: {wildcards.sample_name}"
    threads:
        16
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {params.fqz} --readFilesCommand gunzip -c --outFileNamePrefix {params.out} --outSAMtype BAM SortedByCoordinate --outTmpDir /tmp/TMPDIR/{wildcards.sample_name}"
