# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule strintie_assemble:
    input:
        bam = "results/hisat2/{sample}.bam",
        gtf = "results/genomepy/GRCh38.p13/GRCh38.p13.annotation.gtf"
    output:
        gtf = "results/stringtie/{sample}/transcripts.gtf"
    message:
        "[Stringtie] Assemble transcripts for each sample:"
    shell:
        "stringtie -p {threads} -G {input.gtf} -o {output.gtf} {input.bam}"

rule stringtie_merge:
    input:
        ann = "results/genomepy/GRCh38.p13/GRCh38.p13.annotation.gtf",
        gtf = expand("results/stringtie/{sample}/transcripts.gtf", sample = samples["sample"])
    output:
        gtf = "results/stringtie/merged.gtf"
    message:
        "[Stringtie] Merge transcripts from all samples:"
    shell:
        "stringtie --merge -p {threads} -G {input.ann} -o {output.gtf} {input.gtf}"

rule stringtie_abundance:
    input:
        bam = "results/hisat2/{sample}.bam",
        gtf = "results/stringtie/merged.gtf"
    output:
        gtf = "results/stringtie/{sample}/abundance.gtf"
    message:
        "[Stringtie] Estimate transcript abundances and create table counts for Ballgown:"
    threads:
        4
    shell:
        "stringtie -e -B -p {threads} -G {input.gtf} -o {output.gtf} {input.bam}"
