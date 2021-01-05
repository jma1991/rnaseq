# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

import pandas as pd

class Project:

    def __init__(self, config):
        self.config = config
        self.samples = pd.read_csv(config["samples"])
        self.units = pd.read_csv(config["units"])

    def fastq(self):
        ext = [
            "results/fastq/{sample}/{unit}{read}.fastq.gz"
        ]
        out = []
        for row in self.units.itertuples():
            if pd.isna(row.read2):
                exp = expand(ext, sample = row.sample, unit = row.unit, read = "")
                out.extend(exp)
            else:
                exp = expand(ext, sample = row.sample, unit = row.unit, read = ["_1", "_2"])
                out.extend(exp)
        return out

    def fastq_output(self):
        return [
            self.fastq()
        ]

    ## fastqc

    def fastqc(self):
        ext = [
            "results/fastqc/{sample}/{unit}{read}_fastqc.html",
            "results/fastqc/{sample}/{unit}{read}_fastqc.zip"
        ]
        out = []
        for row in self.units.itertuples():
            if pd.isna(row.read2):
                exp = expand(ext, sample = row.sample, unit = row.unit, read = "")
                out.extend(exp)
            else:
                exp = expand(ext, sample = row.sample, unit = row.unit, read = ["_1", "_2"])
                out.extend(exp)
        return out
        
    def fastqc_output(self):
        return [
            self.fastqc()
        ]

    ## cutadapt

    def cutadapt_single(self):
        ext = [
            "results/cutadapt/{sample}/{unit}.fastq.gz"
        ]
        ind = self.units["read2"].isnull()
        dat = self.units[ind]
        return expand(ext, zip, sample = dat["sample"], unit = dat["unit"])

    def cutadapt_paired(self):
        ext = [
            "results/cutadapt/{sample}/{unit}_1.fastq.gz",
            "results/cutadapt/{sample}/{unit}_2.fastq.gz"
        ]
        ind = self.units["read2"].notnull()
        dat = self.units[ind]
        return expand(ext, zip, sample = dat["sample"], unit = dat["unit"])

    def cutadapt_output(self):
        return [
            self.cutadapt_single(),
            self.cutadapt_paired()
        ]
    
    ## picard

    def picard_collectrnaseqmetrics(self):
        return expand("results/picard/collectrnaseqmetrics/{sample}.txt", sample = self.samples)

    def picard_output(self):
        return self.picard_collectrnaseqmetrics()
    
    ## rseqc

    def rseqc_bam_stat(self):
        return expand("results/rseqc/bam_stat/{sample}", zip, sample = self.samples["sample"])
    
    def rseqc_inner_distance(self):
        return expand("results/rseqc/inner_distance/{sample}/{LB}", zip, sample = self.samples["sample"], LB = self.samples["LB"])
    
    def rseqc_infer_experiment(self):
        return expand("results/rseqc/infer_experiment/{sample}/{LB}", zip, sample = self.samples["sample"], LB = self.samples["LB"])

    def rseqc_junction_annotation(self):
        return expand("results/rseqc/junction_annotation/{sample}/{LB}", zip, sample = self.samples["sample"], LB = self.samples["LB"])
    
    def rseqc_junction_saturation(self):
        return expand("results/rseqc/junction_saturation/{sample}/{LB}", zip, sample = self.samples["sample"], LB = self.samples["LB"])
    
    def rseqc_read_distribution(self):
        return expand("results/rseqc/read_distribution/{sample}/{LB}", zip, sample = self.samples["sample"], LB = self.samples["LB"])
    
    def rseqc_read_duplication(self):
        return expand("results/rseqc/read_duplication/{sample}/{LB}", zip, sample = self.samples["sample"], LB = self.samples["LB"])

    def rseqc_output(self):
        return [
            self.rseqc_bam_stat(),
            self.rseqc_inner_distance(),
            self.rseqc_infer_experiment(),
            self.rseqc_junction_annotation(),
            self.rseqc_junction_saturation(),
            self.rseqc_read_distribution(),
            self.rseqc_read_duplication()
        ]

    ## qualimap

    def qualimap_rnaseq(self):
        return expand("results/qualimap/rnaseq/{sample}/{LB}", sample = self.samples)
    
    def qualimap_output(self):
        return [
            self.qualimap_rnaseq()
        ]

    ## sambamba

    def sambamba_markdup(self):
        return expand("results/sambamba/markdup/{sample}/{LB}/Aligned.sortedByCoord.out.markdup.bam", sample = self.samples)
    
    def sambamba_output(self):
        return [
            self.sambamba_markdup()
        ]

    ## genomepy
    
    def genomepy_install(self):
        ext = [
            "results/genomepy/{genome}/{genome}.annotation.bed.gz",
            "results/genomepy/{genome}/{genome}.annotation.gtf.gz",
            "results/genomepy/{genome}/{genome}.fa",
            "results/genomepy/{genome}/{genome}.fa.fai",
            "results/genomepy/{genome}/{genome}.fa.sizes",
            "results/genomepy/{genome}/{genome}.gaps.bed",
            "results/genomepy/{genome}/README.txt"
        ]
        return expand(ext, genome = self.config["genome"])

    def genomepy_gunzip(self):
        ext = [
            "results/genomepy/{genome}/{genome}.annotation.bed",
            "results/genomepy/{genome}/{genome}.annotation.gtf"
        ]
        return expand(ext, genome = self.config["genome"])
    
    def genomepy_output(self):
        return [
            self.genomepy_install(),
            self.genomepy_gunzip()
        ]

    ## gffread

    def gffread_transcripts(self):
        return expand("results/gffread/{genome}/{genome}.transcripts.fa", genome = self.config["genome"])

    def gffread_tx2gene(self):
        return expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = self.config["genome"])

    def gffread_id2name(self):
        return expand("results/gffread/{genome}/{genome}.id2name.tsv", genome = self.config["genome"])

    def gffread_output(self):
        return [
            self.gffread_transcripts(),
            self.gffread_tx2gene(),
            self.gffread_id2name()
        ]

    ## kallisto

    def kallisto_index(self):
        ext = [
            "results/kallisto/index/{genome}/{genome}.idx"
        ]
        return expand(ext, genome = self.config["genome"])

    def kallisto_quant(self):
        ext = [
            "results/kallisto/quant/{sample}/abundance.h5",
            "results/kallisto/quant/{sample}/abundance.tsv",
            "results/kallisto/quant/{sample}/run_info.json"
        ]
        return expand(ext, zip, sample = self.samples["sample"])

    def kallisto_output(self):
        return [
            self.kallisto_index(),
            self.kallisto_quant()
        ]

    ## star

    def star_index(self):
        itr = [
            "results/star/index/{genome}"
        ]
        return expand(itr, genome = self.genome)

    def star_align(self):
        ext = [
            "results/star/align/{sample}/{LB}/Aligned.sortedByCoord.out.bam",
            "results/star/align/{sample}/{LB}/Log.final.out",
            "results/star/align/{sample}/{LB}/Log.out",
            "results/star/align/{sample}/{LB}/Log.progress.out",
            "results/star/align/{sample}/{LB}/SJ.out.tab"
        ]
        return expand(ext, zip, sample = self.samples["sample"], LB = self.samples["LB"])

    def star_output(self):
        return [
            self.star_index(),
            self.star_align()
        ]

    ## bioconductor-tximport

    def tximport_object(self):
        return "results/tximport/object.rds"

    def tximport_output(self):
        return [
            self.tximport_object()
        ]

    ## bioconductor-rsubread

    def rsubread_object(self):
        return "results/rsubread/object.rds"
    
    def rsubread_output(self):
        return [
            self.rsubread_object()
        ]

    ## bioconductor-deseq2

    def deseq2_object(self):
        return "results/deseq2/object.rds"

    def deseq2_counts(self):
        return "results/deseq2/counts.csv"

    def deseq2_normcounts(self):
        return "results/deseq2/normcounts.csv"

    def deseq2_logcounts(self):
        return "results/deseq2/logcounts.csv"
    
    def deseq2_results(self):
        return expand("results/deseq2/results_{A}_vs_{B}.csv")

    def deseq2_output(self):
        return [
            self.deseq2_object(),
            self.deseq2_counts(),
            self.deseq2_normcounts(),
            self.deseq2_logcounts(),
            self.deseq2_results()
        ]

    ## bioconductor-edger

    def edger_object(self):
        return "results/edger/object.rds"

    def edger_counts(self):
        return "results/edger/counts.csv"

    def edger_normcounts(self):
        return "results/edger/normcounts.csv"

    def edger_logcounts(self):
        return "results/edger/logcounts.csv"

    def edger_output(self):
        return [
            self.edger_object(),
            self.edger_counts(),
            self.edger_normcounts(),
            self.edger_logcounts()
        ]

    ## bioconductor-limma

    def limma_object(self):
        return "results/limma/object.rds"

    def limma_counts(self):
        return "results/limma/counts.csv"

    def limma_normcounts(self):
        return "results/limma/normcounts.csv"

    def limma_logcounts(self):
        return "results/limma/logcounts.csv"

    def limma_output(self):
        return [
            self.limma_object(),
            self.limma_counts(),
            self.limma_normcounts(),
            self.limma_logcounts()
        ]

    ## deeptools

    def deeptools_coverage(self):
        return expand("results/deeptools/coverage/{sample}.bigWig", sample = self.samples)
    
    def deeptools_profile(self):
        return expand("results/deeptools/profile/genes.pdf")

    def deeptools_output(self):
        return [
            self.deeptools_coverage(),
            self.deeptools_profile()
        ]
    
    ## preseq

    def preseq_c_curve(self):
        return expand("results/preseq/c_curve/{sample}.c_curve.txt", sample = self.samples)
    
    def preseq_output(self):
        return self.preseq_c_curve()

    ## multiqc

    def multiqc(self):
        return [
            "results/multiqc/multiqc_data",
            "results/multiqc/multiqc_report.html"
        ]

    def multiqc_output(self):
        return [
            self.multiqc()
        ]
