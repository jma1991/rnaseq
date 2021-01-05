# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

import pandas as pd

class Project:
    
    version = '0.1'

    def __init__(self, config):
        self.config = config
        self.samples = pd.read_csv(config["samples"])
        self.reads = pd.read_csv(config["reads"])

    ## fastq

    def fastq(self):
        ext = [
            "results/fastq/{SM}/{ID}{FQ}.fastq.gz"
        ]
        out = []
        for row in self.reads.itertuples():
            if pd.isna(row.R2):
                exp = expand(ext, SM = row.SM, ID = row.ID, FQ = "")
                out.extend(exp)
            else:
                exp = expand(ext, SM = row.SM, ID = row.ID, FQ = ["_1", "_2"])
                out.extend(exp)
        return out

    def fastq_output(self):
        return [
            self.fastq()
        ]

    ## fastqc

    def fastqc(self):
        ext = [
            "results/fastqc/{SM}/{ID}{FQ}_fastqc.html",
            "results/fastqc/{SM}/{ID}{FQ}_fastqc.zip"
        ]
        out = []
        for row in self.reads.itertuples():
            if pd.isna(row.R2):
                exp = expand(ext, SM = row.SM, ID = row.ID, FQ = "")
                out.extend(exp)
            else:
                exp = expand(ext, SM = row.SM, ID = row.ID, FQ = ["_1", "_2"])
                out.extend(exp)
        return out
        
    def fastqc_output(self):
        return [
            self.fastqc()
        ]

    ## cutadapt

    def cutadapt_single(self):
        ext = [
            "results/cutadapt/{SM}/{ID}.fastq.gz"
        ]
        ind = self.reads["R2"].isnull()
        dat = self.reads[ind]
        return expand(ext, zip, SM = dat["SM"], ID = dat["ID"])

    def cutadapt_paired(self):
        ext = [
            "results/cutadapt/{SM}/{ID}_1.fastq.gz",
            "results/cutadapt/{SM}/{ID}_2.fastq.gz"
        ]
        ind = self.reads["R2"].notnull()
        dat = self.reads[ind]
        return expand(ext, zip, SM = dat["SM"], ID = dat["ID"])

    def cutadapt_output(self):
        return [
            self.cutadapt_single(),
            self.cutadapt_paired()
        ]
    
    ## picard

    def picard_collectrnaseqmetrics(self):
        return expand("results/picard/collectrnaseqmetrics/{SM}.txt", SM = self.samples)

    def picard_output(self):
        return self.picard_collectrnaseqmetrics()
    
    ## rseqc

    def rseqc_bam_stat(self):
        return expand("results/rseqc/bam_stat/{SM}", zip, SM = self.samples["SM"])
    
    def rseqc_inner_distance(self):
        return expand("results/rseqc/inner_distance/{SM}/{LB}", zip, SM = self.samples["SM"], LB = self.samples["LB"])
    
    def rseqc_infer_experiment(self):
        return expand("results/rseqc/infer_experiment/{SM}/{LB}", zip, SM = self.samples["SM"], LB = self.samples["LB"])

    def rseqc_junction_annotation(self):
        return expand("results/rseqc/junction_annotation/{SM}/{LB}", zip, SM = self.samples["SM"], LB = self.samples["LB"])
    
    def rseqc_junction_saturation(self):
        return expand("results/rseqc/junction_saturation/{SM}/{LB}", zip, SM = self.samples["SM"], LB = self.samples["LB"])
    
    def rseqc_read_distribution(self):
        return expand("results/rseqc/read_distribution/{SM}/{LB}", zip, SM = self.samples["SM"], LB = self.samples["LB"])
    
    def rseqc_read_duplication(self):
        return expand("results/rseqc/read_duplication/{SM}/{LB}", zip, SM = self.samples["SM"], LB = self.samples["LB"])

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
        return expand("results/qualimap/rnaseq/{SM}/{LB}", SM = self.samples)
    
    def qualimap_output(self):
        return [
            self.qualimap_rnaseq()
        ]

    ## sambamba

    def sambamba_markdup(self):
        return expand("results/sambamba/markdup/{SM}/{LB}/Aligned.sortedByCoord.out.markdup.bam", SM = self.samples)
    
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

    def gffread_annotation(self):
        return expand("results/gffread/{genome}/{genome}.annotation.tsv", genome = self.config["genome"])

    def gffread_tx2gene(self):
        return expand("results/gffread/{genome}/{genome}.tx2gene.tsv", genome = self.config["genome"])

    def gffread_id2name(self):
        return expand("results/gffread/{genome}/{genome}.id2name.tsv", genome = self.config["genome"])

    def gffread_output(self):
        return [
            self.gffread_transcripts(),
            self.gffread_annotation()
        ]

    ## kallisto

    def kallisto_index(self):
        ext = [
            "results/kallisto/index/{genome}/{genome}.idx"
        ]
        return expand(ext, genome = self.config["genome"])

    def kallisto_quant(self):
        ext = [
            "results/kallisto/quant/{SM}/abundance.h5",
            "results/kallisto/quant/{SM}/abundance.tsv",
            "results/kallisto/quant/{SM}/run_info.json"
        ]
        return expand(ext, zip, SM = self.samples["SM"])

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
            "results/star/align/{SM}/{LB}/Aligned.sortedByCoord.out.bam",
            "results/star/align/{SM}/{LB}/Log.final.out",
            "results/star/align/{SM}/{LB}/Log.out",
            "results/star/align/{SM}/{LB}/Log.progress.out",
            "results/star/align/{SM}/{LB}/SJ.out.tab"
        ]
        return expand(ext, zip, SM = self.samples["SM"], LB = self.samples["LB"])

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
        return expand("results/deeptools/coverage/{SM}.bigWig", SM = self.samples)
    
    def deeptools_profile(self):
        return expand("results/deeptools/profile/genes.pdf")

    def deeptools_output(self):
        return [
            self.deeptools_coverage(),
            self.deeptools_profile()
        ]
    
    ## preseq

    def preseq_c_curve(self):
        return expand("results/preseq/c_curve/{SM}.c_curve.txt", SM = self.samples)
    
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
