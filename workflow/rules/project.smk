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
        return expand("results/rseqc/bam_stat/{sample}.bam_stat.txt", zip, sample = self.samples["sample"])
    
    def rseqc_inner_distance(self):
        ext = [
            "results/rseqc/inner_distance/{sample}.inner_distance.txt",
            "results/rseqc/inner_distance/{sample}.inner_distance_freq.txt",
            "results/rseqc/inner_distance/{sample}.inner_distance_plot.pdf",
            "results/rseqc/inner_distance/{sample}.inner_distance_plot.r"
        ]
        return expand(ext, sample = self.samples["sample"])

    def rseqc_infer_experiment(self):
        ext = [
            "results/rseqc/infer_experiment/{sample}.infer_experiment.txt"
        ]
        return expand(ext, sample = self.samples["sample"])

    def rseqc_junction_annotation(self):
        ext = [
            "results/rseqc/junction_annotation/{sample}.junction.Interact.bed",
            "results/rseqc/junction_annotation/{sample}.junction.bed",
            "results/rseqc/junction_annotation/{sample}.junction.xls",
            "results/rseqc/junction_annotation/{sample}.junction_plot.r",
            "results/rseqc/junction_annotation/{sample}.splice_events.pdf",
            "results/rseqc/junction_annotation/{sample}.splice_junction.pdf"
        ]
        return expand(ext, sample = self.samples["sample"])

    def rseqc_junction_saturation(self):
        ext = [
            "results/rseqc/junction_saturation/{sample}.junctionSaturation_plot.pdf",
            "results/rseqc/junction_saturation/{sample}.junctionSaturation_plot.r"
        ]
        return expand(ext, sample = self.samples["sample"])

    def rseqc_read_distribution(self):
        ext = [
            "results/rseqc/read_distribution/{sample}.read_distribution.txt"
        ]
        return expand(ext, sample = self.samples["sample"])

    def rseqc_read_duplication(self):
        ext = [
            "results/rseqc/read_duplication/{sample}.DupRate_plot.pdf",
            "results/rseqc/read_duplication/{sample}.DupRate_plot.r",
            "results/rseqc/read_duplication/{sample}.pos.DupRate.xls",
            "results/rseqc/read_duplication/{sample}.seq.DupRate.xls"
        ]
        return expand(ext, sample = self.samples["sample"])

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
        return expand("results/qualimap/rnaseq/{sample}", sample = self.samples["sample"])
    
    def qualimap_output(self):
        return [
            self.qualimap_rnaseq()
        ]

    ## sambamba

    def sambamba_markdup(self):
        return expand("results/sambamba/markdup/{sample}/Aligned.sortedByCoord.out.bam", sample = self.samples["sample"])
    
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
        return expand(ext, sample = self.samples["sample"])

    def kallisto_output(self):
        return [
            self.kallisto_index(),
            self.kallisto_quant()
        ]

    ## star

    def star_index(self):
        ext = [
            "results/star/index/{genome}"
        ]
        return expand(ext, genome = self.config["genome"])

    def star_align(self):
        ext = [
            "results/star/align/{sample}/Aligned.sortedByCoord.out.bam",
        ]
        return expand(ext, sample = self.samples["sample"])

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
        return expand("results/deseq2/object.{type}.rds", type = "tximport")

    def deseq2_counts(self):
        return expand("results/deseq2/counts.{type}.csv", type = "tximport")

    def deseq2_normcounts(self):
        return expand("results/deseq2/normcounts.{type}.csv", type = "tximport")

    def deseq2_logcounts(self):
        return expand("results/deseq2/logcounts.{type}.csv", type = "tximport")

    def deseq2_results(self):
        results = []
        contrasts = self.config["contrasts"]
        for contrast, conditions in contrasts.items():
            result = expand("results/deseq2/results_{A}_vs_{B}.{type}.csv", A = conditions["A"], B = conditions["B"], type = "tximport")
            results.extend(result)
        return results

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
        return expand("results/edger/object.{type}.rds", type = "tximport")

    def edger_counts(self):
        return expand("results/edger/counts.{type}.csv", type = "tximport")

    def edger_normcounts(self):
        return expand("results/edger/normcounts.{type}.csv", type = "tximport")

    def edger_logcounts(self):
        return expand("results/edger/logcounts.{type}.csv", type = "tximport")

    def edger_results(self):
        results = []
        contrasts = self.config["contrasts"]
        for contrast, conditions in contrasts.items():
            result = expand("results/edger/results_{A}_vs_{B}.{type}.csv", A = conditions["A"], B = conditions["B"], type = "tximport")
            results.extend(result)
        return results

    def edger_output(self):
        return [
            self.edger_object(),
            self.edger_counts(),
            self.edger_normcounts(),
            self.edger_logcounts(),
            self.edger_results()
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

    ## bioconductor-goseq

    def goseq(self):
        results = []
        contrasts = self.config["contrasts"]
        for contrast, conditions in contrasts.items():
            result = expand("results/{result}/goseq_{A}_vs_{B}.{type}.csv", result = "deseq2", A = conditions["A"], B = conditions["B"], type = "tximport")
            results.extend(result)
        return results

    def goseq_output(self):
        return [
            self.goseq()
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
