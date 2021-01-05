# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

import pandas as pd

wildcard_constraints:
    sample_name = "|".join(pep.sample_table["sample_name"]),
    subsample_name = "|".join(pep.subsample_table["subsample_name"])

def all_single_end(sample_name):
    u = pep.subsample_table.loc[sample_name, ]
    m = [is_single_end(sample_name, x) for x in u["subsample_name"]]
    return all(m)

def all_paired_end(sample_name):
    u = pep.subsample_table.loc[sample_name, ]
    m = [is_paired_end(sample_name, x) for x in u["subsample_name"]]
    return all(m)

def is_single_end(sample_name, subsample_name):
    obj = pep.subsample_table.loc[(sample_name, subsample_name), "read2"]
    if isinstance(obj, pd.core.series.Series):
        raise ValueError()
    return pd.isnull(obj)

def is_paired_end(sample_name, subsample_name):
    obj = pep.subsample_table.loc[(sample_name, subsample_name), "read2"]
    if isinstance(obj, pd.core.series.Series):
        raise ValueError()
    return pd.notnull(obj)

def get_fastqs(wildcards):
    if is_single_end(wildcards.sample_name, wildcards.subsample_name):
        u = pep.subsample_table.loc[(wildcards.sample_name, wildcards.subsample_name), "read1"]
        return {"fq1": u}
    else:
        u = pep.subsample_table.loc[(wildcards.sample_name, wildcards.subsample_name), ["read1", "read2"]].dropna()
        return {"fq1": u.read1, "fq2": u.read2}

def kallisto_quant_index(wildcards):
    x = pep.sample_table.loc[wildcards.sample_name, "genome"]
    return f"results/kallisto/index/{x}/{x}.idx"

def kallisto_quant_fastq(wildcards):
    if all_single_end(wildcards.sample_name):
        return expand("results/cutadapt/{sample_name}/{subsample_name}.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"])
    elif all_paired_end(wildcards.sample_name):
        return expand("results/cutadapt/{sample_name}/{subsample_name}_{read_name}.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"], read_name = ["1", "2"])
    else:
        raise ValueError()

def kallisto_quant_stranded(sample_name):
    d = {"U": "", "F": "--fr-stranded", "R": "--rf-stranded"}
    k = pep.sample_table.loc[sample_name, "stranded"]
    return d.get(k)

def kallisto_quant_arguments(wildcards):
    str = kallisto_quant_stranded(wildcards.sample_name)
    if all_single_end(wildcards.sample_name):
        return f"--single {str} --fragment-length=200 --sd=20"
    elif all_paired_end(wildcards.sample_name):
        return f"{str}"
    else:
        raise ValueError()

def star_align_index(wildcards):
    x = pep.sample_table.loc[wildcards.sample_name, "genome"]
    return f"results/star/index/{x}"

def star_align_fastq(wildcards):
    if all_single_end(wildcards.sample_name):
        return expand("results/cutadapt/{sample_name}/{subsample_name}.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"])
    elif all_paired_end(wildcards.sample_name):
        return expand("results/cutadapt/{sample_name}/{subsample_name}_{read_name}.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"], read_name = ["1", "2"])
    else:
        raise ValueError()

def star_align_reads(wildcards):
    if all_single_end(wildcards.sample_name):
        fqz = expand("results/cutadapt/{sample_name}/{subsample_name}.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"])
        return ",".join(fqz)
    elif all_paired_end(wildcards.sample_name):
        fq1 = expand("results/cutadapt/{sample_name}/{subsample_name}_1.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"])
        fq2 = expand("results/cutadapt/{sample_name}/{subsample_name}_2.fastq.gz", sample_name = wildcards.sample_name, subsample_name = pep.subsample_table.loc[wildcards.sample_name, "subsample_name"])
        return ",".join(fq1) + " " + ",".join(fq2)
    else:
        raise ValueError()

class Output:

    def __init__(self, sample_table, subsample_table):
        self.sample_table = sample_table
        self.subsample_table = subsample_table
        self.contrasts = {"A": ["MEF_KO", "D3_KO", "D5_KO", "D7_KO", "iPS_KO", "ESC_KO"], "B": ["MEF_WT", "D3_WT", "D5_WT", "D7_WT", "iPS_WT", "ESC_WT"]}

    def cutadapt_single(self):
        out = []
        for sample_name, subsample_name in self.subsample_table.index:
            if is_single_end(sample_name, subsample_name):
                itr = [
                    "results/cutadapt/{sample_name}/{subsample_name}.fastq.gz"
                ]
                exp = expand(itr, sample_name = sample_name, subsample_name = subsample_name)
                out.extend(exp)
        return out

    def cutadapt_paired(self):
        out = []
        for sample_name, subsample_name in self.subsample_table.index:
            if is_paired_end(sample_name, subsample_name):
                itr = [
                    "results/cutadapt/{sample_name}/{subsample_name}_1.fastq.gz",
                    "results/cutadapt/{sample_name}/{subsample_name}_2.fastq.gz"
                ]
                exp = expand(itr, sample_name = sample_name, subsample_name = subsample_name)
                out.extend(exp)
        return out

    def genomepy_install(self):
        itr = [
            "results/genomepy/{genome}/{genome}.annotation.bed",
            "results/genomepy/{genome}/{genome}.annotation.gtf",
            "results/genomepy/{genome}/{genome}.fa",
            "results/genomepy/{genome}/{genome}.fa.fai",
            "results/genomepy/{genome}/{genome}.fa.sizes",
            "results/genomepy/{genome}/{genome}.gaps.bed",
            "results/genomepy/{genome}/README.txt"
        ]
        return expand(itr, genome = self.sample_table["genome"].unique())


    def kallisto_index(self):
        itr = [
            "results/kallisto/index/{genome}/{genome}.idx"
        ]
        return expand(itr, genome = self.sample_table["genome"].unique())

    def kallisto_quant(self):
        itr = [
            "results/kallisto/quant/{sample_name}/abundance.h5",
            "results/kallisto/quant/{sample_name}/abundance.tsv"
        ]
        return expand(itr, sample_name = self.sample_table.index)

    def star_index(self):
        itr = [
            "results/star/index/{genome}",
        ]
        return expand(itr, genome = self.sample_table["genome"].unique())

    def star_align(self):
        itr = [
            "results/star/align/{sample_name}/Aligned.sortedByCoord.out.bam",
            "results/star/align/{sample_name}/Log.final.out",
            "results/star/align/{sample_name}/Log.out",
            "results/star/align/{sample_name}/Log.progress.out",
            "results/star/align/{sample_name}/SJ.out.tab"
        ]
        return expand(itr, sample_name = self.sample_table.index)

    def sambamba_markdup(self):
        return expand("results/sambamba/{sample_name}/Aligned.sortedByCoord.out.markdup.bam", sample_name = self.sample_table.index)

    def rsubread_featurecounts(self):
        return "results/rsubread/counts.rds"

    def deseq2_init(self):
        return expand("results/deseq2/dds.{type}.rds", type = ['rsubread', 'tximport'])

    def deseq2_counts(self):
        return expand("results/deseq2/counts.{type}.csv", type = ['rsubread', 'tximport'])

    def deseq2_normcounts(self):
        return expand("results/deseq2/normcounts.{type}.csv", type = ['rsubread', 'tximport'])

    def deseq2_logcounts(self):
        return expand("results/deseq2/logcounts.{type}.csv", type = ['rsubread', 'tximport'])

    def deseq2_results(self):
        return expand(expand("results/deseq2/condition_{A}_vs_{B}.{{type}}.csv", zip, A = self.contrasts['A'], B = self.contrasts['B']), type = ["rsubread", "tximport"])

    def edger_init(self):
        return expand("results/edger/dge.{type}.rds", type = ['rsubread', 'tximport'])

    def edger_counts(self):
        return expand("results/edger/counts.{type}.csv", type = ['rsubread', 'tximport'])

    def edger_normcounts(self):
        return expand("results/edger/normcounts.{type}.csv", type = ['rsubread', 'tximport'])

    def edger_logcounts(self):
        return expand("results/edger/logcounts.{type}.csv", type = ['rsubread', 'tximport'])

    def edger_results(self):
        return expand(expand("results/edger/condition_{A}_vs_{B}.{{type}}.csv", zip, A = self.contrasts['A'], B = self.contrasts['B']), type = ["rsubread", "tximport"])

    def limma_dge(self):
        return expand("results/limma/dge.{type}.rds", type = ['rsubread', 'tximport'])

    def limma_els(self):
        return expand("results/limma/els.{type}.rds", type = ['rsubread', 'tximport'])

    def limma_counts(self):
        return expand("results/limma/counts.{type}.csv", type = ['rsubread', 'tximport'])

    def limma_normcounts(self):
        return expand("results/limma/normcounts.{type}.csv", type = ['rsubread', 'tximport'])

    def limma_logcounts(self):
        return expand("results/limma/logcounts.{type}.csv", type = ['rsubread', 'tximport'])

    def limma_results(self):
        return expand(expand("results/limma/condition_{A}_vs_{B}.{{type}}.csv", zip, A = self.contrasts['A'], B = self.contrasts['B']), type = ["rsubread", "tximport"])

    def plot_corr(self):
        return expand("results/{result}/corr.{type}.pdf", result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_dist(self):
        return expand("results/{result}/dist.{type}.pdf", result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_prcomp(self):
        return expand("results/{result}/prcomp.{type}.pdf", result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_cmdscale(self):
        return expand("results/{result}/cmdscale.{type}.pdf", result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_diff(self):
        return expand(expand("results/{{result}}/condition_{A}_vs_{B}.{{type}}.diff.pdf", zip, A = self.contrasts['A'], B = self.contrasts['B']), result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_pvalue(self):
        return expand(expand("results/{{result}}/condition_{A}_vs_{B}.{{type}}.pvalue.pdf", zip, A = self.contrasts['A'], B = self.contrasts['B']), result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_volcano(self):
        return expand(expand("results/{{result}}/condition_{A}_vs_{B}.{{type}}.volcano.pdf", zip, A = self.contrasts['A'], B = self.contrasts['B']), result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def plot_heatmap(self):
        return expand(expand("results/{{result}}/condition_{A}_vs_{B}.{{type}}.heatmap.pdf", zip, A = self.contrasts['A'], B = self.contrasts['B']), result = ["deseq2", "edger", "limma"], type = ["rsubread", "tximport"])

    def rseqc_read_distribution(self):
        return expand("results/rseqc/{sample_name}.read_distribution.txt", sample_name = self.sample_table.index)

    def rseqc_gene_body_coverage(self):
        ext = [
            "results/rseqc/{sample_name}.geneBodyCoverage.curves.pdf",
            "results/rseqc/{sample_name}.geneBodyCoverage.r",
            "results/rseqc/{sample_name}.geneBodyCoverage.txt"
        ]
        return expand(ext, sample_name = self.sample_table.index)

    def rseqc_inner_distance(self):
        ext = [
            "results/rseqc/{sample_name}.inner_distance.txt",
            "results/rseqc/{sample_name}.inner_distance_freq.txt",
            "results/rseqc/{sample_name}.inner_distance_plot.pdf",
            "results/rseqc/{sample_name}.inner_distance_plot.r"
        ]
        return expand(ext, sample_name = self.sample_table.index)

    def rseqc_read_gc(self):
        ext = [
            "results/rseqc/{sample_name}.GC.xls",
            "results/rseqc/{sample_name}.GC_plot.pdf",
            "results/rseqc/{sample_name}.GC_plot.r"
        ]
        return expand(ext, sample_name = self.sample_table.index)

    def rseqc_read_duplication(self):
        ext = [
            "results/rseqc/{sample_name}.DupRate_plot.pdf",
            "results/rseqc/{sample_name}.DupRate_plot.r",
            "results/rseqc/{sample_name}.pos.DupRate.xls",
            "results/rseqc/{sample_name}.seq.DupRate.xls"
        ]
        return expand(ext, sample_name = self.sample_table.index)

    def rseqc_junction_annotation(self):
        ext = [
            "results/rseqc/{sample_name}.junction.Interact.bed",
            "results/rseqc/{sample_name}.junction.bed",
            "results/rseqc/{sample_name}.junction.xls",
            "results/rseqc/{sample_name}.junction_plot.r",
            "results/rseqc/{sample_name}.splice_events.pdf",
            "results/rseqc/{sample_name}.splice_junction.pdf"
        ]
        return expand(ext, sample_name = self.sample_table.index)

    def rseqc_junction_saturation(self):
        ext = [
            "results/rseqc/{sample_name}.junctionSaturation_plot.pdf",
            "results/rseqc/{sample_name}.junctionSaturation_plot.r"
        ]
        return expand(ext, sample_name = self.sample_table.index)

    def rseqc_infer_experiment(self):
        return expand("results/rseqc/{sample_name}.infer_experiment.txt", sample_name = self.sample_table.index)

    def rseqc_bam_stat(self):
        return expand("results/rseqc/{sample_name}.bam_stat.txt", sample_name = self.sample_table.index)

    def bam_coverage(self):
        return expand("results/deeptools/{sample_name}.RPKM.bigWig", sample_name = self.sample_table.index)
