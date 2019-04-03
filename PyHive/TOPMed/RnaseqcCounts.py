"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class RnaseqcCounts(Singularity):
    """
    Runs the Rnaseqc Counts task in the broadinstitute/gtex_rnaseq singularity image.
    """
    PIPELINE = 'rnaseqccounts'
    CMD_KWARGS = ['gatk_flags']
    CMD_ARGS = ['md_bam_file', 'genes_gtf', 'genome_fasta', 'memory']
    CMD = (
        "python3 /src/run_rnaseqc.py {md_bam_file} {genes_gtf} {genome_fasta} {PREFIX} "
        "--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java "
        "-o {WORKING_DIR} --memory {memory} --rnaseqc_flags noDoC strictMode {ARGS}")

    FILES = {
        'gene_rpkm': "{PREFIX}.gene_rpkm.gct.gz",
        'gene_counts': "{PREFIX}.gene_reads.gct.gz",
        'exon_counts': "{PREFIX}.exon_reads.gct.gz",
        'count_metrics': "{PREFIX}.metrics.tsv",
        'count_outputs': "{PREFIX}.tar.gz"
    }




