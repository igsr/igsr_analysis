"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class RnaseqcCounts(Singularity):
    """
    Runs the Rnaseqc Counts task in the broadinstitute/gtex_rnaseq singularity image.
    """

    CMD_KWARGS = ['gatk_flags']
    CMD_ARGS = ['bam_file', 'genes_gtf', 'genome_fasta', 'memory']
    CMD = (
        "python3 /src/run_rnaseqc.py {bam_file} {genes_gtf} {genome_fasta} {PREFIX} "
        "--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java "
        "--memory {memory} --rnaseqc_flags noDoC strictMode {ARGS}")




