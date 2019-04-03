"""
Created on 21 March 2019

@author: galdam
"""


from Utils.RunSingularity import Singularity


class Star(Singularity):
    """
    Runs the Star alignment task in the broadinstitute/gtex_rnaseq singularity image.
    """

    PIPELINE = 'Star'
    CMD_KWARGS = [
        "outFilterMultimapNmax", "alignSJoverhangMin", "alignSJDBoverhangMin", "outFilterMismatchNmax",
        "outFilterMismatchNoverLmax", "alignIntronMin", "alignIntronMax", "alignMatesGapMax", "outFilterType",
        "outFilterScoreMinOverLread", "outFilterMatchNminOverLread", "limitSjdbInsertNsj", "outSAMstrandField",
        "outFilterIntronMotifs", "alignSoftClipAtReferenceEnds", "quantMode", "outSAMattrRGline", "outSAMattributes",
        "chimSegmentMin", "chimJunctionOverhangMin", "chimOutType", "chimMainSegmentMultNmax", "sjdbFileChrStartEnd",
    ]
    CMD_ARGS = ['star_index', 'num_threads']
    CMD = ("/src/run_STAR.py {star_index} {fastq1} {fastq2} {PREFIX} "
           "--output_dir {WORKING_DIR} --threads {num_threads} {ARGS}")
    FILES = {
        'bam_file': "{PREFIX}.Aligned.sortedByCoord.out.bam",
        'transcriptome_bam': "{PREFIX}.Aligned.toTranscriptome.out.bam",
        'chimeric_junctions': "{PREFIX}.Chimeric.out.junction",
        'chimeric_bam_file': "{PREFIX}.Chimeric.out.sorted.bam",
        'read_counts': "{PREFIX}.ReadsPerGene.out.tab",
        'junctions': "{PREFIX}.SJ.out.tab",
        'junctions_pass1': "{PREFIX}._STARpass1/SJ.out.tab"
        #logs = ["star_out/${prefix}.Log.final.out", "star_out/${prefix}.Log.out", "star_out/${prefix}.Log.progress.out"]
    }

    def get_cmd_args(self):
        """
        Extending the function to unpack the fastq files
        :return:
        """
        # Call super to unpack the majority of the arguments
        options_dict = super().get_cmd_args()

        # Include the fastq files
        fastq1, fastq2 = sorted(self.param_required('fastq'))
        options_dict['fastq1'] = fastq1
        options_dict['fastq2'] = fastq2
        return options_dict
