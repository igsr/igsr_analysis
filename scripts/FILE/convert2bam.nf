#!/usr/bin/env nextflow

/* 
 * Script to convert a .cram file to a bam file
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.queue = 'production-rh74'
params.threads = 10

files = Channel.fromPath( '/hps/nobackup/production/reseq-info/ernesto/HIGHCOV/converted/*.cram' )

process convert2bam {
        /*
        Process to run samtools view to convert the .cram file to .bam format
        It will also create an index for the converted BAM file

        Returns
        -------
        Path to a BAM file
        */
	tag { cram_f }

        publishDir "converted_bam", mode: 'copy', overwrite: true

        memory '15 GB'
        executor 'lsf'
	maxForks 100
        queue "${params.queue}"
        cpus "${params.threads}"
        errorStrategy 'ignore'

        input:
        file cram_f from files

        output:
        file "${cram_f}.bam" into out_bam

        """
        samtools view -b -o ${cram_f}.bam ${cram_f} --threads ${params.threads}
        """
}
