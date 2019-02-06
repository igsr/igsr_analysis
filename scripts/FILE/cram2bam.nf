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
// params defaults for ascp client
params.key_file = '/homes/ernesto/.aspera/connect/etc/asperaweb_id_dsa.openssh' // Private-key file name (id_rsa) for authentication
params.transfer_rate = '900M'
params.port = 33001 // TCP port used for SSH authentication

//print usage
if (params.help) {
    log.info ''
    log.info 'Script to download .cram file/s from the archive and convert them to .bam'
    log.info '-------------------------------------------------------------------------'
    log.info ''
    log.info 'This script uses the ASPERA file transfer service for the download.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow cram2bam.nf --file input.txt'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--file FILE  File with the urls pointing to the .cram files to be converted.'
    log.info '    This file needs to have the following format:'
    log.info '    url,dest,prefix'
    log.info '    era-fasp@fasp.sra.ebi.ac.uk:/vol1/ERZ454/ERZ454001/ERR1457180.cram,/nfs/production/reseq-info/work/ernesto/isgr/SCRATCH/31_01_2019/ERR1457180.cram,ERR1457180'
    log.info ''
    exit 1
}

Channel
    .fromPath(params.file)
    .splitCsv(header:true)
    .map{ row-> tuple(row.url, file(row.dest), row.prefix) }
    .set { paths_ch }

process downloadFile {
	/*
	This process uses ASPERA for downloading the file

	Returns
	-------
	Path to a CRAM file downloaded from the repo
	*/

	memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        set url, file(dest),prefix from paths_ch

	output:
	file(dest) into cramf
	val prefix into prefix1

	"""	
	ascp -i ${params.key_file} -Tr -Q -l ${params.transfer_rate} -P${params.port} -L- ${url} ${dest}
        """
}


process convert2bam {
	/*
	Process to run samtools view to convert the .cram file to .bam format
	It will also create an index for the converted BAM file

	Returns
	-------
	Path to a BAM file
	*/

	publishDir "converted", mode: 'copy', overwrite: true

	memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

	input:
	file f from cramf
	val prefix from prefix1

	output:
	file "${prefix}.bam" into out_bam
	val prefix into prefix2

	"""
	${params.samtools_folder}/samtools view -b -o ${prefix}.bam ${f}
	"""
}

process make_index {
	/*
	Process to create a samtools index on the converted BAM file
	*/

	publishDir "converted", mode: 'copy', overwrite: true

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

	input:
	file out_bam from out_bam
	val prefix from prefix2

	output:
	file "${prefix}.bam.bai" into out_bai

	"""
	${params.samtools_folder}/samtools -c index ${out_bam}
	"""
}