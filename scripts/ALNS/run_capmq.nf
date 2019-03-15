#!/usr/bin/env nextflow

/* 
 * Script to run capmq on a CRAM file
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 * and the piece of C code used for doing the CAPMQ cloned from:
 * https://github.com/mcshane/capmq
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.queue = 'standard'
params.threads = 1

//print usage
if (params.help) {
    log.info ''
    log.info 'Script to run capmq on a CRAM file in order to downgrade the MQ'
    log.info '---------------------------------------------------------------'
    log.info ''
    log.info 'This script will downgrade the MQ for each of the CRAM files. The association between'
    log.info 'each of the CRAM files and the maximum MQ value used for the capping is passed in a file'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run_capmq.nf --file input.txt'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--file FILE  File with the urls of the CRAM files to cap and the MQ cap value to use.'
    log.info '    This file needs to have the following format:'
    log.info '    url\tmqcap'
    log.info '    file1.cram,20'
    log.info '    file2.cram,30'
    log.info ''
    exit 1
}

Channel
    .fromPath(params.file)
    .splitCsv(header:true)
    .map{ row-> tuple(row.url, row.mqcap) }
    .set { paths_ch}


process run_capmq {
        /*
        This process run CAPMQ on a certain CRAM file

        Returns
        -------
	Path to the capped CRAM file, it will also index 
	the CRAM file
        */

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	maxForks 100

        input:
        set url, mqcap from paths_ch

        """
	${params.capmq_folder}/capmq -C${mqcap} -O CRAM ${url} ${url}.capped.cram
	${params.samtools_folder}/samtools index ${url}.capped.cram
        """
}

