#!/usr/bin/env nextflow
/*
* Nextflow script for running TABIX the VCFs that are in a directory 
*
* @author
* Ernesto Lowy <ernesto.lowy@gmail.com>
*
*/

// Define defaults
def defaults = [
    queue: 'production-rh74',
    errorStrategy: 'ignore'
]

// params defaults
params.help = false
params.queue = defaults.queue // lsf queue name
params.errorStrategy = defaults.errorStrategy

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to run TABIX on the VCFs in a directory'
    log.info '------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run_tabix_on_dir.nf --dir \'dir/*.vcf.gz\' '
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --dir DIR  VCF files in a folder.'
    log.info '  --outdir DIR   o put *.tbi files.'
    log.info '  -with-singularity FILE Name of the Singularity image (i.e. variant_filtering.simg).'
    log.info ''
    exit 1
}


Channel.fromPath(params.dir).set { files_ch }
 
process run_tabix {
	/*
	Process to run TABIX on a dir
	*/
	memory '512 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
        errorStrategy "${params.errorStrategy}"

	publishDir "${params.outdir}", mode: 'copy', overwrite: true

	input:
  	file vcf from files_ch
  	
	output:
	file "${vcf}.tbi" into out_ch_tbi
  
	script:
  	"""
	tabix $vcf
  	"""
}