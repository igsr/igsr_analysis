#!/usr/bin/env nextflow
/*
* Nextflow script for selecting a particular variant type ('SNPs'/'INDELs') in the VCFs that are in a directory
* This script relies on BCFTools for this selection. 
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
    log.info 'Pipeline to select SNPs or INDELs from VCFs in a directory'
    log.info '----------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow select_vt.nf --dir \'dir/*.vcf.gz\' --vt snps'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --dir DIR  VCF files in a folder.'
    log.info '  --vt 'snps'/'indels'   Variant type to select.'
    log.info '  --outdir DIR  Name of the output dir, where VCFs with selected'
    log.info '    variants will be saved.'
    log.info '  -with-singularity FILE Name of the Singularity image (i.e. variant_filtering.simg).'
    log.info ''
    exit 1
}


Channel.fromPath(params.dir).set { files_ch }
 
process select_vt {
	/*
	Process to select a particular variant ('snps'/'indels') from VCF
	It will allso generate a .tbi index on the relevant .VCF
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
	file "${vcf.baseName}.${params.vt}.vcf.gz" into out_ch
	file "${vcf.baseName}.${params.vt}.vcf.gz.tbi" into out_ch_tbi
  
	script:
  	"""
	bcftools view -v ${params.vt} $vcf -o ${vcf.baseName}.${params.vt}.vcf.gz -Oz
	tabix ${vcf.baseName}.${params.vt}.vcf.gz
  	"""
}