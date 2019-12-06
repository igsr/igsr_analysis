#!/usr/bin/env nextflow

/* 
 * Script for selecting a single sample VCF from a multisample VCF
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// Define defaults
def defaults = [
    queue: 'production-rh74'
]

// params defaults
params.help = false
params.queue = defaults.queue // lsf queue name

//print usage
if (params.help) {
    log.info ''
    log.info 'Script to select a sample from a multi-sample VCF'
    log.info '-------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow select_sample.nf --vcf multi.vcf.gz --sample NA12878 --prefix <out.NA12878>'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --vcf VCF	VCF containing the desired sample.'
    log.info '	--sample SAMPLE	  Sample that will be selected.'
    log.info '  --prefix PREFIX	  Prefix for output file.'
    log.info ''
    exit 1
}


process selectSAMPLE {
	/*
	Function to split a multi-chromosome VCF into single chromosome VCF

	Returns
	-------
	Returns 2 files per chromosome:
		1) A VCF format file for selected sample
		2) A tabix index for that VCF
	*/
	publishDir 'final_dir', mode: 'move'

	memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	output:
	file "${params.prefix}.${params.sample}.vcf.gz" into chr_vcf
	file "${params.prefix}.${params.sample}.vcf.gz.tbi" into chr_vcf_tbi

	"""
	bcftools view -s ${params.sample} ${params.vcf} -o ${params.prefix}.${params.sample}.vcf.gz -O z
	tabix ${params.prefix}.${params.sample}.vcf.gz
	"""
}
