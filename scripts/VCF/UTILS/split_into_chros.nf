#!/usr/bin/env nextflow

/* 
 * Script for splitting a VCF containing multiple chromosomes into 1 VCF per chromosome.
 * It will only consider the biallelic sites
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
    log.info 'Pipeline to split a VCF into chromosomes'
    log.info '---------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow split_into_chros.nf --vcf callset.vcf.gz --chros chr1,chr2'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --vcf VCF	VCF that will be splited.'
    log.info '	--prefix PREFIX		  String used for output files.'
    log.info '  --chros LIST_OF_CHROS	  List of chromosome-VCFs to generate.'
    log.info ''
    exit 1
}


chr_str = params.chros.toString()
chr_list = Channel.from( chr_str.split(','))

process splitVCF {
	/*
	Function to split a multi-chromosome VCF into single chromosome VCF

	Returns
	-------
	Returns 2 files per chromosome:
		1) A VCF format file for each splitted chromosome
		2) A tabix index for that VCF
	*/
	publishDir 'final_dir', mode: 'move'

	memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	input:
	val chr from chr_list

	output:
	file "${params.prefix}.${chr}.biallelic.vcf.gz*" into chr_vcf

	"""
	bcftools view -m2 -M2 -r ${chr} ${params.vcf} -o ${params.prefix}.${chr}.biallelic.vcf.gz -O z
	${params.tabix_folder}/tabix ${params.prefix}.${chr}.biallelic.vcf.gz
	"""
}
