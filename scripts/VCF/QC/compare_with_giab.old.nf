#!/usr/bin/env nextflow

/* 
 * VCF benchmarking using GIAB (Genome in a Bottle)
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to benchmark a VCF using GIAB'
    log.info '--------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow benchmark_giab.nf --pwd $DB_PASS --output_dir ./out --dbname elowy_hgdp_03102018 --runs runs.txt'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be assessed.'
    log.info '  --chros CHROSTR	  Chromosomes that will be analyzed: chr20 or chr1,chr2 or chr20:10000000-10000100.'
    log.info '  --prefix STR	  String used as prefix for the output VCF file.'
    log.info ''
    exit 1
}

process dropGTPs {
	/*
	This process will drop the genotypes from the initial VCF and will also select the variants and the desired chromosomes.
	Additionally, only the variants with the 'PASS' label in the filter column are considered

	Returns
	-------
	Path to a VCF file containing just the filtered sites
	*/

	output:
	file 'out.sites.vcf.gz' into out_sites_vcf

	"""
	${params.bcftools_folder}/bcftools view -c1 -G ${params.vcf} -f.,PASS -r ${params.chros} -o out.sites.vcf.gz -Oz
	${params.tabix} out.sites.vcf.gz
	"""
}
