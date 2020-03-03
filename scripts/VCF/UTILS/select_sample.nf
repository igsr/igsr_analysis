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
params.sample = false
params.samplel = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Script to select a sample from a multi-sample VCF'
    log.info '-------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow select_sample.nf --vcf multi.vcf.gz --sample NA12878 --samplel sample_list.txt'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --vcf VCF	VCF containing the desired sample.'
    log.info '	--sample SAMPLE	  Sample that will be selected.'
    log.info '  --samplel FILE   File with samples that will be analysed.'
    log.info '  --prefix PREFIX	  Prefix for output file.'
    log.info ''
    exit 1
}

process select_sample {
	/*
	Use BCFTools to produce a VCF for sample in params.sample
	*/
	publishDir 'final_dir', mode: 'move'

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	output:
        file "${params.prefix}.${params.sample}.vcf.gz" into chr_vcf
        file "${params.prefix}.${params.sample}.vcf.gz.tbi" into chr_vcf_tbi
	when:
  	params.sample
	
        script:
	"""
	bcftools view -s ${params.sample} ${params.vcf} -o ${params.prefix}.${params.sample}.vcf.gz -O z
        tabix ${params.prefix}.${params.sample}.vcf.gz
	"""
}

if (params.samplel) {
   Channel.fromPath(params.samplel)
        .splitText()
	.map{it -> it.trim()}
        .set { sample_list }  

	process select_sample_list {
		/*
		Use BCFTools to generate a different VCF for each of the samples
		in params.samplel
		*/

		publishDir 'final_dir', mode: 'move'

		memory '1 GB'
		executor 'lsf'
		queue "${params.queue}"
		cpus 1

		input:
		val s from sample_list

		output:
        	file "${params.prefix}.${s}.vcf.gz" into chr_vcf_s
        	file "${params.prefix}.${s}.vcf.gz.tbi" into chr_vcf_s_tbi
		when:
		params.samplel

		script:
		"""
		bcftools view -c1 -s ${s} ${params.vcf} -o ${params.prefix}.${s}.vcf.gz -O z
		tabix ${params.prefix}.${s}.vcf.gz
		"""
	}		
}


