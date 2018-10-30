#!/usr/bin/env nextflow

/* 
 * VCF annotation Nextflow pipeline script
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// Define defaults
def defaults = [
    queue: 'production-rh7'
]

// params defaults
params.help = false
params.queue = defaults.queue // lsf queue name

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to decorate and validate a certain VCF'
    log.info '------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow decorate.nf --phased_vcf test.phased.vcf.gz --ann_vcf ann.vcf.gz --region 20:1-64444167'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --phased_vcf PHASED_VCF	  Phased VCF that will be decorated.'
    log.info '  --ann_vcf ANN_VCF	  VCF with DP annotations and info on the number of samples that will be used to annotate PHASED_VCF .'
    log.info '  --region chr20:1000-2000  Region from PHASED_VCF that will be decorated .'
    log.info '  --sample_panel SAMPLE_PANEL	   File with a table containing the sample ids and the superpopulation they belong to .'
    log.info '  --pops EAS,EUR,AFR,AMR,SAS	   comma separated list of populations or superpopulations that be analysed .'
    log.info '  --exome EXOME .BED		   format file with coordinates of the exomes .'
    log.info '  --tabix TABIX path to tabix binary .'
    log.info '	--vcf_validator VCF_VALIDATOR	   path to vcf_validator binary .'
    log.info '  --igsr_root IGSR_ROOT folder containing igsr codebase .'
    log.info '  --output_dir OUTPUT_DIR	     output folder where the final decorated VCF will be placed .'
    log.info ''
    exit 1
}


process getAlleleFrequencies {
	/*
	Function to calculate AFs on a VCF file

	Returns
	-------
	Path to a tsv file containing the allele frequencies
	*/

	memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	output:
	file 'out_annotate.txt' into out_Annotate
        
	"""
	perl ${params.igsr_root}/scripts/VCF/ANNOTATION/calculate_allele_frq_from_vcf.pl -vcf ${params.phased_vcf} -sample_panel ${params.sample_panel} -out_file out_annotate.txt -region ${params.region} -tabix ${params.tabix} -pop ${params.pops}
	"""
}

process getDepths {
	/*
	Function to get the DP annotation from the VCF file
	*/

	memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	
	output:
	file 'out_depths.txt' into out_Depths

	"""
	bcftools query -f '%CHROM\\t%POS\\t%INFO/DP\\n' -r ${params.region} ${params.ann_vcf} -o out_depths.txt
	"""
}

process processAF {
	/*
	Function to create a table containing the AFs and the depths for each position
	*/

	input:
	file out_Annotate
	file out_Depths

	output:
	file 'out_processAF.txt' into out_processAF
	
	"""
	python ${params.igsr_root}/scripts/VCF/ANNOTATION/process_AFs.py --region ${params.region} --pops ${params.pops} --outfile out_processAF.txt --af_file ${out_Annotate} --depth_f ${out_Depths}
	"""
}

process getOverlappingVariants {
        /*
        Function to create a table containing the AFs and the depths for each posion
        */

 	input:
  	file out_processAF

	output:
	file  'out_getOverlappingVariants.txt' into out_getOverlappingVariants

        """
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/get_overlapping_variants.py --outfile out_getOverlappingVariants.txt --afs_table ${out_processAF} --exome ${params.exome} --label 'EX_TARGET' 
        """
}

process addAnnotation {
	/*
	Function to add the VariantType annotation to the annotation table
	*/

        input:
        file out_getOverlappingVariants

        output:
        file  'out_addAnnotation.txt' into out_addAnnotation

        """
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/add_annotation.py --outfile out_addAnnotation.txt --file1 ${out_getOverlappingVariants} --header 'VT' --label 'SNP'
        """
}

process addNumberSamples {
	  /*
  	  Function to add the NS (number of samples) annotation to the annotation table
  	  */

	  input:
	  file out_addAnnotation

	  output:
	  file  'out_addNumberSamples.txt' into out_addNumberSamples

	  """
	  python ${params.igsr_root}/scripts/VCF/ANNOTATION/add_number_of_samples.py --outfile out_addNumberSamples.txt --file1 ${out_addAnnotation} --ann_vcf ${params.ann_vcf}
	  """
}

process compressAFmatrix {
	/*
	Function to compress output of 'addNumberSamples'
	*/

	input:
	file out_addNumberSamples

	output:
	file  'AFmatrix.txt.gz' into out_AFmatrix_gz
	
	"""
	bgzip -c ${out_addNumberSamples} > AFmatrix.txt.gz
	"""
}

process runAnnotate {
        /*
        Function to run Tabix on a VCF file and decorate VCF and also
	change chromosome names
        */

	memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file out_AFmatrix_gz

	output:
        file  'out_decorate.vcf.gz' into out_decorate

        """
	${params.tabix}/tabix -f -s1 -b2 -e3 ${out_AFmatrix_gz}
	${params.bcftools}/bcftools annotate -r ${params.region} -a ${out_AFmatrix_gz} -h ${params.igsr_root}/SUPPORTING/annots_26062018.txt --rename-chrs ${params.igsr_root}/SUPPORTING/ensembl2ucsc_chrdict.txt -c CHROM,FROM,TO,REF,ALT,DP,AN,AC,AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF,EX_TARGET,VT,NS ${params.phased_vcf} -o out_decorate.vcf.gz -Oz
        """
}

process runReheader {
	/*
	Function to modify the header on the decorated VCF
	*/

	memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	input:
	file out_decorate

	output:
	file 'out_reheaded.vcf.gz' into out_reheaded

	"""
	${params.bcftools}/bcftools reheader -h ${params.igsr_root}/SUPPORTING/header_26062018.txt ${out_decorate} -o out_reheaded.vcf.gz
	"""
}

process runValidator {
	/*
	Function run the VCF validator, it will also store the output
	of the validator in the final folder
	*/
	publishDir 'results', saveAs:{ filename -> "$filename" }

	memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	input:
	file out_reheaded

	output:
	file "${params.chr}.vcf.validation.txt"

	"""
	zcat ${out_reheaded} | ${params.vcf_validator} 2> ${params.chr}.vcf.validation.txt 
	"""
}

process moveFinalFile {
	/*
	Process to move the final output file to the output folder set in params.output_dir
	*/
	publishDir 'results', saveAs:{ filename -> "$filename" }

	input:
	file out_reheaded

	output:
	file "${params.chr}.GRCh38.phased.vcf.gz"
	
	script:
	"""
	mv ${out_reheaded} ${params.chr}.GRCh38.phased.vcf.gz
	"""
}