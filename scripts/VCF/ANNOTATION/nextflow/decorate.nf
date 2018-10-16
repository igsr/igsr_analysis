#!/usr/bin/env nextflow

/* 
 * VCF annotation pipeline script
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */


process getAlleleFrequencies {
	/*
	Function to calculate AFs on a VCF file

	Returns
	-------
	Path to a tsv file containing the allele frequencies
	*/

	output:	
        file 'out_annotate.txt' into out_Annotate

	"""
	perl ${params.scripts_dir}/calculate_allele_frq_from_vcf.pl -vcf ${params.phased_vcf} -sample_panel ${params.sample_panel} -out_file out_annotate.txt -region ${params.region} -tabix ${params.tabix} -pop ${params.pops}
	"""
}

process getDepths {
	/*
	Function to get the DP annotation from the VCF file
	*/
	
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
	python ${params.scripts_dir}/process_AFs.py --region ${params.region} --pops ${params.pops} --outfile out_processAF.txt --af_file ${out_Annotate} --depth_f ${out_Depths}
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
        python ${params.scripts_dir}/get_overlapping_variants.py --outfile out_getOverlappingVariants.txt --afs_table ${out_processAF} --exome ${params.exome} --label 'EX_TARGET' 
        """
}

process addAnnotation {
	/*
	Function to create a table containing the AFs and the depths for each posion
	*/

        input:
        file out_getOverlappingVariants

        output:
        file  'out_addAnnotation.txt' into out_addAnnotation

        """
        python ${params.scripts_dir}/add_annotation.py --outfile out_addAnnotation.txt --file1 ${out_getOverlappingVariants} --header 'VT' --label 'SNP'
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
	  python ${params.scripts_dir}/add_number_of_samples.py --outfile out_addNumberSamples.txt --file1 ${out_addAnnotation} --ann_vcf ${params.ann_vcf}
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

        input:
        file out_AFmatrix_gz

	output:
        file  'out_decorate.vcf.gz' into out_decorate

        """
        tabix -f -s1 -b2 -e3 ${out_AFmatrix_gz}
	bcftools annotate -r ${params.region} -a ${out_AFmatrix_gz} -h ${params.igsr_root}/SUPPORTING/annots_26062018.txt --rename-chrs ${params.igsr_root}/SUPPORTING/ensembl2ucsc_chrdict.txt -c CHROM,FROM,TO,REF,ALT,DP,AN,AC,AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF,EX_TARGET,VT,NS ${params.phased_vcf} -o out_decorate.vcf.gz -Oz
        """
}

process runReheader {
	/*
	Function to modify the header on the decorated VCF
	*/

	input:
	file out_decorate

	output:
	file 'out_reheaded.vcf.gz' into out_reheaded

	"""
	bcftools reheader -h ${params.igsr_root}/SUPPORTING/header_26062018.txt ${out_decorate} -o out_reheaded.vcf.gz
	"""
}

process runValidator {
	/*
	Function run the VCF validator
	*/

	input:
	file out_reheaded

	"""
	zcat ${out_reheaded} | /nfs/production/reseq-info/work/ernesto/bin/vcf_validator/vcf_validator_linux
	"""
}