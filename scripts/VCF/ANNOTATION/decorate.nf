#!/usr/bin/env nextflow

/*
* VCF annotation Nextflow pipeline script
*
1* @author
* Ernesto Lowy <ernesto.lowy@gmail.com>
*
*/

// Define defaults
def defaults = [
    queue: 'production-rh7',
    errorStrategy: 'ignore'
]

// params defaults
params.help = false
params.queue = defaults.queue // lsf queue name
params.errorStrategy = defaults.errorStrategy 

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
    log.info '  --ifile FILE  File containing the input metadata in order to run the pipeline.'
    log.info '  --phased_vcf PHASED_VCF   Phased VCF that will be decorated.'
    log.info '  --ann_vcf ANN_VCF         VCF with DP annotations and info on the number of samples that will be used to annotate PHASED_VCF .'
    log.info '  --region chr20:1000-2000  Region from PHASED_VCF that will be decorated .'
    log.info '  --sample_panel SAMPLE_PANEL        File with a table containing the sample ids and the superpopulation they belong to .'
    log.info '  --pops EAS,EUR,AFR,AMR,SAS         comma separated list of populations or superpopulations that be analysed .'
    log.info '  --exome EXOME .BED                 format file with coordinates of the exomes .'
    log.info '  --tabix TABIX path to tabix binary .'
    log.info '  --chrnames  chrnotation_dict   path to dictionary containing the correspondence between chr names.'
    log.info '  --vcf_validator VCF_VALIDATOR      path to vcf_validator binary .'
    log.info '  --igsr_root IGSR_ROOT folder containing igsr codebase .'
    log.info '  --output_dir OUTPUT_DIR      output folder where the final decorated VCF will be placed .'
    log.info '  --outprefix OUTPREFIX        String used as basename for output files .'
    log.info '  --header HEADER_FILE         File with new header to use in the final VCF .'
    log.info '  --sample_list                Path to the file with samples present in the VCF to be annotated. These sample ids will be used '
    log.info '                               in the new header .'
    log.info ''
    exit 1
}

Channel
    .fromPath(params.ifile)
    .splitCsv(header:true)
    .map{ row-> tuple(row.phased_vcf, row.ann_vcf, row.region, row.header, row.outprefix, row.outdir) }
    .set { metadata_ch}


process getAlleleFrequencies_snps {
        /*
        Process to calculate AFs on a SNP VCF file

        Returns
        -------
        Path to a tsv file containing the allele frequencies
        */

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch

        output:
        file 'out_annotate.snp.txt' into out_AnnotateSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch1

        """
	echo "Processed region: ${region}"
        ${params.bcftools_folder}/bcftools view -v snps -r ${region} -o out.snp.vcf.gz -Oz ${phased_vcf}
        tabix out.snp.vcf.gz
        perl ${params.igsr_root}/scripts/VCF/ANNOTATION/calculate_allele_frq_from_vcf.pl -vcf out.snp.vcf.gz -sample_panel ${params.sample_panel} -out_file out_annotate.snp.txt -region ${region} -tabix ${params.tabix} -pop ${params.pops}
        """
}

process getAlleleFrequencies_indels {
        /*
        Process to calculate AFs on a INDEL VCF file

        Returns
        -------
        Path to a tsv file containing the allele frequencies
        */

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

	input:
        set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch1

        output:
        file 'out_annotate.indel.txt' into out_AnnotateINDEL
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch2

        """
	echo "Processed region: ${region}"
        ${params.bcftools_folder}/bcftools view -v indels -r ${region} -o out.indels.vcf.gz -Oz ${phased_vcf}
        tabix out.indels.vcf.gz
        perl ${params.igsr_root}/scripts/VCF/ANNOTATION/calculate_allele_frq_from_vcf.pl -vcf out.indels.vcf.gz -sample_panel ${params.sample_panel} -out_file out_annotate.indel.txt -region ${region} -tabix ${params.tabix} -pop ${params.pops}
        """
}

process getDepths_snps {
        /*
        Function to get the DP annotation from the VCF file
        */

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

	input:
        set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch2

        output:
        file 'out_depths.snps.txt' into out_DepthsSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch3

        """
	echo "Processed region: ${region}"
        ${params.bcftools_folder}/bcftools view -v snps -r ${region} -o out.ann.snp.vcf.gz -Oz ${ann_vcf}
        tabix out.ann.snp.vcf.gz
        ${params.bcftools_folder}/bcftools query -f '%CHROM\\t%POS\\t%INFO/DP\\n' -r ${region} out.ann.snp.vcf.gz -o out_depths.snps.txt
        """
}

process getDepths_indels {
        /*
        Function to get the DP annotation from the VCF file
        */

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

	input:
        set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch3

        output:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch4
        file 'out_depths.indels.txt' into out_DepthsINDEL

        """
	echo "Processed region: ${region}"
        ${params.bcftools_folder}/bcftools view -v indels -r ${region} -o out.ann.indel.vcf.gz -Oz ${ann_vcf}
        tabix out.ann.indel.vcf.gz
        ${params.bcftools_folder}/bcftools query -f '%CHROM\\t%POS\\t%INFO/DP\\n' -r ${region} out.ann.indel.vcf.gz -o out_depths.indels.txt
        """
}

process processAF_snps {
        /*
        Function to create a table containing the AFs and the depths for each position
        */

        memory '5 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch4
        file out_AnnotateSNP
        file out_DepthsSNP

        output:
        file 'out_processAFSNP.txt' into out_processAFSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch5

        """
	echo "Processed region: ${region}"
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/process_AFs.py --region ${region} --pops ${params.pops} --outfile out_processAFSNP.txt --af_file ${out_AnnotateSNP} --depth_f ${out_DepthsSNP}
        """
}

process processAF_indels {
        /*
        Function to create a table containing the AFs and the depths for each position
        */

        memory '5 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch5
        file out_AnnotateINDEL
        file out_DepthsINDEL

        output:
	echo "Processed region: ${region}"
        file 'out_processAFINDEL.txt' into out_processAFINDEL
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch6

        """
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/process_AFs.py --region ${region} --pops ${params.pops} --outfile out_processAFINDEL.txt --af_file ${out_AnnotateINDEL} --depth_f ${out_DepthsINDEL}
        """
}


process getOverlappingVariants_snps {
        /*
        Function to create a table containing the AFs and the depths for each posion
        */

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        file out_processAFSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch6

        output:
        file  'out_getOverlappingVariantsSNP.txt' into out_getOverlappingVariantsSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch7

        """
	echo "Processed region: ${region}"
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/get_overlapping_variants.py --outfile out_getOverlappingVariantsSNP.txt --afs_table ${out_processAFSNP} --exome ${params.exome} --label 'EX_TARGET' --bedtools_folder ${params.bedtools_folder}
        """
}

process getOverlappingVariants_indels {
        /*
        Function to create a table containing the AFs and the depths for each posion
        */

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        file out_processAFINDEL
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch7

        output:
        file  'out_getOverlappingVariantsINDEL.txt' into out_getOverlappingVariantsINDEL
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch8

        """
	echo "Processed region: ${region}"
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/get_overlapping_variants.py --outfile out_getOverlappingVariantsINDEL.txt --afs_table ${out_processAFINDEL} --exome ${params.exome} --label 'EX_TARGET' --bedtools_folder ${params.bedtools_folder}
        """
}

process addAnnotation_snps {
        /*
        Function to add the VariantType annotation to the annotation table
        */

        memory '1 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        file out_getOverlappingVariantsSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch8

        output:
        file  'out_addAnnotationSNP.txt' into out_addAnnotationSNP
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch9

        """
	echo "Processed region: ${region}"
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/add_annotation.py --outfile out_addAnnotationSNP.txt --file1 ${out_getOverlappingVariantsSNP} --header 'VT' --label 'SNP'
        """
}

process addAnnotation_indels {
        /*
        Function to add the VariantType annotation to the annotation table
        */

        memory '5 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        file out_getOverlappingVariantsINDEL
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch9

        output:
        file  'out_addAnnotationINDEL.txt' into out_addAnnotationINDEL
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch10

        """
	echo "Processed region: ${region}"
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/add_annotation.py --outfile out_addAnnotationINDEL.txt --file1 ${out_getOverlappingVariantsINDEL} --header 'VT' --label 'INDEL'
        """
}

process mergeAnnotations {
        /*
        Function to merge 2 annotation files
        */

        memory '1 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch10
        file out_addAnnotationSNP
        file out_addAnnotationINDEL

        output:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch11
        file 'out_addAnnotation.txt' into out_addAnnotation

        """
	echo "Processed region: ${region}"
        cat ${out_addAnnotationSNP} ${out_addAnnotationINDEL} > out_addAnnotation.txt
        """
}

process addNumberSamples {
        /*
        Function to add the NS (number of samples) annotation to the annotation table
        */

        memory '1 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch11
        file out_addAnnotation

        output:
        file  'out_addNumberSamples.txt' into out_addNumberSamples
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch12

        """
	echo "Processed region: ${region}"
        python ${params.igsr_root}/scripts/VCF/ANNOTATION/add_number_of_samples.py --outfile out_addNumberSamples.txt --file1 ${out_addAnnotation} --phased_vcf ${phased_vcf}
        """
}

process compressAFmatrix {
        /*
        Function to compress output of 'addNumberSamples'
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
        file out_addNumberSamples
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch12

        output:
        file  'AFmatrix.txt.gz' into out_AFmatrix_gz
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch13

        """
	echo "Processed region: ${region}"
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
	errorStrategy "${params.errorStrategy}"

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch13
        file out_AFmatrix_gz

        output:
        file  'out_decorate.vcf.gz' into out_decorate
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch14

        """
	echo "Processed region: ${region}"
        zcat ${out_AFmatrix_gz} |sort -n -k2,3 | bgzip  > ${out_AFmatrix_gz}.sorted.gz
        ${params.tabix} -f -s1 -b2 -e3 ${out_AFmatrix_gz}.sorted.gz
        ${params.bcftools_folder}/bcftools annotate -r ${region} -a ${out_AFmatrix_gz}.sorted.gz -h ${params.igsr_root}/SUPPORTING/annots_26062018.txt --rename-chrs ${params.chrnames} -c CHROM,FROM,TO,REF,ALT,DP,AN,AC,AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF,EX_TARGET,VT,NS ${phased_vcf} -o out_decorate.vcf.gz -Oz
        """
}

process runReheader {
        /*
        Function to modify the header on the decorated VCF
        This function will also change the chromosome notations
        */

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch14
        file out_decorate

        output:
        file 'out_reheaded.vcf.gz' into out_reheaded
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch15

        """
	echo "Processed region: ${region}"
        ${params.bcftools_folder}/bcftools reheader -h ${header} -s ${params.sample_list} ${out_decorate} -o out_reheaded.vcf.gz
        """
}

process runValidator {
        /*
        Function run the VCF validator, it will also store the output
        of the validator in the final folder
        */
        publishDir "${outdir}", mode: 'copy', overwrite: true

        memory '1 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy "${params.errorStrategy}"

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch15
        file out_reheaded

        output:
        file "${outprefix}.vcf.validation.txt"
	file out_reheaded into out_reheaded1
	set phased_vcf,ann_vcf,region,header,outprefix,outdir into metadata_ch16

        """
	echo "Processed region: ${region}"
        zcat ${out_reheaded} | ${params.vcf_validator} 2> ${outprefix}.vcf.validation.txt
        """
}

process moveFinalFile {
        /*
        Process to move the final output file to the output folder set in params.output_dir
        */
        publishDir "${outdir}", mode: 'copy', overwrite: true

        input:
	set phased_vcf,ann_vcf,region,header,outprefix,outdir from metadata_ch16
        file out_reheaded1

        output:
        file "${outprefix}.GRCh38.phased.vcf.gz*"

        script:
        """
	echo "Processed region: ${region}"
        mv ${out_reheaded1} ${outprefix}.GRCh38.phased.vcf.gz
        ${params.tabix} ${outprefix}.GRCh38.phased.vcf.gz
        """
}