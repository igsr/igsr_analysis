/* 
 * Workflow to filter a VCF
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.split_multiallelics = false
params.threads = 1
params.queue = 'production-rh7'
params.executor = 'local'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to apply a certain fitted Logistic Regression model to filter a certain VCF file'
    log.info '-----------------------------------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow MLfilter_applymodel.nf --vcf VCF --model --cutoff 0.95 --threads 5 --vt snps --split_multiallelics true'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be filtered.'
    log.info '  --model FILE Path to serialized ML fitted model.'
    log.info '  --cutoff FLOAT/FLOATs Cutoff value used in the filtering. It also accepts comma-separated list of floats: 0.95,0.96'
    log.info '  --annotations ANNOTATION_STRING String containing the annotations to filter, for example:'
    log.info '    %CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.'
    log.info '  --vt  VARIANT_TYPE   Type of variant to filter. Poss1ible values are 'snps'/'indels'.'
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info ''
    exit 1
}


chrList=['chr20']
chrChannel=Channel.from( chrList )

chrChannel.into { chrs_splitmultiallelic; chrs_get_variant_annotations; 
		chrs_apply_model; chrs_splitVCF; chrs_reannotatevcf }

//Apply a fitted model obtained after running MLfilter_trainmodel.nf

process get_header {
        /*
        Process to get the header of the unfiltered VCF
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        output:
        file 'header.txt' into header

        """
        bcftools view -h ${params.vcf} > header.txt
        """
}

process modify_header {
        /*
        Process to modify the header of the unfiltered VCF
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file header

        output:
        file 'newheader.txt' into newheader

        """
        #!/usr/bin/env python

        from VCF.VcfUtils import VcfUtils

        vcf_object=VcfUtils(vcf='${params.vcf}')

        vcf_object.add_to_header(header_f='${header}', outfilename='newheader1.txt',
                                 line_ann='##FILTER=<ID=MLFILT,Description="Binary classifier filter">')
        vcf_object.add_to_header(header_f='newheader1.txt', outfilename='newheader.txt',
                                 line_ann='##INFO=<ID=prob_TP,Number=1,Type=Float,Description="Probability of being a True positive">')
        """
}

process splitVCF {
        /*
        This process will select a single chromosome from the VCF
        */
        tag {chr}

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        val chr from chrs_splitVCF

        output:
        file "unfilt.${chr}.vcf.gz" into unfilt_vcf_chr

        """
        bcftools view -r ${chr} ${params.vcf} -o unfilt.${chr}.vcf.gz --threads ${params.threads} -Oz
        """
}

process replace_header {
        /*
        Process to replace header in the unfiltered VCF
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file newheader
        file unfilt_vcf_chr

        output:
        file 'unfilt_reheaded.vcf.gz' into unfilt_vcf_chr_reheaded

        """
        bcftools reheader -h ${newheader} -o 'unfilt_reheaded.vcf.gz' ${unfilt_vcf_chr}
        """
}

process split_multiallelic {
   	   /*
   	   This process is used to split the multiallelic sites into different lines per allele. It uses bcftools norm 
   	   for this
   	    */
	    tag {chr}

	    memory '5 GB'
	    executor 'local'
	    queue "${params.queue}"
	    cpus "${params.threads}"

	    input:
	    val chr from chrs_splitmultiallelic

	    output:
	    file 'out.splitted.vcf.gz' into splitted_vcf
	    file 'out.splitted.vcf.gz.tbi' into splitted_vcf_tbi

	    """
	    bcftools norm -r ${chr} -m -${params.vt} ${params.vcf} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
	    tabix out.splitted.vcf.gz
	    """
}

process get_variant_annotations {
	/*
	Process to get the variant annotations for the selected ${params.vt} from the unfiltered VCF file
	*/
	tag {chr}

	memory '2 GB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

	input:
	val chr from chrs_get_variant_annotations
	file splitted_vcf
	file splitted_vcf_tbi

	output:
	file 'unfilt_annotations.vt.tsv.gz' into unfilt_annotations

	"""
	bcftools view -c1 -r ${chr} -v ${params.vt} ${splitted_vcf} -o out.onlyvariants.vt.vcf.gz -Oz --threads ${params.threads}
	tabix out.onlyvariants.vt.vcf.gz
	bcftools query -H -r ${chr} -f '${params.annotations}' out.onlyvariants.vt.vcf.gz | bgzip -c > unfilt_annotations.vt.tsv.gz
	"""
}

cutoff_values=[0.95,0.96]

process apply_model {
	/*
	Process to read-in the serialized ML model created by running MLfilter_trainmodel.nf
	and to apply this model on the unfiltered VCF
	*/
	tag "Apply model for $chr with $cutoff"

	memory '5 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

	input:
	file unfilt_annotations
        val chr from chrs_apply_model
	each cutoff from cutoff_values

	output:
	file 'predictions.tsv' into predictions
	val cutoff into cutoff_ch

	"""
	#!/usr/bin/env python

	from VCF.VCFfilter.MLclassifier import MLclassifier

	ML_obj=MLclassifier(fitted_model = '${params.model}')

	ML_obj.predict(outprefix="predictions", annotation_f='${unfilt_annotations}', cutoff=${cutoff})
	"""
}

process compress_predictions {
	/*
	Process to compress and index the 'predictions.tsv' file generated by process 'apply_model'
	*/
	tag "Compress predictions for $chr with $cutoff"

	memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

	input:
	file predictions
	val cutoff from cutoff_ch

	output:
	set file('predictions.tsv.gz'), file('predictions.tsv.gz.tbi'), val(cutoff) into predictions_table

	"""
	bgzip -c ${predictions} > 'predictions.tsv.gz'
	tabix -f -s1 -b2 -e2 'predictions.tsv.gz'
	"""
}


reannotate_parameters = unfilt_vcf_chr_reheaded.combine(predictions_table).combine(chrs_reannotatevcf)


process reannotate_vcf {
	/*
	Process to reannotate the unfiltered VCF with the information generated after applying the classifier
	*/
	tag "reannotate_vcf for $chr with $cutoff"

	memory '500 MB'
	executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

	publishDir "results_${chr}", mode: 'copy', overwrite: true

	input:
	set file(unfilt_vcf_chr_reheaded), file(predictions_table), file(predictions_table_tabix), val(cutoff), val(chr) from reannotate_parameters

	output:
	file(output_cutoff)
	
	script:
        output_cutoff="filt.${cutoff}".replace('.', '_')+".vcf.gz"

	"""
	bcftools annotate -a ${predictions_table} ${unfilt_vcf_chr_reheaded} -c CHROM,POS,FILTER,prob_TP -o ${output_cutoff} --threads ${params.threads} -Oz
	"""
}
