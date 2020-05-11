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
params.queue = 'production-rh74'
params.executor = 'local'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to apply a certain fitted Logistic Regression model to filter a certain VCF file'
    log.info '-----------------------------------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow MLfilter_applymodel_wonorm.nf --vcf VCF --model MODEL.sav --cutoff 0.95 --threads 5 --vt snps --prefix achr'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be filtered.'
    log.info '  --model FILE Path to serialized ML fitted model.'
    log.info '  --cutoff FLOAT/FLOATs Cutoff value used in the filtering. It also accepts comma-separated list of floats: 0.95,0.96'
    log.info '  --annotations ANNOTATION_STRING String containing the annotations to filter, for example:'
    log.info '    %CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.'
    log.info '  --prefix STR  results${params.prefix} used to put the VCFs in\n.'
    log.info '  --vt  VARIANT_TYPE   Type of variant to filter. Poss1ible values are 'snps'/'indels'.'
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

//Apply a fitted model obtained after running MLfilter_trainmodel.nf


process select_variants {
        /*
        Process to select the desired variants type (snps/indels)
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        output:
	set file("out.snps.gtps.vcf.gz"), file("out.indels.gtps.vcf.gz") into gtps_vts
        file "out.${params.vt}.vcf.gz" into no_gtps_vts

        """
	bcftools view -v snps ${params.vcf} -o out.snps.gtps.vcf.gz --threads ${params.threads} -Oz
        bcftools view -v indels ${params.vcf} -o out.indels.gtps.vcf.gz --threads ${params.threads} -Oz
	bcftools view -G -v ${params.vt} ${params.vcf} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
        """	 
}

process run_vt_uniq {
        /*
        Process to run vt uniq
        Returns
        -------
        Path to final normalized file
        */

        memory '20 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file no_gtps_vts

        output:
        file "out.normalized.vcf.gz" into out_uniq

        """
        vt uniq ${no_gtps_vts} | bgzip -c > out.normalized.vcf.gz
        """
}

process excludeNonVariants {
        /*
        This process will select the variants on the unfiltered vcf (all chros) for
        the particular type defined by 'params.vt'
        Returns
        -------
        Path to a site-VCF file containing just the variants on a particular chromosome
        */

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_uniq

        output:
        file "out.onlyvariants.vcf.gz" into out_onlyvariants

        """
        bcftools view -c1 ${out_uniq} -o out.onlyvariants.vcf.gz --threads ${params.threads} -Oz
        """
}

// start the application of the model

process get_variant_annotations {
	/*
	Process to get the variant annotations for the selected ${params.vt} from the unfiltered VCF file
	*/
	
	memory '8 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

	input:
	file out_onlyvariants

	output:
	file 'unfilt_annotations.vt.tsv.gz' into unfilt_annotations

	"""
	tabix ${out_onlyvariants}
	bcftools view -v ${params.vt} ${out_onlyvariants} -o out.onlyvariants.vt.vcf.gz -Oz --threads ${params.threads}
	tabix out.onlyvariants.vt.vcf.gz
	bcftools query -H -f '${params.annotations}' out.onlyvariants.vt.vcf.gz | bgzip -c > unfilt_annotations.vt.tsv.gz
	"""
}

cutoff_list = Channel.from( params.cutoff.split(',') )

cutoff_list.set { cutoff_values}

process apply_model {
	/*
	Process to read-in the serialized ML model created by running MLfilter_trainmodel.nf
	and to apply this model on the unfiltered VCF
	*/
	tag "Apply model with $cutoff"

	memory { 5.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	errorStrategy 'retry'
        maxRetries 5

	input:
	file unfilt_annotations
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
	tag "Compress predictions with $cutoff"

	memory '1.GB'
        executor 'lsf'
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

reannotate_parameters = gtps_vts.combine(predictions_table)

process reannotate_vcf {
	/*
	Process to reannotate the unfiltered VCF with the information generated after applying the classifier
	*/
	tag "reannotate_vcf with $cutoff"

	memory { 20.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        errorStrategy 'retry'
        maxRetries 5

	publishDir "results_${params.prefix}"

	exec:
	def selected
	def non_selected

	if ("${params.vt}"=="snps") {
           selected="snps"
           non_selected="indels"
        } else if ("${params.vt}"=="indels") {
           selected="indels"
           non_selected="snps"
        }

	input:
        set file("out.snps.gtps.vcf.gz"), file("out.indels.gtps.vcf.gz"), file(predictions_table), file(predictions_table_tabix), val(cutoff) from reannotate_parameters

        output:
        file(output_cutoff)

        script:
        output_cutoff="filt.${cutoff}".replace('.', '_')+".vcf.gz"
        output_cutoff_tabix="filt.${cutoff}".replace('.', '_')+".vcf.gz.tbi"

	"""
        bcftools annotate -a ${predictions_table} out.${selected}.gtps.vcf.gz -c CHROM,POS,FILTER,prob_TP -o reannotated.vcf.gz --threads ${params.threads} -Oz
        bcftools concat reannotated.vcf.gz out.${non_selected}.gtps.vcf.gz -o out.merged.vcf.gz --threads ${params.threads} -Oz
        bcftools sort -T tmpdir/ out.merged.vcf.gz -o ${output_cutoff} -Oz
        tabix ${output_cutoff}
        """
}