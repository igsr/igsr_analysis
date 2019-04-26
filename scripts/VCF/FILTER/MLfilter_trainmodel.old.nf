/* 
 * Workflow to Train a Logistic Regression ML model in order to be applied in a filtering
 * pipline
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.threads = 1
params.queue = 'production-rh7'
params.rfe = false
params.region = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to train a Logistic Regression binary classifier'
    log.info '---------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow MLfilter_trainmodel.nf --vcf VCF --true VCF --vt snps --annotations ANNOTATION_STRING --threads 5'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be used for training.'
    log.info '  --true VCF  Path to the VCF file containing the gold-standard sites.'
    log.info '  --vt  VARIANT_TYPE   Type of variant to use for training the model. Poss1ible values are 'snps'/'indels'.'
    log.info '  --region chr1,chr2  Comma separated list of chromosomes to analyse. If omitted then all chros will be analysed.'
    log.info '  --annotations ANNOTATION_STRING	String containing the annotations to filter, for example:'
    log.info '	%CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.' 
    log.info '  --rfe BOOL If true, then do Recursive Feature Elimination in order to select the annotations that are more relevant for prediction.'
    log.info '             If false, then train the ML model and skip RFE.'
    log.info '  --no_features INT Number of features that will be selected if params.rfe is true.'
    log.info '  --tmpdir FOLDER What folder to use as tmpdir for bcftools sort.'
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

chrs_splitmultiallelic_withchr=Channel.empty()
chrs_intersecionCallSets=Channel.empty()
chrs_trainModel=Channel.empty()
chrs_rfe=Channel.empty()

if (params.region) {
    log.info '\t--region provided'

    chrList = Channel.from( params.region.split(',') )

    chrList.into { chrs_splitmultiallelic_withchr ; chrs_intersecionCallSets; chrs_trainModel; chrs_rfe}
}

// Normalization
process split_multiallelic {
        /*
        This process will split the multiallelic variants by using BCFTools

        Returns
        -------
        Path to splitted VCF
        */

        memory '2 GB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

        output:
        file "out.splitted.vcf.gz" into out_splitted_nochr

        when:
        !params.region

        """
        bcftools norm -m -any ${params.vcf} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
        """
}


process split_multiallelic_withchr {
        /*
        This process will split the multiallelic variants for a particular region by using BCFTools

        Returns
        -------
        Path to splitted VCF
        */

	tag {"Processing: "+chr}

        memory '2 GB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        val chr from chrs_splitmultiallelic_withchr
	
	when:
	    params.region


        output:
        file "out.splitted.${chr}.vcf.gz" into out_splitted_chr

        """
        bcftools norm -r ${chr} -m -any ${params.vcf} -o out.splitted.${chr}.vcf.gz -Oz --threads ${params.threads}
        """
}

process allelic_primitives {
        /*
        Process to run vcflib vcfallelicprimitives to decompose of MNPs

        Returns
        -------
        Path to decomposed VCF
        */

        memory '9 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
	file out_splitted from out_splitted_chr.mix(out_splitted_nochr)

        output:
        file "out.splitted.decomp.vcf.gz" into out_decomp

        """
        tabix -f ${out_splitted}
        vcfallelicprimitives -k -g ${out_splitted} |bgzip -c > out.splitted.decomp.vcf.gz
        """
}

process select_variants {
        /*
        Process to select the desired variants type (snps/indels)
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_decomp

        output:
        file "out.${params.vt}.vcf.gz" into out_vts

        """
        bcftools view -v ${params.vt} ${out_decomp} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
        """
}


process run_bcftools_sort {
        /*
        Process to run bcftools sort

        Returns
        -------
        Path to sorted VCF
        */

        memory '9 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_vts

        output:
        file "${params.outprefix}.sort.vcf.gz" into out_sort

        """
//        bcftools sort ${out_vts} -o ${params.outprefix}.sort.vcf.gz -Oz -T ${params.tmpdir}
	bcftools sort ${out_vts} -o ${params.outprefix}.sort.vcf.gz -Oz
        """
}

process run_vt_uniq {
        /*
        Process to run vt uniq

        Returns
        -------
        Path to final normalized file
        */

        memory '9 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_sort

        output:
        file "${params.outprefix}.normalized.vcf.gz" into out_uniq

        """
        vt uniq ${out_sort} | bgzip -c > ${params.outprefix}.normalized.vcf.gz
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

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_uniq

        output:
        file "out.onlyvariants.vcf.gz" into out_onlyvariants_chr, out_onlyvariants_nochr

        """
        bcftools view -c1 ${out_uniq} -o out.onlyvariants.vcf.gz --threads ${params.threads} -Oz
        """
}

process intersecionCallSets {
        /*
        Process to find the intersection between out_sites_vts and the Gold standard call set
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_onlyvariants from out_onlyvariants_nochr

        output:
        file 'dir/' into out_intersect_nochr

	when:
        !params.region

        """
        tabix ${out_onlyvariants}
        bcftools isec -c ${params.vt}  -p 'dir/' ${out_onlyvariants} ${params.true}
        """
}

process intersecionCallSets_withchr {
        /*
        Process to find the intersection between out_sites_vts and the Gold standard call set
	for a particular region
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_onlyvariants from out_onlyvariants_chr
        val chr from chrs_intersecionCallSets

        output:
        file 'dir/' into out_intersect_withchr

	when:
        params.region

        """
        tabix ${out_onlyvariants}
        bcftools isec -r ${chr} -c ${params.vt}  -p 'dir/' ${out_onlyvariants} ${params.true}
        """
}

process compressIntersected {
        /*
        Process to compress the files generated by bcftools isec
        */
        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
	file out_intersect from out_intersect_withchr.mix(out_intersect_nochr)

        output:
        file 'FP.vcf.gz' into fp_vcf
        file 'FN.vcf.gz' into fn_vcf
        file 'TP.vcf.gz' into tp_vcf

        """
        bgzip -c ${out_intersect}/0000.vcf > FP.vcf.gz
        bgzip -c ${out_intersect}/0001.vcf > FN.vcf.gz
        bgzip -c ${out_intersect}/0002.vcf > TP.vcf.gz
        """
}

process get_variant_annotations {
        /*
        Process to get the variant annotations for training files
        and for VCF file to annotate (for a single chromosome in this case)
        */

        memory '2 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file tp_vcf
        file fp_vcf

        output:
        file 'TP_annotations.tsv.gz' into tp_annotations_train_nochr, tp_annotations_train_withchr, tp_annotations_rfe_nochr, tp_annotations_rfe_withchr
        file 'FP_annotations.tsv.gz' into fp_annotations_train_nochr, fp_annotations_train_withchr, fp_annotations_rfe_nochr, fp_annotations_rfe_withchr

        """
        bcftools query -H -f '${params.annotations}' ${tp_vcf} | bgzip -c > TP_annotations.tsv.gz
        bcftools query -H -f '${params.annotations}' ${fp_vcf} | bgzip -c > FP_annotations.tsv.gz
        """
}

process train_model {
        /*
        Process that takes TP_annotations.tsv and FP_annotations.tsv created above and will train the Logistic
        Regression binary classifier
        */

        memory '5 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        publishDir "trained_model", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_train_nochr
        file fp_annotations from fp_annotations_train_nochr

        output:
        file 'fitted_logreg_vts.sav' into trained_model_nochr
        file 'fitted_logreg_vts.score' into trained_model_score_nochr

        when:
        !params.rfe && !params.region

        """
        #!/usr/bin/env python

        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.train(outprefix="fitted_logreg_vts",
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}')
        """
}

process train_model_withchr {
        /*
        Process that takes TP_annotations.tsv and FP_annotations.tsv created above and will train the Logistic
        Regression binary classifier
        */

        memory '5 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        publishDir "trained_model_${chr}", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_train_withchr
        file fp_annotations from fp_annotations_train_withchr
        val chr from chrs_trainModel

        output:
        file 'fitted_logreg_vts.sav' into trained_model_withchr
        file 'fitted_logreg_vts.score' into trained_model_score_withchr

        when:
        !params.rfe && params.region

        """
        #!/usr/bin/env python

        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.train(outprefix="fitted_logreg_vts",
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}')
        """
}
 
process rfe {
        /*
        Process to do Recursive Feature Elimination in order to
        select the annotations that are more relevant for prediction'
        */

        memory '5 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        publishDir "selected_feats", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_rfe_nochr
        file fp_annotations from fp_annotations_rfe_nochr

        output:
        file 'selected_feats.txt' into selected_feats_nochr

        when:
        params.rfe && !params.region

        """
        #!/usr/bin/env python
        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.rfe(
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}',
                        n_features='${params.no_features}',
                        outreport="selected_feats.txt")
        """
}

process rfe_with_chr {
        /*
        Process to do Recursive Feature Elimination in order to
        select the annotations that are more relevant for prediction'
        */

        memory '5 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        publishDir "selected_feats_${chr}", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_rfe_withchr
        file fp_annotations from fp_annotations_rfe_withchr
        val chr from chrs_rfe

        output:
        file 'selected_feats.txt' into selected_feats_withchr

        when:
        params.rfe && params.region

        """
        #!/usr/bin/env python
        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.rfe(
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}',
                        n_features='${params.no_features}',
                        outreport="selected_feats.txt")
        """
}
