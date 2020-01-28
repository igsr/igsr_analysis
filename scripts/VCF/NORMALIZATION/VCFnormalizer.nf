/* 
 * Workflow to Normalize a certain VCF
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
    log.info 'Pipeline to Normalize a VCF'
    log.info '---------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow VCFnormalizer.nf --vcf VCF --vt snps --threads 5 --outprefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be normalized.'
    log.info '  --vt  VARIANT_TYPE   Type of variant that will be selected and normalized. Poss1ible values are 'snps'/'indels'.'
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info '  --outprefix OUTPREFIX Prefix for output files.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

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
        file "out.splitted.vcf.gz" into out_splitted

        """
        bcftools norm -m -any ${params.vcf} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
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
	file out_splitted from out_splitted

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

        memory '15 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_vts

        output:
        file "${params.outprefix}.sort.vcf.gz" into out_sort

        """
	mkdir tmpdir
	bcftools sort -T tmpdir/ ${out_vts} -o ${params.outprefix}.sort.vcf.gz -Oz
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

	publishDir "normalized_vcf", mode: 'copy', overwrite: true

        input:
        file out_sort

        output:
        file "${params.outprefix}.normalized.vcf.gz" into out_uniq

        """
        vt uniq ${out_sort} | bgzip -c > ${params.outprefix}.normalized.vcf.gz
        """
}
