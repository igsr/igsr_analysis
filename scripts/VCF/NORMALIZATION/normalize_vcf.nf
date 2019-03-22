/* 
 * Workflow to normalize a VCF
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

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to normalize a VCF'
    log.info '---------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow normalize_vcf.nf --vcf VCF --threads 5'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be used.'
    log.info '  --ref FASTA  Fasta file with the reference.'
    log.info '  --vt SNP/INDEL  Variant type to analyse.'
    log.info '  --outprefix  STRING Prefix used for final file.'
    log.info '  --threads INT	  Number of threads used by the processes in the pipeline. Default: 1'
    log.info ''
    exit 1
}

process split_multiallelic {
	/*
	This process will split the multiallelic variants by using BCFTools

	Returns
	-------
	Path to splitted VCF
	*/

	memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

	output:
	file "out.splitted.vcf.gz" into out_splitted

	"""
	bcftools norm -m -any ${params.vcf} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
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
	file out_splitted

        output:
        file "out.${params.vt}.vcf.gz" into out_vts

        """
        bcftools view -v ${params.vt} ${out_splitted} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
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
	file out_vts

	output:
	file "out.splitted.decomp.vcf.gz" into out_decomp

	"""
	tabix -f ${out_vts}
	vcfallelicprimitives -k -g ${out_vts} |bgzip -c > out.splitted.decomp.vcf.gz
	"""
}

process run_vt_normalize {
        /*
        Process to run vt normalize

        Returns
        -------
        Path to normalized VCF
        */

        memory '9 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_decomp

        output:
        file "${params.outprefix}.norm.vcf.gz" into out_norm

        """
	vt normalize -r ${params.ref} ${out_decomp} | bgzip -c > ${params.outprefix}.norm.vcf.gz 
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
        file out_norm

        output:
        file "${params.outprefix}.sort.vcf.gz" into out_sort

        """
        bcftools sort ${out_norm} | bgzip -c > ${params.outprefix}.sort.vcf.gz
        """
}

process run_vt_uniq {
        /*
        Process to run vt uniq

        Returns
        -------
        Path to final normalized file
        */

        publishDir 'norm_file', saveAs:{ filename -> "$filename" }

        memory '9 GB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file out_sort

        output:
        file "${params.outprefix}.normalized.vcf.gz"

        """
        vt uniq ${out_sort} | bgzip -c > ${params.outprefix}.normalized.vcf.gz
        """
}

