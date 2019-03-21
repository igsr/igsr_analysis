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

process select_variants {
	/*
	Process to select the desired variants type (snps/indels)
	*/

	memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

	output:
	file "out.${params.vt}.vcf.gz" into out_vts
	
	"""
	bcftools view -v ${params.vt} ${params.vcf} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
	"""
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

	input:
	file out_vts

	output:
	file "out.splitted.vcf.gz" into out_splitted

	"""
	bcftools norm -m -any ${out_vts} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
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
        cpus "${params.threads}"

	input:
	file out_splitted

	output:
	file "out.splitted.decomp.vcf.gz" into out_decomp

	"""
	tabix -f ${out_splitted}
	vcfallelicprimitives -k -g ${out_splitted} |bgzip -c > out.splitted.decomp.vcf.gz
	"""
}

process run_vt {
        /*
        Process to run the following mini-piple: vt normalize | sort | uniq

        Returns
        -------
        Path to normalized VCF
        */

	publishDir "norm_result", mode: 'copy', overwrite: true

        memory '2 GB'
        executor 'local'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_decomp

        output:
        file "${params.outprefix}.norm.vcf.gz" into out_norm

        """
	vt normalize -r ${params.ref} ${out_decomp} | bcftools sort - | vt uniq - | bgzip -c > ${params.outprefix}.norm.vcf.gz 
        """
}
