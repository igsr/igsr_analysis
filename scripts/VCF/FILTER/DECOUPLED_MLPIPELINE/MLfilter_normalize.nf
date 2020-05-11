/* 
 * Workflow to normalize a VCF file
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
    log.info 'Pipeline to normalize a certain file'
    log.info '----------------------------------- '
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow MLfilter_normalize.nf --vcf VCF --prefix chr1_out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be normalized.'
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

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

process replace_header {
        /*
        Process to replace header in the unfiltered VCF
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file newheader

        output:
        file 'unfilt_reheaded.vcf.gz' into unfilt_vcf_reheaded

        """
        bcftools reheader -h ${newheader} -o 'unfilt_reheaded.vcf.gz' ${params.vcf}
        """
}

// normalize vcf

process split_multiallelic {
        /*
        This process will split the multiallelic variants by using BCFTools
        Returns
        -------
        Path to splitted VCF
        */

        memory '70 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

	input:
	file unfilt_vcf_reheaded

        output:
        file "out.splitted.vcf.gz" into out_splitted

        """
        bcftools norm -m -any ${unfilt_vcf_reheaded} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
        """
}

process allelic_primitives {
        /*
        Process to run vcflib vcfallelicprimitives to decompose of MNPs

        Returns
        -------
        Path to decomposed VCF
        */

       	memory { 100.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	errorStrategy 'retry'
        maxRetries 5

        input:
	file out_splitted from out_splitted

        output:
        file "out.splitted.decomp.vcf.gz" into out_decomp

        """
        tabix -f ${out_splitted}
        vcfallelicprimitives -k -g ${out_splitted} |bgzip -c > out.splitted.decomp.vcf.gz
        """
}

process sort_out_decomp {
	/*
	Process to bcftools sort the 'out_decomp'
	*/

	memory '15 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        publishDir "res_norm_${params.prefix}", mode: 'copy', overwrite: true

	input:
	file out_decomp

	output:
	file "out_norm.vcf.gz" into out_norm
	
	"""
	mkdir -p tmpdir
        bcftools sort -T tmpdir/ ${out_decomp} -o out_norm.vcf.gz -Oz
	"""
}
