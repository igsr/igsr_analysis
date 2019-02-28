#!/usr/bin/env nextflow

/* 
 * Script to convert a .cram file to a bam file
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.queue = 'production-rh74'
params.threads = 1
// params defaults for ascp client
params.key_file = '/homes/ernesto/.aspera/connect/etc/asperaweb_id_dsa.openssh' // Private-key file name (id_rsa) for authentication
params.transfer_rate = '900M'
params.port = 33001 // TCP port used for SSH authentication
params.wget = true
params.ascp = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Script to download .cram file/s from the archive and convert them to .bam'
    log.info '-------------------------------------------------------------------------'
    log.info ''
    log.info 'This script uses the ASPERA file transfer service for the download.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow cram2bam.nf --file input.txt --cmd wget'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--file FILE  File with the urls pointing to the .cram files to be converted.'
    log.info '    This file needs to have the following format:'
    log.info '    url,dest'
    log.info '    era-fasp@fasp.sra.ebi.ac.uk:/vol1/ERZ454/ERZ454001/ERR1457180.cram,/nfs/production/reseq-info/work/ernesto/isgr/SCRATCH/31_01_2019/ERR1457180.cram'
    log.info '  --wget true|false  If true, then the file will be downloaded using wget'
    log.info '  --ascp  true|false  If true, then the file will be downloaded using aspera'
    log.info ''
    exit 1
}

Channel
    .fromPath(params.file)
    .splitCsv(header:true)
    .map{ row-> tuple(row.url, file(row.dest)) }
    .into { paths_ch1; paths_ch2 }


process downloadFile_byWGET {
        /*
        This process uses WGET for downloading the file

        Returns
        -------
	Path to the CRAM file to be downloaded
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	maxForks 25
	errorStrategy 'ignore'

        input:
        set url, file(dest) from paths_ch1

        output:
        file(dest) into wget_ch
        when:
        params.wget

        """
	wget ${url} -O ${dest}
        """
}

process downloadFile_byASCP {
	/*
	This process uses ASPERA for downloading the file

	Returns
	-------
	Path to the file to be downloaded
	*/

	memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	maxForks 25
	errorStrategy 'ignore' 

        input:
        set url, file(dest) from paths_ch2

	output:
	file(dest) into ascp_ch
	when:
	params.ascp

	"""	
	ascp -i ${params.key_file} -Tr -Q -l ${params.transfer_rate} -P${params.port} -L- ${url} ${dest}
        """
}

process md5 {
        /*
        Calculate md5 for downloaded file

        Returns
        -------
        md5sum of downloaded file
        */

        publishDir "converted", mode: 'copy', overwrite: true

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1
	errorStrategy 'ignore'

        input:
        file down_f from ascp_ch.mix(wget_ch)

        output:
        file down_f into down_f
        file "${down_f}.md5" into md5

        """
        md5sum ${down_f} > ${down_f}.md5
        """
}

process convert2bam {
        /*
        Process to run samtools view to convert the .cram file to .bam format
        It will also create an index for the converted BAM file

        Returns
        -------
        Path to a BAM file
        */

        publishDir "converted", mode: 'copy', overwrite: true

        memory '12 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"
	errorStrategy 'ignore'

        input:
        file down_f from down_f

        output:
        file "${down_f}.bam" into out_bam

        """
        samtools view -b -o ${down_f}.bam ${down_f} --threads ${params.threads}
        """
}

process make_index {
        /*
        Process to create a samtools index on the converted BAM file
        */

        publishDir "converted", mode: 'copy', overwrite: true

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"
	errorStrategy 'ignore'

        input:
        file out_bam from out_bam

        output:
        file "${out_bam}.bai" into out_bai

        """
	samtools index ${out_bam} -@ ${params.threads}
        """
}