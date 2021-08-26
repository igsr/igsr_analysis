#!/usr/bin/env nextflow

/* 
 * Workflow to download a file from the European Nucleotide Archive (ENA). It will also calculate the MD5sum values on the
 * downloaded files. This tool uses two different programs for accessing and downloading the files: wget or ascp (the Aspera command line tool).
 * ascp is generally faster but the choice depends on the user.
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Script to download a file/s from the ENA'
    log.info '----------------------------------------'
    log.info ''
    log.info 'This script uses the ASPERA/WGET for the download.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow downloadFromENA.nf --file input.txt --wget true'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--file FILE  File with the urls pointing to the files to be downloaded.'
    log.info '    This file needs to have the following format:'
    log.info '    url,dest'
    log.info '    /vol1/ERZ454/ERZ454001/ERR1457180.cram,/target/dir/ERR1457180.cram'
    log.info '  --wget true|false  If true, then the file will be downloaded using wget'
    log.info '  --ascp  true|false  If true, then the file will be downloaded using aspera'
    log.info ''
    exit 1
}

pathFile = file(params.file)
if( !pathFile.exists() ) {
  exit 1, "The specified file with paths does not exist: ${params.file}"
}

Channel
    .fromPath(params.file)
    .splitCsv(header:true)
    .map{ row-> tuple(row.url, row.dest) }
    .into { paths_ch1; paths_ch2 }

process downloadFile_byWGET {
        /*
        This process uses WGET for downloading the file

        Returns
        -------
	Path to the CRAM file to be downloaded
        */

        input:
        set url, val(dest) from paths_ch1

        output:
        val(dest) into wget_ch
        when:
        params.wget

        """
	wget ${params.wget_prefix}${url} -O ${dest}
        """
}

process downloadFile_byASCP {
	/*
	This process uses ASPERA for downloading the file

	Returns
	-------
	Path to the file to be downloaded
	*/
        input:
        set url, val(dest) from paths_ch2

	output:
	val(dest) into ascp_ch
	when:
	params.ascp

	"""	
	ascp -i ${params.key_file} -Tr -Q -l ${params.transfer_rate} -P${params.port} -L- ${params.ascp_prefix}${url} ${dest}
        """
}

process md5 {
        /*
        Calculate md5 for downloaded file

        Returns
        -------
        md5sum of downloaded file
        */

        publishDir "MD5DIR", mode: 'copy', overwrite: true

        input:
        file down_f from ascp_ch.mix(wget_ch)

        output:
        file "${down_f}.md5" into md5

        """
        md5sum ${down_f} > ${down_f}.md5
        """
}
