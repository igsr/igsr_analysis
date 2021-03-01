/* 
 * Workflow to Calculate the coverage per base on sample/s
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.queue = 'production-rh74'
params.region = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to Calculate coverage per base for sample/s'
    log.info '----------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow cov_per_sample.nf --file aln_paths.txt'
    log.info ''
    log.info 'Options:'
    log.info '    --help	Show this message and exit.'
    log.info '    --file	File containing the paths to the alignment files. The format of the file will be:'
    log.info '              url'
    log.info '				path/to/file.bam'
    log.info '				The header is necessary.'
    log.info '    --genome      File with genome as required by BEDTools.'
    log.info '    --window      Number of bases for each of the BEDTools generated windows.'
    log.info '    --outprefix   Prefix for output file.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

process make_windows {
    /*
    Process to create genomic windows of a certain width (in bases)
    */
 
    memory '500 MB'
    executor 'lsf'
    cpus 1

    output:
    stdout ivals_ch
 
    """
    bedtools makewindows -g ${params.genome} -w ${params.window} |awk -F'\t' '{print \$1":"\$2"-"\$3}' -
    """
}

ivals_ch.splitText().map{it -> it.trim()}.set{monoival_ch}

process get_cov {
	/*
    Process to run SAMTools depth on params.pos_file and get a
	pos file that will be used later
    */
	tag "Depth for ival: $ival"
	
	memory { 20.GB * task.attempt }
	executor 'lsf'
    queue "${params.queue}"
    cpus 1
	maxForks 1000

	errorStrategy 'retry' 
    maxRetries 5

	input:
		val ival from monoival_ch

	output:
		file("*.cov.gz") into cov_chunks

	exec:
	def match = (ival =~ /(chr.*):(\d+)-(\d+)/)
	match.matches()
	def chrom = match.group(1)
	def start = match.group(2)
	def toadd=10-start.length()
    def nstart='0'*toadd+start
	
	script:
	"""
    samtools depth -a -d 0 -r ${ival} -f ${params.file} |bgzip -c > ${chrom}.${nstart}.cov.gz
    """
}

sorted_covchunks = cov_chunks.collect().sort { a,b ->
     return a.baseName <=> b.baseName
}

process merge_chunks {
	/*
	This process will collect and merge each of the genomic
	chunks generated above
	*/

	memory '500 MB'
    executor 'lsf'
    queue "${params.queue}"
    cpus 1

	input:
		file(cov_f) from sorted_covchunks

	output:
		file("merged.cov.gz") into merged_file

	script:
	"""
	zcat $cov_f |gzip -c > merged.cov.gz
	"""
}

process aggregate_depth {
	/*
	This process will aggregate across all sample-level coverages
	to produce an aggregated number
	*/
	publishDir "results", mode: 'copy', overwrite: true

	memory '500 MB'
    executor 'lsf'
    queue "${params.queue}"
    cpus 1

	input:
        file(merged_file) from merged_file

	output:
		file("out.cov.gz") into agg_file

	"""
	sum_covs.py --ifile ${merged_file} --prefix out
	"""
}


