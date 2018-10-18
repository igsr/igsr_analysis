#!/usr/bin/env nextflow

/*
 * Pipeline to run SEQTK QA
 * 
 *  @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
*/

// Define defaults
def defaults = [
    min_length: 70
]


// params defaults
params.help = false
params.min_length = defaults.min_length

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to run FASTQA on a FASTQ file'
    log.info '---------------------------------'
    log.info 'FASTQA is based on  SEQTK .'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run_seqtk_qa.nf --pwd $DB_PASS --output_dir ./out --dbname elowy_hgdp_03102018 --runs runs.txt'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--pwd PWD	Password for connecting the ReseqTrack DB.'
    log.info '	--dbname DBNAME	 ReseqTrack DB name.'
    log.info '	--host HOSTNAME	 Hostname of the MySQL ReseqTrack DB.'
    log.info '	--user USER		 Username of the MySQL ReseqTrack DB.'
    log.info '	--port PORT		 Port number of the MySQL ReseqTrack DB.'
    log.info '	--min_length MIN_LENGTH     Minimum sequence length.'
    log.info '	--program PROGRAM	      Path to seqtk binary.'
    log.info '	--output_dir OUTPUT_DIR     Directory that will contain the output files.'
    log.info '	--runs RUN_FILE	      File containing runs to be analyzed'
    log.info '	       The format of RUN_FILE is:'
    log.info '	       runID'
    log.info '	       ERR0001'
    log.info '	       ERR0002'
    log.info '	       ERR0003' 
    log.info ''
    exit 1
}

 /*
  * Validate input file
  */
  run_file = file(params.runs)
  if( !run_file.exists() ) exit 1, "Missing file with run ids: ${run_file}"

 /*
  * Create a channel for different run ids
  */
Channel
    .fromPath(run_file)
    .splitCsv(header:true)
    .map{ row-> row.runId }
    .set { runs_ch }

process runFastqSimpleQA {
        /*
	Process to run perl fastq_simple_qa_by_seqtk.pl on the FASTQs in a RESEQTRACK DB
        */

	memory '2 GB'
        executor 'lsf'
        queue 'standard'
        cpus 1

	input:
	val x from runs_ch

	output:
	file out_fastqa
	val x into run_id

	script:
        """
	perl ${params.reseqtrack}/scripts/qc/fastq_simple_qa_by_seqtk.pl -dbhost ${params.host} -dbname ${params.dbname} -dbuser ${params.user} \
	-dbpass ${params.pwd} -dbport ${params.port} -run_id $x -collection_type FASTQ -new_collection_type FQ_OK -min_length ${params.min_length} \
	-output_dir ${params.output_dir} -program ${params.program} -clobber 1> out_fastqa 
        """
}

process moveFinalFile {
	/*
	Process to move the final output file to the output folder set in params.output_dir
	*/
	publishDir "${params.output_dir}", saveAs:{ filename -> "$filename" }

	input:
        val run_id
	file out_fastqa

	output:
	file "${run_id}.seqtk.out"
	
	script:
	"""
	mv ${out_fastqa} ${run_id}.seqtk.out
	"""
}