nextflow.enable.dsl=2
include { GIRAFFE } from "../nf_modules/vg_toolkit.nf"

// params defaults
params.help = false
params.cpus = 1
params.prefix = 'output'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to align FASTQ file/s to a graph using Giraffe'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow alignment_pipeline.nf --fastq ERR3239334.fastq --outprefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '  --ifile File with comma-separated fields: sample,fastq1,fastq2.'
    log.info '  --gbwtname  GBWT index file.'
    log.info '  --graphname  GBWTGraph file.'
    log.info '  --minimizername  Minimizer index file.'
    log.info '  --distname  Distance index file.'
    log.info '  --prefix  Output prefix.'
    exit 1
}

log.info 'Starting the process.....'

// Check input path parameters to see if they exist
checkPathParamList = [
    params.gbwtname, params.graphname, params.minimizername, params.distname
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// check mandatory parameter
if (!params.ifile) exit 1, 'Please specify an input file with --ifile'

Channel.fromPath(params.ifile)
    .splitCsv()
    .map { row -> tuple(row[0], row[1], row[2])}
    .set { value_list}

workflow  {
    main:
        GIRAFFE( value_list, params.gbwtname, params.graphname, params.minimizername, 
        params.distname, params.cpus)
}
