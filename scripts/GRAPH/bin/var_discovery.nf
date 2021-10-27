nextflow.enable.dsl=2
include { GIRAFFE; CHUNK; AUGMENT; SNARLS; PACK; CALL } from "../nf_modules/vg_toolkit.nf"
include { SAVE_FILE } from "../nf_modules/utils.nf"

// params defaults
params.help = false
params.cpus = 1
params.prefix = 'output'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to identify variants using .gam file/s'
    log.info '-----------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow var_discovery.nf --gam input.gam --xgname input.xg --prefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '  --gam   .gam file with alignments.'
    log.info '  --xgname  xg index file.'
    log.info '  --prefix  Output prefix.'
    log.info '  --outdir  Output directory'
    log.info '  --cpus  Number of cpus to use. Default=1.'
    exit 1
}

log.info 'Starting the process.....'

// Check input path parameters to see if they exist
checkPathParamList = [
    params.gam, params.xgname
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// check mandatory parameters
if (!params.gam) exit 1, 'Please specify an input gam file with --gam'
if (!params.xgname) exit 1, 'Please specify a xg index file with --xgname'

workflow  {
    main:
        CHUNK(params.xgname, params.prefix)
        AUGMENT(CHUNK.out.chunkFile, params.gam, CHUNK.out.chunkFile.baseName)
        SNARLS(AUGMENT.out.aug_pgFile, AUGMENT.out.aug_pgFile.baseName)
        PACK(AUGMENT.out.aug_pgFile, AUGMENT.out.aug_gamFile, AUGMENT.out.aug_pgFile.baseName)
        CALL(AUGMENT.out.aug_pgFile, SNARLS.out.snarlsFile, PACK.out.packFile, CHUNK.out.chunkFile.baseName)
        SAVE_FILE(CALL.out.vcfFile, params.outdir, CALL.out.vcfFile.baseName+"_final.vcf", mode='move')
}
