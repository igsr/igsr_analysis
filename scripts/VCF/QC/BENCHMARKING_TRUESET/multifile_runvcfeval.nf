#!/usr/bin/env nextflow

/*
 * VCF benchmarking
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
    log.info 'Wrapper of vcf_eval.nf to analyse multiple VCFs'
    log.info '-----------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow multifile_runvcfeval.nf --vt snps --calc_gtps true --singularity variant_filtering.simg'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --vt  VARIANT_TYPE   Type of variant to benchmark. Possible values are 'snps'/'indels'.'
    log.info '  --calc_gtps BOOL  If true, then calculte the genotype concordance between params.vcf'
    log.info '                    and params.true. If true, the compared VCFs should contain genotype'
    log.info '                    information.'
    log.info '  --singularity FILE  Path to singularity image.'
    log.info ''
    exit 1
}

Channel
    .fromPath(params.file)
    .splitCsv(header:true)
    .map{ row-> tuple(row.vcf, row.true_set, row.bed, row.chr) }
    .set { comp_ch }

process run_vcfeval {
    publishDir 'out/'

    input:
    set val(vcf), val(true_set), val(bed), val(chr) from comp_ch

    output:
    file "${chr}_cmd" into out_file

    script:
    """
    echo nextflow -C /homes/ernesto/lib/igsr_analysis_master/igsr_analysis/scripts/VCF/QC/BENCHMARKING_TRUESET/vcf_eval.config run /homes/ernesto/lib/igsr_analysis_master/igsr_analysis/scripts/VCF/QC/BENCHMARKING_TRUESET/vcf_eval.nf --vcf ${vcf} --true ${true_set} --vt ${params.vt} --chros ${chr} --high_conf_regions ${bed} --calc_gtps ${params.calc_gtps} --with-singularity ${params.singularity} > ${chr}_cmd
    """
}

