from VcfQC import VcfQC


if __name__ == '__main__':
    '''
    vcfQC = VcfQC(vcf="/Users/ernesto/projects/IGSR/13_02_17/data/lc_bams.samtools.20170111.NA20818.vcf.gz",vcftools_folder="/usr/local/bin/")
    vcfQC.subset_vcf(bed="/Users/ernesto/projects/IGSR/REFERENCE/centromeres_locations.bed",outprefix="lc_bams.samtools.20170111.NA20818.nocentromeres",outdir="/Users/ernesto/projects/IGSR/13_02_17/data")
    
    vcfQC = VcfQC(vcf="/Users/ernesto/projects/IGSR/13_02_17/data/lc_bams.samtools.20170111.NA20818.vcf.gz",bcftools_folder="/usr/local/bin/")
    vcfQC.stats(outprefix="lc_bams.samtools.20170111.NA20818",outdir="/Users/ernesto/projects/IGSR/13_02_17/data")
    '''
    
    vcfQC = VcfQC(vcf="/Users/ernesto/projects/IGSR/13_02_17/data/lc_bams.samtools.20170111.NA20818.vcf.gz",picard_folder="/usr/local/bin/")
    vcfQC.calc_concordance(truth_vcf="ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.hg38.vcf.ucsc.reheaded.refcorrected.final.bgz",truth_sample="HG02143",call_sample="HG02143",outprefix="HG02143_LCvsOMNI", intervals="test.txt")
    print("h")