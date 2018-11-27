Phasing pipeline
================
This workflow is used to generate a phased VCF by using Beagle/Shapeit and emulates the analyses followed during the phase 3 of the 1000 genomes project [ `click here <https://www.nature.com/articles/nature15393>`_] 
It can be run after generating and integrated call set by using the https://github.com/igsr/igsr_analysis/blob/master/PyHive/PipeConfig/INTEGRATION/VCFIntegrationGATKUG.pm workflow.
This pipeline is designed to be run for SNPs or INDELs independently or for both variant types together in the same VCF.
# Important: This workflow can only analyze biallelic variants and it will crash if you try to analyze multiallelic sites.

1. Input preparation:

This pipeline will take as input the VCF file containing either SNPs, INDELs or both type of variants together. 
 In order to generate a combined VCF containing both SNPs+INDELs from 1 SNP VCF+ 1 INDEL VCF you can do the following:
bcftools concat /nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/PRODUCTION/SEQUENCING_GENOTYPES/WHOLE_GENOME/SNPS/working_dir/combined.all.vcf.gz.merged.vcf.gz.recalibrated_snps_raw_indels.vcf.gz /nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/PRODUCTION/SEQUENCING_GENOTYPES/WHOLE_GENOME/INDELs/final_dir/all.merged.20181122.indels.filt.vcf.gz -o combined.snps_indels.vcf.gz -Oz

bcftools sort -T /nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/DEVEL/NEWPHASING/PHASING/tmp_sort combined.snps_indels.vcf.gz -o combined.snps_indels.sorted.vcf.gz -Oz
