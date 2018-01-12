initial_file="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/FILTERING/05_2017/LC_DIR/GATK/WHOLE_GENOME/final_dir/FILTERING_ANALYSIS/GIAB/POSTFILTER/INCLUDING_OTHER_TRANCHES/lc_bams.gatk.20170720.filt.NA12878.onlyvariants.exp.vcf.gz"

prefix="lc_bams.gatk.20170720.filt.NA12878.exp"
non_valid_regions="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/exclude_nonvalid.bed"
giab_snps="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/ANALYSIS_READY/NA12878.giab.SNP.non_valid.vcf.gz"
giab_indels="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/ANALYSIS_READY/NA12878.giab.INDEL.non_valid.vcf.gz"
working_dir="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/FILTERING/05_2017/LC_DIR/GATK/WHOLE_GENOME/final_dir/FILTERING_ANALYSIS/GIAB/POSTFILTER/INCLUDING_OTHER_TRANCHES/"

#select only the desired chros and drop genotypes if present and select PASS if filtered, also select only variants
~/bin/bcftools-1.6/bcftools view -c1 -G ${initial_file} -f.,PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX -o ${prefix}.sites.vcf.gz -Oz
tabix ${prefix}.sites.vcf.gz

#exclude non-valid regions
~/bin/bcftools-1.6/bcftools view -T ^${non_valid_regions} ${prefix}.sites.vcf.gz -o ${prefix}.sites.non_valid.vcf.gz -Oz
tabix ${prefix}.sites.non_valid.vcf.gz

#SNPSs
~/bin/bcftools-1.6/bcftools view -v snps ${prefix}.sites.non_valid.vcf.gz -o ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.vcf.gz -O z
tabix SNPS/${prefix}.sites.non_valid.snps.vcf.gz

#calculate the intersection
~/bin/bcftools-1.6/bcftools isec -p ${working_dir}/SNPS/dir/ ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.vcf.gz ${giab_snps}

#INDELs
~/bin/bcftools-1.6/bcftools view -v indels ${prefix}.sites.non_valid.vcf.gz -o ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.vcf.gz -O z
tabix ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.vcf.gz

~/bin/bcftools-1.6/bcftools sort -o ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.sorted.vcf.gz -Oz ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.vcf.gz
tabix ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.sorted.vcf.gz

#calculate the intersection
~/bin/bcftools-1.6/bcftools isec -p ${working_dir}/INDELS/dir/ ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.sorted.vcf.gz ${giab_indels}

