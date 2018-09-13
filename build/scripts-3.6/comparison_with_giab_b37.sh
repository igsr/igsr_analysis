initial_file="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/P3DIR/VARIANT_CALLING/SAMTOOLS/ALL.wgs.samtools_pass_filter.20130502.snps_indels.low_coverage.sites.vcf.gz"

prefix="ALL.wgs.samtools_pass_filter.20130502.snps_indels.low_coverage"
non_valid_regions="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/B37DIR/exclude_nonvalid.bed"
giab_snps="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/B37DIR/ANALYSIS_READY/NA12878.giab.SNP.non_valid.vcf.gz"
giab_indels="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/B37DIR/ANALYSIS_READY/NA12878.giab.INDEL.non_valid.vcf.gz"
working_dir="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/P3DIR/VARIANT_CALLING/SAMTOOLS/GIAB/"

#select only the desired chros and drop genotypes if present and select PASS if filtered
~/bin/bcftools-1.6/bcftools view -c1 -G ${initial_file} -f.,PASS -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X -o ${prefix}.sites.vcf.gz -Oz
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

