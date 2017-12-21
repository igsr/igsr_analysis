initial_file="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/FILTERING/05_2017/LC_DIR/BCFTOOLS/WHOLE_GENOME/final_dir/FILTERING_ANALYSIS/GIAB/PRE_FILTER/lc_bams.bcftools.20170319.NA12878.onlyvariants.vcf.gz"

prefix="lc_bams.bcftools.20170319.NA12878"
non_valid_regions="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/exclude_nonvalid.bed"
high_conf_regions="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/chr20DIR/HIGH_CONF_REGIONS/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
giab_snps="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/ANALYSIS_READY/NA12878.giab.SNP.inghighconf.non_valid.vcf.gz"
giab_indels="/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GIAB/ANALYSIS_READY/NA12878.giab.INDEL.inghighconf.non_valid.vcf.gz"
working_dir="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/FILTERING/05_2017/LC_DIR/BCFTOOLS/WHOLE_GENOME/final_dir/FILTERING_ANALYSIS/GIAB/PRE_FILTER"

#select only the desired chros and drop genotypes if present and select PASS if filtered
~/bin/bcftools-1.6/bcftools view -G ${initial_file} -f.,PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX -o ${prefix}.sites.vcf.gz -Oz
tabix ${prefix}.sites.vcf.gz

#exclude non-valid regions
~/bin/bcftools-1.6/bcftools view -T ^${non_valid_regions} ${prefix}.sites.vcf.gz -o ${prefix}.sites.non_valid.vcf.gz -Oz
tabix ${prefix}.sites.non_valid.vcf.gz

#SNPSs
~/bin/bcftools-1.6/bcftools view -v snps ${prefix}.sites.non_valid.vcf.gz -o ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.vcf.gz -O z
tabix SNPS/${prefix}.sites.non_valid.snps.vcf.gz

#select variants in high-conf regions as defined by GIAB
~/bin/bcftools-1.6/bcftools view -R ${high_conf_regions} ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.vcf.gz -o ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.inhighconf.vcf.gz -Oz
tabix ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.inhighconf.vcf.gz

#calculate the intersection
~/bin/bcftools-1.6/bcftools isec -p ${working_dir}/SNPS/dir/ ${working_dir}/SNPS/${prefix}.sites.non_valid.snps.inhighconf.vcf.gz ${giab_snps}

#INDELs
~/bin/bcftools-1.6/bcftools view -v indels ${prefix}.sites.non_valid.vcf.gz -o ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.vcf.gz -O z
tabix ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.vcf.gz

#select variants in high-conf regions as defined by GIAB
~/bin/bcftools-1.6/bcftools view -R ${high_conf_regions} ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.vcf.gz -o ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.inhighconf.vcf.gz -Oz
~/bin/bcftools-1.6/bcftools sort -o ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.inhighconf.sorted.vcf.gz -Oz ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.inhighconf.vcf.gz
tabix ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.inhighconf.sorted.vcf.gz

#calculate the intersection
~/bin/bcftools-1.6/bcftools isec -p ${working_dir}/INDELS/dir/ ${working_dir}/INDELS/${prefix}.sites.non_valid.indels.inhighconf.sorted.vcf.gz ${giab_indels}

