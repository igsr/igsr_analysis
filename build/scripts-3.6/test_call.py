import subprocess

#cmd="java -jar /homes/ernesto/bin/GATK/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R /nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/FILTERING/DEVEL_PIPELINE/SMALL_INPUT/LC_EX_FREEBAYES/working_dir/chr20/lc_ex_bams.freebayes.20170121_NA12878_chr20.vcf.gz.norm.30000000_40000000.vcf.gz -o lc_ex_bams.freebayes.20170121_NA12878_chr20.vcf.gz.norm.vcf.aprimitives.vcf.gz"

cmd="java -jar /homes/ernesto/bin/GATK//GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R /nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/FILTERING/DEVEL_PIPELINE/SMALL_INPUT/LC_EX_FREEBAYES/working_dir/chr20/lc_ex_bams.freebayes.20170121_NA12878_chr20.vcf.gz.norm.30000000_40000000.vcf.gz | /nfs/software/ensembl/RHEL7/linuxbrew/bin//bgzip -c > lc_ex_bams.freebayes.20170121_NA12878_chr20.vcf.gz.norm.30000000_40000000.vcf.gz.aprimitives.vcf.gz"


'''
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)

if p.returncode!=0:
    raise Exception(stderr)
#    raise Exception("Error running cmd: {0}".format(command))
   '''
try:
    out = subprocess.check_output(cmd, shell=True)
    t = 0, out
except subprocess.CalledProcessError as e:
    print(e.returncode)
    raise

print("h")
