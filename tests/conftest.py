import pytest

def pytest_addoption(parser):
    parser.addoption('--bam', default='data/exampleBAM.bam' ,action='store_true', help='Path to test BAM file')
    parser.addoption('--faix', default='data/canonical_chros.fa.fai' ,action='store_true', help='Path to faix file')
    parser.addoption('--hive_lib', default='~/lib/ensembl-hive_2.4/' ,action='store_true', help='Path folder containing eHive scripts')
    parser.addoption('--vcf', default='data/test.vcf.gz' ,action='store_true', help='Path to vcf file')
    parser.addoption('--region', default='data/region.bed' ,action='store_true', help='BED file with a small region in chr1')
    parser.addoption('--vcf_chr20', default='data/test_chr20.vcf.gz' ,action='store_true', help='Path to vcf file used for testing Beagle')
    parser.addoption('--vcf_gts', default='data/GLs.HG00136.vcf.gz' ,action='store_true', help='Path to vcf file with GTs')
    parser.addoption('--vcf_gts_ucsc', default='data/GLs.HG00136.ucsc.vcf.gz' ,action='store_true', help='Path to vcf file with GTs with UCSC-style chro names')
    parser.addoption('--vt_folder', default='~/bin/vt/' ,action='store_true', help='Path to folder containing vt binary')
    parser.addoption('--chk_indel_folder', default='~/bin/' ,action='store_true', help='Folder with chk_indel_rg binary')
    parser.addoption('--samtools_folder', default='/homes/ernesto/bin/samtools-1.6/bin/', action='store_true', help='Folder with samtools binary')
    parser.addoption('--java_folder', default='/usr/bin/',action='store_true', help='Folder with java binary')
    parser.addoption('--picard_folder', default='~/bin/',action='store_true', help='Folder with Picard jar file')
    parser.addoption('--bedtools_folder', default='/homes/ernesto/bin/bedtools-2.25.0/bin/', action='store_true', help='Folder with bedtools binary')
    parser.addoption('--makeBGLCHUNKS_folder', default='~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/' ,action='store_true', help='Folder with makeBGLCHUNKS binary')
    parser.addoption('--bcftools_folder', default='~/bin/bcftools-1.6/' ,action='store_true', help='Folder containing bcftools binaries')
    parser.addoption('--beagle_jar', default='beagle.08Jun17.d8b.jar' ,action='store_true', help='Name of Beagle jar file')
    parser.addoption('--beagle_folder', default='~/bin/beagle/' ,action='store_true', help='Folder with Beagle jar file')
    parser.addoption('--prepareGenFromBeagle4_folder', default='/homes/ernesto/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/prepareGenFromBeagle4/bin/', action='store_true', help='Folder with prepareGenFromBeagle4 binary')
    parser.addoption('--ligateHAPLOTYPES_folder', default='/homes/ernesto/bin/shapeit2_v2_12/bin/ligateHAPLOTYPES/bin/', action='store_true', help='Folder with ligateHAPLOTYPES binary')
    parser.addoption('--shapeit_folder', default='~/bin/shapeit2_v2_12/bin/' ,action='store_true', help='Folder with SHAPEIT binary')
    parser.addoption('--gatk_folder', default='~/bin/GATK/', action='store_true', help='Path to folder containing the GATK jar file')
    parser.addoption('--bgzip_folder', default='/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/', action='store_true', help='Path to folder containing the Bgzip binary')
    parser.addoption('--vcflib_folder', default='~/bin/vcflib/bin/', action='store_true', help='Path to folder containing the vcflib binaries')
    parser.addoption('--hostname', default='mysql-g1kdcc-public', action='store_true', help='host name')
    parser.addoption('--username', default='g1krw', action='store_true', help='user name')
    parser.addoption('--port', default=4197, action='store_true', help='port')
    parser.addoption('--pwd', default='test', action='store_true', help='password')
    parser.addoption('--db', default='g1k_archive_staging_track', action='store_true', help='database name')
    parser.addoption('--reference', default='data/exampleFASTA.fasta', action='store_true', help='Path to Fasta file')
    parser.addoption('--vcf_ambiguity', default='data/test.amb.vcf.gz', action='store_true', help='Path to VCF file containing the REF or ALT column with some ambiguity codes')
    parser.addoption('--vcflist', default=['data/test.vcf.gz','data/test1.vcf.gz'], action='store_true', help='List with VCF paths')
    parser.addoption('--snptools_folder', default='~/bin/snptools/', action='store_true', help='Folder with SNPTools binaries')
    parser.addoption('--chr_file', default='data/chr_file.txt', action='store_true', help='File with chros for get_chros function')
    parser.addoption('--glm', default='SNP', action='store_true', help='--glm option for GATK UG') 
    parser.addoption('--output_mode', default='EMIT_ALL_SITES', action='store_true', help='--output_mode option for GATK UG')
