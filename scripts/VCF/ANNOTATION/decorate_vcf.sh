#
# Script to decorate a VCF with some annotations. It also changes its header and finally it validates it
# ------
# author: Ernesto Lowy (ernestolowy@gmail.com)
#
# vcf_validator_linux, bcftools need to be in $PATH
# Set the IGSR_ROOT env variable  to the folder where git https://github.com/igsr/igsr_analysis.git has been cloned

usage="USAGE: sh annotate.py <phased.vcf.gz> <ann.vcf.gz> <20:10000000-10050000> <EAS,EUR,AFR,AMR,SAS> <exome> <tabix_bin>"

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ] || [ "$5" == "" ] || [ "$6" == "" ]; then
    echo $usage
    exit 0
fi

region=`echo $3 |sed s/:/./g`

#create folder for annotation table
mkdir ./ann_table
rm -f ./ann_table/*

#create the annotation table
cmd1="python $IGSR_ROOT/scripts/VCF/ANNOTATION/annotate.py --AFcalc $IGSR_ROOT/scripts/VCF/ANNOTATION/ --phased_vcf $1 --sample_panel $IGSR_ROOT/SUPPORTING/integrated_allsamples.20180619.superpopulations.panel --tabix $6 --region $3 --pops $4 --exome $5 --outdir ./ann_table --ann_vcf $2"
echo "[INFO] running Python annotate.py"
$cmd1
echo "[INFO] running Python annotate.py-DONE"

echo "[INFO] compression annotation table"
cmd2="bgzip ann_table/annot_tab2.${region}.txt"
$cmd2
echo "[INFO] compression annotation table-DONE"

echo "[INFO] running tabix on the annotation table"
cmd3="tabix -f -s1 -b2 -e3 ann_table/annot_tab2.${region}.txt.gz"
$cmd3
echo "[INFO] running tabix on the annotation table-DONE"

#decorate the VCF and also change chromosome names
echo "[INFO] running bcftools annotate"
dec_out=$(basename -s .vcf.gz $1).$region".decorated.vcf.gz"
cmd4="bcftools annotate -r $3 -a ann_table/annot_tab2.${region}.txt.gz -h $IGSR_ROOT/SUPPORTING/annots_26062018.txt --rename-chrs $IGSR_ROOT/SUPPORTING/ensembl2ucsc_chrdict.txt -c CHROM,FROM,TO,REF,ALT,DP,AN,AC,AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF,EX_TARGET,VT,NS $1 -o $dec_out -Oz"
$cmd4
echo "[INFO] running bcftools annotate-DONE"

#Now, modify the header on the decorated VCF
echo "[INFO] running bcftools reheader"
reheaded_out=$(basename -s .vcf.gz $dec_out)".reheaded.vcf.gz"
cmd5="bcftools reheader -h $IGSR_ROOT/SUPPORTING/header_26062018.txt $dec_out -o $reheaded_out"
$cmd5
echo "[INFO] running bcftools reheader-DONE"

#validate the VCF
echo "[INFO] running vcf validator"
`zcat $reheaded_out | /nfs/production/reseq-info/work/ernesto/bin/vcf_validator/vcf_validator_linux` 
echo "[INFO] running vcf validator-DONE"
