#!/bin/bash -ue
/homes/ernesto/bin/bcftools-1.6//bcftools view -c1 -G null -f.,PASS -r null -o out.sites.vcf.gz -Oz
/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/tabix/tabix out.sites.vcf.gz
