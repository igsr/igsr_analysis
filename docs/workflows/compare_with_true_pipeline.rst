Benchmarking with a true set
============================

This workflow is used in our project to benchmark the sites identified in our call sets for a certain sample with a reference call set for which the calls have enough quality to be considered true positives.

Additionally, this workflow will also compare the phased genotypes in our call sets with the ones identified in the reference call set.

Dependencies
------------

* Nextflow

 This pipeline uses a workflow management system named Nextflow. This software can be downloaded from:

 https://www.nextflow.io/

* Bgzip and Tabix

 Bgzip and Tabix are part of the HTSlib project, which can be downloaded from:

 https://github.com/samtools/htslib

* BCFTools

 Downloadable from:

 http://www.htslib.org/download/

* IGSR-analysis code base

 The scripts needed to run this workflow can be downloaded by cloning the IGSR-analysis github repo from:

 https://github.com/igsr/igsr_analysis.git

* vcflib

  This library can be downloaded from:
  
  https://github.com/vcflib/vcflib

Using Docker or Singularity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to skip the installation of these dependencies you can pull the Docker Image from::

	docker pull elowy01/igsr_analysis:latest

or build the Singularity image by using::

	singularity build igsr_analysis.simg docker://elowy01/igsr_analysis
  
How to run the pipeline
-----------------------

* You can start your pipeline by doing::

	nextflow run $IGSR_CODEBASE/scripts/VCF/QC/BENCHMARKING_TRUESET/vcf_eval.nf --vt snps --vcf VCF --true TRUE.vcf.gz --filt_str "PASS,." --chros chr20
	--high_conf_regions highconf.file.bed --calc_gtps true -with-singularity igsr_analysis.simg

 Where:
  * ``$IGSR_CODEBASE`` is the folder containing the igsr codebase downloaded from ``https://github.com/igsr/igsr_analysis.git``
  * ``--vcf`` is the VCF that will be benchmarked with the true set. Notice that you will need to create a tabix index of this file before running this pipeline
  * ``--true`` is the path to the VCF containing the true call set
  * ``--chros`` is the chromosome that will be analyzed. i.e.: ``chr20`` 
  * ``--vt`` is the parameter used to set the type of variants that will be analyzed. i.e. ``'snps'/'indels'``
  * ``--calc_gtps`` if true, then the genotype concordance between ``--true`` and ``--vcf`` will be calculated
  * ``--filt_str`` Filter string used by BCFTools in order to subset a certain subset of variants to be analysed. i.e. ``"PASS,."``
  * ``--high_conf_regions`` BED file used to control the genomic regions for which the benchmarking will be done. This parameter is optional
  * ``--calc_gtps`` If 'true' then the genotype concordance between ``--vcf`` and ``--true`` will be calculated
  * ``-with-singularity`` Parameter used to specify the Singularity image that Nextflow will use to run the workflow. This parameter is optional
  

Pipeline output
---------------

This worklow will create a folder name ``results/`` with the following relevant files:

* ``TP_true.vcf.gz``

Will contain the set of sites that were idendified both in our call set and in the true call set

* ``TP.stats``

Are the stats calculated by running ``bcftools stats TP.vcf.gz``

* ``FP.vcf.gz``

Will contain the set of sites identified in our call set and absent in the true call set

* ``FP.stats``

Are the stats calculated by running ``bcftools stats FP.vcf.gz``

* ``FN.vcf.gz``

Will contain the set of sites that were not idendified in our call set and are present in the true call set

* ``FN.stats``

Are the stats calculated by running ``bcftools stats FN.vcf.gz``

* ``TP_true.highconf.vcf.gz``

Will contain the set of sites from the true call set that were idendified both in our call set and in the true set but restricted to
the regions passed with ``params.high_conf_regions``

* ``TP_target.highconf.vcf.gz``

Will contain the set of sites from the target call set that were idendified both in our call set and in the true call set but restricted to
the regions passed with ``params.high_conf_regions``

* ``TP.highconf.stats``

Are the stats calculated for the sites identified both in the true and target call sets

* ``FP.highconf.vcf.gz``

Will contain the set of sites identified in our call set and absent in the true call set but restricted to
the regions passed with ``params.high_conf_regions``

* ``FP.highconf.stats``

Are the stats calculated by running ``bcftools stats FP.highconf.vcf.gz``

* ``FN.highconf.vcf.gz``

Will contain the set of sites that were not idendified in our call set and are present in the true call set but restricted to
the regions passed with ``params.high_conf_regions``

* ``FN.highconf.stats``

Are the stats calculated by running ``bcftools stats FN.highconf.vcf.gz``

* ``GT_concordance.txt``

This file contains the tables produced after comparing the phased genotypes in our call set with the true call set

