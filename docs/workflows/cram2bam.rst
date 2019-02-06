CRAM to BAM conversion
======================

This workflow is used in IGSR for:

1) Downloading a CRAM file from the archive (ENA, IGSR FTP, etc...)
2) Convert it to BAM
3) Create an index for converted BAM.

This workflow relies on the ASPERA service for the fast download of
the data from the archives
   
Dependencies
------------

* Nextflow
 This pipeline uses a workflow management system named Nextflow. This software can be downloaded from:

 https://www.nextflow.io/

* SAMTools
 Downloadable from:

 http://www.htslib.org/download/

* Aspera connect software:
 This ascp client can be obtained from:

 http://asperasoft.com/software/transfer-clients/connect-web-browser-plug-in/

How to run the pipeline
-----------------------

* First, you need to create a ``nexflow.config`` file that can be used by Nextflow to set the required variables. Here goes an example of one of these files::

        params.samtools_folder='~/bin/samtools-1.9/' // folder containin the samtools binary
	// params defaults for ascp client
	params.key_file = '/homes/ernesto/.aspera/connect/etc/asperaweb_id_dsa.openssh' // Private-key file name (id_rsa) for authentication
	params.transfer_rate = '900M'
	params.port = 33001 // TCP port used for SSH authentication

* Then, you can start your pipeline by doing::

	nextflow -C nextflow.config run $IGSR_CODEBASE/scripts/FILE/cram2bam.nf --file input.txt

 Where:
  * ``-C`` option allows you to specify the path to the ``nextflow.config`` file
  * ``$IGSR_CODEBASE`` is the folder containing the igsr codebase downloaded from ``https://github.com/igsr/igsr_analysis.git``
  * ``--file`` File with the urls pointing to the CRAM files to be
    converted. This file shold have a content similar to::

	url,dest,prefix
	era-fasp@fasp.sra.ebi.ac.uk:/vol1/ERZ454/ERZ454001/ERR1457180.cram,/path/in/dest/ERR1457180.cram,ERR1457180

    Where ``url`` points to the location of the file to be downloaded, ``dest`` is
    the path in the local machine where it will be downloaded
    and ``prefix`` is used as the string used in the converted BAM file
    and its respective index

Pipeline output
---------------

This worklow will create a folder name ``converted/`` with 2 output files:

* ``prefix.bam`` 

 BAM file resulting after converting the downloaded CRAM file
* ``prefix.bam.bai``

 The index created after running ``samtools index prefix.bam``
