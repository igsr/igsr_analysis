WES BAM qc pipeline
===================

Workflow used in IGSR to assess the quality of a certain file in the BAM format produced in a Whole Exome Sequencing (WES) experiment.
This workflow consists on running 3 different types of tests:

* Chkindel_rg
  
This test consists on using a simple algorithm to identify runs with unbalanced ratio of short insertion and deletion (greater than 5), which is indicative of low quality data.
The code to run this test can be found at:

https://github.com/lh3/samtools-legacy/blob/master/examples/chk_indel.c 

* VerifyBAMID
  
This test is used to assess sample contamination and sample mix-ups, it uses the VerifyBAMID software.
More info on this useful piece of software can be found at:

https://genome.sph.umich.edu/wiki/VerifyBamID

* Coverage
  
For assessing the coverage we use Picard CollectHsMetrics. 
This software generates a complete report on the depth of coverage for the targeted regions from the WES.
More information on this Picard tool can be found at:

https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics

In order to run this workflow we need to do the following:

1. Preparing the environment

  Modify your ``$PYTHONPATH`` to include the required libraries::

  	 export PYTHONPATH=${ehive_dir}/wrappers/python3/:$PYTHONPATH

  Modify your ``$PERL5LIB`` to include the required libraries::

  	 export PERL5LIB=${ehive_dir}/modules/:${igsr_analysis_dir}/:${PERL5LIB}

  Modify your ``$PATH`` to include the location of the eHive scripts::

  	 export PATH=${ehive_dir}/scripts/:${PATH}

  * Install dependency

    1) Clone repo by doing ``git clone https://github.com/igsr/igsr_analysis.git`` in the desired folder
    2) ``pip install ${igsr_analysis_dir}/dist/igsr_analysis-0.91.dev0.tar.gz``
    3) Modify ``$PYTHONPATH`` to add the folder where your pip installs the Python packages

    And you are ready to go!

  * Conventions used in this section

    * ``${igsr_analysis_dir}`` is the folder where you have cloned ``https://github.com/igsr/igsr_analysis.git``
    * ``${ehive_dir}`` is the folder where you have cloned ``https://github.com/Ensembl/ensembl-hive.git``

2. Databases

  The pipeline uses two databases. They may be on different servers or the
  same server.

  2.1 The ReseqTrack database

    The pipeline queries a ``ReseqTrack`` database to find the VCF that will be
    filtered by the pipeline. It will also add file metadata for the final
    filtered VCF.

    In order to create a ``ReseqTrack`` database use the following 

    commands::

	mysql -h <hostname> -P <portnumber> -u <username> -p???? -e "create database testreseqtrack" # where testreseqtrack 
    	             		      		      	                  		 	     # is the name you want 
												     # to give to the ReseqTrack DB
    	mysql -h <hostname> -P <portnumber> -u <username> -p???? testreseqtrack < $RESEQTRACK/sql/table.sql
    	mysql -h <hostname> -P <portnumber> -u <username> -p???? testreseqtrack < $RESEQTRACK/sql/views.sql

    * Conventions used in this section:
    
     ``$RESEQTRACK`` is the folder where you have cloned ``https://github.com/EMBL-EBI-GCA/reseqtrack.git``

  2.2 The Hive database

    This is database is used by the Hive code to manage the pipeline and job
    submission etc. The pipeline will be created automatically when you run
    the ``init_pipeline.pl`` script.  Write access is needed to this database.

3. Initialise the pipeline
  
  The pipeline is initialised with the hive script ``init_pipeline.pl``. Here is
  an example of how to initialise a pipeline::

     init_pipeline.pl PyHive::PipeConfig::QC::RunBamQCsonWES \
     		      -pipeline_url mysql://g1krw:$DB_PASS@mysql-rs-1kg-prod:4175/hive_dbname \
     		      -db testreseqtrack \
     		      -pwd $DB_PASS \
     		      -hive_force_init 1

  The first argument is the the module that defines this pipeline.  
  Then ``-pipeline_url`` controls the Hive database connection details, in this 
  example::

	 g1krw= username
	 $DB_PASS= password
	 mysql-rs-1kg-prod= hostname
	 4175= Port number
	 hive_dbname= Hive DB name

  Then ``-db`` is the name of the Reseqtrack database name used in the section 2.1
  ``-pwd`` is the ReseqTrack DB password

  The rest of the options are documented in the `PyHive::PipeConfig::QC::RunBamQCsonWES <https://github.com/igsr/igsr_analysis/blob/master/PyHive/PipeConfig/QC/RunBamQCsonWES.pm>`_
  module file. You will probably want to override the defaults for many of
  these options so take a look.

4. Seeding the pipeline

  In order to seed the pipeline with the VCF file that will be analyzed use the hive script 
  ``seed_pipeline.pl``::

	 seed_pipeline.pl \
    	 		  -url mysql://g1krw:$DB_PASS@mysql-rs-1kg-prod:4175/hive_dbname \
    			  -logic_name find_files \
    			  -input_id "{ 'file' => '/path/to/file/input_file.txt' }"

  Where ``-url`` controls the Hive database connection details and ``/path/to/file/input_file.txt`` 
  contains the filename of the VCF to be analyzed. This file must exist in the ReseqTrack database

5. Sync the hive database

  This should always be done before [re]starting a pipeline:

  Run e.g.::

	 beekeeper.pl -url mysql://g1krw:{password}@mysql-g1k:4175/my_hive_db_name -sync

  where ``-url`` are the details of your hive database.  Look at the output from
  ``init_pipeline.pl`` to see what your url is.

6. Run the pipeline

  Run e.g.::

    beekeeper.pl -url mysql://g1krw:{password}@mysql-g1k:4175/my_hive_db_name -loop &

  Note the '&' makes it run in the background.

  Look at the pod for ``beekeeper.pl`` to see the various options.  E.g. you might
  want to use the ``-hive_log_dir`` flag so that all ``output/error`` gets recorded in
  files.

  While the pipeline is running, you can check the 'progress' view of the hive
  database to see the current status.  If a job has failed, check the msg
  view.
