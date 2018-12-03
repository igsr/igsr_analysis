Phasing pipeline
================
This workflow is used to generate a phased VCF by using Beagle/Shapeit and emulates the analyses from the [ `phase 3 <https://www.nature.com/articles/nature15393>`_] of the 1000 genomes project.  
It can be run after generating and integrated call set by using the PyHive::PipeConfig::PHASING.pm <https://github.com/igsr/igsr_analysis/blob/master/PyHive/PipeConfig/INTEGRATION/PHASING.pm>

This pipeline is designed to be run for SNPs or INDELs independently or for both variant types together in the same VCF.
*Important:* This workflow can only analyze biallelic variants and it will crash if you try to analyze multiallelic sites.

1. Input preparation:

This pipeline will take as input the VCF file that contains either SNPs, INDELs or both types together. 
In order to generate a combined VCF containing both SNPs+INDELs from 1 SNP VCF + 1 INDEL VCF you can do the following::

   bcftools concat input1.snps.vcf.gz input2.indels.vcf.gz -o combined.snps_indels.vcf.gz -Oz

After concatenating the SNP with the INDEL file you will need to sort the new VCF::

   bcftools sort -T ./tmp_sort combined.snps_indels.vcf.gz -o combined.snps_indels.sorted.vcf.gz -Oz

2. Preparing the environment

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

3. Databases

  The pipeline uses two databases. They may be on different servers or the
  same server.

  3.1 The ReseqTrack database

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

  3.2 The Hive database

    This is database is used by the Hive code to manage the pipeline and job
    submission etc. The pipeline will be created automatically when you run
    the ``init_pipeline.pl`` script.  Write access is needed to this database.

4. Initialise the pipeline

  The pipeline is initialised with the hive script ``init_pipeline.pl``. Here is
  an example of how to initialise a pipeline::

     init_pipeline.pl PyHive::PipeConfig::INTEGRATION::PHASING \
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

  The rest of the options are documented in the `PyHive::PipeConfig::INTEGRATION::PHASING.pm <https://github.com/igsr/igsr_analysis/blob/master/PyHive/PipeConfig/INTEGRATION/PHASING.pm>`_
  module file. You will probably want to override the defaults for many of
  these options so take a look.

5. Seeding the pipeline

  In order to seed the pipeline with the VCF file that will be analyzed use the hive script
  ``seed_pipeline.pl``::

         seed_pipeline.pl \
                          -url mysql://g1krw:$DB_PASS@mysql-rs-1kg-prod:4175/hive_dbname \
                          -logic_name find_files \
                          -input_id "{ 'file' => '/path/to/file/input_file.txt' }"

  Where ``-url`` controls the Hive database connection details and ``/path/to/file/input_file.txt``
  contains the filename of the VCF to be analyzed.

6. Sync the hive database

  This should always be done before [re]starting a pipeline:

  Run e.g.::

         beekeeper.pl -url mysql://g1krw:{password}@mysql-g1k:4175/my_hive_db_name -sync

  where ``-url`` are the details of your hive database.  Look at the output from
  ``init_pipeline.pl`` to see what your url is.

7. Run the pipeline

  Run e.g.::

    beekeeper.pl -url mysql://g1krw:{password}@mysql-g1k:4175/my_hive_db_name -loop &

  Note the '&' makes it run in the background.

  Look at the pod for ``beekeeper.pl`` to see the various options.  E.g. you might
  want to use the ``-hive_log_dir`` flag so that all ``output/error`` gets recorded in
  files.

  While the pipeline is running, you can check the 'progress' view of the hive
  database to see the current status.  If a job has failed, check the msg
  view.