# igsr_analysis
This repo contains code that is relevant for the analysis (Mapping, BAM qc, Variant Calling, Filtering etc...) of IGSR data 

The International Genome Sample Resource [(IGSR)](http://www.internationalgenome.org/) is a project funded by the [Wellcome Trust](https://wellcome.ac.uk/) created after the finalization of the 1000 Genomes Project in order to maintain and expand the resource. It has the following aims:

* Ensure the future access to and usability of the 1000 Genomes reference data
* Incorporate additional published genomic data on the 1000 Genomes samples
* Expand the data collection to include new populations not represented in the 1000 Genomes Project 

This repository contains code used in the different analyses pipelines that we use in the project. To use this code please follow the steps described below

The International Genome Sample Resource [(IGSR)](http://www.internationalgenome.org/) is a project funded by the [Wellcome Trust](https://wellcome.ac.uk/) created after the finalization of the 1000 Genomes Project in order to maintain and expand the resource. It has the following aims:

* Ensure the future access to and usability of the 1000 Genomes reference data
* Incorporate additional published genomic data on the 1000 Genomes samples
* Expand the data collection to include new populations not represented in the 1000 Genomes Project 

This repository contains code used in the different analyses pipelines that we use in the project. To use this code please follow the steps described below

### Preparing environment
Modify your $PYTHONPATH to include the required libraries:<br>
```export PYTHONPATH=${ehive_dir}/wrappers/python3/:$PYTHONPATH```

Modify your $PERL5LIB to include the required libraries:<br>
```export PERL5LIB=${ehive_dir}/modules/:${igsr_analysis_dir}/:${PERL5LIB}```

Modify your $PATH to include the location of the eHive scripts:
```export PATH=${ehive_dir}/scripts/:${PATH}```

### Install dependency

1) Clone repo by doing ```git clone https://github.com/igsr/igsr_analysis.git``` in the desired folder
2) ```pip install ${igsr_analysis_dir}/dist/igsr_analysis-0.91.dev0.tar.gz```
3) Modify $PYTHONPATH to add the folder where your pip installs the Python packages

And you are ready to go! 

### Conventions used in this README file:

```${igsr_analysis_dir}``` is the folder where you have cloned https://github.com/igsr/igsr_analysis.git<br>
```${ehive_dir}``` is the folder where you have cloned https://github.com/Ensembl/ensembl-hive.git<br>
