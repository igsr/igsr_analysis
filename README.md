# igsr_analysis
This repo contains code that is relevant for the analysis (Mapping, BAM qc, Variant Calling, Filtering etc...) of IGSR data

### Preparing environment
Modify your $PYTHONPATH to include the required libraries:<br>
```export PYTHONPATH=${ehive_dir}/wrappers/python3/:$PYTHONPATH```

Modify your $PERL5LIB to include the required libraries:<br>
```export PERL5LIB=${ehive_dir}/modules/:${igsr_analysis_dir}/:${PERL5LIB}```

Modify your $PATH to include the location of the eHive scripts:
```export PATH=${ehive_dir}/scripts/:${PATH}```

### Install dependency

1) Clone repo by doing ```git clone https://github.com/igsr/igsr_analysis.git``` in the desired folder
2) ```pip install ${igsr_analysis_dir}/dist/igsr_analysis-0.90.tar.gz```

And you are ready to go! 

### Conventions used in this README file:

```${igsr_analysis_dir}``` is the folder where you have cloned https://github.com/igsr/igsr_analysis.git<br>
```${ehive_dir}``` is the folder where you have cloned https://github.com/Ensembl/ensembl-hive.git<br>