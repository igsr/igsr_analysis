# igsr_analysis
This repo contains code that is relevant for the analysis (Mapping, BAM qc, Variant Calling, Filtering etc...) of IGSR data

### Preparing environment
Modify your $PYTHONPATH to include the required libraries:<br>
```export PYTHONPATH=${igsr_analysis_dir}/igsr_analysis:${ehive_dir}/wrappers/python3/:$PYTHONPATH```

Where ```$igsr_analysis_dir``` and ```$ehive_dir``` are the folders where ```https://github.com/igsr/igsr_analysis.git``` and ```https://github.com/Ensembl/ensembl-hive.git``` have been downloaded
