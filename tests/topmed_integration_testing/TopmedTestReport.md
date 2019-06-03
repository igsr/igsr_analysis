Overview
====================
In production, the pipeline follows the following hierarchy:
`eHive -> bsub -> python eHive Runnable -> Singularity -> topmed component`
eHive catches any errors that occur in bsub, eg memory limits, and logs these.
Any errors thrown by singularity or the topmed component will be caught by the python eHive runnable.
As the capture of errors at the bsub level is a well tested and proven feature of eHive, we will only be concered in this testing protocol with exercising errors raised by topmed comonents. Therefore errors triggered by bsub memory limits will not be necessary.
In these tests a python script, `pywrap.py`, will emulate the functionality of the python eHive runnable. It uses the same command to execute singularity and capture logs.

The tests were excetued at: `/hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing`


# Summary of Results

**Index Bam**

Test                    | Status | Deviations | Mitigation
------------------------|--------|------------|------------
Test 1: Functioning Run | PASS   | NA         | NA
Test 2: Missing files   | PASS   | NA         | NA
Test 3: Truncated files | PASS   | NA         | NA

**STAR**

Test                      | Status | Deviations                  | Mitigation
--------------------------|--------|-----------------------------|-----------------------------------------------------------------------------------------------------
Test 1: Functioning Run   | PASS   | Messages posted to stderr   | Messages are not related to an error. Stderr is only read for troubleshooting. Exit code is still 0.
Test 2: Truncated Run 1   | PASS   | Empty BAM files are created | Exit code is non-zero. Files are only kept on a 0 exit code.
Test 3: Truncated Run 2   | PASS   | Empty BAM files are created | Exit code is non-zero. Files are only kept on a 0 exit code.
Test 4: Missing Reference | PASS   | Empty BAM files are created | Exit code is non-zero. Files are only kept on a 0 exit code.

**RSEM**

Test                      | Status | Deviations                | Mitigation
--------------------------|--------|---------------------------|-----------------------------------------------------------------------------------------------------
Test 1: Functioning Run   | PASS   | Messages posted to stderr | Messages are not related to an error. Stderr is only read for troubleshooting. Exit code is still 0.
Test 2: Truncated Bam     | FAIL   | RSEM runs successfully on a truncated BAM file | TBD
Test 3: Missing Reference | PASS   | NA                   | NA

**Mark Duplicates**

Test                    | Status | Deviations | Mitigation
------------------------|--------|------------|------------
Test 1: Functioning Run | PASS   | Messages posted to stderr | Messages are not related to an error. Stderr is only read for troubleshooting. Exit code is still 0.
Test 2: Truncated Bam  | FAIL | Output file is created | Exit code is non-zero. Files are only kept on a 0 exit code.
Test 3: Missing Bam  |  PASS | NA  | NA

**Rnaseqc Counts**

Test                               | Status | Deviations | Mitigation
-----------------------------------|--------|------------|------------
Test 1: Functioning Run            | PASS   | NA         | NA
Test 2: Missing Bam                | PASS   | NA         | NA
Test 3: Missing Gencode Reference  | PASS   | NA         | NA
Test 4: Missing Reference Sequence | PASS   | NA         | NA
Test 5: Missing Output Dir         | PASS   | NA         | NA



Setting up Resources
====================
Setup example fastq:
```
cd /hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing

seqtk sample -s100 /hps/nobackup/production/reseq-info/galdam/galdam_geuvadis_20190405/FastQs/GEUV:NA20812/sequence_read/ERR188137_1.fastq.gz 10000 | gzip > star/subsampleseq_1.fastq.gz
seqtk sample -s100 /hps/nobackup/production/reseq-info/galdam/galdam_geuvadis_20190405/FastQs/GEUV:NA20812/sequence_read/ERR188137_2.fastq.gz 10000 | gzip > star/subsampleseq_2.fastq.gz

cp star/subsampleseq_1.fastq.gz star/subsampleseq.truncated_1.fastq.gz
truncate -s -900 star/subsampleseq.truncated_1.fastq.gz
cp star/subsampleseq_2.fastq.gz star/subsampleseq.truncated_2.fastq.gz
truncate -s -900 star/subsampleseq.truncated_2.fastq.gz
```
Setup example bamfile:
```
samtools view -s 0.05 -b -h GEUV.NA20812.TSI.star.Aligned.sortedByCoord.out.bam > example_bam_file.bam
cp example_bam_file.bam truncated_bam.bam
truncate -s -900 truncated_bam.bam

cp rsem/Aligned.toTranscriptome.bam rsem/truncated.Aligned.toTranscriptome.bam
truncate -s -900 rsem/truncated.Aligned.toTranscriptome.bam
```


# IndexBam Tests
## Setup Environment
```
cd /hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing
s_img=/nfs/production/reseq-info/work/GTEx-Pipeline/singularity_cache/broadinstitute/gtex_rnaseq:V8
workdir=/hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing/index_bam
```

## IndexBam Test 1: Functioning Run

#### Command:
```
TEST="IndexBam.Test1.functioning_run"
bsub -J $TEST -q production-rh74 python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img samtools index $workdir/example_bam_file.bam"
```

#### Expected Results:
* 1) The file `index_bam/example_bam_file.bam.bai` will be created.
* 2) No messages will be recorded in the error log.
* 3) Exit code will be 0

#### Actual Results:
* 1) PASS: The file `example_bam_file.bam.bai` has been created
* 2) PASS: Error log contains no messages
* 3) PASS: Exit code was 0

**Result:** PASS

## IndexBam Test 2: Missing File
#### Command:
```
TEST="IndexBam.Test2.missing_file"
bsub -J $TEST -q production-rh74 python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img samtools index $workdir/not_a_file.bam"
```
#### Expected Results:
* 1) No `not_a_file.bam.bai` file will be created
* 2) Error log will include reference to `not_a_file.bam` no existing
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: No `not_a_file.bam.bai` file has been created
* 2) PASS: Error log includes the message about the missing file:
```
samtools index: failed to open "[...]/index_bam/not_a_file.bam": No such file or directory
```
* 3) PASS: non-zero exit status 1

**Result:** PASS

## IndexBam Test 3: Truncated File
#### Command:
```
TEST="IndexBam.Test3.truncated_file"
bsub -J $TEST -q production-rh74 python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img samtools index $workdir/truncated_bam.bam"
```
#### Expected Results:
* 1) No `truncated_bam.bam.bai` file will be created
* 2) Error log will include a reference to being unable to read the file.
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: No `truncated_bam.bam.bai` file has been created
* 2) PASS: Error log includes the message about being able to read the file:
```
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
[E::bgzf_read] Read block operation failed with error -1 after 0 of 4 bytes
samtools index: failed to create index for "[...]/index_bam/truncated_bam.bam"
```
* 3) PASS: non-zero exit status 1

**Result:** PASS

# STAR Tests
## Setup Environment
```
cd /hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing
s_img=/nfs/production/reseq-info/work/GTEx-Pipeline/singularity_cache/broadinstitute/gtex_rnaseq:V8
workdir=/hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing/star
refs_dir=/nfs/production/reseq-info/work/GTEx-Pipeline/references/GRCh38
MEM=35000
```

## STAR Test 1: Functioning Run
#### Command:
```
TEST="STAR.Test1.functioning_run"
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_STAR.py $refs_dir/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh74 \
$workdir/subsampleseq_1.fastq.gz $workdir/subsampleseq_2.fastq.gz $TEST --output_dir $workdir/$TEST --threads 8"
```

#### Expected Results:
* 1) No error messages will be recorded in the error log
* 2) Files will be created at `star/STAR.Test1.functioning_run/`:
  * `STAR.Test1.functioning_run.Aligned.sortedByCoord.out.bam`
  * `STAR.Test1.functioning_run.Aligned.toTranscriptome.out.bam`
  * `STAR.Test1.functioning_run.Chimeric.out.sorted.bam`
* 3) Exit code will be 0

#### Actual Results:
* 1) DEVIATION: Error log contains message: `[bam_sort_core] merging from 0 files and 8 in-memory blocks...`
  * Message is not a failure, considering this deviation a pass.
* 2) PASS: The files listed were created.
* 3) PASS: Exit code was 0

**Result:** PASS (with deviations)

## STAR Test 2: Truncated Run 1
#### Command:
```
TEST="STAR.Test2.truncated_run_1"
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_STAR.py $refs_dir/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh74 \
$workdir/subsampleseq.truncated_1.fastq.gz $workdir/subsampleseq_2.fastq.gz $TEST --output_dir $workdir/$TEST --threads 8"
```

#### Expected Results:
* 1) Error log includes reference to incomplete file (fastq 1)
* 2) No bam files created at `star/STAR.Test2.truncated_run_1/`
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes:
```
gzip: [...]/star/subsampleseq.truncated_1.fastq.gz: unexpected end of file
EXITING because of FATAL ERROR in reads input: short read sequence line: 1
```
* 2) DEVIATION: The following bam files were created at `star/STAR.Test2.truncated_run_1/`:
    - `STAR.Test2.truncated_run_1.Aligned.out.bam`
    - `STAR.Test2.truncated_run_1.Aligned.toTranscriptome.out.bam`
  - However, as the files are empty, this deviation can be considered a pass.
* 3) PASS: non-zero exit status 1

**Result:** PASS (with deviations)

## STAR Test 3: Truncated Run 2
#### Command:
```
TEST="STAR.Test3.truncated_run_2"
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_STAR.py $refs_dir/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh74 \
$workdir/subsampleseq_1.fastq.gz $workdir/subsampleseq.truncated_2.fastq.gz $TEST --output_dir $workdir/$TEST --threads 8"
```

#### Expected Results:
* 1) Error log includes reference to incomplete file (fastq 2)
* 2) No bam files created at `star/STAR.Test3.truncated_run_2/`
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes:
```
gzip: [...]/star/subsampleseq.truncated_2.fastq.gz: unexpected end of file
EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length
```
* 2) DEVIATION: The following bam files were created at `star/STAR.Test3.truncated_run_2/`:
    - `STAR.Test3.truncated_run_2.Aligned.out.bam`
    - `STAR.Test3.truncated_run_2.Aligned.toTranscriptome.out.bam`
  - However, as the files are empty, this deviation can be considered a pass.
* 3) PASS: non-zero exit status 1

**Result:** PASS (with deviations)

## STAR Test 4: Missing Reference
#### Command:
```
TEST="STAR.Test4.missing_reference"
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_STAR.py $refs_dir/not_a_STAR_reference $workdir/subsampleseq_1.fastq.gz \
$workdir/subsampleseq.truncated_2.fastq.gz $TEST --output_dir $workdir/$TEST --threads 8"
```

#### Expected Results:
* 1) Error log includes reference to the missing reference file
* 2) No bam files created at `star/STAR.Test4.missing_reference/`
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log contains the following:
```
EXITING because of FATAL ERROR: could not open genome file [..]/not_a_STAR_reference/genomeParameters.txt
SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions
```
* 2) DEVIATION: The following bam files were created at `star/STAR.Test4.missing_reference/`:
    - `STAR.Test4.missing_reference.Aligned.out.bam`
    - `STAR.Test4.missing_reference.Aligned.toTranscriptome.out.bam`
  - However, as the files are empty, this deviation can be considered a pass.
* 3) PASS: non-zero exit status 1

**Result:** PASS (with deviations)

# RSEM Tests
## Setup Environment
```
cd /hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing
s_img=/nfs/production/reseq-info/work/GTEx-Pipeline/singularity_cache/broadinstitute/gtex_rnaseq:V8
workdir=/hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing/rsem
refs_dir=/nfs/production/reseq-info/work/GTEx-Pipeline/references/GRCh38
MEM=700
```

## RSEM Test 1: Functioning Run
#### Command:
```
TEST="RSEM.Test1.functioning_run"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_RSEM.py -o $workdir/$TEST --threads 8 $refs_dir/rsem \
$workdir/Aligned.toTranscriptome.bam $TEST"
```

#### Expected Results:
* 1) No message will be included in error log.
* 2) Files wil be created at `rsem/RSEM.Test1.functioning_run/`:
  * `RSEM.Test1.functioning_run.rsem.genes.results`
  * `RSEM.Test1.functioning_run.rsem.isoforms.results`
* 3) Exit code will be 0

#### Actual Results:
* 1) DEVIATION: Error log contains warning messages from Perl but no messages related to errors.
* 2) PASS: The files have been created.
* 3) PASS: Exit code was 0

**Result:** PASS (with deviations)

## RSEM Test 2: Truncated Bam
#### Command:
```
TEST="RSEM.Test2.truncated_bam"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_RSEM.py -o $workdir/$TEST --threads 8 $refs_dir/rsem \
$workdir/truncated.Aligned.toTranscriptome.bam $TEST"
```

#### Expected Results:
* 1) Error log will include an error message relating to the truncated file
* 2) Results files will not be created at: `rsem/RSEM.Test2.truncated_bam/`
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes the message
```
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated.
```
* 2) FAIL: Files *were* created at `rsem/RSEM.Test2.truncated_bam/`:
  - `RSEM.Test2.truncated_bam.rsem.genes.results`
  - `RSEM.Test2.truncated_bam.rsem.isoforms.results`
* 3) FAIL: Error code *is* 0

**Results:** Fail

## RSEM Test 3: Missing Reference
#### Command:
```
TEST="RSEM.Test3.missing_reference"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -n 8 -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img /src/run_RSEM.py -o $workdir/$TEST --threads 8 $refs_dir/not_an_rsem_reference \
$workdir/truncated.Aligned.toTranscriptome.bam $TEST"
```
#### Expected Results:
* 1) Error log will include an error message relating to the missing file
* 2) Results files will not be created at: `rsem/RSEM.Test3.missing_reference`
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes the message:
`Cannot open [..]/not_an_rsem_reference/rsem_reference.grp! It may not exist.`
* 2) PASS: No results files were created in the output location
* 3) PASS: non-zero exit status 1

**Results:** Pass

# MarkDuplicates Tests
## Setup Environment
```
cd /hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing
s_img=/nfs/production/reseq-info/work/GTEx-Pipeline/singularity_cache/broadinstitute/gtex_rnaseq:V8
workdir=/hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing/mark_duplicates/
MEM=3500
```

## MarkDuplicates Test 1: Functioning Run
#### Command:
```
TEST="MarkDuplicates.Test1.functioning_run"
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img python3 -u /src/run_MarkDuplicates.py $workdir/example_bam_file.bam markdups -o $workdir --memory 5"
```

#### Expected Results:
* 1) Error log contains no error messages
* 2) The file `example_bam_file.md.bam` has been created
* 3) Exit code will be 0

#### Actual Results:
* 1) PASS: The error log includes the output from picard. There are no errors listed. eg:
```
[Fri May 31 15:52:45 UTC 2019] picard.sam.markduplicates.MarkDuplicates done. Elapsed time: 0.84 minutes.
Runtime.totalMemory()=3769159680
```
* 2) PASS: The file `example_bam_file.md.bam` exists.
* 3) PASS: Exit code was 0

**Result:** PASS (with deviations)

## MarkDuplicates Test 2: Truncated Bam
#### Command:
```
TEST="MarkDuplicates.Test2.truncated_bam"
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img python3 -u /src/run_MarkDuplicates.py \
$workdir/truncated_bam.bam markdups -o $workdir --memory 5"
```

#### Expected Results:
* 1) Error log contains reference to the truncated file
* 2) No file created at `truncated_bam.md.bam`
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes the message:
```
Exception in thread "main" htsjdk.samtools.FileTruncatedException: Premature end of file: [..]/mark_duplicates/truncated_bam.bam
```
* 2) FAIL: The file `truncated_bam.md.bam` *does* exist. This is a failing condition.
* 3) PASS: non-zero exit status 1

**Result:** FAIL

## MarkDuplicates Test 3: Missing Bam
#### Command:
```
TEST="MarkDuplicates.Test3.missing_bam"
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --pwd $workdir $s_img python3 -u /src/run_MarkDuplicates.py \
$workdir/not_a_bam_file.bam markdups -o $workdir --memory 5"
```

#### Expected Results:
* 1) Error log contains a reference to the missing file
* 2) The file `not_a_bam_file.md.bam` will not be created.
* 3) Exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log contains the following line:
```
Exception in thread "main" htsjdk.samtools.SAMException:
Cannot read non-existent file: [..]/not_a_bam_file.bam
```
* 2) PASS: The file `not_a_bam_file.md.bam` has not been created.
* 3) PASS: non-zero exit status 1

**Result:** PASS

# RnaseqcCounts Tests
## Setup Environment
```
cd /hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing
s_img=/nfs/production/reseq-info/work/GTEx-Pipeline/singularity_cache/broadinstitute/gtex_rnaseq:V8
workdir=/hps/nobackup/production/reseq-info/galdam/Topmed_Integration_Testing/rnaseq_counts
refs_dir=/nfs/production/reseq-info/work/GTEx-Pipeline/references/GRCh38
MEM=3500
```

## RnaseqcCounts Test 1: Functioning Run
```
TEST="RnaseqcCounts.Test1.functioning_run"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --workdir $workdir $s_img python3 /src/run_rnaseqc.py $workdir/example_bam_file.md.bam \
$refs_dir/gencode/gencode.v29.GRCh38.ERCC.genes.gtf $refs_dir/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta rnaseqccounts \
--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -o $workdir/$TEST --memory 3 --rnaseqc_flags noDoC strictMode"
```

#### Expected Results:
* 1) Error log will contain no messages
* 2) Files will be created at `rnaseq_counts/RnaseqcCounts.Test1.functioning_run`
* 3) Exit status return code will be 0

#### Actual Results:
* 1) PASS: Error log is empty
* 2) PASS: Files were created
* 3) PASS: Exit code was 0

**Results:** PASS

## RnaseqcCounts Test 2: Missing Bam
#### Command:
```
TEST="RnaseqcCounts.Test2.missing_bam"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --workdir $workdir $s_img python3 /src/run_rnaseqc.py $workdir/not_a_bam_file.bam \
$refs_dir/gencode/gencode.v29.GRCh38.ERCC.genes.gtf $refs_dir/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta rnaseqccounts \
--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -o $workdir/$TEST --memory 7 --rnaseqc_flags noDoC strictMode"
```

#### Expected Results:
* 1) Error log will contain message referencing the missing file
* 2) No output files will be created at `rnaseq_counts/RnaseqcCounts.Test2.missing_bam`
* 3) The exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes the message:
```
org.broadinstitute.sting.utils.exceptions.UserException$CouldNotReadInputFile: Couldn't read file [..]/not_a_bam_file.bam because java.io.FileNotFoundException: [..]/not_a_bam_file.bam (No such file or directory)
```
* 2) PASS: Temporary files *do* exist at `rnaseq_counts/RnaseqcCounts.Test2.missing_bam` but none of the output files have been generated
* 3) PASS: Exit code returned: `non-zero exit status 1`

**Results:** PASS (with deviation)

## RnaseqcCounts Test 3: Missing Gencode Reference
#### Command:
```
TEST="RnaseqcCounts.Test3.missing_gencode_reference"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --workdir $workdir $s_img python3 /src/run_rnaseqc.py $workdir/example_bam_file.md.bam \
$refs_dir/gencode/not_a_gencode_refernce.gtf $refs_dir/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta rnaseqccounts \
--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -o $workdir/$TEST --memory 7 --rnaseqc_flags noDoC strictMode"
```

#### Expected Results:
* 1) Error log will reference the missing gencode reference
* 2) No output files generated at `rnaseq_counts/RnaseqcCounts.Test3.missing_gencode_reference/`
* 3) The exit code will be non-zero

#### Actual Results:
* 1) PASS: Error log includes the line: `java.io.FileNotFoundException: /[..]/not_a_gencode_refernce.gtf (No such file or directory)`
* 2) PASS: No output files generated. Only temporary files in the directory.
* 3) PASS: Exit code returned: `non-zero exit status 1`

**Results:** PASS

## RnaseqcCounts Test 4: Missing Reference Sequence
#### Command:
```
TEST="RnaseqcCounts.Test4.missing_reference_seq"
mkdir $workdir/$TEST
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --workdir $workdir $s_img python3 /src/run_rnaseqc.py $workdir/example_bam_file.md.bam \
$refs_dir/gencode/gencode.v29.GRCh38.ERCC.genes.gtf $refs_dir/not_a_genome_reference.fasta rnaseqccounts \
--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -o $workdir/$TEST --memory 7 --rnaseqc_flags noDoC strictMode"
```
#### Expected Results:
* 1) Error log will reference the missing sequence reference file
* 2) No output files generated at `rnaseq_counts/RnaseqcCounts.Test4.missing_reference_seq`
* 3) The exit code will be non-zero

#### Actual Results:
* 1) PASS: The error log includes: `java.io.FileNotFoundException: /[..]/not_a_genome_reference.fasta.fai (No such file or directory)`
* 2) PASS: No output files generated. Only temporary files in the directory.
* 3) PASS: Exit code returned: `non-zero exit status 1`

**Results:** PASS

## RnaseqcCounts Test 5: Missing Output Dir
#### Command:
```
TEST="RnaseqcCounts.Test6.missing_output_dir"
rm -r $workdir/not_an_output_dir
bsub -J $TEST -q production-rh74 -M$MEM -R"select[mem>$MEM] rusage[mem=$MEM]" python3 pywrap.py $TEST \
"singularity exec --workdir $workdir $s_img python3 /src/run_rnaseqc.py $workdir/example_bam_file.md.bam \
$refs_dir/gencode/gencode.v29.GRCh38.ERCC.genes.gtf $refs_dir/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta rnaseqccounts \
--java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -o $workdir/not_an_output_dir --memory 5 --rnaseqc_flags noDoC strictMode"
```

#### Expected Results:
* 1) The error message will reference the missing directory
* 2) The exit code will be non-zero

#### Actual Results:
* 1) PASS: The standard error log includes the message:
```FileNotFoundError: [Errno 2] No such file or directory: '[..]/rnaseq_counts/not_an_output_dir'```
* 2) PASS: Exit code returned: `non-zero exit status 1`

**Results:** PASS
