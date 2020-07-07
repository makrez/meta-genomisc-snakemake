# Analysis Pipeline for Metagenomic assembly

## Pipeline Description

This Pipeline performs metagenomics assembly based on pooled reads from several
samples.

It performs the follwing steps:

1) QC: Software trimmomatic. Adapterclip for Nextera Libraries.

2) Host removal with bowtie2

Here, two different

3) Metaassembly of the genomes with metaSPAdes.


## Files

### Config file

This file is in JSON format.

**Tool specifications**

* {tool}_version: the specific version used for the tool

* {tool}_threads: the number of threads that are used

* {tool}_hours: the maximum number of hours before the job gets killed

* {tool}_mem_mb: the maximum number of RAM (in megabytes) to be used by this tool

* short_sh_commands_threads: Number of threads used for short bash commands. Usually used when one-liners are called.

* short_sh_commands_hours: Number of hours used for short bash commands.

**Data specifications**

* DataFolder: Folder where the data is located

* extension: File extension (only tested with .gz files, however, it could be fastq.gz or anything.gz)

* Lanes: Number of lanes that the samples were sequenced. # //TODO: This may produce errors if e.g. sample_1_L1_[...] and sample_2_L2_[...] are present, but not sample_1_L2_[...] and/or sample_2_L1_[...]

* mates: Paired-end identifier (e.g. R1, R2)

* defaultIlluminaExt: "001" //TODO: If this is not present, it will produce " __ " and throw errors

**Reference Files and other Resources**

* Host reference: The location of the file reference_list.txt (see Reference_files and resources)

**Test parameters**

These parameters are used for all tools for testing purposes with a small sample data set.

* testing_threads: Number of threads used for ALL tools in the testing case

* testing_hours: Number of hours used for ALL tools in the testing case

* testing_mem_mb: Number of memory in MB used for ALL tools in the testing case

### Reference files

The reference files have to be added manually to the pipeline.

The fasta files have to be put into the folder `analysis/02_filter_host/reference` and
a text file has to be generated. The structure of the .txt file looks as follows:

```
ref1.fasta,ref2.fasta # add as many as needed
```

Once the fasta files are in the folder, the easiest way to create `reference_list.txt`
is to run the bash command:

```
ls | grep 'fasta' | paste -sd ',' > reference_list.txt
```

## Creating the test environment

The following commands should produce a test folder. Inside this folder,
the same logic is reproduced but the data folder and the configuration about
threads, memory and running time are changed.

```
chmod + x sync_test.sh # requires rsync
./sync_test.sh
```

This test can be run on the IBU cluster with the command:

```
cd test;
module load Utils/snakemake
/opt/cluster/software/Conda/miniconda/3/bin/snakemake --printshellcmds --drmaa " --partition=phshort --ntasks=1 --mem={resources.mem_mb} --cpus-per-task={threads} --time={resources.hours}:0 --mail-type=END,FAIL " --latency-wait 300 --jobs 1 --jobname <jobname>_{jobid}
```

If this test finishes successfully, the whole script can be run.
