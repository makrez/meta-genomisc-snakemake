#!/bin/bash
module load Utils/snakemake
#/opt/cluster/software/Conda/miniconda/3/bin/snakemake --printshellcmds --verbose -p --drmaa " --partition=pall --ntasks=1 --mem={resources.mem_mb} --cpus-per-task={threads} --time={resources.hours}:0 --mail-type=END,FAIL " --latency-wait 300 --jobs 2 --jobname p532_{jobid}
/opt/cluster/software/Conda/miniconda/3/bin/snakemake --printshellcmds --verbose --drmaa -np #
