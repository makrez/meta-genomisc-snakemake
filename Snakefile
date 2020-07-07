import os, string, re

configfile: "config/config.json"

### Config Variables
### ----------------------------------------------------------------------------

DataFolder = config['DataFolder']
fastx_extension = config["extension"]

# Lanes
num_lanes = int(config["Lanes"])
lanes=[]
for i in range(1,num_lanes+1):
    lanes.append(f"L{i}")

# Samples
files=os.listdir(DataFolder)
samples = set([re.sub(r"_L._R._001\.fastq\.gz", "", file) for file in files])

# Full extension: Files are {sample}{extension}
def get_extensions():
    out = []
    for sample in samples:
        for file in files:
            if sample in file:
                tmp = re.sub(sample, "", file)
                out.append(tmp)
    return out

extensions = set(get_extensions())
extensions = list(filter(None, extensions))

# Mate pairs (e.g. "R1" + "R2")
mates = list(set(config["mates"].values()))

# Default Illumina Extension (usually "001")
illumina_ext = config["defaultIlluminaExt"]

wildcard_constraints:
    lane = "|".join(lanes),
    sample = "|".join(samples)

## Check if Reference file exists
if os.path.isfile(config["host_reference"]):
    print("Host reference file exists")
else:
	sys.exit("Host reference missing")

print(config['samtools']['samtools_version'])
print(int(config['bowtie2']['bowtie2_ref']['bowtie2_threads']))
print(int(config['bowtie2']['bowtie2_map']['bowtie2_threads']))


## Pipeline with pooling pool_reads
#------------------------------------------------------------------------------

#rule all:
#    input: "analysis/04_meta-assembly/scaffolds.fasta"

#include: "rules/organise_data"
#include: "rules/QC"
#include: "rules/filter_host"
#include: "rules/pool_reads"
#include: "rules/meta_assembly"

#-------------------------------------------------------------------------------

## Pipeline with individual metagenomics assembly (one meta-assembly per sample)
#-------------------------------------------------------------------------------

rule all:
    input: expand("analysis/04_meta-assembly/{sample}/scaffolds.fasta", sample=samples)
include: "rules/organise_data"
include: "rules/QC"
include: "rules/filter_host"
include: "rules/pool_reads_individual"
include: "rules/meta_assembly_individual"
