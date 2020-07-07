rule meta_assembly:
  input:
    FORWARD_POOLED="results/02_filter_host/{sample}/{sample}_FP.fastq.gz",
    REVERSE_POOLED="results/02_filter_host/{sample}/{sample}_RP.fastq.gz"

  output:
    "results/04_meta-assembly/{sample}/scaffolds.fasta"

  params:
    spades=config["spades"]["spades_version"]

  threads:
    int(config['spades']['spades_threads'])

  resources:
    mem_mb = int(config['spades']['spades_mem_mb']),
    hours = int(config['spades']['spades_hours']),
    mem_gb = int(config['spades']['spades_mem_gb'])

  shell:
    " module add vital-it/7 ;"
    " module add UHTS/Assembler/SPAdes/{params.spades} ;"
    " srun spades.py "
    "  --meta " # Metagenomics flag
    "  -t {threads} "
    "  -m {resources.mem_gb} "
    "  -1 {input.FORWARD_POOLED} "
    "  -2 {input.REVERSE_POOLED} "
    "  -o results/04_meta-assembly/{wildcards.sample}/ ;"
