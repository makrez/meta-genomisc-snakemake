rule meta_assembly:
  input:
    FORWARD_POOLED="analysis/03_pool_reads/pooled_reads_FP.fastq.gz",
    REVERSE_POOLED="analysis/03_pool_reads/pooled_reads_RP.fastq.gz"

  output:
    "analysis/04_meta-assembly/scaffolds.fasta"

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
    " srun spades.py --continue "
    " -o analysis/04_meta-assembly/ ;"
    " #srun spades.py "
    "  #--meta " # Metagenomics flag
    "  #-t {threads} "
    "  #-m {resources.mem_gb} "
    "  #-1 {input.FORWARD_POOLED} "
    "  #-2 {input.REVERSE_POOLED} "
    "  #-o analysis/04_meta-assembly/ ;"
