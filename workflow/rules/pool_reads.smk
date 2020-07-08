rule pool_reads:
  input:
    FORWARD=expand("results/02_filter_host/{sample}/{sample}_FP.fastq.gz",
                    sample = samples) if config["remove_host"]["remove_host"] \
                    else \
                    expand("results/01_QC/{sample}/{sample}_FP.fastq.gz",
                           sample = samples),
    REVERSE=expand("results/02_filter_host/{sample}/{sample}_RP.fastq.gz",
                    sample = samples) if config["remove_host"]["remove_host"] \
                    else \
                    expand("results/01_QC/{sample}/{sample}_RP.fastq.gz",
                           sample = samples),

  threads:
    int(config['short_sh_commands_threads'])

  resources:
    mem_mb = int(config['spades']['spades_mem_mb']),
    hours = int(config['short_sh_commands_hours'])

  output:
    FORWARD_POOLED="results/03_pool_reads/pooled_reads_FP.fastq.gz",
    REVERSE_POOLED="results/03_pool_reads/pooled_reads_RP.fastq.gz"

  shell:
    "srun /bin/cat {input.FORWARD} >> {output.FORWARD_POOLED} ;"
    "srun /bin/cat {input.REVERSE} >> {output.REVERSE_POOLED} ;"
