rule pool_reads:
  input:
    FORWARD=expand("analysis/02_filter_host/{sample}/{sample}_{lane}_FP.fastq.gz",
                    sample = samples,
                    lane = lanes),
    REVERSE=expand("analysis/02_filter_host/{sample}/{sample}_{lane}_RP.fastq.gz",
                    sample = samples,
                    lane = lanes)
  threads:
    int(config['short_sh_commands_threads'])

  resources:
    mem_mb = int(config['spades']['spades_mem_mb']),
    hours = int(config['short_sh_commands_hours'])

  output:
    FORWARD_POOLED="analysis/03_pool_reads/pooled_reads_FP.fastq.gz",
    REVERSE_POOLED="analysis/03_pool_reads/pooled_reads_RP.fastq.gz"

  shell:
    "srun /bin/cat {input.FORWARD} >> {output.FORWARD_POOLED} ;"
    "srun /bin/cat {input.REVERSE} >> {output.REVERSE_POOLED} ;"
