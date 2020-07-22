rule meta_assembly:
  input:
    FORWARD="results/02_filter_host/{sample}/{sample}_FP.fastq.gz" \
            if config["remove_host"]["remove_host"] else \
            "results/01_QC/{sample}/{sample}_FP.fastq.gz",
    REVERSE="results/02_filter_host/{sample}/{sample}_RP.fastq.gz" \
            if config["remove_host"]["remove_host"] else \
            "results/01_QC/{sample}/{sample}_RP.fastq.gz",

  output:
    "results/04_meta-assembly/{sample}_assembly_finished.txt"

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['spades']['spades_threads'])

  resources:
    mem_mb = int(config['spades']['spades_mem_mb']),
    hours = int(config['spades']['spades_hours']),
    mem_gb = int(config['spades']['spades_mem_gb'])

  run:
    if (config["assembler"] == "spades"):
    # Spades Assembler will be run
      shell(
        " set +u ;"
        " source {params.conda_profile} ;"
        " conda activate spades ;"
        " srun spades.py "
        "  --meta "
        "  -t {threads} "
        "  -m {resources.mem_gb} "
        "  -1 {input.FORWARD} "
        "  -2 {input.REVERSE} "
        "  -o results/04_meta-assembly/spades/{wildcards.sample} ;"
        "srun /bin/touch {output} ;"
       )
    elif (config["assembler"] == "megahit"):
        # Megahit assembler will be run
        shell(
          " mkdir -p results/04_meta-assembly/megahit/ ;"
          " module add UHTS/Assembler/megahit/1.1.4 ;"
          " srun megahit "
          "  -1 {input.FORWARD} "
          "  -2 {input.REVERSE} "
          "  -o results/04_meta-assembly/megahit/{wildcards.sample} ;"
          "  /bin/touch {output} ; "
         )
    else:
        print("Assembler not recognized. \
              Use either the strings 'spades' or 'megahit' in the config file")
