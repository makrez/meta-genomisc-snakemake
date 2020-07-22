rule meta_assembly:
  input:
    FORWARD="results/03_pool_reads/pooled_reads_FP.fastq.gz",
    REVERSE="results/03_pool_reads/pooled_reads_RP.fastq.gz"

  output:
    "results/04_meta-assembly/pooled_assembly_finished.txt"

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['meta_assembly']['meta_assembly_threads'])

  resources:
    mem_mb = int(config['meta_assembly']['meta_assembly_mem_mb']),
    hours = int(config['meta_assembly']['meta_assembly_hours']),
    mem_gb = int(config['meta_assembly']['meta_assembly_mem_gb'])

  run:
    if (config["assembler"] == "spades"):
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
        "  -o results/04_meta-assembly/spades ;"
        "srun /bin/touch {output} ;"
       )
    elif (config["assembler"] == "megahit"):
        shell(
          " mkdir -p results/04_meta-assembly/;"
          " module add UHTS/Assembler/megahit/1.1.4 ;"
          " srun megahit "
          "  -1 {input.FORWARD} "
          "  -2 {input.REVERSE} "
          "  -o results/04_meta-assembly/megahit ;"
          " /bin/touch {output} ; "
         )
    else:
        print("Assembler not recognized. \
              Use either the strings 'spades' or 'megahit' in the config file")
