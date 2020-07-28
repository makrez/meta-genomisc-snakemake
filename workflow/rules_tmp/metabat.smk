rule contig_depth:
  input:
    BAM = "results/05_map_to_contigs/{sample}/{sample}.sorted.bam",

  output:
    DEPTH = "results/06_metabat/{sample}_depth.txt",

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['metabat2']['metabat2_threads'])

  resources:
    mem_mb = int(config['metabat2']['metabat2_mem_mb']),
    hours = int(config['metabat2']['metabat2_hours'])

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate metabat2 ;"
    " srun jgi_summarize_bam_contig_depths "
    "  --outputDepth {output.DEPTH} {input.BAM} ;"

#-------------------------------------------------------------------------------

checkpoint metabat2:
  input:
    DEPTH = "results/06_metabat/{sample}_depth.txt",

  output:
    #LINK = "results/06_metabat/{sample}_metabat_finished.txt",
    FILES = directory("results/06_metabat/{sample}")

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['metabat2']['metabat2_threads'])

  resources:
    mem_mb = int(config['metabat2']['metabat2_mem_mb']),
    hours = int(config['metabat2']['metabat2_hours'])

  run:
    if (config["assembler"] == "spades"):
      shell(
        " set +u ;"
        " source {params.conda_profile} ;"
        " conda activate metabat2 ;"
        " srun metabat2 "
        "  -i results/04_meta-assembly/spades/contigs.fasta "
        "  -a {input.DEPTH} "
        "  -t {threads} "
        "  -o results/06_metabat/{wildcards.sample}/{wildcards.sample} ;"
        " #/bin/touch {output} ;"
        )

    if (config["assembler"] == "megahit"):
      shell(
        " set +u ;"
        " source {params.conda_profile} ;"
        " conda activate metabat2 ;"
        " srun metabat2 "
        "  -i results/04_meta-assembly/megahit/final.contigs.fa "
        "  -a {input.DEPTH} "
        "  -t {threads} "
        "  -o results/06_metabat/{wildcards.sample}/{wildcards.sample} ;"
        " #/bin/touch {output} ;"
        )

#-------------------------------------------------------------------------------
