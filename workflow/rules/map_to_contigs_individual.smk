rule create_reference:
  input:
    link = "results/04_meta-assembly/{sample}_assembly_finished.txt",

  output:
    indexed_reference = "results/04_meta-assembly/spades/{sample}/ref.1.bt2",
    log_ref = "results/04_meta-assembly/spades/{sample}/log_ref.txt"

  params:
    bowtie2=config['bowtie2']['bowtie2_version']

  threads:
    int(config['bowtie2']['bowtie2_ref']['bowtie2_threads'])

  resources:
    mem_mb = int(config['bowtie2']['bowtie2_ref']['bowtie2_mem_mb']),
    hours = int(config['bowtie2']['bowtie2_ref']['bowtie2_hours'])

  shell:
    " module add UHTS/Aligner/bowtie2/{params.bowtie2} ;"
    " srun bowtie2-build --threads {threads} "
    "  -c results/04_meta-assembly/spades/{wildcards.sample}/contigs.fasta "
    "  results/02_filter_host/reference/ref ;"
    " srun /bin/touch {output.log_ref}; "

#-------------------------------------------------------------------------------

rule run_bowtie2:
  input:
    R1= "results/01_QC/{sample}/{sample}_FP.fastq.gz",
    R2= "results/01_QC/{sample}/{sample}_RP.fastq.gz",
    ref =  "results/04_meta-assembly/spades/{sample}/ref.1.bt2",

  output:
    temp("results/05_map_to_contigs/{sample}/{sample}.sam")

  params:
    bowtie2=config['bowtie2']['bowtie2_version']

  threads:
    int(config['bowtie2']['bowtie2_map']['bowtie2_threads'])

  resources:
    mem_mb = int(config['bowtie2']['bowtie2_map']['bowtie2_mem_mb']),
    hours = int(config['bowtie2']['bowtie2_map']['bowtie2_hours'])

  shell:
    " module add UHTS/Aligner/bowtie2/{params.bowtie2} ;"
    " srun bowtie2 -x results/02_filter_host/reference/ref " #//# TODO: change ref
    "  -1 {input.R1} "
    "  -2 {input.R2} "
    "  -p {threads} "
    "  -S {output} ;"

#-------------------------------------------------------------------------------

rule sam_to_bam:
  input:
    "results/05_map_to_contigs/{sample}/{sample}.sam"

  output:
    temp("results/05_map_to_contigs/{sample}/{sample}.bam")

  params:
    samtools=config["samtools"]["samtools_version"]

  threads:
    int(config['samtools']['samtools_threads'])

  resources:
    mem_mb = int(config['samtools']['samtools_mem_mb']),
    hours = int(config['samtools']['samtools_hours'])

  shell:
    " /bin/echo 'Running SAM to BAM' ;" # TODO: Why does this not work?
    " module add UHTS/Analysis/samtools/{params.samtools} ;"
    " srun samtools view -@ {threads} -b -S -o {output} {input} ;"

#-------------------------------------------------------------------------------

rule sort_bam:
  input:
    "results/05_map_to_contigs/{sample}/{sample}.bam"

  output:
    "results/05_map_to_contigs/{sample}/{sample}.sorted.bam"

  params:
    samtools=config["samtools"]["samtools_version"]

  threads:
    int(config['samtools']['samtools_threads'])

  resources:
    mem_mb = int(config['samtools']['samtools_mem_mb']),
    hours = int(config['samtools']['samtools_hours'])

  shell:
    " module add UHTS/Analysis/samtools/{params.samtools} ;"
    " /bin/echo Sorting BAM ;"
    " srun samtools sort -@ {threads} -o {output} {input}; "
