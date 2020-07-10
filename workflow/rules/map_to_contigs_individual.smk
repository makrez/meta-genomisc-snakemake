rule create_reference:
  input:
    link = "results/04_meta-assembly/{sample}_assembly_finished.txt",

  output:
    indexed_reference = "results/05_map_to_contigs/{sample}/contigs.fasta.ann",
    log_ref = "results/05_map_to_contigs/{sample}/log_ref.txt"

  params:
    bwa=config['bwa']['bwa_version']

  threads:
    int(config['bwa']['bwa_threads'])

  resources:
    mem_mb = int(config['bwa']['bwa_mem_mb']),
    hours = int(config['bwa']['bwa_hours'])

  shell:
    " module add UHTS/Aligner/bwa/{params.bwa} ;"
    " awk '/^>/ {{printf(\"\\n%s\\n\",$0);next; }} {{ printf(\"%s\",$0);}}  "
    "  END {{printf(\"\\n\");}}' "
    "  < results/04_meta-assembly/spades/{wildcards.sample}/contigs.fasta | "
    "  sed '1{{/^$/d}}' "
    "  > results/05_map_to_contigs/{wildcards.sample}/contigs.fasta  ;"
    " srun bwa index results/05_map_to_contigs/{wildcards.sample}/contigs.fasta ;"
    " srun /bin/touch {output.log_ref}; "

#-------------------------------------------------------------------------------

rule bwa:
  input:
    R1= "results/01_QC/{sample}/{sample}_FP.fastq.gz",
    R2= "results/01_QC/{sample}/{sample}_RP.fastq.gz",
    ref = "results/05_map_to_contigs/{sample}/contigs.fasta.ann",

  output:
    "results/05_map_to_contigs/{sample}/{sample}.sam"

  params:
    bwa=config['bwa']['bwa_version']

  threads:
    int(config['bwa']['bwa_threads'])

  resources:
    mem_mb = int(config['bwa']['bwa_mem_mb']),
    hours = int(config['bwa']['bwa_hours'])

  shell:
    " module add UHTS/Aligner/bwa/{params.bwa} ;"
    " srun bwa mem "
    "  -t {threads} "
    "  results/05_map_to_contigs/{wildcards.sample}/contigs.fasta "
    "  {input.R1} "
    "  {input.R2} "
    "  > {output} ;"

#-------------------------------------------------------------------------------

rule sam_to_bam:
  input:
    "results/05_map_to_contigs/{sample}/{sample}.sam"

  output:
    "results/05_map_to_contigs/{sample}/{sample}.bam"

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
    " srun samtools view -@ {threads} -b -o {output} {input} ;"

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
