rule create_reference:
  input:
    reference = config["host_reference"],
    link = "results/01_QC/link_rule.txt"

  output:
    indexed_reference = "results/02_filter_host/reference/ref.1.bt2",
    log_ref = "results/02_filter_host/reference/log_ref.txt"

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
    "  -c {input.reference} results/02_filter_host/reference/ref ;"
    " srun /bin/touch {output.log_ref}; "

#-------------------------------------------------------------------------------

rule run_bowtie2:
  input:
    R1= "results/01_QC/{sample}/{sample}_FP.fastq.gz",
    R2= "results/01_QC/{sample}/{sample}_RP.fastq.gz",
    ref = "results/02_filter_host/reference/log_ref.txt"

  output:
    temp("results/02_filter_host/{sample}/{sample}.sam")

  params:
    bowtie2=config['bowtie2']['bowtie2_version']

  threads:
    int(config['bowtie2']['bowtie2_map']['bowtie2_threads'])

  resources:
    mem_mb = int(config['bowtie2']['bowtie2_map']['bowtie2_mem_mb']),
    hours = int(config['bowtie2']['bowtie2_map']['bowtie2_hours'])

  shell:
    " module add UHTS/Aligner/bowtie2/{params.bowtie2} ;"
    " srun bowtie2 -x results/02_filter_host/reference/ref "
    "  -1 {input.R1} "
    "  -2 {input.R2} "
    "  -p {threads} "
    "  -S {output} ;"

#-------------------------------------------------------------------------------

rule sam_to_bam:
  input:
    "results/02_filter_host/{sample}/{sample}.sam"

  output:
    temp("results/02_filter_host/{sample}/{sample}.bam")

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
    "results/02_filter_host/{sample}/{sample}.bam"

  output:
    "results/02_filter_host/{sample}/{sample}.sorted.bam"

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

#-------------------------------------------------------------------------------

rule unmapped_bam_to_fastq:
  input:
    "results/02_filter_host/{sample}/{sample}.sorted.bam"

  output:
    FORWARD="results/02_filter_host/{sample}/{sample}_FP.fastq.gz",
    REVERSE="results/02_filter_host/{sample}/{sample}_RP.fastq.gz"

  params:
    samtools=config["samtools"]["samtools_version"]

  threads:
    int(config['samtools']['samtools_threads'])

  resources:
    mem_mb = int(config['samtools']['samtools_mem_mb']),
    hours = int(config['samtools']['samtools_hours'])

  shell:
    " module add UHTS/Analysis/samtools/{params.samtools} ;"
    " /bin/echo extracting unmapped reads for {wildcards.sample} and writing it to fastq file ;"
    " srun samtools view -@ {threads} -u -f 0x4 {input} | "
    " samtools bam2fq -N -O -@ {threads} "
    "  -0 /dev/null -s /dev/null "
    "  -1 {output.FORWARD} "
    "  -2 {output.REVERSE} - ;"
