rule trim:
  input:
    R1 = f"{DataFolder}" + "{sample}" + config['mates']['mate1'] + f"{fastx_extension}",
    R2 = f"{DataFolder}" + "{sample}" + config['mates']['mate2'] + f"{fastx_extension}"

  output:
    FP = "results/01_QC/{sample}/{sample}_FP.fastq.gz",
    RP = "results/01_QC/{sample}/{sample}_RP.fastq.gz",
    FU = "results/01_QC/{sample}/{sample}_FU.fastq.gz",
    RU = "results/01_QC/{sample}/{sample}_RU.fastq.gz"

  params:
    trimmomatic=config['trimmomatic']['trimmomatic_version']

  threads:
    int(config['trimmomatic']['trimmomatic_threads'])

  resources:
    mem_mb = int(config['trimmomatic']['trimmomatic_mem_mb']),
    hours = int(config['trimmomatic']['trimmomatic_hours'])

  shell:
    " module add UHTS/Analysis/trimmomatic/{params.trimmomatic} ;"
    " srun trimmomatic PE -threads {threads} "
    " -phred33 "
    " {input.R1} {input.R2} "
    " {output.FP} {output.FU} "
    " {output.RP} {output.RU} "
    " ILLUMINACLIP:resources/nextera_adapters.fasta:2:30:10  "
    " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30;"

#-------------------------------------------------------------------------------

rule fastqc_run:
  input:
    R1="results/01_QC/{sample}/{sample}_FP.fastq.gz",
    R2="results/01_QC/{sample}/{sample}_RP.fastq.gz"

  output:
    zip1 = "results/01_QC/result_fastqc/{sample}_FP_fastqc.zip",
    zip2 = "results/01_QC/result_fastqc/{sample}_RP_fastqc.zip",
    link_rule = "results/01_QC/result_fastqc/{sample}_RP_fastqc.txt"

  params:
    fastqc=config['fastqc']['fastqc_version']

  threads:
    int(config['fastqc']['fastqc_threads'])

  resources:
    mem_mb = int(config['fastqc']['fastqc_mem_mb']),
    hours = int(config['fastqc']['fastqc_hours'])

  shell:
    " module add UHTS/Quality_control/fastqc/{params.fastqc}; "
    "srun fastqc {input.R1} -t {threads} "
    "  -o results/01_QC/result_fastqc/; "
    "srun fastqc {input.R2} -t {threads} "
    " -o results/01_QC/result_fastqc/; "
    "srun /bin/echo FINISHED {wildcards.sample} > {output.link_rule}; "

#-------------------------------------------------------------------------------

rule create_link_file:
  input:
    expand("results/01_QC/result_fastqc/{sample}_RP_fastqc.txt",
            sample = samples)

  params:
    trimmomatic=config['trimmomatic']['trimmomatic_version']

  threads:
    int(config['short_sh_commands_threads'])

  resources:
    mem_mb = 16000,
    hours = int(config['short_sh_commands_hours'])

  output:
    "results/01_QC/link_rule.txt"

  shell:
    "srun /bin/cat {input} > {output}"
