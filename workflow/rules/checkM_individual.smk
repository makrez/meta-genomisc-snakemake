rule checkM:
  input:
    LINK = "results/07_prokka/summary/{sample}_summary.csv"

  output:
    QUAL= directory("results/08_checkM/{sample}/"),
    OUT = "results/08_checkM/{sample}/{sample}.link",
    LINK = "results/08_checkM/{sample}.link"

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",
    version = "results/10_report/conda_software_versions.txt" #//# TODO: Change numbering for report

  threads:
    int(config['checkM']['checkM_threads'])

  resources:
    mem_mb = int(config['checkM']['checkM_mem_mb']),
    hours = int(config['checkM']['checkM_hours'])

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate checkM ;"
    " checkm | head -n2 | tail -n -1 | sed -E 's/\.|://g' >> {params.version} ;"
    " srun checkm "
    "  lineage_wf "
    "  -t {threads} "
    "  -x fa "
    "  results/06_metabat/{wildcards.sample}/ "
    "  {output.QUAL} > {output.OUT} ;"
    " /bin/touch {output.LINK} ;"

#-------------------------------------------------------------------------------
#
# rule convert_prokka_summary:
#   input:
#     TXT = "results/07_prokka/{sample}/{id}/{sample}_{id}.txt",
#
#   output:
#     CSV = "results/07_prokka/{sample}/{id}/{sample}_{id}.csv"
#
#   threads:
#     int(config['short_sh_commands_threads'])
#
#   resources:
#     mem_mb = int(config['short_commands_mb']),
#     hours = int(config['short_sh_commands_hours'])
#
#   shell:
#     # remove trailing white spaces
#     " srun /bin/sed 's/ *$//g' < {input.TXT} | "
#     "  /bin/sed 's/ /_/g' | /bin/sed 's/:_/ /g' | "
#     "  /bin/awk -v var={wildcards.sample} -v genome={wildcards.id} "
#     "  '{{print $1\",\"$2\",\"var\",\"genome}}' > {output.CSV} ;"

#-------------------------------------------------------------------------------

    # rule move_genomes:
    #   input:
    #     TXT = "results/07_prokka/{sample}/{id}/{sample}_{id}.txt",
    #     GENOME = "results/07_prokka/{sample}/{id}/{sample}_{id}.fna",
    #
    #   output:
    #     DEST = "results/07_prokka/genomes/{sample}/{sample}_{id}.fna"
    #
    #   threads:
    #     int(config['short_sh_commands_threads'])
    #
    #   resources:
    #     mem_mb = int(config['short_commands_mb']),
    #     hours = int(config['short_sh_commands_hours'])
    #
    #   run:
    #     shell(" cp {input.GENOME} {output.DEST}")
    #     # dest_dir = "results/07_prokka/genomes/{wildcards.sample}/"
    #     # print(f'{input.GENOME}')
    #     # shutil.copy(f'{input.GENOME}', dest_dir)

#-------------------------------------------------------------------------------
#
# def aggregate_prokka(wildcards):
#     checkpoint_output = checkpoints.metabat2.get(**wildcards).output[0]
#     ivals = glob_wildcards(os.path.join(checkpoint_output,
#                         "{sample}.{id}.fa")).id
#     print("ivals={}".format(ivals))
#     return expand("results/07_prokka/{sample}/{id}/{sample}_{id}.csv",
#            sample=wildcards,
#            id=ivals)
#
# rule concatenate_prokka:
#   input:
#     aggregate_prokka
#
#   output:
#     "results/07_prokka/summary/{sample}_summary.csv"
#
#   threads:
#     int(config['short_sh_commands_threads'])
#
#   resources:
#     mem_mb = int(config['short_commands_mb']),
#     hours = int(config['short_sh_commands_hours'])
#
#   shell:
#     " srun /bin/cat {input} > {output} "
#
#
# # #-------------------------------------------------------------------------------
# def aggregate_gtdb(wildcards):
#     checkpoint_output = checkpoints.metabat2.get(**wildcards).output[0]
#     ivals = glob_wildcards(os.path.join(checkpoint_output,
#                         "{sample}.{id}.fa")).id
#     print("ivals={}".format(ivals))
#     return expand("results/07_prokka/genomes/{sample}/{sample}_{id}.fna",
#            sample=wildcards,
#            id=ivals)
#
# rule gtdb:
#   input:
#     aggregate_gtdb
#
#   output:
#     LINK = "results/07_prokka/genomes/{sample}/gtdb_link.txt",
#     DIR = directory("results/07_prokka/taxonomy/{sample}/")
#
#   params:
#     conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",
#     version = "results/10_report/conda_software_versions.txt"
#
#   threads:
#     int(config['gtdb']['gtdb_threads'])
#
#   resources:
#     mem_mb = int(config['gtdb']['gtdb_mem_mb']),
#     hours = int(config['gtdb']['gtdb_hours'])
#
#   shell:
#     " set +u ;"
#     " source {params.conda_profile} ;"
#     " conda activate gtdbtk ;"
#     " gtdbtk --version >> {params.version} ;"
#     " srun gtdbtk classify_wf "
#     "  --genome_dir results/07_prokka/genomes/{wildcards.sample}/"
#     "  --out_dir results/07_prokka/taxonomy/{wildcards.sample}/ "
#     "  --cpus {threads} ;"
#     " srun /bin/touch {output.LINK} ;"

#-------------------------------------------------------------------------------
