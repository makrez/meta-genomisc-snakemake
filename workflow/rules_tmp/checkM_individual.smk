rule checkM:
  input:
    LINK = "results/06_metabat/{sample}_metabat_finished.txt"

  output:
    QUAL= directory("results/08_checkM/{sample}/"),
    OUT = "results/08_checkM/{sample}/{sample}.link",
    LINK = "results/08_checkM/{sample}.link",
    RESULT = "results/08_checkM/{sample}/{sample}_result.txt",

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
    "  --file {output.RESULT} "
    "  --tab_table "
    "  results/06_metabat/{wildcards.sample}/ "
    "  {output.QUAL} > {output.OUT} ;"
    " /bin/touch {output.LINK} ;"

#-------------------------------------------------------------------------------

rule parse_output:
  input:
    expand("results/08_checkM/{sample}/{sample}_result.txt", sample= samples)

  output:
    "results/08_checkM/checkm_summary.txt"

  threads:
    int(config['short_commands_threads'])

  resources:
    mem_mb = int(config['short_commands_mb']),
    hours = int(config['short_commands_hours'])

  shell:
    " /bin/cat {input} > {output} "
