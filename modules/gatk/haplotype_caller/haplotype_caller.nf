nextflow.preview.dsl = 2

params.contamination
params.make_gvcf = false
params.make_bamout
params.hc_scatter
params.gcs_project_for_requester_pays = false


process GATK_HAPLOTYPE_CALLER {
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"

    input:
    tuple path(input_bam), path(input_bam_index)
    tuple path(ref_dict), path(ref_fasta), path(ref_fasta_index)
    path(interval_list)

    output:
    tuple path("${output_file_name}"), path("${output_file_name}.tbi")
    path("${vcf_basename}.bamout.bam") optional true

    script:
    // TODO: Derive these using the file extension utils in NXF
    vcf_basename
    bamout_arg

    """
    set -e

    ${gatk_path} --java-options "-Xmx${command_mem_gb}G ${java_opt}" \
          HaplotypeCaller \
          -R ${ref_fasta} \
          -I ${input_bam} \
          -L ${interval_list} \
          -O ${output_filename} \
          -contamination ${params.contamination ? params.contamination : 0} \
          -G StandardAnnotation \
          -G StandardHCAnnotation \
          ${params.make_gvcf ? "-G AS_StandardAnnotation" : ""} \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
          ${params.make_gvcf ? "-ERC GVCF" : ""} \
          ${params.gcs_project_for_requester_pays ? "--gcs-project-for-requester-pays ${params.gcs_project_for_requester_pays}" : ""} \
          ${bamout_arg}
    """
}
