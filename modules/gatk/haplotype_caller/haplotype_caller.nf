//===========================
// Module parameters
//===========================

params.contamination
params.make_gvcf
params.make_bamout
params.hc_scatter
params.gcs_project_for_requester_pays = false

//===========================
// Process definition
//===========================

process GATK_HAPLOTYPE_CALLER {

//---------------------------
//  directives
//---------------------------


//---------------------------
    input:
//---------------------------
    path(input_bam)
    path(input_bam_index)
    path(interval_list)
    val(output_file_name)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_index)


//---------------------------
    output:
//---------------------------
    path("${output_file_name}")
    path("${output_file_name}.tbi")
    // NOTE: In NXF we can specify an optional output. See Line-227 of haplotypeCaller_gvcf_gatk4.wdl
    path("${vcf_basename}.bamout.bam") optional true

//---------------------------
    script:
//---------------------------
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
