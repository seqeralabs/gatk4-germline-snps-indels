nextflow.enable.dsl = 2

params.contamination
params.make_gvcf = false
params.make_bamout = false
params.memory = '10'
params.java_opts = ''
params.hc_scatter

process GATK_HAPLOTYPE_CALLER {
    container = "broadinstitute/gatk:4.1.8.1"
    memory "${params.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(input_bam), path(input_bam_index)
    tuple path(ref_dict), path(ref_fasta), path(ref_fasta_index)
    path(interval_list)

    output:
    tuple path("${output_file_name}"), path("${output_file_name}.tbi")
    path("${vcf_basename}.bamout.bam") optional true

    script:
    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = input_bam.getBaseName() + output_suffix
    bamout_arg = make_bamout ? "-bamout ~{vcf_basename}.bamout.bam" : ""

    """
    set -e

    /gatk/gatk --java-options "-Xmx${params.memory}G ${params.java_opts}" \
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
          ${bamout_arg}
    """
}
