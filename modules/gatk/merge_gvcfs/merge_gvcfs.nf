nextflow.preview.dsl = 2

params.make_gvcf = false
params.memory = '10'
params.gatk_path

process GATK_MERGE_GVCFS {
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    memory "${params.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(input_vcf), path(input_vcf_index)


    output:
    tuple path("${output_filename}"), path("${output_filename}.tbi")

    script:
    input_str = input_vcf.reduce("") { a, b -> a + " --INPUT ${b}" }.join(' ')
    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = input_vcf.getBaseName() + output_suffix

    """
    set -e

    ${gatk_path} --java-options "-Xmx${$params.memory}G"  \
      MergeVcfs \
      ${input_str} \
      --OUTPUT ${output_filename}
    
    """
}
