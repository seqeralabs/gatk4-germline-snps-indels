nextflow.enable.dsl = 2

params.make_gvcf = false
params.memory = '10'


process GATK_MERGE_GVCFS {
    container = "broadinstitute/gatk:4.1.8.1"
    memory "${params.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(input_vcf), path(input_vcf_index)


    output:
    tuple path("${output_filename}"), path("${output_filename}.tbi")

    script:
    input_str = input_vcf.reduce("") { list_of_vcfs, new_vcf -> list_of_vcfs + " --INPUT ${new_vcf}" }.join(' ')
    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = input_vcf.getBaseName() + output_suffix

    """
    set -e

    /gatk/gatk --java-options "-Xmx${$params.memory}G"  \
        MergeVcfs \
        ${input_str} \
        --OUTPUT ${output_filename}
    
    """
}


workflow test {
    input_vcf_ch = Channel.value(["./test_data/*vcf", "./test_data/*tbi"])
    GATK_MERGE_GVCFS(input_vcf_ch)

}
