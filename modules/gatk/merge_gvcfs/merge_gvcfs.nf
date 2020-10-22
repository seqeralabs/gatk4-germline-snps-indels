nextflow.enable.dsl = 2



params_map = [
        make_gvcf: params.gatk_merge_gvcfs_make_gvcf = true,
        memory: params.gatk_merge_gvcfs_memory = '10'
]


process GATK_MERGE_GVCFS {
    container = "broadinstitute/gatk:4.1.8.1"
    memory "${params_map.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(input_vcf), path(input_vcf_index)


    output:
    tuple path("${output_filename}"), path("${output_filename}.tbi")

    script:
    input_str = input_vcf
            .reduce("") { list_of_vcfs, new_vcf -> list_of_vcfs + " --INPUT ${new_vcf}" }
            .join(' ')
    output_suffix = params_map.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = input_vcf.getBaseName() + output_suffix

    """
    set -e

    /gatk/gatk --java-options "-Xmx${$params_map.memory}G"  \
          MergeVcfs \
          ${input_str} \
          --OUTPUT ${output_filename}
    
    """
}


workflow test {
    input_vcf_ch = Channel.fromPath(["${baseDir}/test_data/*.vcf.gz",
                                     "${baseDir}/test_data/*.tbi"]).buffer(size: 2)

    GATK_MERGE_GVCFS(input_vcf_ch)

}
