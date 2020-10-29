nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
params.java_opts = ""


process GATK_MERGE_VCFS {
    tag "${sampleId}"

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus


    input:

    tuple val(sampleId),
            path(input_vcfs_to_merge),
            path(inputs_vcf_indices)


    output:
    tuple   val(sampleId),
            path("${sampleId}.merged.vcf")

    script:

    input_vcfs_params = input_vcfs_to_merge
            .sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --INPUT ")

    """
    set -e

    ${params.gatk_path} --java-options "${params.java_opts}"  \
                        MergeVcfs \
                        --INPUT ${input_vcfs_params} \
                        --OUTPUT ${sampleId}.merged.vcf
    """
}

workflow test {
   
}
