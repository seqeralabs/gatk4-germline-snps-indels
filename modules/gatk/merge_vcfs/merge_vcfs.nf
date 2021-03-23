nextflow.enable.dsl = 2

params.gatk_path = "gatk"
params.java_opts = "-Xmx4G"


process GATK_MERGE_VCFS {
    tag "${sampleId}"
    memory "16 GB"
    cpus 16


    input:
    tuple val(sampleId),
            path(input_vcfs_to_merge),
            path(inputs_vcf_indices)


    output:
    tuple val(sampleId), path("${sampleId}.merged.vcf")

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

    stub:

    """
    touch ${sampleId}.merged.vcf
    """

}

