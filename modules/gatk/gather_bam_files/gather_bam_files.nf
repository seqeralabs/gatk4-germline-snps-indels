nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
params.java_opts = ""
params.compression_level = 5

process GATK_GATHER_BAM_FILES {
    tag "${sampleId}"

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus

    input:
    tuple val(sampleId), path(input_recalibrated_bams)

    output:
    tuple val(sampleId),
            path("${sampleId}.recal.merged.bam"),
            path("${sampleId}.recal.merged.bai"),
            path("${sampleId}.recal.merged.bam.md5")

    script:

    inputs_bams_to_merge = input_recalibrated_bams
            .sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --INPUT ")

    """
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level}  ${params.java_opts} " \
                 GatherBamFiles \
                --INPUT ${inputs_bams_to_merge} \
                --OUTPUT "${sampleId}.recal.merged.bam" \
                --CREATE_INDEX true \
                --CREATE_MD5_FILE true
    """
}

workflow test {

}
