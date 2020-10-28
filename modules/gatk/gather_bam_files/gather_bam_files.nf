
process gatherBamFiles {
    tag "${sampleId}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

    input:
    tuple val(sampleId), path(input_recalibrated_bams)

    output:
    tuple   val(sampleId),
            path("${sampleId}.recal.merged.bam"),
            path("${sampleId}.recal.merged.bai"),
            path("${sampleId}.recal.merged.bam.md5")

    script:

    gatk_path = "/gatk/gatk"
    inputs_bams_to_merge = input_recalibrated_bams
            .sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --INPUT ")

    """
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level}  \
    ${params.java_opt_gatherbams} " \
    GatherBamFiles \
    --INPUT ${inputs_bams_to_merge} \
    --OUTPUT "${sampleId}.recal.merged.bam" \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true
    """
}
