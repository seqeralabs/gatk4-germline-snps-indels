process markDuplicates {
    tag "${sampleId}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

    input:
    val(sampleId)
    path(input_mapped_merged_bam)

    output:
    val(sampleId)
    path("*_merged.deduped.bam")
    path("*_merged.deduped.metrics.txt")

    script:
    gatk_path = "/gatk/gatk"

    """
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_markdups}" \
    MarkDuplicates \
    --INPUT ${input_mapped_merged_bam} \
    --OUTPUT ${sampleId}_merged.deduped.bam \
    --METRICS_FILE ${sampleId}_merged.deduped.metrics.txt \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER "queryname" \
    --CREATE_MD5_FILE true
    """

}
