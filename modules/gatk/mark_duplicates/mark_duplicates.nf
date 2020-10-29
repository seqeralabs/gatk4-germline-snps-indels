nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
params.java_opts = ""
params.compression_level = 5

process GATK_MARK_DUPLICATES {
    tag "${sampleId}"

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus

    input:
    val(sampleId)
    path(input_mapped_merged_bam)

    output:
    val(sampleId)
    path("*_merged.deduped.bam")
    path("*_merged.deduped.metrics.txt")

    script:

    """
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opts}" \
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

workflow test {
   
}
