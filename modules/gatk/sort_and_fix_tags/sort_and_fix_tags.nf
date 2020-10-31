nextflow.enable.dsl = 2

params.gatk_path = "/gatk/gatk"
params.java_opts_sort = ""
params.java_opts_fix = ""
params.compression_level = 5


process GATK_SORT_AND_FIX_TAGS {
    tag "${sampleId}"

    container = "broadinstitute/gatk:4.1.8.1"
    memory "16 GB"
    cpus 16

    input:
    val(sampleId)
    path(input_mapped_merged_marked_bam)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    output:
    tuple val(sampleId),
            path("${sampleId}.mapped.merged.duplicate_marked.sorted.bam"),
            path("${sampleId}.mapped.merged.duplicate_marked.sorted.bai"),
            path("${sampleId}.mapped.merged.duplicate_marked.sorted.bam.md5")


    script:

    """
    set -o pipefail

    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opts_sort}" \
     SortSam \
        --INPUT ${input_mapped_merged_marked_bam} \
        --OUTPUT /dev/stdout \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX false \
        --CREATE_MD5_FILE false \
    | \
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opts_fix}" \
    SetNmMdAndUqTags \
        --INPUT /dev/stdin \
        --OUTPUT ${sampleId}.mapped.merged.duplicate_marked.sorted.bam \
        --CREATE_INDEX true \
        --CREATE_MD5_FILE true \
        --REFERENCE_SEQUENCE ${ref_fasta}
    """
}
