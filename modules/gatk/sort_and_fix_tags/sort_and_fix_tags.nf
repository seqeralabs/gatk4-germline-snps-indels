params.gatk_path = "gatk"
params.java_opts_sort = "-Xms4000m"
params.java_opts_fix = "-Xms500m"
params.compression_level = 5


process GATK_SORT_AND_FIX_TAGS {
    tag "${sampleId}"
    label 'gatk4_container'

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


    stub:

    """
    touch ${sampleId}.mapped.merged.duplicate_marked.sorted.bai
    touch ${sampleId}.mapped.merged.duplicate_marked.sorted.bam
    touch ${sampleId}.mapped.merged.duplicate_marked.sorted.bam.md5
    
    """
}
