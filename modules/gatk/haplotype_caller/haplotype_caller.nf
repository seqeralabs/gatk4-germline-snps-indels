
process haplotypeCaller {
    tag "${sampleId}_${interval_chunk_name}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

    input:

    tuple val(sampleId),
            path(input_recal_merged_bam),
            path(input_recal_merged_bai),
            path(input_recal_merged_md5),
            val(scatter_id),
            val(interval_chunk_name),
            path(interval_list_file)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)


    output:

    tuple val(sampleId),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf"),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf.idx")


    script:

    gatk_path = "/gatk/gatk"

    """
    set -e
    ${gatk_path} --java-options "${params.java_opt_haplotype}" \
    HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${input_recal_merged_bam} \
    --output "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf" \
    -contamination 0 \
    -ERC GVCF \
    -L ${interval_list_file}
    """
}

