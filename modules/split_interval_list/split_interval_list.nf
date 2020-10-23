nextflow.enable.dsl = 2

params.localization_optional = true
params.scatter_count = 2
params.scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"


process GATK_SPLIT_INTERVAL_LIST {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    path(interval_list)

    val(sample_names_unique_done)


    output:
    path("scatterDir")

    script:

    """
    gatk --java-options -Xms3g SplitIntervals \
                                -L ${interval_list} \
                                -O scatterDir \
                                -scatter ${params.scatter_count} \
                                -R ${ref_fasta} \
                                -mode ${params.scatter_mode} \
                                --interval-merging-rule OVERLAPPING_ONLY
    """

}


