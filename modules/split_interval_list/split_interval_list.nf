// module parameters
params.localization_optional = true
// NOTE: scatter_count needs to be provided at the call site. See Line-130-131 in JointGenotyping.wdl
params.scatter_count = 2
params.scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
params.gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"

process SPLIT_INTERVAL_LIST {
    // NOTE: We could extract all of the following directives as parameters
    memory '3.75 GB'
    disk '10 GB'

    input:
    file(interval_list)
    file(ref_fasta)
    file(ref_fasta_index)
    file(ref_dict)
    // NOTE: This value comes from a successful completion of CHECK_SAMPLES_UNIQUE process
    val(sample_names_unique_done)


    output:
    // NOTE: Instead of passing an Array[File], as in WDL, we can simply put the dir path in the channel
    path("scatterDir")

    script:

    """
    gatk --java-options -Xms3g SplitIntervals \
                                -L ${interval_list} \
                                -O scatterDir \
                                -scatter ${scatter_count} \
                                -R ${ref_fasta} \
                                -mode ${scatter_mode} \
                                --interval-merging-rule OVERLAPPING_ONLY
    """

}
