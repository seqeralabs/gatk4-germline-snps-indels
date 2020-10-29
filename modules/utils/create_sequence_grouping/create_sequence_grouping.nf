nextflow.enable.dsl = 2

params.java_opts = ""

process UTILS_CREATE_SEQUENCE_GROUPING {
    tag { "Create sequence grouping" }

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory "16 GB"
    cpus 16

    input:
    path(ref_dict)

    output:
    path("*txt")
    path("*unmapped.txt")

    script:
    """
    sequenceGrouping.py ${ref_dict} "sequence_grouping.txt" "sequence_grouping_with_unmapped.txt"
    """
}
