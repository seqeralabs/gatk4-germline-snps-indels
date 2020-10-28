process createSequenceGrouping {
    tag { "Create sequence grouping" }

    memory '16 GB'
    cpus 16

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"

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
