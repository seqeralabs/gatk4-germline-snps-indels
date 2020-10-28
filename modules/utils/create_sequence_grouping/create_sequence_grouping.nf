nextflow.enable.dsl = 2

params.container = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
params.memory = '16'
params.cpus = 16
// FIXME
params.java_opts = "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"

process UTILS_CREATE_SEQUENCE_GROUPING {
    tag { "Create sequence grouping" }

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus

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
