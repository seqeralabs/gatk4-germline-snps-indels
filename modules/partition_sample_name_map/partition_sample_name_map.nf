nextflow.enable.dsl = 2

params.line_limit

process UTILS_PARTITION_SAMPLE_NAME_MAP {
    container "us.gcr.io/broad-gotc-prod/python:2.7"

    input:
    path(sample_name_map)

    output:
    path("partition_*")

    shell:
    '''
    cut -f 2 !{sample_name_map} > sample_paths
    split -l !{params.line_limit} -d sample_paths partition_

    # Let the OS catch up with creation of files for glob command
    sleep 1
    '''

}
