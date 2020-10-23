nextflow.enable.dsl = 2

params.output_file_name

process UTILS_GATHER_PICARD_METRICS {
    container "us.gcr.io/broad-gotc-prod/python:2.7"

    input:
    path(metrics_files)


    output:
    path("${params.output_file_name}")


    shell:

    metrics_files_str = {sep=' ' metrics_files}

    '''
    # Don't use this task to gather tens of thousands of files.
    # Cromwell can't handle it.

    # This cannot gather metrics with histograms

    head -n 7 !{metrics_files[0]} > !{output_file_name}

    for metrics_file in !{metrics_files_str}; do
        sed -n '1,7d;p' $metrics_file | grep -v '^$' >> !{output_file_name}
    done
    '''


}
