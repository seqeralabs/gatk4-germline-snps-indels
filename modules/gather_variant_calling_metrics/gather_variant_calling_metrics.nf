nextflow.enable.dsl = 2

params.output_prefix

process GATK_GATHER_VARIANT_CALLING_METRICS {
    container "us.gcr.io/broad-gotc-prod/gatk4-joint-genotyping:1.3.0-1527875152"

    input:
    path(input_details)
    path(input_summaries)


    output:
    path("${params.output_prefix}.variant_calling_detail_metrics")
    path("${params.output_prefix}.variant_calling_summary_metrics")

    shell:

    '''
    set -euo pipefail

    input_details_fofn=!{write_lines(input_details)}
    input_summaries_fofn=!{write_lines(input_summaries)}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat $input_details_fofn | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
    sleep 1
    ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
    echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat $input_summaries_fofn | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
    sleep 1
    ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
    echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=$(cat $input_details_fofn | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("--INPUT metrics/%s ", $1)}')

    /usr/gitc/gatk --java-options -Xms2g AccumulateVariantCallingMetrics \
                                         $INPUT \
                                         --OUTPUT ~{output_prefix}
    '''

}
