nextflow.enable.dsl = 2

params.output_filename

process GATK_GATHER_TRANCHES {
    container "us.gcr.io/broad-gotc-prod/gatk4-joint-genotyping:1.3.0-1527875152"

    input:
    path(tranches)


    output:
    path("${params.output_filename}")


    shell:
//    NOTE: This assumes that gsutil is installed
    tranches_lines = write_lines(tranches)

    '''
    set -euo pipefail

    tranches_fofn=!{tranches_lines}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf tranches
    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat $tranches_fofn | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I tranches/; do
    sleep 1
    ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
    echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

    /usr/gitc/gatk --java-options -Xms6g GatherTranches \
                                        --input inputs.list \
                                        --output !{output_filename}
    '''

}
