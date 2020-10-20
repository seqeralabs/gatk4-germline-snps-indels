// module parameters
params.sample_num_threshold = 50

process CHECK_SAMPLES_UNIQUE {
    container "us.gcr.io/broad-gotc-prod/python:2.7"

    // NOTE: WDL preemptible is not applicable. See https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    // NOTE: We could extract all of the following directives as parameters
    memory '1 GB'
    disk '10 GB'

    input:
    file(sample_name_map)


    output:
    // NOTE: In WDL sample_names_unique_done serves the same purpose. See CheckSamplesUnique and SplitIntervalList
    val true


    shell:
    '''
    set -euo pipefail
    if [[ $(cut -f 1 !{sample_name_map} | wc -l) -ne $(cut -f 1 !{sample_name_map} | sort | uniq | wc -l) ]]
    then
      echo "Samples in the sample_name_map are not unique" 1>&2
      exit 1
    elif [[ $(cut -f 1 !{sample_name_map} | wc -l) -lt !{params.sample_num_threshold} ]]
    then
      echo "There are fewer than !{params.sample_num_threshold} samples in the sample_name_map" 1>&2
      echo "Having fewer than !{params.sample_num_threshold} samples means there likely isn't enough data to complete joint calling" 1>&2
      exit 1
    else
      echo true
    fi
    '''

}
