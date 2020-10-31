nextflow.enable.dsl = 2

params.gitc_path = "/usr/gitc"


process BWA_GET_BWA_VERSION {
    tag "BWA version"

    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"

    output:
    stdout

    script:

    """

    ${params.gitc_path}/bwa 2>&1  \
    | grep -e '^Version'  \
    | sed 's/Version: //' \
    | tr -d '\n'
    """
}
