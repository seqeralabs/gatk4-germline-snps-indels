process GET_BWA_VERSION {
    tag "BWA version"

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"

    output:
    stdout

    script:

    bwa_path = "/usr/gitc"
    """

    ${bwa_path}/bwa 2>&1  \
    | grep -e '^Version'  \
    | sed 's/Version: //' \
    | tr -d '\n'
    """
}
