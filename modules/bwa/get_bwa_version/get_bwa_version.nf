nextflow.enable.dsl = 2

params.gitc_path = "/usr/gitc"


process BWA_GET_BWA_VERSION {
    tag "BWA version"
    label "gitc_container"

    output:
    stdout

    script:

    """

    ${params.gitc_path}/bwa 2>&1  \
    | grep -e '^Version'  \
    | sed 's/Version: //' \
    | tr -d '\n'
    """

    stub:
    """
    echo 0.7.15-r1140
    """
}
