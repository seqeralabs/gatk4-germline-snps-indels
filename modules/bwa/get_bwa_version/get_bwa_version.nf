nextflow.enable.dsl = 2

params.gitc_path = "/usr/gitc"
params.bwa_path = "${params.gitc_path}/bwa"

process BWA_GET_BWA_VERSION {
    tag "BWA version"
    memory "16 GB"
    cpus 16

    output:
    stdout

    script:

    """

    ${params.bwa_path} 2>&1  \
    | grep -e '^Version'  \
    | sed 's/Version: //' \
    | tr -d '\n'
    """

    stub:
    """
    echo 0.7.15-r1140
    """
}
