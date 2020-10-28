nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
// FIXME
params.java_opts = "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"


process GATK_GATHER_BQSR_REPORTS {
    tag "${sampleId}"

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus

    input:
    tuple val(sampleId), path(input_bqsr_reports)

    output:
    tuple val(sampleId),
            path("${sampleId}.recal_data.csv")

    script:
    input_bqsr_params = input_bqsr_reports.join(" -I ")

    """
    ${params.gatk_path} --java-options ${params.java_opts} \
                        GatherBQSRReports \
                        -I ${input_bqsr_params} \
                        -O "${sampleId}.recal_data.csv"
    """
}
