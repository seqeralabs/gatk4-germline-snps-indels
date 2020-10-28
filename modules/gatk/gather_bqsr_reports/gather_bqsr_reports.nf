
process gatherBqsrReports {
    tag "${sampleId}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

    input:
    tuple val(sampleId), path(input_bqsr_reports)

    output:
    tuple val(sampleId),
            path("${sampleId}.recal_data.csv")

    script:
    gatk_path = "/gatk/gatk"
    input_bqsr_params = input_bqsr_reports.join(" -I ")

    """
    ${gatk_path} --java-options ${params.java_opt_bqsrreport} \
    GatherBQSRReports \
    -I ${input_bqsr_params} \
    -O "${sampleId}.recal_data.csv"
    """
}
