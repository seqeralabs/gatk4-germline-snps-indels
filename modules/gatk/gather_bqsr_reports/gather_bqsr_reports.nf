nextflow.enable.dsl = 2

params.gatk_path = "gatk"
params.java_opts = "-Xms3000m"


process GATK_GATHER_BQSR_REPORTS {
    tag "${sampleId}"
    memory "16 GB"
    cpus 16

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

    stub:

    """
    touch "${sampleId}.recal_data.csv"
    """
}
