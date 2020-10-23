nextflow.enable.dsl = 2

process GATK_GATHER_VCFS {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    output:
    path("${params.output_vcf_name}")
    path("${params.output_vcf_name}.tbi")

    script:

    input_vcfs_str = input_vcfs.collect {" --input "}

    """
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g GatherVcfsCloud \
                               --ignore-safety-checks \
                               --gather-type BLOCK \
                               --input ${input_vcfs_str} \
                               --output ${params.output_vcf_name}

    tabix ${params.output_vcf_name}
   """

}
