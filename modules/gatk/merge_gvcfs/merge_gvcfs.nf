nextflow.enable.dsl = 2

params.make_gvcf = false
params.memory = '10'
params.container = "broadinstitute/gatk:4.1.8.1"
params.tool_path = "/gatk/gatk"

process GATK_MERGE_GVCFS {
    container params.container
    memory "${params.memory}GB"

    input:
    path(input_vcf)


    output:
    tuple path("${output_filename}"), path("${output_filename}.tbi")

    script:
    input_str =  input_vcf.findAll { file(it).getExtension() == "gz"}
                  	.join(" --INPUT ")



    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = "merged" + output_suffix


    """
    set -e

    ${params.tool_path} --java-options "-Xmx${params.memory}G"  \
                              MergeVcfs \
                             --INPUT ${input_str} \
                             --OUTPUT ${output_filename}
    
    """
}


workflow test {

    params.vcf_file_pattern = "${baseDir}/test_data/*vcf*"
   
    input_vcf_ch = Channel.fromPath(params.vcf_file_pattern).collect()

    GATK_MERGE_GVCFS(input_vcf_ch)

}
