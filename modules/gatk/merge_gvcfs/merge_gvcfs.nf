nextflow.enable.dsl = 2

params.make_gvcf = false
params.gatk_merge_gvcfs_memory = '10'


process GATK_MERGE_GVCFS {
    container = "broadinstitute/gatk:4.1.8.1"
    memory "${params.gatk_merge_gvcfs_memory}GB"

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

    /gatk/gatk --java-options "-Xmx${params.gatk_merge_gvcfs_memory}G"  \
                              MergeVcfs \
                             --INPUT ${input_str} \
                             --OUTPUT ${output_filename}
    
    """
}


workflow test {
    input_vcf_ch = Channel.fromPath("${baseDir}/test_data/*vcf*").collect()

    GATK_MERGE_GVCFS(input_vcf_ch)

}
