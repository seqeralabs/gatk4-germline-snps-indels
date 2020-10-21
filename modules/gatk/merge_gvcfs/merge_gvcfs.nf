//===========================
// Module parameters
//===========================

params.make_gvcf = false
params.memory = '10'

//===========================
// Process definition
//===========================

process GATK_MERGE_GVCFS {


//---------------------------
//  directives
//---------------------------
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    memory "${params.memory}GB"
    // NOTE: WDL pre-emptible is not applicable. See https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    // NOTE: Instead of WDL preemptive, we can rely on the excutor-independent `errorStrategy` and `maxRetries` directives in NXF
    errorStrategy
    maxRetries

//---------------------------
    input:
//---------------------------
    tuple path(input_vcf), path(input_vcf_index)


//---------------------------
    output:
//---------------------------
    tuple path("${output_filename}"), path("${output_filename}.tbi")

//---------------------------
    script:
//---------------------------
    input_str = input_vcf.reduce("") { a, b -> a + " --INPUT ${b}" }.join(' ')
    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = input_vcf.getBaseName() + output_suffix

    """
    set -e

    ${gatk_path} --java-options "-Xmx${$params.memory}G"  \
      MergeVcfs \
      --INPUT ${input_str} \
      --OUTPUT ${output_filename}
    
    """
}
