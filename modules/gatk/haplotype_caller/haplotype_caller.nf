nextflow.enable.dsl = 2

import java.nio.file.Paths

params.container = "broadinstitute/gatk:4.1.8.1"
params.tool_path = "/gatk/gatk"
params.memory = '4'
params.java_opts = "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"

params.contamination = 0
params.make_gvcf = false
params.make_bamout = false

process GATK_HAPLOTYPE_CALLER {
    container params.container
    memory "${params.gatk_haplotype_caller_memory}GB"


    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    tuple path(input_bam), path(input_bam_index)
    path(interval_list)


    output:
    tuple path("${output_vcf}"), path("${output_vcf}.tbi")
    path("${bam_basename}.bamout.bam") optional true

    script:
    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    bam_basename = input_bam.getBaseName()
    output_vcf = input_bam.getBaseName() + output_suffix
    bamout_arg = params.make_bamout ? "-bamout ${bam_basename}.bamout.bam" : ""

    """
    set -e

    /gatk/gatk --java-options "-Xmx${params.gatk_haplotype_caller_memory}G ${params.gatk_haplotype_caller_java_opts}" \
          HaplotypeCaller \
          -R ${ref_fasta} \
          -I ${input_bam} \
          -O ${output_vcf} \
          -contamination ${params.gatk_haplotype_caller_contamination} \
          ${params.make_gvcf ? "-G AS_StandardAnnotation" : ""} \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
          ${params.make_gvcf ? "-ERC GVCF" : ""} \
          ${bamout_arg}
    """


//    -L ${interval_list} \
}


workflow test {


    params.ref_fasta_file = "${baseDir}/test_data/Homo_sapiens_assembly38.fasta"
    params.ref_fasta_fai_file = "${baseDir}/test_data/Homo_sapiens_assembly38.fasta.fai"
    params.ref_dict_file = "${baseDir}/test_data/Homo_sapiens_assembly38.dict"
    params.bam_file_pattern = "${baseDir}/test_data/*bam*"
    params.interval_list = "${baseDir}/test_data/test-intervals.hg38.list"


    ref_fasta_ch = Channel.value([Paths.get(params.ref_fasta_file),
                                  Paths.get(params.ref_fasta_fai_file)])

    ref_dict_ch = Channel.value(Paths.get(params.ref_dict_file))

    input_bam_ch = Channel.fromPath(params.bam_file_pattern).toSortedList().flatten().collate(2)

    interval_list_ch = Channel.value(Paths.get())


    GATK_HAPLOTYPE_CALLER(
            ref_fasta_ch,
            ref_dict_ch,
            input_bam_ch,
            interval_list_ch
    )

}




