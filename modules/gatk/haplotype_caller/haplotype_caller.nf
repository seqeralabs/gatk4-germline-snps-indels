nextflow.enable.dsl = 2

import java.nio.file.Paths


params.contamination = 0
//params.make_gvcf = false
//params.make_bamout = false

params.memory = '4'
params.gatk_java_opts = ''

process GATK_HAPLOTYPE_CALLER {
    container = "broadinstitute/gatk:4.1.8.1"
    memory "${params.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    tuple path(input_bam), path(input_bam_index)
    path(interval_list)


    output:
    tuple path("${output_file_name}"), path("${output_file_name}.tbi")
    path("${vcf_basename}.bamout.bam") optional true

    script:
    output_suffix = params.make_gvcf ? ".g.vcf.gz" : ".vcf.gz"
    output_filename = input_bam.getBaseName() + output_suffix
    bamout_arg = params.make_bamout ? "-bamout ${input_bam.getBaseName()}.bamout.bam" : ""

    """
    set -e

    /gatk/gatk --java-options "-Xmx${params.memory}G ${params.gatk_java_opts}" \
          HaplotypeCaller \
          -R ${ref_fasta} \
          -I ${input_bam} \
          -L ${interval_list} \
          -O ${output_filename} \
          -contamination ${params.contamination} \
          -G StandardAnnotation \
          -G StandardHCAnnotation \
          ${params.make_gvcf ? "-G AS_StandardAnnotation" : ""} \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
          ${params.make_gvcf ? "-ERC GVCF" : ""} \
          ${bamout_arg}
    """
}


workflow test {

    params.make_gvcf = true
    params.make_bamout = true

    ref_fasta_ch = Channel.value([Paths.get("./test_data/Homo_sapiens_assembly38.fasta"), Paths.get("./test_data/Homo_sapiens_assembly38.fasta.fai")])
    ref_dict_ch = Channel.value(Paths.get("./test_data/Homo_sapiens_assembly38.dict"))
    input_bam_ch = Channel.fromPath(["./test_data/*bam", "./test_data/*bai"])
    interval_list_ch = Channel.value(Paths.get("./test_data/*intervals.txt"))

    GATK_HAPLOTYPE_CALLER(
            ref_fasta_ch,
            ref_dict_ch,
            input_bam_ch,
            interval_list_ch
    )

}
