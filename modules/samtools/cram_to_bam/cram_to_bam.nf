nextflow.enable.dsl = 2

import java.nio.file.Paths

params.memory = '4'
params.container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
params.tool_path = "samtools"

process SAMTOOLS_CRAM_TO_BAM {
    container = params.container
    memory "${params.samtools_cram_to_bam_memory}GB"

    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    path(input_cram)


    output:
    tuple path("${sample_name}.bam"), path("${sample_name}.bai")


    script:
    sample_name = input_cram.getBaseName()

    """
    set -e
    set -o pipefail
    
    ${params.tool_path} view -h -T ${ref_fasta} ${input_cram}
    ${params.tool_path} view -b -o ${sample_name}.bam -
    ${params.tool_path} index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai
    
    """

}

workflow test {

    params.ref_fasta_file = "${baseDir}/test_data/Homo_sapiens_assembly38.fasta"
    params.ref_fasta_fai_file = "${baseDir}/test_data/Homo_sapiens_assembly38.fasta.fai"
    params.ref_dict_file = "${baseDir}/test_data/Homo_sapiens_assembly38.dict"
    params.cram_file_pattern = "${baseDir}/test_data/*cram"

    ref_fasta_ch = Channel.value([Paths.get(params.ref_fasta_file),
                                  Paths.get(params.ref_fasta_fai_file)])

    ref_dict_ch = Channel.value(Paths.get(params.ref_dict_file))

    input_cram_ch = Channel.fromPath(params.cram_file_pattern)

    SAMTOOLS_CRAM_TO_BAM(
            ref_fasta_ch,
            ref_dict_ch,
            input_cram_ch
    )

}
