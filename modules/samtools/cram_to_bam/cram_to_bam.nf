nextflow.enable.dsl = 2

import java.nio.file.Paths

params.samtools_cram_to_bam_memory = '4'

process SAMTOOLS_CRAM_TO_BAM {
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
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
    
    samtools view -h -T ${ref_fasta} ${input_cram}
    samtools view -b -o ${sample_name}.bam -
    samtools index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai
    
    """

}

workflow test {

    ref_fasta_ch = Channel.value([Paths.get("${baseDir}/test_data/Homo_sapiens_assembly38.fasta"),
                                  Paths.get("${baseDir}/test_data/Homo_sapiens_assembly38.fasta.fai")])

    ref_dict_ch = Channel.value(Paths.get("${baseDir}/test_data/Homo_sapiens_assembly38.dict"))

    input_cram_ch = Channel.fromPath("${baseDir}/test_data/*cram")

    SAMTOOLS_CRAM_TO_BAM(
            ref_fasta_ch,
            ref_dict_ch,
            input_cram_ch
    )

}
