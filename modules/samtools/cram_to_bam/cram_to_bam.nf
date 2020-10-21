nextflow.enable.dsl = 2

import java.nio.file.Paths

params.memory = '4'


ref_fasta_ch = Channel.value([Paths.get("./test_data/NC000962_3.fasta"), Paths.get("/test_data/NC000962_3.fasta.fai")])

ref_dict_ch = Channel.value(Paths.get("/test_data/NC000962_3.dict"))

input_cram_ch = Channel.fromPath("./test_data/*cram")

process SAMTOOLS_CRAM_TO_BAM {
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    memory "${params.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(ref_fasta)
    path(ref_fasta_index)
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
//    ref_fasta_ch = Channel.value([Paths.get("./test_data/NC000962_3.fasta"), Paths.get("/test_data/NC000962_3.fasta.fai")])
    ref_fasta = "./test_data/NC000962_3.fasta"
    ref_fasta_fai = "./test_data/NC000962_3.fasta.fai"

//    ref_dict_ch = Channel.value(Paths.get("/test_data/NC000962_3.dict"))
    ref_dict = "./test_data/NC000962_3.dict"

    input_cram_ch = Channel.fromPath("./test_data/*cram")

    SAMTOOLS_CRAM_TO_BAM(ref_fasta, ref_fasta_fai, ref_dict, input_cram_ch)

}
