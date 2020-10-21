nextflow.preview.dsl = 2

params.memory = '10'

process SAMTOOLS_CRAM_TO_BAM {
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    memory "${params.memory}GB"
    errorStrategy 'retry'
    maxRetries 3

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
    set -o pipefall
    
    samtools view -h -T ${ref_fasta} ${input_cram}
    samtools view -b -o ${sample_name}.bam -
    samtools index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai
    
    """

}

workflow test {
    ref_fasta_ch = Channel.value(["/test_data/*fasta", "/test_data/*fasta.fai"])
    ref_dict_ch = Channel.value("/test_data/*dict")
    input_cram_ch = Channel.fromPath("/test_data/*cram")
    SAMTOOLS_CRAM_TO_BAM(ref_fasta_ch, ref_dict_ch, input_cram_ch)

}
