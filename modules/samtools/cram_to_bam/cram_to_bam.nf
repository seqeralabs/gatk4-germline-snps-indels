//===========================
// Module parameters
//===========================
params.samtools_path

//===========================
// Process definition
//===========================

process SAMTOOLS_CRAM_TO_BAM {

//---------------------------
// directives
//---------------------------
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"


//---------------------------
    input:
//---------------------------
    path(ref_fasta)
    path(ref_fasta_index)
    path(ref_dict)
    path(input_cram)


//---------------------------
    output:
//---------------------------
    path("${sample_name}.bam")
    path("${sample_name}.bai")


//---------------------------
    script:
//---------------------------
    sample_name = input_cram.getBaseName()

    """
    set -e
    set -o pipefall
    
    ${params.samtools_path} view -h -T ${ref_fasta} ${input_cram}
    ${params.samtools_path} view -b -o ${sample_name}.bam -
    ${params.samtools_path} index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai
    
    """

}
