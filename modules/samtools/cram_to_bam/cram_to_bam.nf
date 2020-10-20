//===========================
// Module parameters
//===========================



//===========================
// Process definition
//===========================

process SAMTOOLS_CRAM_TO_BAM {

//---------------------------
// directives
//---------------------------
    // TODO: We could add the samtools docker image
    container = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    // TODO: Allocate these resources dynamically. See Line-135-137 of haplotypecaller_gvcf_gatk4.wdl
    //    memory
    //    disk
    // NOTE: WDL pre-emptible is not applicable. See https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    // NOTE: Instead of WDL preemptive, we can rely on the excutor-independent `errorStrategy` and `maxRetries` directives in NXF


//---------------------------
    input:
//---------------------------
    path(ref_fasta)
    path(ref_fasta_index)
    path(ref_dict)
    path(input_cram)
    path(sample_name)

//---------------------------
    output:
//---------------------------
    path("${sample_name}.bam")
    path("${sample_name}.bai")


//---------------------------
    script:
//---------------------------
    """
    set -e
    set -o pipefall
    
    ${samtools_path} view -h -T ${ref_fasta} ${input_cram}
    ${samtools_path} view -b -o ${sample_name}.bam -
    ${samtools_path} index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai
    
    """

}
