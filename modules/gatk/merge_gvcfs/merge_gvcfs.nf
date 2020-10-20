//===========================
// Module parameters
//===========================

params.contamination
params.output_filename
params.gatk_path
params.machine_mem_gb
params.command_mem_gb

//===========================
// Process definition
//===========================

process GATK_MERGE_GVCFS {


//---------------------------
// directives
//---------------------------
    container
    memory
    disk
    // NOTE: WDL pre-emptible is not applicable. See https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    // NOTE: Instead of WDL preemptive, we can rely on the excutor-independent `errorStrategy` and `maxRetries` directives in NXF
    errorStrategy
    maxRetries

//---------------------------
    input:
//---------------------------
    path(input_vcfs)
    path(input_vcfs_indexes)


//---------------------------
    output:
//---------------------------

//---------------------------
    script:
//---------------------------

    """
    set -e

    ${gatk_path} --java-options "-Xmx${command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ${sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
    
"""
}
