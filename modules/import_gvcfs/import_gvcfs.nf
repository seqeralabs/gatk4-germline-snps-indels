nextflow.enable.dsl = 2

params.workspace_dir_name
params.batch_size
params.sample_name_map

process GATK_IMPORT_GVCFS {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    path(interval)


    output:
    path("${params.workspace_dir_name}.tar")

    script:
    """
    set -euo pipefail

    rm -rf ${params.workspace_dir_name}

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.

    gatk --java-options -Xms8g \
                        GenomicsDBImport \
                        --genomicsdb-workspace-path ${params.workspace_dir_name} \
                        --batch-size ${params.batch_size} \
                        -L ${interval} \
                        --sample-name-map ${params.sample_name_map} \
                        --reader-threads 5 \
                        --merge-input-intervals \
                        --consolidate

    tar -cf ${params.workspace_dir_name}.tar ${params.workspace_dir_name}
    """


}



