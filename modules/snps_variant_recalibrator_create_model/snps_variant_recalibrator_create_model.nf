nextflow.enable.dsl = 2

params.model_report_filename
pararms.recalibration_filename
params.tranches_filename
params.use_allele_specific_annotations

process GATK_SNPS_VARIANT_RECALIBRATOR_CREATE_MODEL {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    path(sites_only_variant_filtered_vcf)
    path(sites_only_variant_filtered_vcf_index)
    path(hapmap_resource_vcf)
    path(omni_resource_vcf)
    path(one_thousand_genomes_resource_vcf)
    path(dbsnp_resource_vcf)
    path(hapmap_resource_vcf_index)
    path(omni_resource_vcf_index)
    path(one_thousand_genomes_resource_vcf_index)
    path(dbsnp_resource_vcf_index)

    output:
    path("${params.model_report_filename}")


    script:
    tranche_str = recalibration_tranche_values.collect().join(" -tranche ")
    an_str = recalibration_annotation_values.collect().join(" -an ")

    """
    set -euo pipefail

    gatk --java-options -Xms100g VariantRecalibrator \
                                 -V ${sites_only_variant_filtered_vcf} \
                                 -O ${recalibration_filename} \
                                 --tranches-file ${tranches_filename} \
                                 --trust-all-polymorphic \
                                 -tranche ${tranche_str} \
                                 -an ${an_str} \
                                 ${params.use_allele_specific_annotations ? '--use-allele-specific-annotations' : ''} \
                                 -mode SNP \
                                 --sample-every-Nth-variant ${downsampleFactor} \
                                 --output-model ${model_report_filename} \
                                 --max-gaussians ${max_gaussians} \
                                 -resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
                                 -resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
                                 -resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
                                 -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}
    """

}
