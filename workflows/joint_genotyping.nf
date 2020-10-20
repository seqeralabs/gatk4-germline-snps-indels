nextflow.preview.dsl = 2

// Import modules
include { CHECK_SAMPLES_UNIQUE } from "../modules/check_samples_unique/check_samples_unique"
include { SPLIT_INTERVAL_LIST } from "../modules/split_interval_list/split_interval_list"
include { IMPORT_GVCFS } from "../modules/import_gvcfs/import_gvcfs"
include { GENOTYPE_GVCFS } from "../modules/genotype_gvcfs/genotype_gvcfs"
include { GNARLY_GENOTYPER } from "../modules/gnarly_genotyper/gnarly_genotyper"
include { HARD_FILTER_AND_MAKE_SITES_ONLY_VCF } from "../modules/hard_filter_and_make_sites_only_vcf/hard_filter_and_make_sites_only_vcf"
include { INDELS_VARIANT_RECALIBRATOR } from "../modules/indels_variant_recalibrator/indels_variant_recalibrator"
include { SNPS_VARIANT_RECALIBRATOR_CREATE_MODEL } from "../modules/snps_variant_recalibrator_create_model/snps_variant_recalibrator_create_model"
include { SNPS_VARIANT_RECALIBRATOR } from "../modules/snps_variant_recalibrator/snps_variant_recalibrator"
include { GATHER_TRANCHES } from "../modules/gather_tranches/gather_tranches"
include { APPLY_RECALIBRATION } from "../modules/apply_recalibration/apply_recalibration"
include { GATHER_VCFS } from "../modules/gather_vcfs/gather_vcfs"
include { SELECT_FINGERPRINT_SITE_VARIANTS } from "../modules/select_fingerprint_site_variants/select_fingerprint_site_variants"
include { COLLECT_VARIANT_CALLING_METRICS } from "../modules/collect_variant_calling_metrics/collect_variant_calling_metrics"
include { GATHER_VARIANT_CALLING_METRICS } from "../modules/gather_variant_calling_metrics/gather_variant_calling_metrics"
include { CROSS_CHECK_FINGERPRINT } from "../modules/cross_check_fingerprint/cross_check_fingerprint"
include { GATHER_PICARD_METRICS } from "../modules/gather_picard_metrics/gather_picard_metrics"
include { GET_FINGERPRINTING_INTERVAL_INDICES } from "../modules/get_fingerprinting_interval_indices/get_fingerprinting_interval_indices"
include { PARTITION_SAMPLE_NAME_MAP } from "../modules/partition_sample_name_map/partition_sample_name_map"

// Parameters


// Workflow
workflow joint_genotyping {
    take:

    main:

}
