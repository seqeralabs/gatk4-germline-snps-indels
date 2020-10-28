//------------------
params.GATK_GENOMICS_DB_IMPORT = [
        java_opts: "-Xmx4G"
]
include { GATK_GENOMICS_DB_IMPORT } from "../modules/gatk/genomics_db_import/genomics_db_import.nf" addParams(*: params.GATK_GENOMICS_DB_IMPORT)


nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
params.java_opts = ""


process GATK_GENOMICS_DB_IMPORT {

}
