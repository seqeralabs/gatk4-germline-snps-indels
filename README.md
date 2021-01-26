# gatk4-germline-snps-indels
Workflow for germline short variant discovery using GATK4 written as per the Nextflow DSL2 best-practices.

This pipeline was developed with significant contributions by [Diamond Age Data Science](https://diamondage.com/). 


# Usage with Tower

This pipeline is readily executable with `Tower`. 

For example, if you wish to run this pipeline with `AWS Batch` please follow the steps outlined [here for setting up the AWS enviroment](https://help.tower.nf/docs/compute-envs/aws-batch/) and [here for launching the pipeline](https://help.tower.nf/docs/launch/overview/)  on the official Tower setup docs.


# Stub-run (mock-run) for the pipeline

For testing the workflow, you can enable the `stub-run` feature on tower while launching the pipeline.

For testing locally, you can use the following command 

```
nextflow run 'https://github.com/seqeralabs/gatk4-germline-snps-indels' \
		 -params-file params.yaml \
		 -main-script gatk4-germline-snps-indels.nf \
		 -latest \
		 -stub-run
```

# Folder structure

The code organization is as follows.

```
$ tree -L 2
.
├── LICENSE
├── README.md
├── codespaces.md
├── containers
│   ├── build.sh
│   └── gatk4
├── gatk4-germline-snps-indels.nf
├── modules
│   ├── bwa
│   ├── gatk
│   ├── picard
│   └── utils
├── nextflow.config
├── params.yaml
├── test_data
│   ├── manifest-3samples.txt
├── test_params.yaml
└── workflows
    ├── format_conversion
    ├── preprocessing_mapping
    ├── quality_recalibration
    └── variant_discovery



```

## Salient features of the design

- The `workflows` folder containers the (sub-)workflows in its own unique folder

``` 
workflows/
└── format_conversion
    ├── format_conversion.nf
    ├── nextflow.config
    └── test_params.yaml
```

- The top-level file (in this case `gatk4-germline-snps-indels.nf`) is the primary workflow. It could be thought of as a counter-part for the top-level `main.nf`


- The `modules` are designed around the idea of a specific tool (for example `gatk`) as the core abstraction, which contains wrappers for various sub-commands for `gatk`.

``` 
modules/gatk/
│ 
├── paired_fastq_to_unmapped_bam
│   ├── README.md
│   ├── nextflow.config
│   ├── paired_fastq_to_unmapped_bam.nf
│   ├── test_data
│   │   ├── SRR047732_1.10000_reads.filt.fastq.gz
│   │   └── SRR047732_2.10000_reads.filt.fastq.gz
│   └── test_params.yaml
```


- We have used default parameters for each module, making it executable in a standalone manner. In the  `modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf`, the default parameters have been initialized with default values.

```
params.gatk_path = "/gatk/gatk"
params.java_opts = ""

```

- We have named the `process` within these modules following the `TOOLNAME_PROCESSNAME`pattern. For example, in the  `modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf`, the name of the wrapped process is `GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM`. This is helpful when we are stitching these processes together in a workflow file.

- At the workflow level, we have used `namespaced` parameters to override the defaults for a module as shown below

```
// NOTE: The param override must be above the inclusion of the module.

params.GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM = [
        java_opts: "-Xms3000m"

include { GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM } from "../../modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf" addParams (params.GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM)

```

This allows us to `namespace` the parameters for every module either from (i) the workflow file (ii) `nextflow.config` (iii) or the `params.yaml` file.

Here's how a `yaml` file would look like, with this pattern

``` yaml
GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM:
    java_opts:
        "-Xms3000m"
```

- Use the `modules/utils` module to keep the ad-hoc utilities and name them with the same pattern, i.e. `UTILS_PROCESSNAME`

- You can use the `containers` folder to build the individual containers.

# Parameters

It is our recommendation that all pipeline related parameters like input/output location, tool oriented options should be in the `yaml` file and all executor/container level optimization and overrides should be in the `config` file.

An example `params.yaml` is shown below 


``` yaml
fasta:
  "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
dbSNP_vcf:
  "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
```


An example `nextflow.config` file is shown below

``` nextflow

manifest {
    name = 'GATK4 Germline SNP and Indel Analysis'
    description = 'Workflow for germline short variant discovery using GATK4'
    version = '0.0.1'
    mainScript = 'gatk4-germline-snps-indels.nf'
    defaultBranch = 'master'
    homePage = 'https://github.com/seqeralabs/gatk4-germline-snps-indels'
    nextflowVersion = '>=20.07.1'
}


process.errorStrategy = 'retry'
process.maxRetries = 3


profiles {
    dockerhub {
        process {
            errorStrategy = 'retry'
            ext.registry = 'https://hub.docker.com/r/seqeralabs/'
        }

        docker {
            enabled = true
            fixOwnership = true
        }

    }
}
```


# Testing and Optimization

- `Module` level tests

These are conceptually similar to the traditional `Unit test` and the `modules` should be designed to be testable individually. For a reference please take a look at [`/modules/gatk/paired_fastq_to_unmapped_bam/README.md`](/modules/gatk/paired_fastq_to_unmapped_bam/README.md). 

- `Workflow` level tests

These are conceptually similar to the `Integration and End-to-End tests` and the same strategy for the `module-level` tests can be used for the `(sub)workflow-level` tests as well. 

