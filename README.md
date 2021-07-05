# gatk4-germline-snps-indels
Workflow for germline short variant discovery using GATK4 written as per the Nextflow DSL2 best-practices.

This pipeline was developed with significant contributions by [Diamond Age Data Science](https://diamondage.com/), you can read more about the infrastructure usage on the [AWS blogpost](https://aws.amazon.com/blogs/industries/running-gatk-workflows-on-aws-a-user-friendly-solution/).

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with Docker containers making installation trivial and results highly reproducible.


# Quickstart

To run this pipeline locally on your machine, install [Nextflow](https://nextflow.io) and use the following command:

```
nextflow run 'https://github.com/seqeralabs/gatk4-germline-snps-indels' \
		 -params-file params/params.yaml \
		 -profile docker
		 -latest
```

# Usage with Nextflow Tower

Alternatively, this pipeline is readily executable with Nextflow Tower. It is availble in the Community Showcase workspace at [tower.nf](https://tower.nf). While launching the workflow, please make sure (i) to specify fill the correct parameters via the `params` field and (ii) to `Config profile` field to `docker`.

# Folder structure

The code organization is as follows.

```
$ tree -L 2
.
├── README.md
├── codespaces.md
├── conf
│   └── optimized_processes.config
├── containers
│   ├── Readme.md
│   ├── build.sh
│   └── gatk4
├── main.nf
├── modules
│   ├── bwa
│   ├── gatk
│   └── picard
├── nextflow.config
├── params
│   ├── params.yaml
│   └── test_params.yaml
├── resources
├── test_data
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


- We have used default parameters for each module, making it executable in a standalone manner. In the  `modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf`, the parameters have been initialized with default values.

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
    mainScript = 'main.nf'
    defaultBranch = 'master'
    homePage = 'https://github.com/seqeralabs/gatk4-germline-snps-indels'
    nextflowVersion = '>=21.04.1'
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

