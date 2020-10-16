## Creating a codespaces environment for Nextflow DSL2 development

### Install Nextflow

`wget -qO- https://get.nextflow.io | bash`
`sudo mv nextflow /usr/local/bin/`

### Install Nextflow syntax highlighting extension

Select **Extensions** from the left menu bar.
Search for **Nextflow** and select the extension described *Nextflow language support* from Nextflow.
Install the Nextflow extension.

### Install Docker extension

Select **Extensions** from the left menu bar.
Search for **Docker** and select the from Microsoft.
Install the Docker extension.

### Connect to a container registry of your choice.

Select the Docker icon from the left menu bar.
Connect to your registry (Azure, Dockerhub, AWS ECR, Google GCR etc).
This can be a public or private registry. You may want to use different registies depending on the deployment platform and the location of your compute environment. 

For this exercise, we will use Dockerhub for development and testing. Later for production scaling can change the registry to an ECR in the same region as our compute to run with AWS Batch.


### Create the nextflow configuration file.

In the terminal enter:

    $ code nextflow.config
 
Enter the following details into the Nextflow configuration:
```
manifest {
    name = 'GATK4 Germline SNP and Indel Analysis'
    description = 'Workflow for germline short variant discovery using GATK4'
    version = '0.0.1'
    author = 'Evan Floden <evan@seqera.io>'
    mainScript = 'gatk4-germline-snps-indels.nf'
    defaultBranch = 'master'
    homePage = 'https://github.com/seqeralabs/gatk4-germline-snps-indels'
    nextflowVersion = '20.07.0'
}

dockerhub {
  docker.enabled = true
  docker.fixOwnership = true
  process.ext.registry = 'https://hub.docker.com/r/seqeralabs/' }
}
```
Where process.ext.registry is the base name of the registry.


### Create first container

In the terminal enter:

    $ mkdir containers/gatk
    $ code containers/gatk/conda.env

```

```





