# Custom containers

This pipeline depends on  `gatk4` and `genomes-in-the-cloud` containers. A sample build setup is provided for `gatk4` which could be a good starting point for building your own custom containers.

Moreover, the `/containers/build.sh` file contains a basic `bash` script which builds and pushes all the containers in the `/containers` folder to the `quay.io` docker registry.  A sample has been provided for AWS ECR  as well.



