#!/bin/bash
set -uex

# NOTE: Make sure you've set the environment correctly and are logged in to the registry.

QUAY_NAMESPACE="seqeralabs"

for container_dir in $(find * -type d); do
  echo "Building $container_dir ..."
  cd $container_dir
  CONTAINER_TAG=0.0.1
  CONTAINER_NAME=quay.io/$QUAY_NAMESPACE/$container_dir:$CONTAINER_TAG
  docker build -t $CONTAINER_NAME .
  CONTAINER_ID=$(docker run -d $CONTAINER_NAME)
  docker commit $CONTAINER_ID $CONTAINER_NAME
  docker push quay.io/$QUAY_NAMESPACE/$container_dir
  docker stop $CONTAINER_ID
  cd ..
done

#---------------------------------------

# NOTE: An example of building and pushing multiple containers to ECR registry.
# TODO: Please modify this script for your own use-case.

#AWS_REGISTRY="ecr/MY_REGISTRY"
#AWS_REGION="MY_AWS_REGION"
#for container_dir in $(find * -type d); do
#  echo "Building $container_dir ..."
#  cd $container_dir
#  aws ecr create-repository \
#    --repository-name $container_dir \
#    --image-scanning-configuration scanOnPush=false \
#    --region $MY_AWS_REGION
#  docker build -t $AWS_REGISTRY/$container_dir:0.0.1 .
#  docker push $AWS_REGISTRY/$container_dir:0.0.1
#  cd ..
#done
