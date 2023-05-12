#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project1
mkdir -p "$WORKDIR"
cd "$WORKDIR"

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project1.tar.gz
tar -xvf project1.tar.gz
rm project1.tar.gz
