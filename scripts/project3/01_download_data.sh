#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project3
mkdir -p "$WORKDIR"
cd "$WORKDIR"

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3.tar.gz
tar -xvf project3.tar.gz
rm project3.tar.gz
