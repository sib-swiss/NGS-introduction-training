#!/usr/bin/env bash

PROJDIR=/config/project/projects/project3
mkdir -p "$PROJDIR"
cd "$PROJDIR"

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3.tar.gz
tar -xvf project3.tar.gz
rm project3.tar.gz
