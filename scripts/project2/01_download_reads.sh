#!/usr/bin/env bash

PROJDIR=/config/workspace/projects/
mkdir -p "$PROJDIR"
cd "$PROJDIR"

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project2.tar.gz
tar -xvf project2.tar.gz
rm project2.tar.gz
