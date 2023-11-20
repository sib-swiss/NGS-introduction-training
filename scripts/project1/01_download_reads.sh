#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1
mkdir -p "$PROJDIR"
cd "$PROJDIR"

wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project1.tar.gz
tar -xvf project1.tar.gz
rm project1.tar.gz
