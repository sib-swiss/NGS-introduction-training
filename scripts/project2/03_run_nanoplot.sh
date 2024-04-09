#!/usr/bin/env bash

cd ~/workspace/projects/project2

mkdir -p nanoplot

NanoPlot \
--fastq_rich reads/Cell_2.fastq.gz \
--outdir nanoplot/Cell_2

NanoPlot \
--fastq_rich reads/EV_2.fastq.gz \
--outdir nanoplot/EV_2