


The last part of this course will be project-based-learning. This means that you will work in groups on a single question. We will split up into groups of four people.

## The projects

There are two different projects you can choose from:

1. Evaluate two aligners to align short-read RNA-seq data
2. Evaluate two aligners to align long-read genomic data

## Roles & organisation

For both projects, there are three important components:

1. What makes an aligner suitable for the job, i.e. on what features should it be evaluated?
2. Run and visualise the alignment
3. The evaluation itself

We will be working in groups of four people. At the start of the project, 2 people will focus on component 1 (evaluation features), the other two on component 2 (run the alignment). Divide the roles amongst yourselves. Components 1 and 2 should eventually lead to the evaluation itself, but the evaluation will probably be performed iteratively (e.g. going back and forth between what to evaluate and run code). Make a simple workplan (roles, moments of means of communication), and send it to the instructors before the end of the day.

## Presentations

We will discuss the results of each group based on presentations. Each group will have 10 minutes for its presentation (~ 10 ppt slides), followed by 5 minutes discussion.

## Material

### Project 1: Short-read RNA-seq

hisat2 (splice-aware) vs bwa mem (splice unaware?)

Reads are from [this paper](https://www.nature.com/articles/s41467-019-10601-6), and here's the [SRA page](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7821918). Reads have been downloaded, and can be found at `/data/reads/toxoplasma`.

Compare the two aligners. In particular, you could try answering the following questions:

* Is there a difference in computational speed?
* What are the alignment rates?
* How do the aligners handle splicing?
* How are spliced alignments stored in the SAM file?
* Do you see differences in soft clipping?
* What would be the effect of the aligner if you would be measuring gene expression? (To investigate this you'll need to run e.g. [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)).

??? hint "Example code hisat2"
    Everything in between `<>` should be replaced with specific arguments

    ```sh
    hisat2-build <reference_sequence_fasta> <index_basename>

    hisat2 \
    -x <index_basename> \
    -1 <mate1.fastq.gz> \
    -2 <mate2.fastq.gz> \
    -p <threads> \
    | samtools sort \
    | samtools view -bh \
    > <alignment_file.bam>

    ```

??? hint
    Have a look at IGV on a particular gene, e.g. Nbr1

### Project 2: Long-read genome sequencing

https://github.com/Magdoll/cDNA_Cupcake/wiki/Best-practice-for-aligning-Iso-Seq-to-reference-genome:-minimap2,-deSALT,-GMAP,-STAR,-BLAT#refstar
https://www.ncbi.nlm.nih.gov/sra/?term=SRR5286960

compare hisat2 and minimap2

* Have a look at the quality report. What are the average read lengths?
* What is the average read quality? What kind of accuracy would you expect?
* Is there a difference in computational speed?
* What are the differences in alignment rates?
