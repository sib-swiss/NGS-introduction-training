
The last part of this course will consist of project-based-learning. This means that you will work in groups on a single question. We will split up into groups of five people.

!!! note "If online"
    If the course takes place online, we will use break-out rooms to communicate within groups. Please stay in the break-out room during the day, also if you are working individually. 

## Roles & organisation

Project based learning is about learning by doing, but also about *peer instruction*. This means that you will be both a learner and a teacher. There will be differences in levels among participants, but because of that, some will learn efficiently from people that have just learned, and others will teach and increase their understanding.

Each project has **tasks** and **questions**. By performing the tasks, you should be able to answer the questions. At the start of the project, make sure that each of you gets a task assigned. You should consider the tasks and questions as a guidance. If interesting questions pop up during the project, you are **encouraged** to work on those. Also, you don't have to perform all the tasks and answer all the questions.

In the afternoon of day 2, you will divide the initial tasks, and start on the project. On day 3, you can work on the project in the morning and in the first part of the afternoon. We will conclude the projects with a **10-minute presentation** of each group.

## Project 1: Short-read RNA-seq

**Aim:** Compare `hisat2` (splice-aware) with `bwa mem` (splice unaware) while aligning a mouse RNA-seq dataset.

In this project you will be working with data from:

Singhania A, Graham CM, Gabryšová L, Moreira-Teixeira L, Stavropoulos E, Pitt JM, et al (2019). *Transcriptional profiling unveils type I and II interferon networks in blood and tissues across diseases*. Nat Commun. 10:1–21. [https://doi.org/10.1038/s41467-019-10601-6](https://doi.org/10.1038/s41467-019-10601-6)

Here's the [SRA page](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7821918). Reads have been downloaded, and can be found at `/data/reads/toxoplasma`.

### Tasks:

* Do a QC on the data with `fastqc`
* Check which options to use, and run `bwa-mem`
* Check which options to use, and run `hisat2`
* Find out how to measure computational speed, and measure it
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Compare the bam files of the two aligners in IGV
* Run `featurecounts` on both alignments
* Compare the count matrices in `R` or `python`

### Questions:

* Is there a difference in computational speed?
* What are the alignment rates?
* How do the aligners handle splicing?
* How are spliced alignments stored in the SAM file?
* Do you see differences in soft clipping?
* What would be the effect of the aligner if you would be measuring gene expression? (To investigate this you'll need to run e.g. [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)).

!!! hint "Example code hisat2"
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

!!! hint "Spliced alignments"
    Have a look at IGV on a particular gene, e.g. Nbr1

## Project 2: Long-read genome sequencing

**Aim**: Align long reads from RNA-seq data to a reference genome.

In this project, you will be working with data from:

Clark, M. B. et al (2020). *Long-read sequencing reveals the complex splicing profile of the psychiatric risk gene CACNA1C in human brain*. Molecular Psychiatry, 25(1), 37–47. [https://doi.org/10.1038/s41380-019-0583-1](https://doi.org/10.1038/s41380-019-0583-1).

Here you can find the [SRA page](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5286960). Reads have been downloaded, and can be found at `/data/reads/CACNA1C`.

It is Oxford Nanopore Technology sequencing data of amplicons of the gene CACNA1C. It is primarily used to discover new splice variants. In this project, we will align a few of the samples to the reference genome, and assess the quality of reads and the alignment.

### Tasks:

* Perform QC with `fastqc`
* Perform QC with `NanoPlot`
* Align with `minimap2` with default parameters
* Figure how you should set parameters `-x` and `-G`
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Compare different samples in read quality, alignment rates, depth, etc.
* **Bonus** see if you can identify different splice variants using [FLAIR](https://github.com/BrooksLabUCSC/flair) (see hints below).

### Questions:

* Have a look at the quality report. What are the average read lengths? Is that expected?
* What is the average read quality? What kind of accuracy would you expect?
* Note any differences between `fastqc` and `NanoPlot`? How is that compared to the publication?
* Check out the options `-x` and `-G` of `minimap2`. Are the defaults appropriate?
* You might consider using `-x map-ont` or `-x splice`. Do you see differences in the alignment in e.g. IGV?
* How are spliced alignments stored in the SAM file with the different settings of `-x`?
* How deep is the gene sequenced?

!!! hint "Intron sizes"
    Check out the the intron sizes of CACNA1C in e.g. IGV or UCSC genome browser. How does that relate to the parameter `-G`?

!!! hint "Accuracy from quality scores"
    Find the equation to calculate error probability from quality score on [Wikipedia](https://en.wikipedia.org/wiki/Phred_quality_score).

!!! hint "Comparing `fastqc` and `Nanoplot`"
    For comparing `fastqc` and `NanoPlot`, check out this [blog](https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/) of the author of NanoPlot, and this [thread](https://github.com/wdecoster/NanoPlot/issues/191).

!!! hint "Running FLAIR"
    FLAIR is a set of python scripts that can be used to identify and quantify (new) isoforms based on alignment files of long-read sequencing data. You can basically follow the pipeline as described [here](https://github.com/BrooksLabUCSC/flair).

    To use FLAIR first clone the git repository, and activate the (pre-configured) conda environment:

    ```sh
    cd
    git clone https://github.com/BrooksLabUCSC/flair.git
    conda activate flair_env
    ```

    After that, generate a FLAIR working directory and convert your `.bam` alignment file to `bed12` format.

    ```sh
    mkdir ~/flair_output

    python3 ~/flair/bin/bam2Bed12.py \
    -i alignment.bam \
    > alignment.bed12
    ```

    Generate (and run) a shell script to run the modules `flair.py correct`, `flair.py collapse` and `flair.py quantify`. To do this, carefully follow the manual at https://github.com/BrooksLabUCSC/flair. Structure and document your script(s), so you can easily re-run the analysis.

    Files you will need are:

    * Reference genome (chromosome 12 only): `/data/references/GRCh38.p13.chr12.fa`
    * GTF: `/data/references/Homo_sapiens.GRCh38.100.gtf`
    * Reads manifest: `/data/reads/lrrnaseq/batch_combined/reads_manifest.tsv`
