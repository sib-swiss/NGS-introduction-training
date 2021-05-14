
The last part of this course will consist of project-based-learning. This means that you will work in groups on a single question. We will split up into groups of five people.

!!! note "If working with Docker"
    If you are working with Docker, I assume you are working independently and therefore can not work in a group. However, you can test your skills with these real biological datasets. Realize that the datasets and calculations are (much) bigger compared to the exercises, so check if your computer is up for it. You'll probably need around 4 cores, 16G of RAM and 50G of harddisk.

!!! note "If online"
    If the course takes place online, we will use break-out rooms to communicate within groups. Please stay in the break-out room during the day, also if you are working individually.

## Roles & organisation

Project based learning is about learning by doing, but also about *peer instruction*. This means that you will be both a learner and a teacher. There will be differences in levels among participants, but because of that, some will learn efficiently from people that have just learned, and others will teach and increase their understanding.

Each project has **tasks** and **questions**. By performing the tasks, you should be able to answer the questions. You should consider the tasks and questions as a guidance. If interesting questions pop up during the project, you are **encouraged** to work on those. Also, you don't have to perform all the tasks and answer all the questions.

In the afternoon of day 1, you will start on the project. On day 3, you can work on the project in the morning and in the first part of the afternoon. We will conclude the projects with a **10-minute presentation** of each group.

## Working directories

Each group has access to a shared working directory. It is mounted in the root directory (`/`). Make a soft link in your home directory:

```sh
cd ~
ln -s /group_work/<group name> ./
```

Now you can find your group directory at `~/<group name>`. Use this to share files.

!!! warning
    Do not remove the soft link with `rm -r`, this will delete the entire source directory. If you want to remove only the softlink, use `rm` (without `-r`), or `unlink`. More info [here](https://linuxize.com/post/how-to-remove-symbolic-links-in-linux/).

## :fontawesome-solid-seedling: Project 1: Short-read RNA-seq data of *Arabidopsis thaliana* grown in space

**Aim:** Compare `hisat2` (splice-aware) with `bowtie2` (splice unaware) while aligning an Arabidopsis RNA-seq dataset.

The analysis of this dataset is reported in this paper:

Vandenbrink JP, Herranz R, Poehlman WL, Alex Feltus F, Villacampa A, Ciska M, et al. (2019) *RNA-seq analyses of Arabidopsis thaliana seedlings after exposure to blue-light phototropic stimuli in microgravity *. Am J Bot. 106:1466–76.

Reads can be found on the NASA repository [here](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-251/). You can download data with `wget`. Here is an example for downloading the forward and reverse reads for sample 131:

```sh
wget -O sample_131_R1.fastq.gz \
https://genelab-data.ndc.nasa.gov/genelab/static/media/dataset/GLDS-251_rna-seq_13JUN2017HiSeq_Run_Sample_131_UMISS_Hoeksema_ACAGTG_L001_R1_001.fastq.gz?version=1

wget -O sample_131_R2.fastq.gz \
https://genelab-data.ndc.nasa.gov/genelab/static/media/dataset/GLDS-251_rna-seq_13JUN2017HiSeq_Run_Sample_131_UMISS_Hoeksema_ACAGTG_L001_R2_001.fastq.gz?version=1

```

!!! note "Use of `-O`"
    The name of the file is very long we use the option `-O` here to give it a shorter name.

Download the reference genome sequence like this:

```sh
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

You can download the gtf like this:

```sh
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/Araport11_GTF_genes_transposons.Mar202021.noChr.gtf.gz
gunzip Araport11_GTF_genes_transposons.Mar202021.noChr.gtf.gz
```

### Tasks:

* Check out the project page, and download one or two samples that interest you (download both the forward and reverse reads from the same sample).
* Do a QC on the data with `fastqc`
* Trim adapters and low quality bases with `cutadapt` (the adapter sequences are the same as in the exercises).
* Check which options to use, and align with `bowtie2`
* Check which options to use, and align with `hisat2`
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Compare the bam files of the two aligners in IGV. For this, download only a part of the bam file (e.g. the region `1:22145-42561`).
* Compare different samples in read quality, alignment rates, depth, etc.
* Run `featureCounts` on both alignments
* Compare the count matrices in `R` or `python` (Rstudio server is running on the same machine. Approach it with your credentials and username `rstudio`)

### Questions:

* What are the alignment rates?
* How do the aligners handle splicing?
* How are spliced alignments stored in the SAM file? (have a look at the CIGAR string)
* Do you see differences in soft clipping?
* What would be the effect of the aligner if you would be measuring gene expression? (To investigate this you'll need to run [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)).

!!! hint "Run your processes on multiple cores!"
    We are now doing computations on a full genome, with full transcriptomic data. This is quite a bit more than we have used during the exercises. Therefore, computations take longer. However, most tools support parallel processing, in which you can specify how many cores you want to use to run in parallel. Your environment contains **four** cores, so this is also the maximum number of processes you can specify. Below you can find the options used in each command to specify multi-core processing.

    | command       	| option    	|
    |---------------	|-----------	|
    | `bowtie2-build` 	| `--threads` 	|
    | `hisat2-build`  	| `--threads` 	|
    | `fastqc`       	| `--threads` 	|
    | `cutadapt`      	| `--cores`   	|
    | `bowtie2`       	| `--threads` 	|
    | `hisat2`        	| `--threads` 	|
    | `featureCounts`  	| `-T`        	|

!!! hint "Example code `hisat2` and `featureCounts`"
    Everything in between `<>` should be replaced with specific arguments.

    Here's an example for `hisat2`:

    ```sh
    hisat2-build <reference_sequence_fasta> <index_basename>

    hisat2 \
    -x <index_basename> \
    -1 <foward_reads.fastq.gz> \
    -2 <reverse_reads.fastq.gz> \
    -p <threads> \
    | samtools sort \
    | samtools view -bh \
    > <alignment_file.bam>
    ```

    Example code `featureCounts`:

    ```sh
    featureCounts \
    -p \
    -T 2 \
    -a <annotations.gtf> \
    -o <output.counts.txt> \
    <bowtie2_alignment.bam> <hisat2_alignment.bam>
    ```

!!! hint "Reading in the count data in R"
    You can read in the count data table, and compare the log2 counts of the two aligners like this:

    ```r
    cts <- read.delim('project_work/project3/counts/counts.txt', comment.char = '#')
    plot(log2(cts$..alignments.SRR7822040.chr5.bt2.bam), log2(cts$..alignments.SRR7822040.chr5.hs2.bam))
    ```

## :fontawesome-solid-brain: Project 2: Long-read genome sequencing

**Aim**: Align long reads from RNA-seq data to a reference genome.

In this project, you will be working with data from:

Clark, M. B. et al (2020). *Long-read sequencing reveals the complex splicing profile of the psychiatric risk gene CACNA1C in human brain*. Molecular Psychiatry, 25(1), 37–47. [https://doi.org/10.1038/s41380-019-0583-1](https://doi.org/10.1038/s41380-019-0583-1).

Here you can find the [BioProject page](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB34660).

It is Oxford Nanopore Technology sequencing data of amplicons of the gene CACNA1C. It is primarily used to discover new splice variants. In this project, we will align a few of the samples to the reference genome, and assess the quality of reads and the alignment.

Download the human reference genome like this:

```sh
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### Tasks:

* Check out the BioProject, and download two samples that interest you.
* Perform QC with `fastqc`
* Perform QC with `NanoPlot`
* Align with `minimap2` with default parameters
* Figure how you should set parameters `-x` and `-G`
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Compare different samples in read quality, alignment rates, depth, etc.

### Questions:

* Have a look at the quality report. What are the average read lengths? Is that expected?
* What is the average read quality? What kind of accuracy would you expect?
* Note any differences between `fastqc` and `NanoPlot`? How is that compared to the publication?
* Check out the options `-x` and `-G` of `minimap2`. Are the defaults appropriate?
* You might consider using `-x map-ont` or `-x splice`. Do you see differences in the alignment in e.g. IGV?
* How are spliced alignments stored in the SAM file with the different settings of `-x`?
* How deep is the gene sequenced?
* Do you already see evidence for splice variants in the alignments?

!!! hint "Downloading from SRA"
    ```sh
    prefetch [SRR number]
    fastq-dump --gzip [SRR number]
    ```

!!! hint "Accuracy from quality scores"
    Find the equation to calculate error probability from quality score on [Wikipedia](https://en.wikipedia.org/wiki/Phred_quality_score).

!!! hint "Comparing `fastqc` and `Nanoplot`"
    For comparing `fastqc` and `NanoPlot`, check out this [blog](https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/) of the author of NanoPlot, and this [thread](https://github.com/wdecoster/NanoPlot/issues/191).

!!! hint "Running `minimap2`"
    Here's an example command for `minimap2`:

    ```sh
    minimap2 \
    -a \
    -x [PARAMETER] \
    -G [PARAMETER] \
    [REFERENCE].fa \
    [FASTQFILE].fastq.gz \
    | samtools sort \
    | samtools view -bh > [OUTPUT].bam
    ```

!!! hint "Intron sizes"
    Check out the the intron sizes of CACNA1C in e.g. IGV or UCSC genome browser. How does that relate to the parameter `-G`?

!!! hint "More resources"
    Need e.g. a gtf file? Here's the [ensembl page](https://www.ensembl.org/Homo_sapiens/Info/Index)

## :fontawesome-solid-disease: Project 3: Short-read RNA-seq of mice.

**Aim:** Compare `hisat2` (splice-aware) with `bowtie2` (splice unaware) while aligning a mouse RNA-seq dataset.

In this project you will be working with data from:

Singhania A, Graham CM, Gabryšová L, Moreira-Teixeira L, Stavropoulos E, Pitt JM, et al (2019). *Transcriptional profiling unveils type I and II interferon networks in blood and tissues across diseases*. Nat Commun. 10:1–21. [https://doi.org/10.1038/s41467-019-10601-6](https://doi.org/10.1038/s41467-019-10601-6)

Here's the [BioProject page](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA490485). We'll be working with a single sample of this project found [here](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7822040). Since the mouse genome is rather large, we have prepared reads for you that originate from chromosome 5. Use those for the project. Download them like this:

```sh
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/SRR7822040.chr5_R1.fastq.gz
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/SRR7822040.chr5_R2.fastq.gz
```

Download the mouse reference of chromosome 5 like this:

```sh
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.chr5.fa.gz
```

And download the gtf file for chromosome 5 like this:

```sh
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3/Mus_musculus.GRCm38.102.chr5.gtf.gz
gunzip Mus_musculus.GRCm38.102.chr5.gtf.gz
```

### Tasks:

* Check out the project page, and download one or two samples that interest you (download both the forward and reverse reads from the same sample).
* Do a QC on the data with `fastqc`
* Trim adapters and low quality bases with `cutadapt` (the adapter sequences are the same as in the exercises).
* Check which options to use, and align with `bowtie2`
* Check which options to use, and align with `hisat2`
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Compare the bam files of the two aligners in IGV. For this, download only a part of the bam file (e.g. the region `5:32592000-32999545`).
* Compare different samples in read quality, alignment rates, depth, etc.
* Run `featureCounts` on both alignments
* Compare the count matrices in `R` or `python` (Rstudio server is running on the same machine. Approach it with your credentials and username `rstudio`)

### Questions:

* Check the description at the SRA sample page. What kind of sample is this?
* How does the quality of the reads look? Anything special about the overrepresented sequences? (Hint: [blast](https://blast.ncbi.nlm.nih.gov/) some overrepresented sequences, and see what they are)
* Did trimming improve the QC results? What could be the cause of the warnings/errors in the `fastqc` reports?
* What are the alignment rates?
* How do the aligners handle splicing?
* How are spliced alignments stored in the SAM file?
* Do you see differences in soft clipping?
* What would be the effect of the aligner if you would be measuring gene expression? (To investigate this you'll need to run [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)).

!!! hint "Run your processes on multiple cores!"
    We are now doing computations on a full genome, with full transcriptomic data. This is quite a bit more than we have used during the exercises. Therefore, computations take longer. However, most tools support parallel processing, in which you can specify how many cores you want to use to run in parallel. Your environment contains **four** cores, so this is also the maximum number of processes you can specify. Below you can find the options used in each command to specify multi-core processing.

    | command       	| option    	|
    |---------------	|-----------	|
    | `bowtie2-build` 	| `--threads` 	|
    | `hisat2-build`  	| `--threads` 	|
    | `fastqc`       	| `--threads` 	|
    | `cutadapt`      	| `--cores`   	|
    | `bowtie2`       	| `--threads` 	|
    | `hisat2`        	| `--threads` 	|
    | `featureCounts`  	| `-T`        	|

!!! hint "Example code `hisat2` and `featureCounts`"
    Everything in between `<>` should be replaced with specific arguments.

    Here's an example for `hisat2`:

    ```sh
    hisat2-build <reference_sequence_fasta> <index_basename>

    hisat2 \
    -x <index_basename> \
    -1 <foward_reads.fastq.gz> \
    -2 <reverse_reads.fastq.gz> \
    -p <threads> \
    | samtools sort \
    | samtools view -bh \
    > <alignment_file.bam>
    ```

    Example code `featureCounts`:

    ```sh
    featureCounts \
    -p \
    -T 2 \
    -a <annotations.gtf> \
    -o <output.counts.txt> \
    <bowtie2_alignment.bam> <hisat2_alignment.bam>
    ```

!!! hint "Spliced alignments"
    Have a look at IGV on a particular gene, e.g. Pisd

!!! hint "Reading in the count data in R"
    You can read in the count data table, and compare the log2 counts of the two aligners like this:

    ```r
    cts <- read.delim('project_work/project3/counts/counts.txt', comment.char = '#')
    plot(log2(cts$..alignments.SRR7822040.chr5.bt2.bam), log2(cts$..alignments.SRR7822040.chr5.hs2.bam))
    ```
