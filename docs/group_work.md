
The last part of this course will consist of project-based-learning. This means that you will work in groups on a single question. We will split up into groups of five people.

!!! note "If working with Docker"
    If you are working with Docker, I assume you are working independently and therefore can not work in a group. However, you can test your skills with these real biological datasets. Realize that the datasets and calculations are (much) bigger compared to the exercises, so check if your computer is up for it. You'll probably need around 4 cores, 16G of RAM and 50G of harddisk.

!!! note "If online"
    If the course takes place online, we will use break-out rooms to communicate within groups. Please stay in the break-out room during the day, also if you are working individually.

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/07_group_work.pdf){: .md-button }

## Roles & organisation

Project based learning is about learning by doing, but also about *peer instruction*. This means that you will be both a learner and a teacher. There will be differences in levels among participants, but because of that, some will learn efficiently from people that have just learned, and others will teach and increase their understanding.

Each project has **tasks** and **questions**. By performing the tasks, you should be able to answer the questions. You should consider the tasks and questions as a guidance. If interesting questions pop up during the project, you are **encouraged** to work on those. Also, you don't have to perform all the tasks and answer all the questions.

In the afternoon of day 1, you will start on the project. On day 3, you can work on the project in the morning and in the first part of the afternoon. We will conclude the projects with a **10-minute presentation** of each group.

## Working directories

Each group has access to a shared working directory. It is mounted in the root directory (`/`). Make a soft link in your home directory:

```sh
cd ~/project
ln -s /group_work/GROUP_NAME/ ./
# replace [GROUP_NAME] with your group directory
```

Now you can find your group directory at `~/<group name>`. Use this to share files.

!!! warning
    Do not remove the soft link with `rm -r`, this will delete the entire source directory. If you want to remove only the softlink, use `rm` (without `-r`), or `unlink`. More info [here](https://linuxize.com/post/how-to-remove-symbolic-links-in-linux/).


## Project 1: Variant analysis of human data

**Aim**: Find variants on chromosome 20 from three samples

In this project you will be working with Illumina reads from three samples: a father, mother and a child. You will perform quality control, align the reads, mark duplicates, detect variants and visualize them. 

You can get the data by running these commands:

```sh
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project1.tar.gz
tar -xvf project1.tar.gz
rm project1.tar.gz
```

### Tasks

!!! warning "Important!"
    **Stick to the principles for reproducible analysis** described [here](day1/reproducibility.md)

* Download the required data
* Do a QC on the data with `fastqc`
* Trim adapters and low quality bases with `fastp`. Make sure to include the option `--detect_adapter_for_pe`. To prevent overwriting `fastp.html`, specify a report filename for each sample with the option `--html`. 
* After trimming the adapters, run `fastqc` again to see whether all adapters are gone.
* Create an index for bowtie2. At the same time create a fasta index (`.fai` file) with `samtools faidx`. 
* Check which options to use, and align with `bowtie2`. At the same time add readgroups to the aligned reads (see hints below). Make sure you end up with an indexed and sorted bam file. 
* Mark duplicates on the individual bam files with `gatk MarkDuplicates` (see hints below).
* Merge the three bam files with `samtools merge`. Index the bam file afterwards. 
* Run `freebayes` to call variants. Only call variants on the region `chr20:10018000-10220000` by specifying the `-r` option. 
* Load your alignments together with the vcf containing the variants in IGV. Check out e.g. `chr20:10,026,397-10,026,638`. 
* Run `multiqc` to get an overall quality report.

### Questions

* Have a look at the quality of the reads. Are there any adapters in there? Did adapter trimming change that? How is the base quality? Could you improve that?
* How many duplicates were in the different samples (hint: use `samtools flagstat`)? Why is it important to remove them for variant analysis?
* Why did you add read groups to the bam files? Where is this information added in the bam file? 
* Are there variants that look spurious? What could be the cause of that? What information in the vcf can you use to evaluate variant quality? 
* There are two high quality variants in `chr20:10,026,397-10,026,638`. What are the genotypes of the three samples according to freebayes? Is this according to what you see in the alignments? If the alternative alleles are present in the same individual, are they in phase or in repulsion? Note: you can also load vcf files in IGV. 

### Hints

You can add readgroups to the alignment file with `bowtie2` with the options `--rg-id` and `--rg`, e.g. (`$SAMPLE` is a variable containing a sample identifier):

```sh
bowtie2 \
-x ref.fa \
-1 r1.fastq.gz \
-2 r2.fastq.gz \
--rg-id $SAMPLE \
--rg SM:$SAMPLE \
```

To run `gatk MarkDuplicates` you will only need to specify `--INPUT` and `--OUTPUT`, e.g.:

```sh
gatk MarkDuplicates \
--INPUT sample.bam \
--OUTPUT sample.md.bam \
--METRICS_FILE sample.metrics.txt 
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

### Tasks

!!! warning "Important!"
    **Stick to the principles for reproducible analysis** described [here](day1/reproducibility.md)

* Check out the BioProject, and download two samples that interest you.
* Perform QC with `fastqc`
* Perform QC with `NanoPlot`
* Align with `minimap2` with default parameters
* Figure how you should set parameters `-x` and `-G`
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Compare different samples in read quality, alignment rates, depth, etc.

### Questions

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

**Aim:** Generate a count matrix to estimate differential gene expression. 

In this project you will be working with data from:

Singhania A, Graham CM, Gabryšová L, Moreira-Teixeira L, Stavropoulos E, Pitt JM, et al (2019). *Transcriptional profiling unveils type I and II interferon networks in blood and tissues across diseases*. Nat Commun. 10:1–21. [https://doi.org/10.1038/s41467-019-10601-6](https://doi.org/10.1038/s41467-019-10601-6)

Here's the [BioProject page](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA490485). Since the mouse genome is rather large, we have prepared reads for you that originate from chromosome 5. Use those for the project. Download them like this:

```sh
wget https://ngs-introduction-training.s3.eu-central-1.amazonaws.com/project3.tar.gz
tar -xvf project3.tar.gz
rm project3.tar.gz
```

### Tasks

!!! warning "Important!"
    **Stick to the principles for reproducible analysis** described [here](day1/reproducibility.md)

* Download the tar file, and find out what's in the data folder
* Do a QC on the fastq files with `fastqc`
* Trim adapters and low quality bases with `fastp`
* After trimming the adapters, run `fastqc` again to see whether all adapters are gone.
* Check which options to use, and align with `hisat2`
* Evaluate the alignment quality (e.g. alignment rates, mapping quality)
* Have a look at the alignments in IGV, e.g. check out `Sparcl1`. For this, you can use the built-in genome (*Mouse (mm10)*). Do you see any evidence for differential splicing?
* Run `featureCounts` on both alignments. Have a look at the option `-Q`. For further suggestions, see the hints below. 
* Compare the count matrices in `R` (find a script to get started [here](https://github.com/sib-swiss/NGS-introduction-training/blob/main/scripts/project3/07_run_DESeq2.R); Rstudio server is running on the same machine. Approach it with your credentials and username `rstudio`)

### Questions

* Check the description at the SRA sample page. What kind of sample is this?
* How does the quality of the reads look? Anything special about the overrepresented sequences? (Hint: [blast](https://blast.ncbi.nlm.nih.gov/) some overrepresented sequences, and see what they are)
* Did trimming improve the QC results? What could be the cause of the warnings/errors in the `fastqc` reports?
* What are the alignment rates?
* How are spliced alignments stored in the SAM file?
* Are there any differences between the treatments in the percentage of assigned alignments by `featureCounts`? What is the cause of this? 
* Can you find any genes that seem to be differentially expressed? 
* What is the effect of setting the option `-Q` in `featureCounts`?

### Hints

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

Here's some example code for `hisat2` and `featureCounts`. Everything in between `<>` should be replaced with specific arguments.

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
*.bam
```

