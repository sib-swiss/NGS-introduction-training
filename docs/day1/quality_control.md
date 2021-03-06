## Learning outcomes

**After having completed this chapter you will be able to:**

- Find information about a sequence run on the Sequence Read Archive
- Run `fastqc` on sequence reads and interpret the results
- Trim adapters and low quality bases using `cutadapt`

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/quality_control.pdf){: .md-button }

* `fastqc` command line [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)
* `cutadapt` [manual](https://cutadapt.readthedocs.io/en/stable/)
* Unix command line [E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

## Exercises

### Download and evaluate an E. coli dataset

**Exercise:** If you haven't already done so, create a directory called `workdir` in your home directory and make the directory your current directory.

!!! note "If working with Docker"
    If you have mounted your local directory to `/root/workdir`, this directory should already exist.

??? done "Answer"
    ```sh
    cd ~
    mkdir workdir
    cd workdir
    ```

Check out the dataset at [SRA](https://www.ncbi.nlm.nih.gov/sra/?term=SRR519926).

**Exercise:** Browse around the SRA entry and answer these questions:

**A.** Is the dataset paired-end or single end?

**B.** Which instrument was used for sequencing?

**C.** What is the read length?

**D.** How many reads do we have?

??? done "Answers"
    A. paired-end

    B. Illumina MiSeq

    C. 2 x 251 bp

    D. 400596

Make a directory `reads` in `~/workdir` and download the reads from the SRA database using `prefetch` and `fastq-dump` from [SRA-Tools](https://ncbi.github.io/sra-tools/) into the `reads` directory:

!!! note "Running `sra-tools` for the first time"
    If you run `sra-tools` for the first time, you have to set a config file. We'll be just using a minimal file for now. To do that run:

    ```sh
    vdb-config --interactive
    ```
    When the GUI pops up, press the key ++x++. The config file will be created in `~/.ncbi/user-settings.mkfg`.

```sh
mkdir reads
cd reads
prefetch SRR519926
fastq-dump --split-files SRR519926
```

**Exercise:** Check whether the download was successful by counting the number of reads in the fastq files and compare it to the SRA entry.

!!! tip "Tip"
    A read in a fastq file consists of four lines (more on that at [file types](../day2/file_types.md)). Use Google to figure out how to count the number of reads in a fastq file.

??? done "Answer"
    e.g. from [this](https://www.biostars.org/p/139006/) thread on Biostars:

    ```sh
    ## forward read
    echo $(cat SRR519926_1.fastq | wc -l)/4 | bc

    ## reverse read
    echo $(cat SRR519926_2.fastq | wc -l)/4 | bc
    ```

### Run fastqc

**Exercise:** Run fastqc on the fastq files.

!!! tip "Tip"
    `fastqc` accepts multiple files as input, so you can use a [wildcard](https://en.wikipedia.org/wiki/Glob_(programming)) to run `fastqc` on all the files in one line of code. Use it like this: `*.fastq`.  

??? done "Answer"
    ```sh
    fastqc *.fastq
    ```

**Exercise:** Download the html files to your local computer, and view the results. How is the quality? Where are the problems?

??? done "Answer"
    There seems to be:

    * Low quality towards the 3' end (per base sequence quality)
    * Full sequence reads with low quality (per sequence quality scores)
    * Adapters in the sequences (adapter content)

    We can probably fix most of these issues by trimming.

### Trim the reads

We will use [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) for trimming adapters and low quality bases from our reads. The most used adapters for Illumina are TruSeq adapters. To run `cutadapt` you need to specify the adapter sequences with options `-a` (or `--adapter`) and `-A`. A reference for the adapter sequences can be found [here](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html).

**Exercise:** The script below will trim the sequence reads. However, some parts are missing. We want to:

* trim the reads based on a base quality of 10 both on the 3' and 5' end of the reads,
* keep only reads with a read length not shorter than 25 base pairs.

Fill in the missing options and execute the script to trim the data.

!!! hint
    Check out the helper of `cutadapt` with:
    ```sh
    cutadapt --help
    ```

```sh
##!/usr/bin/env bash

TRIMMED_DIR=trimmed_data
READS_DIR=reads

mkdir $TRIMMED_DIR

cutadapt \
--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
[QUALITY CUTOFF OPTION] \
[MINIMUM LENGTH OPTION] \
--output $TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
--paired-output $TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
$READS_DIR/SRR519926_1.fastq \
$READS_DIR/SRR519926_2.fastq
```

??? done "Answer"
    Your script should look like this:

    ```sh
    ##!/usr/bin/env bash

    TRIMMED_DIR=trimmed_data
    READS_DIR=reads

    mkdir $TRIMMED_DIR

    cutadapt \
    --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --quality-cutoff 10,10 \
    --minimum-length 25 \
    --output $TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
    --paired-output  $TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
    $READS_DIR/SRR519926_1.fastq \
    $READS_DIR/SRR519926_2.fastq
    ```

!!! note "The use of `\`"
    In the script above you see that we're using `\` at the end of many lines. We use it to tell bash to ignore the newlines. If we would not do it, the `cutadapt` command would become a very long line, and the script would become very difficult to read. It is in general good practice to put every option of a long command on a newline in your script and use `\` to ignore the newlines when executing.

**Exercise:** Run `fastqc` on the trimmed fastq files and answer these questions:

**A.** Has the quality improved?

**B.** How many reads do we have left?


??? done "Answers"
    Running `fastqc`:
    ```
    cd ~/workdir/trimmed_data
    fastqc paired_trimmed*.fastq
    ```

    A. Yes, low quality 3' end, per sequence quality and adapter sequences have improved.

    B. 315904
