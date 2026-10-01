# Test 2 — Lab Notebook

Source: [Official course template](https://github.com/mestato/EPP622_2026/wiki/Test-2-Wiki-template-rendered)

**Instructions:** Copy each the template to your own github wiki. You are welcome to break into multiple pages if that makes it easier for you, just keep all the sections.

**Every command must be copy-pasteable and complete.** Not `ran fastp on the samples`. The actual command, with the actual arguments, in a formatted code block.

**Each command needs a comment.** What are you running and why? No need to be too long, but something like "Running FastQC to check the quality of the reads."


Your Name: Bryan Nsoh

Server: Sphinx

Repository URL: https://github.com/Bryan-Nsoh/EPP622-2026-SNPCalling

## Setup

Create directory and symbolically link provided files.

```bash
# Note that all pasted code is run from .sh files

# Add path that all subsequent code will reuse to bashrc
export READS_DIR="/nfs/home/bnsoh/test2/reads"
export RESULTS_DIR="/nfs/home/bnsoh/test2/results"

# make symlinks of the data from the course project to my test2 folder
#!/usr/bin/env bash
for file in /lustre/isaac24/proj/UTK0505/test2/reads/*.fastq.gz
do
        ln -s $file /nfs/home/bnsoh/test2/reads
        echo "copied $file to here"
done
```

## Step 1: Quality Control (FastQC + MultiQC)

### Commands

Create results directories that will be reused downstream. Run fastqc and multiqc sequrntially on the symlinked data. 

```
#!/usr/bin/env bash

mkdir -p "$RESULTS_DIR"/01_fastqc/
mkdir -p "$RESULTS_DIR"/01_multiqc/

for file in "$READS_DIR"/*.fastq.gz
do
        fastqc "$file" -o "$RESULTS_DIR"/01_fastqc/
done

multiqc "$RESULTS_DIR"/01_fastqc/ -o "$RESULTS_DIR"/01_multiqc/
```

### Results

MultiQC report path and name: `/nfs/home/bnsoh/test2/results/01_multiqc/multiqc_report.html`

| Sample | Total reads | Quality issues? | Average Genome Coverage |
| --- | --- | --- | --- |
| A | \- | \- | \- |
| B | \- | \- | \- |
| C | \- | \- | \- |
| D | \- | \- | \- |

## Step 2: Read Trimming and Filtering (fastp)

WARNING: Because these reads are simulated, fastp fails to identify that they need poly-G tail trimming. Turn that on with a flag.

ADVICE: Its easier to compare the reads pre and post trimming if you run multiqc on the results of fastqc of the trimmed files (and do not use the .json files from fastp).

### Commands

Comments 

```
commands
```
### Results

MultiQC report path and name: `_______`

| Sample | Reads before | Reads after | Mean length before | Mean length after | Average Genome Coverage After |
| --- | --- | --- | --- | --- | --- |
| A | \- | \- | \- | \- | \- |
| B | \- | \- | \- | \- | \- |
| C | \- | \- | \- | \- | \- |
| D | \- | \- | \- | \- | \- |

## Step 3: Alignment (bwa-mem2 + samtools)

### Commands
Comments 

```
commands
```


### Results
MultiQC report path and name: `_______`

| Sample | % mapped | % properly paired | Mean depth |
| --- | --- | --- | --- |
| A | \- | \- | \- |
| B | \- | \- | \- |
| C | \- | \- | \- |
| D | \- | \- | \- |

## Step 4: Variant Calling (bcftools mpileup + call)

### Commands

Comments 

```
commands
```

### Results
bcftools stats report path and name: `_______`

| Sample | Raw records | SNP records | Indel records | Multiallelic records |
| --- | --- | --- | --- | --- |
| A | \- | \- | \- | \- |
| B | \- | \- | \- | \- |
| C | \- | \- | \- | \- |
| D | \- | \- | \- | \- |

## Step 5: Filtering (bcftools filter)

### Commands

I am going to try a gentle filtering and more stingent filtering.

Gentle: exclude variants with a QUAL less than 20 or less than 5 total read depth across all samples.

Stringent: exclude variants with a QUAL less than 30 or less than 10 total read depth across all samples or a site level mapping quality (root mean square) of less than 40

```bash
#filtering gentle
#filtering stringent
```

```bash
# bcftools stats on gentle
# bcftools stats on stringent
```

### Results

bcftools stats report path and name for lenient filtering: `_______`

bcftools stats report path and name: for stringent filtering: `_______`

| Sample | SNP records raw | SNP records gentle | SNP records stringent |
| --- | --- | --- | --- |
| A | \- | \- | \- |
| B | \- | \- | \- |
| C | \- | \- | \- |
| D | \- | \- | \- |

## Step 6: Evaluation Against the Truth Set

### Commands

Link to the provided evaluate_snps_v4.py and truth_snps.vcf. Find the precision, recall and F1 for gently filtered and stringently filtered SNP sets.

Example of how to run evaluate_snps_v4.py:

```bash
python3 evaluate_snps_v4.py \
    --truth truth_snps.vcf --calls <yourvcffile> \
    --all-samples --tag strict --outdir outdir
```

```
code
```

## Directory structure

From `this folder` that holds all my analysis, here is the directory structure (paste the output of: `tree .`):

```
#output of tree command:
```

### Results

| Sample | Filter | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | unfiltered | \- | \- | \- | \- | \- | \- |
| A | gentle | \- | \- | \- | \- | \- | \- |
| A | strict | \- | \- | \- | \- | \- | \- |
| B | unfiltered | \- | \- | \- | \- | \- | \- |
| B | gentle | \- | \- | \- | \- | \- | \- |
| B | strict | \- | \- | \- | \- | \- | \- |
| C | unfiltered | \- | \- | \- | \- | \- | \- |
| C | gentle | \- | \- | \- | \- | \- | \- |
| C | strict | \- | \- | \- | \- | \- | \- |
| D | unfiltered | \- | \- | \- | \- | \- | \- |
| D | gentle | \- | \- | \- | \- | \- | \- |
| D | strict | \- | \- | \- | \- | \- | \- |

## Assignment Analysis - Questions to answer

Minimum one paragraph each.

> How did samples differ from each other in read count versus in read quality? Why does each matter for SNP calling?

> What are the relative merits and drawbacks of filtering? Interpret the above chart in terms of the original quality issues, and how those yielded different outcomes. When is filtering a good idea vs when it is it a bad idea? Think specifically about false positives vs false negatives. 

> Pick **two false negatives** and **two false positives**, from any sample. For each, tell me the sample, tell me the chromosome and coordinate location, put in a screenshot of IGV or JBrowse, and explain the specific cause of the FN or FP call.

## AI Analysis - Pipeline Audit

This is an evaluation of the pipeline you generated at the beginning of the learning module. Provide full sentences and thoughtful rationale.

Model used:

Date generated:

Exact prompt:

Minimum one paragraph each.

> Three things it got right.

> Three things you don't understand and would want to learn more about before using. Or do you understand all of it? If so be ready to answer questions about it.

> How specifically would you assess the results of this pipeline for accuracy?

> What did the labs teach you that turned out to matter for evaluating the AI-generated code?

> What is the value of learning about a computational genomics task before using AI for it — or is there none?
