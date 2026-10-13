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
#!/usr/bin/env bash
set -euo pipefail

# Note that all pasted code is run from .sh files

# Add path that all subsequent code will reuse to bashrc
export READS_DIR="/nfs/home/bnsoh/test2/reads"
export RESULTS_DIR="/nfs/home/bnsoh/test2/results"
export REFERENCE_DIR="/nfs/home/bnsoh/test2/reference"

mkdir -p "$READS_DIR" "$RESULTS_DIR" "$REFERENCE_DIR"

# resuable function to link a file indemnpotently so we dont duplicate and make a mess
link_file() {
    local file="$1" target="$2/${1##*/}"
    [[ -L "$target" && "$(readlink "$target")" == "$file" ]] && return 0
    ln -s "$file" "$target"
}

# make symlinks of the data from the course project to my test2 folder
for file in /lustre/isaac24/proj/UTK0505/test2/reads/*.fastq.gz; do
    link_file "$file" "$READS_DIR"
done

# Link reference genome and indexes for bwamem
for file in /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa*; do
    link_file "$file" "$REFERENCE_DIR"
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

MultiQC gives a rounded number for total reads and does not report total bases per sample. To find out those numbers use this:

```
for f in *.fastq.gz; do
  echo -ne "$f\t"
  zcat "$f" | awk 'NR%4==2{r++; b+=length($0)} END{printf "%d\t%d\t%.1f\n", r, b, b/r}'
done | column -t
```

This gives the total reads, then the total bases, then the average read length for each file.

### Results

MultiQC report path and name: `/nfs/home/bnsoh/test2/results/01_multiqc/multiqc_report.html`

One sample has adapter and polyG issues, one sample has low quality, one sample has low coverage, and one is problem free. Label each sample below with its quality issues.
For each sample, report R1 + R2.

| Sample | Total reads | Total bases | Quality issues? |
| --- | --- | --- | --- |
| A | 1609300 | 241395000 | no obvious issues |
| B | 321860 | 48279000 | Too few unique reads comparatively |
| C | 1609300 | 241395000 | Lots of overrepresented sequences perhaps polyG tails |
| D | 1609300 | 241395000 | Generally very poor quality scores |


## Step 2: Read Trimming and Filtering (fastp)

WARNING: Because these reads are simulated, fastp fails to identify that they need poly-G tail trimming. Turn that on with a flag.

ADVICE: Its easier to compare the reads pre and post trimming if you run multiqc on the results of fastqc of the trimmed files (and do not use the .json files from fastp).

### Commands

Comments 

```
#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$RESULTS_DIR/02_fastp"
mkdir -p "$RESULTS_DIR/02_fastp_fastqc"
mkdir -p "$RESULTS_DIR/02_fastp_multiqc"

# Clean read pairs together and save separate reports per sample
for sample in A B C D; do
    fastp \
        -i "$READS_DIR/Sample_${sample}_R1.fastq.gz" \
        -I "$READS_DIR/Sample_${sample}_R2.fastq.gz" \
        -o "$RESULTS_DIR/02_fastp/Sample_${sample}_R1.trimmed.fastq.gz" \
        -O "$RESULTS_DIR/02_fastp/Sample_${sample}_R2.trimmed.fastq.gz" \
        --trim_poly_g \
        --html "$RESULTS_DIR/02_fastp/Sample_${sample}.fastp.html" \
        --json "$RESULTS_DIR/02_fastp/Sample_${sample}.fastp.json"
done

# Check cleaned reads and combine their FastQC reports
fastqc "$RESULTS_DIR"/02_fastp/*.trimmed.fastq.gz \
    -o "$RESULTS_DIR/02_fastp_fastqc"

multiqc "$RESULTS_DIR/02_fastp_fastqc" \
    -o "$RESULTS_DIR/02_fastp_multiqc"

```

MultiQC gives a rounded number for total reads and does not report total bases per sample. To find out those numbers use this:

```
for f in *.fastq.gz; do
  echo -ne "$f\t"
  zcat "$f" | awk 'NR%4==2{r++; b+=length($0)} END{printf "%d\t%d\t%.1f\n", r, b, b/r}'
done | column -t
```

This gives the total reads, then the total bases, then the average read length for each file.

### Results

MultiQC report path and name: `/nfs/home/bnsoh/test2/results/02_fastp_multiqc/multiqc_report.html`

For each sample, report R1 + R2. Use the reference genome size above to calculate coverage.

| Sample | Reads before | Reads after | Mean length before | Mean length after |
| --- | --- | --- | --- | --- |
| A | 1609300 | 1609300 | 150.0 | 150.0 |
| B | 321860 | 321860 | 150.0 | 150.0 |
| C | 1609300 | 1609300 | 150.0 | 143.8 |
| D | 1609300 | 1609300 | 150.0 | 150.0 |


## Step 3: Alignment (bwa-mem2 + samtools)

### Commands
Comments 

```
#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$RESULTS_DIR/03_alignment"
mkdir -p "$RESULTS_DIR/03_multiqc"

for sample in A B C D; do
    # first align
    bwa-mem2 mem \
        "$REFERENCE_DIR/sacCer3.fa" \
        "$RESULTS_DIR/02_fastp/Sample_${sample}_R1.trimmed.fastq.gz" \
        "$RESULTS_DIR/02_fastp/Sample_${sample}_R2.trimmed.fastq.gz" \
        > "$RESULTS_DIR/03_alignment/Sample_${sample}.sam"

    # then convert sam to sorted BAM and sort by reference position
    samtools sort \
        -o "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" \
        "$RESULTS_DIR/03_alignment/Sample_${sample}.sam"

    # index the BAM for access by genomic region
    samtools index "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam"

    # summarize alignment counts
    samtools flagstat "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" \
        > "$RESULTS_DIR/03_alignment/Sample_${sample}.flagstat.txt"
done

multiqc "$RESULTS_DIR/03_alignment" -o "$RESULTS_DIR/03_multiqc"

```


### Results

MultiQC report path and name: `/nfs/home/bnsoh/test2/results/03_multiqc/multiqc_report.html`

| Sample | % mapped | % properly paired | Mean depth |
| --- | --- | --- | --- |
| A | 100.0% | 100.0% | 20.00 |
| B | 100.0% | 100.0% | 4.00 |
| C | 98.0% | 98.0% | 16.46 |
| D | 100.0% | 100.0% | 19.13 |


## Step 4: Variant Calling (bcftools mpileup + call)

### Commands

Comments 

```
#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$RESULTS_DIR/04_variants" "$RESULTS_DIR/04_multiqc"

for sample in A B C D; do
    bcftools mpileup \
        -f "$REFERENCE_DIR/sacCer3.fa" \
        -a FORMAT/DP,FORMAT/AD \
        -Ou \
        "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" |
    bcftools call \
        -m -v \
        -Ov \
        -o "$RESULTS_DIR/04_variants/Sample_${sample}.vcf"
    bcftools stats "$RESULTS_DIR/04_variants/Sample_${sample}.vcf" \
        > "$RESULTS_DIR/04_variants/Sample_${sample}.stats.txt"
done

multiqc "$RESULTS_DIR/04_variants" -o "$RESULTS_DIR/04_multiqc" --force

```

### Results

bcftools stats report path and name: `/nfs/home/bnsoh/test2/results/04_multiqc/multiqc_report.html`

## Step 5: Filtering (bcftools filter)

### Commands

I am going to try a gentle filtering and more stingent filtering. (WARNING: These are different filters than what we have seen in class.)

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

bcftools stats report path and name for gentle filtering: `/nfs/home/bnsoh/test2/results/05_multiqc/multiqc_report.html`

bcftools stats report path and name: for stringent filtering: `/nfs/home/bnsoh/test2/results/05_multiqc/multiqc_report.html`. I placed both reports in one file.

| VCF | number of SNPs |
| --- | --- |
| A raw | 4839 |
| A gentle | 4825 |
| A stringent | 4704 |
| B raw | 3930 |
| B gentle | 1524 |
| B stringent | 54 |
| C raw | 4819 |
| C gentle | 4779 |
| C stringent | 4525 |
| D raw | 4739 |
| D gentle | 4492 |
| D stringent | 4208 |


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

From `/nfs/home/bnsoh/test2` that holds all my analysis, here is the directory structure (paste the output of: `tree .`):

```
.
├── reads
│   ├── Sample_A_R1.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_A_R1.fastq.gz
│   ├── Sample_A_R2.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_A_R2.fastq.gz
│   ├── Sample_B_R1.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_B_R1.fastq.gz
│   ├── Sample_B_R2.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_B_R2.fastq.gz
│   ├── Sample_C_R1.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_C_R1.fastq.gz
│   ├── Sample_C_R2.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_C_R2.fastq.gz
│   ├── Sample_D_R1.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_D_R1.fastq.gz
│   └── Sample_D_R2.fastq.gz -> /lustre/isaac24/proj/UTK0505/test2/reads/Sample_D_R2.fastq.gz
├── reference
│   ├── sacCer3.fa -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa
│   ├── sacCer3.fa.0123 -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa.0123
│   ├── sacCer3.fa.amb -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa.amb
│   ├── sacCer3.fa.ann -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa.ann
│   ├── sacCer3.fa.bwt.2bit.64 -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa.bwt.2bit.64
│   ├── sacCer3.fa.fai -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa.fai
│   └── sacCer3.fa.pac -> /lustre/isaac24/proj/UTK0505/test2/ref/sacCer3.fa.pac
├── results
│   ├── 01_fastqc
│   │   ├── Sample_A_R1_fastqc.html
│   │   ├── Sample_A_R1_fastqc.zip
│   │   ├── Sample_A_R2_fastqc.html
│   │   ├── Sample_A_R2_fastqc.zip
│   │   ├── Sample_B_R1_fastqc.html
│   │   ├── Sample_B_R1_fastqc.zip
│   │   ├── Sample_B_R2_fastqc.html
│   │   ├── Sample_B_R2_fastqc.zip
│   │   ├── Sample_C_R1_fastqc.html
│   │   ├── Sample_C_R1_fastqc.zip
│   │   ├── Sample_C_R2_fastqc.html
│   │   ├── Sample_C_R2_fastqc.zip
│   │   ├── Sample_D_R1_fastqc.html
│   │   ├── Sample_D_R1_fastqc.zip
│   │   ├── Sample_D_R2_fastqc.html
│   │   └── Sample_D_R2_fastqc.zip
│   ├── 01_multiqc
│   │   ├── multiqc_data
│   │   │   ├── fastqc-status-check-heatmap.txt
│   │   │   ├── fastqc_adapter_content_plot.txt
│   │   │   ├── fastqc_overrepresented_sequences_plot.txt
│   │   │   ├── fastqc_per_base_n_content_plot.txt
│   │   │   ├── fastqc_per_base_sequence_quality_plot.txt
│   │   │   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   │   │   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   │   │   ├── fastqc_per_sequence_quality_scores_plot.txt
│   │   │   ├── fastqc_sequence_counts_plot.txt
│   │   │   ├── fastqc_sequence_duplication_levels_plot.txt
│   │   │   ├── fastqc_top_overrepresented_sequences_table.txt
│   │   │   ├── llms-full.txt
│   │   │   ├── multiqc.log
│   │   │   ├── multiqc.parquet
│   │   │   ├── multiqc_citations.txt
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_fastqc.txt
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc_software_versions.txt
│   │   │   └── multiqc_sources.txt
│   │   └── multiqc_report.html
│   ├── 02_fastp
│   │   ├── Sample_A.fastp.html
│   │   ├── Sample_A.fastp.json
│   │   ├── Sample_A_R1.trimmed.fastq.gz
│   │   ├── Sample_A_R2.trimmed.fastq.gz
│   │   ├── Sample_B.fastp.html
│   │   ├── Sample_B.fastp.json
│   │   ├── Sample_B_R1.trimmed.fastq.gz
│   │   ├── Sample_B_R2.trimmed.fastq.gz
│   │   ├── Sample_C.fastp.html
│   │   ├── Sample_C.fastp.json
│   │   ├── Sample_C_R1.trimmed.fastq.gz
│   │   ├── Sample_C_R2.trimmed.fastq.gz
│   │   ├── Sample_D.fastp.html
│   │   ├── Sample_D.fastp.json
│   │   ├── Sample_D_R1.trimmed.fastq.gz
│   │   └── Sample_D_R2.trimmed.fastq.gz
│   ├── 02_fastp_fastqc
│   │   ├── Sample_A_R1.trimmed_fastqc.html
│   │   ├── Sample_A_R1.trimmed_fastqc.zip
│   │   ├── Sample_A_R2.trimmed_fastqc.html
│   │   ├── Sample_A_R2.trimmed_fastqc.zip
│   │   ├── Sample_B_R1.trimmed_fastqc.html
│   │   ├── Sample_B_R1.trimmed_fastqc.zip
│   │   ├── Sample_B_R2.trimmed_fastqc.html
│   │   ├── Sample_B_R2.trimmed_fastqc.zip
│   │   ├── Sample_C_R1.trimmed_fastqc.html
│   │   ├── Sample_C_R1.trimmed_fastqc.zip
│   │   ├── Sample_C_R2.trimmed_fastqc.html
│   │   ├── Sample_C_R2.trimmed_fastqc.zip
│   │   ├── Sample_D_R1.trimmed_fastqc.html
│   │   ├── Sample_D_R1.trimmed_fastqc.zip
│   │   ├── Sample_D_R2.trimmed_fastqc.html
│   │   └── Sample_D_R2.trimmed_fastqc.zip
│   ├── 02_fastp_multiqc
│   │   ├── multiqc_data
│   │   │   ├── fastqc-status-check-heatmap.txt
│   │   │   ├── fastqc_adapter_content_plot.txt
│   │   │   ├── fastqc_overrepresented_sequences_plot.txt
│   │   │   ├── fastqc_per_base_n_content_plot.txt
│   │   │   ├── fastqc_per_base_sequence_quality_plot.txt
│   │   │   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   │   │   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   │   │   ├── fastqc_per_sequence_quality_scores_plot.txt
│   │   │   ├── fastqc_sequence_counts_plot.txt
│   │   │   ├── fastqc_sequence_duplication_levels_plot.txt
│   │   │   ├── fastqc_sequence_length_distribution_plot.txt
│   │   │   ├── fastqc_top_overrepresented_sequences_table.txt
│   │   │   ├── llms-full.txt
│   │   │   ├── multiqc.log
│   │   │   ├── multiqc.parquet
│   │   │   ├── multiqc_citations.txt
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_fastqc.txt
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc_software_versions.txt
│   │   │   └── multiqc_sources.txt
│   │   └── multiqc_report.html
│   ├── 03_alignment
│   │   ├── Sample_A.flagstat.txt
│   │   ├── Sample_A.sam
│   │   ├── Sample_A.sorted.bam
│   │   ├── Sample_A.sorted.bam.bai
│   │   ├── Sample_B.flagstat.txt
│   │   ├── Sample_B.sam
│   │   ├── Sample_B.sorted.bam
│   │   ├── Sample_B.sorted.bam.bai
│   │   ├── Sample_C.flagstat.txt
│   │   ├── Sample_C.sam
│   │   ├── Sample_C.sorted.bam
│   │   ├── Sample_C.sorted.bam.bai
│   │   ├── Sample_D.flagstat.txt
│   │   ├── Sample_D.sam
│   │   ├── Sample_D.sorted.bam
│   │   ├── Sample_D.sorted.bam.bai
│   │   └── mean_depth_mqc.tsv
│   ├── 03_multiqc
│   │   ├── multiqc_data
│   │   │   ├── llms-full.txt
│   │   │   ├── multiqc.log
│   │   │   ├── multiqc.parquet
│   │   │   ├── multiqc_citations.txt
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc_mean_depth.txt
│   │   │   ├── multiqc_samtools_flagstat.txt
│   │   │   ├── multiqc_sources.txt
│   │   │   ├── samtools-flagstat-pct-table.txt
│   │   │   └── samtools-flagstat-table.txt
│   │   └── multiqc_report.html
│   ├── 04_multiqc
│   │   ├── multiqc_data
│   │   │   ├── bcftools-stats-subtypes.txt
│   │   │   ├── bcftools_stats_indel-lengths.txt
│   │   │   ├── bcftools_stats_variant_depths.txt
│   │   │   ├── bcftools_stats_vqc_Count_Indels.txt
│   │   │   ├── bcftools_stats_vqc_Count_SNP.txt
│   │   │   ├── bcftools_stats_vqc_Count_Transitions.txt
│   │   │   ├── bcftools_stats_vqc_Count_Transversions.txt
│   │   │   ├── llms-full.txt
│   │   │   ├── multiqc.log
│   │   │   ├── multiqc.parquet
│   │   │   ├── multiqc_bcftools_stats.txt
│   │   │   ├── multiqc_citations.txt
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc_software_versions.txt
│   │   │   └── multiqc_sources.txt
│   │   └── multiqc_report.html
│   ├── 04_variants
│   │   ├── Sample_A.stats.txt
│   │   ├── Sample_A.vcf
│   │   ├── Sample_B.stats.txt
│   │   ├── Sample_B.vcf
│   │   ├── Sample_C.stats.txt
│   │   ├── Sample_C.vcf
│   │   ├── Sample_D.stats.txt
│   │   └── Sample_D.vcf
│   ├── 05_filtering
│   │   ├── Sample_A.gentle.stats.txt
│   │   ├── Sample_A.gentle.vcf
│   │   ├── Sample_A.stringent.stats.txt
│   │   ├── Sample_A.stringent.vcf
│   │   ├── Sample_B.gentle.stats.txt
│   │   ├── Sample_B.gentle.vcf
│   │   ├── Sample_B.stringent.stats.txt
│   │   ├── Sample_B.stringent.vcf
│   │   ├── Sample_C.gentle.stats.txt
│   │   ├── Sample_C.gentle.vcf
│   │   ├── Sample_C.stringent.stats.txt
│   │   ├── Sample_C.stringent.vcf
│   │   ├── Sample_D.gentle.stats.txt
│   │   ├── Sample_D.gentle.vcf
│   │   ├── Sample_D.stringent.stats.txt
│   │   └── Sample_D.stringent.vcf
│   ├── 05_multiqc
│   │   ├── multiqc_data
│   │   │   ├── bcftools-stats-subtypes.txt
│   │   │   ├── bcftools_stats_indel-lengths.txt
│   │   │   ├── bcftools_stats_variant_depths.txt
│   │   │   ├── bcftools_stats_vqc_Count_Indels.txt
│   │   │   ├── bcftools_stats_vqc_Count_SNP.txt
│   │   │   ├── bcftools_stats_vqc_Count_Transitions.txt
│   │   │   ├── bcftools_stats_vqc_Count_Transversions.txt
│   │   │   ├── llms-full.txt
│   │   │   ├── multiqc.log
│   │   │   ├── multiqc.parquet
│   │   │   ├── multiqc_bcftools_stats.txt
│   │   │   ├── multiqc_citations.txt
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc_software_versions.txt
│   │   │   └── multiqc_sources.txt
│   │   └── multiqc_report.html
│   └── 06_evaluation
│       ├── Sample_A.gentle
│       │   ├── Sample_A_gentle_false_negatives.tsv
│       │   ├── Sample_A_gentle_false_positives.tsv
│       │   ├── Sample_A_gentle_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_A.raw
│       │   ├── Sample_A_raw_false_negatives.tsv
│       │   ├── Sample_A_raw_false_positives.tsv
│       │   ├── Sample_A_raw_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_A.stringent
│       │   ├── Sample_A_stringent_false_negatives.tsv
│       │   ├── Sample_A_stringent_false_positives.tsv
│       │   ├── Sample_A_stringent_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_B.gentle
│       │   ├── Sample_B_gentle_false_negatives.tsv
│       │   ├── Sample_B_gentle_false_positives.tsv
│       │   ├── Sample_B_gentle_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_B.raw
│       │   ├── Sample_B_raw_false_negatives.tsv
│       │   ├── Sample_B_raw_false_positives.tsv
│       │   ├── Sample_B_raw_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_B.stringent
│       │   ├── Sample_B_stringent_false_negatives.tsv
│       │   ├── Sample_B_stringent_false_positives.tsv
│       │   ├── Sample_B_stringent_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_C.gentle
│       │   ├── Sample_C_gentle_false_negatives.tsv
│       │   ├── Sample_C_gentle_false_positives.tsv
│       │   ├── Sample_C_gentle_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_C.raw
│       │   ├── Sample_C_raw_false_negatives.tsv
│       │   ├── Sample_C_raw_false_positives.tsv
│       │   ├── Sample_C_raw_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_C.stringent
│       │   ├── Sample_C_stringent_false_negatives.tsv
│       │   ├── Sample_C_stringent_false_positives.tsv
│       │   ├── Sample_C_stringent_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_D.gentle
│       │   ├── Sample_D_gentle_false_negatives.tsv
│       │   ├── Sample_D_gentle_false_positives.tsv
│       │   ├── Sample_D_gentle_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_D.raw
│       │   ├── Sample_D_raw_false_negatives.tsv
│       │   ├── Sample_D_raw_false_positives.tsv
│       │   ├── Sample_D_raw_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       ├── Sample_D.stringent
│       │   ├── Sample_D_stringent_false_negatives.tsv
│       │   ├── Sample_D_stringent_false_positives.tsv
│       │   ├── Sample_D_stringent_genotype_confusion.tsv
│       │   ├── evaluation.log
│       │   └── summary_all_runs.tsv
│       └── results.md
└── scripts
    ├── 00_copy_data.sh
    ├── 01_fastqc_multiqc.sh
    ├── 01_multiqc_table.sh
    ├── 02_fastp_multiqc.sh
    ├── 03_bwamem_multiqc.sh
    ├── 04_variant_calling.sh
    ├── 05_filtering.sh
    └── 06_evaluation.sh
```

### Results

| Sample | Filter | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | raw | 4839 | 0 | 161 | 1.0000 | 0.9678 | 0.9836 |
| A | gentle | 4825 | 0 | 175 | 1.0000 | 0.9650 | 0.9822 |
| A | strict | 4704 | 0 | 296 | 1.0000 | 0.9408 | 0.9695 |
| B | raw | 3728 | 202 | 1272 | 0.9486 | 0.7456 | 0.8349 |
| B | gentle | 1524 | 0 | 3476 | 1.0000 | 0.3048 | 0.4672 |
| B | strict | 54 | 0 | 4946 | 1.0000 | 0.0108 | 0.0214 |
| C | raw | 4819 | 0 | 181 | 1.0000 | 0.9638 | 0.9816 |
| C | gentle | 4779 | 0 | 221 | 1.0000 | 0.9558 | 0.9774 |
| C | strict | 4525 | 0 | 475 | 1.0000 | 0.9050 | 0.9501 |
| D | raw | 4718 | 21 | 282 | 0.9956 | 0.9436 | 0.9689 |
| D | gentle | 4491 | 1 | 509 | 0.9998 | 0.8982 | 0.9463 |
| D | strict | 4208 | 0 | 792 | 1.0000 | 0.8416 | 0.9140 |


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
