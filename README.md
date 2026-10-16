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
        --detect_adapter_for_pe \
        --html "$RESULTS_DIR/02_fastp/Sample_${sample}.fastp.html" \
        --json "$RESULTS_DIR/02_fastp/Sample_${sample}.fastp.json"
done

# Check cleaned reads and combine their FastQC reports
fastqc "$RESULTS_DIR"/02_fastp/*.trimmed.fastq.gz \
    -o "$RESULTS_DIR/02_fastp_fastqc"

multiqc "$RESULTS_DIR/02_fastp_fastqc" \
    -o "$RESULTS_DIR/02_fastp_multiqc" --force

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
| C | 1609300 | 1577258 | 150.0 | 125.9 |
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

    # Then convert sam to sorted BAM and sort by reference position
    samtools sort \
        -o "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" \
        "$RESULTS_DIR/03_alignment/Sample_${sample}.sam"

    # Index the BAM for access by genomic region
    samtools index "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam"

    # Summarize alignment counts
    samtools flagstat "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" \
        > "$RESULTS_DIR/03_alignment/Sample_${sample}.flagstat.txt"
done

# Calculate mean depth across all reference positions, including uncovered bases
printf '# plot_type: table\n# section_name: Mean depth\nSample\tMean depth\n' \
    > "$RESULTS_DIR/03_alignment/mean_depth_mqc.tsv"
for sample in A B C D; do
    samtools depth -aa "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" |
        awk -v sample="$sample" '{ total += $3; positions++ }
            END { if (!positions) exit 1; printf "Sample_%s\t%.2f\n", sample, total / positions }' \
        >> "$RESULTS_DIR/03_alignment/mean_depth_mqc.tsv"
done

multiqc "$RESULTS_DIR/03_alignment" -o "$RESULTS_DIR/03_multiqc" --force

```


### Results

MultiQC report path and name: `/nfs/home/bnsoh/test2/results/03_multiqc/multiqc_report.html`

| Sample | % mapped | % properly paired | Mean depth |
| --- | --- | --- | --- |
| A | 100.00% | 100.00% | 20.00 |
| B | 100.00% | 100.00% | 4.00 |
| C | 100.00% | 100.00% | 16.44 |
| D | 100.00% | 100.00% | 19.13 |


## Step 4: Variant Calling (bcftools mpileup + call)

### Commands

Comments 

```
#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$RESULTS_DIR/04_variants" "$RESULTS_DIR/04_multiqc"

# Call all four samples together
bcftools mpileup \
    -f "$REFERENCE_DIR/sacCer3.fa" \
    -a FORMAT/DP,FORMAT/AD \
    -Ou \
    "$RESULTS_DIR"/03_alignment/Sample_{A,B,C,D}.sorted.bam |
bcftools call -m -v -Ov \
    -o "$RESULTS_DIR/04_variants/joint.unnamed.vcf"

# Replace BAM-path sample names with Sample_A through Sample_D
for sample in A B C D; do
    printf '%s Sample_%s\n' "$RESULTS_DIR/03_alignment/Sample_${sample}.sorted.bam" "$sample"
done > "$RESULTS_DIR/04_variants/sample_names.txt"
bcftools reheader -s "$RESULTS_DIR/04_variants/sample_names.txt" \
    "$RESULTS_DIR/04_variants/joint.unnamed.vcf" \
    -o "$RESULTS_DIR/04_variants/joint.named.vcf"
bcftools view -v snps "$RESULTS_DIR/04_variants/joint.named.vcf" \
    -Ov -o "$RESULTS_DIR/04_variants/joint.vcf"
rm "$RESULTS_DIR/04_variants/joint.unnamed.vcf" "$RESULTS_DIR/04_variants/joint.named.vcf"

bcftools stats -s - "$RESULTS_DIR/04_variants/joint.vcf" \
    > "$RESULTS_DIR/04_variants/joint.stats.txt"
multiqc "$RESULTS_DIR/04_variants" -o "$RESULTS_DIR/04_multiqc" --force

```

### Results

bcftools stats report path and name: `/nfs/home/bnsoh/test2/results/04_multiqc/multiqc_report.html`

| Stats file | Path |
| --- | --- |
| raw | `/nfs/home/bnsoh/test2/results/04_variants/joint.stats.txt` |

## Step 5: Filtering (bcftools filter)

### Commands

I am going to try a gentle filtering and more stingent filtering. (WARNING: These are different filters than what we have seen in class.)

Gentle: exclude variants with a QUAL less than 20 or less than 5 total read depth across all samples.

Stringent: exclude variants with a QUAL less than 30 or less than 10 total read depth across all samples or a site level mapping quality (root mean square) of less than 40

```bash
#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$RESULTS_DIR/05_filtering" "$RESULTS_DIR/05_multiqc"

# Filter shared sites using total depth across all samples
bcftools filter -i 'QUAL>=20 && INFO/DP>=5' \
    "$RESULTS_DIR/04_variants/joint.vcf" \
    -Ov -o "$RESULTS_DIR/05_filtering/joint.gentle.vcf"

# INFO/MQ is average mapping quality, not RMS
bcftools filter -i 'QUAL>=30 && INFO/DP>=10 && INFO/MQ>=40' \
    "$RESULTS_DIR/04_variants/joint.vcf" \
    -Ov -o "$RESULTS_DIR/05_filtering/joint.stringent.vcf"

for filter in gentle stringent; do
    bcftools stats -s - "$RESULTS_DIR/05_filtering/joint.${filter}.vcf" \
        > "$RESULTS_DIR/05_filtering/joint.${filter}.stats.txt"
done

multiqc "$RESULTS_DIR/05_filtering" -o "$RESULTS_DIR/05_multiqc" --force

```
 
### Results

bcftools stats report path and name for gentle filtering: `/nfs/home/bnsoh/test2/results/05_multiqc/multiqc_report.html`

bcftools stats report path and name: for stringent filtering: `/nfs/home/bnsoh/test2/results/05_multiqc/multiqc_report.html`. I placed both reports in one file.

| Stats file | Path |
| --- | --- |
| gentle | `/nfs/home/bnsoh/test2/results/05_filtering/joint.gentle.stats.txt` |
| stringent | `/nfs/home/bnsoh/test2/results/05_filtering/joint.stringent.stats.txt` |

| VCF | number of SNPs |
| --- | --- |
| raw | 4898 |
| gentle | 4862 |
| stringent | 4749 |


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
#!/usr/bin/env bash
set -euo pipefail

CHECK_DIR="/lustre/isaac24/proj/UTK0505/test2/check"
RESULTS_DIR="/nfs/home/bnsoh/test2/results"

mkdir -p "$RESULTS_DIR/06_evaluation"
table="$RESULTS_DIR/06_evaluation/results.md"
printf '| Sample | Filter | TP | FP | FN | Precision | Recall | F1 |\n| --- | --- | --- | --- | --- | --- | --- | --- |\n' > "$table"

for sample in A B C D; do
    for filter in raw gentle stringent; do
        if [[ "$filter" == raw ]]; then
            calls="$RESULTS_DIR/04_variants/joint.vcf"
        else
            calls="$RESULTS_DIR/05_filtering/joint.${filter}.vcf"
        fi

        mkdir -p "$RESULTS_DIR/06_evaluation/Sample_${sample}.${filter}"
        rm -f "$RESULTS_DIR/06_evaluation/Sample_${sample}.${filter}/summary_all_runs.tsv"
        python3 "$CHECK_DIR/evaluate_snps_v4.py" \
            --truth "$CHECK_DIR/truth_snps.vcf" \
            --calls "$calls" \
            --sample "Sample_${sample}" \
            --tag "$filter" \
            --depth-vcf "$RESULTS_DIR/04_variants/joint.vcf" \
            --outdir "$RESULTS_DIR/06_evaluation/Sample_${sample}.${filter}"
        # dump the results into a convenient .md table
        awk -F '\t' -v sample="$sample" -v filter="$filter" '
            END {
                if (filter == "stringent") filter = "strict"
                printf "| %s | %s | %s | %s | %s | %s | %s | %s |\n", sample, filter, $7, $8, $9, $10, $11, $12
            }' "$RESULTS_DIR/06_evaluation/Sample_${sample}.${filter}/summary_all_runs.tsv" >> "$table"
    done
done

cat "$table"

```

## Directory structure

From `/nfs/home/bnsoh/test2` that holds all my analysis, here is the directory structure (paste the output of: `tree .`):

```
$ tree -I archive
.
├── adapter_rerun.log
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
│   │   │   ├── fastqc_per_base_n_content_plot.txt
│   │   │   ├── fastqc_per_base_sequence_quality_plot.txt
│   │   │   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   │   │   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   │   │   ├── fastqc_per_sequence_quality_scores_plot.txt
│   │   │   ├── fastqc_sequence_counts_plot.txt
│   │   │   ├── fastqc_sequence_duplication_levels_plot.txt
│   │   │   ├── fastqc_sequence_length_distribution_plot.txt
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
│   │   │   ├── bcftools-stats-sequencing-depth.txt
│   │   │   ├── bcftools-stats-singletons.txt
│   │   │   ├── bcftools-stats-sites.txt
│   │   │   ├── bcftools-stats-subtypes.txt
│   │   │   ├── bcftools-stats-tstv.txt
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
│   │   ├── joint.stats.txt
│   │   ├── joint.vcf
│   │   └── sample_names.txt
│   ├── 05_filtering
│   │   ├── joint.gentle.stats.txt
│   │   ├── joint.gentle.vcf
│   │   ├── joint.stringent.stats.txt
│   │   └── joint.stringent.vcf
│   ├── 05_multiqc
│   │   ├── multiqc_data
│   │   │   ├── bcftools-stats-sequencing-depth_joint_gentle.txt
│   │   │   ├── bcftools-stats-sequencing-depth_joint_stringent.txt
│   │   │   ├── bcftools-stats-singletons_joint_gentle.txt
│   │   │   ├── bcftools-stats-singletons_joint_stringent.txt
│   │   │   ├── bcftools-stats-sites_joint_gentle.txt
│   │   │   ├── bcftools-stats-sites_joint_stringent.txt
│   │   │   ├── bcftools-stats-subtypes.txt
│   │   │   ├── bcftools-stats-tstv_joint_gentle.txt
│   │   │   ├── bcftools-stats-tstv_joint_stringent.txt
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
│       │   └── summary_all_runs.tsv
│       ├── Sample_A.raw
│       │   ├── Sample_A_raw_false_negatives.tsv
│       │   ├── Sample_A_raw_false_positives.tsv
│       │   ├── Sample_A_raw_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_A.stringent
│       │   ├── Sample_A_stringent_false_negatives.tsv
│       │   ├── Sample_A_stringent_false_positives.tsv
│       │   ├── Sample_A_stringent_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_B.gentle
│       │   ├── Sample_B_gentle_false_negatives.tsv
│       │   ├── Sample_B_gentle_false_positives.tsv
│       │   ├── Sample_B_gentle_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_B.raw
│       │   ├── Sample_B_raw_false_negatives.tsv
│       │   ├── Sample_B_raw_false_positives.tsv
│       │   ├── Sample_B_raw_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_B.stringent
│       │   ├── Sample_B_stringent_false_negatives.tsv
│       │   ├── Sample_B_stringent_false_positives.tsv
│       │   ├── Sample_B_stringent_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_C.gentle
│       │   ├── Sample_C_gentle_false_negatives.tsv
│       │   ├── Sample_C_gentle_false_positives.tsv
│       │   ├── Sample_C_gentle_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_C.raw
│       │   ├── Sample_C_raw_false_negatives.tsv
│       │   ├── Sample_C_raw_false_positives.tsv
│       │   ├── Sample_C_raw_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_C.stringent
│       │   ├── Sample_C_stringent_false_negatives.tsv
│       │   ├── Sample_C_stringent_false_positives.tsv
│       │   ├── Sample_C_stringent_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_D.gentle
│       │   ├── Sample_D_gentle_false_negatives.tsv
│       │   ├── Sample_D_gentle_false_positives.tsv
│       │   ├── Sample_D_gentle_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_D.raw
│       │   ├── Sample_D_raw_false_negatives.tsv
│       │   ├── Sample_D_raw_false_positives.tsv
│       │   ├── Sample_D_raw_genotype_confusion.tsv
│       │   └── summary_all_runs.tsv
│       ├── Sample_D.stringent
│       │   ├── Sample_D_stringent_false_negatives.tsv
│       │   ├── Sample_D_stringent_false_positives.tsv
│       │   ├── Sample_D_stringent_genotype_confusion.tsv
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
    ├── 06_evaluation.sh
    └── rerun_adapter_correction.sh

33 directories, 245 files
```

### Results

| Sample | Filter | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | raw | 4858 | 0 | 142 | 1.0000 | 0.9716 | 0.9856 |
| A | gentle | 4855 | 0 | 145 | 1.0000 | 0.9710 | 0.9853 |
| A | strict | 4747 | 0 | 253 | 1.0000 | 0.9494 | 0.9740 |
| B | raw | 4426 | 22 | 574 | 0.9951 | 0.8852 | 0.9369 |
| B | gentle | 4422 | 0 | 578 | 1.0000 | 0.8844 | 0.9387 |
| B | strict | 4327 | 0 | 673 | 1.0000 | 0.8654 | 0.9278 |
| C | raw | 4862 | 0 | 138 | 1.0000 | 0.9724 | 0.9860 |
| C | gentle | 4857 | 0 | 143 | 1.0000 | 0.9714 | 0.9855 |
| C | strict | 4747 | 0 | 253 | 1.0000 | 0.9494 | 0.9740 |
| D | raw | 4844 | 11 | 156 | 0.9977 | 0.9688 | 0.9831 |
| D | gentle | 4841 | 1 | 159 | 0.9998 | 0.9682 | 0.9837 |
| D | strict | 4736 | 0 | 264 | 1.0000 | 0.9472 | 0.9729 |


## Assignment Analysis - Questions to answer

Minimum one paragraph each.

> How did samples differ from each other in read count versus in read quality? Why does each matter for SNP calling?

There generally didn't seem to be a correlation between read counts and quality. Sample B generally had very low read counts but relatively high quality scores whereas sample D had as high read counts as A and C but much lower quality. For SNP calling, having a minimum read counts are important to distinguish a potential SNP statistically from noise. Since any statistical procedure like the ones underlying SNP calling need a minimum number of samples regardless of quality to make a call. Although generally high quality reads may offset the minimum number of reads needed to confidently call a SNP. Judging from the final variant calling results though, read count seems much more impactful than read quality because the F1 scores of Sample D do not suffer much despite its generally low read quality. B on the other hand seems to have pretty low recall and thus F1 in spite of the higher quality of its relatively few reads.


> What are the relative merits and drawbacks of filtering? Interpret the above chart in terms of the original quality issues, and how those yielded different outcomes. When is filtering a good idea vs when it is it a bad idea? Think specifically about false positives vs false negatives.

The main effect of filtering seems to be to reduce the rate of False positives. Which is predicated on the presence of false positives in your data and the degree to which the filters capture them. Filtering may help for data with lots of false positives. However for B the low prevalence of false positives and tendency of the filtering algorithm to overassign false negatives mean that filtering is ultimately detrimental in all cases. For samples with low read depth like B, the filtering thresholds remove what little signal there is whereas doing best effort calling on the poor data yields a good F1 score of 0.83. 

> Pick **two false negatives** and **two false positives**, from any sample. For each, tell me the sample, tell me the chromosome and coordinate location, put in a screenshot of IGV or JBrowse, and explain the specific cause of the FN or FP call.

| Example data source | Calling method |
| --- | --- |
| `archive/before_pe_adapter_20261015_132613/` | Original separate-sample calling |

### False negative 1

| Sample | Filter | Chromosome | Position (1-based) | REF → ALT |
| --- | --- | --- | --- | --- |
| A | stringent | I (`NC_001133.9`) | 205684 | T → C |

Truth genotype: `0/1`. Raw call: `GT=0/1`, `QUAL=52.38`, `INFO/DP=12`, `INFO/MQ=23`. Stringent call: absent.

![False negative 1: sample A, 205684](https://raw.githubusercontent.com/Bryan-Nsoh/EPP622-2026-SNPCalling/main/images/A_FN.png)

**Explanation:**

This is a case of False negative after overenthusiastic filtering. This sample clearly as a high read depth with the SNP appearing in a decent number of reads. The variant is correctly called without any filtering but was excluded because it did not meet the read mapping quality threshold of 40(INFO/MQ>=40).

### False negative 2

| Sample | Filter | Chromosome | Position (1-based) | REF → ALT |
| --- | --- | --- | --- | --- |
| C | raw | I (`NC_001133.9`) | 160949 | G → A |

Truth genotype: `0/1`. Raw call: absent.

![False negative 2: sample C, 160949](https://raw.githubusercontent.com/Bryan-Nsoh/EPP622-2026-SNPCalling/main/images/C_FN.png)

**Explanation:**

Unlike the case above which has the alternative allele showing up in many reads the alternative allele, A here only shows up in a single read. So there is weak support for the underlying allele. Not sure what could cause this perhaps more reads would help. Though you would think this many reads should be sufficient.

### False positive 1

| Sample | Filter | Chromosome | Position (1-based) | REF → ALT |
| --- | --- | --- | --- | --- |
| B | raw | II (`NC_001134.8`) | 14817 | T → A |

Called genotype: `1/1`. `QUAL=5.8`, `FORMAT/DP=1` (evaluator values). Truth variant: absent.

![False positive 1: sample B, 14817](https://raw.githubusercontent.com/Bryan-Nsoh/EPP622-2026-SNPCalling/main/images/B_FP.png)

**Explanation:**

The extremely low read depth and poor quality QUAL=5.8 make this unlikely to be a SNP and in general make any calling based on this data highly uncertain. The cause of this false positive could largely be due to lack of filtering based on quality and read depth

### False positive 2

| Sample | Filter | Chromosome | Position (1-based) | REF → ALT |
| --- | --- | --- | --- | --- |
| D | raw | III (`NC_001135.5`) | 66244 | T → G |

Called genotype: `0/1`. `QUAL=3.7`, `FORMAT/DP=13` (evaluator values). Truth variant: absent.

![False positive 2: sample D, 66244](https://raw.githubusercontent.com/Bryan-Nsoh/EPP622-2026-SNPCalling/main/images/D_FP.png)

**Explanation:**

The ratio of alternative allele to the reference nucleotide is 3 to 10. this combined with the relatively low quality scores of QUAL=3.7 reduce the likelihood of this being a SNP. Since the raw calls lacked additional quality filtering however this was wrongly called.



## AI Analysis - Pipeline Audit

This is an evaluation of the pipeline you generated at the beginning of the learning module. Provide full sentences and thoughtful rationale.

Model used: OpenAI Codex, GPT 5.6 sol

Date generated: September 22, 2026

Exact prompt: Write me a complete short-read SNP calling pipeline for a haploid fungal genome, starting from paired FASTQ and ending at a filtered VCF.


Minimum one paragraph each.

> Three things it got right.

- The full selection of tools needed to run the full pipeline
- It considered reusability, making use of variables to hold paths, foldernames etc so the pipeline could be run by someone with a different setup with a little tweaking
- I noticed AIs always use `set -euo pipefail` at the start of every bash script. When i searched online i learned the command causes the entire pipeline to crash and print any error instead of failing silently or inscrutably then proceeding with potentially erroneous outputs to the end. Errors without error messages are the most misleading and the absolute worst to debug. I think all bash scripts should use this and have started using this myself.

> Three things you don't understand and would want to learn more about before using. Or do you understand all of it? If so be ready to answer questions about it.

- The AI used GATK-4 which is a genome analysis toolkit to index the reference genome, to call variants, mark PCR duplicates etc. i hadnt heard of this software before i guess its an alternative to the tools we use but given that we havent discussed GATK in class and i hadn't heard of it, every part of the pipeline that used GATK-4 is difficult to understand
- Some of the bash the AI used is a little too advanced/opaque for me for instance `RG="@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tLB:${SAMPLE}_lib1\\tPU:${SAMPLE}_unit1"`. i could probably figure out what this does if i really tried but it would be painful. There is almost certainly a more human-legible way to achieve whatever is going on here.
- Usage lifecycle of the pipeline could also be more well thought out. Managing and organizing the inputs/outputs appeared to have received some thought but seems like a secondary concern. 

These issues could all be resolved with a better prompt or more specific requests, but you wouldn't know what to ask for in the first place if you didn't understand the tools and their purpose

> How specifically would you assess the results of this pipeline for accuracy?

For intermediary steps I would use multiqc and its html viewer to look at the data before and after filtering specifically looking for read depth phred scores and examining the effects of the thresholds used for filtering reads based on these measures. I would also reuse the python tool provided to obtain F1 scores and precision recall metrics for variant calling ideally running this against another pipeline i know and trust. If the F1 scores are good after a little tuning of the thresholds used for VCF filtering i would trust it more. 

> What did the labs teach you that turned out to matter for evaluating the AI-generated code?

AI written bash can be rather inscrutable especially when written for unfamiliar CLIs. What helped me understand and evaluate the AI code better was familiarity with the CLIs from having worked with them myself.

> What is the value of learning about a computational genomics task before using AI for it — or is there none?

Using AI for execution without understanding could lead to alot of wasted effort and unnecessary confusion. AI tends to be helpful when you can ask very precise questions of it and have the ability to judge its output and the only way to get to that level is to spend the time learning and practicing.
