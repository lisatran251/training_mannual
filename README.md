# Bioinformatics Training Session

## Overview
This training session is designed to introduce participants to essential bioinformatics tools and workflows using Conda for environment management. Participants will learn how to process sequencing data, perform alignment, and call variants using popular bioinformatics software.

## Setup Instructions

### Step 1: Create and Activate the Conda Environment
To begin, create a new Conda environment specifically for this training session:

```bash
conda create -n bioenv
```
*Creates a new Conda environment named `bioenv`.*

Next, initialize Conda:

```bash
conda init
```
*Initializes Conda for your shell, enabling the use of Conda commands.*

Activate the newly created environment:

```bash
conda activate bioenv
```
*Activates the `bioenv` environment, allowing you to use the installed packages.*

### Step 2: Install Required Packages
Install the necessary Python packages and bioinformatics tools:

```bash
conda install pandas matplotlib seaborn
```
*Installs data manipulation and visualization libraries: Pandas, Matplotlib, and Seaborn.*

```bash
conda install -c bioconda bwa fastqc fastp samtools bcftools sra-tools pysam
```
*Installs bioinformatics tools from the Bioconda channel: BWA (alignment), FastQC (quality control), Fastp (trimming), Samtools (file manipulation), BCFtools (variant calling), SRA-tools (data retrieval), and Pysam (Python interface for SAM/BAM files).*

## Data Processing Workflow

### Step 3: Download and Prepare Data
Use the SRA Toolkit to download the sequencing data:

```bash
prefetch SRR396637
```
*Downloads the SRA file for the specified accession number (SRR396637) from the Sequence Read Archive.*

Convert the downloaded SRA file into FASTQ format:

```bash
fastq-dump --split-files SRR396637/SRR396637.sra
```
*Converts the SRA file into FASTQ format, splitting paired-end reads into separate files.*

### Step 4: Quality Control and Trimming
Perform quality control and trimming of the FASTQ files using `fastp`:

```bash
fastp -i SRR396637_1.fastq -I SRR396637_2.fastq -o trimmed_SRR396637_R1.fastq -O trimmed_SRR396637_R2.fastq -q 20 -l 50 --cut_mean_quality 20 --cut_right
```
*Trims the input FASTQ files (`SRR396637_1.fastq` and `SRR396637_2.fastq`), producing trimmed output files. Parameters:*
- `-q 20`: Minimum quality score for bases to keep.
- `-l 50`: Minimum length of reads after trimming.
- `--cut_mean_quality 20`: Trims reads based on the mean quality score.
- `--cut_right`: Trims the right end of the reads.

### Step 5: Alignment
Index the reference genome:

```bash
bwa index reference.fna
```
*Indexes the reference genome file (`reference.fna`) for efficient alignment.*

Align the trimmed reads to the reference genome:

```bash
bwa mem reference.fna trimmed_SRR396637_R1.fastq trimmed_SRR396637_R2.fastq > aligned_reads.sam
```
*Aligns the trimmed FASTQ files to the reference genome, outputting the results in SAM format.*

### Step 6: Convert and Sort SAM to BAM
Convert the SAM file to BAM format and sort it:

```bash
samtools view -Sb aligned_reads.sam | samtools sort -o aligned_reads_sorted.bam
```
*Converts the SAM file to BAM format (`-Sb`), then sorts the BAM file, outputting the sorted file as `aligned_reads_sorted.bam`.*

### Step 7: Variant Calling
Call variants using `bcftools`:

```bash
bcftools mpileup -f reference.fna aligned_reads_sorted.bam | bcftools call -mv -Oz -o variants.vcf
```
*Generates a variant call format (VCF) file from the sorted BAM file. Parameters:*
- `-f reference.fna`: Specifies the reference genome for variant calling.
- `-mv`: Calls variants and outputs both SNPs and indels.
- `-Oz`: Compresses the output VCF file using bgzip.

### Step 8: Visualization
To visualize the results, run the following Python script:

```bash
python3 visualization.py
```
*Executes a Python script (`visualization.py`) that may analyze or visualize the results of the bioinformatics workflow.*

## Conclusion
By following this workflow, participants will gain hands-on experience with bioinformatics tools and techniques essential for analyzing sequencing data. For any questions or further assistance, please feel free to reach out.
