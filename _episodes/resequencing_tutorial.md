---
layout: page
title: Resequencing and Variant Calling Tutorial
published: true
---

> ## Paper Reading
> Please read this paper before class:
> [Van der Auwera & O'Connor (2020) Genomics in the Cloud (GATK Best Practices). O'Reilly Media](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/)
>
> Also see the GATK Best Practices:
> [GATK Best Practices for Germline SNPs & Indels](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932)
{: .prereq}

---

## Overview: Resequencing and Variant Calling

**Resequencing** (also called whole-genome resequencing, WGS, or whole-exome sequencing, WES) aligns short reads from an individual's genome back to a reference genome to discover **genetic variants** — differences between the sample and the reference.

### Types of Variants

| Variant Type | Size | Description |
|-------------|------|-------------|
| SNP | 1 bp | Single nucleotide polymorphism |
| InDel | 1–50 bp | Insertion or deletion |
| CNV | kb–Mb | Copy number variation |
| SV | > 50 bp | Structural variant (inversion, translocation) |

### Applications

- Population genetics / GWAS
- Disease gene discovery (human genetics)
- Crop breeding (QTL, marker-assisted selection)
- Cancer genomics (somatic mutation detection)
- Phylogenomics (SNP-based trees)

### Resequencing vs. RNA-Seq

| Feature | Resequencing (WGS) | RNA-Seq |
|---------|-------------------|---------|
| Input | Genomic DNA | RNA (cDNA) |
| Coverage | Uniform (genome) | Non-uniform (expressed genes) |
| Aligner | BWA-MEM2, Bowtie2 | HISAT2, STAR (splice-aware) |
| Output | Variants (VCF) | Expression counts |

---

## The Variant Calling Workflow

<img src="{{ page.root }}/fig/reseq_workflow.svg" alt="Resequencing and Variant Calling Workflow" style="max-width:520px; width:100%;"/>

---

## 0. Prerequisites

> **Disk Space:** This tutorial requires at least **30 GB of free disk space**. Each sample generates ~5 GB of intermediate files (BAM, markdup BAM, recal BAM, GVCF). Check your available space before starting:
>
> ```bash
> df -h .
> ```
>
> If disk space is limited, follow the cleanup tips marked with the label **[Cleanup]** after each major step.
{: .callout}

---

## 1. Environment Setup

### Create Conda Environment

```bash
conda create -n reseq -c bioconda -c conda-forge python=3.11
conda activate reseq

conda install -c bioconda -c conda-forge fastqc fastp bwa-mem2 samtools
conda install -c bioconda -c conda-forge openjdk=17 picard gatk4 bcftools
conda install -c bioconda -c conda-forge snpeff plink tabix
pip install multiqc

# Clean conda package cache to free disk space
conda clean --all -y
```

### Verify Installations

```bash
fastp --version
bwa-mem2 version
samtools --version
gatk --version
bcftools --version
```

### Create Working Directory

```bash
mkdir -p ~/bch709/reseq
cd ~/bch709/reseq
```

---

## 2. Quality Control

### Download Example Data

For this tutorial, we use publicly available *Arabidopsis thaliana* WGS datasets from the 1001 Genomes Project. We download **two samples** so that we can demonstrate joint genotyping later.

> **Note:** The SRA Toolkit (`fasterq-dump`, `fastq-dump`) is known to cause segmentation faults on WSL (Windows Subsystem for Linux). Downloading FASTQ files directly from ENA avoids this issue entirely.
{: .callout}

**Option A — Download from ENA (European Nucleotide Archive):**

```bash
cd ~/bch709/reseq

# Sample 1: Col-0 re-sequencing (SRR519585)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519585/SRR519585_1.fastq.gz -O sample1_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519585/SRR519585_2.fastq.gz -O sample1_R2.fastq.gz

# Sample 2: TBO-01 re-sequencing (SRR519586)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519586/SRR519586_1.fastq.gz -O sample2_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519586/SRR519586_2.fastq.gz -O sample2_R2.fastq.gz

ls -lh
```

**Option B — Download from Dropbox mirror (if ENA is slow or fails):**

```bash
cd ~/bch709/reseq

# Sample 1
wget "https://www.dropbox.com/scl/fi/8fy6hrhczh7v23ojqox5r/wgs_R1.fastq.gz?rlkey=nmanzksfbtm76xer4jesq2vuy&dl=1" -O sample1_R1.fastq.gz
wget "https://www.dropbox.com/scl/fi/kt29bt8c29vf9i294xeek/wgs_R2.fastq.gz?rlkey=8xm8lk839jqcz2qfjlzx4gqle&dl=1" -O sample1_R2.fastq.gz

# Sample 2 — download from ENA (no Dropbox mirror)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519586/SRR519586_1.fastq.gz -O sample2_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519586/SRR519586_2.fastq.gz -O sample2_R2.fastq.gz

ls -lh
```

> **Dropbox direct download tip:** Add `&dl=1` to the end of a Dropbox share URL, and wrap the URL in double quotes to prevent the shell from interpreting `&` and `?` characters.
{: .callout}

> You can find ENA download links for any SRA accession at [ENA Browser](https://www.ebi.ac.uk/ena/browser/) or [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra).

### Run FastQC

```bash
# -t 4 : use 4 threads
fastqc -t 4 sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz

# Aggregate all QC reports into one
multiqc .
```

### Key Metrics for WGS

| Metric | What to Check |
|--------|--------------|
| Per-base quality | Should be Q30+ across most of the read |
| GC content | Should match species genome GC% |
| Duplication | High duplication may indicate low-complexity library |
| Adapter content | Remove if present |

---

## 3. Read Trimming

**fastp** performs adapter removal, quality trimming, and QC in a single pass — faster and simpler than running separate tools.

| Option | Description |
|--------|-------------|
| `--in1` / `--in2` | Forward / reverse input reads |
| `--out1` / `--out2` | Forward / reverse output reads |
| `--detect_adapter_for_pe` | Auto-detect adapters for paired-end data |
| `--qualified_quality_phred 20` | Minimum base quality threshold (Q20) |
| `--length_required 50` | Discard reads shorter than 50 bp |
| `--thread 4` | Number of threads |
| `--html` / `--json` | QC report outputs (JSON is used by MultiQC) |

```bash
mkdir -p trim

# Trim sample 1
fastp \
  --in1 sample1_R1.fastq.gz \
  --in2 sample1_R2.fastq.gz \
  --out1 trim/sample1_R1_trimmed.fq.gz \
  --out2 trim/sample1_R2_trimmed.fq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 4 \
  --html trim/sample1_fastp_report.html \
  --json trim/sample1_fastp_report.json

# Trim sample 2
fastp \
  --in1 sample2_R1.fastq.gz \
  --in2 sample2_R2.fastq.gz \
  --out1 trim/sample2_R1_trimmed.fq.gz \
  --out2 trim/sample2_R2_trimmed.fq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 4 \
  --html trim/sample2_fastp_report.html \
  --json trim/sample2_fastp_report.json

multiqc trim/ -n trim_report
```

> **[Cleanup]** After trimming, you can remove the original FASTQ files to save disk space:
> ```bash
> rm -f sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz
> ```
{: .callout}

---

## 4. Reference Genome Preparation

### Download Reference

```bash
# Arabidopsis TAIR10 reference from Ensembl Plants
wget https://ftp.ensemblgenomes.org/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz --no-check-certificate
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa reference.fasta
```

> **Note:** The Ensembl TAIR10 reference uses chromosome names `1`, `2`, `3`, `4`, `5`, `Mt`, `Pt` (not `Chr1`, `Chr2`, etc.). Keep this in mind when specifying genomic regions in later steps.
{: .callout}

### Create BWA-MEM2 Index

```bash
# Generates .0123, .amb, .ann, .bwt.2bit.64, .pac
bwa-mem2 index reference.fasta
```

### Create FASTA Index and Sequence Dictionary (required by GATK)

```bash
# .fai index for samtools/GATK
samtools faidx reference.fasta

# .dict sequence dictionary for GATK
picard CreateSequenceDictionary R=reference.fasta O=reference.dict
```

---

## 5. Alignment with BWA-MEM2

BWA-MEM2 uses the Burrows-Wheeler Transform (BWT) + FM-index for short-read alignment.

### Why BWA-MEM2 for DNA (not HISAT2)?

| Feature | BWA-MEM2 | HISAT2/STAR |
|---------|---------|------------|
| Splice-aware | No | Yes |
| DNA reads | Optimal | Not designed for |
| RNA reads | Suboptimal | Optimal |

### Align Reads

The `@RG` (read group) tag is **required** by GATK.

| @RG Field | Description |
|-----------|-------------|
| `ID` | Read group ID |
| `SM` | Sample name (GATK uses this to identify samples) |
| `PL` | Sequencing platform (e.g., ILLUMINA) |
| `LB` | Library name |
| `PU` | Platform unit (e.g., flowcell-barcode.lane) |

```bash
# Align sample 1
bwa-mem2 mem \
  -t 4 \
  -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  reference.fasta \
  trim/sample1_R1_trimmed.fq.gz \
  trim/sample1_R2_trimmed.fq.gz \
  | samtools sort -@ 4 -o sample1.bam

samtools index sample1.bam

# Align sample 2
bwa-mem2 mem \
  -t 4 \
  -R "@RG\tID:sample2\tSM:sample2\tPL:ILLUMINA\tLB:lib2\tPU:unit2" \
  reference.fasta \
  trim/sample2_R1_trimmed.fq.gz \
  trim/sample2_R2_trimmed.fq.gz \
  | samtools sort -@ 4 -o sample2.bam

samtools index sample2.bam
```

### Alignment Statistics

```bash
samtools flagstat sample1.bam
samtools flagstat sample2.bam

samtools stats sample1.bam > sample1.stats
samtools stats sample2.bam > sample2.stats
multiqc . -n alignment_report
```

| Metric | What to Expect |
|--------|---------------|
| Mapping rate | > 95% for matched reference |
| Properly paired | > 90% |
| Average depth | 10–30x for population WGS; 30–60x for clinical |

> **[Cleanup]** After alignment, you can remove the trimmed FASTQ files:
> ```bash
> rm -rf trim/
> ```
{: .callout}

---

## 6. Mark Duplicates

PCR amplification creates duplicate reads. For variant calling, these **must be marked** (not just discarded — GATK handles them).

### Run Picard MarkDuplicates

| Option | Description |
|--------|-------------|
| `I=` | Input BAM |
| `O=` | Output BAM with duplicates flagged |
| `M=` | Duplication metrics file |
| `VALIDATION_STRINGENCY=SILENT` | Suppress warnings on BAM format |

```bash
# Mark duplicates for sample 1
picard MarkDuplicates \
  I=sample1.bam \
  O=sample1.markdup.bam \
  M=sample1.markdup.metrics \
  VALIDATION_STRINGENCY=SILENT

samtools index sample1.markdup.bam

# Mark duplicates for sample 2
picard MarkDuplicates \
  I=sample2.bam \
  O=sample2.markdup.bam \
  M=sample2.markdup.metrics \
  VALIDATION_STRINGENCY=SILENT

samtools index sample2.markdup.bam
```

### Check Duplication Rate

```bash
cat sample1.markdup.metrics
cat sample2.markdup.metrics
multiqc . -n markdup_report
```

> **Note:** Duplication rates > 30–40% may indicate problems with input DNA quality or library complexity. For very high duplication, use a PCR-free library prep.

> **[Cleanup]** After marking duplicates, remove the original sorted BAM files:
> ```bash
> rm -f sample1.bam sample1.bam.bai sample2.bam sample2.bam.bai
> ```
{: .callout}

---

## 7. Base Quality Score Recalibration (BQSR)

GATK BQSR corrects systematic errors in base quality scores from the sequencer. It requires a **known variants VCF** (e.g., dbSNP).

### Download Known Variants

BQSR requires a VCF of known polymorphic sites to distinguish true variants from sequencing errors.

> **Important:** The 1001 Genomes multi-sample VCF is very large (~5 GB compressed, ~80 GB decompressed) and contains formatting inconsistencies that cause GATK errors. For BQSR, we only need the **variant positions** (sites-only), not the genotype data. The commands below extract a sites-only VCF that is much smaller and works correctly with GATK.
{: .callout}

```bash
# For Arabidopsis: download from the 1001 Genomes Project
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz

# Extract sites-only VCF (removes genotype columns — BQSR only needs positions)
# This also fixes malformed lines in the original VCF.
zcat 1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
  | awk 'BEGIN{OFS="\t"} /^##/{print; next} /^#CHROM/{print $1,$2,$3,$4,$5,$6,$7,$8; next} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
  | bgzip > known_sites.vcf.gz

tabix -p vcf known_sites.vcf.gz

# Remove the large original file
rm -f 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

# For human (hg38): download dbSNP from GATK resource bundle
# wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
```

> **No known variants available?** For less-characterized species, skip BQSR or use an iterative "bootstrap" approach: call variants with HaplotypeCaller, use those as known sites, then re-run BQSR. Two rounds are usually sufficient.
{: .callout}

### Step 1: Compute Recalibration Table

| Option | Description |
|--------|-------------|
| `-I` | Input BAM (with duplicates marked) |
| `-R` | Reference genome |
| `--known-sites` | Known polymorphic sites VCF |
| `-O` | Output recalibration table |

```bash
# Sample 1
gatk BaseRecalibrator \
  -I sample1.markdup.bam \
  -R reference.fasta \
  --known-sites known_sites.vcf.gz \
  -O sample1.recal.table

# Sample 2
gatk BaseRecalibrator \
  -I sample2.markdup.bam \
  -R reference.fasta \
  --known-sites known_sites.vcf.gz \
  -O sample2.recal.table
```

### Step 2: Apply Recalibration

```bash
# Sample 1
gatk ApplyBQSR \
  -I sample1.markdup.bam \
  -R reference.fasta \
  --bqsr-recal-file sample1.recal.table \
  -O sample1.recal.bam

samtools index sample1.recal.bam

# Sample 2
gatk ApplyBQSR \
  -I sample2.markdup.bam \
  -R reference.fasta \
  --bqsr-recal-file sample2.recal.table \
  -O sample2.recal.bam

samtools index sample2.recal.bam
```

> **[Cleanup]** After BQSR, remove the markdup BAM files:
> ```bash
> rm -f sample1.markdup.bam sample1.markdup.bam.bai sample2.markdup.bam sample2.markdup.bam.bai
> ```
{: .callout}

---

## 8. Variant Calling with GATK HaplotypeCaller

GATK HaplotypeCaller performs local assembly of haplotypes to call SNPs and indels simultaneously.

### Call Variants (per-sample GVCF mode)

Using **GVCF mode** allows joint genotyping across multiple samples — recommended for multi-sample projects.

| Option | Description |
|--------|-------------|
| `-R` | Reference genome |
| `-I` | Input recalibrated BAM |
| `-O` | Output GVCF (genomic VCF) |
| `-ERC GVCF` | Emit reference confidence — produces GVCF instead of plain VCF |
| `--native-pair-hmm-threads` | Threads for PairHMM calculation |

```bash
# Sample 1
gatk HaplotypeCaller \
  -R reference.fasta \
  -I sample1.recal.bam \
  -O sample1.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads 4

# Sample 2
gatk HaplotypeCaller \
  -R reference.fasta \
  -I sample2.recal.bam \
  -O sample2.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads 4
```

### Joint Genotyping

```bash
# Combine GVCFs from both samples
gatk CombineGVCFs \
  -R reference.fasta \
  -V sample1.g.vcf.gz \
  -V sample2.g.vcf.gz \
  -O cohort.g.vcf.gz

# Joint genotyping across all samples
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V cohort.g.vcf.gz \
  -O cohort.vcf.gz
```

> **[Cleanup]** After joint genotyping, you can remove individual GVCFs and recal BAMs:
> ```bash
> rm -f sample1.g.vcf.gz sample1.g.vcf.gz.tbi sample2.g.vcf.gz sample2.g.vcf.gz.tbi
> rm -f cohort.g.vcf.gz cohort.g.vcf.gz.tbi
> rm -f sample1.recal.bam sample1.recal.bai sample2.recal.bam sample2.recal.bai
> ```
{: .callout}

---

## 9. Variant Filtering

### Separate SNPs and Indels

```bash
# Extract SNPs only
gatk SelectVariants \
  -R reference.fasta \
  -V cohort.vcf.gz \
  --select-type-to-include SNP \
  -O cohort.snps.vcf.gz

# Extract Indels only
gatk SelectVariants \
  -R reference.fasta \
  -V cohort.vcf.gz \
  --select-type-to-include INDEL \
  -O cohort.indels.vcf.gz
```

### Hard Filter SNPs (GATK Best Practices)

```bash
gatk VariantFiltration \
  -R reference.fasta \
  -V cohort.snps.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O cohort.snps.filtered.vcf.gz
```

### Hard Filter Indels

```bash
gatk VariantFiltration \
  -R reference.fasta \
  -V cohort.indels.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O cohort.indels.filtered.vcf.gz
```

### Key GATK Filters Explained

| Filter | Field | Threshold (SNP) | Threshold (Indel) | Meaning |
|--------|-------|-----------------|-------------------|---------|
| QD | Quality by Depth | < 2.0 | < 2.0 | Normalized variant quality |
| FS | Fisher Strand Bias | > 60.0 | > 200.0 | Strand bias test |
| MQ | RMS Mapping Quality | < 40.0 | — | Average mapping quality |
| MQRankSum | MQ Rank Sum | < -12.5 | — | Comparison of MQ for ref/alt |
| ReadPosRankSum | Read Position RS | < -8.0 | < -20.0 | Position bias within reads |

---

## 10. VCF Format

### VCF File Structure

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO              FORMAT     SAMPLE1
1       12345   .       A    T    220   PASS    DP=45;AF=0.5;...  GT:AD:DP   0/1:22,23:45
```

### Key VCF Fields

| Field | Description |
|-------|-------------|
| CHROM | Chromosome |
| POS | Position (1-based) |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| QUAL | Variant quality score |
| FILTER | PASS or filter names |
| INFO | Variant annotations |
| FORMAT | Sample field format |
| GT | Genotype (0/0=hom-ref, 0/1=het, 1/1=hom-alt) |

### Basic VCF Manipulation with BCFtools

| Option | Description |
|--------|-------------|
| `-f PASS` | Keep only variants that passed all filters |
| `-O z` | Output as compressed VCF (`.vcf.gz`) |
| `-o` | Output file name |

```bash
# Count variants
bcftools stats cohort.snps.filtered.vcf.gz | grep "^SN"

# Keep only PASS variants
bcftools view -f PASS -O z -o cohort.snps.pass.vcf.gz cohort.snps.filtered.vcf.gz

# Extract specific genomic region
# Note: Ensembl TAIR10 uses "1", "2", ... (not "Chr1", "Chr2")
bcftools view cohort.snps.pass.vcf.gz 1:100000-200000

# Get variant summary statistics
bcftools stats cohort.snps.pass.vcf.gz > stats.txt
multiqc . -n vcf_report
```

---

## 11. Variant Annotation with SnpEff

SnpEff predicts the functional effect of each variant (missense, nonsense, synonymous, etc.).

### Setup SnpEff Database

```bash
# List available databases
snpEff databases | grep -i arabidopsis

# Download database (example: Arabidopsis TAIR10)
# Note: the exact database name may vary by SnpEff version.
# Use the output of the command above to find the correct name.
snpEff download athalianaTair10
```

### Annotate Variants

```bash
snpEff -v athalianaTair10 cohort.snps.pass.vcf.gz > cohort.snps.annotated.vcf

# snpEff also generates snpEff_summary.html and snpEff_genes.txt
```

### Variant Effect Categories

| Effect | Description |
|--------|-------------|
| HIGH | Frameshift, stop gained/lost, splice site |
| MODERATE | Missense, in-frame InDel |
| LOW | Synonymous, splice region |
| MODIFIER | Intergenic, intronic, upstream/downstream |

---

## 12. Population-Level Analysis

For multiple samples, use population genomics tools:

### PCA (Principal Component Analysis)

| Option | Description |
|--------|-------------|
| `--vcf` | Input VCF file |
| `--make-bed` | Output PLINK binary format (`.bed`/`.bim`/`.fam`) |
| `--pca 10` | Compute top 10 principal components |
| `--allow-extra-chr` | Allow non-human chromosome names (required for plants: Arabidopsis uses `1`–`5`, `Mt`, `Pt`) |

```bash
# Convert VCF to PLINK binary format
plink --vcf cohort.snps.pass.vcf.gz --make-bed --out cohort --allow-extra-chr

# Run PCA (top 10 components)
plink --bfile cohort --pca 10 --out cohort.pca --allow-extra-chr
```

### Linkage Disequilibrium Pruning

Before running GWAS or population structure analyses, prune variants in LD to avoid biasing results.

| `--indep-pairwise` Parameter | Description |
|------------------------------|-------------|
| `50` | Window size (number of SNPs) |
| `10` | Step size (SNPs to shift window) |
| `0.2` | r² threshold (remove one of a pair if r² > 0.2) |

```bash
# Identify LD-pruned variant set
plink --bfile cohort --indep-pairwise 50 10 0.2 --out cohort.ld --allow-extra-chr

# Extract pruned variants
plink --bfile cohort --extract cohort.ld.prune.in --make-bed --out cohort.pruned --allow-extra-chr
```

### GWAS

For GWAS workflows, see the [GWAS tutorial]({{site.baseurl}}/episodes/GWAS/).

---

## 13. Full Workflow Summary

| Step | Tool | Input | Output |
|------|------|-------|--------|
| QC | FastQC + MultiQC | FASTQ | HTML report |
| Trim | fastp | FASTQ | Trimmed FASTQ |
| Align | BWA-MEM2 | FASTQ + Reference | BAM |
| Sort/Index | SAMtools | BAM | Sorted BAM |
| Mark duplicates | Picard | BAM | Markdup BAM |
| BQSR | GATK | BAM + known VCF | Recal BAM |
| Call variants | GATK HaplotypeCaller | BAM | GVCF |
| Joint genotype | GATK GenotypeGVCFs | GVCFs | VCF |
| Filter | GATK VariantFiltration | VCF | Filtered VCF |
| Annotate | SnpEff | VCF | Annotated VCF |

---

## 14. Cleanup

```bash
conda deactivate
# Optional: remove environment when done
conda env remove --name reseq
```

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| `No space left on device` | Intermediate files filling disk | Follow **[Cleanup]** tips after each step; run `conda clean --all -y` |
| `libcrypto.so.1.0.0: cannot open` | OpenSSL version mismatch in conda | `ln -s ${CONDA_PREFIX}/lib/libcrypto.so.3 ${CONDA_PREFIX}/lib/libcrypto.so.1.0.0` (check which `libcrypto.so.*` exists first) |
| GATK `VCF is malformed` on known variants | Multi-sample VCF has inconsistent genotype columns | Use the sites-only extraction command in Section 7 |
| `Contig 'X' is not defined in the header` | Chromosome naming mismatch between VCF and reference | Verify chromosome names match: `head reference.fasta.fai` vs `zcat file.vcf.gz \| grep -v '^#' \| cut -f1 \| sort -u` |
| BWA-MEM2 `Unexpected end of file` | Index files corrupted (often from disk-full during indexing) | Delete index files (`reference.fasta.0123`, `.bwt.2bit.64`, etc.) and re-run `bwa-mem2 index` |
| `fasterq-dump` segfault on WSL | Known SRA Toolkit bug on WSL | Download FASTQ directly from ENA instead |

---

## References

| Resource | Link |
|---------|------|
| GATK Best Practices | [broadinstitute.github.io/gatk](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) |
| BWA-MEM2 paper | [Md et al. 2019, iScience](https://www.sciencedirect.com/science/article/pii/S2589004219310017) |
| fastp paper | [Chen et al. 2018, Bioinformatics](https://doi.org/10.1093/bioinformatics/bty560) |
| SAM flag decoder | [Picard SAM Flags](https://broadinstitute.github.io/picard/explain-flags.html) |
| VCF specification | [samtools.github.io/hts-specs](https://samtools.github.io/hts-specs/VCFv4.2.pdf) |
| SnpEff documentation | [pcingola.github.io/SnpEff](https://pcingola.github.io/SnpEff/) |
| BCFtools manual | [samtools.github.io/bcftools](https://samtools.github.io/bcftools/bcftools.html) |
