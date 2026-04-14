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
| Aligner | BWA-MEM, Bowtie2 | HISAT2, STAR (splice-aware) |
| Output | Variants (VCF) | Expression counts |

---

## The Variant Calling Workflow

```
Raw FASTQ (DNA reads)
    │
    ├── FastQC → MultiQC (QC)
    ├── Trim Galore / fastp (trimming)
    │
    ├── BWA-MEM2 (alignment to reference genome)
    │       │
    │       └── SAMtools sort + index
    │
    ├── Picard MarkDuplicates (remove PCR duplicates)
    │
    ├── GATK BaseRecalibrator + ApplyBQSR (BQSR)
    │
    ├── GATK HaplotypeCaller (variant calling → GVCF)
    │
    ├── GATK GenotypeGVCFs (joint genotyping)
    │
    ├── GATK VariantFiltration (hard filtering)
    │
    └── SnpEff / VEP (variant annotation)
```

---

## 1. Environment Setup

### Create Conda Environment

```bash
conda create -n reseq -c bioconda -c conda-forge python=3.11 \
  fastqc trim-galore bwa-mem2 samtools picard gatk4 \
  snpeff multiqc bcftools plink
conda activate reseq
```

### Verify Installations

```bash
bwa-mem2 version
samtools --version
gatk --version
```

### Create Working Directory

```bash
mkdir -p ~/bch709/reseq
cd ~/bch709/reseq
```

---

## 2. Quality Control

### Download Example Data

For this tutorial, we use a publicly available *Arabidopsis thaliana* WGS dataset from the 1001 Genomes Project. Download using SRA tools:

```bash
conda install -c bioconda sra-tools
cd ~/bch709/reseq

# Arabidopsis accession Col-0 re-sequencing (SRR519585)
fasterq-dump --split-files SRR519585 -O .
mv SRR519585_1.fastq wgs_R1.fastq
mv SRR519585_2.fastq wgs_R2.fastq
gzip wgs_R1.fastq wgs_R2.fastq
ls -lh
```

> Replace `SRR519585` with any SRA accession from your experiment. Browse datasets at [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra).

### Run FastQC

```bash
fastqc -t 4 wgs_R1.fastq.gz wgs_R2.fastq.gz
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

```bash
trim_galore \
  --paired \
  --cores 4 \
  --fastqc \
  --gzip \
  -o trim \
  wgs_R1.fastq.gz wgs_R2.fastq.gz

multiqc --dirs ~/bch709/reseq --filename trim
```

---

## 4. Reference Genome Preparation

### Download Reference

```bash
# Arabidopsis TAIR10 reference (example)
wget https://ftp.ensemblgenomes.org/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz --no-check-certificate
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa reference.fasta
```

### Create BWA-MEM2 Index

BWA-MEM2 is the faster successor to BWA-MEM:

```bash
bwa-mem2 index reference.fasta
```

### Create FASTA Index and Sequence Dictionary (required by GATK)

```bash
samtools faidx reference.fasta
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

The `@RG` (read group) tag is **required** by GATK:

```bash
bwa-mem2 mem \
  -t 8 \
  -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  reference.fasta \
  trim/wgs_R1_val_1.fq.gz \
  trim/wgs_R2_val_2.fq.gz \
  | samtools sort -@ 4 -o sample1.bam

samtools index sample1.bam
samtools flagstat sample1.bam
```

### Alignment Statistics

```bash
samtools stats sample1.bam > sample1.stats
multiqc .
```

| Metric | What to Expect |
|--------|---------------|
| Mapping rate | > 95% for matched reference |
| Properly paired | > 90% |
| Average depth | 10–30x for population WGS; 30–60x for clinical |

---

## 6. Mark Duplicates

PCR amplification creates duplicate reads. For variant calling, these **must be marked** (not just discarded — GATK handles them).

### Run Picard MarkDuplicates

```bash
picard MarkDuplicates \
  I=sample1.bam \
  O=sample1.markdup.bam \
  M=sample1.markdup.metrics \
  VALIDATION_STRINGENCY=SILENT

samtools index sample1.markdup.bam
```

### Check Duplication Rate

```bash
cat sample1.markdup.metrics
multiqc .
```

> **Note:** Duplication rates > 30–40% may indicate problems with input DNA quality or library complexity. For very high duplication, use a PCR-free library prep.

---

## 7. Base Quality Score Recalibration (BQSR)

GATK BQSR corrects systematic errors in base quality scores from the sequencer. It requires a **known variants VCF** (e.g., dbSNP).

### Download Known Variants

BQSR requires a VCF of known polymorphic sites to distinguish true variants from sequencing errors.

```bash
# For Arabidopsis: download from the 1001 Genomes Project
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
mv 1001genomes_snp-short-indel_only_ACGTN.vcf.gz known_variants.vcf.gz
bgzip -d known_variants.vcf.gz
gatk IndexFeatureFile -I known_variants.vcf

# For human (hg38): download dbSNP from GATK resource bundle
# wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
```

> **No known variants available?** For less-characterized species, skip BQSR or use an iterative "bootstrap" approach: call variants with HaplotypeCaller, use those as known sites, then re-run BQSR.

### Step 1: Compute Recalibration Table

```bash
gatk BaseRecalibrator \
  -I sample1.markdup.bam \
  -R reference.fasta \
  --known-sites known_variants.vcf \
  -O sample1.recal.table \
  --verbosity ERROR
```

### Step 2: Apply Recalibration

```bash
gatk ApplyBQSR \
  -I sample1.markdup.bam \
  -R reference.fasta \
  --bqsr-recal-file sample1.recal.table \
  -O sample1.recal.bam

samtools index sample1.recal.bam
```

> **Note for model organisms without known variants:** The BQSR "bootstrap" approach: (1) Call raw variants with HaplotypeCaller → hard-filter → use as known sites → re-run BQSR → call final variants. Two rounds are usually sufficient.

---

## 8. Variant Calling with GATK HaplotypeCaller

GATK HaplotypeCaller performs local assembly of haplotypes to call SNPs and indels simultaneously.

### Call Variants (per-sample GVCF mode)

Using **GVCF mode** allows joint genotyping across multiple samples — recommended for multi-sample projects:

```bash
gatk HaplotypeCaller \
  -R reference.fasta \
  -I sample1.recal.bam \
  -O sample1.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads 4
```

### Joint Genotyping (for multiple samples)

```bash
# Combine GVCFs (if multiple samples)
gatk CombineGVCFs \
  -R reference.fasta \
  -V sample1.g.vcf.gz \
  -V sample2.g.vcf.gz \
  -O cohort.g.vcf.gz

# Genotype
gatk GenotypeGVCFs \
  -R reference.fasta \
  -V cohort.g.vcf.gz \
  -O cohort.vcf.gz
```

---

## 9. Variant Filtering

### Separate SNPs and Indels

```bash
# Extract SNPs
gatk SelectVariants \
  -R reference.fasta \
  -V cohort.vcf.gz \
  --select-type-to-include SNP \
  -O cohort.snps.vcf.gz

# Extract Indels
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

| Filter | Field | Threshold (SNP) | Meaning |
|--------|-------|----------------|---------|
| QD | Quality by Depth | < 2.0 | Normalized variant quality |
| FS | Fisher Strand Bias | > 60.0 | Strand bias test |
| MQ | RMS Mapping Quality | < 40.0 | Average mapping quality |
| MQRankSum | MQ Rank Sum | < -12.5 | Comparison of MQ for ref/alt |
| ReadPosRankSum | Read Position RS | < -8.0 | Position bias within reads |

---

## 10. VCF Format

### VCF File Structure

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO              FORMAT     SAMPLE1
Chr1    12345   .       A    T    220   PASS    DP=45;AF=0.5;...  GT:AD:DP   0/1:22,23:45
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

```bash
# Count variants
bcftools stats cohort.snps.filtered.vcf.gz | grep "^SN"

# Filter PASS variants only
bcftools view -f PASS cohort.snps.filtered.vcf.gz -O z -o cohort.snps.pass.vcf.gz

# Extract specific region
bcftools view cohort.snps.pass.vcf.gz Chr1:100000-200000

# Get variant summary
bcftools stats cohort.snps.pass.vcf.gz > stats.txt
multiqc .
```

---

## 11. Variant Annotation with SnpEff

SnpEff predicts the functional effect of each variant (missense, nonsense, synonymous, etc.).

### Setup SnpEff Database

```bash
# List available databases
snpEff databases | grep -i arabidopsis

# Download database (example: Arabidopsis TAIR10)
snpEff download athalianaTair10
```

### Annotate Variants

```bash
snpEff \
  -v athalianaTair10 \
  cohort.snps.pass.vcf.gz \
  > cohort.snps.annotated.vcf

# View annotation summary
cat snpEff_summary.html
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

```bash
# Convert VCF to PLINK format
plink --vcf cohort.snps.pass.vcf.gz \
  --make-bed \
  --out cohort \
  --allow-extra-chr

# Run PCA
plink --bfile cohort \
  --pca 10 \
  --out cohort.pca \
  --allow-extra-chr
```

### Linkage Disequilibrium Pruning

Before running GWAS or population structure analyses, prune variants in LD to avoid biasing results:

```bash
plink --bfile cohort \
  --indep-pairwise 50 10 0.2 \
  --out cohort.ld \
  --allow-extra-chr

plink --bfile cohort \
  --extract cohort.ld.prune.in \
  --make-bed \
  --out cohort.pruned \
  --allow-extra-chr
```

### GWAS

For GWAS workflows, see the [GWAS tutorial]({{site.baseurl}}/episodes/GWAS/).

---

## 13. Full Workflow Summary

| Step | Tool | Input | Output |
|------|------|-------|--------|
| QC | FastQC + MultiQC | FASTQ | HTML report |
| Trim | Trim Galore | FASTQ | Trimmed FASTQ |
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

## References

| Resource | Link |
|---------|------|
| GATK Best Practices | [broadinstitute.github.io/gatk](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) |
| BWA-MEM2 paper | [Md et al. 2019, iScience](https://www.sciencedirect.com/science/article/pii/S2589004219310017) |
| SAM flag decoder | [Picard SAM Flags](https://broadinstitute.github.io/picard/explain-flags.html) |
| VCF specification | [samtools.github.io/hts-specs](https://samtools.github.io/hts-specs/VCFv4.2.pdf) |
| SnpEff documentation | [pcingola.github.io/SnpEff](https://pcingola.github.io/SnpEff/) |
| BCFtools manual | [samtools.github.io/bcftools](https://samtools.github.io/bcftools/bcftools.html) |
