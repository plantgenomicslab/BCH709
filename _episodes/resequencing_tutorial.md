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
pip install 'multiqc<1.34'

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

Expected output:

```
fastp 1.3.3
2.2.1
samtools 1.23.1
Using htslib 1.23.1
The Genome Analysis Toolkit (GATK) v4.6.2.0
bcftools 1.23.1
Using htslib 1.23.1
```

> **Common Error — `libcrypto.so.1.0.0: cannot open shared object file`**
>
> ```
> samtools: error while loading shared libraries: libcrypto.so.1.0.0:
> cannot open shared object file: No such file or directory
> ```
> **Cause:** OpenSSL version mismatch between samtools/bcftools and conda environment.
> **Fix:** Check which `libcrypto` file exists and create a symlink:
> ```bash
> ls ${CONDA_PREFIX}/lib/libcrypto.so.*
> # Use the actual version you see (often libcrypto.so.3 or libcrypto.so.1.1)
> ln -s ${CONDA_PREFIX}/lib/libcrypto.so.3 ${CONDA_PREFIX}/lib/libcrypto.so.1.0.0
> ```
{: .callout}

> **Common Error — Conda install hangs or fails with "Solving environment"**
>
> **Cause:** Default conda solver is slow and can get stuck resolving dependencies.
> **Fix:** Use mamba or the libmamba solver:
> ```bash
> conda install -n base -c conda-forge mamba
> # Then use `mamba install` instead of `conda install`
> ```
{: .callout}

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

> **Common Error — Dropbox download gives tiny HTML file instead of FASTQ**
>
> Symptom: downloaded file is only a few KB and running `file` on it shows `HTML document`.
> **Cause:** Missed the `&dl=1` at the end of the URL, or forgot to quote the URL.
> **Fix:**
> ```bash
> # WRONG — shell interprets &
> wget https://www.dropbox.com/.../file.fastq.gz?rlkey=...&dl=1
> # RIGHT — quote the URL
> wget "https://www.dropbox.com/.../file.fastq.gz?rlkey=...&dl=1" -O file.fastq.gz
> ```
{: .callout}

> **Common Error — `fasterq-dump: segmentation fault`**
>
> ```
> fasterq-dump: Segmentation fault (core dumped)
> ```
> **Cause:** Known SRA Toolkit bug on WSL (Windows Subsystem for Linux).
> **Fix:** Skip the SRA Toolkit — download FASTQ directly from ENA using `wget` as shown above.
{: .callout}

> You can find ENA download links for any SRA accession at [ENA Browser](https://www.ebi.ac.uk/ena/browser/) or [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra).

### Run FastQC

```bash
# -t 4 : use 4 threads
fastqc -t 4 sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz

# Aggregate all QC reports into one
multiqc . --exclude rsem --exclude gatk
```

Expected output (MultiQC summary):

```
/// MultiQC v1.34
       file_search | Search path: /path/to/reseq
            fastqc | Found 4 reports
     write_results | Data        : multiqc_data
     write_results | Report      : multiqc_report.html
           multiqc | MultiQC complete
```

After running, you should see four `*_fastqc.html`, four `*_fastqc.zip`, plus
`multiqc_report.html` and a `multiqc_data/` directory.

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

Expected output (fastp tail of stderr for sample 1):

```
Filtering result:
reads passed filter: 53051286
reads failed due to low quality: 4575222
reads failed due to too many N: 64736
reads failed due to too short: 1145000
reads failed due to adapter dimer: 0
reads with adapter trimmed: 3182403
bases trimmed due to adapters: 71080529

Duplication rate: 11.0013%

Insert size peak (evaluated by paired-end reads): 119

JSON report: trim/sample1_fastp_report.json
HTML report: trim/sample1_fastp_report.html

fastp v1.3.3, time used: 253 seconds
```

Numbers will vary slightly with each fastp version, but you should see >95% of
reads passing the filter and an insert-size peak in the 100-300 bp range for
typical paired-end Illumina libraries.

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

Expected output (`cat reference.fasta.fai`):

```
1	30427671	3	79	80
2	19698289	30812838	79	80
3	23459830	50760476	79	80
4	18585056	74517269	79	80
5	26975502	93337582	79	80
Pt	154478	120654551	79	80
Mt	367808	120810989	70	71
```

The reference has chromosomes `1`-`5`, plastid (`Pt`), and mitochondrion (`Mt`).

> **Common Error — `reference.dict already exists`**
>
> ```
> Exception in thread "main" picard.PicardException: /path/reference.dict already exists.
> Delete this file and try again, or specify a different output file.
> ```
> **Cause:** Running CreateSequenceDictionary a second time without removing the old `.dict`.
> **Fix:**
> ```bash
> rm -f reference.dict
> picard CreateSequenceDictionary R=reference.fasta O=reference.dict
> ```
{: .callout}

> **Common Error — BWA-MEM2 index files missing**
>
> ```
> ERROR! Can't open index file reference.fasta.0123
> ```
> **Cause:** Index not built, or index was interrupted (disk full, killed process).
> **Fix:** Delete any partial index files and re-run:
> ```bash
> rm -f reference.fasta.0123 reference.fasta.amb reference.fasta.ann \
>       reference.fasta.bwt.2bit.64 reference.fasta.pac
> bwa-mem2 index reference.fasta
> ```
{: .callout}

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

> **Common Error — `BAM file ... does not have any read groups`**
>
> ```
> A USER ERROR has occurred: SAM/BAM/CRAM file <file> does not have any read groups defined
> ```
> **Cause:** Forgot the `-R "@RG\t..."` flag in `bwa-mem2 mem`. GATK requires a read group.
> **Fix:** Re-run alignment with the `-R` flag. Use **literal `\t`** (not actual Tab) — the shell converts it correctly when wrapped in double quotes.
{: .callout}

> **Common Error — `samtools sort` fails with "No space left on device"**
>
> ```
> samtools sort: couldn't allocate memory for bam_mem
> [E::bgzf_flush] File write failed (wrong size)
> ```
> **Cause:** Temporary sort files fill the disk. BWA-MEM2 + sort can use 2-3x the final BAM size in temp space.
> **Fix:** Free up space first (`df -h .`, delete old BAMs), or redirect temp files to a larger disk:
> ```bash
> bwa-mem2 mem ... | samtools sort -T /path/to/larger/disk/tmp -@ 4 -o sample1.bam
> ```
{: .callout}

> **Common Error — Pipe fails silently, resulting BAM is truncated**
>
> Symptom: `samtools quickcheck sample1.bam` reports EOF missing, but no error was shown during alignment.
> **Cause:** `bwa-mem2` crashed inside the pipe; `samtools sort` still produces a (truncated) BAM.
> **Fix:** Add `set -o pipefail` before the command, or check exit codes:
> ```bash
> set -o pipefail
> bwa-mem2 mem ... | samtools sort -@ 4 -o sample1.bam
> echo "Pipeline exit: $?"  # Should be 0
> samtools quickcheck sample1.bam && echo OK
> ```
{: .callout}

### Alignment Statistics

```bash
samtools flagstat sample1.bam
samtools flagstat sample2.bam

samtools stats sample1.bam > sample1.stats
samtools stats sample2.bam > sample2.stats
multiqc . -n alignment_report --exclude rsem --exclude gatk
```

Expected output (`samtools flagstat sample1.bam`):

```
53076354 + 0 in total (QC-passed reads + QC-failed reads)
53051286 + 0 primary
0 + 0 secondary
25068 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
52354774 + 0 mapped (98.64% : N/A)
52329706 + 0 primary mapped (98.64% : N/A)
53051286 + 0 paired in sequencing
26525643 + 0 read1
26525643 + 0 read2
48411760 + 0 properly paired (91.25% : N/A)
52276894 + 0 with itself and mate mapped
52812 + 0 singletons (0.10% : N/A)
490848 + 0 with mate mapped to a different chr
129225 + 0 with mate mapped to a different chr (mapQ>=5)
```

Look for **mapped > 95%** and **properly paired > 90%** as a sanity check that
the reference matches your reads.

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

> **Common Error — `Exception in thread "main" java.lang.OutOfMemoryError: Java heap space`**
>
> ```
> Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
>     at picard.sam.markduplicates.MarkDuplicates.doWork(MarkDuplicates.java:...)
> ```
> **Cause:** Picard's default heap size (~2 GB) is too small for large BAMs.
> **Fix:** Increase Java heap memory:
> ```bash
> picard -Xmx8g MarkDuplicates I=sample1.bam O=sample1.markdup.bam M=sample1.markdup.metrics
> ```
{: .callout}

> **Common Error — MarkDuplicates fails with `SAM validation error`**
>
> ```
> SAM validation error: ERROR: Record ..., MAPQ should be 0 for unmapped read
> ```
> **Cause:** BWA-MEM2 occasionally emits non-standard fields that strict validators reject.
> **Fix:** We already use `VALIDATION_STRINGENCY=SILENT` to suppress these warnings. Ensure that option is in your command.
{: .callout}

### Check Duplication Rate

```bash
cat sample1.markdup.metrics
cat sample2.markdup.metrics
multiqc . -n markdup_report --exclude rsem --exclude gatk
```

Expected output (sample1.markdup.metrics — `LIBRARY` line, tab-delimited):

```
LIBRARY  UNPAIRED_READS_EXAMINED  READ_PAIRS_EXAMINED  SECONDARY_OR_SUPPLEMENTARY_RDS  UNMAPPED_READS  UNPAIRED_READ_DUPLICATES  READ_PAIR_DUPLICATES  READ_PAIR_OPTICAL_DUPLICATES  PERCENT_DUPLICATION  ESTIMATED_LIBRARY_SIZE
lib1     52812                    26138447             25068                            721580          16888                     7358725               0                              0.281567             37224722
```

Note `PERCENT_DUPLICATION = 0.281567` (~28%) — high but typical for older
publicly available 1001 Genomes samples. PCR-free libraries from modern
sequencing should be < 10%.

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

> **Common Error — `VCF file is malformed: there are N genotypes while the header requires M genotypes`**
>
> ```
> A USER ERROR has occurred: Error while trying to create index for known_variants.vcf.
> Error was: htsjdk.tribble.TribbleException: The provided VCF file is malformed at
> approximately line number 1202192: there are 460 genotypes while the header requires
> that 1135 genotypes be present for all records at 1:13410579
> ```
> **Cause:** The 1001 Genomes multi-sample VCF has inconsistent genotype columns across rows.
> **Fix:** Use the sites-only extraction shown above (strips genotype columns and reheaders).
{: .callout}

> **Common Error — `A USER ERROR has occurred: Bad input: ... Contig 'X' is not defined`**
>
> ```
> A USER ERROR has occurred: Bad input: The contig 'Chr1' is not defined in the
> sequence dictionary for reference reference.fasta
> ```
> **Cause:** Chromosome naming mismatch — the VCF has `Chr1` but the reference has `1`.
> **Fix:** Rename chromosomes in the VCF to match the reference:
> ```bash
> # Check names in both
> head reference.fasta.fai
> zcat known_sites.vcf.gz | grep -v '^#' | cut -f1 | sort -u
> # Rename with bcftools if needed
> echo "Chr1 1" > chr_map.txt
> bcftools annotate --rename-chrs chr_map.txt known_sites.vcf.gz -Oz -o known_sites.renamed.vcf.gz
> ```
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

> **Common Error — `ApplyBQSR` crashes with `No space left on device`**
>
> ```
> htsjdk.samtools.SAMException: Exception writing BAM index file
>     Caused by: java.io.IOException: No space left on device
> ```
> **Cause:** Disk filled up during BAM writing. BQSR doubles the size of the markdup BAM temporarily.
> **Fix:** Free up space, then re-run. Check with `df -h .` — you need at least 2x the markdup BAM size free.
> ```bash
> # Remove failed output and retry
> rm -f sample1.recal.bam sample1.recal.bai
> gatk ApplyBQSR -I sample1.markdup.bam -R reference.fasta \
>     --bqsr-recal-file sample1.recal.table -O sample1.recal.bam
> ```
{: .callout}

> **Common Error — `samtools quickcheck` reports EOF missing**
>
> ```
> sample1.recal.bam was missing EOF block when one should be present.
> ```
> **Cause:** BAM file is truncated — usually the writing process was killed or disk filled up.
> **Fix:** Always validate BAMs after writing; re-run if truncated:
> ```bash
> samtools quickcheck sample1.recal.bam && echo OK || echo TRUNCATED
> ```
{: .callout}

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

> **Common Error — `An index is required but was not found`**
>
> ```
> A USER ERROR has occurred: An index is required but was not found for file
> sample1.g.vcf.gz. Support for unindexed block-compressed files has been
> temporarily disabled. Try running IndexFeatureFile on the input.
> ```
> **Cause:** HaplotypeCaller was interrupted before writing the `.tbi` index file.
> **Fix:**
> ```bash
> gatk IndexFeatureFile -I sample1.g.vcf.gz
> ```
{: .callout}

> **Common Error — `Premature end of file` on GVCF**
>
> ```
> htsjdk.samtools.FileTruncatedException: Premature end of file: sample1.g.vcf.gz
> ```
> **Cause:** HaplotypeCaller was killed or crashed before finishing; the GVCF is truncated.
> **Fix:** Delete the truncated file and re-run HaplotypeCaller:
> ```bash
> rm -f sample1.g.vcf.gz sample1.g.vcf.gz.tbi
> gatk HaplotypeCaller -R reference.fasta -I sample1.recal.bam \
>     -O sample1.g.vcf.gz -ERC GVCF --native-pair-hmm-threads 4
> ```
> **Tip:** HaplotypeCaller is slow (30-60 min for full genome). To test the pipeline on a small region first, use `-L 1:1-5000000`.
{: .callout}

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

Expected output (last lines of `gatk GenotypeGVCFs`):

```
INFO  ProgressMeter -            Mt:313550              6.1              32351894        5274149.1
INFO  ProgressMeter - Traversal complete. Processed 32351894 total variants in 6.1 minutes.
INFO  GenotypeGVCFs - Shutting down engine
[org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs] done. Elapsed time: 6.16 minutes.
```

Run `bcftools stats cohort.vcf.gz | grep "^SN"` to count variants. For a
single-sample test on this dataset you should see ~225K SNPs and ~50K indels.

> **Single-sample alternative:** If you only have one sample, skip `CombineGVCFs` and run `GenotypeGVCFs` directly on the single GVCF:
> ```bash
> gatk GenotypeGVCFs -R reference.fasta -V sample1.g.vcf.gz -O cohort.vcf.gz
> ```
> Note: population-level analyses (PCA, LD pruning in Section 12) require 2+ samples.
{: .callout}

> **Common Error — `CombineGVCFs` fails with `Input files ... have incompatible contigs`**
>
> ```
> A USER ERROR has occurred: Input files sample1.g.vcf.gz and reference.fasta have
> incompatible contigs: No overlapping contigs found.
> ```
> **Cause:** The GVCF and reference were built against different assemblies.
> **Fix:** Ensure all GVCFs used the same reference used in the current step. Rebuild GVCFs if necessary.
{: .callout}

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

# Index the filtered VCF (required for region queries)
tabix -p vcf cohort.snps.pass.vcf.gz

# Extract specific genomic region
# Note: Ensembl TAIR10 uses "1", "2", ... (not "Chr1", "Chr2")
bcftools view cohort.snps.pass.vcf.gz 1:100000-200000

# Get variant summary statistics
bcftools stats cohort.snps.pass.vcf.gz > stats.txt
multiqc . -n vcf_report --exclude rsem --exclude gatk
```

Expected output (`bcftools stats cohort.snps.pass.vcf.gz | grep "^SN"`):

```
SN	0	number of samples:	1
SN	0	number of records:	205962
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	205962
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	253
SN	0	number of multiallelic SNP sites:	71
```

A region query like `bcftools view cohort.snps.pass.vcf.gz 1:100000-200000`
should return ~19 PASS SNPs in that 100-kb window.

> **Common Error — `Could not retrieve index file`**
>
> ```
> [E::idx_find_and_load] Could not retrieve index file for 'cohort.snps.pass.vcf.gz'
> Failed to read from cohort.snps.pass.vcf.gz: could not load index
> ```
> **Cause:** Region queries (e.g., `1:100000-200000`) require a `.tbi` tabix index.
> **Fix:**
> ```bash
> tabix -p vcf cohort.snps.pass.vcf.gz
> ```
{: .callout}

> **Common Error — Region query returns nothing**
>
> ```bash
> $ bcftools view cohort.snps.pass.vcf.gz Chr1:100000-200000
> # Only header shown, no variants
> ```
> **Cause:** Chromosome name mismatch — Ensembl TAIR10 uses `1`, not `Chr1`.
> **Fix:** Check names in the VCF and use the exact name:
> ```bash
> zcat cohort.snps.pass.vcf.gz | grep -v '^#' | cut -f1 | sort -u
> bcftools view cohort.snps.pass.vcf.gz 1:100000-200000
> ```
{: .callout}

---

## 11. Variant Annotation with SnpEff

SnpEff predicts the functional effect of each variant (missense, nonsense, synonymous, etc.).

### Setup SnpEff Database

```bash
# List available databases
snpEff databases | grep -i arabidopsis

# Download database
# For SnpEff 5.x, the database name is "Arabidopsis_thaliana"
# Use the output of the command above to confirm the exact name.
snpEff download Arabidopsis_thaliana
```

### Annotate Variants

```bash
snpEff Arabidopsis_thaliana cohort.snps.pass.vcf.gz > cohort.snps.annotated.vcf

# snpEff also generates snpEff_summary.html and snpEff_genes.txt
```

Expected output (files created):

```
cohort.snps.annotated.vcf   # ~166 MB for ~206K SNPs
snpEff_summary.html         # ~330 KB
snpEff_genes.txt            # ~1.9 MB
```

Open `snpEff_summary.html` in a browser to see variant-effect distribution
(MODIFIER, LOW, MODERATE, HIGH) per chromosome, transcript-biotype counts,
and Ts/Tv ratios.

> **Common Error — `ERROR: Genome 'athalianaTair10' not found`**
>
> ```
> java.lang.RuntimeException: ERROR: Genome 'athalianaTair10' not found.
>     Please check your data directory or config file.
> ```
> **Cause:** Database name differs between SnpEff versions. Older docs use `athalianaTair10`, but SnpEff 5.x uses `Arabidopsis_thaliana`.
> **Fix:** Check the exact name with:
> ```bash
> snpEff databases | grep -i "arabidopsis"
> ```
{: .callout}

> **Common Error — stderr mixed into annotated VCF**
>
> Symptom: `cohort.snps.annotated.vcf` starts with timestamps and log messages instead of `##fileformat=VCFv4.2`.
> **Cause:** Using `-v` (verbose) flag prints log messages to stderr, and if you redirect both streams (`&>`) they get mixed into the output.
> **Fix:** Use plain `>` to redirect only stdout, or redirect stderr separately:
> ```bash
> snpEff Arabidopsis_thaliana cohort.snps.pass.vcf.gz > cohort.snps.annotated.vcf 2> snpeff.log
> ```
{: .callout}

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

> **Common Error — `Error: Invalid chromosome code`**
>
> ```
> Error: Invalid chromosome code 'Mt' on line 1 of .bim file.
> (This is disallowed for human analyses; use --allow-extra-chr to permit it.)
> ```
> **Cause:** Arabidopsis uses chromosome names (`Mt`, `Pt`) that don't exist in PLINK's default human chromosome list.
> **Fix:** Add `--allow-extra-chr` to every PLINK command.
{: .callout}

> **Common Error — `At least 2 people required for pairwise analysis`**
>
> ```
> Error: At least 2 people required for pairwise analysis.
> ```
> **Cause:** PCA and LD analyses need multiple samples (founders). Single-sample data cannot be used for population analyses.
> **Fix:** Process at least 2 samples (e.g., sample1 + sample2 in this tutorial) before running PCA.
{: .callout}

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

> **Common Error — `Skipping --indep-pairwise since there are less than two founders`**
>
> ```
> Warning: Skipping --indep-pairwise since there are less than two founders.
> (--make-founders may come in handy here.)
> ```
> **Cause:** Same as PCA — requires at least 2 samples. With a single sample, LD pruning is meaningless.
> **Fix:** Add a second sample, or (for testing) declare the single sample as founder: `plink --bfile cohort --make-founders ...` — but the pruning output will be uninformative.
{: .callout}

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

## Troubleshooting — Quick Reference

Each section above has detailed error boxes. This table is a quick index; see the relevant section for the full fix.

| Problem | Section | Quick Fix |
|---------|---------|-----------|
| `No space left on device` | All | `df -h .`; delete intermediate files following `[Cleanup]` tips |
| `libcrypto.so.1.0.0: cannot open` | 1. Setup | `ln -s ${CONDA_PREFIX}/lib/libcrypto.so.3 ${CONDA_PREFIX}/lib/libcrypto.so.1.0.0` |
| Dropbox download gives HTML file | 2. QC | Quote URL and add `&dl=1` |
| `fasterq-dump: Segmentation fault` | 2. QC | Use ENA download via `wget` |
| `reference.dict already exists` | 4. Reference | `rm reference.dict` and re-run |
| BWA-MEM2 `Can't open index file` | 4–5. Index/Align | Delete partial index and re-run `bwa-mem2 index` |
| `BAM does not have read groups` | 5. Alignment | Add `-R "@RG\tID:...\tSM:..."` to `bwa-mem2 mem` |
| `samtools quickcheck` reports truncation | 5. Alignment, 7. BQSR | File write was interrupted — delete and re-run |
| `OutOfMemoryError: Java heap space` | 6. MarkDup, 8. HC | Increase heap: `picard -Xmx8g` / `gatk --java-options "-Xmx8g"` |
| VCF `malformed: N genotypes while header requires M` | 7. BQSR | Extract sites-only VCF (see Section 7) |
| `Contig 'X' is not defined` | 7. BQSR, 10. bcftools | Match chromosome names: `head reference.fasta.fai` vs VCF chromosomes |
| GVCF `Premature end of file` | 8. HC | Delete truncated GVCF; re-run HaplotypeCaller |
| `An index is required but was not found` | 8. HC | `gatk IndexFeatureFile -I sample1.g.vcf.gz` |
| `Could not retrieve index file` (bcftools) | 10. bcftools | `tabix -p vcf file.vcf.gz` |
| Region query returns no results | 10. bcftools | Check chromosome name: Ensembl TAIR10 uses `1`, not `Chr1` |
| SnpEff `Genome 'X' not found` | 11. SnpEff | `snpEff databases \| grep -i arabidopsis` to find correct name |
| PLINK `Invalid chromosome code 'Mt'` | 12. PLINK | Add `--allow-extra-chr` to every PLINK command |
| PLINK `At least 2 people required` | 12. PLINK | Need 2+ samples for PCA / LD analyses |

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
