---
layout: page
title: ChIP-Seq Tutorial
published: true
---

> ## Paper Reading
> Please read these papers before class:
> - [Kharchenko et al. (2008) Design and analysis of ChIP-seq experiments for DNA-binding proteins. *Nature Biotechnology* 26:1351–1359](https://www.nature.com/articles/nbt.1508)
> - [Landt et al. (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. *Genome Research* 22:1813–1831](https://genome.cshlp.org/content/22/9/1813.long)
{: .prereq}

---

## Overview: ChIP-Sequencing

**ChIP-Seq** (Chromatin Immunoprecipitation followed by Sequencing) identifies genomic regions bound by a protein of interest (transcription factor, histone modification) across the entire genome.

### How ChIP-Seq Works

1. **Crosslink** — Fix protein-DNA interactions with formaldehyde
2. **Shear** — Fragment chromatin by sonication or MNase digestion
3. **Immunoprecipitate** — Pull down target protein with a specific antibody
4. **Reverse crosslinks** — Release DNA from protein
5. **Sequence** — Sequence the enriched DNA fragments

### Types of ChIP-Seq Targets

| Target Type | Examples | Peak Width | Notes |
|------------|---------|-----------|-------|
| Transcription factor (TF) | CTCF, p53, MYC | Narrow (100–500 bp) | Sharp peaks |
| Active promoter | H3K4me3, H3K9ac | Narrow | Near TSS |
| Active enhancer | H3K4me1, H3K27ac | Broad (1–10 kb) | Distal regulatory |
| Repressive | H3K27me3, H3K9me3 | Very broad (10–100 kb) | Heterochromatin |
| Elongating RNAPII | H3K36me3 | Gene body | Transcribed regions |

### ChIP-Seq Controls

| Control Type | Description | When to Use |
|-------------|-------------|-------------|
| Input | Chromatin before IP | Always recommended |
| IgG | Non-specific antibody IP | Alternative to input |
| Mock IP | No antibody | Less common |

---

## The ChIP-Seq Workflow

```
Raw FASTQ (ChIP + Input/Control)
    │
    ├── FastQC → MultiQC (QC)
    ├── Trim Galore (trimming)
    │
    ├── Minimap2 (alignment to reference genome)
    │       │
    │       └── SAMtools sort + index
    │
    ├── Picard MarkDuplicates (remove PCR duplicates)
    │
    ├── deepTools bamCoverage (generate bigWig tracks)
    │
    ├── MACS3 (peak calling)
    │       │
    │       ├── Narrow peaks (TF, H3K4me3)
    │       └── Broad peaks (H3K27me3, H3K36me3)
    │
    ├── deepTools plotFingerprint (IP quality)
    ├── deepTools plotCorrelation (replicate concordance)
    ├── deepTools computeMatrix + plotHeatmap (signal profiles)
    │
    ├── HOMER / MEME-ChIP (motif analysis)
    │
    └── ChIPseeker / DiffBind (peak annotation, differential binding)
```

---

## 1. Environment Setup

### Create Conda Environment

```bash
conda create -n chipseq python=3.11
conda activate chipseq
```

### Install Tools

```bash
conda install -c bioconda -c conda-forge \
  fastqc trim-galore minimap2 samtools picard deeptools \
  macs3 homer bedtools multiqc r-base bioconductor-chipseeker
```

### Verify Installations

```bash
minimap2 --version
macs3 --version
deeptools --version
```

### Create Working Directory

```bash
mkdir -p ~/bch709/chipseq
cd ~/bch709/chipseq
```

---

## 2. Experimental Design Considerations

### Replicates

| Type | Minimum | Recommended |
|------|---------|------------|
| Biological replicates | 2 | 3+ |
| Technical replicates | Not needed | — |

**Replicate concordance targets (ENCODE):**
- IDR (Irreproducibility Discovery Rate) < 0.05
- Pearson correlation of read counts: > 0.9 (same antibody)

### Sequencing Depth

| Target | Minimum Reads | Notes |
|--------|--------------|-------|
| TF ChIP-Seq | 20 million | Narrow peaks |
| Histone ChIP-Seq | 40–50 million | Broad marks |
| Input control | Match ChIP depth | Or >1.5× ChIP |

### Antibody Validation

Use only antibodies validated for ChIP:
- Check manufacturer's ChIP validation data
- Run western blot to confirm specificity
- Use ENCODE validated antibodies when possible

---

## 3. Quality Control

### Download Example Data

```bash
cd ~/bch709/chipseq

# Example: CTCF ChIP-Seq (human K562, from ENCODE)
# ChIP replicate 1
wget https://www.encodeproject.org/files/ENCFF001NQP/@@download/ENCFF001NQP.fastq.gz -O chip_R1.fastq.gz
# Input control
wget https://www.encodeproject.org/files/ENCFF001NQQ/@@download/ENCFF001NQQ.fastq.gz -O input.fastq.gz
```

> For ENCODE data, browse experiments at [encodeproject.org](https://www.encodeproject.org) and download FASTQ files directly.

### Run FastQC

```bash
fastqc -t 4 chip_R1.fastq.gz input.fastq.gz
multiqc .
```

### ChIP-Seq Specific QC Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| Mapping rate | % reads aligned | > 80% |
| Non-redundant fraction (NRF) | Unique reads / total | > 0.8 |
| PCR Bottleneck Coefficient (PBC) | M1/Md | > 0.9 |
| FRiP (Fraction of Reads in Peaks) | Reads in peaks / total | > 1% TF; > 0.3% histone |

---

## 4. Read Trimming

```bash
trim_galore \
  --cores 4 \
  --fastqc \
  --gzip \
  -o trim \
  chip_R1.fastq.gz

# Also trim the input control
trim_galore \
  --cores 4 \
  --fastqc \
  --gzip \
  -o trim \
  input.fastq.gz

multiqc --dirs ~/bch709/chipseq --filename trim
```

---

## 5. Reference Genome Preparation

### Download Reference

```bash
# Human hg38 (example: chr1 only for tutorial)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
gunzip chr1.fa.gz
mv chr1.fa reference.fasta

# Or download the full genome
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

### Index the Reference for Minimap2

Minimap2 does not require a pre-built index step (it builds the index on the fly), but you can pre-compute and save it to speed up repeated alignments:

```bash
minimap2 -d reference.mmi reference.fasta
```

> **Why Minimap2 for ChIP-Seq?**
> Minimap2 (`-ax sr` mode) handles short-read DNA alignment efficiently and is suitable for both short (< 100 bp) and longer Illumina reads. Unlike HISAT2/STAR, it is not splice-aware, making it appropriate for ChIP-Seq where reads come from genomic DNA. It is also widely used for long-read (PacBio/ONT) data.

---

## 6. Alignment with Minimap2

Minimap2 uses the `-ax sr` preset for short paired-end Illumina reads. Output is piped directly to `samtools sort` to avoid intermediate SAM files.

### Align ChIP Sample

```bash
minimap2 \
  -ax sr \
  -t 8 \
  reference.mmi \
  trim/chip_R1_trimmed.fq.gz \
  2> chip.minimap2.log \
  | samtools sort -@ 4 -o chip.bam

samtools index chip.bam
samtools flagstat chip.bam
```

### Align Input Control

```bash
minimap2 \
  -ax sr \
  -t 8 \
  reference.mmi \
  trim/input_trimmed.fq.gz \
  2> input.minimap2.log \
  | samtools sort -@ 4 -o input.bam

samtools index input.bam
```

> **Paired-end ChIP-Seq reads?** Use both FASTQ files:
> ```bash
> minimap2 -ax sr -t 8 reference.mmi \
>   trim/chip_R1_trimmed.fq.gz trim/chip_R2_trimmed.fq.gz \
>   2> chip.minimap2.log | samtools sort -@ 4 -o chip.bam
> ```

### Filter Low-Quality and Unmapped Reads

```bash
# Keep only properly mapped, high-quality reads (MAPQ >= 30)
samtools view -b -q 30 -F 4 chip.bam > chip.filt.bam
samtools index chip.filt.bam

samtools view -b -q 30 -F 4 input.bam > input.filt.bam
samtools index input.filt.bam
```

---

## 7. Mark and Remove Duplicates

For ChIP-Seq, **PCR duplicates must be removed** (unlike RNA-Seq):

```bash
picard MarkDuplicates \
  I=chip.filt.bam \
  O=chip.dedup.bam \
  M=chip.dedup.metrics \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=SILENT

samtools index chip.dedup.bam

picard MarkDuplicates \
  I=input.filt.bam \
  O=input.dedup.bam \
  M=input.dedup.metrics \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=SILENT

samtools index input.dedup.bam
```

---

## 8. Generate BigWig Signal Tracks

BigWig files are used for visualizing ChIP-Seq signal in genome browsers.

### Normalize by Sequencing Depth (RPKM)

```bash
bamCoverage \
  --bam chip.dedup.bam \
  --outFileName chip.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --numberOfProcessors 8

bamCoverage \
  --bam input.dedup.bam \
  --outFileName input.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --numberOfProcessors 8
```

### ChIP vs. Input Ratio (log2 fold change)

```bash
bamCompare \
  --bamfile1 chip.dedup.bam \
  --bamfile2 input.dedup.bam \
  --outFileName chip_vs_input.bw \
  --normalizeUsing RPKM \
  --operation log2 \
  --binSize 10 \
  --numberOfProcessors 8
```

### Visualize in IGV

Load `chip.bw`, `input.bw`, and `chip_vs_input.bw` into [IGV](https://software.broadinstitute.org/software/igv/) for visual inspection.

---

## 9. IP Quality Assessment with deepTools

### Fingerprint Plot

The fingerprint plot shows how well the ChIP enrichment worked. A perfect ChIP shows a steep curve; input should be diagonal.

```bash
plotFingerprint \
  --bamfiles chip.dedup.bam input.dedup.bam \
  --labels ChIP Input \
  --plotFile fingerprint.png \
  --outRawCounts fingerprint.tab \
  --numberOfProcessors 8
```

> **Interpreting the fingerprint plot:** A successful ChIP enrichment shows a steep curve (reads concentrated at a few genomic loci). The input control should be near-diagonal (reads evenly distributed). Poor enrichment = ChIP curve resembles the input.

### Replicate Correlation

```bash
# Bin the genome and count reads in each bin
multiBamSummary bins \
  --bamfiles chip_rep1.dedup.bam chip_rep2.dedup.bam input.dedup.bam \
  --labels Rep1 Rep2 Input \
  --outFileName multibam.npz \
  --numberOfProcessors 8

# Plot correlation heatmap
plotCorrelation \
  --corData multibam.npz \
  --corMethod pearson \
  --skipZeros \
  --plotType heatmap \
  --colorMap RdYlBu \
  --plotFile correlation.png \
  --outFileCorMatrix correlation.txt
```

**Target:** Pearson correlation between biological replicates > 0.9

---

## 10. Peak Calling with MACS3

MACS3 (Model-based Analysis of ChIP-Seq) is the standard tool for identifying genomic regions enriched in ChIP-Seq.

### Narrow Peak Calling (TF, H3K4me3, H3K9ac)

```bash
macs3 callpeak \
  -t chip.dedup.bam \
  -c input.dedup.bam \
  --format BAM \
  --gsize hs \
  --name chip_narrow \
  --outdir macs3_narrow \
  --qvalue 0.05 \
  2>&1 | tee macs3_narrow.log
```

### Broad Peak Calling (H3K27me3, H3K36me3, H3K9me3)

```bash
macs3 callpeak \
  -t chip.dedup.bam \
  -c input.dedup.bam \
  --format BAM \
  --gsize hs \
  --name chip_broad \
  --outdir macs3_broad \
  --broad \
  --broad-cutoff 0.1 \
  2>&1 | tee macs3_broad.log
```

### MACS3 Output Files

| File | Description |
|------|-------------|
| `_peaks.narrowPeak` | Peak coordinates + statistics |
| `_peaks.xls` | Detailed peak table |
| `_summits.bed` | Peak summit positions |
| `_model.r` | Fragment size model (run in R) |
| `_control_lambda.bdg` | Control lambda values |
| `_treat_pileup.bdg` | ChIP pileup signal |

### Interpret NarrowPeak Format

```
Col 1: Chromosome
Col 2: Start (0-based)
Col 3: End
Col 4: Peak name
Col 5: Score
Col 6: Strand
Col 7: Signal value (fold enrichment)
Col 8: -log10(p-value)
Col 9: -log10(q-value)
Col 10: Summit offset from start
```

### Count and Filter Peaks

```bash
# Count total peaks
wc -l macs3_narrow/chip_narrow_peaks.narrowPeak

# View top peaks by fold enrichment (column 7)
sort -k7,7rn macs3_narrow/chip_narrow_peaks.narrowPeak | head -20

# Filter by FDR q-value < 0.01 (column 9 = -log10 qvalue, so > 2)
awk '$9 > 2' macs3_narrow/chip_narrow_peaks.narrowPeak > chip_narrow_filtered.bed
```

---

## 11. Signal Profiles and Heatmaps

Visualize ChIP signal around features (e.g., TSS, peak centers).

### Download Gene Annotation

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
# Or use a BED/GTF file with TSS coordinates
```

### Compute Signal Matrix Around TSS

```bash
computeMatrix reference-point \
  --scoreFileName chip.bw \
  --regionsFileName TSS.bed \
  --referencePoint TSS \
  --upstream 3000 \
  --downstream 3000 \
  --binSize 50 \
  --numberOfProcessors 8 \
  -o tss_matrix.gz

plotHeatmap \
  --matrixFile tss_matrix.gz \
  --outFileName tss_heatmap.png \
  --colorMap Blues \
  --whatToShow 'heatmap and colorbar' \
  --zMin 0 --zMax 5
```

### Profile Plot

```bash
plotProfile \
  --matrixFile tss_matrix.gz \
  --outFileName tss_profile.png \
  --plotTitle "ChIP Signal at TSS"
```

---

## 12. Motif Analysis

Motif analysis identifies DNA sequence motifs enriched in ChIP-Seq peaks — typically the binding motif of the immunoprecipitated protein.

### HOMER Motif Analysis

HOMER requires genome sequences to be installed before running motif analysis:

```bash
# Install HOMER genome (run once)
perl /path/to/homer/configureHomer.pl -install hg38
# Or using the conda-installed HOMER:
configureHomer.pl -install hg38
```

```bash
# Prepare peak file (convert narrowPeak to BED)
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
  chip_narrow_filtered.bed > peaks_homer.bed

# Find motifs (de novo + known)
findMotifsGenome.pl \
  peaks_homer.bed \
  hg38 \
  motif_output/ \
  -size 200 \
  -mask \
  -p 8
```

### HOMER Output

| File | Description |
|------|-------------|
| `homerResults.html` | De novo discovered motifs |
| `knownResults.html` | Matches to known motifs database |
| `homerResults/motif1.motif` | PWM for top motif |

---

## 13. Peak Annotation with ChIPseeker (R)

ChIPseeker annotates peaks relative to genomic features (promoter, exon, intron, etc.).

### Install ChIPseeker

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "clusterProfiler"))
```

### Annotate Peaks

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peaks
peaks <- readPeakFile("chip_narrow_filtered.bed")

# Annotate
peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")

# Plot annotation pie chart
plotAnnoPie(peakAnno)

# Plot distance to TSS
plotDistToTSS(peakAnno, title="Distribution of ChIP peaks\nrelative to TSS")

# Gene ontology of nearby genes
genes <- seq2gene(peaks, tssRegion=c(-1000, 1000),
                   flankDistance=3000, TxDb=txdb)
pathway <- enrichPathway(genes)
dotplot(pathway)
```

---

## 14. Differential Binding Analysis with DiffBind (R)

Compare ChIP-Seq signal between two conditions (e.g., treated vs. untreated).

The sample sheet CSV must contain these columns:

| SampleID | Condition | Replicate | bamReads | Peaks | PeakCaller |
|----------|-----------|-----------|----------|-------|------------|
| chip_cond1_rep1 | Condition1 | 1 | chip_cond1_rep1.dedup.bam | peaks_cond1_rep1.narrowPeak | macs |
| chip_cond1_rep2 | Condition1 | 2 | chip_cond1_rep2.dedup.bam | peaks_cond1_rep2.narrowPeak | macs |
| chip_cond2_rep1 | Condition2 | 1 | chip_cond2_rep1.dedup.bam | peaks_cond2_rep1.narrowPeak | macs |
| chip_cond2_rep2 | Condition2 | 2 | chip_cond2_rep2.dedup.bam | peaks_cond2_rep2.narrowPeak | macs |

```r
library(DiffBind)

# Load sample sheet
samples <- read.csv("samplesheet.csv")
dba_obj <- dba(sampleSheet=samples)

# Count reads in peaks
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps=TRUE)

# Set contrast
dba_obj <- dba.contrast(dba_obj, categories=DBA_CONDITION)

# Run differential analysis (DESeq2 by default)
dba_obj <- dba.analyze(dba_obj)

# Get significant peaks
db_peaks <- dba.report(dba_obj, th=0.05)
```

---

## 15. IDR Analysis (Irreproducibility Discovery Rate)

IDR measures reproducibility of peak calls across replicates. ENCODE requires IDR < 0.05.

```bash
# Install IDR via conda (preferred)
conda install -c bioconda idr

# Sort peaks by score (column 5)
sort -k5,5rn chip_rep1_peaks.narrowPeak > rep1_sorted.bed
sort -k5,5rn chip_rep2_peaks.narrowPeak > rep2_sorted.bed

# Run IDR
idr \
  --samples rep1_sorted.bed rep2_sorted.bed \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file idr_peaks.txt \
  --plot \
  --log-output-file idr.log
```

---

## 16. Full Workflow Summary

| Step | Tool | Input | Output |
|------|------|-------|--------|
| QC | FastQC + MultiQC | FASTQ | HTML report |
| Trim | Trim Galore | FASTQ | Trimmed FASTQ |
| Align | Minimap2 | FASTQ | BAM |
| Sort/Index | SAMtools | BAM | Sorted BAM |
| Filter | SAMtools | BAM | Filtered BAM |
| Dedup | Picard | BAM | Deduplicated BAM |
| Signal tracks | deepTools bamCoverage | BAM | BigWig |
| IP quality | deepTools plotFingerprint | BAM | Plot |
| Replicate QC | deepTools plotCorrelation | BAM | Correlation matrix |
| Peak calling | MACS3 | BAM | BED/narrowPeak |
| Heatmaps | deepTools computeMatrix/plotHeatmap | BigWig + BED | Heatmap PNG |
| Motif analysis | HOMER | BED | Motif report HTML |
| Annotation | ChIPseeker | BED | Annotated peaks |
| Differential | DiffBind | BAM + peaks | Differential peaks |
| Reproducibility | IDR | Peak files | IDR peaks |

---

## 17. Cleanup

```bash
conda deactivate
conda env remove --name chipseq
```

---

## References

| Resource | Link |
|---------|------|
| ENCODE ChIP-Seq standards | [Landt et al. 2012, Genome Research](https://genome.cshlp.org/content/22/9/1813.long) |
| MACS3 paper | [Zhang et al. 2008, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) |
| deepTools paper | [Ramírez et al. 2016, Nucleic Acids Research](https://academic.oup.com/nar/article/44/W1/W160/2499308) |
| ChIPseeker paper | [Yu et al. 2015, Bioinformatics](https://academic.oup.com/bioinformatics/article/31/14/2382/255379) |
| HOMER documentation | [homer.ucsd.edu](http://homer.ucsd.edu/homer/chipseq/) |
| IDR framework | [Li et al. 2011, Ann. Applied Statistics](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Measuring-reproducibility-of-high-throughput-experiments/10.1214/11-AOAS466.full) |
| DiffBind paper | [Ross-Innes et al. 2012, Nature](https://www.nature.com/articles/nature10730) |
