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
    ├── fastp (trimming)
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
conda create -n chipseq -c bioconda -c conda-forge python=3.11 -y
conda activate chipseq

# Core alignment + QC tools (via conda)
conda install -c bioconda -c conda-forge fastqc fastp minimap2 samtools -y
conda install -c bioconda -c conda-forge openjdk=17 picard bedtools -y

# MACS3 and deepTools install cleanest via pip (conda has dependency conflicts)
pip install macs3 deeptools 'multiqc<1.34'

# HOMER (optional, for motif analysis — large download)
conda install -c bioconda homer -y

# IDR (optional, for replicate reproducibility analysis)
# NOTE: the pip package named "idr" is unrelated; install from GitHub:
pip install git+https://github.com/nboley/idr.git
```

> **Note on conda errors:** if you see `TypeError: 'NoneType' object is not iterable` from `conda install`, your conda solver is too old. Either upgrade conda (`conda update -n base conda`) or add `--solver=libmamba` to each install command.
{: .callout}

> **Note:** R packages (`ChIPseeker`, `DiffBind`) are installed separately within R (see Sections 13–14).
{: .callout}

> **libcrypto fix (if samtools reports missing library):**
> ```bash
> # Check which libcrypto is present first
> ls ${CONDA_PREFIX}/lib/libcrypto.so.*
> # Symlink it (use the actual version you see, often .so.3 on modern installs)
> ln -sf ${CONDA_PREFIX}/lib/libcrypto.so.3 ${CONDA_PREFIX}/lib/libcrypto.so.1.0.0
> ```
{: .callout}

### Verify Installations

```bash
fastp --version
minimap2 --version
samtools --version
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
| Input control | Match ChIP depth | Or >1.5x ChIP |

### Antibody Validation

Use only antibodies validated for ChIP:
- Check manufacturer's ChIP validation data
- Run western blot to confirm specificity
- Use ENCODE validated antibodies when possible

---

## 3. Quality Control

### Download Example Data

We use CTCF ChIP-Seq on human K562 cells (ENCODE experiment `ENCSR000DWE`) together with its matching input control (`ENCSR000DWA`). All three files are single-end 36 bp Illumina reads aligned to hg19.

```bash
cd ~/bch709/chipseq

# CTCF ChIP-Seq on human K562 (ENCODE ENCSR000DWE)
wget https://www.encodeproject.org/files/ENCFF001HTP/@@download/ENCFF001HTP.fastq.gz -O chip.fastq.gz

# Matching input control for K562 (ENCODE ENCSR000DWA)
wget https://www.encodeproject.org/files/ENCFF001HTT/@@download/ENCFF001HTT.fastq.gz -O input.fastq.gz

ls -lh *.fastq.gz   # expect ~910MB and ~690MB respectively
```

> **Verify the download:** ENCODE's S3 URLs are large (~1 GB each). If `wget` gets interrupted, the file will be truncated and later steps (fastp, minimap2) fail with `unexpected end of file`. Always validate:
> ```bash
> gunzip -t chip.fastq.gz && echo "chip OK"   # should print "chip OK"
> gunzip -t input.fastq.gz && echo "input OK"
> ```
> If the test fails, resume with `wget -c <url>` (the `-c` flag continues from where it stopped).
{: .callout}

> For ENCODE data, browse experiments at [encodeproject.org](https://www.encodeproject.org) and download FASTQ files directly. Each experiment page lists its "Controls" — always pair ChIP files with their matching input, not a random control.

### Run FastQC

```bash
fastqc -t 4 chip.fastq.gz input.fastq.gz
multiqc . --exclude rsem --exclude gatk
```

> Expected output:
> ```
> Started analysis of chip.fastq.gz
> Approx 5% complete for chip.fastq.gz
> ...
> Analysis complete for chip.fastq.gz
> Started analysis of input.fastq.gz
> ...
> Analysis complete for input.fastq.gz
> ```
> Produces `chip_fastqc.html`, `chip_fastqc.zip`, `input_fastqc.html`, `input_fastqc.zip`, then `multiqc_report.html`.
{: .solution}

### ChIP-Seq Specific QC Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| Mapping rate | % reads aligned | > 80% |
| Non-redundant fraction (NRF) | Unique reads / total | > 0.8 |
| PCR Bottleneck Coefficient (PBC) | M1/Md | > 0.9 |
| FRiP (Fraction of Reads in Peaks) | Reads in peaks / total | > 1% TF; > 0.3% histone |

---

## 4. Read Trimming

**fastp** performs adapter detection, quality trimming, and QC reporting in a single pass.

| Option | Description |
|--------|-------------|
| `--in1` | Input FASTQ (single-end) or forward reads (paired-end) |
| `--out1` | Output trimmed FASTQ |
| `--detect_adapter_for_pe` | Auto-detect adapters (use for paired-end) |
| `--qualified_quality_phred 20` | Minimum base quality threshold (Q20) |
| `--length_required 25` | Discard reads shorter than 25 bp |
| `--thread 4` | Number of threads |

```bash
mkdir -p trim

# Trim ChIP sample
fastp \
  --in1 chip.fastq.gz \
  --out1 trim/chip_trimmed.fq.gz \
  --qualified_quality_phred 20 \
  --length_required 25 \
  --thread 4 \
  --html trim/chip_fastp.html \
  --json trim/chip_fastp.json

# Trim Input control
fastp \
  --in1 input.fastq.gz \
  --out1 trim/input_trimmed.fq.gz \
  --qualified_quality_phred 20 \
  --length_required 25 \
  --thread 4 \
  --html trim/input_fastp.html \
  --json trim/input_fastp.json

multiqc trim/ -n trim_report
```

> For **paired-end** data, add `--in2`, `--out2`, and `--detect_adapter_for_pe` options.
{: .callout}

---

## 5. Reference Genome Preparation

### Download Reference

Our example FASTQs (ENCFF001HTP / ENCFF001HTT) are from early ENCODE and were aligned to **hg19**. For full-genome mapping we download the full assembly; for a quick tutorial run, chr19 alone (~58 MB) is enough.

```bash
# Option A — full hg19 genome (~3 GB after decompression, recommended)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
mv hg19.fa reference.fasta

# Option B — chr19 only, for a fast tutorial run (small disk / short time)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
gunzip chr19.fa.gz
mv chr19.fa reference.fasta
```

> **Using Option B (single chromosome)?** Only reads that map to chr19 will align (~2–6% of total reads). This is fine for testing pipeline commands but won't give biologically meaningful peak counts. For real analysis use Option A.
{: .callout}

### Create FASTA Index

```bash
samtools faidx reference.fasta
```

### Index the Reference for Minimap2

Minimap2 can build the index on the fly, but pre-computing it speeds up repeated alignments:

```bash
# -d : save index to file
minimap2 -d reference.mmi reference.fasta
```

> **Why Minimap2 for ChIP-Seq?**
> Minimap2 (`-ax sr` mode) handles short-read DNA alignment efficiently and is suitable for both short (< 100 bp) and longer Illumina reads. Unlike HISAT2/STAR, it is not splice-aware, making it appropriate for ChIP-Seq where reads come from genomic DNA.

---

## 6. Alignment with Minimap2

Minimap2 uses the `-ax sr` preset for short Illumina reads.

| Option | Description |
|--------|-------------|
| `-ax sr` | Short-read alignment preset |
| `-t 4` | Number of threads |
| `reference.mmi` | Pre-built minimap2 index |

### Align ChIP Sample

```bash
set -o pipefail   # so any failure in the pipe is caught

minimap2 -ax sr -t 4 reference.mmi trim/chip_trimmed.fq.gz \
  2> chip.minimap2.log \
  | samtools sort -@ 4 -o chip.bam

samtools index chip.bam
samtools quickcheck chip.bam && echo "chip BAM OK"
samtools flagstat chip.bam
```

> Expected output (counts vary with depth/library; mapping rate near 55–60% is typical for these old 36 bp ENCODE files vs full hg19):
> ```
> chip BAM OK
> 19956973 + 0 in total (QC-passed reads + QC-failed reads)
> 19956935 + 0 primary
> 38 + 0 supplementary
> 10928167 + 0 mapped (54.76% : N/A)
> ```
{: .solution}

### Align Input Control

```bash
minimap2 -ax sr -t 4 reference.mmi trim/input_trimmed.fq.gz \
  2> input.minimap2.log \
  | samtools sort -@ 4 -o input.bam

samtools index input.bam
```

> **Paired-end ChIP-Seq reads?** Supply both FASTQ files:
> ```bash
> minimap2 -ax sr -t 4 reference.mmi \
>   trim/chip_R1_trimmed.fq.gz trim/chip_R2_trimmed.fq.gz \
>   2> chip.minimap2.log | samtools sort -@ 4 -o chip.bam
> ```

### Filter Low-Quality and Unmapped Reads

| Option | Description |
|--------|-------------|
| `-b` | Output BAM format |
| `-q 30` | Minimum mapping quality (MAPQ >= 30) |
| `-F 4` | Exclude unmapped reads |

```bash
samtools view -b -q 30 -F 4 chip.bam > chip.filt.bam
samtools index chip.filt.bam

samtools view -b -q 30 -F 4 input.bam > input.filt.bam
samtools index input.filt.bam
```

---

## 7. Mark and Remove Duplicates

For ChIP-Seq, **PCR duplicates must be removed** (unlike variant calling where they are only marked).

| Option | Description |
|--------|-------------|
| `I=` | Input BAM |
| `O=` | Output BAM |
| `M=` | Duplication metrics file |
| `REMOVE_DUPLICATES=true` | Remove (not just flag) duplicate reads |
| `VALIDATION_STRINGENCY=SILENT` | Suppress warnings on BAM format |

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

| Option | Description |
|--------|-------------|
| `--bam` | Input BAM file |
| `--outFileName` | Output BigWig file |
| `--normalizeUsing RPKM` | Normalize by reads per kilobase per million |
| `--binSize 10` | Resolution in base pairs |
| `--numberOfProcessors 4` | Number of threads |

```bash
bamCoverage \
  --bam chip.dedup.bam \
  --outFileName chip.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --numberOfProcessors 4

bamCoverage \
  --bam input.dedup.bam \
  --outFileName input.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --numberOfProcessors 4
```

### ChIP vs. Input Ratio (log2 fold change)

```bash
bamCompare \
  --bamfile1 chip.dedup.bam \
  --bamfile2 input.dedup.bam \
  --outFileName chip_vs_input.bw \
  --scaleFactorsMethod None \
  --normalizeUsing RPKM \
  --operation log2 \
  --binSize 10 \
  --numberOfProcessors 4
```

> **Common Error — `--normalizeUsing RPKM` is only valid if you also use `--scaleFactorsMethod None`!**
>
> **Cause:** `bamCompare` applies scaling in two places by default; combining both with RPKM is a conflict.
> **Fix:** Always pass `--scaleFactorsMethod None` when using `--normalizeUsing RPKM` (shown above).
{: .callout}

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
  --numberOfProcessors 4
```

> **Interpreting the fingerprint plot:** A successful ChIP enrichment shows a steep curve (reads concentrated at a few genomic loci). The input control should be near-diagonal (reads evenly distributed). Poor enrichment = ChIP curve resembles the input.

### Replicate Correlation

> ⚠️ **Demonstration only — this tutorial uses a single ChIP replicate.**
> The canonical block below shows what you'd run with two true biological replicates (`chip_rep1.dedup.bam`, `chip_rep2.dedup.bam`). **Do not run it as-is** on this dataset — those files do not exist. A runnable pseudo-replicate version follows.
{: .callout}

**Canonical command (two true biological replicates — read-only example):**

```bash
# Bin the genome and count reads in each bin
multiBamSummary bins \
  --bamfiles chip_rep1.dedup.bam chip_rep2.dedup.bam input.dedup.bam \
  --labels Rep1 Rep2 Input \
  --outFileName multibam.npz \
  --numberOfProcessors 4

# Plot correlation heatmap
plotCorrelation \
  --corData multibam.npz \
  --corMethod pearson \
  --skipZeros \
  --whatToPlot heatmap \
  --colorMap RdYlBu \
  --plotFile correlation.png \
  --outFileCorMatrix correlation.txt
```

**Runnable on this single-replicate dataset (pseudo-replicates):**

To still see `multiBamSummary` and `plotCorrelation` execute end-to-end, split `chip.dedup.bam` into two halves and treat them as pseudo-replicates. This is **not** a substitute for true biological replicates — it only demonstrates the commands.

```bash
# Split chip.dedup.bam into two halves by alternating reads
samtools view -h chip.dedup.bam \
  | awk 'BEGIN{n=0} /^@/ {print > "h1.sam"; print > "h2.sam"; next}
                       {if (n%2==0) print > "h1.sam"; else print > "h2.sam"; n++}'
samtools sort -@ 4 -o chip_pseudo1.dedup.bam h1.sam && samtools index chip_pseudo1.dedup.bam
samtools sort -@ 4 -o chip_pseudo2.dedup.bam h2.sam && samtools index chip_pseudo2.dedup.bam
rm -f h1.sam h2.sam

# Bin the genome and count reads in each bin
multiBamSummary bins \
  --bamfiles chip_pseudo1.dedup.bam chip_pseudo2.dedup.bam input.dedup.bam \
  --labels Pseudo1 Pseudo2 Input \
  --outFileName multibam.npz \
  --numberOfProcessors 4

# Plot correlation heatmap
plotCorrelation \
  --corData multibam.npz \
  --corMethod pearson \
  --skipZeros \
  --whatToPlot heatmap \
  --colorMap RdYlBu \
  --plotFile correlation.png \
  --outFileCorMatrix correlation.txt
```

> Expected output (tab-separated `correlation.txt`, pseudo-replicate split — Pseudo1↔Pseudo2 Pearson is artificially close to 1.0 because both halves come from the same library):
> ```
> #plotCorrelation --outFileCorMatrix
>             'Input'   'Pseudo1'  'Pseudo2'
> 'Input'      1.0000   -0.2314    -0.2298
> 'Pseudo1'   -0.2314    1.0000     0.9987
> 'Pseudo2'   -0.2298    0.9987     1.0000
> ```
> Also produces `correlation.png` heatmap. For real biological replicates with deeper coverage, expect Rep1↔Rep2 Pearson > 0.9 and ChIP↔Input near 0.
{: .solution}

> **Common Error — `error: the following arguments are required: --whatToPlot/-p`**
>
> Older documentation often says `--plotType` but current `plotCorrelation` uses `--whatToPlot` (valid values: `heatmap`, `scatterplot`).
{: .callout}

**Target:** Pearson correlation between biological replicates > 0.9

---

## 10. Peak Calling with MACS3

MACS3 (Model-based Analysis of ChIP-Seq) is the standard tool for identifying enriched genomic regions.

### MACS3 Key Options

| Option | Description |
|--------|-------------|
| `-t` | Treatment (ChIP) BAM file |
| `-c` | Control (Input) BAM file |
| `--format BAM` | Input file format |
| `--gsize hs` | Effective genome size (`hs`=human, `mm`=mouse, `ce`=*C. elegans*, `dm`=*Drosophila*) |
| `--name` | Output file prefix |
| `--outdir` | Output directory |
| `--qvalue 0.05` | FDR threshold for peak calling |
| `--broad` | Call broad peaks (for histone marks) |
| `--broad-cutoff 0.1` | FDR threshold for broad peaks |

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
  --nomodel --extsize 147 \
  2>&1 | tee macs3_narrow.log
```

> Expected output (tail of log; counts will vary):
> ```
> #1 tag size = 36.0
> #1 total tags in treatment: 10604
> #1 tags after filtering in treatment: 10597
> #2 Skipped...
> #2 Use 147 as fragment length
> #4 Write peak in narrowPeak format file... macs3_narrow/chip_narrow_peaks.narrowPeak
> Done!
> ```
> Result: `macs3_narrow/chip_narrow_peaks.narrowPeak` with 9170 peaks for the example data.
{: .solution}

> **Why `--nomodel --extsize 147`?** The default model fails on the 36 bp ENCODE example data with "Total number of paired peaks: 0" (see callout below). The fix is shipped here in the main command so the tutorial runs end-to-end.
{: .callout}

> **Common Error — `MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only find 0`**
>
> **Cause:** MACS3 tries to build a fragment-size model from strand-paired peaks. With low-depth data (or a single chromosome / small reference), it can't find enough peaks to model.
> **Fix:** Bypass the model and set a fixed fragment size (147 bp is a reasonable default for nucleosome-associated reads):
> ```bash
> macs3 callpeak -t chip.dedup.bam -c input.dedup.bam --format BAM --gsize hs \
>     --name chip_narrow --outdir macs3_narrow --qvalue 0.05 \
>     --nomodel --extsize 147
> ```
{: .callout}

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
  --nomodel --extsize 147 \
  2>&1 | tee macs3_broad.log
```

> Expected output: `macs3_broad/chip_broad_peaks.broadPeak` (9162 broad peaks for the example data) plus `chip_broad_peaks.gappedPeak`.
{: .solution}

### MACS3 Output Files

| File | Description |
|------|-------------|
| `_peaks.narrowPeak` | Peak coordinates + statistics |
| `_peaks.xls` | Detailed peak table |
| `_summits.bed` | Peak summit positions |
| `_model.r` | Fragment size model (run in R) |
| `_control_lambda.bdg` | Control lambda values |
| `_treat_pileup.bdg` | ChIP pileup signal |

### NarrowPeak Format

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Start (0-based) |
| 3 | End |
| 4 | Peak name |
| 5 | Score |
| 6 | Strand |
| 7 | Signal value (fold enrichment) |
| 8 | -log10(p-value) |
| 9 | -log10(q-value) |
| 10 | Summit offset from start |

### Count and Filter Peaks

```bash
# Count total peaks
wc -l macs3_narrow/chip_narrow_peaks.narrowPeak

# View top peaks by fold enrichment (column 7)
sort -k7,7rn macs3_narrow/chip_narrow_peaks.narrowPeak | head -20

# Filter by FDR q-value < 0.01 (column 9 = -log10 qvalue, so > 2)
awk '$9 > 2' macs3_narrow/chip_narrow_peaks.narrowPeak > chip_narrow_filtered.bed
```

> Expected output:
> ```
> 9170 macs3_narrow/chip_narrow_peaks.narrowPeak
>
> # Top peaks by fold enrichment:
> chrX  17963258  17963516  chip_narrow_peak_8759  171  .  6.99626  26.6055  17.1159  129
> chr1  225668620 225668878  chip_narrow_peak_763   145  .  5.99679  22.4886  14.5672  129
> chr21  47645076  47645336  chip_narrow_peak_4926  145  .  5.99679  22.4886  14.5672  130
>
> # After q<0.01 filter:
> 228 chip_narrow_filtered.bed
> ```
{: .solution}

---

## 11. Signal Profiles and Heatmaps

Visualize ChIP signal around features (e.g., TSS, peak centers).

### Download Gene Annotation and Build TSS BED

```bash
# Download refGene table (UCSC format: tab-separated text)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip refGene.txt.gz

# Build a BED file of transcription start sites (TSS = txStart on +, txEnd on -)
awk 'BEGIN{OFS="\t"} {
    if ($4=="+") print $3, $5, $5+1, $2, "0", $4;
    else         print $3, $6-1, $6, $2, "0", $4
}' refGene.txt | sort -k1,1 -k2,2n | uniq > TSS.bed

head TSS.bed
wc -l TSS.bed
```

> Expected output:
> ```
> chr1   11868  11869  NR_148357  0  +
> chr1   11873  11874  NR_046018  0  +
> chr1   17435  17436  NR_106918  0  -
> chr1   17435  17436  NR_107062  0  -
> chr1   17435  17436  NR_107063  0  -
> ...
> 81407 TSS.bed
> ```
{: .solution}

The resulting `TSS.bed` has one row per transcript, with coordinates pointing to the exact TSS base.

### Compute Signal Matrix Around TSS

| Option | Description |
|--------|-------------|
| `--scoreFileName` | BigWig signal file |
| `--regionsFileName` | BED file with genomic regions |
| `--referencePoint TSS` | Anchor point for the plot |
| `--upstream` / `--downstream` | Window size around reference point |
| `--binSize 50` | Resolution in base pairs |

```bash
computeMatrix reference-point \
  --scoreFileName chip.bw \
  --regionsFileName TSS.bed \
  --referencePoint TSS \
  --upstream 3000 \
  --downstream 3000 \
  --binSize 50 \
  --numberOfProcessors 4 \
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

> **Match the genome build to your alignment.** This tutorial aligned reads to **hg19** (Section 5), so peak coordinates are in hg19. Install and search against **hg19** — using `hg38` here would extract sequences at the wrong positions and yield meaningless motifs.
{: .callout}

```bash
# Install HOMER genome (run once) — must match the build used for alignment
configureHomer.pl -install hg19
```

**`findMotifsGenome.pl` usage:** `findMotifsGenome.pl <peaks.bed> <genome> <output_dir/> [options]`

| Option | Description |
|--------|-------------|
| `-size 200` | Region size around peak center for motif search |
| `-mask` | Mask repeat sequences |
| `-p 4` | Number of threads |

```bash
# Prepare peak file (convert narrowPeak to BED6)
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
  chip_narrow_filtered.bed > peaks_homer.bed

# Find motifs (de novo + known)
findMotifsGenome.pl peaks_homer.bed hg19 motif_output/ -size 200 -mask -p 4
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

### Install DiffBind

```r
BiocManager::install("DiffBind")
```

### Sample Sheet Format

The sample sheet CSV must contain these columns:

| SampleID | Condition | Replicate | bamReads | Peaks | PeakCaller |
|----------|-----------|-----------|----------|-------|------------|
| chip_cond1_rep1 | Condition1 | 1 | chip_cond1_rep1.dedup.bam | peaks_cond1_rep1.narrowPeak | macs |
| chip_cond1_rep2 | Condition1 | 2 | chip_cond1_rep2.dedup.bam | peaks_cond1_rep2.narrowPeak | macs |
| chip_cond2_rep1 | Condition2 | 1 | chip_cond2_rep1.dedup.bam | peaks_cond2_rep1.narrowPeak | macs |
| chip_cond2_rep2 | Condition2 | 2 | chip_cond2_rep2.dedup.bam | peaks_cond2_rep2.narrowPeak | macs |

### Run DiffBind

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

> ⚠️ **Demonstration only — IDR requires two true biological replicates.**
> This tutorial called peaks once (`macs3_narrow/chip_narrow_peaks.narrowPeak`) on a single-replicate dataset, so `chip_rep1_peaks.narrowPeak` / `chip_rep2_peaks.narrowPeak` **do not exist**. The block below is the canonical command you would run if you had two replicate-level peak calls — **do not run it as-is** on this dataset. Unlike correlation analysis, IDR cannot be meaningfully approximated with pseudo-replicates from one library: its statistical model fundamentally needs two independent biological replicates.
{: .callout}

| Option | Description |
|--------|-------------|
| `--samples` | Two replicate peak files |
| `--input-file-type narrowPeak` | Input file format |
| `--rank p.value` | Column to rank peaks by |
| `--output-file` | Output IDR results |
| `--plot` | Generate IDR diagnostic plot |

**Canonical command (two true biological replicates — read-only example):**

```bash
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

> **To run IDR for real:** start from two independent biological-replicate ChIP samples (e.g., download a second ENCODE replicate paired with the same input), repeat Sections 4–10 to produce `chip_rep1_peaks.narrowPeak` and `chip_rep2_peaks.narrowPeak`, then run the block above. ENCODE's IDR pipeline is documented at <https://github.com/ENCODE-DCC/chip-seq-pipeline2>.
{: .callout}

> **Common Error — `AttributeError: module 'numpy' has no attribute 'int'`**
>
> IDR 2.0.3 is not compatible with NumPy ≥ 1.20 (where `np.int` was removed).
> **Fix:** Downgrade NumPy in the chipseq env, or use the bioconda build (which patches this):
> ```bash
> pip install "numpy<1.20"
> # OR use the bioconda-packaged idr:
> conda install -c bioconda idr
> ```
{: .callout}

> **Common Error — `idr: Image Data Repository access library` (wrong package)**
>
> `pip install idr` pulls a completely unrelated package (for microscopy images). Always install from the Kundaje/Boley repo:
> ```bash
> pip uninstall idr -y
> pip install git+https://github.com/nboley/idr.git
> ```
{: .callout}

---

## 16. Full Workflow Summary

| Step | Tool | Input | Output |
|------|------|-------|--------|
| QC | FastQC + MultiQC | FASTQ | HTML report |
| Trim | fastp | FASTQ | Trimmed FASTQ |
| Align | Minimap2 | FASTQ | BAM |
| Sort/Index | SAMtools | BAM | Sorted BAM |
| Filter | SAMtools | BAM | Filtered BAM (MAPQ >= 30) |
| Dedup | Picard | BAM | Deduplicated BAM |
| Signal tracks | deepTools bamCoverage | BAM | BigWig |
| IP quality | deepTools plotFingerprint | BAM | Plot |
| Replicate QC | deepTools plotCorrelation | BAM | Correlation matrix |
| Peak calling | MACS3 | BAM | BED/narrowPeak |
| Heatmaps | deepTools computeMatrix/plotHeatmap | BigWig + BED | Heatmap PNG |
| Motif analysis | HOMER | BED | Motif report HTML |
| Annotation | ChIPseeker (R) | BED | Annotated peaks |
| Differential | DiffBind (R) | BAM + peaks | Differential peaks |
| Reproducibility | IDR | Peak files | IDR peaks |

---

## 17. Troubleshooting — Quick Reference

| Problem | Section | Quick Fix |
|---------|---------|-----------|
| `TypeError: 'NoneType' object is not iterable` from conda | 1. Setup | Old conda solver — `conda update -n base conda` or add `--solver=libmamba` |
| `libcrypto.so.1.0.0: cannot open` | 1. Setup | `ln -sf ${CONDA_PREFIX}/lib/libcrypto.so.3 ${CONDA_PREFIX}/lib/libcrypto.so.1.0.0` |
| `fastp: igzip: unexpected eof` | 3. QC | FASTQ download was truncated — check with `gunzip -t`, resume with `wget -c` |
| `samtools quickcheck` EOF missing | 6. Align | BAM write was interrupted — delete and re-run |
| MACS3 `Total number of paired peaks: 0` | 10. MACS3 | Low depth — add `--nomodel --extsize 147` |
| `--normalizeUsing RPKM is only valid if you also use --scaleFactorsMethod None` | 8. bamCompare | Add `--scaleFactorsMethod None` |
| `plotCorrelation: error: --whatToPlot required` | 9. deepTools | Use `--whatToPlot heatmap` (not `--plotType`) |
| IDR `module 'numpy' has no attribute 'int'` | 15. IDR | `pip install "numpy<1.20"` or use bioconda idr |
| IDR `Image Data Repository access library` installed | 1. Setup | Wrong package — install from `git+https://github.com/nboley/idr.git` |
| HOMER `findMotifsGenome.pl: command not found` | 12. HOMER | `conda install -c bioconda homer` and `configureHomer.pl -install hg19` |

---

## 18. Cleanup

```bash
# Remove large intermediate files when done
rm -f chip.bam chip.bam.bai chip.filt.bam chip.filt.bam.bai
rm -f input.bam input.bam.bai input.filt.bam input.filt.bam.bai
rm -rf trim/

conda deactivate
# Optional: remove environment entirely
conda env remove --name chipseq
```

---

## References

| Resource | Link |
|---------|------|
| ENCODE ChIP-Seq standards | [Landt et al. 2012, Genome Research](https://genome.cshlp.org/content/22/9/1813.long) |
| MACS3 paper | [Zhang et al. 2008, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) |
| deepTools paper | [Ramirez et al. 2016, Nucleic Acids Research](https://academic.oup.com/nar/article/44/W1/W160/2499308) |
| fastp paper | [Chen et al. 2018, Bioinformatics](https://doi.org/10.1093/bioinformatics/bty560) |
| ChIPseeker paper | [Yu et al. 2015, Bioinformatics](https://academic.oup.com/bioinformatics/article/31/14/2382/255379) |
| HOMER documentation | [homer.ucsd.edu](http://homer.ucsd.edu/homer/chipseq/) |
| IDR framework | [Li et al. 2011, Ann. Applied Statistics](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Measuring-reproducibility-of-high-throughput-experiments/10.1214/11-AOAS466.full) |
| DiffBind paper | [Ross-Innes et al. 2012, Nature](https://www.nature.com/articles/nature10730) |
