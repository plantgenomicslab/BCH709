---
layout: page
title: Hi-C Tutorial — 3D Genome Organization
published: true
---

> ## Paper Reading
> Please read these papers before class:
> - [Lieberman-Aiden et al. (2009) Comprehensive mapping of long-range interactions reveals folding principles of the human genome. *Science* 326:289–293](https://www.science.org/doi/10.1126/science.1181369)
> - [Rao et al. (2014) A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell* 159:1665–1680](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4)
> - [Dixon et al. (2012) Topological domains in mammalian genomes identified by analysis of chromatin interactions. *Nature* 485:376–380](https://www.nature.com/articles/nature11082)
{: .prereq}

---

## Overview: Hi-C and 3D Genome Organization

**Hi-C** (High-throughput Chromosome Conformation Capture) is a genome-wide method for mapping chromatin contacts — which genomic regions are physically close to each other in the nucleus, even if they are far apart on the linear chromosome.

### The Central Dogma of 3D Genomics

Linear DNA is hierarchically folded in the nucleus:

```
DNA double helix (2 nm)
    └── Nucleosome fiber (10 nm)
            └── Chromatin fiber (30 nm)
                    └── Loops (kb–Mb scale)
                            └── TADs / Compartments (100 kb–Mb scale)
                                    └── Chromosomal territories
```

### Key 3D Genome Features

| Feature | Size | Description |
|---------|------|-------------|
| Loops | 10–300 kb | Direct chromatin contacts, often anchored by CTCF |
| TADs (Topologically Associating Domains) | 100 kb–2 Mb | Insulated domains of preferential interaction |
| Compartments | 1–10 Mb | A (active/euchromatin) and B (inactive/heterochromatin) |
| Chromosomal territories | Chromosome scale | Distinct nuclear positions per chromosome |

### How Hi-C Works

1. **Crosslink** — Fix chromatin contacts with formaldehyde
2. **Digest** — Cut DNA with a restriction enzyme (e.g., HindIII, MboI)
3. **Ligation** — Fill ends with biotinylated nucleotides and ligate nearby fragments
4. **Shear** — Fragment ligated DNA
5. **Pull down** — Enrich biotinylated ligation junctions with streptavidin beads
6. **Sequence** — Paired-end sequencing; each read pair represents a chromatin contact

### Hi-C Read Pairs

Each paired-end read represents a **chimeric ligation junction**: one read from each of two genomic loci that were spatially close in the nucleus.

---

## The Hi-C Analysis Workflow

```
Raw paired-end FASTQ (Hi-C reads)
    │
    ├── FastQC → MultiQC (QC)
    │
    ├── HiC-Pro (full pipeline):
    │       ├── Trim and align (Bowtie2, used internally by HiC-Pro)
    │       ├── Parse valid ligation junctions
    │       ├── Remove duplicates
    │       ├── Build contact matrices (multiple resolutions)
    │       └── ICE normalization
    │
    ├── Convert to .cool / .hic format
    │
    ├── Visualization (HiGlass, Juicebox, cooltools)
    │
    ├── A/B Compartment analysis (PCA on contact matrix)
    │
    ├── TAD calling (Insulation Score or Arrowhead)
    │
    └── Loop calling (HiCCUPS / MUSTACHE)
```

---

## 1. Environment Setup

### Create Conda Environment

```bash
conda create -n hic python=3.9
conda activate hic
```

### Install Tools

```bash
# Core Hi-C processing
# HiC-Pro uses Bowtie2 internally — install both together
conda install -c bioconda -c conda-forge \
  hic-pro bowtie2 samtools fastqc trim-galore multiqc

# Analysis and visualization tools
conda install -c bioconda -c conda-forge \
  cooler cooltools hicexplorer pyBigWig

# Optional: Juicer tools for .hic format and HiCCUPS loop calling (JAR-based)
wget https://github.com/aidenlab/juicer/releases/latest/download/juicer_tools.jar
```

### Verify Installations

```bash
HiC-Pro --version
cooler --help
```

### Create Working Directory

```bash
mkdir -p ~/bch709/hic
cd ~/bch709/hic
```

---

## 2. Experimental Design Considerations

### Library Preparation Methods

| Method | Enzyme | Fragment Size | Notes |
|--------|--------|--------------|-------|
| Hi-C (original) | HindIII | ~500 bp | Standard |
| in situ Hi-C | MboI / DpnII | ~100–300 bp | Higher resolution, standard now |
| Micro-C | MNase | Nucleosome-level | Single nucleosome resolution |
| SPRITE | No ligation | — | Multi-way contacts |
| HiChIP | MboI + ChIP | — | Protein-centric contacts |

### Sequencing Depth

| Resolution | Recommended Depth | Notes |
|-----------|-----------------|-------|
| 1 Mb compartments | 50–100 M read pairs | Shallow |
| 100 kb TADs | 200–500 M read pairs | Standard |
| 10 kb loops | 1–3 B read pairs | Deep Hi-C |
| Micro-C loops | > 2 B read pairs | Very deep |

### Biological Replicates

- Minimum 2 replicates
- Correlation of contact matrices between replicates should be > 0.95 (genome-wide insulation)

---

## 3. Quality Control

### Download Example Data

```bash
cd ~/bch709/hic

# Example: in situ Hi-C of human IMR90 fibroblasts (subsampled)
# (Replace with actual SRA accession for your data)
fasterq-dump --split-files SRR1658693 -O .
mv SRR1658693_1.fastq hic_R1.fastq
mv SRR1658693_2.fastq hic_R2.fastq
gzip hic_R1.fastq hic_R2.fastq
```

### Run FastQC

```bash
fastqc -t 4 hic_R1.fastq.gz hic_R2.fastq.gz
multiqc . --exclude rsem --exclude gatk
```

### Hi-C Specific QC Metrics

| Metric | Good Range | Description |
|--------|-----------|-------------|
| Valid pairs | > 30–40% of total | Properly ligated read pairs |
| Cis/Trans ratio | Cis > 70–80% | Most contacts should be within chromosomes |
| Short-range vs. long-range | Varies | Long-range (> 20 kb) enrichment indicates good data |
| Duplication rate | < 30% | Higher may indicate low complexity |

---

## 4. Hi-C Data Processing with HiC-Pro

HiC-Pro is a comprehensive pipeline for Hi-C data processing, from raw reads to contact matrices.

### Prepare Reference Genome

```bash
# Download hg38 (or your organism's genome)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build Bowtie2 index (HiC-Pro uses Bowtie2 internally)
bowtie2-build --threads 8 hg38.fa hg38

# Download chromosome size file (required for contact matrix)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

### Generate Restriction Fragment File

HiC-Pro needs to know where restriction enzyme sites are in the reference:

```bash
# For MboI (GATC sites) — in situ Hi-C
digest_genome.py \
  -r GATC \
  -o hg38_MboI.bed \
  hg38.fa

# For HindIII (AAGCTT sites) — original Hi-C
digest_genome.py \
  -r AAGCTT \
  -o hg38_HindIII.bed \
  hg38.fa
```

### Configure HiC-Pro

Create a configuration file `hicpro_config.txt`:

```ini
# Genome
GENOME_SIZE = /path/to/hg38.chrom.sizes
BOWTIE2_IDX_PATH = /path/to/bowtie2_index/
REFERENCE_GENOME = hg38

# Reads
PAIR1_EXT = _R1
PAIR2_EXT = _R2

# Trimming (optional)
MIN_MAPQ = 10

# Restriction enzyme
GENOME_FRAGMENT = /path/to/hg38_MboI.bed
LIGATION_SITE = GATCGATC
MIN_FRAG_SIZE = 100
MAX_FRAG_SIZE = 100000
MIN_INSERT_SIZE = 100
MAX_INSERT_SIZE = 600

# Contact matrix
BIN_SIZE = 1000000 500000 250000 100000 40000

# Output
N_CPU = 8
SORT_RAM = 768000
```

### Run HiC-Pro

```bash
# Organize input files
mkdir -p rawdata/sample1
cp hic_R1.fastq.gz hic_R2.fastq.gz rawdata/sample1/

# Run HiC-Pro
HiC-Pro \
  -i rawdata \
  -o hicpro_output \
  -c hicpro_config.txt
```

### HiC-Pro Output

```
hicpro_output/
├── logs/                     # Log files for each step
├── stat/                     # QC statistics
│   └── sample1/
│       ├── sample1_allValidPairs.mergestat
│       └── *.stat
└── hic_results/
    ├── data/sample1/         # Valid pairs file
    └── matrix/sample1/       # Contact matrices
        ├── raw/              # Raw counts at each resolution
        └── iced/             # ICE-normalized matrices
```

### Check HiC-Pro Statistics

```bash
cat hicpro_output/stat/sample1/*mergestat

# Key metrics to check:
# - Total pairs
# - Valid pairs (ligation junctions)
# - Cis pairs (same chromosome)
# - Trans pairs (different chromosomes)
# - Intra/inter-chromosomal split
```

---

## 5. Convert to .cool Format

`.cool` is a modern, efficient format for Hi-C contact matrices (HDF5-based).

```bash
# Convert HiC-Pro output to .cool
hicpro2higlass.sh \
  -i hicpro_output/hic_results/data/sample1/sample1.allValidPairs \
  -r 40000 \
  -c /path/to/hg38.chrom.sizes \
  -p 8 \
  -o sample1.cool

# Or convert from HiC-Pro matrix format
cooler load \
  --format hicpro \
  /path/to/hg38.chrom.sizes:40000 \
  hicpro_output/hic_results/matrix/sample1/raw/40000/sample1_40000.matrix \
  sample1_40000.cool
```

### Create Multi-Resolution .mcool

```bash
cooler zoomify \
  --nproc 8 \
  --resolutions 1000,5000,10000,25000,50000,100000,250000,500000,1000000 \
  sample1_40000.cool \
  -o sample1.mcool
```

---

## 6. Visualization

### HiGlass (Web-based, interactive)

[HiGlass](https://higlass.io) is a web tool for interactive Hi-C visualization. Load your `.mcool` file at [app.higlass.io](https://app.higlass.io).

### Juicebox

[Juicebox](https://www.aidenlab.org/juicebox/) is the most commonly used desktop viewer. It requires `.hic` format.

#### Convert to .hic format

```bash
# Using Juicer tools
java -jar juicer_tools.jar pre \
  hicpro_output/hic_results/data/sample1/sample1.allValidPairs \
  sample1.hic \
  hg38
```

### Quick Plot with cooltools

```python
# Python visualization with cooltools
import cooler
import cooltools
import numpy as np
import matplotlib.pyplot as plt

clr = cooler.Cooler("sample1.mcool::resolutions/100000")

# Plot a region
region = ("chr1", 0, 50_000_000)
mat = clr.matrix(balance=True).fetch(region)

plt.figure(figsize=(8, 8))
plt.imshow(np.log1p(mat), cmap="YlOrRd", vmax=2)
plt.colorbar(label="log1p(balanced contacts)")
plt.title("Hi-C contact map: chr1 0–50 Mb (100 kb resolution)")
plt.savefig("hic_map.png", dpi=150, bbox_inches="tight")
```

---

## 7. A/B Compartment Analysis

**A compartments** = transcriptionally active, euchromatin
**B compartments** = transcriptionally inactive, heterochromatin

A/B compartments are identified by PCA on the Hi-C contact matrix (or on the correlation matrix of contact patterns).

### Compute Compartments with cooltools

```bash
# Compute eigenvectors at 100 kb resolution
cooltools eigs-cis \
  sample1.mcool::resolutions/100000 \
  --phasing-track H3K27ac_chip.bw \
  --out-prefix sample1_compartments \
  --n-eigs 3

# Output:
# sample1_compartments.cis.vecs.tsv  — eigenvectors per bin
# sample1_compartments.cis.lam.txt   — eigenvalues
```

### Interpret Compartments

```bash
# The sign of the first eigenvector (E1) defines A/B:
# Positive E1 → A compartment (active)
# Negative E1 → B compartment (inactive)
# Sign must be oriented by a mark (H3K27ac, GC content, gene density)

# Summarize compartment genome coverage
awk '$6 > 0 {sum += $3-$2} END {print "A:", sum}' sample1_compartments.cis.vecs.tsv
awk '$6 < 0 {sum += $3-$2} END {print "B:", sum}' sample1_compartments.cis.vecs.tsv
```

---

## 8. TAD Calling

**TADs (Topologically Associating Domains)** are regions of preferential self-interaction, separated by insulating boundaries.

### Method 1: Insulation Score (cooltools)

The insulation score measures local self-interaction; dips in the score correspond to TAD boundaries.

```bash
cooltools insulation \
  sample1.mcool::resolutions/40000 \
  --window-bp 500000 \
  --out sample1_insulation.tsv

# Output columns: chr, start, end, ..., insulation_score, is_boundary
```

### Method 2: HiCExplorer TAD Calling

```bash
# Find TADs using the insulation score approach
hicFindTADs \
  --matrix sample1.mcool::resolutions/40000 \
  --outPrefix sample1_TADs \
  --correctForMultipleTesting fdr \
  --threshold 0.05 \
  --minDepth 120000 \
  --maxDepth 2000000 \
  --step 40000

# Output:
# sample1_TADs_domains.bed       — TAD coordinates
# sample1_TADs_boundaries.bed    — Boundary coordinates
# sample1_TADs_score.bw          — TAD boundary score (bigWig)
```

### Visualize TADs

```bash
# Plot contact matrix with TAD boundaries
hicPlotMatrix \
  --matrix sample1.mcool::resolutions/40000 \
  --outFileName hic_chr1_TADs.png \
  --region chr1:10000000-50000000 \
  --TAD sample1_TADs_domains.bed \
  --colorMap YlOrRd \
  --title "Hi-C chr1 with TADs" \
  --dpi 300
```

### TAD Statistics

```bash
# Count TADs
wc -l sample1_TADs_domains.bed

# TAD size distribution
awk '{print $3-$2}' sample1_TADs_domains.bed | sort -n | \
  awk 'BEGIN{sum=0; n=0} {sum+=$1; n++} END{print "Mean TAD size:", sum/n, "bp"}'
```

---

## 9. Loop Calling

Chromatin loops are specific point-to-point contacts, often anchored by CTCF/cohesin.

### HiCCUPS (Juicer Tools)

```bash
java -Xmx16g -jar juicer_tools.jar hiccups \
  --norms KR \
  --cpu \
  -r 5000,10000,25000 \
  sample1.hic \
  sample1_loops/

# Output:
# merged_loops.bedpe — Paired genomic coordinates of loops
```

### MUSTACHE (Alternative Loop Caller)

```bash
python -m mustache \
  -f sample1.mcool \
  -r 10000 \
  -o sample1_mustache_loops.tsv \
  -p 8
```

### Analyze Loops

```bash
# Count loops
wc -l sample1_loops/merged_loops.bedpe

# Loop size distribution
awk '{print $5-$2}' sample1_loops/merged_loops.bedpe | sort -n | \
  awk 'BEGIN{sum=0; n=0} {sum+=$1; n++} END{print "Mean loop size:", sum/n, "bp"}'

# Loops near CTCF peaks?
bedtools pairtopair \
  -a sample1_loops/merged_loops.bedpe \
  -b ctcf_peaks.bed \
  -type both
```

---

## 10. Differential Hi-C Analysis

Compare 3D genome organization between two conditions (e.g., WT vs. KO, treated vs. untreated).

### Using multiHiCcompare (R)

```r
library(multiHiCcompare)

# Load contact matrices
# hic1 <- read_table("condition1.txt")
# hic2 <- read_table("condition2.txt")

# Create HICEXP object
hicexp <- make_hicexp(
  data_list = list(cond1_rep1, cond1_rep2, cond2_rep1, cond2_rep2),
  groups = c(1, 1, 2, 2)
)

# Normalize
hicexp <- fastlo(hicexp)

# Test for differential interactions
hicexp <- hic_exactTest(hicexp)

# Get significant differential interactions
results <- results(hicexp, alpha=0.05, logFC_cutoff=1)
```

---

## 11. Integrating Hi-C with Other Data

### Correlate TADs with Gene Expression (RNA-Seq)

```bash
# Find genes within each TAD
bedtools intersect \
  -a sample1_TADs_domains.bed \
  -b genes.bed \
  -wa -wb > genes_in_TADs.txt
```

### Correlate Compartments with ChIP-Seq

```python
import cooltools
import pyBigWig

# Load compartment calls
import pandas as pd
comp = pd.read_table("sample1_compartments.cis.vecs.tsv")

# Load H3K27ac signal
bw = pyBigWig.open("H3K27ac.bw")
comp["H3K27ac"] = comp.apply(
    lambda row: bw.stats(row["chrom"], row["start"], row["end"], type="mean")[0],
    axis=1
)

# Correlate eigenvector with H3K27ac
from scipy import stats
r, p = stats.pearsonr(comp["E1"].dropna(), comp["H3K27ac"].dropna())
print(f"Pearson r = {r:.3f}, p = {p:.2e}")
```

### Overlap Loops with CTCF ChIP-Seq

```bash
# Extract loop anchors
awk '{print $1"\t"$2"\t"$3}' sample1_loops/merged_loops.bedpe > anchors_left.bed
awk '{print $4"\t"$5"\t"$6}' sample1_loops/merged_loops.bedpe > anchors_right.bed
cat anchors_left.bed anchors_right.bed | sort -k1,1 -k2,2n > all_anchors.bed

# Overlap with CTCF peaks
bedtools intersect \
  -a all_anchors.bed \
  -b ctcf_peaks.narrowPeak \
  -u > anchors_with_ctcf.bed

echo "Anchors total: $(wc -l < all_anchors.bed)"
echo "Anchors with CTCF: $(wc -l < anchors_with_ctcf.bed)"
```

---

## 12. Full Workflow Summary

| Step | Tool | Input | Output |
|------|------|-------|--------|
| QC | FastQC + MultiQC | FASTQ | HTML report |
| Full pipeline | HiC-Pro | FASTQ | Valid pairs, matrices |
| Format conversion | cooler | HiC-Pro output | .cool / .mcool |
| Format conversion | Juicer tools | Valid pairs | .hic |
| Visualization | HiGlass / Juicebox | .mcool / .hic | Interactive contact map |
| Compartments | cooltools eigs-cis | .cool | A/B compartment calls |
| TADs | HiCExplorer / cooltools | .cool | TAD boundaries |
| Loops | HiCCUPS / MUSTACHE | .hic / .cool | Loop anchors |
| Differential | multiHiCcompare | Count matrices | Differential interactions |
| Integration | BEDtools | BED + ChIP/RNA | Overlap tables |

---

## 13. Cleanup

```bash
conda deactivate
conda env remove --name hic
```

---

## Glossary

| Term | Definition |
|------|-----------|
| Contact matrix | Matrix where entry (i,j) = number of read pairs linking bins i and j |
| ICE normalization | Iterative Correction and Eigenvector decomposition — balances matrix biases |
| TAD | Topologically Associating Domain — genomic region of preferential self-interaction |
| Loop | Direct chromatin contact between two genomic loci (point-to-point) |
| Compartment A | Active, gene-rich, open chromatin |
| Compartment B | Inactive, gene-poor, heterochromatin |
| Insulation score | Local measure of cross-boundary interaction — dips mark TAD boundaries |
| CTCF | CCCTC-binding factor — key insulator/loop anchor protein |
| Cohesin | Protein complex that extrudes chromatin loops |

---

## References

| Resource | Link |
|---------|------|
| Lieberman-Aiden et al. (2009) Original Hi-C paper | [Science](https://www.science.org/doi/10.1126/science.1181369) |
| Rao et al. (2014) In situ Hi-C / loop calling | [Cell](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4) |
| Dixon et al. (2012) TADs | [Nature](https://www.nature.com/articles/nature11082) |
| HiC-Pro paper | [Servant et al. 2015, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x) |
| cooler / cooltools | [Abdennur & Mirny 2020, Bioinformatics](https://academic.oup.com/bioinformatics/article/36/1/311/5530598) |
| HiCExplorer | [Ramírez et al. 2018, Nature Communications](https://www.nature.com/articles/s41467-017-02525-w) |
| Juicer / HiCCUPS | [Rao et al. 2014, Cell](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4) |
| MUSTACHE loop caller | [Roayaei Ardakany et al. 2020, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02167-0) |
| HiGlass visualization | [Kerpedjiev et al. 2018, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1486-1) |
