---
layout: page
title: Vibe Coding in Life Sciences
published: true
---

{% include gh_variables.html %}

> ## Source
> This content is adapted from ["When Vibe Coding Meets Life Science"](https://www.linkedin.com/pulse/when-vibe-coding-meets-life-science-gozde-eskici-phd)
> by Gozde Eskici, Ph.D. (The Second Translation newsletter, April 14, 2025)
{: .callout}

## What is Vibe Coding?

On February 2nd, 2025, **Andrej Karpathy**, one of the most grounded and trusted voices in AI, posted:

> "There's a new kind of coding I call 'vibe coding', where you fully give in to the vibes, embrace exponentials, and forget that the code even exists. It's possible because the LLMs (e.g. Cursor Composer w Sonnet) are getting too good. Also I just talk to Composer with SuperWhisper so I barely even touch the keyboard... I'm building a project or webapp, but it's not really coding – I just see stuff, say stuff, run stuff, and copy paste stuff, and it mostly works."

And just like that, a new term entered the tech lexicon.

**IBM** followed with a crisp definition:

> "Vibe coding is a fresh take in coding where users express their intention using plain speech and the AI transforms that thinking into executable code."

## Vibe Coding in Biotech

In biotech, coding isn't about building websites—it's about running genome pipelines, training disease models, or scripting CRISPR screens. Historically, that's meant technical depth, time, and a bioinformatics team.

But what if a scientist could just say:

- "Compare these RNA-seq datasets."
- "Predict disease progression from this clinical data."
- "Simulate how this protein binds ligands."

...and the AI handles the rest?

Thanks to LLMs, Biopython, and Colab-powered interfaces, we're now close. The act of building has become more conversational, more iterative—more "vibey."

## Why This Matters

Biotech has long been bottlenecked by translation—the gap between idea and execution. Vibe coding changes that by:

- **Speeding up iteration** on experiments and product ideas
- **Lowering the technical barrier** for scientists-turned-founders
- **Enabling more people** to prototype, test, and scale ideas
- **Unlocking leaner, faster teams**—especially at the seed stage

In short: the next great biotech startup might not hire an engineer first. They might just vibe their way to an MVP.

## 8 Tools Rewriting the Rules of Life Sciences

### 1. Superbio.ai – The No-Code AI Marketplace for Life Sciences

Founded by Berke Buyukkucak and Ronjon Nag, Superbio lets you run cutting-edge AI tools for:
- Drug discovery
- Protein design
- Literature review

No code needed. Their GPT-powered Co-pilot helps you find the right tool.

**Link:** [superbio.ai](https://superbio.ai)

---

### 2. Recursion's LOWE – An LLM-Orchestrated Wet Lab

Recursion's internal tool, LOWE (unveiled by Chris Gibson), is a glimpse into AI-native drug discovery:
- Say the assay, LOWE designs + executes it via robotics
- Taps Recursion's proprietary phenomics & chemistry stack
- Scientists use AI tools via natural language

**Link:** [recursion.com](https://recursion.com)

---

### 3. DrBioRight 2.0 – Talk to Your Cancer Proteomics Atlas

Built at MD Anderson by the Han Liang Lab, this LLM-powered chatbot lets researchers ask things like:

> "Which proteins in pathway X are altered in this tumor?"

No code needed—just real answers and even plots.

**Publication:** [Nature Communications (2025)](https://www.nature.com/articles/s41467-025-56650-0)

**Link:** [drbioright.org](https://drbioright.org)

---

### 4. BioChatter – Build Your Own Bio-AI Assistant (Open Source)

From EMBL-EBI, this Python toolkit lets labs build custom AI assistants:
- Connects to APIs, databases, bio tools
- No coding needed to use
- Fully open-source + on-prem ready

**Link:** [biochatter.org](https://biochatter.org)

---

### 5. OLAF – Conversational Bioinformatics OS

From Weill Cornell by Dylan Riffle et al., OLAF lets you say:

> "Analyze this RNA-seq file"

And it:
- Writes + runs the code
- Returns results and visuals
- Keeps it transparent and inspectable

**Publication:** [ArXiv](https://arxiv.org/abs/2503.12465)

---

### 6. TinyBio (acquired by Seqera) – The "ChatGPT for Scientists"

A platform started by Sasha Dagayev and Vishal Patel in 2022 with real-time code execution:
- Load your data, describe the task
- Supports 50+ bio libraries (Scanpy, Seurat, etc.)
- Self-heals errors and annotates code for you

**Link:** [tinybio.cloud](https://tinybio.cloud)

---

### 7. Scispot (Scibot) – Your Lab's AI Analyst

YC-backed startup embedding AI into lab workflows. When brothers Guru and Satya watched their mother wait too long for a life-saving treatment, they saw firsthand how the slow pace of science could cost lives. That frustration sparked a bold question: "What if labs ran on intelligence instead of inefficiency?"

Their AI assistant, Scibot, makes lab data conversational:
- "Summarize this week's PCR results" → done
- Flags issues, builds dashboards—no Python needed

**Link:** [scispot.com](https://scispot.com)

---

### 8. Synthace – Conversational Automation for Wet Labs

A digital experiment pioneer now powered by LLMs:
- Describe experiments in plain English
- AI generates protocols + sends them to robots

**Link:** [synthace.com](https://synthace.com)

---

## TL;DR

> ## Key Takeaways
> - **Vibe coding** lets users build through intent, not syntax
> - In **biotech**, that means less friction, faster feedback, and broader access
> - This shift expands the pool of technical founders and unlocks whitespace across bio-AI tooling, platforms, and infrastructure
>
> The tools above don't just "assist" scientists—they enable them to do more with less code and more creativity. If vibe coding is the new interface for software, these teams are translating it into biology.
{: .callout}

---

# BCH709 Bioinformatics Vibe Coding Lab Materials

> ## Lab Overview
> **Audience**: BCH709 Genome Informatics graduate students
> **Goal**: Experience how prompt specificity transforms code quality and output
> **Core Lesson**: The more specific your prompt, the closer the result to what you actually need
> **Structure**: 2 Examples (Python 1, R 1) + 2 Homework Assignments (Python 1, R 1)
{: .callout}

## The Vibe Coding Workflow

Vibe coding is a programming approach where you describe a task in natural language, an AI generates the code, you run it and inspect the output, then refine the prompt and iterate.

```
Natural-language prompt → AI generates code → Execute → Inspect results → Revise prompt → Repeat
```

**The key to vibe coding is "saying exactly what you want."**
A prompt controls not only the code but also the **execution environment** — and without that, reproducibility breaks down.

---

# Step 0: Brainstorming and Environment Setup

> ## Lesson for this step
> Before writing any code, ask the AI about possible approaches and required tools first.
{: .callout}

---

## Step 0A. Brainstorming Prompt (For Students — Paste Directly into AI)

Copy and paste the prompt below into the AI first. **This is a strategy question, not a code request.**

~~~
I am a beginner student in BCH709.
Before writing any code, brainstorm the approaches and libraries I need for the following analyses.

Analysis 1 (Python, GFF3 analysis):
- Input: MGI.gff3.gz (GFF3, gzip), chrom.sizes (TSV: chrom, length_bp)
- Goal: Count genes, exons (preventing isoform overcounting), snRNAs, and lncRNAs per chromosome; compute density
- Output: TSV table + dropped_seqids.txt (QC artifact)

Analysis 2 (R, RNA-Seq TPM analysis):
- Input: iceplant_TPM_DT_ZT.tab.gz (gzip TSV, gene_id + sample TPM values)
- Goal: Select top 200 genes by CV (with mean >= 1 filter), save log2(TPM+1) heatmap
- Output: TSV + heatmap PNG

Requirements:
1) Break each analysis into functional units (input parsing, statistical computation, visualization, file output, QC).
2) For each functional unit, suggest 1–2 candidate libraries.
3) Pick one recommended combination for beginners and explain why.
4) List the exact conda-forge package names for that combination.
~~~

**Why this prompt works:**
- Students don't get stuck trying to recall package names they've never heard of
- The AI produces a structured "function → package" mapping that justifies each choice
- The conda install command follows naturally from the output

---

## Step 0B. Environment Setup Prompt (For Students)

Once the brainstorming results are in, use the following prompt to get the install commands:

~~~
Using the library combination you just recommended, generate conda environment creation commands.

Conditions:
- Python environment name: bch709-python
- R environment name: bch709-R
- Pin Python 3.11, R 4.3
- Include import/library verification tests after installation
- Present the commands in copy-paste order so a beginner can just run them one by one
~~~

---

## Step 0C. Environment Creation (Beginner-Friendly Commands)

### For Analysis 1: `bch709-python`

```bash
# Create environment
conda create -n bch709-python -y python=3.11
conda activate bch709-python

# Install packages
conda install -c conda-forge -y \
  pandas numpy matplotlib seaborn biopython tqdm

# Verify installation
python -c "import pandas, numpy, matplotlib, Bio; print('bch709-python OK')"

# Create working directories
mkdir -p data results
```

### For Analysis 2: `bch709-R`

```bash
# Create environment
conda create -n bch709-R -y -c conda-forge \
  r-base=4.3 r-data.table r-ggplot2 r-pheatmap r-viridislite r-scales
conda activate bch709-R

# Verify installation
R -q -e 'library(data.table); library(ggplot2); library(pheatmap); cat("bch709-R OK\n")'

# Create working directories
mkdir -p results
```

### Data Downloads

```bash
# GFF3 (Example 1, Homework 1)
curl -L -o data/MGI.gff3.gz http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz

# mRNA FASTA (Homework 1)
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz

# Ice plant TPM (Example 2, Homework 2)
# https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling
# iceplant_TPM_DT_ZT.tab.gz
```

---

### How to Include the Environment in Every Prompt

From this point forward, **start every prompt** with a sentence like this:

~~~
Write Python code that runs in the bch709-python conda environment.
Assume pandas, numpy, matplotlib, and biopython are installed.
~~~

Or for R:

~~~
Write R code that runs in the bch709-R conda environment.
You may use data.table, ggplot2, and pheatmap.
~~~

> ## Important
> **If you don't specify the environment, the AI will assume an arbitrary one.**
{: .callout}

---

# Part 1: Vibe Coding Examples

---

## Example 1 (Python): Per-Chromosome Feature Counts from MGI GFF3 + QC

### Background

Extract chromosome-level feature counts (genes, exons, snRNAs, lncRNAs) from the MGI GFF3 file and cross-reference against an external chrom.sizes file to verify data integrity.

**GFF3 file structure (9 tab-separated columns):**

```
chr1  MGI  gene  3214482  3671498  .  -  .  ID=MGI:MGI:1918911;Name=Xkr4;biotype=protein_coding
```

| Column | Content |
|--------|---------|
| 1 | chromosome (seqid) |
| 2 | source |
| 3 | feature type (gene, mRNA, exon, snRNA, lnc_RNA, etc.) |
| 4 | start position |
| 5 | end position |
| 6 | score |
| 7 | strand (+/−) |
| 8 | phase |
| 9 | attributes (key=value pairs, semicolon-delimited) |

**chrom.sizes file structure (2 tab-separated columns):**
```
chr1    195465000
chr2    182105000
chrX    171031299
```

### Critical Design Decisions (Definitions That Must Appear in Your Prompt)

**1. Preventing exon overcounting:**
- When a gene has multiple transcript isoforms in GFF3, exon features are repeated per transcript
- **Definition**: Count the number of unique (start, end, strand) exon intervals per chromosome

**2. Chromosome length source:**
- Do NOT estimate lengths from the GFF — use an **external chrom.sizes file** as the reference
- Any GFF seqid not present in chrom.sizes is excluded and logged to `dropped_seqids.txt`

**3. snRNA/lncRNA definition:**
- Count lines where the feature type is `snRNA`, `lnc_RNA`, or `lncRNA` (type-based, first-pass definition)

---

### Stage 1: Vague Prompt

**Prompt:**
> "Extract the gene, exon, snRNA, and lncRNA counts per chromosome from the GFF3 file."

**AI-generated code:**
```python
import gzip, re
from collections import defaultdict

counts = defaultdict(lambda: defaultdict(int))

with gzip.open("data/MGI.gff3.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom = fields[0]
        ftype = fields[2]
        if ftype in ("gene", "exon", "snRNA", "lncRNA"):
            counts[chrom][ftype] += 1

for chrom in sorted(counts):
    print(chrom, dict(counts[chrom]))
```

**Problems:**
- No chromosome length reference — scaffolds and contigs are included alongside primary chromosomes
- Exons are double-counted across isoforms
- May miss `lnc_RNA` (the standard GFF3 spelling)
- No tracking of dropped seqids → no QC possible
- No file output

---

### Stage 2: Specifying Output Format

**Prompt:**
~~~
Write Python code that runs in the bch709-python conda environment.
pandas and numpy are installed.

Input:
- GFF3 gzip: data/MGI.gff3.gz (9 columns, tab-separated)
- Chrom sizes: data/chrom.sizes (TSV: chrom, length_bp)

Task:
- Only include chromosomes present in chrom.sizes
- Save GFF seqids not found in chrom.sizes to results/dropped_seqids.txt
- Count genes per chromosome (type == "gene")
- Count exons per chromosome: to prevent isoform overcounting, define count as number of unique (start, end, strand) intervals
- Count snRNA/lncRNA: lines where type is "snRNA", "lnc_RNA", or "lncRNA"

Output:
- results/chr_feature_counts.tsv
- Columns: chrom, chrom_length_bp, n_gene, n_exon_unique, n_snRNA, n_lncRNA
- Include header
- Print top 5 rows to console
~~~

**AI-generated code:**
```python
from pathlib import Path
import gzip
import pandas as pd
from collections import defaultdict

DATA_GFF  = Path("data/MGI.gff3.gz")
DATA_SIZES = Path("data/chrom.sizes")
OUT_TSV   = Path("results/chr_feature_counts.tsv")
OUT_DROP  = Path("results/dropped_seqids.txt")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# Load chrom sizes
sizes = pd.read_csv(DATA_SIZES, sep="\t", header=None, names=["chrom","length_bp"])
chrom_set = set(sizes["chrom"])
chrom_len = dict(zip(sizes["chrom"], sizes["length_bp"]))

n_gene = defaultdict(int)
exon_intervals = defaultdict(set)
n_snRNA = defaultdict(int)
n_lncRNA = defaultdict(int)
dropped_seqids = set()

SNRNA_TYPES  = {"snRNA"}
LNCRNA_TYPES = {"lnc_RNA", "lncRNA"}

with gzip.open(DATA_GFF, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqid, _, ftype, start, end, _, strand, _, _ = fields

        if seqid not in chrom_set:
            dropped_seqids.add(seqid)
            continue

        if ftype == "gene":
            n_gene[seqid] += 1
        elif ftype == "exon":
            exon_intervals[seqid].add((int(start), int(end), strand))
        elif ftype in SNRNA_TYPES:
            n_snRNA[seqid] += 1
        elif ftype in LNCRNA_TYPES:
            n_lncRNA[seqid] += 1

rows = []
for chrom in sizes["chrom"]:
    rows.append({
        "chrom": chrom,
        "chrom_length_bp": int(chrom_len[chrom]),
        "n_gene": n_gene.get(chrom, 0),
        "n_exon_unique": len(exon_intervals.get(chrom, set())),
        "n_snRNA": n_snRNA.get(chrom, 0),
        "n_lncRNA": n_lncRNA.get(chrom, 0),
    })

df = pd.DataFrame(rows)
df.to_csv(OUT_TSV, sep="\t", index=False)
OUT_DROP.write_text("\n".join(sorted(dropped_seqids)) + "\n")

print("Saved:", OUT_TSV)
print("Dropped seqids:", len(dropped_seqids))
print(df.head(5).to_string(index=False))
```

**Improvements:** TSV output, unique exon definition, dropped seqid QC
**Still missing:** Density not computed, no sorting, unclear handling of chromosomes with zero features

---

### Stage 3: Adding Detailed Constraints

**Prompt:**
~~~
Extend the Stage 2 results.

Additional output requirements:
- Add density columns:
  - gene_per_Mb = n_gene / (chrom_length_bp / 1e6)
  - exon_unique_per_Mb, snRNA_per_Mb, lncRNA_per_Mb (same formula)
- Sort by gene_per_Mb descending
- Round to 4 decimal places
- Include all chromosomes from chrom.sizes, filling counts with 0 where no features exist
- Save unique dropped seqids sorted alphabetically to dropped_seqids.txt
- Print to console: number of dropped seqids, number of excluded feature lines, top 5 rows of the result table
~~~

**AI-generated code (final):**
```python
from pathlib import Path
import gzip
import pandas as pd
from collections import defaultdict

DATA_GFF  = Path("data/MGI.gff3.gz")
DATA_SIZES = Path("data/chrom.sizes")
OUT_TSV   = Path("results/chr_feature_counts.tsv")
OUT_DROP  = Path("results/dropped_seqids.txt")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# 1) Chrom sizes
sizes = pd.read_csv(DATA_SIZES, sep="\t", header=None, names=["chrom","length_bp"])
chrom_set = set(sizes["chrom"])
chrom_len = dict(zip(sizes["chrom"], sizes["length_bp"]))

# 2) Counters
n_gene = defaultdict(int)
exon_intervals = defaultdict(set)
n_snRNA = defaultdict(int)
n_lncRNA = defaultdict(int)
dropped_seqids = set()
dropped_lines = 0

SNRNA_TYPES  = {"snRNA"}
LNCRNA_TYPES = {"lnc_RNA", "lncRNA"}

with gzip.open(DATA_GFF, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqid, _, ftype, start, end, _, strand, _, _ = fields

        if seqid not in chrom_set:
            dropped_seqids.add(seqid)
            dropped_lines += 1
            continue

        if ftype == "gene":
            n_gene[seqid] += 1
        elif ftype == "exon":
            exon_intervals[seqid].add((int(start), int(end), strand))
        elif ftype in SNRNA_TYPES:
            n_snRNA[seqid] += 1
        elif ftype in LNCRNA_TYPES:
            n_lncRNA[seqid] += 1

# 3) Build result table (include all chroms from chrom.sizes; fill 0 where no features)
rows = []
for chrom in sizes["chrom"]:
    L = float(chrom_len[chrom])
    g  = n_gene.get(chrom, 0)
    ex = len(exon_intervals.get(chrom, set()))
    sn = n_snRNA.get(chrom, 0)
    ln = n_lncRNA.get(chrom, 0)
    Mb = L / 1e6 if L > 0 else 1

    rows.append({
        "chrom": chrom,
        "chrom_length_bp": int(L),
        "n_gene": g,
        "n_exon_unique": ex,
        "n_snRNA": sn,
        "n_lncRNA": ln,
        "gene_per_Mb":        round(g  / Mb, 4),
        "exon_unique_per_Mb": round(ex / Mb, 4),
        "snRNA_per_Mb":       round(sn / Mb, 4),
        "lncRNA_per_Mb":      round(ln / Mb, 4),
    })

df = pd.DataFrame(rows).sort_values("gene_per_Mb", ascending=False)
df.to_csv(OUT_TSV, sep="\t", index=False)
OUT_DROP.write_text("\n".join(sorted(dropped_seqids)) + "\n")

print(f"Saved: {OUT_TSV}")
print(f"Saved: {OUT_DROP}")
print(f"Dropped seqids: {len(dropped_seqids)}")
print(f"Dropped feature lines: {dropped_lines}")
print(df.head(5).to_string(index=False))
```

**Sample output:**
```
Saved: results/chr_feature_counts.tsv
Saved: results/dropped_seqids.txt
Dropped seqids: 47
Dropped feature lines: 1823

 chrom  chrom_length_bp  n_gene  n_exon_unique  n_snRNA  n_lncRNA  gene_per_Mb  ...
 chr11      122082543     2847        42156        12        45      23.3224  ...
 chr19       61431566     1892        31245         8        32      30.8024  ...
 chr17       94987271     2456        38912        10        38      25.8563  ...
  chr1      195465000     3215        52341        18        67      16.4488  ...
  chr2      182105000     2987        48562        15        58      16.4027  ...
```

---

### Example 1 Comparison Summary

| Aspect | Stage 1 (Vague) | Stage 2 (Format specified) | Stage 3 (Detailed constraints) |
|--------|-----------------|---------------------------|-------------------------------|
| Chromosome scope | Everything (incl. scaffolds) | chrom.sizes only | chrom.sizes only + zero-fill |
| Exon definition | Duplicate-counted | Unique interval | Unique interval |
| QC artifact | None | dropped_seqids.txt | Dropped count + line count + file |
| Density | None | None | 4 per_Mb columns |
| Sorting | None | None | gene_per_Mb descending |
| Reusability | Low | Medium | High (publication-ready) |

### QC Interpretation Points (Questions Students Must Answer)

1. **What seqids ended up in dropped_seqids.txt?** (Alternative contigs? Unplaced scaffolds? Mitochondrial?)
2. **What fraction of total genes were dropped?** Could this fraction affect the analysis conclusions?
3. **If the prompt had NOT specified using chrom.sizes**, what errors could have occurred **silently**?

> ## Key Message
> `dropped_seqids.txt` is not a failure log — it is a **declaration of what this analysis accepts as data**.
{: .callout}

---

## Example 2 (R): Top 200 Variable Genes from Ice Plant TPM + Heatmap

### Background

Analyze time-course (Zeitgeber Time) TPM expression data from ice plant (*Mesembryanthemum crystallinum*). To keep computation and interpretation manageable, **extract only the top 200 genes by coefficient of variation (CV)** and visualize their expression patterns.

```bash
# Data source
# https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling
# iceplant_TPM_DT_ZT.tab.gz
```

**Data structure:**

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier (e.g., Mc01g0010) |
| DT_ZT2_rep1 through DT_ZT22_rep3 | 6 time points × 3 replicates = 18 samples, TPM values |

- DT = Day Treatment, ZT = Zeitgeber Time
- ZT2, ZT6, ZT10, ZT14, ZT18, ZT22 (4-hour intervals)

### Critical Design Decisions

**1. mean >= 1 filter**: Genes with near-zero TPM can produce infinite or artificially extreme CV values
**2. CV = sd / mean**: Measures relative variability independent of absolute expression level
**3. log2(TPM + 1) transformation**: Reduces the influence of extreme values in heatmap visualization

---

### Stage 1: Vague Prompt

**Prompt:**
> "Pick the 200 most variable genes from the TPM data and draw a heatmap."

**AI-generated code:**
```r
data <- read.table("iceplant_TPM_DT_ZT.tab", header=TRUE, sep="\t", row.names=1)
vars <- apply(data, 1, var)
top200 <- head(sort(vars, decreasing=TRUE), 200)
heatmap(as.matrix(data[names(top200), ]))
```

**Problems:**
- Uses variance instead of CV → biased toward highly expressed genes
- No mean >= 1 filter → noise genes creep into the top ranks
- No file output, no image size/resolution control
- No log2 transformation → heatmap color scale is dominated by extreme values

---

### Stage 2: Specifying Output Format

**Prompt:**
~~~
Write R code that runs in the bch709-R conda environment.
data.table and pheatmap are installed.

Input:
- Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz (gzip TSV)
- First column: gene_id; remaining columns: numeric TPM

Task:
- Compute mean_tpm and sd_tpm for each gene
- CV = sd_tpm / mean_tpm
- Keep only genes with mean_tpm >= 1
- Select top 200 by CV descending

Output:
- Save to results/iceplant_cv_top200.tsv
- Columns: gene_id, mean_tpm, sd_tpm, cv
- Round to 4 decimal places
- Print top 10 to console
~~~

**AI-generated code:**
```r
library(data.table)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz"
dt <- fread(f)
stopifnot("gene_id" %in% names(dt))

sample_cols <- setdiff(names(dt), "gene_id")
mat <- as.matrix(dt[, ..sample_cols])
mode(mat) <- "numeric"

mean_tpm <- rowMeans(mat, na.rm = TRUE)
sd_tpm   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- mean_tpm >= 1
res <- data.table(
  gene_id  = dt$gene_id[keep],
  mean_tpm = mean_tpm[keep],
  sd_tpm   = sd_tpm[keep]
)
res[, cv := sd_tpm / mean_tpm]
setorder(res, -cv)

top200 <- res[1:min(200, .N)]
top200[, `:=`(
  mean_tpm = round(mean_tpm, 4),
  sd_tpm   = round(sd_tpm, 4),
  cv       = round(cv, 4)
)]

fwrite(top200, "results/iceplant_cv_top200.tsv", sep = "\t")
print(top200[1:min(10, .N)])
cat("Saved: results/iceplant_cv_top200.tsv\n")
```

**Improvements:** Explicit CV definition, mean >= 1 filter, TSV output
**Still missing:** No heatmap generated, no plot specifications

---

### Stage 3: Adding Detailed Constraints

**Prompt:**
~~~
Using the top 200 genes from Stage 2, generate a heatmap.

Heatmap data:
- Values: log2(TPM + 1)
- Rows: gene_id (maintain CV descending order, cluster_rows = FALSE)
- Columns: original sample column order (cluster_cols = FALSE)

Plot specifications:
- Save to: results/iceplant_cv_top200_heatmap.png
- Size: 1800 × 1200 px, dpi 200
- Rotate x-axis labels 90 degrees
- Title: "Ice plant log2(TPM+1), CV top200 (mean>=1)"
- Also print top 10 summary table to console
~~~

**AI-generated code (final):**
```r
library(data.table)
library(pheatmap)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz"
dt <- fread(f)
stopifnot("gene_id" %in% names(dt))

sample_cols <- setdiff(names(dt), "gene_id")
mat <- as.matrix(dt[, ..sample_cols])
mode(mat) <- "numeric"

mean_tpm <- rowMeans(mat, na.rm = TRUE)
sd_tpm   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- mean_tpm >= 1
res <- data.table(
  gene_id  = dt$gene_id[keep],
  mean_tpm = mean_tpm[keep],
  sd_tpm   = sd_tpm[keep]
)
res[, cv := sd_tpm / mean_tpm]
setorder(res, -cv)
top200 <- res[1:min(200, .N)]

# Summary TSV
top200_out <- copy(top200)
top200_out[, `:=`(
  mean_tpm = round(mean_tpm, 4),
  sd_tpm   = round(sd_tpm, 4),
  cv       = round(cv, 4)
)]
fwrite(top200_out, "results/iceplant_cv_top200.tsv", sep = "\t")

# Heatmap matrix
idx <- match(top200$gene_id, dt$gene_id)
submat <- mat[idx, , drop = FALSE]
submat <- log2(submat + 1)
rownames(submat) <- top200$gene_id

# Save (1800×1200 px, dpi 200)
png("results/iceplant_cv_top200_heatmap.png",
    width = 1800, height = 1200, res = 200)
pheatmap(
  submat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 6,
  fontsize_row = 3,
  angle_col = 90,
  main = "Ice plant log2(TPM+1), CV top200 (mean>=1)"
)
dev.off()

print(top200_out[1:min(10, .N)])
cat("Saved: results/iceplant_cv_top200.tsv\n")
cat("Saved: results/iceplant_cv_top200_heatmap.png\n")
```

---

### Example 2 Comparison Summary

| Aspect | Stage 1 (Vague) | Stage 2 (Format specified) | Stage 3 (Detailed constraints) |
|--------|-----------------|---------------------------|-------------------------------|
| Variability metric | Variance | CV (sd/mean) | CV (sd/mean) |
| Filtering | None | mean >= 1 | mean >= 1 |
| Data transformation | None | None | log2(TPM+1) |
| File output | None | TSV | TSV + PNG (size/dpi specified) |
| Reusability | Low | Medium | High |

### Interpretation Points

- Without the `mean >= 1` filter, genes with near-zero expression but a single spike dominate the top ranks
- CV captures **relative variability independent of absolute expression** → enables fair comparison across expression levels
- `log2(TPM+1)` transformation equalizes the heatmap color distribution, making expression patterns visible

---

# Part 2: Homework Assignments

---

## Homework 1 (Python): mRNA FASTA Analysis + GC Distribution Graph Page

### Problem Description

Extract sequence information from the UCSC human mRNA FASTA file (mrna.fa.gz), analyze GC content distribution, and produce a graph and an HTML report page.

### Objective

**Write a single prompt using vibe coding that produces the desired result in one shot.** The goal is to get all outputs correct from a single, well-crafted prompt.

### Input Data

```bash
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
```

**FASTA file structure:**
```
>BC032353 /gb=BC032353 /gi=34783011 /ug=Hs.444613 /len=1873
CGGCAGCGGCTGCGGGGAGGATGGCGGCGACGGCGACTTTTAAATATATTGG...
```

### Expected Output

1. **`results/mrna_metrics.tsv`**: accession, length, gc_content (4 decimal places, sorted by gc_content descending, top 20 only)
2. **`results/gc_content_distribution.png`**: Histogram + density curve (1600×900 px, dpi 200)
   - Bins: 0.00 to 1.00, step 0.01
   - Bar color: light gray; border: dark gray
   - Density curve: blue, line width 2.0
   - Mean: red dashed vertical line; Median: green dashed vertical line
   - Caption showing n, mean, median, sd
3. **`results/gc_content_distribution.svg`**: Same graph in SVG format
4. **`docs/gc_content_distribution.html`**: Title + summary statistics + embedded PNG image

### Chart Selection Discussion (In-Class Activity)

Use the following prompt to have the AI compare chart types before making a graph:

~~~
I want to show the distribution of gc_content values.
Compare histogram, density curve, and ECDF — which is most appropriate?

My situation:
- Continuous values between 0 and 1, very large sample size (tens of thousands)
- I want to convey distribution shape, central tendency (mean, median), and tails
- Must be reproducible at publication-quality level

For each option, give pros and cons, then recommend your top pick with justification.
~~~

Reference: [R Graph Gallery — Distribution section](https://r-graph-gallery.com/)

### Grading Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Prompt quality | 40% | Did you specify input format, accession parsing rule, gzip handling, and output specs (filename, columns, decimals, sorting, graph size/color/font)? |
| Code correctness | 40% | Correct gzip streaming parse, accurate computation, all 4 files generated |
| Result interpretation | 20% | Explain 3 possible biological/technical reasons why the top-20 GC content values are high |

---

## Homework 2 (R): Z-Score Clustering of CV Top 200 Genes + Pattern Visualization

### Problem Description

Using the 200 genes from the `results/iceplant_cv_top200.tsv` file generated in Example 2:

1. **Z-score normalize** (row-wise)
2. **Hierarchical clustering** (ward.D2 method, euclidean distance)
3. **Cut tree at k=4** to assign clusters
4. **Save cluster mean expression pattern** line plots
5. **Save cluster assignment** table as TSV

### Objective

Practice writing prompts that **precisely control R visualization output**.

### Input Data

```bash
# TSV from Example 2 + original TPM data
# results/iceplant_cv_top200.tsv (used as gene_id list)
# Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz (original TPM)
```

### Expected Output

**1. Clustered Heatmap** — `results/cv_top200_cluster_heatmap.pdf`

- Data: Mean TPM per ZT time point → row-wise Z-score normalization
  - Z-score = (value − row_mean) / row_sd
- Rows: gene_id (hierarchical clustering, ward.D2, euclidean)
- Columns: ZT time points (ZT2, ZT6, ZT10, ZT14, ZT18, ZT22) — fixed in chronological order
- Colors: blue (low) → white (0) → red (high)
- k=4 cutree result displayed as a color bar to the left of rows
- PDF size: 8 × 12 inches

**2. Cluster Mean Pattern Line Plot** — `results/cluster_patterns.pdf`

- 2×2 panel layout (facet_wrap)
- Each panel: Mean TPM ± SD across ZT for genes in that cluster
- X-axis: ZT (2, 6, 10, 14, 18, 22)
- Y-axis: Mean TPM
- Panel title format: "Cluster 1 (n=XX genes)"
- Include error bars
- PDF size: 10 × 8 inches

**3. Assignment Table** — `results/cluster_assignment.tsv`

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| cluster | 1–4 |
| peak_ZT | ZT with highest mean TPM |
| trough_ZT | ZT with lowest mean TPM |
| amplitude | max(mean_TPM) − min(mean_TPM), 2 decimal places |

- Sort by cluster number ascending, then amplitude descending within each cluster

### Prompt-Writing Hints

Definitions that must appear in your prompt:

~~~
Analysis procedure:
1. From iceplant_TPM_DT_ZT.tab.gz, extract only the top 200 genes (by gene_id list)
2. Parse ZT time point from sample names (regex: digits after "ZT")
3. Average the 3 replicates per ZT → gene × ZT matrix (200 rows × 6 columns)
4. Z-score normalize: for each row (gene), compute (value - mean) / sd
5. Hierarchical clustering: dist(euclidean) → hclust(ward.D2)
6. cutree(k=4) to assign 4 clusters

For the heatmap:
- Use pheatmap or ComplexHeatmap
- Show cluster assignment as annotation_row color bar
- Columns (ZT) in chronological order (cluster_cols = FALSE)

For the line plot:
- Use original TPM values (no log transformation) for mean ± SD
- ggplot2 facet_wrap(~cluster, ncol=2)
~~~

### Grading Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Prompt quality | 40% | Did you specify Z-score definition, clustering method (ward.D2, euclidean), k=4, replicate-to-mean procedure, output filenames and specs? |
| Code correctness | 40% | Accurate Z-score normalization, distance computation, clustering, cutree, summary statistics, all 3 files correctly generated |
| Result interpretation | 20% | Interpret each of the 4 cluster expression patterns in the context of CAM photosynthesis (2 sentences per cluster, e.g., "Cluster 2 peaks during nighttime (ZT14–ZT22) → candidate CAM-related genes") |

---

# Appendix: Effective Vibe Coding Prompt Template

~~~
[Environment]: Write [language] code that runs in the [env name] conda environment.
[Installed packages] are available.

**Input specification:**
- File: [filename] ([format], [delimiter], [gzip?])
- Structure: [column descriptions, special parsing rules]
- Additional inputs: [chrom.sizes, gene lists, or other reference files]

**Analysis conditions:**
- [Filter criteria (e.g., mean >= 1)]
- [Computation method (e.g., CV = sd/mean)]
- [Definitions (e.g., exon count = unique intervals only)]

**Output 1 — Table:**
- Filename: [filename]
- Columns: [list column names]
- Decimal places: [number]
- Sorting: [criterion, direction]
- Filter: [top N]

**Output 2 — Plot:**
- Filename: [filename]
- Size: [px or inches], dpi: [value]
- Colors: [specify explicitly]
- Axes/labels/legend: [specify explicitly]

**Output 3 — QC:**
- [Dropped items filename]
- [Summary information to print to console]
~~~

---

# Input/Output Prompt Checklist

### When specifying input:
- Filename and path
- File format (GFF3, TSV, FASTA, etc.)
- Compression (gzip or not)
- Delimiter (tab, comma, space)
- Header presence
- Data structure (column names, what rows represent)
- Special structures (e.g., GFF3 attribute parsing rules)
- External reference files (chrom.sizes, etc.)

### When specifying output:
- File format (TSV, CSV, PDF, PNG, SVG, HTML)
- Filename
- Column names and order
- Decimal places
- Sorting criterion (ascending/descending)
- Filter conditions (top N, minimum threshold, etc.)
- Plot: size, resolution, colors, font, legend position, axis range
- QC artifacts (dropped items, summary statistics)

### When specifying analysis definitions:
- Metric definitions (CV = sd/mean, Z-score = (x−mean)/sd)
- Filter rules (mean >= 1)
- Deduplication handling (unique intervals, etc.)
- Transformation methods (log2(TPM+1))
- Clustering parameters (method, distance metric, k)

---

> ## Data Sources
> - MGI GFF3: [http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz](http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz)
> - Human mRNA: [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz)
> - Ice Plant TPM: [https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling](https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling)
{: .callout}
