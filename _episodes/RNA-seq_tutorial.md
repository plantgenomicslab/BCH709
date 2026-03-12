---
layout: page
title: RNA-Seq Tutorial
published: true
---

> ## Paper Reading
> Please read this paper before class:
> [Conesa et al. (2016) A survey of best practices for RNA-seq data analysis. *Genome Biology* 17:13](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8)
{: .prereq}

---

## Overview: RNA Sequencing

![RNA Sequencing]({{site.baseurl}}/fig/rnaseq.png)

RNA-Seq (RNA Sequencing) is the standard method for measuring genome-wide gene expression. Key characteristics:

1. The transcriptome is spatially and temporally dynamic
2. Data comes from functional units (coding regions)
3. Only a small fraction of the genome is transcribed at any given time

> ## Introduction
> Sequence-based assays of transcriptomes (RNA-seq) are in wide use because of their favorable properties for quantification, transcript discovery, and splice isoform identification, as well as adaptability for numerous more specialized measurements. RNA-Seq studies present challenges shared with prior methods such as microarrays and SAGE tagging, and also new ones specific to high-throughput sequencing platforms.
>
> RNA-Seq experiments are diverse in their aims and design goals, currently including multiple types of RNA isolated from whole cells or from specific sub-cellular compartments or biochemical classes, such as total polyA+ RNA, polysomal RNA, nuclear ribosome-depleted RNA, various size fractions of RNA and a host of others.
{: .prereq}

---

## The RNA-Seq Workflow

![RNA Sequencing Workflow]({{site.baseurl}}/fig/rnaseq_workflow.png)

### Seven Stages of a Data Science Project
1. Define the question of interest
2. Get the data
3. Clean the data
4. Explore the data
5. Fit statistical models
6. Communicate the results
7. Make your analysis reproducible

---

## 1. Experimental Design

### Sample Information
Document the following for each sample:

| Item | Description |
|------|-------------|
| Material type | Tissue, cell line, primary cell type, etc. |
| Treatments | TALENs, CRISPR, drug treatment, etc. |
| Subcellular fraction | Whole cell, nuclear, cytoplasmic, etc. |
| Input amount | Bulk vs. low-input (affects reproducibility) |
| Lot/catalog # | For commercial cell lines |
| Culture protocol | If cells were propagated in vitro |
| QC/phenotyping | Purity confirmation methods |

### RNA Information
Key properties to report:

- RNA type: Total RNA, Poly-A(+), Poly-A(-)
- Size fraction: >200 nt (long RNA) vs. <200 nt (small RNA)
- rRNA depletion: RiboMinus, RiboZero (note kit used)

### Library Preparation Protocol
Document all steps:

- RNA isolation method
- Size selection method
- rRNA removal method
- Oligo-dT selection method
- DNase I treatment

### Replicate Strategy

- **Biological replicates:** Minimum 3; recommended 6–12 (Schurch et al., 2016)
- **Technical replicates:** Not required for RNA-Seq (Marioni et al., 2008)
- **Power analysis:** Use tools like Scotty or RnaSeqSampleSize before starting

> ## Reading Materials
> - [Auer & Doerge (2010) Statistical Design and Analysis of RNA Sequencing Data. *Genetics* 185:405–416](https://www.genetics.org/content/185/2/405)
> - [Busby et al. (2013) Scotty: a web tool for designing RNA-Seq experiments. *Bioinformatics* 29:656–657](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582267/)
> - [Marioni et al. (2008) RNA-seq: an assessment of technical reproducibility. *Genome Research*](https://genome.cshlp.org/content/18/9/1509.full.html)
> - [Schurch et al. (2016) How many biological replicates are needed in an RNA-seq experiment? *RNA* 22:839](https://rnajournal.cshlp.org/content/22/6/839.long)
> - [Zhao et al. (2018) RnaSeqSampleSize: real data based sample size estimation. *BMC Bioinformatics* 19:191](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5)
{: .prereq}

---

## 2. Sequencing Parameters

### Key Parameters to Report

| Parameter | Examples |
|-----------|---------|
| Platform | Illumina NovaSeq, PacBio, Oxford Nanopore |
| Format | Single-end (SE) or Paired-end (PE) |
| Read length | 100 bp, 150 bp, 250 bp |
| Barcode placement | Standard or custom |
| Custom primers | Sequence and position |

![RNA Library Prep]({{site.baseurl}}/fig/sequencing.png)

### Recommended Sequencing Depth

| Library Type | Minimum Aligned Read Pairs |
|-------------|--------------------------|
| Long RNA-Seq (polyA+) | 30 million |
| RAMPAGE | 20 million |
| Small RNA-Seq | 30 million |

### Quantitative Standards (Spike-ins)

ERCC spike-in controls are highly recommended for calibrating quantification, sensitivity, and linearity.

Report:
- Stage of addition (before polyA selection, at cDNA synthesis, or prior to sequencing)
- FASTA file of spike-in sequences
- Source (ERCC, home-made, etc.)
- Concentration of each spike-in

---

## 3. Data Files and FASTQ Format

### What is FASTQ?

FASTQ is a text-based format storing both nucleotide sequences and per-base quality scores.

Each record has exactly **4 lines**:

```
@A00261:180:HL7GCDSXX:2:1101:30572:1047/2          # Line 1: header
AAAATACATTGATGACCATCTAAAGTCTACGGCGTAT...            # Line 2: sequence
+                                                   # Line 3: separator
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF...              # Line 4: quality scores
```

### Illumina Read Header Format

```
@Instrument:RunID:FlowcellID:Lane:Tile:X:Y/ReadNum
```

| Field | Example | Meaning |
|-------|---------|---------|
| Instrument | A00261 | Instrument name |
| RunID | 180 | Run ID |
| FlowcellID | HL7GCDSXX | Flowcell ID |
| Lane | 2 | Flowcell lane |
| Tile | 1101 | Tile number within lane |
| X | 30572 | X-coordinate of cluster |
| Y | 1047 | Y-coordinate of cluster |
| ReadNum | /2 | Pair member (1 or 2) |

### Phred Quality Scores

![FASTQ Quality]({{site.baseurl}}/fig/quality.png)

Quality scores represent the probability of an incorrect base call:

| Phred Score | Error Probability | Accuracy |
|------------|-----------------|---------|
| Q10 | 1 in 10 | 90% |
| Q20 | 1 in 100 | 99% |
| Q30 | 1 in 1,000 | 99.9% |
| Q40 | 1 in 10,000 | 99.99% |

---

## 4. Hands-On: Setting Up the Environment

### Create Conda Environment

```bash
conda create -n rnaseq python=3.11
conda activate rnaseq
```

### Install Tools

```bash
conda install -c bioconda -c conda-forge \
  fastqc trim-galore hisat2 samtools subread multiqc
```

### Create Working Directory

```bash
mkdir -p ~/bch709/rnaseq
cd ~/bch709/rnaseq
```

### Download Example Data

```bash
wget https://www.dropbox.com/s/y7yehmfze1l6cgz/pair1.fastq.gz
wget https://www.dropbox.com/s/xsrth6icapyr4p0/pair2.fastq.gz
ls -lh
```

### Inspect a FASTQ File

```bash
zcat pair2.fastq.gz | head -12
```

---

## 5. Quality Control

### What to Check

| Metric | Tool |
|--------|------|
| Number of reads | FastQC |
| Per-base sequence quality | FastQC |
| Per-sequence GC content | FastQC |
| Per-base N content | FastQC |
| Sequence length distribution | FastQC |
| Adapter content | FastQC |
| Duplication levels | FastQC |

[FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### Run FastQC

```bash
fastqc -t 4 pair1.fastq.gz pair2.fastq.gz
```

### Aggregate QC Reports with MultiQC

![MultiQC]({{site.baseurl}}/fig/multiqc.png)

[MultiQC](https://multiqc.info/) aggregates results from multiple tools into a single interactive report.

```bash
multiqc .
```

### Open Reports in Your Browser

MultiQC generates an HTML report in your working directory. Open it directly:

**macOS / Linux:**
```bash
open multiqc_report.html
```

**Windows (WSL):**
```bash
explorer.exe multiqc_report.html
```

---

## 6. Read Trimming

### Why Trim?

- Remove adapter sequences
- Trim low-quality bases from read ends
- Remove reads that are too short (< 18 nt)
- Improve downstream alignment accuracy

![Trimming]({{site.baseurl}}/fig/trim.png)

### Common Trimming Tools

| Tool | Link |
|------|------|
| Trim Galore | [link](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) |
| fastp | [link](https://github.com/OpenGene/fastp) |
| Cutadapt | [link](https://github.com/marcelm/cutadapt/) |
| Trimmomatic | [link](http://www.usadellab.org/cms/?page=trimmomatic) |
| Skewer | [link](https://github.com/relipmoc/skewer) |

### Run Trim Galore

```bash
cd ~/bch709/rnaseq

trim_galore \
  --paired \
  --three_prime_clip_R1 5 \
  --three_prime_clip_R2 5 \
  --cores 4 \
  --max_n 40 \
  --fastqc \
  --gzip \
  -o trim \
  pair1.fastq.gz pair2.fastq.gz
```

### Post-Trimming QC Report

```bash
multiqc --dirs ~/bch709/rnaseq --filename trim
```

---

## 7. Alignment (Mapping)

![Mapping]({{site.baseurl}}/fig/mapping.png)

For RNA-Seq, reads must be aligned to a **reference genome** using a **splice-aware** aligner that can span intron-exon junctions.

### Common RNA-Seq Aligners

| Aligner | Algorithm | Link |
|---------|-----------|------|
| HISAT2 | Graph FM index (BWT-based) | [paper](https://www.nature.com/articles/s41587-019-0201-4) |
| STAR | Suffix arrays | [paper](https://academic.oup.com/bioinformatics/article/29/1/15/272537) |
| GSNAP | SNP/indel-aware | [paper](https://dx.doi.org/10.1007/978-1-4939-3578-9_15) |

> [Baruzzo et al. (2017) Simulation-based comprehensive benchmarking of RNA-seq aligners. *Nature Methods* 14:135](https://www.nature.com/articles/nmeth.4106)

![Algorithm Overview]({{site.baseurl}}/fig/algorithm.png)

### Download Reference Files

```bash
cd ~/bch709/rnaseq

wget https://www.dropbox.com/s/851ob9e3ktxhyxz/bch709.fasta
wget https://www.dropbox.com/s/e9dvdkrl9dta4qg/bch709.gtf
```

### Build HISAT2 Index

```bash
hisat2-build bch709.fasta bch709
```

### Align Reads

The `--dta` flag is required when using HISAT2 output with featureCounts or StringTie for downstream quantification. Pipe directly to `samtools sort` to skip writing the intermediate SAM file to disk.

```bash
hisat2 \
  -x bch709 \
  --threads 4 \
  --dta \
  -1 trim/pair1_val_1.fq.gz \
  -2 trim/pair2_val_2.fq.gz \
  --summary-file alignment.txt \
  | samtools sort -@ 4 -o align_sort.bam

samtools index align_sort.bam
cat alignment.txt
```

---

## 8. SAM/BAM Format

### SAM Record Format

```
<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [TAGS]
```

![SAM Format]({{site.baseurl}}/fig/sam1.png)

### SAM Flags

SAM flags are bitwise integers encoding alignment properties (strand, pairing, etc.):

![SAM Flags]({{site.baseurl}}/fig/flag.jpg)

Convert a decimal flag to binary to interpret it:

```bash
python3 -c "print(bin(163))"
# or
echo 'obase=2; 163' | bc
```

Check what a flag means: [SAM Flag Decoder](https://broadinstitute.github.io/picard/explain-flags.html)

### SAM Optional Tags

![SAM Tags]({{site.baseurl}}/fig/samtag.png)

Full SAM specification: [samtools.github.io/hts-specs](http://samtools.github.io/hts-specs/)

---

## 9. BAM Processing with SAMtools

SAMtools provides utilities for manipulating SAM/BAM files: sorting, indexing, merging, and statistics.

> If you aligned without piping (produced a SAM file), convert and sort it:
> ```bash
> samtools view -Sb align.sam | samtools sort -@ 4 -o align_sort.bam
> samtools index align_sort.bam
> ```

```bash
# Compute alignment statistics
samtools stats align_sort.bam > align_sort.bam.stat
cat align_sort.bam.stat
```

### File Size Comparison

| Format | Size | Notes |
|--------|------|-------|
| SAM (align.sam) | ~903 MB | Text format; avoid writing to disk if possible |
| BAM (align.bam) | ~166 MB | Binary, ~5× smaller |

### Visualize Alignments

```bash
# Terminal viewer
COLUMNS=150 samtools tview -d t align_sort.bam bch709.fasta
```

GUI Viewers:
- [IGV (Integrative Genomics Viewer)](https://software.broadinstitute.org/software/igv/)
- [Tablet](https://ics.hutton.ac.uk/tablet/)

### Alignment QC

```bash
multiqc --dirs ~/bch709/rnaseq --filename align
```

### Key Alignment Metrics

| Metric | Tool |
|--------|------|
| Mapping rate, pairing | samtools stats |
| Insert size distribution | samtools stats |
| Gene body coverage | RSeQC |
| Strand specificity | RSeQC |
| Biotype counts | QoRTs |
| Sequencing saturation | QoRTs |

Tools: [QoRTs](https://hartleys.github.io/QoRTs/), [RSeQC](http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/), [Qualimap](http://qualimap.bioinfo.cipf.es/)

---

## 10. Read Quantification

![Gene Counts]({{site.baseurl}}/fig/genecount.png)

Reads overlapping annotated gene features are counted as a proxy for gene expression.

![Count Methods]({{site.baseurl}}/fig/count.png)

### Quantification Tools

| Tool | Approach |
|------|---------|
| featureCounts | Fast, summarizes reads per gene/exon |
| HTSeq | Python-based, flexible |
| RSEM | Transcript-level, EM-based |
| Salmon | Quasi-mapping, very fast |
| Kallisto | Pseudoalignment |

### Run featureCounts

The `-p` flag is for paired-end reads. Add `-T 4` to use multiple threads. The `-s` flag sets strandedness (0=unstranded, 1=forward, 2=reverse).

```bash
featureCounts \
  -p \
  -T 4 \
  -a bch709.gtf \
  -o counts.txt \
  align_sort.bam

# View the count matrix (skip the first commented header line)
grep -v "^#" counts.txt | head
```

### PCR Duplicates

For RNA-Seq, PCR duplicates are generally **not removed** because many identical reads reflect true high-abundance transcripts.

- **Do NOT** computationally deduplicate standard RNA-Seq
- Use **UMIs** during library prep if deduplication is needed

### Multi-Mapping Reads

| Strategy | Tool |
|----------|------|
| Discard multi-mappers | featureCounts, HTSeq |
| Distribute counts | Cufflinks |
| Probabilistic (EM) | RSEM |
| Prioritize features | Rcount |

> ## Reference
> - Fu et al. (2018) Elimination of PCR duplicates in RNA-seq using UMIs. *BMC Genomics* 19:531
> - Parekh et al. (2016) The impact of amplification on differential expression. *Scientific Reports* 6:25533
{: .challenge}

---

## 11. Differential Expression Analysis

After generating a count matrix, the next step is testing for statistically significant differences between conditions.

### Common Tools

| Tool | Model | Notes |
|------|-------|-------|
| DESeq2 | Negative binomial | Recommended for small n |
| edgeR | Negative binomial | GLM-based |
| Limma-Voom | Normal (after voom) | Good for large studies |

### Downstream: Functional Analysis

- **Gene Ontology (GO) enrichment** — hypergeometric test
- **GSEA** — gene set enrichment analysis
- **Pathway analysis** — KEGG, Reactome

---

## 12. Full Workflow Summary

```
Raw FASTQ
    │
    ├── FastQC → MultiQC (QC report)
    │
    ├── Trim Galore (adapter/quality trimming)
    │       │
    │       └── FastQC → MultiQC (post-trim QC)
    │
    ├── HISAT2 (alignment to genome)
    │       │
    │       └── SAMtools (sort, index, stats)
    │
    ├── featureCounts (read quantification)
    │
    └── DESeq2 / edgeR (differential expression)
            │
            └── GO / GSEA (functional enrichment)
```

---

## 13. Cleanup

```bash
conda deactivate
# Optional: remove the environment when done
conda env remove --name rnaseq
```

---

## References

| Resource | Link |
|---------|------|
| Conesa et al. (2016) Best practices for RNA-seq | [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8) |
| ENCODE RNA-Seq Standards | [Hong et al. 2016](https://academic.oup.com/database/article-lookup/doi/10.1093/database/baw001) |
| QC Fail blog | [sequencing.qcfail.com](https://sequencing.qcfail.com/) |
| Conda documentation | [docs.conda.io](https://docs.conda.io/en/latest/) |
| BioConda | [bioconda.github.io](https://bioconda.github.io/) |
