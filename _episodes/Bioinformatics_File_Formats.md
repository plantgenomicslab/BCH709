---
layout: page
title: Bioinformatics File Formats
published: true
---

{% include gh_variables.html %}

## Bioinformatics File Formats

Understanding common bioinformatics file formats is essential for working with genomic data.

### SAM/BAM Format

SAM (Sequence Alignment Map) is a text format for storing sequence alignments. BAM is its binary, compressed version.

> ## SAM Format Structure
> A SAM file consists of:
> 1. **Header section** (lines starting with @)
> 2. **Alignment section** (tab-delimited fields)
>
> **Header lines:**
> - `@HD` - Header line
> - `@SQ` - Reference sequence dictionary
> - `@RG` - Read group
> - `@PG` - Program used
{: .keypoints}

**Alignment fields (11 mandatory):**

| Col | Field | Description |
|-----|-------|-------------|
| 1 | QNAME | Query name |
| 2 | FLAG | Bitwise flag |
| 3 | RNAME | Reference name |
| 4 | POS | Position |
| 5 | MAPQ | Mapping quality |
| 6 | CIGAR | CIGAR string |
| 7 | RNEXT | Mate reference name |
| 8 | PNEXT | Mate position |
| 9 | TLEN | Template length |
| 10 | SEQ | Sequence |
| 11 | QUAL | Quality string |

> ## Working with SAM/BAM
> ```bash
> # View SAM file
> $ less alignment.sam
>
> # Convert SAM to BAM
> $ samtools view -bS alignment.sam > alignment.bam
>
> # Sort BAM file
> $ samtools sort alignment.bam -o alignment.sorted.bam
>
> # Index BAM file
> $ samtools index alignment.sorted.bam
>
> # View BAM file
> $ samtools view alignment.bam | head
>
> # View specific region
> $ samtools view alignment.sorted.bam chr1:1000-2000
>
> # Get statistics
> $ samtools flagstat alignment.bam
> $ samtools stats alignment.bam
> ```
{: .keypoints}

### BED Format

BED (Browser Extensible Data) format is used to define genomic regions.

> ## BED Format Structure
> BED files are tab-delimited with at least 3 columns.
{: .keypoints}

| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome |
| 2 | chromStart | Start position (0-based) |
| 3 | chromEnd | End position |
| 4 | name | Feature name (optional) |
| 5 | score | Score (optional) |
| 6 | strand | + or - (optional) |

**Example:**
```
chr1    1000    2000    gene1    100    +
chr1    3000    4000    gene2    200    -
chr2    5000    6000    gene3    150    +
```

**Important:** BED uses 0-based, half-open coordinates

> ## Working with BED Files
> ```bash
> # Sort BED file
> $ sort -k1,1 -k2,2n input.bed > sorted.bed
>
> # Merge overlapping intervals
> $ bedtools merge -i sorted.bed > merged.bed
>
> # Find intersections
> $ bedtools intersect -a file1.bed -b file2.bed > common.bed
>
> # Subtract regions
> $ bedtools subtract -a file1.bed -b file2.bed > unique.bed
>
> # Get flanking regions
> $ bedtools flank -i genes.bed -g genome.txt -b 1000 > flanks.bed
> ```
{: .keypoints}

### VCF Format

VCF (Variant Call Format) stores genetic variation data.

> ## VCF Format Structure
> VCF files have:
> 1. **Meta-information lines** (starting with ##)
> 2. **Header line** (starting with #CHROM)
> 3. **Data lines** (one per variant)
{: .keypoints}

**Fixed columns:**

| Col | Field | Description |
|-----|-------|-------------|
| 1 | CHROM | Chromosome |
| 2 | POS | Position (1-based) |
| 3 | ID | Variant ID |
| 4 | REF | Reference allele |
| 5 | ALT | Alternate allele(s) |
| 6 | QUAL | Quality score |
| 7 | FILTER | Filter status |
| 8 | INFO | Additional info |
| 9 | FORMAT | Genotype format |
| 10+ | SAMPLE | Sample genotypes |

**Example:**
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
chr1    100     rs123   A       G       30      PASS    DP=50   GT:DP   0/1:50
```

> ## Working with VCF Files
> ```bash
> # View VCF file
> $ bcftools view variants.vcf | head
>
> # Compress and index
> $ bgzip variants.vcf
> $ tabix -p vcf variants.vcf.gz
>
> # Filter variants
> $ bcftools filter -i 'QUAL>30' variants.vcf.gz > filtered.vcf
>
> # Extract specific region
> $ bcftools view variants.vcf.gz chr1:1000-2000 > region.vcf
>
> # Statistics
> $ bcftools stats variants.vcf.gz > stats.txt
>
> # Convert to table
> $ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' variants.vcf.gz
> ```
{: .keypoints}


## Common Bioinformatics Tools

### samtools

samtools is a suite of programs for interacting with SAM/BAM files.

> ## Essential samtools Commands
> ```bash
> # Convert SAM to BAM
> $ samtools view -bS input.sam > output.bam
>
> # Sort BAM file
> $ samtools sort input.bam -o sorted.bam
>
> # Index BAM file (required for many operations)
> $ samtools index sorted.bam
>
> # View alignment statistics
> $ samtools flagstat sorted.bam
>
> # Calculate depth
> $ samtools depth sorted.bam > depth.txt
>
> # Extract reads from region
> $ samtools view sorted.bam chr1:1000-2000 > region.sam
>
> # Extract unmapped reads
> $ samtools view -f 4 sorted.bam > unmapped.sam
>
> # Extract properly paired reads
> $ samtools view -f 2 sorted.bam > proper_pairs.sam
>
> # Merge multiple BAM files
> $ samtools merge merged.bam file1.bam file2.bam file3.bam
>
> # Create FASTA index
> $ samtools faidx reference.fa
>
> # Extract sequence from FASTA
> $ samtools faidx reference.fa chr1:1000-2000
> ```
{: .keypoints}

### bedtools

bedtools is a powerful suite for genomic arithmetic operations.

> ## Essential bedtools Commands
> ```bash
> # Find overlapping features
> $ bedtools intersect -a genes.bed -b peaks.bed > overlaps.bed
>
> # Count overlaps
> $ bedtools intersect -a genes.bed -b peaks.bed -c > counts.bed
>
> # Find features NOT overlapping
> $ bedtools intersect -a genes.bed -b peaks.bed -v > no_overlap.bed
>
> # Merge overlapping intervals
> $ bedtools merge -i sorted.bed > merged.bed
>
> # Calculate coverage
> $ bedtools coverage -a genes.bed -b reads.bam > coverage.bed
>
> # Get closest feature
> $ bedtools closest -a query.bed -b reference.bed > closest.bed
>
> # Generate genome windows
> $ bedtools makewindows -g genome.txt -w 1000 > windows.bed
>
> # Get FASTA sequences for BED regions
> $ bedtools getfasta -fi genome.fa -bed regions.bed > sequences.fa
>
> # Shuffle features randomly
> $ bedtools shuffle -i features.bed -g genome.txt > shuffled.bed
>
> # Compute Jaccard statistic
> $ bedtools jaccard -a file1.bed -b file2.bed
> ```
{: .keypoints}

> ## Combining Tools in Pipelines
> ```bash
> # Example: Find genes with mapped reads and count
> $ bedtools intersect -a genes.bed -b aligned.bam -c | \
>     awk '$4 > 10' | \
>     sort -k4,4rn > highly_expressed.bed
>
> # Example: Get sequences of peaks
> $ bedtools sort -i peaks.bed | \
>     bedtools merge -i - | \
>     bedtools getfasta -fi genome.fa -bed - > peak_sequences.fa
>
> # Example: Calculate mapping statistics per gene
> $ samtools view -F 4 aligned.bam | \
>     bedtools bamtobed -i stdin | \
>     bedtools intersect -a genes.bed -b stdin -c > gene_counts.bed
> ```
{: .challenge}
