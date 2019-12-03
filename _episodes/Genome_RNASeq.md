---
layout: page
title: Genome RNASeq
published: true
---



## RNA reads alignment
The alignment process consists of choosing an appropriate reference genome to map our reads against and performing the read alignment using one of several splice-aware alignment tools such as STAR or HISAT2. The choice of aligner is often a personal preference and also dependent on the computational resources that are available to you.


## Environment
```bash
conda create -n genomernaseq starseqr bioconductor-deseq2 samtools gffread star
conda activate genomernaseq
```
## Location
```
/data/gpfs/assoc/bch709/<YOURID>/genomernaseq
```


## Reads preparation
```bash
/data/gpfs/assoc/bch709/Course_material/RNASeq_trimmed_fastq
```

## Genome annotation preparation
```bash
/data/gpfs/assoc/bch709/Course_material/genome/bch709_assembly.fasta
/data/gpfs/assoc/bch709/Course_material/genome/bch709.gff
```

## GTF format
The Gene Transfer Format (GTF) is a widely used format for storing gene annotations. You can obtain GTF files easily from the UCSC table browser and Ensembl


Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
2. source - name of the program that generated this feature, or the data source (database or project name)
3. feature - feature type name, e.g. Gene, Variation, Similarity
4. start - Start position of the feature, with sequence numbering starting at 1.
5. end - End position of the feature, with sequence numbering starting at 1.
6. score - A floating point value.
7. strand - defined as + (forward) or - (reverse).
8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
9. attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

## GFF3 to GTF
```bash
gffread --help

gffread bch709.gff -T -o bch709.gtf
```

## STAR Aligner
We will align our reads to the reference genome using STAR (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments.

![star]({{site.baseurl}}/fig/STAR.png)

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs)

## Genome sequence preparation

Now that you have the genome and annotation files, you will create a genome index using the following script:
```bash
STAR --runMode genomeGenerate --genomeDir /data/gpfs/assoc/bch709/<YOURID>/genomernaseq/reference --genomeFastaFiles /data/gpfs/assoc/bch709/<YOURID>/genomernaseq/bch709_assembly.fasta --sjdbGTFfile  /data/gpfs/assoc/bch709/<YOURID>/genomernaseq/bch709.gtf --runThreadN 8 --genomeSAindexNbases 10
```


