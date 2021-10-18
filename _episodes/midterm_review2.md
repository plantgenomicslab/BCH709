---
layout: page
title: midterm review
published: true
---


## How does RNA sequencing allow the determination of a transcriptome for an entire genome? RNASeq provides a good approach for bioinformatics in one example. (min 300, max. 800 words) (20 points)


## What is the best sequencing platform for RNA-Seq? (Please describe Next generation and 3rd generation sequencing). Please justify your answer. (min 300, max. 800 words) (20 points)



## Describe the fastq and fasta format. Describe how fastq and fasta are used in one example each. Explain why fastq uses ASCII codes. (min 300, max. 800 words) (20 points)



### FASTA
FASTA stores a variable number of sequence records, and for each record it stores the sequence itself, and a sequence ID. Each record starts with a header line whose first character is >, followed by the sequence ID. The next lines of a record contain the actual sequence.

The Wikipedia artice gives several examples for peptide sequences, but since FASTQ and SAM are used exclusively (?) for nucleotide sequences, here’s a nucleotide example:

```
>Mus_musculus_tRNA-Ala-AGC-1-1 (chr13.trna34-AlaAGC)
GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGcCCTGGGTTCGATCC
CCAGCACCTCCA
>Mus_musculus_tRNA-Ala-AGC-10-1 (chr13.trna457-AlaAGC)
GGGGGATTAGCTCAAATGGTAGAGCGCTCGCTTAGCATGCAAGAGGtAGTGGGATCGATG
CCCACATCCTCCA
```

The ID can be in any arbitrary format, although several conventions exist.

In the context of nucleotide sequences, FASTA is mostly used to store reference data; that is, data extracted from a curated database; the above is adapted from GtRNAdb (a database of tRNA sequences).

### FASTQ
FASTQ was conceived to solve a specific problem of FASTA files: when sequencing, the confidence in a given base call (that is, the identity of a nucleotide) varies. This is expressed in the Phred quality score. FASTA had no standardised way of encoding this. By contrast, a FASTQ record contains a sequence of quality scores for each nucleotide.

A FASTQ record has the following format:

A line starting with @, containing the sequence ID.

One or more lines that contain the sequence.

A new line starting with the character +, and being either empty or repeating the sequence ID.

One or more lines that contain the quality scores.

Here’s an example of a FASTQ file with two records:
```
@071112_SLXA-EAS1_s_7:5:1:817:345
GGGTGATGGCCGCTGCCGATGGCGTC
AAATCCCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIII
IIII9IG9IC
@071112_SLXA-EAS1_s_7:5:1:801:338
GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBI
```
FASTQ files are mostly used to store short-read data from high-throughput sequencing experiments. As a consequence, the sequence and quality scores are usually put into a single line each, and indeed many tools assume that each record in a FASTQ file is exactly four lines long, even though this isn’t guaranteed.

In FASTQ files, quality scores are encoded into a compact form, which uses only 1 byte per quality value. In this encoding, the quality score is represented as the character with an ASCII code equal to its value + 33

As for FASTA, the format of the sequence ID isn’t standardised, but different producers of FASTQ use fixed notations that follow strict conventions.



### IUPAC code

[Johnson, 2021, An extended IUPAC nomenclature code for polymorphic nucleic acids](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/)

![IUPAC](https://www.researchgate.net/profile/Amarda-Shehu/publication/264056973/figure/fig12/AS:340683517906961@1458236688590/IUPAC-code-is-adapted-from-100.png)



## Please create a "midterm" Conda environment, activate and install the following software. After installation, please export installed packages in the "midterm" environment by 'conda env export'  and upload file here. (20 points)


```
conda create -n midterm
conda activate midterm
conda install -c bioconda  -c conda-forge trim-galore hisat2 subread seqkit
conda env export >> midterm.yaml
```




## Download the below file and describe the command that you use by answering the following questions.
```
ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
```

```bash
mkdir -p midterm/Q4

cd midterm/Q4

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
```

## Please check how many CDSs are in this file?
```bash
 zcat Homo_sapiens.GRCh38.cds.all.fa.gz  | egrep -c ">"
```
## How many bp in this file?
```bash
zcat Homo_sapiens.GRCh38.cds.all.fa.gz  | grep -v ">" | wc -c
```

## How long is this CDS? ENST00000570495.5
```bash
seqkit  grep -p ENST00000570495.5 Homo_sapiens.GRCh38.cds.all.fa.gz > ENST00000570495.fasta
grep  -v "^>" ENST00000570495.fasta | wc

```
```bash
seqkit fx2tab Homo_sapiens.GRCh38.cds.all.fa.gz > Homo_sapiens.GRCh38.cds.all.fa.tab
grep  ENST00000570495 Homo_sapiens.GRCh38.cds.all.fa.tab
grep  ENST00000570495 Homo_sapiens.GRCh38.cds.all.fa.tab | cut -f 2
grep  ENST00000570495 Homo_sapiens.GRCh38.cds.all.fa.tab | cut -f 2 | wc
 ```



## Fastqc might have an error, please skip fastqc.

We have paired end RNA-Seq fastq files for DT (drought treatment) and WT (wildtype) for certain plant. Please process these adaptors and quality trimming, mapping to  bch709.fasta by Hisat2 and quantify featurecount. 

Copy the your commands that you used for this question and paste it as text (Text entry) (please use history command)
How many reads were in the fastq files. (Text entry)
How was the mapping (alignment rate) for DT and WT? (Text entry)
How were the read counts on Chr4_gene_1834 in DT and WT samples? (Text entry)
Run MultiQC under rnaseq_test environment, and download & upload MultiQC results.  (File upload)

Download below file (GTF file for annotation)
https://www.dropbox.com/s/e9dvdkrl9dta4qg/bch709.gtf (Links to an external site.)

Copy below files (Reads and reference genome)
/data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/bch709.fasta

/data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/DT1_R1.fastq.gz

/data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/DT1_R2.fastq.gz

/data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/WT1_R1.fastq.gz

/data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/WT1_R2.fastq.gz


```bash
mkdir -p ~/midterm/Q6

cd ~/midterm/Q6

wget https://www.dropbox.com/s/e9dvdkrl9dta4qg/bch709.gtf

cp /data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/bch709.fasta .

cp /data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/DT1_R1.fastq.gz .

cp /data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/DT1_R2.fastq.gz .

cp /data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/WT1_R1.fastq.gz .

cp /data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq/WT1_R2.fastq.gz .
```



```bash
trim_galore  WT1_R1.fastq.gz WT1_R2.fastq.gz --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 2  --max_n 40  --gzip -o trim  --fastqc

trim_galore DT1_R1.fastq.gz DT1_R2.fastq.gz --paired   --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 2  --max_n 40  --gzip -o trim  --fastqc

ls -algh trim

```

```bash
hisat2-build bch709.fasta bch709

hisat2 -x bch709 --threads 2 -1 trim/WT1_R1_val_1.fq.gz  -2 trim/WT1_R2_val_2.fq.gz  -S align_WT1.sam 2> align_WT1.txt

samtools view -Sb align_WT1.sam > align_WT1.bam

hisat2 -x bch709 --threads 2 -1 trim/DT1_R1_val_1.fq.gz  -2 trim/DT1_R2_val_2.fq.gz  -S align_DT.sam 2> align_DT1.txt

samtools view -Sb align_DT1.sam > align_DT1.bam

```


```bash
featureCounts -p  -a bch709.gtf align_WT1.bam align_DT1.bam -o counts.txt
```


How many reads were in the fastq files. (Text entry)

```bash
seqkit stats  DT1_R1.fastq.gz DT1_R2.fastq.gz WT1_R1.fastq.gz WT1_R2.fastq.gz
```

How was the mapping (alignment rate) for DT and WT? (Text entry)
```bash
cat align_WT1.txt
cat align_DT1.txt
```

How were the read counts on Chr4_gene_1834 in DT and WT samples? (Text entry)

```bash
egrep "Chr4_gene_1834" counts.txt
```

Run MultiQC under rnaseq_test environment, and download & upload MultiQC results.  (File upload)

```bash
 conda install -c bioconda -c conda-forge multiqc
```


