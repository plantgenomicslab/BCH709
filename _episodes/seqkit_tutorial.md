---
layout: page
title: seqkit_tutorial
published: true
---


# Tutorial
## seqkit

### Create and Installation

```bash
conda create -n seqkit -c bioconda seqkit bedops csvtk
```

### Environment activate
```bash 
conda activate seqkit
```

### Export environment
```bash
conda env export > seqkit.yml
```

### Conda deactivate
```bash
conda deactivate
```

### Conda environment remove
```bash
conda remove --name seqkit --all
```

### Reinstall environment

```bash
conda env create -f seqkit.yml
```

### Environment activate
```bash 
conda activate seqkit
```

### use seqkit

```bash
seqkit --help
```

### make homework folder
```
mkdir /data/gpfs/assoc/bch709-4/${USER}
rm -rf ~/bch709
ln -s /data/gpfs/assoc/bch709-4/${USER} ~/bch709
cd ~/bch709
```

### Usage

```bash
SeqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation

Version: 2.0.0

Author: Wei Shen <shenwei356@gmail.com>

Documents  : http://bioinf.shenwei.me/seqkit
Source code: https://github.com/shenwei356/seqkit
Please cite: https://doi.org/10.1371/journal.pone.0163962

Available Commands:
  amplicon        extract amplicon (or specific region around it) via primer(s)
  bam             monitoring and online histograms of BAM record features
  common          find common sequences of multiple files by id/name/sequence
  completion      generate the autocompletion script for the specified shell
  concat          concatenate sequences with same ID from multiple files
  convert         convert FASTQ quality encoding between Sanger, Solexa and Illumina
  duplicate       duplicate sequences N times
  faidx           create FASTA index file and extract subsequence
  fish            look for short sequences in larger sequences using local alignment
  fq2fa           convert FASTQ to FASTA
  fx2tab          convert FASTA/Q to tabular format (and length, GC content, average quality...)
  genautocomplete generate shell autocompletion script (bash|zsh|fish|powershell)
  grep            search sequences by ID/name/sequence/sequence motifs, mismatch allowed
  head            print first N FASTA/Q records
  head-genome     print sequences of the first genome with common prefixes in name
  help            Help about any command
  locate          locate subsequences/motifs, mismatch allowed
  mutate          edit sequence (point mutation, insertion, deletion)
  pair            match up paired-end reads from two fastq files
  range           print FASTA/Q records in a range (start:end)
  rename          rename duplicated IDs
  replace         replace name/sequence by regular expression
  restart         reset start position for circular genome
  rmdup           remove duplicated sequences by ID/name/sequence
  sample          sample sequences by number or proportion
  sana            sanitize broken single line FASTQ files
  scat            real time recursive concatenation and streaming of fastx files
  seq             transform sequences (extract ID, filter by length, remove gaps...)
  shuffle         shuffle sequences
  sliding         extract subsequences in sliding windows
  sort            sort sequences by id/name/sequence/length
  split           split sequences into files by id/seq region/size/parts (mainly for FASTA)
  split2          split sequences into files by size/parts (FASTA, PE/SE FASTQ)
  stats           simple statistics of FASTA/Q files
  subseq          get subsequences by region/gtf/bed, including flanking sequences
  tab2fx          convert tabular format to FASTA/Q format
  translate       translate DNA/RNA to protein sequence (supporting ambiguous bases)
  version         print version information and check for update
  watch           monitoring and online histograms of sequence features

Flags:
      --alphabet-guess-seq-length int   length of sequence prefix of the first FASTA record based on which seqkit guesses the sequence type (0 for whole seq) (default 10000)
  -h, --help                            help for seqkit
      --id-ncbi                         FASTA head is NCBI-style, e.g. >gi|110645304|ref|NC_002516.2| Pseud...
      --id-regexp string                regular expression for parsing ID (default "^(\\S+)\\s?")
      --infile-list string              file of input files list (one file per line), if given, they are appended to files from cli arguments
  -w, --line-width int                  line width when outputing FASTA format (0 for no wrap) (default 60)
  -o, --out-file string                 out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
      --quiet                           be quiet and do not show extra information
  -t, --seq-type string                 sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
  -j, --threads int                     number of CPUs. can also set with environment variable SEQKIT_THREADS) (default 4)
```
### Current folder
```bash
cd ~/bch709
```


### Datasets

Test fastq file
- [`fastq`](https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/data/test.fastq.gz)

```bash
wget https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/data/test.fastq.gz
ll -tr
```


Datasets from [The miRBase Sequence Database -- Release 21](ftp://mirbase.org/pub/mirbase/21/)

- [`hairpin.fa`](https://mirbase.org/download/hairpin.fa)
```bash
wget https://mirbase.org/download/hairpin.fa
ll -tr
```
- [`mature.fa.gz`](https://mirbase.org/download/mature.fa)
```bash
wget https://mirbase.org/download/mature.fa
ll -tr
```
- [`miRNA.diff`](https://mirbase.org/download/miRNA.diff)
```bash
wget https://mirbase.org/download/miRNA.diff
ll -tr
```
Human genome from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/) (For `seqkit subseq`)

https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/


- [`GRCh38_latest_genomic.fna.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)
```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
ls -algh
```

- [`GRCh38_latest_genomic.gtf.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz)
```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz
ls -algh
```

- `GRCh38_latest_genomic.bed.gz` 
create with command
```bash
 zcat GRCh38_latest_genomic.gtf.gz | gtf2bed --do-not-sort | gzip -c > GRCh38_latest_genomic.bed.gz
```

Only DNA and gtf/bed data of Chr1 were used:

- `chr1.fa.gz`
```bash
seqkit grep -p NC_000001.11 GRCh38_latest_genomic.fna.gz -o chr1.fa.gz
ls -algh
```
- `chr1.gtf.gz`
```bash
zcat GRCh38_latest_genomic.gtf.gz | grep "NC_000001.11" | gzip -c > chr1.gtf.gz
ls -algh
```
- `chr1.bed.gz`
```bash
zcat GRCh38_latest_genomic.bed.gz| grep "NC_000001.11" | gzip -c > chr1.bed.gz
ls -algh
````

## seq

Usage

```bash
transform sequences (revserse, complement, extract ID...)

Usage:
  seqkit seq [flags]

Flags:
  -p, --complement                complement sequence (blank for Protein sequence)
      --dna2rna                   DNA to RNA
  -G, --gap-letter string         gap letters (default "- ")
  -l, --lower-case                print sequences in lower case
  -n, --name                      only print names
  -i, --only-id                   print ID instead of full head
  -q, --qual                      only print qualities
  -g, --remove-gaps               remove gaps
  -r, --reverse                   reverse sequence)
      --rna2dna                   RNA to DNA
  -s, --seq                       only print sequences
  -u, --upper-case                print sequences in upper case
  -v, --validate-seq              validate bases according to the alphabet
  -V, --validate-seq-length int   length of sequence to validate (0 for whole seq) (default 10000)

```

Examples

1. Read and print

- From file:
```bash
cat hairpin.fa
seqkit seq hairpin.fa
```

```bash
zcat test.fastq.gz
seqkit seq test.fastq.gz
```

- From stdin:
```bash
zcat hairpin.fa.gz | seqkit seq
```

2. Sequence types

- By default, `seqkit seq` automatically detect the sequence type
```bash
 echo -e ">seq\nacgtryswkmbdhvACGTRYSWKMBDHV" | seqkit stat
```
```bash
echo -e ">seq\nACGUN ACGUN" | seqkit stat
```
```bash
echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit stat
```

3. Only print names

- Full name:
```bash
seqkit seq hairpin.fa -n
```

- Only ID:
```bash
seqkit seq hairpin.fa -n -i
```

- Custom ID region by regular expression (this could be applied to all subcommands):
```bash
seqkit seq hairpin.fa -n -i --id-regexp "^[^\s]+\s([^\s]+)\s"
```

Please check regular expression in [Regex](https://regex101.com/)

4. Only print seq (global flag `-w` defines the output line width, 0 for no wrap)

```bash
seqkit seq hairpin.fa -s -w 0
```

5. Reverse comlement sequence
```bash
seqkit seq hairpin.fa.gz -r -p
```

6. Remove gaps and to lower/upper case
```bash
echo -e ">seq\nACGT-ACTGC-ACC" | seqkit seq -g -u
```


7. RNA to DNA
```bash
echo -e ">seq\nUCAUAUGCUUGUCUCAAAGAUUA" | seqkit seq --rna2dna
```


## subseq

Usage

```bash
get subsequences by region/gtf/bed, including flanking sequences.

Recommendation: use plain FASTA file, so seqkit could utilize FASTA index.

The definition of region is 1-based and with some custom design.

Examples:

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n

Usage:
  seqkit subseq [flags]

Flags:
      --bed string        by BED file
      --chr value         select limited sequence with sequence IDs (multiple value supported, case ignored) (default [])
  -d, --down-stream int   down stream length
      --feature value     select limited feature types (multiple value supported, case ignored, only works with GTF) (default [])
      --gtf string        by GTF (version 2.2) file
  -f, --only-flank        only return up/down stream sequence
  -r, --region string     by region. e.g 1:12 for first 12 bases, -12:-1 for last 12 bases, 13:-1 for cutting first 12 bases. type "seqkit subseq -h" for more examples
  -u, --up-stream int     up stream length

```

Examples

***Recommendation: use plain FASTA file, so seqkit could utilize FASTA index.***

1. First 12 bases
```bash
cat hairpin.fa | seqkit subseq -r 1:12
```
1. Last 12 bases
```bash
cat hairpin.fa | seqkit subseq -r -12:-1
```
1. Subsequences without first and last 12 bases
```bash
cat hairpin.fa | seqkit subseq -r 13:-13
```


1. Get subsequence by GTF file


#### Human genome example:
***AVOID loading all data from GRCh38_latest_genomic.gtf.gz, the uncompressed data are so big and may exhaust your RAM.***

We could specify chromesomes and features.

```bash
seqkit subseq --gtf GRCh38_latest_genomic.gtf.gz --chr NC_000001.11 --feature CDS  GRCh38_latest_genomic.fna.gz > chr1.gtf.cds.fa
```

```bash
seqkit stat chr1.gtf.cds.fa
```

5. Get subsequences by BED file.
```bash
seqkit subseq --bed GRCh38_latest_genomic.bed.gz --chr NC_000001.11 GRCh38_latest_genomic.fna.gz >  chr1.bed.gz.fa
```

We may need to remove duplicated sequences

```bash
seqkit subseq --bed GRCh38_latest_genomic.bed.gz --chr NC_000001.11  GRCh38_latest_genomic.fna.gz | seqkit rmdup > chr1.bed.rmdup.fa
```

```bash
seqkit stat chr1.gz.*.fa
```

## sliding

- sliding sequences, circular genome supported

Usage:
  seqkit sliding [flags]

Flags:
  -C, --circular-genome   circular genome
  -s, --step int        step size
  -W, --window int      window size

```

Examples

1. General use

```bash
echo -e ">seq\nACGTacgtNN" | seqkit sliding -s 3 -W 6
```


2. Circular genome
```bash
echo -e ">seq\nACGTacgtNN" | seqkit sliding -s 3 -W 6 -C
```

3. Generate GC content for ploting
```bash
cat hairpin.fa | seqkit fx2tab | head -n 1 | seqkit tab2fx | seqkit sliding -s 5 -W 30 | seqkit fx2tab -n -g
```

## stat

Usage
simple statistics of FASTA files

```
Usage:
  seqkit stat [flags]
```

Examples

1. General use
```bash
seqkit stat *.f*{a,q}.gz
```

## fq2fa

covert FASTQ to FASTA

Usage:
  seqkit fq2fa [flags]

```bash
seqkit fq2fa test.fastq.gz -o test_.fa.gz
zcat test_.fa.gz
zcat test.fastq.gz
```

## fx2tab & tab2fx

Usage (fx2tab)

```bash
covert FASTA/Q to tabular format, and provide various information,
like sequence length, GC content/GC skew.

Usage:
  seqkit fx2tab [flags]

Flags:
  -B, --base-content value   print base content. (case ignored, multiple values supported) e.g. -B AT -B N (default [])
  -g, --gc                   print GC content
  -G, --gc-skew              print GC-Skew
  -H, --header-line          print header line
  -l, --length               print sequence length
  -n, --name                 only print names (no sequences and qualities)
  -i, --only-id              print ID instead of full head

Usage (tab2fx)
covert tabular format (first two/three columns) to FASTA/Q format

Usage:
  seqkit tab2fx [flags]

Flags:
  -p, --comment-line-prefix value   comment line prefix (default [#,//])


```

Examples

1. Default output

```bash
seqkit fx2tab hairpin.fa| head -n 2
```

1. Print sequence length, GC content, and only print names (no sequences),
we could also print title line by flag `-T`.

```bash
seqkit fx2tab hairpin.fa -l -g -n -i -H | head -n 4 | csvtk -t -C '&' pretty
```

1. Use fx2tab and tab2fx in pipe
```bash
cat hairpin.fa | seqkit fx2tab | seqkit tab2fx
zcat test.fastq.gz | seqkit fx2tab | seqkit tab2fx
```

1. Sort sequences by length (use `seqkit sort -l`)
```bash
cat hairpin.fa | seqkit fx2tab -l | sort -t"`echo -e '\t'`" -n -k4,4 | seqkit tab2fx
```

```bash
seqkit sort -l hairpin.fa
```

1. Sorting or filtering by GC (or other base by -flag `-B`) content could also achieved in similar way.

1. Get first 1000 sequences

```bash
seqkit fx2tab hairpin.fa | head -n 1000 | seqkit tab2fx

seqkit fx2tab test.fastq.gz | head -n 1000 | seqkit tab2fx
```

***Extension***

After converting FASTA to tabular format with `seqkit fx2tab`,
it could be handled with CSV/TSV tools,
 e.g. [csvtk](https://github.com/shenwei356/csvtkt), a cross-platform, efficient and practical CSV/TSV toolkit

- `csvtk grep` could be used to filter sequences (similar with `seqkit grep`)
- `csvtk inter` computates intersection of multiple files. It could achieve similar function
as `seqkit common -n` along with shell.
- `csvtk join` joins multiple CSV/TSV files by multiple IDs.
- [csv_melt](https://github.com/shenwei356/datakit/blob/master/csv_melt)
provides melt function, could be used in preparation of data for ploting.




## grep

Usage

```
search sequences by pattern(s) of name or sequence motifs

Usage:
  seqkit grep [flags]

Flags:
  -n, --by-name               match by full name instead of just id
  -s, --by-seq                match by seq
  -d, --degenerate            pattern/motif contains degenerate base
      --delete-matched        delete matched pattern to speedup
  -i, --ignore-case           ignore case
  -v, --invert-match          invert the sense of matching, to select non-matching records
  -p, --pattern value         search pattern (multiple values supported) (default [])
  -f, --pattern-file string   pattern file
  -r, --use-regexp            patterns are regular expression

```

Examples

1. Extract human hairpins (i.e. sequences with name starting with `hsa`)
```bash
cat hairpin.fa | seqkit grep -r -p ^hsa
```

1. Remove human and mice hairpins.
```bash
cat hairpin.fa| seqkit grep -r -p ^hsa -p ^mmu -v
```

1. Extract new entries by information from miRNA.diff.gz

    1. Get IDs of new entries.
```bash
cat miRNA.diff | grep ^# -v | grep NEW | cut -f 2 > list
more  list ##q to out
```
2. Extract by ID list file
```bash
cat hairpin.fa | seqkit grep -f list > new.fa
```

1. Extract sequences starting with AGGCG
```bash
cat hairpin.fa | seqkit grep -s -r -i -p ^aggcg
```

1. Extract sequences with TTSAA (AgsI digest site) in SEQUENCE. Base S stands for C or G.
```bash
cat hairpin.fa | seqkit grep -s -d -i -p UUSAA
```
It's equal to but simpler than:
```bash
cat hairpin.fa| seqkit grep -s -r -i -p UU[CG]AA
```

## locate

Usage: locate subsequences/motifs
```
Motifs could be EITHER plain sequence containing "ACTGN" OR regular
expression like "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" for ORFs.
Degenerate bases like "RYMM.." are also supported by flag -d.

By default, motifs are treated as regular expression.
When flag -d given, regular expression may be wrong.
For example: "\w" will be wrongly converted to "\[AT]".

Usage:
  seqkit locate [flags]

Flags:
  -d, --degenerate                pattern/motif contains degenerate base
  -i, --ignore-case               ignore case
  -P, --only-positive-strand      only search at positive strand
  -p, --pattern value             search pattern/motif (multiple values supported) (default [])
  -f, --pattern-file string       pattern/motif file (FASTA format)
  -V, --validate-seq-length int   length of sequence to validate (0 for whole seq) (default 10000)

```

Examples

1. Locate ORFs.
```bash
cat hairpin.fa | seqkit locate -i -r -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)"
```

1. Locate Motif.
```bash
cat hairpin.fa | seqkit locate -i -d -p AUGGACUN
```

***Notice that `seqkit grep` only searches in positive strand, but `seqkit loate` could recognize both strand***


## rmdup

Usage

```
remove duplicated sequences by id/name/sequence

Usage:
  seqkit rmdup [flags]

Flags:
    -n, --by-name                by full name instead of just id
    -s, --by-seq                 by seq
    -D, --dup-num-file string    file to save number and list of duplicated seqs
    -d, --dup-seqs-file string   file to save duplicated seqs
    -i, --ignore-case            ignore case
    -m, --md5                    use MD5 instead of original seqs to reduce memory usage when comparing by seqs

```

Examples

Similar to `common`.

1. General use
```bash
cat hairpin.fa | seqkit rmdup -s -o clean.fa.gz
zcat test.fastq.gz | seqkit rmdup -s -o clean.fa.gz
```
1. Save duplicated sequences to file
```bash
cat hairpin.fa | seqkit rmdup -s -i -o clean.fa.gz -d duplicated.fa.gz -D duplicated.detail.txt
cat duplicated.detail.txt 
```

## common

Usage

```
find common sequences of multiple files by id/name/sequence

Usage:
  seqkit common [flags]

Flags:
    -n, --by-name       match by full name instead of just id
    -s, --by-seq        match by sequence
    -i, --ignore-case   ignore case
    -m, --md5           use MD5 instead of original seqs to reduce memory usage when comparing by seqs

```

Examples
1. Create file
```bash
echo -e ">1234\nATGGTATTAGATA" >> file1.fa
echo -e ">2234\nATGGTATTTGATA" >> file1.fa
echo -e ">3234_A\nATGGTATTCGATA" >> file1.fa
echo -e ">1234_B\nATGGTATTCGATA" >> file2.fa
echo -e ">4234\nATGGTATTAGATA" >> file2.fa
echo -e ">2234\nATGGTATTGGGATA" >> file2.fa
cat file1.fa
cat file2.fa


1. By ID (default)
```bash
seqkit common file*.fa -o common.fasta
```
1. By full name
```bash
seqkit common file*.fa -n -o common.fasta
```

1. By sequence
```bash
seqkit common file*.fa -s -i -o common.fasta
```

## split

Usage

```
split sequences into files by name ID, subsequence of given region,
part size or number of parts.

The definition of region is 1-based and with some custom design.

Examples:

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n

Usage:
  seqkit split [flags]

Flags:
Flags:
  -i, --by-id              split squences according to sequence ID
  -p, --by-part int        split squences into N parts
  -r, --by-region string   split squences according to subsequence of given region. e.g 1:12 for first 12 bases, -12:-1 for last 12 bases. type "seqkit split -h" for more examples
  -s, --by-size int        split squences into multi parts with N sequences
  -d, --dry-run            dry run, just print message and no files will be created.
  -f, --force              overwrite output directory
  -k, --keep-temp          keep tempory FASTA and .fai file when using 2-pass mode
  -m, --md5                use MD5 instead of region sequence in output file when using flag -r (--by-region)
  -O, --out-dir string     output directory (default value is infile.split)
  -2, --two-pass           two-pass mode read files twice to lower memory usage. (only for FASTA format)

```

Examples

1. Split sequences into parts with at most 10000 sequences

```bash
seqkit split hairpin.fa -s 10000
```

1. Split sequences into 4 parts

```bash
seqkit split hairpin.fa -p 4
```
    ***To reduce memory usage when spliting big file, we should alwasy use flag `--two-pass`***

```bash
seqkit split hairpin.fa -p 4 -2
```


1. Split sequences by species. i.e. by custom IDs (first three letters)
```bash
seqkit split hairpin.fa -i --id-regexp "^([\w]+)\-" -2
```

1. Split sequences by sequence region (for example, sequence barcode)
```bash
seqkit split hairpin.fa -r 1:3 -2
```
        [INFO] split by region: 1:3
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] read sequence IDs and sequence region from FASTA file ...
        [INFO] create and read FASTA index ...
        [INFO] write 463 sequences to file: hairpin.region_1:3_AUG.fa.gz
        [INFO] write 349 sequences to file: hairpin.region_1:3_ACU.fa.gz
        [INFO] write 311 sequences to file: hairpin.region_1:3_CGG.fa.gz

    Sequence suffix could be defined as `-r -12:-1`

## sample

Usage

```
sample sequences by number or proportion.

Usage:
  seqkit sample [flags]

Flags:
  -n, --number int         sample by number (result may not exactly match)
  -p, --proportion float   sample by proportion
  -s, --rand-seed int      rand seed for shuffle (default 11)
  -2, --two-pass           2-pass mode read files twice to lower memory usage. Not allowed when reading from stdin

```

Examples

1. Sample by proportion
```bash
cat hairpin.fa | seqkit sample -p 0.1 -o sample.fa.gz
```
1. Sample by number
```bash
cat hairpin.fa | seqkit sample -n 1000 -o sample.fa.gz
```
    ***To reduce memory usage when spliting big file, we could use flag `--two-pass`***

    ***We can also use `seqkit sample -p` followed with `seqkit head -n`:***
```bash
cat hairpin.fa | seqkit sample -p 0.1 | seqkit head -n 1000 -o sample.fa.gz
```

1. Set rand seed to reproduce the result
```
cat hairpin.fa | seqkit sample -p 0.1 -s 11
```

1. Most of the time, we could shuffle after sampling
```bash
cat hairpin.fa | seqkit sample -p 0.1 | seqkit shuffle -o sample.fa.gz
```

Note that when sampling on FASTQ files, make sure using same random seed by
flag `-s` (`--rand-seed`)

## head

Usage

```
print first N FASTA/Q records

Usage:
  seqkit head [flags]

Flags:
  -n, --number int   print first N FASTA/Q records (default 10)

```

Examples

1. FASTA
```bash
seqkit head -n 1 hairpin.fa
```
        >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
        UAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA

1. FASTQ
```bash
seqkit head -n 1 test.fastq.gz
```


## replace

Usage replace name/sequence/by regular expression.

Note that the replacement supports capture variables.
e.g. $1 represents the text of the first submatch.
ATTENTION: use SINGLE quote NOT double quotes in *nix OS.

Examples: Adding space to all bases.
```bash
seqkit head -n 1 hairpin.fa | seqkit replace -p "(.)" -r '$1 ' -s
```

more on: http://shenwei356.github.io/seqkit/usage/#replace

Usage:
  seqkit replace [flags]

Flags:
  -s, --by-seq               replace seq
  -i, --ignore-case          ignore case
  -p, --pattern string       search regular expression
  -r, --replacement string   replacement. supporting capture variables.  e.g. $1 represents the text of the first submatch. ATTENTION: use SINGLE quote NOT double quotes in *nix OS or use the \ escape character. record number is also supported by "{NR}"

```

Examples

1. Remove descriptions
```bash
echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p " .+"
```

1. Replace "-" with "="
```bash
echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p "\-" -r '='
```

1. Remove gaps in sequences.
```bash
echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p " |-" -s
```

1. Add space to every base. **ATTENTION: use SINGLE quote NOT double quotes in *nix OS**
```bash
echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p "(.)" -r '$1 ' -s
```

1. Transpose sequence with [csvtk](https://github.com/shenwei356/csvtk)
```bash
echo -e ">seq1\nACTGACGT\n>seq2\nactgccgt" | seqkit replace -p "(.)" -r     "\$1 " -s | seqkit seq -s -u | csvtk space2tab | csvtk -t transpose
```


1. Rename with number of record
```bash
echo -e ">abc\nACTG\n>123\nATTT" |  seqkit replace -p .+ -r "seq_{NR}"
```

## shuffle

Usage

```
shuffle sequences.

By default, all records will be readed into memory.
For FASTA format, use flag -2 (--two-pass) to reduce memory usage. FASTQ not
supported.

Firstly, seqkit reads the sequence IDs. If the file is not plain FASTA file,
seqkit will write the sequences to tempory files, and create FASTA index.

Secondly, seqkit shuffles sequence IDs and extract sequences by FASTA index.

Usage:
  seqkit shuffle [flags]

Flags:
  -k, --keep-temp       keep tempory FASTA and .fai file when using 2-pass mode
  -s, --rand-seed int   rand seed for shuffle (default 23)
  -2, --two-pass        two-pass mode read files twice to lower memory usage. (only for FASTA format)

```

Examples

1. General use.
```bash
seqkit shuffle hairpin.fa > shuffled.fa
```

1. ***For big genome, you'd better use two-pass mode*** so seqkit could use
   FASTA index to reduce memory usage
```bash
time seqkit shuffle -2 GRCh38_latest_genomic.fna.gz > shuffle.fa
```

## sort

Usage

```
sort sequences by id/name/sequence/length.

By default, all records will be readed into memory.
For FASTA format, use flag -2 (--two-pass) to reduce memory usage. FASTQ not
supported.

Firstly, seqkit reads the sequence head and length information.
If the file is not plain FASTA file,
seqkit will write the sequences to tempory files, and create FASTA index.

Secondly, seqkit sort sequence by head and length information
and extract sequences by FASTA index.

Usage:
  seqkit sort [flags]

Flags:
  -l, --by-length               by sequence length
  -n, --by-name                 by full name instead of just id
  -s, --by-seq                  by sequence
  -i, --ignore-case             ignore case
  -k, --keep-temp               keep tempory FASTA and .fai file when using 2-pass mode
  -r, --reverse                 reverse the result
  -L, --seq-prefix-length int   length of sequence prefix on which seqkit sorts by sequences (0 for whole sequence) (default 10000)
  -2, --two-pass                two-pass mode read files twice to lower memory usage. (only for FASTA format)

```

Examples

***For FASTA format, use flag -2 (--two-pass) to reduce memory usage***

1. sort by ID
```bash
echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet
```


1. sort by ID, ignoring case.
```bash
echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -i
```

1. sort by seq, ignoring case.
```bash
echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -s -i
```

1. sort by sequence length
```bash
echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAAnnn\n>seq3\nacgt" | seqkit sort --quiet -l
```

### Quick glance

1. Sequence number
```bash
seqkit stat hairpin.fa
```

1. First 10 bases
```bash
cat hairpin.fa | seqkit subseq -r 1:10 | seqkit sort -s | seqkit seq -s | head -n 10

### Repeated hairpin sequences

We may want to check how may identical hairpins among different species there are.
`seqkit rmdup` could remove duplicated sequences by sequence content,
and save the replicates to another file (here is `duplicated.fa.gz`),
as well as replicating details (`duplicated.detail.txt`,
1th column is the repeated number,
2nd column contains sequence IDs seperated by comma).
```bash
seqkit rmdup -s -i hairpin.fa -o clean.fa.gz -d duplicated.fa.gz -D duplicated.detail.txt
```

The result shows the most conserved miRNAs among different species,
`mir-29b`, `mir-125`, `mir-19b-1` and `mir-20a`.
And the `dre-miR-430c` has the most multicopies in *Danio rerio*.

### Hairpins in different species

1. Before spliting by species, let's take a look at the sequence names.
```bash
seqkit seq hairpin.fa.gz -n | head -n 3
```

The first three letters (e.g. `cel`) are the abbreviation of species names. So we could split hairpins by the first letters by defining custom sequence ID parsing regular expression `^([\w]+)\-`.

By default, `seqkit` takes the first non-space letters as sequence ID.

For example,

    |   FASTA head                                                  |     ID                                            |
    |:--------------------------------------------------------------|:--------------------------------------------------|
    | >123456 gene name                                             | 123456                                            |
    | >longname                                                     | longname                                          |
    | >gi&#124;110645304&#124;ref&#124;NC_002516.2&#124; Pseudomona | gi&#124;110645304&#124;ref&#124;NC_002516.2&#124; |

    But for some sequences from NCBI,
    e.g. `>gi|110645304|ref|NC_002516.2| Pseudomona`, the ID is `NC_002516.2`.
    In this case, we could set sequence ID parsing regular expression by flag
    `--id-regexp "\|([^\|]+)\| "` or just use flag `--id-ncbi`. If you want
    the `gi` number, then use `--id-regexp "^gi\|([^\|]+)\|"`.

1. Split sequences by species.
A custom ID parsing regular expression is used, `^([\w]+)\-`.
```bash
seqkit split hairpin.fa -i --id-regexp "^([\w]+)\-" --two-pass
```
    ***To reduce memory usage when splitting big file, we should always use flag `--two-pass`***

2. Species with most miRNA hairpins. Third column is the sequences number.
```bash
cd hairpin.fa.split/;
seqkit stat hairpin.part_* | csvtk space2tab | csvtk -t sort -k num_seqs:nr | csvtk -t pretty| more
 
Here, a CSV/TSV tool [csvtk](https://github.com/shenwei356/csvtk) is used to sort and view the result.


## Softwares

1. [seqkit](https://github.com/shenwei356/seqkit). (Go).
   Version [v2](https://github.com/shenwei356/seqkit/releases/tag/v0.3.1.1).
   Compiled with Go 1.7rc5.
1. [fasta_utilities](https://github.com/jimhester/fasta_utilities). (Perl).
   Version [3dcc0bc](https://github.com/jimhester/fasta_utilities/tree/3dcc0bc6bf1e97839476221c26984b1789482579).
   Lots of dependencies to install.
1. [fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). (Perl).
   Version [0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2).
   Can't handle multi-line FASTA files.
1. [seqmagick](http://seqmagick.readthedocs.io/en/latest/index.html#installation). (Python).
   Version 0.6.1
1. [seqtk](https://github.com/lh3/seqtk). (C).
   Version [1.1-r92-dirty](https://github.com/lh3/seqtk/tree/fb85aad4ce1fc7b3d4543623418a1ae88fe1cea6).

Not used:

1. [pyfaidx](https://github.com/mdshw5/pyfaidx). (Python).
   Version [0.4.7.1](https://pypi.python.org/packages/source/p/pyfaidx/pyfaidx-0.4.7.1.tar.gz#md5=f33604a3550c2fa115ac7d33b952127d). *Not used, because it exhausted my memory (10G) when computing reverse-complement on a 5GB fasta file of 250 bp.*

A Python script [memusg](https://github.com/shenwei356/memusg) was used to compute running time and peak memory usage of a process.

## Features

Categories          |Features               |seqkit  |fasta_utilities|fastx_toolkit|pyfaidx|seqmagick|seqtk
:-------------------|:----------------------|:------:|:-------------:|:-----------:|:-----:|:-------:|:---:
**Formats support** |Multi-line FASTA       |Yes     |Yes            |--           |Yes    |Yes      |Yes
                    |FASTQ                  |Yes     |Yes            |Yes          |--     |Yes      |Yes
                    |Multi-line  FASTQ      |Yes     |Yes            |--           |--     |Yes      |Yes
                    |Validating sequences   |Yes     |--             |Yes          |Yes    |--       |--
                    |Supporting RNA         |Yes     |Yes            |--           |--     |Yes      |Yes
**Functions**       |Searching by motifs    |Yes     |Yes            |--           |--     |Yes      |--
                    |Sampling               |Yes     |--             |--           |--     |Yes      |Yes
                    |Extracting sub-sequence|Yes     |Yes            |--           |Yes    |Yes      |Yes
                    |Removing duplicates    |Yes     |--             |--           |--     |Partly   |--
                    |Splitting              |Yes     |Yes            |--           |Partly |--       |--
                    |Splitting by seq       |Yes     |--             |Yes          |Yes    |--       |--
                    |Shuffling              |Yes     |--             |--           |--     |--       |--
                    |Sorting                |Yes     |Yes            |--           |--     |Yes      |--
                    |Locating motifs        |Yes     |--             |--           |--     |--       |--
                    |Common sequences       |Yes     |--             |--           |--     |--       |--
                    |Cleaning bases         |Yes     |Yes            |Yes          |Yes    |--       |--
                    |Transcription          |Yes     |Yes            |Yes          |Yes    |Yes      |Yes
                    |Translation            |--      |Yes            |Yes          |Yes    |Yes      |--
                    |Filtering by size      |Indirect|Yes            |--           |Yes    |Yes      |--
                    |Renaming header        |Yes     |Yes            |--           |--     |Yes      |Yes
**Other features**  |Cross-platform         |Yes     |Partly         |Partly       |Yes    |Yes      |Yes
                    |Reading STDIN          |Yes     |Yes            |Yes          |--     |Yes      |Yes
                    |Reading gzipped file   |Yes     |Yes            |--           |--     |Yes      |Yes
                    |Writing gzip file      |Yes     |--             |--           |--     |Yes      |--

**Note 2**: See [usage](http://shenwei356.github.io/seqkit/usage/) for detailed options of seqkit.


