---
layout: page
title: seqkit_tutorial
published: true
---


# Tutorial
## seqkit

#### Installation

```bash
conda install -c bioconda seqkit bedops csvtk
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

### Datasets

Test fastq file
- [`fastq`](https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/data/test.fastq.gz)

Datasets from [The miRBase Sequence Database -- Release 21](ftp://mirbase.org/pub/mirbase/21/)

- [`hairpin.fa.gz`](ftp://mirbase.org/pub/mirbase/21/hairpin.fa.gz)
- [`mature.fa.gz`](ftp://mirbase.org/pub/mirbase/21/mature.fa.gz)
- [`miRNA.diff.gz`](ftp://mirbase.org/pub/mirbase/21/miRNA.diff.gz)

Human genome from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/) (For `seqkit subseq`)

https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/


- [`GRCh38_latest_genomic.fna.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)
- [`GRCh38_latest_genomic.gtf.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz)
- `GRCh38_latest_genomic.bed.gz` is converted from `GRCh38_latest_genomic.gtf.gz` by [`gtf2bed`](http://bedops.readthedocs.org/en/latest/content/reference/file-management/conversion/gtf2bed.html?highlight=gtf2bed)
with command
```bash
        zcat GRCh38_latest_genomic.gtf.gz | gtf2bed --do-not-sort | gzip -c > GRCh38_latest_genomic.bed.gz
```

Only DNA and gtf/bed data of Chr1 were used:

- `chr1.fa.gz`
```bash
seqkit grep -p NC_000001.11 GRCh38_latest_genomic.fna.gz -o chr1.fa.gz
```
- `chr1.gtf.gz`
```bash
zcat GRCh38_latest_genomic.gtf.gz | grep "NC_000001.11" | gzip -c > chr1.gtf.gz
```
- `chr1.bed.gz`
```bash
zcat GRCh38_latest_genomic.bed.gz| grep "NC_000001.11" | gzip -c > chr1.bed.gz
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

seqkit seq hairpin.fa.gz
            >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
            UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
            UAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
```
```bash

seqkit seq test.fastq.gz
            @HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG
            TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTG
            CCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGT
            ATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGA
            TCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT
            +
            HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHH
            HHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIII
            HIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHG
            IHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH
```
    - From stdin:

            zcat hairpin.fa.gz | seqkit seq


1. Sequence types

    - By default, `seqkit seq` automatically detect the sequence type

            $ echo -e ">seq\nacgtryswkmbdhvACGTRYSWKMBDHV" | seqkit stat
            file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTA   DNA          1       28       28       28       28

            $ echo -e ">seq\nACGUN ACGUN" | seqkit stat
            file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTA   RNA          1       11       11       11       11

            $ echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit stat
            file  format  type     num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTA   Protein         1       23       23       23       23

            $ echo -e "@read\nACTGCN\n+\n@IICCG" | seqkit stat
            file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTQ   DNA          1        6        6        6        6

    - You can also set sequence type by flag `-t` (`--seq-type`).
      But this only take effect on subcommands `seq` and `locate`.

            $ echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit seq -t dna
            [INFO] when flag -t (--seq-type) given, flag -v (--validate-seq) is automatically switched on
            [ERRO] error when parsing seq: seq (invalid DNAredundant letter: e)


1. Only print names

    - Full name:

            $ seqkit seq hairpin.fa.gz -n
            cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
            cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop
            cel-mir-1 MI0000003 Caenorhabditis elegans miR-1 stem-loop

    - Only ID:

            $ seqkit seq hairpin.fa.gz -n -i
            cel-let-7
            cel-lin-4
            cel-mir-1

    - Custom ID region by regular expression (this could be applied to all subcommands):

            $ seqkit seq hairpin.fa.gz -n -i --id-regexp "^[^\s]+\s([^\s]+)\s"
            MI0000001
            MI0000002
            MI0000003

Please check regular expression in [Regex](https://regex101.com/)

1. Only print seq (global flag `-w` defines the output line width, 0 for no wrap)

        $ seqkit seq hairpin.fa.gz -s -w 0
        UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAACUAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
        AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGAGUGUACUAUUGAUGCUUCACACCUGGGCUCUCCGGGUACCAGGACGGUUUGAGCAGAU
        AAAGUGACCGUACCGAGCUGCAUACUUCCUUACAUGCCCAUACUAUAUCAUAAAUGGAUAUGGAAUGUAAAGAAGUAUGUAGAACGGGGUGGUAGU

1. Reverse comlement sequence

        $ seqkit seq hairpin.fa.gz -r -p
        >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        UCGAAGAGUUCUGUCUCCGGUAAGGUAGAAAAUUGCAUAGUUCACCGGUGGUAAUAUUCC
        AAACUAUACAACCUACUACCUCACCGGAUCCACAGUGUA

1. Remove gaps and to lower/upper case

        $ echo -e ">seq\nACGT-ACTGC-ACC" | seqkit seq -g -u
        >seq
        ACGTACTGCACC

1. RNA to DNA

        $ echo -e ">seq\nUCAUAUGCUUGUCUCAAAGAUUA" | seqkit seq --rna2dna
        >seq
        TCATATGCTTGTCTCAAAGATTA


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

        $ zcat hairpin.fa.gz | seqkit subseq -r 1:12

1. Last 12 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r -12:-1

1. Subsequences without first and last 12 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r 13:-13



1. Get subsequence by GTF file

        $ cat t.fa
        >seq
        actgACTGactgn
        $ cat t.gtf
        seq     test    CDS     5       8       .       .       .       gene_id "A"; transcript_id "";
        seq     test    CDS     5       8       .       -       .       gene_id "B"; transcript_id "";
        $ seqkit

        $ seqkit subseq --gtf t.gtf t.fa
        >seq_5:8:. A
        ACTG
        >seq_5:8:- B
        CAGT

    Human genome example:

    ***AVOID loading all data from GRCh38_latest_genomic.gtf.gz,
    the uncompressed data are so big and may exhaust your RAM.***

    We could specify chromesomes and features.

        $ seqkit subseq --gtf GRCh38_latest_genomic.gtf.gz --chr NC_000001.11 --feature cds  GRCh38_latest_genomic.fna.gz > chr1.gtf.cds.fa

        $ seqkit stat chr1.gtf.cds.fa
        file             format  type  num_seqs    sum_len  min_len  avg_len  max_len
        chr1.gtf.cds.fa  FASTA   DNA     65,012  9,842,274        1    151.4   12,045

1. Get CDS and 3bp up-stream sequences

        $ seqkit subseq --gtf t.gtf t.fa -u 3
        >seq_5:8:._us:3 A
        ctgACTG
        >seq_5:8:-_us:3 B
        agtCAGT

1. Get 3bp up-stream sequences of CDS, not including CDS

        $ seqkit subseq --gtf t.gtf t.fa -u 3 -f
        >seq_5:8:._usf:3 A
        ctg
        >seq_5:8:-_usf:3 B
        agt

1. Get subsequences by BED file.

    ***AVOID loading all data from GRCh38_latest_genomic.gtf.gz,
    the uncompressed data are so big and may exhaust your RAM.***

        $  seqkit subseq --bed GRCh38_latest_genomic.bed.gz --chr NC_000001.11 GRCh38_latest_genomic.fna.gz >  chr1.bed.gz.fa

    We may need to remove duplicated sequences

        $ seqkit subseq --bed GRCh38_latest_genomic.bed.gz --chr NC_000001.11  GRCh38_latest_genomic.fna.gz | seqkit rmdup > chr1.bed.rmdup.fa
        [INFO] 141060 duplicated records removed

    Summary:

        $ seqkit stat chr1.gz.*.fa
```bash

file               format  type  num_seqs        sum_len  min_len   avg_len    max_len
chr1.bed.gz.fa     FASTA   DNA    372,161  1,303,704,051        1   3,503.1  1,551,949
chr1.bed.rmdup.fa  FASTA   DNA     63,548    674,967,858        1  10,621.4  1,551,949
```

## sliding

Usage

```
sliding sequences, circular genome supported

Usage:
  seqkit sliding [flags]

Flags:
  -C, --circular-genome   circular genome
  -s, --step int        step size
  -W, --window int      window size

```

Examples

1. General use

        $ echo -e ">seq\nACGTacgtNN" | seqkit sliding -s 3 -W 6
        >seq_sliding:1-6
        ACGTac
        >seq_sliding:4-9
        TacgtN

2. Circular genome

        $ echo -e ">seq\nACGTacgtNN" | seqkit sliding -s 3 -W 6 -C
        >seq_sliding:1-6
        ACGTac
        >seq_sliding:4-9
        TacgtN
        >seq_sliding:7-2
        gtNNAC
        >seq_sliding:10-5
        NACGTa

3. Generate GC content for ploting

        $ zcat hairpin.fa.gz | seqkit fx2tab | head -n 1 | seqkit tab2fx | seqkit sliding -s 5 -W 30 | seqkit fx2tab -n -g
        cel-let-7_sliding:1-30          50.00
        cel-let-7_sliding:6-35          46.67
        cel-let-7_sliding:11-40         43.33
        cel-let-7_sliding:16-45         36.67
        cel-let-7_sliding:21-50         33.33
        cel-let-7_sliding:26-55         40.00
        ...

## stat

Usage

```
simple statistics of FASTA files

Usage:
  seqkit stat [flags]

```

Eexamples

1. General use

        $ seqkit stat *.f{a,q}.gz
        file           format  type  num_seqs    sum_len  min_len  avg_len  max_len
        hairpin.fa.gz  FASTA   RNA     28,645  2,949,871       39      103    2,354
        mature.fa.gz   FASTA   RNA     35,828    781,222       15     21.8       34
        test.fastq.gz  FASTQ   DNA      2,500    567,516      226      227      229
        reads_2.fq.gz  FASTQ   DNA      2,500    560,002      223      224      225

## fq2fa

Usage

```
covert FASTQ to FASTA

Usage:
  seqkit fq2fa [flags]

```

Examples

    seqkit fq2fa test.fastq.gz -o test_.fa.gz


## fx2tab & tab2fx

Usage (fx2tab)

```
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

```

Usage (tab2fx)

```
covert tabular format (first two/three columns) to FASTA/Q format

Usage:
  seqkit tab2fx [flags]

Flags:
  -p, --comment-line-prefix value   comment line prefix (default [#,//])


```

Examples

1. Default output

        $ seqkit fx2tab hairpin.fa.gz | head -n 2
        cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop      UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAACUAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
        cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop      AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGAGUGUACUAUUGAUGCUUCACACCUGGGCUCUCCGGGUACCAGGACGGUUUGAGCAGAU


1. Print sequence length, GC content, and only print names (no sequences),
we could also print title line by flag `-T`.

        $ seqkit fx2tab hairpin.fa.gz -l -g -n -i -H | head -n 4 | csvtk -t -C '&' pretty
        #name       seq   qual   length   GC
        cel-let-7                99       43.43
        cel-lin-4                94       54.26
        cel-mir-1                96       40.62

1. Use fx2tab and tab2fx in pipe

        $ zcat hairpin.fa.gz | seqkit fx2tab | seqkit tab2fx

        $ zcat test.fastq.gz | seqkit fx2tab | seqkit tab2fx

1. Sort sequences by length (use `seqkit sort -l`)

        $ zcat hairpin.fa.gz | seqkit fx2tab -l | sort -t"`echo -e '\t'`" -n -k4,4 | seqkit tab2fx
        >cin-mir-4129 MI0015684 Ciona intestinalis miR-4129 stem-loop
        UUCGUUAUUGGAAGACCUUAGUCCGUUAAUAAAGGCAUC
        >mmu-mir-7228 MI0023723 Mus musculus miR-7228 stem-loop
        UGGCGACCUGAACAGAUGUCGCAGUGUUCGGUCUCCAGU
        >cin-mir-4103 MI0015657 Ciona intestinalis miR-4103 stem-loop
        ACCACGGGUCUGUGACGUAGCAGCGCUGCGGGUCCGCUGU

        $ seqkit sort -l hairpin.fa.gz

    Sorting or filtering by GC (or other base by -flag `-B`) content could also achieved in similar way.

1. Get first 1000 sequences

        $ seqkit fx2tab hairpin.fa.gz | head -n 1000 | seqkit tab2fx

        $ seqkit fx2tab test.fastq.gz | head -n 1000 | seqkit tab2fx

**Extension**

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

        $ zcat hairpin.fa.gz | seqkit grep -r -p ^hsa
        >hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop
        UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAU
        ACAAUCUACUGUCUUUCCUA
        >hsa-let-7a-2 MI0000061 Homo sapiens let-7a-2 stem-loop
        AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACAUCAAGGGAGAUAACUGUACAGCCU
        CCUAGCUUUCCU

1. Remove human and mice hairpins.

        $ zcat hairpin.fa.gz | seqkit grep -r -p ^hsa -p ^mmu -v

1. Extract new entries by information from miRNA.diff.gz

    1. Get IDs of new entries.

            $ zcat miRNA.diff.gz | grep ^# -v | grep NEW | cut -f 2 > list
            $ more list
            cfa-mir-486
            cfa-mir-339-1
            pmi-let-7


    2. Extract by ID list file

            $ zcat hairpin.fa.gz | seqkit grep -f list > new.fa

1. Extract sequences starting with AGGCG

        $ zcat hairpin.fa.gz | seqkit grep -s -r -i -p ^aggcg

1. Extract sequences with TTSAA (AgsI digest site) in SEQUENCE. Base S stands for C or G.

        $ zcat hairpin.fa.gz | seqkit grep -s -d -i -p UUSAA

    It's equal to but simpler than:

        $ zcat hairpin.fa.gz | seqkit grep -s -r -i -p UU[CG]AA


## locate

Usage

```
locate subsequences/motifs

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

        $ zcat hairpin.fa.gz | seqkit locate -i -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)"
        seqID   patternName     pattern strand  start   end     matched
        cel-lin-4       A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        +  136      AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGA
        cel-mir-1       A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        +  54       95      AUGGAUAUGGAAUGUAAAGAAGUAUGUAGAACGGGGUGGUAG
        cel-mir-1       A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        -  43       51      AUGAUAUAG

1. Locate Motif.

        $ zcat hairpin.fa.gz | seqkit locate -i -d -p AUGGACUN
        seqID         patternName   pattern    strand   start   end   matched
        cel-mir-58a   AUGGACUN      AUGGACUN   +        81      88    AUGGACUG
        ath-MIR163    AUGGACUN      AUGGACUN   -        122     129   AUGGACUC

    Notice that `seqkit grep` only searches in positive strand, but `seqkit loate` could recognize both strand


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

        $ zcat hairpin.fa.gz | seqkit rmdup -s -o clean.fa.gz
        [INFO] 2226 duplicated records removed

        $ zcat test.fastq.gz | seqkit rmdup -s -o clean.fa.gz
        [INFO] 0 duplicated records removed

1. Save duplicated sequences to file

        $ zcat hairpin.fa.gz | seqkit rmdup -s -i -m -o clean.fa.gz -d duplicated.fa.gz -D duplicated.detail.txt

        $ cat duplicated.detail.txt   # here is not the entire list
        3       hsa-mir-424, mml-mir-424, ppy-mir-424
        3       hsa-mir-342, mml-mir-342, ppy-mir-342
        2       ngi-mir-932, nlo-mir-932
        2       ssc-mir-9784-1, ssc-mir-9784-2

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

1. By ID (default)

        seqkit common file*.fa -o common.fasta

1. By full name

        seqkit common file*.fa -n -o common.fasta

1. By sequence

        seqkit common file*.fa -s -i -o common.fasta

1. By sequence (***for large sequences***)

        seqkit common file*.fa -s -i -o common.fasta --md5


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

        $ seqkit split hairpin.fa.gz -s 10000
        [INFO] split into 10000 seqs per file
        [INFO] write 10000 sequences to file: hairpin.fa.part_001.gz
        [INFO] write 10000 sequences to file: hairpin.fa.part_002.gz
        [INFO] write 8645 sequences to file: hairpin.fa.part_003.gz

1. Split sequences into 4 parts

        $ seqkit split hairpin.fa.gz -p 4
        [INFO] split into 4 parts
        [INFO] read sequences ...
        [INFO] read 28645 sequences
        [INFO] write 7162 sequences to file: hairpin.fa.part_001.gz
        [INFO] write 7162 sequences to file: hairpin.fa.part_002.gz
        [INFO] write 7162 sequences to file: hairpin.fa.part_003.gz
        [INFO] write 7159 sequences to file: hairpin.fa.part_004.gz


    ***To reduce memory usage when spliting big file, we should alwasy use flag `--two-pass`***

        $ seqkit split hairpin.fa.gz -p 4 -2
        [INFO] split into 4 parts
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] create and read FASTA index ...
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 28645 sequences loaded
        [INFO] write 7162 sequences to file: hairpin.part_001.fa.gz
        [INFO] write 7162 sequences to file: hairpin.part_002.fa.gz
        [INFO] write 7162 sequences to file: hairpin.part_003.fa.gz
        [INFO] write 7159 sequences to file: hairpin.part_004.fa.gz

1. Split sequences by species. i.e. by custom IDs (first three letters)

        $ seqkit split hairpin.fa.gz -i --id-regexp "^([\w]+)\-" -2
        [INFO] split by ID. idRegexp: ^([\w]+)\-
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] create and read FASTA index ...
        [INFO] create FASTA index for hairpin.fa.gz.fa
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 28645 sequences loaded
        [INFO] write 48 sequences to file: hairpin.id_cca.fa.gz
        [INFO] write 3 sequences to file: hairpin.id_hci.fa.gz
        [INFO] write 106 sequences to file: hairpin.id_str.fa.gz
        [INFO] write 1 sequences to file: hairpin.id_bkv.fa.gz
        ...

1. Split sequences by sequence region (for example, sequence barcode)

        $ seqkit split hairpin.fa.gz -r 1:3 -2
        [INFO] split by region: 1:3
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] read sequence IDs and sequence region from FASTA file ...
        [INFO] create and read FASTA index ...
        [INFO] write 463 sequences to file: hairpin.region_1:3_AUG.fa.gz
        [INFO] write 349 sequences to file: hairpin.region_1:3_ACU.fa.gz
        [INFO] write 311 sequences to file: hairpin.region_1:3_CGG.fa.gz

    **If region is too long, we could use falg `--md5`**,
    i.e. use MD5 instead of region sequence in output file.

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

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 -o sample.fa.gz
        [INFO] sample by proportion
        [INFO] 2814 sequences outputed

1. Sample by number

        $ zcat hairpin.fa.gz | seqkit sample -n 1000 -o sample.fa.gz
        [INFO] sample by number
        [INFO] 949 sequences outputed

    ***To reduce memory usage when spliting big file, we could use flag `--two-pass`***

    ***We can also use `seqkit sample -p` followed with `seqkit head -n`:***

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 | seqkit head -n 1000 -o sample.fa.gz

1. Set rand seed to reproduce the result

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 -s 11

1. Most of the time, we could shuffle after sampling

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 | seqkit shuffle -o sample.fa.gz

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

        $ seqkit head -n 1 hairpin.fa.gz
        >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
        UAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA

1. FASTQ

        $ seqkit head -n 1 test.fastq.gz
```
        @HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG
        TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTG
        CCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGT
        ATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGA
        TCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT
        +
        HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHH
        HHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIII
        HIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHG
        IHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH

```


## replace

Usage

```bash 
replace name/sequence/by regular expression.

Note that the replacement supports capture variables.
e.g. $1 represents the text of the first submatch.
ATTENTION: use SINGLE quote NOT double quotes in *nix OS.

Examples: Adding space to all bases.

    seqkit replace -p "(.)" -r '$1 ' -s

Or use the \ escape character.

    seqkit replace -p "(.)" -r "\$1 " -s

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

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p " .+"
        >seq1
        ACGT-ACGT

1. Replace "-" with "="

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p "\-" -r '='
        >seq1 abc=123
        ACGT-ACGT

1. Remove gaps in sequences.

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p " |-" -s
        >seq1 abc-123
        ACGTACGT

1. Add space to every base. **ATTENTION: use SINGLE quote NOT double quotes in *nix OS**

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p "(.)" -r '$1 ' -s
        >seq1 abc-123
        A C G T - A C G T

1. Transpose sequence with [csvtk](https://github.com/shenwei356/csvtk)

        $ echo -e ">seq1\nACTGACGT\n>seq2\nactgccgt" | seqkit replace -p "(.)" -r     "\$1 " -s | seqkit seq -s -u | csvtk space2tab | csvtk -t transpose
        A       A
        C       C
        T       T
        G       G
        A       C
        C       C
        G       G
        T       T

1. Rename with number of record

        echo -e ">abc\nACTG\n>123\nATTT" |  seqkit replace -p .+ -r "seq_{NR}"
        >seq_1
        ACTG
        >seq_2
        ATTT


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

        $ seqkit shuffle hairpin.fa.gz > shuffled.fa
        [INFO] read sequences ...
        [INFO] 28645 sequences loaded
        [INFO] shuffle ...
        [INFO] output ...

1. ***For big genome, you'd better use two-pass mode*** so seqkit could use
   FASTA index to reduce memory usage

        $ time seqkit shuffle -2 GRCh38_latest_genomic.fna.gz > shuffle.fa
        [INFO] create and read FASTA index ...
        [INFO] create FASTA index for GRCh38_latest_genomic.fna.gz
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 194 sequences loaded
        [INFO] shuffle ...
        [INFO] output ...

        real    0m35.080s
        user    0m45.521s
        sys     0m3.411s

Note that when sampling on FASTQ files, make sure using same random seed by
flag `-s` (`--rand-seed`) for read 1 and 2 files.

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

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet
        >SEQ2
        acgtnAAAA
        >seq1
        ACGTNcccc

1. sort by ID, ignoring case.

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -i
        >seq1
        ACGTNcccc
        >SEQ2
        acgtnAAAA

1. sort by seq, ignoring case.

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -s -i
        >SEQ2
        acgtnAAAA
        >seq1
        ACGTNcccc

1. sort by sequence length

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAAnnn\n>seq3\nacgt" | seqkit sort --quiet -l
        >seq3
        acgt
        >seq1
        ACGTNcccc
        >SEQ2
        acgtnAAAAnnn


### Quick glance

1. Sequence number

        $ seqkit stat hairpin.fa.gz
        file           format  type  num_seqs    sum_len  min_len  avg_len  max_len
        hairpin.fa.gz  FASTA   RNA     28,645  2,949,871       39      103    2,354

1. First 10 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r 1:10 | seqkit sort -s | seqkit seq -s | head -n 10
        AAAAAAAAAA
        AAAAAAAAAA
        AAAAAAAAAG
        AAAAAAAAAG
        AAAAAAAAAG
        AAAAAAAAAU
        AAAAAAAAGG
        AAAAAAACAU
        AAAAAAACGA
        AAAAAAAUUA

    hmm, nothing special, non-coding RNA~


### Repeated hairpin sequences

We may want to check how may identical hairpins among different species there are.
`seqkit rmdup` could remove duplicated sequences by sequence content,
and save the replicates to another file (here is `duplicated.fa.gz`),
as well as replicating details (`duplicated.detail.txt`,
1th column is the repeated number,
2nd column contains sequence IDs seperated by comma).

    $ seqkit rmdup -s -i hairpin.fa.gz -o clean.fa.gz -d duplicated.fa.gz -D duplicated.detail.txt

    $ head -n 5 duplicated.detail.txt
    18      dre-mir-430c-1, dre-mir-430c-2, dre-mir-430c-3, dre-mir-430c-4, dre-mir-430c-5, dre-mir-430c-6, dre-mir-430c-7, dre-mir-430c-8, dre-mir-430c-9, dre-mir-430c-10, dre-mir-430c-11, dre-mir-430c-12, dre-mir-430c-13, dre-mir-430c-14, dre-mir-430c-15, dre-mir-430c-16, dre-mir-430c-17, dre-mir-430c-18
    16      hsa-mir-29b-2, mmu-mir-29b-2, rno-mir-29b-2, ptr-mir-29b-2, ggo-mir-29b-2, ppy-mir-29b-2, sla-mir-29b, mne-mir-29b, ppa-mir-29b-2, bta-mir-29b-2, mml-mir-29b-2, eca-mir-29b-2, aja-mir-29b, oar-mir-29b-1, oar-mir-29b-2, rno-mir-29b-3
    15      dme-mir-125, dps-mir-125, dan-mir-125, der-mir-125, dgr-mir-125-1, dgr-mir-125-2, dmo-mir-125, dpe-mir-125-2, dpe-mir-125-1, dpe-mir-125-3, dse-mir-125, dsi-mir-125, dvi-mir-125, dwi-mir-125, dya-mir-125
    13      hsa-mir-19b-1, ggo-mir-19b-1, age-mir-19b-1, ppa-mir-19b-1, ppy-mir-19b-1, ptr-mir-19b-1, mml-mir-19b-1, sla-mir-19b-1, lla-mir-19b-1, mne-mir-19b-1, bta-mir-19b, oar-mir-19b, chi-mir-19b
    13      hsa-mir-20a, ssc-mir-20a, ggo-mir-20a, age-mir-20, ppa-mir-20, ppy-mir-20a, ptr-mir-20a, mml-mir-20a, sla-mir-20, lla-mir-20, mne-mir-20, bta-mir-20a, eca-mir-20a

The result shows the most conserved miRNAs among different species,
`mir-29b`, `mir-125`, `mir-19b-1` and `mir-20a`.
And the `dre-miR-430c` has the most multicopies in *Danio rerio*.

### Hairpins in different species

1. Before spliting by species, let's take a look at the sequence names.

        $ seqkit seq hairpin.fa.gz -n | head -n 3
        cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop
        cel-mir-1 MI0000003 Caenorhabditis elegans miR-1 stem-loop

    The first three letters (e.g. `cel`) are the abbreviation of species names.
    So we could split hairpins by the first letters by defining custom
    sequence ID parsing regular expression `^([\w]+)\-`.

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

        $ seqkit split hairpin.fa.gz -i --id-regexp "^([\w]+)\-" --two-pass

    ***To reduce memory usage when splitting big file, we should always use flag `--two-pass`***

2. Species with most miRNA hairpins. Third column is the sequences number.

        $ cd hairpin.fa.gz.split/;
        $ seqkit stat hairpin.id_* | csvtk space2tab | csvtk -t sort -k num_seqs:nr | csvtk -t pretty| more
        file                     format   type   num_seqs   sum_len   min_len   avg_len   max_len
        hairpin.id_GRCh38_latest_genomic.fna.gzsta     FASTA    RNA    1,881      154,242   82        82        82
        hairpin.id_mmu.fasta     FASTA    RNA    1,193      107,370   90        90        90
        hairpin.id_bta.fasta     FASTA    RNA    808        61,408    76        76        76
        hairpin.id_gga.fasta     FASTA    RNA    740        42,180    57        57        57
        hairpin.id_eca.fasta     FASTA    RNA    715        89,375    125       125       125
        hairpin.id_mtr.fasta     FASTA    RNA    672        231,840   345       345       345

    Here, a CSV/TSV tool [csvtk](https://github.com/shenwei356/csvtk)
    is used to sort and view the result.

For human miRNA hairpins

1. Length distribution.
 `seqkit fx2tab` could show extra information like sequence length, GC content.
 A distribution ploting script is used, (
 [plot_distribution.py](https://github.com/shenwei356/bio_scripts/blob/master/plot/plot_distribution.py) )

        $ seqkit fx2tab hairpin.id_GRCh38_latest_genomic.fna.gz.gz -l | cut -f 3  | plot_distribution.py -o hairpin.id_GRCh38_latest_genomic.fna.gz.gz.lendist.png

    ![hairpin.id_GRCh38_latest_genomic.fna.gz.gz.lendist.png](/files/hairpin/hairpin.id_GRCh38_latest_genomic.fna.gz.gz.lendist.png)


## Bacteria genome

### Dataset

[Pseudomonas aeruginosa PAO1](http://www.ncbi.nlm.nih.gov/nuccore/110645304),
files:

-  Genbank file [`PAO1.gb`](/files/PAO1/PAO1.gb)
-  Genome FASTA file [`PAO1.fasta`](/files/PAO1/PAO1.fasta)
-  GTF file [`PAO1.gtf`](/files/PAO1/PAO1.gtf) was created with [`extract_features_from_genbank_file.py`](https://github.com/shenwei356/bio_scripts/blob/master/file_formats/extract_features_from_genbank_file.py), by

        extract_features_from_genbank_file.py  PAO1.gb -t . -f gtf > PAO1.gtf


### Motif distribution

Motifs

    $ cat motifs.fa
    >GTAGCGS
    GTAGCGS
    >GGWGKTCG
    GGWGKTCG

1. Sliding. Remember flag `--id-ncbi`, do you?
  By the way, do not be scared by the long flag `--circle-genome`, `--step`
  and so on. They have short ones, `-c`, `-s`

        $ seqkit sliding --id-ncbi --circular-genome --step 20000 --window 200000 PAO1.fasta -o PAO1.fasta.sliding.fa

        $ seqkit stat PAO1.fasta.sliding.fa
        file                   format  type  num_seqs     sum_len  min_len  avg_len  max_len
        PAO1.fasta.sliding.fa  FASTA   DNA        314  62,800,000  200,000  200,000  200,000

1. Locating motifs

        $ seqkit locate --id-ncbi --ignore-case --degenerate --pattern-file motifs.fa  PAO1.fasta.sliding.fa -o  PAO1.fasta.sliding.fa.motifs.tsv

1. Ploting distribution ([plot_motif_distribution.R](/files/PAO1/plot_motif_distribution.R))

        # preproccess
        $ perl -ne 'if (/_sliding:(\d+)-(\d+)\t(.+)/) {$loc= $1 + 100000; print "$loc\t$3\n";} else {print}' PAO1.fasta.sliding.fa.motifs.tsv  > PAO1.fasta.sliding.fa.motifs.tsv2

        # plot
        $ ./plot_motif_distribution.R

    Result

    ![motif_distribution.png](files/PAO1/motif_distribution.png)


### Find multicopy genes

1. Get all CDS sequences

        $ seqkit subseq --id-ncbi --gtf PAO1.gtf --feature cds PAO1.fasta -o PAO1.cds.fasta

        $ seqkit stat *.fasta
        file            format  type  num_seqs    sum_len    min_len    avg_len    max_len
        PAO1.cds.fasta  FASTA   DNA      5,572  5,593,306         72    1,003.8     16,884
        PAO1.fasta      FASTA   DNA          1  6,264,404  6,264,404  6,264,404  6,264,404

1. Get duplicated sequences

        $ seqkit rmdup --by-seq --ignore-case PAO1.cds.fasta -o PAO1.cds.uniq.fasta --dup-seqs-file PAO1.cds.dup.fasta --dup-num-file PAO1.cds.dup.text

        $ cat PAO1.cds.dup.text
        6       NC_002516.2_500104:501120:-, NC_002516.2_2556948:2557964:+, NC_002516.2_3043750:3044766:-, NC_002516.2_3842274:3843290:-, NC_002516.2_4473623:4474639:+, NC_002516.2_5382796:5383812:-
        2       NC_002516.2_2073555:2075438:+, NC_002516.2_4716660:4718543:+
        2       NC_002516.2_2072935:2073558:+, NC_002516.2_4716040:4716663:+
        2       NC_002516.2_2075452:2076288:+, NC_002516.2_4718557:4719393:+

### Flanking sequences

1. Get CDS and 1000 bp upstream sequence

        $ seqkit subseq --id-ncbi --gtf PAO1.gtf --feature cds PAO1.fasta --up-stream 1000

1. Get 1000 bp upstream sequence of CDS, *NOT* including CDS.

        $ seqkit subseq --id-ncbi --gtf PAO1.gtf --feature cds PAO1.fasta --up-stream 1000 --only-flank

<div id="disqus_thread"></div>
<script>
/**
* RECOMMENDED CONFIGURATION VARIABLES: EDIT AND UNCOMMENT THE SECTION BELOW TO INSERT DYNAMIC VALUES FROM YOUR PLATFORM OR CMS.
* LEARN WHY DEFINING THESE VARIABLES IS IMPORTANT: https://disqus.com/admin/universalcode/#configuration-variables
*/
/*
var disqus_config = function () {
this.page.url = PAGE_URL; // Replace PAGE_URL with your page's canonical URL variable
this.page.identifier = PAGE_IDENTIFIER; // Replace PAGE_IDENTIFIER with your page's unique identifier variable
};
*/
(function() { // DON'T EDIT BELOW THIS LINE
var d = document, s = d.createElement('script');

s.src = '//fastakit.disqus.com/embed.js';

s.setAttribute('data-timestamp', +new Date());
(d.head || d.body).appendChild(s);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>



    # Benchmark


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

A Python script [memusg](https://github.com/shenwei356/memusg) was used
to compute running time and peak memory usage of a process.

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

## Datasets

All test data is available here: [seqkit-benchmark-data.tar.gz](http://app.shenwei.me/data/seqkit/seqkit-benchmark-data.tar.gz)  (2.2G)

### dataset_A.fa - large number of short sequences

Dataset A is reference genomes DNA sequences of gastrointestinal tract from
[NIH Human Microbiome Project](http://hmpdacc.org/):
[`Gastrointestinal_tract.nuc.fsa`](http://downloads.hmpdacc.org/data/reference_genomes/body_sites/Gastrointestinal_tract.nuc.fsa) (FASTA format, ~2.7G).

### dataset_B.fa - small number of large sequences

Dataset B is Human genome from [ensembl](http://uswest.ensembl.org/info/data/ftp/index.html).

- Genome DNA:  [`GRCh38_latest_genomic.fna.gz`](ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/GRCh38_latest_genomic.fna.gz) (Gzipped FASTA file, ~900M)
. Decompress it and rename to dataset_B.fa (~2.9G).
- GTF file:  [`GRCh38_latest_genomic.gtf.gz`](ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/GRCh38_latest_genomic.gtf.gz) (~44M)
- BED file: `Homo_sapiens.GRCh38.84.bed.gz` was converted from `GRCh38_latest_genomic.gtf.gz` by  [`gtf2bed`](http://bedops.readthedocs.org/en/latest/content/reference/file-management/conversion/gtf2bed.html?highlight=gtf2bed)  with command

        $ zcat GRCh38_latest_genomic.gtf.gz | gtf2bed --do-not-sort | gzip -c > Homo_sapiens.GRCh38.84.bed.gz


### dataset_C.fq – Illumina single end reads (SE100)

Dataset C is Illumina single end (SE 100bp) reads file (~2.2G).


Summary

    $ seqkit stat *.fa
    file          format  type   num_seqs        sum_len  min_len       avg_len      max_len
    dataset_A.fa  FASTA   DNA      67,748  2,807,643,808       56      41,442.5    5,976,145
    dataset_B.fa  FASTA   DNA         194  3,099,750,718      970  15,978,096.5  248,956,422
    dataset_C.fq  FASTQ   DNA   9,186,045    918,604,500      100           100          100
    
### Sequence ID list

Parts of sequences IDs was sampled and shuffled from original data.
They were used in test of extracting sequences by ID list.

Commands:

    $ seqkit sample -p 0.3  dataset_A.fa | seqkit seq --name --only-id | shuf > ids_A.txt
    $ seqkit sample -p 0.3  dataset_B.fa | seqkit seq --name --only-id | shuf > ids_B.txt    
    $ seqkit sample -p 0.03 dataset_C.fq | seqkit seq --name --only-id | shuf > ids_C.txt

Numbers:

    $ wc -l ids*.txt
        20138 ids_A.txt
        58 ids_B.txt
    2754516 ids_C.txt

### BED file

Only BED data of chromosome 19 was used in test of subsequence with BED file:

    $ zcat GRCh38_latest_genomic.bed.gz| grep -E "^19" | gzip -c > chr19.bed.gz



