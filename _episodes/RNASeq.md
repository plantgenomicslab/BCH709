---
layout: page
title: Compile and Software installation (II)
published: true
---

{% include gh_variables.html %}

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8


## RNA Sequencing
![RNA Sequencing]({{{site.baseurl}}/fig/rnaseq.png)

1. The transcriptome is spatially and temporally dynamic
2. Data comes from functional units (coding regions)
3. Only a tiny fraction of the genome

>## Introduction
>Sequence based assays of transcriptomes (RNA-seq) are in wide use because of their favorable properties for quantification, transcript discovery and splice isoform identification, as well as adaptability for numerous more specialized measurements. RNA-Seq studies present some challenges that are shared with prior methods such as microarrays and SAGE tagging, and they also present new ones that are specific to high-throughput sequencing platforms and the data they produce. This document is part of an ongoing effort to provide the community with standards and guidelines that will be updated as RNASeq matures and to highlight unmet challenges. The intent is to revise this document periodically to capture new advances and increasingly consolidate standards and best practices.
>
>RNA-Seq experiments are diverse in their aims and design goals, currently including multiple types of RNA isolated from whole cells or from specific sub-cellular compartments or biochemical classes, such as total polyA+ RNA, polysomal RNA, nuclear ribosome-depleted RNA, various size fractions of RNA and a host of others. The goals of individual experiments range from major transcriptome “discovery” that seeks to define and quantify all RNA species in a starting RNA sample to experiments that simply need to detect significant changes in the more abundant RNA classes across many samples.  
{: .prereq}


![RNA Sequencing workflow]({{{site.baseurl}}/fig/rnaseq_workflow.png)

### Seven stages to data science
1. Define the question of interest
2. Get the data
3. Clean the data
4. Explore the data
5. Fit statistical models
6. Communicate the results
7. Make your analysis reproducible

### What do we need to prepare ?
##### Sample Information
a. What kind of material it is should be noted: Tissue, cell line, primary cell type, etc…
b. It’s ontology term (a DCC wrangler will work with you to obtain this)
c. If any treatments or genetic modifications (TALENs, CRISPR, etc…) were done to the sample
prior to RNA isolation.
d. If it’s a subcellular fraction or derived from another sample. If derived from another sample,
that relationship should be noted.
e. Some sense of sample abundance: RNA-Seq data from “bulk” vs. 10,000 cell equivalents can
give very different results, with lower input samples typically being less reproducible. Having
a sense of the amount of starting material here is useful.
f. If you received a batch of primary or immortalized cells, the lot #, cat # and supplier should be
noted.
g. If cells were cultured out, the protocol and methods used to propagate the cells should be
noted.
h. If any cell phenotyping or other characterizations were done to confirm it’s identify, purity,
etc.. those methods should be noted.




##### RNA Information: 
RNAs come in all shapes and sizes. Some of the key properties to report are:
a. Total RNA, Poly-A(+) RNA, Poly-A(-) RNA
b. Size of the RNA fraction: we typically have a + 200 and – 200 cutoff, but there is a wide
range, i.e. microRNA-sized, etc…
c. If the RNA was treated with Ribosomal RNA depletion kits (RiboMinus, RiboZero): please
note the kit used.

##### Protocols: 
There are several methods used to isolate RNAs with that work fine for the purposes of RNA-Seq. For all the ENCODE libraries that we make, we provide a document that lists in detail:
a. The RNA isolation methods,
b. Methods of size selections
c. Methods of rRNA removal
d. Methods of oligo-dT selections
e. Methods of DNAse I treatments

#### Experimental Design
- Balanced design
- Technical replicates not necessary (Marioni et al., 2008)
- Biological replicates: 6 - 12 (Schurch et al., 2016)
- Power analysis

```
[Paul L. Auer and R. W. Doerge "Statistical Design and Analysis of RNA Sequencing Data" Genetics June 1, 2010 vol. 185 no.2 405-416](https://www.genetics.org/content/185/2/405)
pBusby, Michele A., et al. "Scotty: a web tool for designing RNA-Seq experiments to measure differential gene expression." Bioinformatics 29.5 (2013): 656-657](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582267/)
[Marioni, John C., et al. "RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays." Genome research (2008)](https://genome.cshlp.org/content/18/9/1509.full.html)
[Schurch, Nicholas J., et al. "How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?." Rna (2016)](https://rnajournal.cshlp.org/content/22/6/839.long)
[Zhao, Shilin, et al. "RnaSeqSampleSize: real data based sample size estimation for RNA sequencing." BMC bioinformatics 19.1 (2018): 191](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5)

```

#### Replicate number
In all cases, experiments should be performed with two or more biological replicates, unless there is a
compelling reason why this is impractical or wasteful (e.g. overlapping time points with high temporal
resolution). A biological replicate is defined as an independent growth of cells/tissue and subsequent
analysis. Technical replicates from the same RNA library are not required, except to evaluate cases
where biological variability is abnormally high. In such instances, separating technical and biological
variation is critical. In general, detecting and quantifying low prevalence RNAs is inherently more variable
than high abundance RNAs. As part of the ENCODE pipeline, annotated transcript and genes are
quantified using RSEM and the values are made available for downstream correlation analysis. Replicate
concordance: the gene level quantification should have a Spearman correlation of >0.9 between
isogenic replicates and >0.8 between anisogenic replicates.


#### RNA extraction
- Sample processing and storage
- Total RNA/mRNA/small RNA
- DNAse treatment
- Quantity & quality
- RIN values (Strong effect)
- Batch effect
- Extraction method bias (GC bias)

```
Romero, Irene Gallego, et al. "RNA-seq: impact of RNA degradation on transcript quantification." BMC biology 12.1 (2014): 42
Kim, Young-Kook, et al. "Short structured RNAs with low GC content are selectively lost during extraction from a small number of cells." Molecular cell 46.6 (2012): 893-89500481-9).
```

##### RNA Quantification and Quality Control: When working with bulk samples, throughout the various
steps we periodically assess the quality and quantity of the RNA. This is typically done on a
BioAnalyzer. Points to check are:
a. Total RNA
b. After oligo-dT size selections
c. After rRNA-depletions
d. After library construction

#### Library prep
- PolyA selection
- rRNA depletion
- Size selection
- PCR amplification (See section PCR duplicates)
- Stranded (directional) libraries
   - Accurately identify sense/antisense transcript
   - Resolve overlapping genes
- Exome capture
- Library normalisation
- Batch effect

![RNA library]({{{site.baseurl}}/fig/library.png)

![RNA Sequencing tool]({{{site.baseurl}}/fig/frag.png)  


##### Sequencing: 
There are several sequencing platforms and technologies out there being used. It is important to provide the following pieces of information:
a. Platform: Illumina, PacBio, Oxford Nanopore, etc…
b. Format: Single-end, Pair-end,
c. Read Length: 101 bases, 125 bases, etc…
d. Unusual barcode placement and sequence: Some protocols introduce barcodes in noncustomary places. If you are going to deliver a FASTQ file that will contain the barcode
sequences in it or other molecular markers – you will need to report both the position in the
read(s) where they are and their sequence(s).
e. Please provide the sequence of any custom primers that were used to sequence the library

![RNA library]({{{site.baseurl}}/fig/sequencing.png)



##### Sequencing depth.
The amount of sequencing needed for a given sample is determined by the goals of the experiment and
the nature of the RNA sample. Experiments whose purpose is to evaluate the similarity between the
transcriptional profiles of two polyA+ samples may require only modest depths of sequencing.
Experiments whose purpose is discovery of novel transcribed elements and strong quantification of
known transcript isoforms requires more extensive sequencing.
• Each Long RNA-Seq library must have a minimum of 30 million aligned reads/mate-pairs.
• Each RAMPAGE library must have a minimum of 20 million aligned reads/mate-pairs.
• Each small RNA-Seq library must have a minimum of 30 million aligned reads/mate-pairs.

##### Quantitative Standards (spike-ins).
It is highly desirable to include a ladder of RNA spike-ins to calibrate quantification, sensitivity, coverage
and linearity. Information about the spikes should include the stage of sample preparation that the spiked
controls were added, as the point of entry affects use of spike data in the output. In general, introducing
spike-ins as early in the process as possible is the goal, with more elaborate uses of different spikes at
different steps being optional (e.g. before poly A+ selection, at the time of cDNA synthesis, or just prior to
sequencing). Different spike-in controls are needed for each of the RNA types being analyzed (e.g. long
RNAs require different quantitative controls from short RNAs). Such standards are not yet available for all
RNA types. Information about quantified standards should also include:
a) A FASTA (or other standard format) file containing the sequences of each spike in.
b) Source of the spike-ins (home-made, Ambion, etc..)
c) The concentration of each of the spike-ins in the pool used. 

[Hong et al., 2016, Principles of metadata organization at the ENCODE data coordination center.](https://academic.oup.com/database/article-lookup/doi/10.1093/database/baw001)



![RNA Sequencing tool]({{{site.baseurl}}/fig/rna_software.png)


#### QC FAIL?
https://sequencing.qcfail.com/



### Fastq format
FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores.


The format is similar to fasta though there are differences in syntax as well as integration of quality scores. Each sequence requires at least 4 lines:

1. The first line is the sequence header which starts with an ‘@’ (not a ‘>’!).
Everything from the leading ‘@’ to the first whitespace character is considered the sequence identifier.
Everything after the first space is considered the sequence description
2. The second line is the sequence.
3. The third line starts with ‘+’ and can have the same sequence identifier appended (but usually doesn’t anymore).
4. The fourth line are the quality scores

The FastQ sequence identifier generally adheres to a particular format, all of which is information related to the sequencer and its position on the flowcell. The sequence description also follows a particular format and holds information regarding sample information.


```
$ mkdir rnaseq

$ cd rnaseq/

$ pwd

$ wget https://www.dropbox.com/s/pfu2r922esr1ccf/paired2.fastq.gz

$ wget https://www.dropbox.com/s/o642h0ca6djdbib/paired1.fastq.gz

$ ls -algh

$ zcat paired2.fastq.gz | head
 ```

```
@A00261:180:HL7GCDSXX:2:1101:30572:1047/2
AAAATACATTGATGACCATCTAAAGTCTACGGCGTATGCGACTGATGAAGTATATTGCACCACCTGAGGGTGATGCTAATACTACTGTTGACGATAATGCTGATCTTCTTGCTAAGCTTAATATTGTTGGTGTTGAACCTAATGTTGGTG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF
@A00261:180:HL7GCDSXX:2:1101:21088:1094/2
ATCTCACATCGTTCCCTCAAGATTCTGAATTTTGGCAGCTCATTGCATTCTGTGCCGGCACTGGTGGTTCGATGCTTGTCATTGGTTCTGCTGCTGGTGTAGCCTTCATGGGGATGGAGAAAGTCGATTTCTTTTGGTATTTCCGAAAGG
+
FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00261:180:HL7GCDSXX:2:1101:21251:1125/2
CGGTGGAAAAGGAAACAGCTTTGGAAGGTTGATTCCATTACAGATTCGATTCGAAACTATGGTTCAGATTTCCGATCTTCCACGGGATTTGACAGAGGAGGTGCTCTCTAGGATTCCGGTGACATCTATGAGAGCAGTGAGATTTACTTG
```

- A00261  : Instrument name
- 180 : run ID
- HL7GCDSXX : Flowcell ID
- 2 : Flowcell lane
- 1101 : tile number within the flowcell lane
- 30572 : X-coordinate of the cluster within tile
- 1047 : Y-coordinate of the cluster within tile
- /2 : member of a pair 1 or 2 (Paired end reads only)


#### Quality Scores
Quality scores are a way to assign confidence to a particular base within a read. Some sequencers have their own proprietary quality encoding but most have adopted Phred-33 encoding. Each quality score represents the probability of an incorrect basecall at that position.

##### Phred Quality Score Encoding
Quality scores started as numbers (0-40) but have since changed to an ASCII encoding to reduce filesize and make working with this format a bit easier, however they still hold the same information. ASCII codes are assigned based on the formula found below. This table can serve as a lookup as you progress through your analysis.

##### Quality Score Interpretation
Once you know what each quality score represents you can then use this chart to understand the confidence in a particular base.


![FASTQ quality]({{{site.baseurl}}/fig/quality.png)


#### Conda enviroment

```
$ conda create -n rnaseq


$ conda activate rnaseq


```



### Read QC
- Number of reads
- Per base sequence quality
- Per sequence quality score
- Per base sequence content
- Per sequence GC content
- Per base N content
- Sequence length distribution
- Sequence duplication levels
- Overrepresented sequences
- Adapter content
- Kmer content

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```
$ conda install fastqc
```
***if you have error, please add below first***
```
$ conda config --add channels conda-forge
$ conda config --add channels defaults
$ conda config --add channels r
$ conda config --add channels bioconda
```

##### Run fastqc
```
$ fastqc --help
$ fastqc -t <YOUR CPU COUNT> paired1.fastq.gz  paired2.fastq.gz

```

### How to make a report?
![MultiQC]({{{site.baseurl}}/fig/multiqc.png)
[MultiQC](https://multiqc.info/)
```
$ conda install multiqc
$ multiqc --help
$ multiqc .
$ cp -r multiqc* <YOUR DESKTOP FOLDER>
```
### Trim the reads
- Trim IF necessary
   - Synthetic bases can be an issue for SNP calling
   - Insert size distribution may be more important for assemblers
- Trim/Clip/Filter reads
- Remove adapter sequences
- Trim reads by quality
- Sliding window trimming
- Filter by min/max read length
- Remove reads less than ~18nt
- Demultiplexing/Splitting

![Trimming]({{{site.baseurl}}/fig/trim.png)

[Cutadapt](https://github.com/marcelm/cutadapt/)
[fastp](https://github.com/OpenGene/fastp)
[Skewer](https://github.com/relipmoc/skewer)
[Prinseq](http://prinseq.sourceforge.net/)
[Trimmomatics](http://www.usadellab.org/cms/?page=trimmomatic)
[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

##### Install Trim Galore
```
$ conda install trim-galore
```

##### Run trimming
```
$ trim_galore --help

$ trim_galore --paired   --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 2  --max_n 40  --gzip -o trim paired1.fastq.gz paired2.fastq.gz 

$ multiqc .
```

### Align the reads (mapping)
![mapping]({{{site.baseurl}}/fig/mapping.png)
 - Aligning reads back to a reference sequence
 - Mapping to genome vs transcriptome
 - Splice-aware alignment (genome)

[STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537)
[HISAT2](https://www.nature.com/articles/s41587-019-0201-4)
[GSNAP](https://dx.doi.org/10.1007/978-1-4939-3578-9_15)
[Bowtie2](https://www.nature.com/articles/nmeth.1923)
[Novoalign](http://www.novocraft.com/products/novoalign/)

[Baruzzo, Giacomo, et al. "Simulation-based comprehensive benchmarking of RNA-seq aligners." Nature methods 14.2 (2017): 135](https://www.nature.com/articles/nmeth.4106)

![algorithm]({{{site.baseurl}}/fig/algorithm.png)
![banana]({{{site.baseurl}}/fig/banana.png)
![banana]({{{site.baseurl}}/fig/BackwardMatching.png)

#### Install HISAT2 (graph FM index, spin off version Burrows-Wheeler Transform)
```
$ conda install hisat2
```
##### Download reference sequence
```
$ wget https://www.dropbox.com/s/0onch14nnxx9b94/bch709.fasta
```
#### HISAT2 indexing
```
$ hisat2-build --help
$ hisat2-build bch709.fasta bch709
```

#### HISAT2 mapping

```
$ hisat2 -x bch709 --threads <YOUR CPU COUNT> -1 trim/paired1_val_1.fq.gz -2 trim/paired2_val_2.fq.gz  -S align.sam
```

#### SAM file format

Check result
```
$ head align.sam
```
#### SAM file

```
<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [<TAG>:<VTYPE>:<VALUE> [...]]
```
![sam1]({{{site.baseurl}}/fig/sam1.png)

##### SAM flag
![flag]({{{site.baseurl}}/fig/flag.jpg)

##### Demical to binary
```
$ echo 'obase=163;10' | bc
```
If you don't have bc, please install through conda
```
$ conda install bc
```



##### SAM tag
There are a bunch of predefined tags, please see the SAM manual for more information. For the tags used in this example:

Any tags that start with X? are reserved fields for end users: XT:A:M, XN:i:2, XM:i:0, XO:i:0, XG:i:0

![samtag]({{{site.baseurl}}/fig/samtag.png)

##### More information is below.
http://samtools.github.io/hts-specs/


#### SAMtools
SAM (Sequence Alignment/Map) format is a generic format for storing large nucleotide sequence alignments. SAM aims to be a format that:

- Is flexible enough to store all the alignment information generated by various alignment programs;
- Is simple enough to be easily generated by alignment programs or converted from existing alignment formats;
- Is compact in file size;
- Allows most of operations on the alignment to work on a stream without loading the whole alignment into memory;
- Allows the file to be indexed by genomic position to efficiently retrieve all reads aligning to a locus.

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. http://samtools.sourceforge.net/
```
$ samtools view -Sb align.sam > align.bam
$ samtools sort align.bam  align_sort
$ samtools index align_sort.bam
```

##### BAM file
A BAM file (.bam) is the binary version of a SAM file. A SAM file (.sam) is a tab-delimited text file that contains sequence alignment data. 

|SAM/BAM|size|
|-----|----|
|align.sam | 903M|
|align.bam | 166M|


#### Alignment visualization
```
samtools tview align_sort.bam bch709.fasta
```
![tview]({{{site.baseurl}}/fig/tview.png)

[IGV](https://software.broadinstitute.org/software/igv/)
SAMTOOLS Tview
[Tablet](https://ics.hutton.ac.uk/tablet/)


#### Alignment QC
- Number of reads mapped/unmapped/paired etc
- Uniquely mapped
- Insert size distribution
- Coverage
- Gene body coverage
- Biotype counts / Chromosome counts
- Counts by region: gene/intron/non-genic
- Sequencing saturation
- Strand specificity

samtools > stats 
bamtools > stats 
[QoRTs](https://hartleys.github.io/QoRTs/)
[RSeQC](http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/)
[Qualimap](http://qualimap.bioinfo.cipf.es/)


#### Quantification • Counts
![genecount]({{{site.baseurl}}/fig/genecount.png)
- Read counts = gene expression
- Reads can be quantified on any feature (gene, transcript, exon etc)
- Intersection on gene models
- Gene/Transcript level

![count]({{{site.baseurl}}/fig/count.png)

[featureCounts](http://subread.sourceforge.net/)
[HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/)
[RSEM](https://deweylab.github.io/RSEM/)
[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)
[Rcount](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu680)

```
conda install -c bioconda subread
conda install -c bioconda rsem
```
#### Quantification method
- PCR duplicates
 - Ignore for RNA-Seq data
 - Computational deduplication (Don't!)
 - Use PCR-free library-prep kits
 - Use UMIs during library-prep

- Multi-mapping
 - Added (BEDTools multicov)
 - Discard (featureCounts, HTSeq)
 - Distribute counts (Cufflinks)
 - Rescue
   - Probabilistic assignment (Rcount, Cufflinks)
   - Prioritise features (Rcount)
   - Probabilistic assignment with EM (RSEM)

 Fu, Yu, et al. "Elimination of PCR duplicates in RNA-seq and small RNA-seq using unique molecular identifiers." BMC genomics 19.1 (2018): 531
 Parekh, Swati, et al. "The impact of amplification on differential expression analyses by RNA-seq." Scientific reports 6 (2016): 25533
 Klepikova, Anna V., et al. "Effect of method of deduplication on estimation of differential gene expression using RNA-seq." PeerJ 5 (2017): e3091


#### Differential expression
DESeq2
edgeR (Neg-binom > GLM > Test)
Limma-Voom (Neg-binom > Voom-transform > LM > Test)

```
conda install -c bioconda bioconductor-deseq2
```

#### Functional analysis • GO
Gene enrichment analysis (Hypergeometric test)
Gene set enrichment analysis (GSEA)
Gene ontology / Reactome databases

#### Conda deactivate
```
$ conda deactivate
$ conda env remove --name rnaseq
```

Conesa, Ana, et al. "A survey of best practices for RNA-seq data analysis." Genome biology 17.1 (2016): 13


>## HOME WORK due next Tuesday
>1. Create your enviroment name rnaseq, install following tools
>- star
>- bioconductor-limma
>- fastqc
>- rsem
>- subread
>- hisat2
>- samtools
>- bowtie2
>- bioconductor-deseq2
>- trim-galore 
>- multiqc
>- trinity 
>- qorts
>2. Then export your enviroment to rnaseq.yaml (conda env export ......)
>3. copy contents of rnaseq.yaml and paste to [this site](https://docs.google.com/forms/d/e/1FAIpQLSe0faV4UHKEQ8CnqJ-iXOTAmKM2IxN6g7ZNSSmtcw5FuPwmWA/viewform?usp=sf_link).
{: .solution}



## Last class question
1. I can't paste (Ubuntu)
![paste]({{{site.baseurl}}/fig/paste.png)
![paste2]({{{site.baseurl}}/fig/paste2.png)

2. ChIP-Seq
Maybe we need a special session.


## Check your CPUs and Memory

```
$ lscpu
$ free
$ htop
```

If your have error ***command not found***
Please check below site.
https://command-not-found.com/

```Mac
brew install htop util-linux
```

```Ubuntu
sudo apt install htop util-linux
```

### Prompt Customization

```
echo '###BCH709 ' >> ~/.bashrc

echo 'export PS1="\[\033[38;5;2m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;33m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;166m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;4m\]\\$\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bashrc

echo "alias ls='ls --color=auto'" >> ~/.bashrc

source ~/.bashrc
```
More information is [here](https://plantgenomicslab.github.io/BCH709/bash/index.html)

```
echo '###BCH709 ' >> ~/.bash_profile

echo 'export PS1="\[\033[38;5;2m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;33m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;166m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;4m\]\\$\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bash_profile

echo "alias ls='ls --color=auto'" >> ~/.bash_profile

source ~/.bash_profile
```


## Software

| Software | Version | Manual | Available for | Description |
| -------- | ------------ | ------ | ------------- | ----------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.7 | [Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)| Linux, MacOS, Windows | Quality control tool for high throughput sequence data. |
| [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) | 2.1.0| [Link](https://ccb.jhu.edu/software/hisat2/index.shtml) | Linux, MacOS, Windows | Mapping RNA sequences against genome |
| [BWA](http://bio-bwa.sourceforge.net/) | 0.7.17 | [Link](http://bio-bwa.sourceforge.net/bwa.shtml) | Linux, MacOS | Mapping DNA sequences against reference genome. |

## QuickStart Software Installation Instructions
​
First, I would recommend to put all the software under `~/bch709/bin` folder.
Please use `mkdir -p ~/bch709/bin` 


### Download example
Second please make example folder below `~/bch709/`
Such as `~/bch709/example`  and download below files.

#### Download reads example 
```
https://www.dropbox.com/s/p11a5cw5lb5y0un/test.fq
```

#### Download genome example
```
https://www.dropbox.com/s/q2srdymrp76jevm/test_genome.fasta
```

## Current location
Please check your current location and move to   `~/bch709/example`


## FastQC Source Code Installation
Please download and unzip
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
> ## MacOS
> 
>To install unzip: 
>
>~~~
>$ brew install unzip
>~~~
>{: .bash}
> . 
{: .solution}

> ## On Ubuntu systems
> 
>To install unzip: 
>
>~~~
>$ sudo apt install unzip
>~~~
>{: .bash}
> . 
{: .solution}

### Please check  INSTALL.txt & README.md
```
less <FILENAME>
```

### Run FastQC

**Test your installation by running:**

```
./fastqc 
```

### Which error did you get?



### JAVA Installation
```bash
sudo apt install default-jre
```

```bash
brew cask install java
```


### Run FastQC
```
./fastqc <your input>
```

### Result
- HTML file

### Link Your Directory (Short Cut in WINDOWS)
```bash
$ mkdir /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop 
$ ln -s /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop ~/bch709/results
```

```bash
$ mkdir ~/Desktop/BCH709_Desktop
$ ln -s ~/Desktop/BCH709_Desktop ~/bch709/results
```

### Check your results
unzip your results




## HISAT2
https://ccb.jhu.edu/software/hisat2/index.shtml

![hisat2]({{{site.baseurl}}/fig/hisat2.png)


### Check HISAT2 manual
```
less MANUAL
```

### compile your HISAT2
use make

### How to run HISAT2?
- need to build index first `hisat2-build`
```bash
./hisat2-build <YOUR_GENOME_SEQUENCE> <YOUR_GENOME_INDEX>
```
- map your reads by `hisat2`
```bash
./hisat2 -x <YOUR_GENOME_INDEX> -U <SEQUENCING_READS> -S output.sam
```

## Different way to install software?
## Github
### What is Git and why use it?
Git is an open source distributed version control system that facilitates GitHub activities on your laptop or desktop. Version control using Git is the most reasonable way to keep track of changes in code, manuscripts, presentations, and data analysis projects.  

### Why Version Control?
Version control is essential in creating any project that takes longer than 5 minutes to complete. Even if your memory is longer than 5 minutes, next month you are not likely to be able to retrace your steps.  
![github-workflow]({{{site.baseurl}}/fig/git_overview.png)

![github]({{{site.baseurl}}/fig/github_dri.png)

## HISAT2 installation from Github
```
https://github.com/DaehwanKimLab/hisat2

```
![github-workflow]({{{site.baseurl}}/fig/hisat2_github.png)
![github-workflow]({{{site.baseurl}}/fig/hisat2_github2.png)

### Download
```bash

$ cd ~/bch709/bin

$ git clone https://github.com/DaehwanKimLab/hisat2.git

$ cd hisat2

$ ls -algh

$ less MANUAL
```
![github-workflow]({{{site.baseurl}}/fig/hisat2_github3.png)

### Compile
```bash
$ make -j <YOUR CPU>
```

## HISAT2 installation from binary
![github-workflow]({{{site.baseurl}}/fig/hisat2_binary.png)


### BWA
​​
> ## BWA Source Code Installation
>
> If you prefer to install from source, follow the instructions below:
>
> ~~~
> $ cd ~/bch709/bin
> $ curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
> $ tar xvf bwa-0.7.17.tar.bz2
> $ cd bwa-0.7.17
> $ make
> ~~~
> {: .bash}
{: .solution}
​
### Please check README.md
```
less <FILENAME>
```
#### Tell me your error


>## How to solve it?
>~~~
>sudo apt install zlib1g-dev
>~~~
> {: .bash}
>~~~
>brew install zlib
>~~~
> {: .bash}
{: .solution}

**Test your installation by running:**
​
~~~
$ ./bwa
~~~
{: .bash}

### BWA Github Installation
***Seach BWA on Google***


### How to run BWA?
```
bwa index <YOUR_GENOME_SEQUENCE>
bwa mem  <YOUR_GENOME_SEQUENCE> <SEQUENCING_READS>
```


## Conda?
- Dependencies is one of the main reasons to use Conda.
Sometimes, install a package is not as straight forward as you think. Imagine a case like this: You want to install package Matplotlib, when installing, it asks you to install Numpy, and Scipy, because Matplotlib need these Numpy and Scipy to work. They are called the dependencies of Matplotlib. For Numpy and Scipy, they may have their own dependencies. These require even more packages.
 
- Conda provide a solution for this situation: when you install package Matplotlib, it will automatically install all the dependencies like Numpy and Scipy. So you don’t have to install them one by one, manually. This can save you great amount of time.
 
- The other advantage of conda, is that conda can have multiple environments for different projects. As mentioned at the very beginning, it can have two separate environments of different versions of software.
Using conda environment on BioHPC

### Installing Packages Using Conda
>
>Conda is a package manager, which helps you find and install packages such as numpy or scipy. It also serves as an environment manager, and allows you to have multiple isolated environments for different projects on a single machine. Each environment has its own installation directories, that doesn’t share packages with other environments.
>
>For example, you need python 2.7 and Biopython 1.60 in project A, while you also work on another project B, which needs python 3.5 and Biopython 1.68. You can use conda to create two separate environments for each project, and you can switch between different versions of packages easily to run your project code.
{: .callout}

### Anaconda or Miniconda?  
- Anaconda includes both Python and conda, and additionally bundles a suite of other pre-installed packages geared toward scientific computing. Because of the size of this bundle, expect the installation to consume several gigabytes of disk space.

- Miniconda gives you the Python interpreter itself, along with a command-line tool called conda which operates as a cross-platform package manager geared toward Python packages, similar in spirit to the apt or yum tools that Linux users might be familiar with.

### Miniconda3

Miniconda is a package manager that simplifies the installation process. Please first install miniconda3 (installation instructions below), and then proceed to the installation of individual tools. 
​​
### Install Miniconda
Visit the [miniconda](https://docs.conda.io/en/latest/miniconda.html) page and get ready to download the installer of your choice/system.

There are several different env in this world.
***https://repo.anaconda.com/miniconda/***
![package]({{site.baseurl}}/fig/slide_package.png)

> ## MacOS
> 
>To install miniconda3, type:
>
>~~~
>$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
>$ bash Miniconda3-latest-MacOSX-x86_64.sh
>~~~
>{: .bash}
> Then, follow the instructions that you are prompted with on the screen to install Miniconda3. 
{: .solution}

> ## Ubuntu
> 
>To install miniconda3, type:
>
>~~~
>$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
>$ bash Miniconda3-latest-Linux-x86_64.sh
>~~~
>{: .bash}
> Then, follow the instructions that you are prompted with on the screen to install Miniconda3. 
{: .solution}

![conda1]({{site.baseurl}}/fig/conda_excute.png)
![conda2]({{site.baseurl}}/fig/conda_excute2.png)


### Reload your enviroment
#### Ubuntu
```bash
$ source ~/.bashrc
```

#### Ubuntu
```bash
$ source ~/.bash_profile
```
​
### Initialize Miniconda3

```bash
$ conda init
```



### Create new conda environment

#### Create a conda environment named test with latest anaconda package.
```bash 
$ conda create -n bch709 python=3
```
#### Alternatively you can specify python version
```bash
$ conda create -n snowdeer_env python=2.7.16
```

*Usually, the conda environment is installed in your home directory on computer, /home/\<your username\>. The newly created environment <test> will be installed in the directory /home/wyim/miniconda3/envs/test*
 


#### Use the environment you just created
Activate your environment:
```bash  
$ conda activate bch709
```
It will show your environment name at the beginning of the prompt.

![conda3]({{site.baseurl}}/fig/conda.png)

### Install packages in the conda environment

Install from default conda channel
You can search if your package is in the default source from Anaconda collection. Besides the 200 pre-built Anaconda packages, it contains over 600 extra scientific and analytic packages. All the dependencies will also be installed automatically.
``` 
$ conda search <package>
$ conda install <package>
```
### Install from conda-forge channel (example: hisat2)
Conda channels are the remote repository that conda takes to search or download the packages. If you want to install a package that is not in the default Anaconda channel, you can tell conda which channel containing the package, so that conda can find and install.
Conda-forge is a GitHub community-led conda channel, containing general packages which are not in the default Anaconda source. All the packages from conda-forge is listed at https://bioconda.github.io/conda-recipe_index.html

```bash
$ conda search hisat2
$ conda search -c bioconda hisat2
```

### Install from bioconda channel (example: hisat2)
Bioconda is another channel of conda, focusing on bioinformatics software. Instead of adding “-c” to search a channel only one time, “add channels” tells Conda to always search in this channel, so you don’t need to specify the channel every time. Remember to add channel in this order, so that Bioconda channel has the highest priority. Channel orders will be explained in next part.
```bash
 $ conda config --add channels conda-forge
 $ conda config --add channels defaults
 $ conda config --add channels r
 $ conda config --add channels bioconda
```
Adding channels will not generate any command line output.
Then, you can install Stringtie from the Bioconda channel
```bash   
 $ conda install hisat2
```
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/hisat2/README.html)


All the bioconda packages can be found here: https://bioconda.github.io/conda-recipe_index.html


### Install R and R packages
The Conda package manager is not limited to Python. R and R packages are well supported by a conda channel maintained by the developers of Conda. The R-essentials include all the popular R packages with all of their dependencies. The command below opens R channel by “-c r”, and then install the r-essentials using R channel.
```bash
$ conda install -c r r-essentials
```

#### Update R packages
```bash 
$ conda update -c r r-essentials
$ conda update -c r r-<package name>
```

### More conda commands:
#### See all available environments
You can check the list of all separate environments, and it will show * at your current environment. In the figure below, it shows root, since I’m not in any conda environment.
```bash   
$ conda env list
```

#### List all package installed
This will show all the packages and versions you’ve installed.
```bash
$ conda list
```
#### Update packages or conda itself
This will update to the newest version of one package, or conda itself.
update package
```bash
$ conda update <package>
```
update package in env
```bash
conda update  --name <ENV_name> <package>
```
update conda itself
```bash
$ conda update -n test --all
$ conda update -n bch709 --all
```

#### Uninstall package from the environment
```   
$ conda uninstall <package name>
```

#### Exit current environment:
You can exit, when you finish your work in the current environment.
```bash   
$ conda deactivate
```

#### Remove environment
```bash   
$ conda env remove --name bch709
```
When you finish your project, you might want to remove the environment. However, it is not recommended because you might want to update some work in this project in the future.

#### Enviroment export
```bash
conda env export  --name <ENVIRONMENT> --file <outputfilename>.yaml
```
#### Envrioment import 
```bash
conda env create --file <outputfilename>.yaml  
``` 


### HPC clusters
>This exercise mainly deals with using HPC clusters for large scale data (Next Generation Sequencing analysis, Genome annotation, evolutionary studies etc.). These clusters have several processors with large amounts of RAM (compared to typical desktop/laptop), which makes it ideal for running programs that are computationally intensive. The operating system of these clusters are primarily UNIX and are mainly operated via command line. All the commands that you have learned in the previous exercises can be used on HPC.
>
>Pronghorn High Performance Computing offers shared cluster computing infrastructure for researchers and students at UNR. Brief descriptions for the available resources can be found here: https://www.unr.edu/research-computing/hpc. To begin with, you need to request permission for accessing these resources either through your department or through your advisor. All workshop attendees will have their account setup on HPC class education cluster and they can use their UNR NetID and the password for logging-in. You should have already received a confirmation email about your account creation with instructions on how to connect to the cluster. In this exercise we will specifically teach you how to connect to a remote server (HPC), transfer files in and out of the server, and running programs by requesting resources.
>
>
You can log onto its front-end/job-submission system (pronghorn.rc.unr.edu) using your UNR NetID and password. Logging into HPC class requires an SSH client if you are >using Windows but Mac/Linux have these built into their OS. There are several available for download for the Windows platform.
>
>```
>ssh <YOURID>@pronghorn.rc.unr.edu
>```
{: .prereq}


There are a number of ways to transfer data to and from HPC clusters. Which you should use depends on several factors, including the ease of use for you personally, connection speed and bandwidth, and the size and number of files which you intend to transfer. Most common options include scp, rsync (command line) and SCP and SFTP clients (GUI). scp (secure copy) is a simple way of transferring files between two machines that use the SSH (Secure SHell) protocol. You may use scp to connect to any system where you have SSH (login) access. scp is available as a protocol choice in some graphical file transfer programs and also as a command line program on most Linux, UNIX, and Mac OS X systems. scp can copy single files, but will also recursively copy directory contents if given a directory name. scp can be used as follows:

- to a remote system from local
```
scp sourcefile username@pronghorn.rc.unr.edu:somedirectory/
```
- from a remote system to local
```
scp username@pronghorn.rc.unr.edu:somedirectory/sourcefile destinationfile
```
- recursive directory copy to a remote system from local
```
scp -r SourceDirectory/ username@pronghorn.rc.unr.edu:somedirectory/
```

***rsync*** is a fast and extraordinarily versatile file copying tool. It can synchronize file trees across local disks, directories or across a network

- Synchronize a local directory with the remote server directory
```
rsync -avhP path/to/SourceDirectory username@pronghorn.rc.unr.eduu:somedirectory/
```

- Synchronize a remote directory with the local directory
```
rsync -avhP username@hpronghorn.rc.unr.edu:SourceDirectory/ path/to/Destination/
```



>## HOME WORK
>1. Create your enviroment name bch709, install following tools
>- bwa
>- trinity
>- hisat2
>- samtools
>- bowtie2
>- canu
>- trim-galore 
>- multiqc 
>2. Then export your enviroment to bch.yaml (conda env export ......)
>3. copy contents of bch.yaml to [this site](https://docs.google.com/forms/d/e/1FAIpQLSe0faV4UHKEQ8CnqJ-iXOTAmKM2IxN6g7ZNSSmtcw5FuPwmWA/viewform?usp=sf_link).
{: .solution}

### Reference:

- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

