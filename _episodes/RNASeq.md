---
layout: page
title: RNA-Seq
published: true
---

{% include gh_variables.html %}

>## Paper reading
>Please read this paper
>https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8
{: .prereq}


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

## Seven stages to data science
1. Define the question of interest
2. Get the data
3. Clean the data
4. Explore the data
5. Fit statistical models
6. Communicate the results
7. Make your analysis reproducible

## What do we need to prepare ?
### Sample Information
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




### RNA Information: 
RNAs come in all shapes and sizes. Some of the key properties to report are:
a. Total RNA, Poly-A(+) RNA, Poly-A(-) RNA
b. Size of the RNA fraction: we typically have a + 200 and – 200 cutoff, but there is a wide
range, i.e. microRNA-sized, etc…
c. If the RNA was treated with Ribosomal RNA depletion kits (RiboMinus, RiboZero): please
note the kit used.

### Protocols: 
There are several methods used to isolate RNAs with that work fine for the purposes of RNA-Seq. For all the ENCODE libraries that we make, we provide a document that lists in detail:
a. The RNA isolation methods,
b. Methods of size selections
c. Methods of rRNA removal
d. Methods of oligo-dT selections
e. Methods of DNAse I treatments

### Experimental Design
- Balanced design
- Technical replicates not necessary (Marioni et al., 2008)
- Biological replicates: 6 - 12 (Schurch et al., 2016)
- Power analysis

>## Reading materials
>[Paul L. Auer and R. W. Doerge "Statistical Design and Analysis of RNA Sequencing Data" Genetics June 1, 2010 vol. 185 no.2 405-416](https://www.genetics.org/content/185/2/405)  
>[Busby, Michele A., et al. "Scotty: a web tool for designing RNA-Seq experiments to measure differential gene expression." Bioinformatics 29.5 (2013): 656-657](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582267/)  
>[Marioni, John C., et al. "RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays." Genome research (2008)](https://genome.cshlp.org/content/18/9/1509.full.html)  
>[Schurch, Nicholas J., et al. "How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?." Rna (2016)](https://rnajournal.cshlp.org/content/22/6/839.long)  
>[Zhao, Shilin, et al. "RnaSeqSampleSize: real data based sample size estimation for RNA sequencing." BMC bioinformatics 19.1 (2018): 191](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5)  
{: .prereq}


### Replicate number
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


### RNA extraction
- Sample processing and storage
- Total RNA/mRNA/small RNA
- DNAse treatment
- Quantity & quality
- RIN values (Strong effect)
- Batch effect
- Extraction method bias (GC bias)

>## Reading materials
>Romero, Irene Gallego, et al. "RNA-seq: impact of RNA degradation on transcript quantification." BMC biology 12.1 (2014): 42
>Kim, Young-Kook, et al. "Short structured RNAs with low GC content are selectively lost during extraction from a small number of cells." Molecular cell 46.6 (2012): 893-89500481-9).
{: .prereq}

### RNA Quantification and Quality Control: When working with bulk samples, throughout the various
steps we periodically assess the quality and quantity of the RNA. This is typically done on a
BioAnalyzer. Points to check are:
a. Total RNA
b. After oligo-dT size selections
c. After rRNA-depletions
d. After library construction

### Library prep
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


### Sequencing: 
There are several sequencing platforms and technologies out there being used. It is important to provide the following pieces of information:
a. Platform: Illumina, PacBio, Oxford Nanopore, etc…
b. Format: Single-end, Pair-end,
c. Read Length: 101 bases, 125 bases, etc…
d. Unusual barcode placement and sequence: Some protocols introduce barcodes in noncustomary places. If you are going to deliver a FASTQ file that will contain the barcode
sequences in it or other molecular markers – you will need to report both the position in the
read(s) where they are and their sequence(s).
e. Please provide the sequence of any custom primers that were used to sequence the library

![RNA library]({{{site.baseurl}}/fig/sequencing.png)



### Sequencing depth.
The amount of sequencing needed for a given sample is determined by the goals of the experiment and
the nature of the RNA sample. Experiments whose purpose is to evaluate the similarity between the
transcriptional profiles of two polyA+ samples may require only modest depths of sequencing.
Experiments whose purpose is discovery of novel transcribed elements and strong quantification of
known transcript isoforms requires more extensive sequencing.  
• Each Long RNA-Seq library must have a minimum of 30 million aligned reads/mate-pairs.  
• Each RAMPAGE library must have a minimum of 20 million aligned reads/mate-pairs.  
• Each small RNA-Seq library must have a minimum of 30 million aligned reads/mate-pairs.  

### Quantitative Standards (spike-ins).
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



![RNA Sequencing tool]({{{site.baseurl}}/fig/rnasoftware.png)


### QC FAIL?
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


```bash
pwd

cd ~/

mkdir -p bch709/rnaseq

cd bch709/rnaseq/

pwd

wget https://nevada.box.com/shared/static/g9dc1dn6h23u6tktv0uttdw2q1ad3xuz.gz -O pair1.fastq.gz

wget https://nevada.box.com/shared/static/y320e7atipagawwvl4plkclmeh86a1vg.gz -O pair2.fastq.gz
ls -algh

zcat pair2.fastq.gz | head
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

![Fastq_file]({{{site.baseurl}}/fig/fastq.png)

![basequality]({{{site.baseurl}}/fig/basequality.png)




### Quality Scores
Quality scores are a way to assign confidence to a particular base within a read. Some sequencers have their own proprietary quality encoding but most have adopted Phred-33 encoding. Each quality score represents the probability of an incorrect basecall at that position.

### Phred Quality Score Encoding
Quality scores started as numbers (0-40) but have since changed to an ASCII encoding to reduce filesize and make working with this format a bit easier, however they still hold the same information. ASCII codes are assigned based on the formula found below. This table can serve as a lookup as you progress through your analysis.

### Quality Score Interpretation
Once you know what each quality score represents you can then use this chart to understand the confidence in a particular base.


![FASTQ quality]({{{site.baseurl}}/fig/quality.png)


### Conda enviroment

```bash
conda create -n rnaseq_test


conda activate rnaseq_test


```

### Reads QC
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

```bash
conda search fastqc
conda install -c conda-forge -c bioconda fastqc
```


### Run fastqc
```bash
fastqc --help
fastqc -t <YOUR CPU COUNT> pair1.fastq.gz  pair2.fastq.gz

```
### Download folder setting
*in your desktop*

#### Windows
```bash
mkdir /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop 
ln -s /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop ~/bch709
```

#### MacOS
```bash
mkdir ~/Desktop/BCH709_Desktop 
ln -s ~/Desktop/BCH709_Desktop ~/bch709
```

### Download results
*in your desktop*
```bash
scp <YOURID>@pronghorn.rc.unr.edu:~/bch709/rnaseq/*_fastqc.html ~/bch709
cd  ~/bch709
```

#### open current location in Windows
```bash
explorer.exe .
```

#### MacOS open curretn directory
```bash
open .
```

### Practice downloading
```bash
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/10x.jpg
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/Global_vs_Local.png
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/SNPS.png
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/Multi_QC_Results.png
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/HMM.jpeg
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/hicmovie.gif
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/fig/difference-between-global-and-local.html
https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/Homo_sapiens.GRCh38.cds.all.fa.gz
```

### How to make a report?
![MultiQC]({{{site.baseurl}}/fig/multiqc.png)
[MultiQC](https://multiqc.info/)
```bash
conda search multiqc 
conda install -c conda-forge -c bioconda multiqc
multiqc --help
multiqc .
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

### Install Trim Galore
```bash
conda install -c bioconda -c conda-forge trim-galore
```

### Run trimming
```bash
trim_galore --help

trim_galore --paired   --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 2  --max_n 40  --gzip -o trim pair1.fastq.gz pair2.fastq.gz --fastqc

multiqc .
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

### Install HISAT2 (graph FM index, spin off version Burrows-Wheeler Transform)
```bash
conda install -c conda-forge -c bioconda hisat2
```
### Download reference sequence
```bash
wget https://nevada.box.com/shared/static/5v14j6gjt16c7k5d42b7g51csjmmj16l.fasta -O bch709.fasta

```
### HISAT2 indexing
```bash
hisat2-build --help
hisat2-build bch709.fasta bch709
```

### HISAT2 mapping
```bash
hisat2 -x bch709 --threads <YOUR CPU COUNT> -1 trim/paired1_val_1.fq.gz -2 trim/paired2_val_2.fq.gz  -S align.sam 2> summarymetrics.txt
```

### SAM file format

Check result
```bash
head align.sam
```
### SAM file

```
<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [<TAG>:<VTYPE>:<VALUE> [...]]
```
![sam1]({{{site.baseurl}}/fig/sam1.png)  

### SAM flag
![flag]({{{site.baseurl}}/fig/flag.jpg)  

### Demical to binary
```bash
echo 'obase=163;10' | bc
```
If you don't have bc, please install through conda
```bash
conda install -c conda-forge bc
```



### SAM tag
There are a bunch of predefined tags, please see the SAM manual for more information. For the tags used in this example:

Any tags that start with X? are reserved fields for end users: XT:A:M, XN:i:2, XM:i:0, XO:i:0, XG:i:0

![samtag]({{{site.baseurl}}/fig/samtag.png)

### More information is below.
http://samtools.github.io/hts-specs/


### SAMtools
SAM (Sequence Alignment/Map) format is a generic format for storing large nucleotide sequence alignments. SAM aims to be a format that:

- Is flexible enough to store all the alignment information generated by various alignment programs;
- Is simple enough to be easily generated by alignment programs or converted from existing alignment formats;
- Is compact in file size;
- Allows most of operations on the alignment to work on a stream without loading the whole alignment into memory;
- Allows the file to be indexed by genomic position to efficiently retrieve all reads aligning to a locus.

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. http://samtools.sourceforge.net/
```bash
samtools view -Sb align.sam > align.bam

samtools sort align.bam  -o align_sort.bam

samtools index align_sort.bam
```

### BAM file
A BAM file (.bam) is the binary version of a SAM file. A SAM file (.sam) is a tab-delimited text file that contains sequence alignment data. 

|SAM/BAM|size|
|-----|----|
|align.sam | 903M|
|align.bam | 166M|


### Alignment visualization
```bash
samtools tview align_sort.bam bch709.fasta
```
![tview]({{{site.baseurl}}/fig/tview.png)  

[IGV](https://software.broadinstitute.org/software/igv/)  
SAMTOOLS Tview  
[Tablet](https://ics.hutton.ac.uk/tablet/)    


### Alignment QC
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


### Quantification • Counts
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

```bash
conda install -c conda-forge -c bioconda subread
conda install -c conda-forge -c bioconda rsem
```

```bash
wget https://www.dropbox.com/s/e9dvdkrl9dta4qg/bch709.gtf

featureCounts -p  -a bch709.gtf align_sort.bam -o counts.txt
```


### Quantification method
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

>## Reference
>Fu, Yu, et al. "Elimination of PCR duplicates in RNA-seq and small RNA-seq using unique molecular identifiers." BMC genomics 19.1 (2018): 531
>Parekh, Swati, et al. "The impact of amplification on differential expression analyses by RNA-seq." Scientific reports 6 (2016): 25533
>Klepikova, Anna V., et al. "Effect of method of deduplication on estimation of differential gene expression using RNA-seq." PeerJ 5 (2017): e3091
{: .challenge}

### Differential expression
DESeq2
edgeR (Neg-binom > GLM > Test)
Limma-Voom (Neg-binom > Voom-transform > LM > Test)

```bash
conda install -c conda-forge -c bioconda bioconductor-deseq2
```

### Functional analysis • GO
Gene enrichment analysis (Hypergeometric test)
Gene set enrichment analysis (GSEA)
Gene ontology / Reactome databases

>## Reading material
>Conesa, Ana, et al. "A survey of best practices for RNA-seq data analysis." Genome biology 17.1 (2016): 13
{: .challenge}

### Please revisit
Introduction to R
https://learn.datacamp.com/courses/free-introduction-to-r




## Transcriptome Assembly
***De novo*** assembly
![denovo]({{{site.baseurl}}/fig/denovo.png)

### Assembly?
![assembly]({{{site.baseurl}}/fig/assembly.png)  
![assembly1]({{{site.baseurl}}/fig/assembly1.png)  

### Sequencing coverage
![coverage]({{{site.baseurl}}/fig/coverage.png)  
![averagecoverage]({{{site.baseurl}}/fig/averagecoverage.png)  

### Assembly law
![assemblylaw]({{{site.baseurl}}/fig/assemblylaw.png)  
![assemblylaw1]({{{site.baseurl}}/fig/assemblylaw1.png)  

### Overlap graph
![overlapgraph]({{{site.baseurl}}/fig/overlapgraph.png)  
![greedy1]({{{site.baseurl}}/fig/greedy1.png)  
![greedy2]({{{site.baseurl}}/fig/greedy2.png)  
![greedy3]({{{site.baseurl}}/fig/greedy3.png)  
![greedy4]({{{site.baseurl}}/fig/greedy4.png)  
![greedy5]({{{site.baseurl}}/fig/greedy5.png)  
![greedy6]({{{site.baseurl}}/fig/greedy6.png)  
![greedy7]({{{site.baseurl}}/fig/greedy7.png)  
![greedy8]({{{site.baseurl}}/fig/greedy8.png)  
![greedy9]({{{site.baseurl}}/fig/greedy9.png)  

### K-mer: substring of length k
k-mers are subsequences of length ***k*** contained within a biological sequence.
![kmer]({{{site.baseurl}}/fig/kmer.png)
### De bruijn
![hamiltonian_Eulerian]({{{site.baseurl}}/fig/hamiltonian_Eulerian.png)
![DBG]({{{site.baseurl}}/fig/DBG.png)  
![dbg1]({{{site.baseurl}}/fig/dbg2.png)  


### OLC vs De bruijn
![assemblyalgorithm]({{{site.baseurl}}/fig/assemblyalgorithm.png)  
![olcdbg]({{{site.baseurl}}/fig/olcdbg.png)  
![realdbg]({{{site.baseurl}}/fig/realdbg.png)  
![realdbg2]({{{site.baseurl}}/fig/realdbg2.png)  
![realdbg3]({{{site.baseurl}}/fig/realdbg3.png)  

### Transcriptome assembler
![trinity]({{{site.baseurl}}/fig/trinity.png)  


#### Transcriptome assembly error

![transcriptome_error]({{{site.baseurl}}/fig/transcriptome_error.png)  

# Running Trinity

Trinity is run via the script: 'Trinity' found in the base installation directory.

Usage info is as follows:

```bash
conda create -n transcriptome_assembly

conda activate transcriptome_assembly

conda install -c bioconda trinity

```

     ###############################################################################
     #
     #     ______  ____   ____  ____   ____  ______  __ __
     #    |      ||    \ |    ||    \ |    ||      ||  |  |
     #    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
     #    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
     #      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
     #      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
     #      |__|  |__|\_||____||__|__||____|  |__|  |____/
     #
     ###############################################################################
     #
     # Required:
     #
     #  --seqType <string>      :type of reads: ( fa, or fq )
     #
     #  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
     #                            provided in Gb of RAM, ie.  '--max_memory 10G'
     #
     #  If paired reads:
     #      --left  <string>    :left reads, one or more file names (separated by commas, not spaces)
     #      --right <string>    :right reads, one or more file names (separated by commas, not spaces)
     #
     #  Or, if unpaired reads:
     #      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
     #
     #  Or,
     #      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
     #                                   ex.
     #                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
     #                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
     #                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
     #                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
     #
     #                      # if single-end instead of paired-end, then leave the 4th column above empty.
     #
     ####################################
     ##  Misc:  #########################
     #
     #  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
     #                                   if paired: RF or FR,
     #                                   if single: F or R.   (dUTP method = RF)
     #                                   See web documentation.
     #
     #  --CPU <int>                     :number of CPUs to use, default: 2
     #  --min_contig_length <int>       :minimum assembled contig length to report
     #                                   (def=200)
     #
     #  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
     #
     #  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
     #                                   (see genome-guided param section under --show_full_usage_info)
     #
     #  --jaccard_clip                  :option, set if you have paired reads and
     #                                   you expect high gene density with UTR
     #                                   overlap (use FASTQ input file format
     #                                   for reads).
     #                                   (note: jaccard_clip is an expensive
     #                                   operation, so avoid using it unless
     #                                   necessary due to finding excessive fusion
     #                                   transcripts w/o it.)
     #
     #  --trimmomatic                   :run Trimmomatic to quality trim reads
     #                                        see '--quality_trimming_params' under full usage info for tailored settings.
     #
     #
     #  --no_normalize_reads            :Do *not* run in silico normalization of reads. Defaults to max. read coverage of 50.
     #                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
     #                                       (note, as of Sept 21, 2016, normalization is on by default)
     #
     #
     #
     #  --output <string>               :name of directory for output (will be
     #                                   created if it doesn't already exist)
     #                                   default( your current working directory: "/Users/bhaas/GITHUB/trinityrnaseq/trinity_out_dir"
     #                                    note: must include 'trinity' in the name as a safety precaution! )
     #
     #  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
     #
     #  --cite                          :show the Trinity literature citation
     #
     #  --version                       :reports Trinity version (BLEEDING_EDGE) and exits.
     #
     #  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
     #
     #
     ###############################################################################
     #
     #  *Note, a typical Trinity command might be:
     #
     #        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
     #
     #
     #    and for Genome-guided Trinity:
     #
     #        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
     #                --genome_guided_max_intron 10000 --CPU 6
     #
     #     see: /Users/bhaas/GITHUB/trinityrnaseq/sample_data/test_Trinity_Assembly/
     #          for sample data and 'runMe.sh' for example Trinity execution
     #
     #     For more details, visit: http://trinityrnaseq.github.io
     #
     ###############################################################################


<a name='strand_specific_assembly'>
Trinity performs best with strand-specific data, in which case sense and antisense transcripts can be resolved.  For protocols on strand-specific RNA-Seq, see: [Borodina T, Adjaye J, Sultan M. A strand-specific library preparation protocol for RNA sequencing. Methods Enzymol. 2011;500:79-98. PubMed PMID: 21943893](http://www.ncbi.nlm.nih.gov/pubmed/21943893).


If you have strand-specific data, specify the library type.  There are four library types:

- Paired reads:
    * *RF*: first read (/1) of fragment pair is sequenced as anti-sense (reverse(*R*)), and second read (/2) is in the sense strand (forward(*F*)); typical of the dUTP/UDG sequencing method.
    * *FR*: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the antisense strand (reverse)

- Unpaired (single) reads:
    * *F*: the single read is in the sense (forward) orientation
    * *R*: the single read is in the antisense (reverse) orientation

By setting the *--SS_lib_type* parameter to one of the above, you are indicating that the reads are strand-specific.  By default, reads are treated as not strand-specific.

![strand specific specification](https://raw.githubusercontent.com/wiki/trinityrnaseq/trinityrnaseq/images/strand_specificity.jpg)

Other important considerations:

- Whether you use Fastq or Fasta formatted input files, be sure to keep the reads oriented as they are reported by Illumina, if the data are strand-specific. This is because, Trinity will properly orient the sequences according to the specified library type.  If the data are not strand-specific, no worries because the reads will be parsed in both orientations.

- If you have both paired and unpaired data, and the data are NOT strand-specific, you can combine the unpaired data with the left reads of the paired fragments.  Be sure that the unpaired reads have a /1 as a suffix to the accession value similarly to the left fragment reads.  The right fragment reads should all have /2 as the accession suffix.  Then, run Trinity using the --left and --right parameters as if all the data were paired.

- If you have multiple paired-end library fragment sizes, set the '--group_pairs_distance' according to the larger insert library.  Pairings that exceed that distance will be treated as if they were unpaired by the Butterfly process.

- by setting the '--CPU option', you are indicating the maximum number of threads to be used by processes within Trinity. Note that Inchworm alone will be internally capped at 6 threads, since performance will not improve for this step beyond that setting)


<a name="typical_trinity_command_line"></a>
## Typical Trinity Command Line

A typical Trinity command for assembling non-strand-specific RNA-seq data would be like so, running the entire process on a single high-memory server (aim for \~1G RAM per \~1M \~76 base Illumina paired reads, but often *much* less memory is required):

Run Trinity like so:

     Trinity --seqType fq --max_memory 50G \
             --left reads_1.fq.gz  --right reads_2.fq.gz --CPU 6

If you have multiple sets of fastq files, such as corresponding to multiple tissue types or conditions, etc., you can indicate them to Trinity like so:

     Trinity --seqType fq --max_memory 50G  \
             --left condA_1.fq.gz,condB_1.fq.gz,condC_1.fq.gz \
             --right condA_2.fq.gz,condB_2.fq.gz,condC_2.fq.gz \
             --CPU 6

or better yet, create a 'samples.txt' file that describes the data like so:


     #      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
     #                                   ex.
     #                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
     #                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
     #                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
     #                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq

This same samples file can then be used later on with other downstream analysis steps, including expression quantification and differential expression analysis.


>note that fastq files can be gzip-compressed as shown above, in which case they should require a '.gz' extension.



<a name="typical_options"></a>


### Run trimming
```
$ trim_galore --paired   --three_prime_clip_R1 10 --three_prime_clip_R2 10 --cores 2  --max_n 40  --gzip -o trimmed_fastq fastq/WTD1_1.fastq.gz fastq/WTD1_2.fastq.gz 
.
.
.
.
.

$ multiqc . -n rnaseq_data
```


***PLEASE CHECK YOUR MULTIQC***

### Trinity run
```
htop
mkdir -p /data/gpfs/assoc/bch709-2/YOURID/rnaseq/transcriptome_assembly

cd /data/gpfs/assoc/bch709-2/YOURID/rnaseq/transcriptome_assembly

Trinity --seqType fq --max_memory 6G --left ..... --right ....

```

Paper need to read
https://academic.oup.com/gigascience/article/8/5/giz039/5488105
https://www.nature.com/articles/nbt.1883
https://www.nature.com/articles/nprot.2013.084


>## Assignment
>1. de Bruijn graph construction (10 pts)
> - Draw (by hand) the de Bruijn graph for the following reads using k=3 (assume all reads are from the forward strand, no sequencing errors)  
> - Please provide the assembled sequence from this reads.
>ATG  
>AGT  
>CAT  
>GTA  
>GTT  
>TGT  
>TAG  
>TTA  
>TAC  
>
>2. Trimming Practice
> - Make `Trim` folder under `/data/gpfs/assoc/bch709-2/YOURID/BCH709_assignment`
> - Change directory to the folder `Trim`
> - Copy all fastq.gz files in `/data/gpfs/assoc/bch709-2/Course_materials/RNASeq_raw_fastq` to `Trim` folder
> - Run `Trim-Galore` and process `MultiQC`
> - Generate trim report output from `MultiQC` and upload `html` file to WebCanvas.
> *** IF YOU USE "LOOP" FOR THIS JOB, YOU WILL GET ADDITIONAL 10 POINTS ***
{: .solution}

### Reference:

- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/
