---
layout: page
title: Final Review
published: true
---



### K-mer spectrum
Analyzing k-mer frequencies in whole-genome sequencing data is becoming a common method for estimating genome size.  Genome size refers to the amount of haploid nuclear DNA of an organism and is typically measured in picograms or megabases (where 1 pg is equivalent to 978 Mb). Assuming that each k-mer is sequenced on average with C copies (k-mer coverage) and N denotes the number of genomic k-mers in the reads, the relationship N = C x (G–k + 1) allows to estimate Genome Size with G ≈ N/C as G ≫ k. Both C and N can be statistically inferred from a k-mer frequency histogram (or k-mer distribution in short), which summarizes how many distinct k-mers occur at a specific frequency within a given whole-genome sequencing data set. 


![kmer2]({{site.baseurl}}/fig/kmer2.png)
![genomescope]({{site.baseurl}}/fig/genomescope.png)




## Long read sequencing

Single Molecule, Real-Time (SMRT) Sequencing is the core technology powering our long-read sequencing platforms. This innovative approach was the first of its kind and is now a proven technology used in all fields of life science.


### Long Reads (PacBio)
With reads tens of kilobases in length you can readily assemble complete genomes and sequence full-length transcripts.

![pacbio]({{site.baseurl}}/fig/pacbio2.jpg)
A 20 kb size-selected human library using the SMRTbell Express Template Prep Kit 2.0 on a Sequel II System (2.0 Chemistry, Sequel II System Software v8.0, 30-hour movie).

SMRT Sequencing enables simultaneous collection of data from millions of wells using the natural process of DNA replication to sequence long fragments of native DNA or RNA.  
![pacbio]({{site.baseurl}}/fig/pacbio.png)

### PacBio Error
![pacbio]({{site.baseurl}}/fig/pacbio_error.png)
![hgap]({{site.baseurl}}/fig/HGAP.png)

### Repeat in Genome
![bandage]({{site.baseurl}}/fig/bandage.png)


## Genome Annotation

After the sections of DNA sequence have been assembled into a complete genome sequence we need to identify where the genes and key features are. We have our aligned and assembled genome sequence but how do we identify where the genes and other functional regions of the genome are located.

- Annotation involves marking where the genes start and stop in the DNA sequence and also where other relevant and interesting regions are in the sequence.

- Although genome annotation pipelines can differ from one another, for example, some elements can be manual while others have to be automated, they all share a core set of features.

- They are generally divided into two distinct phases: gene prediction and manual annotation.

### Gene prediction
There are two types of gene prediction: 
Ab initio – this technique relies on signals within the DNA sequence. It is an automated process whereby a computer is given instructions for finding genes in the sequence and is then left to find them. The computer looks for common sequences known to be found at the start and end of genes such as promoter sequences (where proteins bind that switch on genes), start codons (where the code for the gene product, RNA ?or protein, starts) and stop codons (where the code for the gene product ends).  

### Evidence-based
This technique relies on evidence beyond the DNA sequence. It involves gathering various pieces of genetic information from the transcript sequence (mRNA), and known protein sequences of the genome. With these pieces of evidence it is then possible to get an idea of the original DNA sequence by working backwards through transcription? and translation? (reverse transcription/translation). For example, if you have the protein sequence it is possible to work out the family of possible DNA sequences it could be derived from by working out which amino acids? make up the protein and then which combination of codons could code for those amino acids and so on, until you get to the DNA sequence. 
The information taken from these two prediction methods is then combined and lined up with the sequenced genome.

### Key Point
- Gene annotation is one of the core mechanisms through which we decipher the information that is contained in genome sequences.

- Gene annotation is complicated by the existence of 'transcriptional complexity', which includes extensive alternative splicing and transcriptional events outside of protein-coding genes.

- The annotation strategy for a given genome will depend on what it is hoped to achieve, as well as the resources available.

- The availability of next-generation data sets has transformed gene annotation pipelines in recent years, although their incorporation is rarely straightforward.

- Even human gene annotation is far from complete: transcripts are missing and existing models are truncated. Most importantly, 'functional annotation' — the description of what transcripts actually do — remains far from comprehensive.

- Efforts are now under way to integrate gene annotation pipelines with projects that seek to describe regulatory sequences, such as promoter and enhancer elements.

- Gene annotation is producing increasingly complex resources. This can present a challenge to usability, most notably in a clinical context, and annotation projects must find ways to resolve such problems.

### The core annotation workflows for different gene types.
![nrg.2016.119-f2]({{site.baseurl}}/fig/nrg.2016.119-f2.jpg)


## Reads Count

In RNA-seq gene expression data analysis, we come across various expression units such as RPM, RPKM, FPKM and raw reads counts. Most of the times it's difficult to understand basic underlying methodology to calculate these units from mapped sequence data.

I have seen a lot of post of such normalization questions and their confusion among readers. Hence, I attempted here to explain these units in the much simpler way (avoided complex mathematical expressions).
Note: To use biopython, you need to install it.

## <span style="color:#33a8ff">Why different normalized expression units?</span> ##

The expression units provide a digital measure of the abundance of transcripts. Normalized expression units are necessary to remove technical biases in sequenced data such as depth of sequencing (more sequencing depth produces more read count for gene expressed at same level) and gene length (differences in gene length generate unequal reads count for genes expressed at the same level; longer the gene more the read count).


## <span style="color:#33a8ff"> Gene expression units and calculation </span> ##

 **<span style="color:#060606">RPM (Reads per million mapped reads) </span>**

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;RPM&space;=&space;\frac{Number&space;\&space;of&space;\&space;reads&space;\&space;mapped&space;\&space;to&space;\&space;gene&space;\times&space;10^6}{Total&space;\&space;number&space;\&space;of&space;\&space;mapped&space;\&space;reads}" />

For example, You have sequenced one library with 5 million(M) reads. Among them, total 4 M matched to the genome sequence and 5000 reads matched to a given gene.

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;\large&space;RPM&space;=&space;\frac{5000&space;\times&space;10^6}{4&space;\times&space;10^6}&space;=&space;1250" />

Notes:

 - RPM does not consider the transcript length normalization.
 - RPM Suitable for sequencing protocols where reads are generated irrespective of gene length


**<span style="color:#060606">RPKM (Reads per kilo base per million mapped reads)</span>**

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;RPKM&space;=&space;\frac{Number&space;\&space;of&space;\&space;reads&space;\&space;mapped&space;\&space;to&space;\&space;gene&space;\times&space;10^3&space;\times&space;10^6}{Total&space;\&space;number&space;\&space;of&space;\&space;mapped&space;\&space;reads&space;\times&space;gene&space;\&space;length&space;\&space;in&space;\&space;bp}" />

 Here, 10^3 normalizes for gene length and 10^6 for sequencing depth factor.

FPKM (Fragments per kilo base per million mapped reads) is analogous to RPKM and used especially in paired-end RNA-seq experiments. In paired-end RNA-seq experiments, two (left and right) reads are sequenced from same DNA fragment. When we map paired-end data, both reads or only one read with high quality from a fragment can map to reference sequence. To avoid confusion or multiple counting, the fragments to which both or single read mapped is counted and represented for FPKM calculation.

For example, You have sequenced one library with 5 M reads. Among them, total 4 M matched to the genome sequence and 5000 reads matched to a given gene with a length of 2000 bp.

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;RPKM&space;=&space;\frac{5000&space;\times&space;10^3&space;\times&space;10^6}{4&space;\times&space;10^6&space;\times&space;2000}&space;=&space;625" />

Notes:

 - RPKM considers the gene length for normalization
 - RPKM is suitable for sequencing protocols where reads sequencing depends on gene length
 - Used in single-end RNA-seq experiments (FPKM for paired-end RNA-seq data)


**<span style="color:#060606">TPM (Transcript per million)</span>**

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;TPM&space;=&space;\frac{Number&space;\&space;of&space;\&space;reads&space;\&space;mapped&space;\&space;to&space;\&space;gene&space;\times&space;read&space;\&space;length&space;\times&space;10^6}{Total&space;\&space;number&space;\&space;of&space;\&space;transcripts&space;\&space;sampled&space;\times&space;gene&space;\&space;length&space;\&space;in&space;\&space;bp}" />

Here, read length refers to the average number of nucleotides mapped to a gene.

For example, You have sequenced one library with 5M 100 bp reads. Among them, total 4M matched to the genome sequence and 5000 reads matched to a given gene with a length of 2000 bp. There were 10K transcripts were sampled from a genome sequence i.e. reads mapped to 10K genes. Suppose all 100 bp mapped from 5000 reads.

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;TPM&space;=&space;\frac{5000&space;\times&space;100&space;\times&space;10^6}{10000&space;\times&space;2000}&space;=&space;25000" />

Notes:

 - TPM considers the gene length for normalization
 - TPM proposed as an alternative to RPKM due to inaccuracy in RPKM measurement (Wagner et al., 2012)
 - TPM is suitable for sequencing protocols where reads sequencing depends on gene length

**<span style="color:#060606">Relationship between RPKM and TPM,</span>**

 <img src="https://latex.codecogs.com/gif.latex?\bg_green&space;RPKM&space;=&space;TPM&space;\times&space;\frac{10^3&space;\times&space;total&space;\&space;number&space;\&space;of&space;\&space;transcripts&space;\&space;sampled}{read&space;\&space;length&space;\times&space;total&space;\&space;number&space;\&space;of&space;\&space;mapped&space;\&space;reads}" />


References:

 - Mortazavi A, Williams BA, McCue K, Schaeffer L, Wold B. Mapping and quantifying mammalian transcriptomes by RNA-Seq. Nature methods. 2008 Jul 1;5(7):621-8.
 - Wagner GP, Kin K, Lynch VJ. Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. Theory in biosciences. 2012 Dec 1;131(4):281-5.


**<span style="color:#33a8ff">How to cite?</span>**

Bedre, R. Bioinformatics data analysis and visualization toolkit. GitHub repository, <a href="https://github.com/reneshbedre/bioinfokit">https://github.com/reneshbedre/bioinfokit</a>


<span style="color:#9e9696">If you have any questions, comments or recommendations, please email me at 
<b>reneshbe@gmail.com</b></span>

<span style="color:#9e9696"><i> Last updated: November 25, 2019</i> </span>




## RNA-Seq with Genome
The alignment process consists of choosing an appropriate reference genome to map our reads against and performing the read alignment using one of several splice-aware alignment tools such as STAR or HISAT2. The choice of aligner is often a personal preference and also dependent on the computational resources that are available to you.

## BLAST (Basic Local Alignment Search Tool) 
BLAST is a popular program for searching biosequences against databases. BLAST was developed and is maintained by a group at the National Center for Biotechnology Information (NCBI). Salient characteristics of BLAST are:

### Local alignments
BLAST tries to find patches of regional similarity, rather than trying to find the best alignment between your entire query and an entire database sequence.
### Ungapped alignments
Alignments generated with BLAST do not contain gaps. BLAST's speed and statistical model depend on this, but in theory it reduces sensitivity. However, BLAST will report multiple local alignments between your query and a database sequence.

### Explicit statistical theory
BLAST is based on an explicit statistical theory developed by Samuel Karlin and Steven Altschul (PNAS 87:2284-2268. 1990) The original theory was later extended to cover multiple weak matches between query and database entry PNAS 90:5873. 1993).

CAUTION: the repetitive nature of many biological sequences (particularly naive translations of DNA/RNA) violates assumptions made in the Karlin & Altschul theory. While the P values provided by BLAST are a good rule-of-thumb for initial identification of promising matches, care should be taken to ensure that matches are not due simply to biased amino acid composition.

CAUTION: The databases are contaminated with numerous artifacts. The intelligent use of filters can reduce problems from these sources. Remember that the statistical theory only covers the likelihood of finding a match by chance under particular assumptions; it does not guarantee biological importance.


## Gene Ontology
Gene Ontology project is a major bioinformatics initiative Gene ontology is an annotation system The project provides the controlled and consistent vocabulary of terms and gene product annotations, i.e. terms occur only once, and there is a dictionary of allowed words
GO describes how gene products behave in a cellular context A consistent description of gene products attributes in terms of their associated biological processes, cellular components and molecular functions in a species-independent manner Each GO term consists of a unique alphanumerical identifier, a common name, synonyms (if applicable), and a definition Each term is assigned to one of the three ontologies Terms have a textual definition When a term has multiple meanings depending on species, the GO uses a "sensu" tag to differentiate among them (trichome differentiation (sensu Magnoliophyta) 
