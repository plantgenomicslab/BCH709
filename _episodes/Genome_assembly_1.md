---
layout: page
title: Genome assembly
published: true
---

>## HOMEWORK (10/31/19) due by 11/05
>### Check Genome Size by Illumina Reads
>
>```bash
>cd /data/gpfs/assoc/bch709/<YOUR_ID>/
>mkdir -p  Genome_assembly/Illumina
>cd Genome_assembly/Illumina
>```
>
>### Create Preprocessing Env
>```bash
>conda create -n preprocessing python=3
>conda install -c bioconda trim-galore jellyfish multiqc
>```
>
>### Reads Download
>```
>https://www.dropbox.com/s/ax38m9wra44lsgi/WGS_R1.fq.gz
>https://www.dropbox.com/s/kp7et2du5c2v385/WGS_R2.fq.gz
>```
>
>
>### Reads Trimming
>
>```bash
>#!/bin/bash
>#SBATCH --job-name=Trim
>#SBATCH --cpus-per-task=32
>#SBATCH --time=2:00:00
>#SBATCH --mem=100g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=<YOUR ID>@unr.edu
>#SBATCH -o trim.out # STDOUT
>#SBATCH -e trim.err # STDERR
>
>trim_galore --paired   --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 16  --max_n 40  -o trimmed_fastq <READ_R1> <READ_R2>
>fastqc <READ_R1> <READ_R2>
>multiqc . -n WGS_Illumina
>```
>
>### K-mer counting
>```bash
>mkdir kmer
>cd kmer
>```
>
>```bash
>#!/bin/bash
>#SBATCH --job-name=Trim
>#SBATCH --cpus-per-task=32
>#SBATCH --time=2:00:00
>#SBATCH --mem=100g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=<YOUR ID>@unr.edu
>#SBATCH -o trim.out # STDOUT
>#SBATCH -e trim.err # STDERR
>gunzip trimmed_fastq/*.gz
>jellyfish count -C -m 21 -s 1000000000 -t 10 trimmed_fastq/<trim_galore output>  -o reads.jf
>jellyfish histo -t 10 reads.jf > reads.histo
>```
>
>
>
>
>### Count Reads Number in file
>
>```bash
>echo $(zcat WGS_R1.fq.gz |wc -l)/4 | bc
>```
>
>### Advanced approach
>
>```bash
>for i in `ls *.fq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done
>For all gzip compressed fastq files, display the number of reads since 4 lines = 1 reads
>```
>
> ***Keep all results at Pronghorn***
>
{: .callout}


Meeting schedule [Link](https://docs.google.com/spreadsheets/d/1c4RzQle8AZPRdayYW5Ov3b16uWKMyUVXl8-iNnuCDSI/edit?usp=sharing)


## Genome assembly

![Genome assembly]({{site.baseurl}}/fig/genome_assembly.png)

After DNA sequencing is complete, the fragments of DNA that come out of the machine are all jumbled up. Like a jigsaw puzzle we need to take the pieces of the genome and put them back together.

### The human genome required a significant technological push
25✕ larger than any previously sequenced genome
1. Construction of genetic and physical maps of human and mouse genomes
2. Sequencing of yeast and worm genomes
3. Pilot projects to test the feasibility and cost effectiveness of large-scale sequencing
![Genome assembly]({{site.baseurl}}/fig/genome_history.png)


![human_genome]({{site.baseurl}}/fig/human_genome.png)

![human_genome2]({{site.baseurl}}/fig/human_genome2.png)

![geneticmap]({{site.baseurl}}/fig/geneticmap.png)


### What’s the challenge?
- The technology of DNA sequencing? is not 100 per cent accurate and therefore there are likely to be errors in the DNA sequence that is produced.  

- So, to account for the errors that could potentially occur, each base in the genome? is sequenced a number of times over, this is called coverage. For example, 30 times (30-fold) coverage means each base? is sequenced 30 times.  

- Effectively, the more times you sequence, or `read`, the same section of DNA, the more confidence you have that the final sequence is correct.
30- to 50-fold coverage is currently the standard used when sequencing human genomes to a high level of accuracy.  

- During the Human Genome Project? coverage was only between 5- and 10-fold and used a different sequencing technology to those used today.  

- Coverage has increased because of a few reasons:
    Although most current sequencing techniques are now faster than they were during the Human Genome Project, some sequencing technologies have a higher error rate.
    Some sequencing technologies deal with shorter reads of DNA which means that gaps are more likely to occur when the genome is assembled. Having a higher coverage reduces the likelihood of there being gaps in the final assembled sequence.
    It is also much cheaper to carry out sequencing to a higher coverage than it was at the time of the Human Genome Project.

- High coverage means that after sequencing DNA we have lots and lots of pieces of DNA sequence (reads).  

- To put this into perspective, once a human genome has been fully sequenced we have around 100 gigabases (100,000,000,000 bases) of sequence data.  

- Like the pieces of a jigsaw puzzle, these DNA reads are jumbled up so we need to piece them together and put them in the correct order to assemble the genome sequence.  

![assembly_depth]({{site.baseurl}}/fig/assembly_depth.png)
### Coverage calculation

Example: I know that the genome I am sequencing is 10 Mbases. I want a 50x coverage to do a good assembly. I am ordering 125 bp Illumina reads. How many reads do I need?

- (125xN)/10e+6=50
- N=(50x10e+6)/125=4e+6 (4 million reads)

### *De novo* assembly
- *De novo* sequencing is when the genome of an organism is sequenced for the first time.  

- In *de novo* assembly there is no existing reference genome sequence for that species to use as a template for the assembly of its genome sequence.  
- If you know that the new species? is very similar to another species that does have a reference genome, it is possible to assemble the sequence using a similar genome as a guide.  
- To help assemble a de novo sequence a physical gene? map can be developed before sequencing to highlight the **landmarks** so the scientists know where sections of DNA are located in relation to each other.  
- Producing a gene map can be an expensive process, so some assembly programmes rely on data consisting of a mix of single and paired-end reads 
- Single reads are where one end or the whole of a fragment of DNA is sequenced. These sequences can then be joined together by finding overlapping regions in the sequence to create the full DNA sequence.  

- Paired-end reads are where both ends of a fragment of DNA are sequenced. The distance between paired-end reads can be anywhere between 200 base pairs? and several thousand. The key advantage of paired-end reads is that scientists know how far apart the two ends are. This makes it easier to assemble them into a continuous DNA sequence. Paired-end reads are particularly useful when assembling a de novo sequence as they provide long-range information that you wouldn’t otherwise have in the absence of a gene map.  

- Assembly of a de novo sequence begins with a large number of short sections or “reads” of DNA.  
- These reads are compared to each other and those sharing the same DNA sequence are grouped together.  
- From here they are assembled into progressively larger sections to form long contiguous (together in sequence) sequences called “contigs”.  
- These contigs can then be grouped together with information taken from other technologies to provide clues for how to stitch the contigs together and roughly how far apart to place them, even if the sequence in between is still unknown. This is called “scaffolding”.  
- The assembly can be further refined by ordering the individual scaffolds into chromosomes?. A physical gene map is a useful tool for doing this.
- The resulting assembly is then fed on to the next stage of the process – annotation, which identifies where the genes and other features in the sequence start and stop.  
- The assembly of a genome is a computer-intensive job. It usually takes around 20 hours per gigabase of sequence for genome assembly programmes to stitch together an organism’s genome sequence from the reads of DNA sequence generated by the sequencing machines.  
- So, with the 100 gigabases of sequence data we have after sequencing a human genome, it will take 2,000 hours or around 83 days to assemble the complete sequence.  

![assembly_real.png]({{site.baseurl}}/fig/assembly_real.png)

Genome assembly refers to the process of taking a large number of short DNA sequences and putting them back together to create a representation of the original chromosomes from which the DNA originated. De novo genome assemblies assume no prior knowledge of the source DNA sequence length, layout or composition. In a genome sequencing project, the DNA of the target organism is broken up into millions of small pieces and read on a sequencing machine. These “reads” vary from 20 to 1000 nucleotide base pairs (bp) in length depending on the sequencing method used. Typically for Illumina type short read sequencing, reads of length 36 - 150 bp are produced. These reads can be either `single ended` as described above or `paired end`. A good summary of other types of DNA sequencing can be found below.

![Genome assembly]({{site.baseurl}}/fig/genome_assembly2.png)



Paired end reads are produced when the fragment size used in the sequencing process is much longer (typically 250 - 500 bp long) and the ends of the fragment are read in towards the middle. This produces two “paired” reads. One from the left hand end of a fragment and one from the right with a known separation distance between them. (The known separation distance is actually a distribution with a mean and standard deviation as not all original fragments are of the same length.) This extra information contained in the paired end reads can be useful for helping to tie pieces of sequence together during the assembly process.

The goal of a sequence assembler is to produce long contiguous pieces of sequence (contigs) from these reads. The contigs are sometimes then ordered and oriented in relation to one another to form scaffolds. The distances between pairs of a set of paired end reads is useful information for this purpose.


The mechanisms used by assembly software are varied but the most common type for short reads is assembly by de Bruijn graph. See this document for an explanation of the de Bruijn graph genome assembler "SPAdes.”

Genome assembly is a very difficult computational problem, made more difficult because many genomes contain large numbers of identical sequences, known as repeats. These repeats can be thousands of nucleotides long, and some occur in thousands of different locations, especially in the large genomes of plants and animals.

### Why do we want to assemble an organism’s DNA?
Determining the DNA sequence of an organism is useful in fundamental research into why and how they live, as well as in applied subjects. Because of the importance of DNA to living things, knowledge of a DNA sequence may be useful in practically any biological research. For example, in medicine it can be used to identify, diagnose and potentially develop treatments for genetic diseases. Similarly, research into pathogens may lead to treatments for contagious diseases.

### What do we need to do?

- Put the pieces together in the correct order to construct the complete genome sequence and identify any areas of interest.  
- This is done using processes called alignment and assembly:  
    Alignment is when the new DNA sequence is compared to existing DNA sequences to find any similarities or discrepancies between them and then arranged to show these features. Alignment is a vital part of assembly.

    Assembly involves taking a large number of DNA reads, looking for areas in which they overlap with each other and then gradually piecing together the ‘jigsaw’. It is an attempt to reconstruct the original genome. This is primarily carried out for de novo sequences?.


### The protocol in a nutshell:
- Obtain sequence read file(s) from sequencing machine(s).  
- Look at the reads - get an understanding of what you’ve got and what the quality is like.  
- Raw data cleanup/quality trimming if necessary.  
- Choose an appropriate assembly parameter set.  
- Assemble the data into contigs/scaffolds.  
- Examine the output of the assembly and assess assembly quality.  




![Genome assembly]({{site.baseurl}}/fig/genome_assembly3.png)

**Flowchart of de novo assembly protocol.**

## Long read sequencing

Single Molecule, Real-Time (SMRT) Sequencing is the core technology powering our long-read sequencing platforms. This innovative approach was the first of its kind and is now a proven technology used in all fields of life science.


### Long Reads (PacBio)
With reads tens of kilobases in length you can readily assemble complete genomes and sequence full-length transcripts.

![pacbio]({{site.baseurl}}/fig/pacbio2.jpg)
A 20 kb size-selected human library using the SMRTbell Express Template Prep Kit 2.0 on a Sequel II System (2.0 Chemistry, Sequel II System Software v8.0, 30-hour movie).

SMRT Sequencing enables simultaneous collection of data from millions of wells using the natural process of DNA replication to sequence long fragments of native DNA or RNA.  
![pacbio]({{site.baseurl}}/fig/pacbio.png)

## PacBio Error
![pacbio]({{site.baseurl}}/fig/pacbio_error.png)
![hgap]({{site.baseurl}}/fig/HGAP.png)


### Long Reads Assembly

![longreads_assembly]({{site.baseurl}}/fig/longreads_assembly.png)

### OLC Algorithm

![olc]({{site.baseurl}}/fig/olc.png) 

### OLC example
![olc1]({{site.baseurl}}/fig/olc1.png)  
![olc2]({{site.baseurl}}/fig/olc2.png) 
![olc3]({{site.baseurl}}/fig/olc3.png)  

### Long reads publications
- Mitsuhashi, Satomi et al. (2019) Long-read sequencing for rare human genetic diseases. Journal of human genetics
- Wenger, Aaron M et al. (2019) Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome. Nature biotechnology
- Eichler, Evan E et al. (2019) Genetic Variation, Comparative Genomics, and the Diagnosis of Disease. The New England journal of medicine
- Mantere, Tuomo et al. (2019) Long-Read Sequencing Emerging in Medical Genetics Frontiers in genetics
- Wang, Bo et al. (2019) Reviving the Transcriptome Studies: An Insight into the Emergence of Single-molecule Transcriptome Sequencing Frontiers in genetics
- Pollard, Martin O et al. (2018) Long reads: their purpose and place. Human molecular genetics
- Sedlazeck, Fritz J et al. (2018) Accurate detection of complex structural variations using single-molecule sequencing. Nature methods
- Ardui, Simon et al. (2018) Single molecule real-time (SMRT) sequencing comes of age: applications and utilities for medical diagnostics. Nucleic acids research
- Nakano, Kazuma et al. (2017) Advantages of genome sequencing by long-read sequencer using SMRT technology in medical area. Human cell
- Chaisson, Mark J P et al. (2015) Genetic variation and the de novo assembly of human genomes. Nature reviews. Genetics
- Rhoads, Anthony et al. (2015) PacBio sequencing and its applications. Genomics, proteomics & bioinformatics
- Berlin, Konstantin et al. (2015) Assembling large genomes with single-molecule sequencing and locality-sensitive hashing. Nature biotechnology
- Huddleston, John et al. (2014) Reconstructing complex regions of genomes using long-read sequencing technology. Genome research
- Koren, Sergey et al. (2013) Reducing assembly complexity of microbial genomes with single-molecule sequencing. Genome biology
- Travers, Kevin J et al. (2010) A flexible and efficient template format for circular consensus sequencing and SNP detection. Nucleic acids research
- Flusberg, Benjamin A et al. (2010) Direct detection of DNA methylation during single-molecule, real-time sequencing. Nature methods
- Eid, John et al. (2009) Real-time DNA sequencing from single polymerase molecules. Science

## Homework
1. Create Conda enviroment named "genomeassembly"  
2. Please install following software  
```
conda install -c bioconda spades canu pacbio_falcon samtools minimap2 multiqc 
conda install -c r r-ggplot2 r-stringr r-scales r-argparse
conda install -c conda-forge nano
```




## Repeat in Genome
![repeat]({{site.baseurl}}/fig/repeat.png)
![repeat2]({{site.baseurl}}/fig/repeat2.png)
![repeat3]({{site.baseurl}}/fig/repeat3.png)
![repeat4]({{site.baseurl}}/fig/repeat4.png)
![repeat5]({{site.baseurl}}/fig/repeat5.png)
![repeat6]({{site.baseurl}}/fig/repeat6.png)
Fig: E.coli genome (4.6Mbp)


### String Graph Algorithm
![stringgraph]({{site.baseurl}}/fig/string_graph.png)



## Hybrid Assembly
![hybrid_assembly]({{site.baseurl}}/fig/hybrid_assembly.png)

## Genome assembly complexity

### GC-content
- Regions of low or high GC-content have a lower coverage (Illumina, not PacBio)

### Secondary structure
- Regions that are tightly bound get less coverage

### Ploidy level
- On higher ploidy levels you potentially have more alleles present

### Size of organism
- Hard to extract enough DNA from small organisms

### Pooled individuals
- Increases the variability of the DNA (more alleles)

### Inhibiting compounds
- Lower coverage and shorter fragments

### Presence of additional genomes/contamination
- Lower coverage of what you actually are interested in, potentially chimeric assemblies

### Heterozygosity
![hetero]({{site.baseurl}}/fig/hetero.jpg)

### Haplotype
![haplotype]({{site.baseurl}}/fig/SNPS.png)


## Contiguity
- Fewer contigs
- Longer contigs

### Metrics
- Number of contigs
- Average contig length
- Median contig length
- Maximum contig length
- “N50”, “NG50”, “D50”

## N50
![N50]({{site.baseurl}}/fig/N50.png)
![N502]({{site.baseurl}}/fig/N502.jpg)

## N50 example

|Contig Length|Cumulative Sum|  
| --- | --- |  
|100|100|  
|200|300|  
|230|530|  
|400|930|  
|750|1680|  
|852|2532|  
|950|3482|  
|990|4472|  
|1020|5492|  
|1278|6770|  
|1280|8050|  
|1290|9340|    

## Completeness : Total size
Proportion of the original genome represented by the assembly
Can be between 0 and 1

![size]({{site.baseurl}}/fig/size.png)

## Completeness: core genes

- Proportion of coding sequences can be estimated based on known core genes thought to be present in a wide variety of organisms.  

- Assumes that the proportion of assembled genes is equal to the proportion of assembled core genes.
![coregene]({{site.baseurl}}/fig/coregenes.png)

## Correctness
Proportion of the assembly that is free from errors  
Errors include
1. Mis-joins
2. Repeat compressions
3. Unnecessary duplications
4. Indels / SNPs caused by assembler


## Optical mapping
![bionano3]({{site.baseurl}}/fig/bionano3.jpg)
![bionano2]({{site.baseurl}}/fig/bionano2.jpg)
![bionano]({{site.baseurl}}/fig/bionano.jpg)

## 10x Genomics
![10x]({{site.baseurl}}/fig/10x.jpg)

## Assembly results
![assembly_results]({{site.baseurl}}/fig/assembly_results.png)

## Dotplot
![genome_plot2]({{site.baseurl}}/fig/dotplot2.png)
![genome_plot]({{site.baseurl}}/fig/genome_plot.png)


## Why assemble genomes
- Produce a reference for new species
- Genomic variation can be studied by comparing against the reference
- A wide variety of molecular biological tools become available or more effective
- Coding sequences can be studied in the context of their related non-coding (eg regulatory) sequences
- High level genome structure (number, arrangement of genes and repeats) can be studied

## Assembly is very challenging (“impossible”) because
- sequencing bias under represents certain regions
- Reads are short relative to genome size
- Repeats create tangled hubs in the assembly graph
- Sequencing errors cause detours and bubbles in the assembly graph

## Flow Cytometry
![flowcytometry]({{site.baseurl}}/fig/flowcytometry.png)


## K-mer spectrum
![kmer2]({{site.baseurl}}/fig/kmer2.png)
![genomescope]({{site.baseurl}}/fig/genomescope.png)


## Check Genome Size by Illumina Reads

```bash
cd /data/gpfs/assoc/bch709/<YOUR_ID>/
mkdir Genome_assembly/Illumina
cd Genome_assembly/Illumina
```

### Create Preprocessing Env
```bash
conda create -n preprocessing python=3
conda install -c bioconda trim-galore jellyfish multiqc 
```

### Reads Download
```
https://www.dropbox.com/s/ax38m9wra44lsgi/WGS_R1.fq.gz
https://www.dropbox.com/s/kp7et2du5c2v385/WGS_R2.fq.gz
```


### Reads Trimming

```bash
#!/bin/bash
#SBATCH --job-name=Trim
#SBATCH --cpus-per-task=32
#SBATCH --time=2:00:00
#SBATCH --mem=100g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o trim.out # STDOUT
#SBATCH -e trim.err # STDERR

trim_galore --paired   --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 16  --max_n 40  -o trimmed_fastq <READ_R1> <READ_R2>
fastqc <READ_R1> <READ_R2>
multiqc . -n WGS_Illumina
```

### K-mer counting
```bash
mkdir kmer
cd kmer
```

```bash
#!/bin/bash
#SBATCH --job-name=Trim
#SBATCH --cpus-per-task=32
#SBATCH --time=2:00:00
#SBATCH --mem=100g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o trim.out # STDOUT
#SBATCH -e trim.err # STDERR

jellyfish count -C -m 21 -s 1000000000 -t 10 <trim_galore output>  -o reads.jf
jellyfish histo -t 10 reads.jf > reads.histo
```




### Count Reads Number in file

```bash
echo $(zcat WGS_R1.fq.gz |wc -l)/4 | bc
```

### Advanced approach

```bash
for i in `ls *.fq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done
For all gzip compressed fastq files, display the number of reads since 4 lines = 1 reads
```


### Upload to GenomeScope
http://qb.cshl.edu/genomescope/genomescope2.0






### Genome assembly Spades
```bash
mkdir Spades
cd Spades
conda activate genomeassembly

```

```bash
#!/bin/bash
#SBATCH --job-name=Spades
#SBATCH --cpus-per-task=64
#SBATCH --time=2:00:00
#SBATCH --mem=140g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o Spades.out # STDOUT
#SBATCH -e Spades.err # STDERR

spades.py -k 21,33,55,77 --careful -1 <trim_galore output> -2 <trim_galore output> -o spades_output --memory 140 --threads 64
```
