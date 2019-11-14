---
layout: page
title: Genome assembly
published: true
---
# Meeting schedule
Please feel free to fill up below link to review our course.
Meeting schedule [Link](https://docs.google.com/spreadsheets/d/1c4RzQle8AZPRdayYW5Ov3b16uWKMyUVXl8-iNnuCDSI/edit?usp=sharing)


>## Homework 11/09/2019 Due by 11/14/2019
>
>## PacBio reads preprocessing env
>```bash
>conda activate preprocessing
>```
>### Reads download
>```bash
>cd /data/gpfs/assoc/bch709/<YOURID>/genomeassembly/
>cd PacBio
>```
>** if this is not there, please use 'mkdir'
>
>### Download below reads
>```bash
>https://www.dropbox.com/s/1e4a4jpp3eszt6u/BCH709_Pacbio_02.fastq.gz
>https://www.dropbox.com/s/au2528hpm8vvr4c/BCH709_Pacbio_01.fastq.gz
>```
>
>### Check PacBio reads statistics *Submit below job*
>```bash
>#!/bin/bash
>#SBATCH --job-name=PacBio_stat
>#SBATCH --cpus-per-task=64
>#SBATCH --time=2:00:00
>#SBATCH --mem=140g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=<YOUR ID>@unr.edu
>#SBATCH -o PacBio_stat.out # STDOUT
>#SBATCH -e PacBio_stat.err # STDERR
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>
>NanoStat --fastq BCH709_Pacbio_02.fastq.gz > BCH709_Pacbio_02.stat.txt
>NanoPlot -t 2 --fastq  BCH709_Pacbio_02.fastq.gz --maxlength 25000 --plots hex dot pauvre -o pacbio_stat
>
>
>NanoStat --fastq BCH709_Pacbio_01.fastq.gz > BCH709_Pacbio_01.stat.txt
>NanoPlot -t 2 --fastq  BCH709_Pacbio_01.fastq.gz --maxlength 25000 --plots hex dot pauvre -o pacbio_stat
>```
>
>
>### Transfer your result
>```bash
>*.png *.html *.txt
>```
>
>
>
>
>### Genome assembly Spades
>```bash
>mkdir Spades_pacbio
>cd Spades_pacbio
>conda activate genomeassembly
>
>```
>### Submit below job
>
>```bash
>#!/bin/bash
>#SBATCH --job-name=Spades_pacbio
>#SBATCH --cpus-per-task=64
>#SBATCH --time=2:00:00
>#SBATCH --mem=140g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=<YOUR ID>@unr.edu
>#SBATCH -o Spades.out # STDOUT
>#SBATCH -e Spades.err # STDERR
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>zcat  <LOCATION_BCH709_Pacbio_02.fastq.gz> <LOCATION_BCH709_Pacbio_01.fastq.gz> >> merged_pacbio.fastq
>spades.py -k 21,33,55,77 --careful -1 <Illumina_trim_galore output_reads1> -2 <Illumina_trim_galore output_reads2> --pacbio merged_pacbio.fastq -o spades_output --memory 140 --threads 64
>```
>
>
>
>## Canu assembly
>```
>mkdir canu_pacbio
>cd canu_pacbio
>conda activate genomeassembly
>```
>### Submit below job
>```bash
>#!/bin/bash
>#SBATCH --job-name=Canu_pacbio
>#SBATCH --cpus-per-task=64
>#SBATCH --time=2:00:00
>#SBATCH --mem=140g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=<YOUR ID>@unr.edu
>#SBATCH -o Canu.out # STDOUT
>#SBATCH -e Canu.err # STDERR
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>
>canu -p bch709 -d canu_outdir genomeSize=11m -pacbio-raw <LOCATION_BCH709_Pacbio_02.fastq.gz> <LOCATION_BCH709_Pacbio_01.fastq.gz>   corMemory=186 corThreads=64 batMemory=186  ovbMemory=24 ovbThreads=12 corOutCoverage=120  ovsMemory=32-186 maxMemory=186 ovsThreads=20 gridOptions='--time=12-00:00:00 -p cpu-s2-core-0 -A cpu-s2-bch709-0'
>```
{: .callout}



>## Homework 11/05/2019 due by 11/07
>Please calculate N50 of `before_rr.fasta` and check how many reads in "BCH709_0001.fastq.gz" send to `wyim@unr.edu`
>
>### N50 calculation
>```bash
>conda activate genomeassembly
>conda install -c bioconda assembly-stats
>assembly-stats before_rr.fasta
>```
>
>### Reads download
>```bash
>mkdir PacBio
>cd PacBio
>```
>
>```bash
>https://www.dropbox.com/s/1e4a4jpp3eszt6u/BCH709_Pacbio_02.fastq.gz
>https://www.dropbox.com/s/au2528hpm8vvr4c/BCH709_Pacbio_01.fastq.gz
>```
>Please count how many reads are there. (`zcat`, `wc`, `bc`, etcs)
{: .solution}


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
> echo $(zcat WGS_R1.fq.gz |wc -l)/4 | bc
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
{: .solution}




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

## Why assemble genomes
- Produce a reference for new species
- Genomic variation can be studied by comparing against the reference
- A wide variety of molecular biological tools become available or more effective
- Coding sequences can be studied in the context of their related non-coding (eg regulatory) sequences
- High level genome structure (number, arrangement of genes and repeats) can be studied

## Assembly is very challenging (“impossible”) because
- Sequencing bias under represents certain regions
- Reads are short relative to genome size
- Repeats create tangled hubs in the assembly graph
- Sequencing errors cause detours and bubbles in the assembly graph

## Prerequisite
### Flow Cytometry
![flowcytometry]({{site.baseurl}}/fig/flowcytometry.png)


### K-mer spectrum
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
gunzip ../trimmed_fastq/*.gz
jellyfish count -C -m 21 -s 1000000000 -t 10 <trim_galore output>  -o reads.jf
jellyfish histo -t 10 reads.jf > reads_jf.histo
```

### for multiqc
```bash
jellyfish:
  fn: '*_jf.hist'
```

### Count Reads Number in file

```bash
echo $(zcat WGS_R1.fq.gz |wc -l)/4 | bc

echo $(cat WGS_R1.fq |wc -l)/4 | bc
```

### Advanced approach

```bash
for i in `ls *.fq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done
***For all gzip compressed fastq files, display the number of reads since 4 lines = 1 reads***
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



###Spade
![spades]({{site.baseurl}}/fig/spades.jpg)
![spades2]({{site.baseurl}}/fig/spades2.jpg)
## Log
```bash
Command line: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades.py -k21,33,55,77     --careful       -1      /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz    -2      /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz    -o      /data/gpfs/assoc/bch709/spiderman/gee/spades_output       --memory        140     --threads       64

System information:
  SPAdes version: 3.13.1
  Python version: 3.7.3
  OS: Linux-3.10.0-957.27.2.el7.x86_64-x86_64-with-centos-7.6.1810-Core

Output dir: /data/gpfs/assoc/bch709/spiderman/gee/spades_output
Mode: read error correction and assembling
Debug mode is turned OFF

Dataset parameters:
  Multi-cell mode (you should set '--sc' flag if input data was obtained with MDA (single-cell) technology or --meta flag if processing metagenomic dataset)
  Reads:
    Library number: 1, library type: paired-end
      orientation: fr
      left reads: ['/data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz']
      right reads: ['/data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz']
      interlaced reads: not specified
      single reads: not specified
      merged reads: not specified
Read error correction parameters:
  Iterations: 1
  PHRED offset will be auto-detected
  Corrected reads will be compressed
Assembly parameters:
  k: [21, 33, 55, 77]
  Repeat resolution is enabled
  Mismatch careful mode is turned ON
  MismatchCorrector will be used
  Coverage cutoff is turned OFF
Other parameters:
  Dir for temp files: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/tmp
  Threads: 64
  Memory limit (in Gb): 140


======= SPAdes pipeline started. Log can be found here: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/spades.log


===== Read error correction started.


== Running read error correction tool: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-hammer /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/configs/config.info

  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  75)   Starting BayesHammer, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  76)   Loading config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/configs/config.info
  0:00:00.001     4M / 4M    INFO    General                 (main.cpp                  :  78)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.001     4M / 4M    INFO    General                 (memory_limit.cpp          :  49)   Memory limit set to 140 Gb
  0:00:00.001     4M / 4M    INFO    General                 (main.cpp                  :  86)   Trying to determine PHRED offset
  0:00:00.038     4M / 4M    INFO    General                 (main.cpp                  :  92)   Determined value is 33
  0:00:00.038     4M / 4M    INFO    General                 (hammer_tools.cpp          :  36)   Hamming graph threshold tau=1, k=21, subkmer positions = [ 0 10 ]
  0:00:00.038     4M / 4M    INFO    General                 (main.cpp                  : 113)   Size of aux. kmer data 24 bytes
     === ITERATION 0 begins ===
  0:00:00.042     4M / 4M    INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:00.042     4M / 4M    INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:00.043     4M / 4M    INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:00.043     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45829 Gb
  0:00:00.043     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 131072
  0:00:02.300    17G / 17G   INFO   K-mer Splitting          (kmer_data.cpp             :  97)   Processing /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz
  0:00:19.373    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             : 107)   Processed 3022711 reads
  0:00:19.373    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             :  97)   Processing /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
  0:00:37.483    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             : 107)   Processed 6045422 reads
  0:00:37.483    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             : 112)   Total 6045422 reads processed
  0:00:39.173   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:43.628   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 318799406 kmers in total.
  0:00:43.628   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:48.548   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:58.437   320M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:01:00.516   320M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 147835808 bytes occupied (3.70981 bits per kmer).
  0:01:00.518   320M / 18G   INFO   K-mer Counting           (kmer_data.cpp             : 356)   Arranging kmers in hash map order
  0:02:50.247     5G / 18G   INFO    General                 (main.cpp                  : 148)   Clustering Hamming graph.
  0:05:15.609     5G / 18G   INFO    General                 (main.cpp                  : 155)   Extracting clusters
  0:06:20.894     5G / 18G   INFO    General                 (main.cpp                  : 167)   Clustering done. Total clusters: 47999941
  0:06:20.900     2G / 18G   INFO   K-mer Counting           (kmer_data.cpp             : 376)   Collecting K-mer information, this takes a while.
  0:06:23.367     9G / 18G   INFO   K-mer Counting           (kmer_data.cpp             : 382)   Processing /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz
  0:06:41.328     9G / 18G   INFO   K-mer Counting           (kmer_data.cpp             : 382)   Processing /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
  0:06:59.453     9G / 18G   INFO   K-mer Counting           (kmer_data.cpp             : 389)   Collection done, postprocessing.
  0:07:00.669     9G / 18G   INFO   K-mer Counting           (kmer_data.cpp             : 403)   There are 318799406 kmers in total. Among them 268145204 (84.1109%) are singletons.
  0:07:00.669     9G / 18G   INFO    General                 (main.cpp                  : 173)   Subclustering Hamming graph
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 649)   Subclustering done. Total 11739 non-read kmers were generated.
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 650)   Subclustering statistics:
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 651)     Total singleton hamming clusters: 31640728. Among them 6970 (0.0220286%) are good
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 652)     Total singleton subclusters: 14379. Among them 2710 (18.8469%) are good
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 653)     Total non-singleton subcluster centers: 20505272. Among them 19729791 (96.2181%) are good
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 654)     Average size of non-trivial subcluster: 14.0052 kmers
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 655)     Average number of sub-clusters per non-singleton cluster: 1.25432  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 656)     Total solid k-mers: 19739471
  0:15:33.177     9G / 18G   INFO   Hamming Subclustering    (kmer_cluster.cpp          : 657)     Substitution probabilities: [4,4]((0.944751,0.0179573,0.0182179,0.0190737),(0.0183827,0.9447,0.0178155,0.0191018),(0.0189441,0.0177016,0.94483,0.0185242),(0.0190054,0.0181679,0.0179184,0.944908))
  0:15:33.233     9G / 18G   INFO    General                 (main.cpp                  : 178)   Finished clustering.
  0:15:33.233     9G / 18G   INFO    General                 (main.cpp                  : 197)   Starting solid k-mers expansion in 32 threads.
  0:16:01.042     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 0 produced 1807215 new k-mers.
  0:16:28.728     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 1 produced 152036 new k-mers.
  0:16:56.455     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 2 produced 13249 new k-mers.
  0:17:24.110     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 3 produced 1204 new k-mers.
  0:17:51.837     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 4 produced 93 new k-mers.
  0:18:19.634     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 5 produced 20 new k-mers.
  0:18:47.540     9G / 18G   INFO    General                 (main.cpp                  : 218)   Solid k-mers iteration 6 produced 7 new k-mers.
  0:18:47.540     9G / 18G   INFO    General                 (main.cpp                  : 222)   Solid k-mers finalized
  0:18:47.540     9G / 18G   INFO    General                 (hammer_tools.cpp          : 220)   Starting read correction in 32 threads.
  0:18:47.540     9G / 18G   INFO    General                 (hammer_tools.cpp          : 233)   Correcting pair of reads: /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz and /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
  0:19:07.696    12G / 18G   INFO    General                 (hammer_tools.cpp          : 168)   Prepared batch 0 of 3022711 reads.
  0:19:28.112    12G / 18G   INFO    General                 (hammer_tools.cpp          : 175)   Processed batch 0
  0:19:33.861    12G / 18G   INFO    General                 (hammer_tools.cpp          : 185)   Written batch 0
  0:19:35.020     9G / 18G   INFO    General                 (hammer_tools.cpp          : 274)   Correction done. Changed 10322067 bases in 4896755 reads.
  0:19:35.020     9G / 18G   INFO    General                 (hammer_tools.cpp          : 275)   Failed to correct 299 bases out of 782307432.
  0:19:35.053   128M / 18G   INFO    General                 (main.cpp                  : 255)   Saving corrected dataset description to /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/corrected.yaml
  0:19:35.054   128M / 18G   INFO    General                 (main.cpp                  : 262)   All done. Exiting.

== Compressing corrected reads (with pigz)

== Dataset description file was created: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/corrected.yaml


===== Read error correction finished.


===== Assembling started.


== Running assembler: K21

  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K21/configs/config.info
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K21/configs/careful_mode.info
  0:00:00.000     4M / 4M    INFO    General                 (memory_limit.cpp          :  49)   Memory limit set to 140 Gb
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  87)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  88)   Maximum k-mer length: 128
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  89)   Assembling dataset (/data/gpfs/assoc/bch709/spiderman/gee/spades_output/dataset.info) with K=21
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  90)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  52)   SPAdes started
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  59)   Starting from stage: construction
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  66)   Two-step RR enabled: 0
  0:00:00.000     4M / 4M    INFO   StageManager             (stage.cpp                 : 132)   STAGE == de Bruijn graph construction
  0:00:00.001     4M / 4M    INFO    General                 (read_converter.hpp        :  77)   Converting reads to binary format for library #0 (takes a while)
  0:00:00.001     4M / 4M    INFO    General                 (read_converter.hpp        :  78)   Converting paired reads
  0:00:00.224    80M / 80M   INFO    General                 (binary_converter.hpp      :  93)   16384 reads processed
  0:00:00.400    92M / 92M   INFO    General                 (binary_converter.hpp      :  93)   32768 reads processed
  0:00:00.752   112M / 112M  INFO    General                 (binary_converter.hpp      :  93)   65536 reads processed
  0:00:01.451   160M / 160M  INFO    General                 (binary_converter.hpp      :  93)   131072 reads processed
  0:00:02.845   252M / 252M  INFO    General                 (binary_converter.hpp      :  93)   262144 reads processed
  0:00:06.716   360M / 360M  INFO    General                 (binary_converter.hpp      :  93)   524288 reads processed
  0:00:14.404   424M / 424M  INFO    General                 (binary_converter.hpp      :  93)   1048576 reads processed
  0:00:29.053   436M / 436M  INFO    General                 (binary_converter.hpp      :  93)   2097152 reads processed
  0:00:43.429   336M / 436M  INFO    General                 (binary_converter.hpp      : 117)   3022500 reads written
  0:00:43.649     4M / 436M  INFO    General                 (read_converter.hpp        :  87)   Converting single reads
  0:00:43.742   132M / 436M  INFO    General                 (binary_converter.hpp      : 117)   211 reads written
  0:00:43.743     4M / 436M  INFO    General                 (read_converter.hpp        :  95)   Converting merged reads
  0:00:43.831   132M / 436M  INFO    General                 (binary_converter.hpp      : 117)   0 reads written
  0:00:43.836     4M / 436M  INFO    General                 (construction.cpp          : 111)   Max read length 130
  0:00:43.836     4M / 436M  INFO    General                 (construction.cpp          : 117)   Average read length 129.408
  0:00:43.836     4M / 436M  INFO    General                 (stage.cpp                 : 101)   PROCEDURE == k+1-mer counting
  0:00:43.837     4M / 436M  INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 1024 files using 32 threads. This might take a while.
  0:00:43.838     4M / 436M  INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:43.838     4M / 436M  INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45829 Gb
  0:00:43.838     4M / 436M  INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 65536
  0:00:52.146    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 289)   Processed 12090422 reads
  0:00:52.146    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 295)   Adding contigs from previous K
  0:00:53.617   132M / 18G   INFO    General                 (kmer_splitters.hpp        : 308)   Used 12090422 reads
  0:00:53.617   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:54.511   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 13905356 kmers in total.
  0:00:54.511   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:55.190   128M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Extension index construction
  0:00:55.190   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:55.190   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:55.190   128M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:55.190   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45703 Gb
  0:00:55.190   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 131072
  0:00:58.373    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 380)   Processed 13905356 kmers
  0:00:58.373    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 385)   Used 13905356 kmers.
  0:01:00.118   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:01:00.573   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 13723215 kmers in total.
  0:01:00.573   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:01:00.954   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:01:01.253   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:01:01.347   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 6371792 bytes occupied (3.71446 bits per kmer).
  0:01:01.356   144M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build:  99)   Building k-mer extensions from k+1-mers
  0:01:01.625   144M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 103)   Building k-mer extensions from k+1-mers finished.
  0:01:01.627   144M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Early tip clipping
  0:01:01.627   144M / 18G   INFO    General                 (construction.cpp          : 253)   Early tip clipper length bound set as (RL - K)
  0:01:01.627   144M / 18G   INFO   Early tip clipping       (early_simplification.hpp  : 181)   Early tip clipping
  0:01:06.382   144M / 18G   INFO   Early tip clipping       (early_simplification.hpp  : 184)   1510915 22-mers were removed by early tip clipper
  0:01:06.382   144M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Condensing graph
  0:01:06.393   144M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 355)   Extracting unbranching paths
  0:01:06.791   160M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 374)   Extracting unbranching paths finished. 471430 sequences extracted
  0:01:07.038   160M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 310)   Collecting perfect loops
  0:01:07.152   160M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 343)   Collecting perfect loops finished. 0 loops collected
  0:01:07.280   232M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Filling coverage indices (PHM)
  0:01:07.280   232M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:01:07.280   232M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:01:07.468   280M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 6464608 bytes occupied (3.7192 bits per kmer).
  0:01:07.482   336M / 18G   INFO    General                 (construction.cpp          : 388)   Collecting k-mer coverage information from reads, this takes a while.
  0:01:14.159   336M / 18G   INFO    General                 (construction.cpp          : 508)   Filling coverage and flanking coverage from PHM
  0:01:14.413   336M / 18G   INFO    General                 (construction.cpp          : 464)   Processed 942734 edges
  0:01:14.456   264M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == EC Threshold Finding
  0:01:14.459   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 181)   Kmer coverage valley at: 18
  0:01:14.460   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 201)   K-mer histogram maximum: 56
  0:01:14.460   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 237)   Estimated median coverage: 57. Coverage mad: 8.8956
  0:01:14.460   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 259)   Fitting coverage model
  0:01:14.547   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 2
  0:01:14.774   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 4
  0:01:15.594   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 8
  0:01:17.226   264M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 16
  0:01:19.130   200M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 309)   Fitted mean coverage: 57.5973. Fitted coverage std. dev: 7.63442
  0:01:19.133   200M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 334)   Probability of erroneous kmer at valley: 0.999423
  0:01:19.133   200M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 358)   Preliminary threshold calculated as: 36
  0:01:19.133   200M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 362)   Threshold adjusted to: 36
  0:01:19.133   200M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 375)   Estimated genome size (ignoring repeats): 10414853
  0:01:19.133   196M / 18G   INFO    General                 (genomic_info_filler.cpp   : 112)   Mean coverage was calculated as 57.5973
  0:01:19.133   196M / 18G   INFO    General                 (genomic_info_filler.cpp   : 127)   EC coverage threshold value was calculated as 36
  0:01:19.133   196M / 18G   INFO    General                 (genomic_info_filler.cpp   : 128)   Trusted kmer low bound: 0
  0:01:19.133   196M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Raw Simplification
  0:01:19.133   196M / 18G   INFO    General                 (simplification.cpp        : 128)   PROCEDURE == InitialCleaning
  0:01:19.133   196M / 18G   INFO    General                 (graph_simplification.hpp  : 669)   Flanking coverage based disconnection disabled
  0:01:19.133   196M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Self conjugate edge remover
  0:01:19.155   196M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Self conjugate edge remover triggered 0 times
  0:01:19.155   196M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification
  0:01:19.155   196M / 18G   INFO    General                 (simplification.cpp        : 357)   Graph simplification started
  0:01:19.155   196M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:01:19.155   196M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 1
  0:01:19.155   196M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:19.186   180M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 258 times
  0:01:19.186   180M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:22.193   252M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 31922 times
  0:01:22.193   252M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:24.254   276M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 100239 times
  0:01:24.254   276M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 2
  0:01:24.254   276M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:24.283   276M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 94 times
  0:01:24.283   276M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:25.837   216M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 9261 times
  0:01:25.837   216M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:25.938   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 6090 times
  0:01:25.938   212M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 3
  0:01:25.938   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:25.943   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:25.943   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.269   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1276 times
  0:01:26.269   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.288   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 1393 times
  0:01:26.288   212M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 4
  0:01:26.288   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.290   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.290   212M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.367   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 252 times
  0:01:26.367   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.373   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 473 times
  0:01:26.373   208M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 5
  0:01:26.373   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.373   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.373   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.410   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 103 times
  0:01:26.410   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.413   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 225 times
  0:01:26.413   208M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 6
  0:01:26.413   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.413   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.413   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.436   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 71 times
  0:01:26.436   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.438   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 137 times
  0:01:26.438   208M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 7
  0:01:26.438   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.438   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.438   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.456   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 48 times
  0:01:26.456   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.457   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 117 times
  0:01:26.457   208M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 8
  0:01:26.457   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.457   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.457   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.470   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 35 times
  0:01:26.470   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.471   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 58 times
  0:01:26.471   208M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 9
  0:01:26.471   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.471   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.471   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.475   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 14 times
  0:01:26.475   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.476   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 48 times
  0:01:26.476   208M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 10
  0:01:26.476   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.476   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.476   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.481   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 16 times
  0:01:26.481   208M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.482   204M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 76 times
  0:01:26.482   204M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 11
  0:01:26.482   204M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.486   204M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 2 times
  0:01:26.486   204M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.686   252M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 691 times
  0:01:26.687   252M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.696   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:01:26.696   248M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 12
  0:01:26.696   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.699   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.699   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.699   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:26.699   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:26.699   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:01:26.699   248M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification Cleanup
  0:01:26.699   248M / 18G   INFO    General                 (simplification.cpp        : 196)   PROCEDURE == Post simplification
  0:01:26.699   248M / 18G   INFO    General                 (graph_simplification.hpp  : 458)   Disconnection of relatively low covered edges disabled
  0:01:26.699   248M / 18G   INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:01:26.699   248M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:01:26.699   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.703   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.703   248M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:26.863   244M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 87 times
  0:01:26.863   244M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:26.868   244M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:26.868   244M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:27.019   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 10 times
  0:01:27.019   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:27.019   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:27.019   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:27.019   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:27.019   240M / 18G   INFO    General                 (simplification.cpp        : 330)   Disrupting self-conjugate edges
  0:01:27.072   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Removing isolated edges
  0:01:27.076   240M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Removing isolated edges triggered 0 times
  0:01:27.076   240M / 18G   INFO    General                 (simplification.cpp        : 470)   Counting average coverage
  0:01:27.092   240M / 18G   INFO    General                 (simplification.cpp        : 476)   Average coverage = 64.7081
  0:01:27.092   240M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Contig Output
  0:01:27.092   240M / 18G   INFO    General                 (contig_output.hpp         :  22)   Outputting contigs to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K21/simplified_contigs.fasta
  0:01:27.662   240M / 18G   INFO    General                 (launch.hpp                : 151)   SPAdes finished
  0:01:27.853   128M / 18G   INFO    General                 (main.cpp                  : 109)   Assembling time: 0 hours 1 minutes 27 seconds
Max read length detected as 130

== Running assembler: K33

  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K33/configs/config.info
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K33/configs/careful_mode.info
  0:00:00.000     4M / 4M    INFO    General                 (memory_limit.cpp          :  49)   Memory limit set to 140 Gb
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  87)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  88)   Maximum k-mer length: 128
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  89)   Assembling dataset (/data/gpfs/assoc/bch709/spiderman/gee/spades_output/dataset.info) with K=33
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  90)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  52)   SPAdes started
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  59)   Starting from stage: construction
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  66)   Two-step RR enabled: 0
  0:00:00.000     4M / 4M    INFO   StageManager             (stage.cpp                 : 132)   STAGE == de Bruijn graph construction
  0:00:00.000     4M / 4M    INFO    General                 (read_converter.hpp        :  59)   Binary reads detected
  0:00:00.002     4M / 4M    INFO    General                 (construction.cpp          : 111)   Max read length 130
  0:00:00.002     4M / 4M    INFO    General                 (construction.cpp          : 117)   Average read length 129.408
  0:00:00.002     4M / 4M    INFO    General                 (stage.cpp                 : 101)   PROCEDURE == k+1-mer counting
  0:00:00.003     4M / 4M    INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 1024 files using 32 threads. This might take a while.
  0:00:00.004     4M / 4M    INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:00.004     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45829 Gb
  0:00:00.004     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 32768
  0:00:08.184    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 289)   Processed 12090422 reads
  0:00:08.184    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 295)   Adding contigs from previous K
  0:00:12.153   132M / 18G   INFO    General                 (kmer_splitters.hpp        : 308)   Used 12090422 reads
  0:00:12.153   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:13.121   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 15869046 kmers in total.
  0:00:13.121   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:13.954   128M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Extension index construction
  0:00:13.955   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:13.955   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:13.955   128M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:13.955   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45703 Gb
  0:00:13.955   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 65536
  0:00:17.375    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 380)   Processed 15869046 kmers
  0:00:17.375    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 385)   Used 15869046 kmers.
  0:00:19.295   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:19.837   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 15720358 kmers in total.
  0:00:19.837   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:20.485   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:20.950   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:00:21.173   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 7297672 bytes occupied (3.71374 bits per kmer).
  0:00:21.184   144M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build:  99)   Building k-mer extensions from k+1-mers
  0:00:21.553   144M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 103)   Building k-mer extensions from k+1-mers finished.
  0:00:21.555   144M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Early tip clipping
  0:00:21.555   144M / 18G   INFO    General                 (construction.cpp          : 253)   Early tip clipper length bound set as (RL - K)
  0:00:21.555   144M / 18G   INFO   Early tip clipping       (early_simplification.hpp  : 181)   Early tip clipping
  0:00:31.948   144M / 18G   INFO   Early tip clipping       (early_simplification.hpp  : 184)   2455865 34-mers were removed by early tip clipper
  0:00:31.948   144M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Condensing graph
  0:00:31.958   144M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 355)   Extracting unbranching paths
  0:00:32.405   160M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 374)   Extracting unbranching paths finished. 398618 sequences extracted
  0:00:32.693   160M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 310)   Collecting perfect loops
  0:00:32.821   160M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 343)   Collecting perfect loops finished. 0 loops collected
  0:00:32.930   216M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Filling coverage indices (PHM)
  0:00:32.930   216M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:32.930   216M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:33.140   280M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 7374928 bytes occupied (3.71789 bits per kmer).
  0:00:33.157   344M / 18G   INFO    General                 (construction.cpp          : 388)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:39.359   344M / 18G   INFO    General                 (construction.cpp          : 508)   Filling coverage and flanking coverage from PHM
  0:00:39.646   344M / 18G   INFO    General                 (construction.cpp          : 464)   Processed 797195 edges
  0:00:39.704   256M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == EC Threshold Finding
  0:00:39.705   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 181)   Kmer coverage valley at: 16
  0:00:39.705   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 201)   K-mer histogram maximum: 49
  0:00:39.705   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 237)   Estimated median coverage: 51. Coverage mad: 8.8956
  0:00:39.705   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 259)   Fitting coverage model
  0:00:39.774   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 2
  0:00:39.950   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 4
  0:00:40.499   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 8
  0:00:41.109   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 309)   Fitted mean coverage: 51.0559. Fitted coverage std. dev: 7.2284
  0:00:41.111   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 334)   Probability of erroneous kmer at valley: 0.997188
  0:00:41.111   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 358)   Preliminary threshold calculated as: 30
  0:00:41.111   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 362)   Threshold adjusted to: 30
  0:00:41.111   256M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 375)   Estimated genome size (ignoring repeats): 10731880
  0:00:41.111   256M / 18G   INFO    General                 (genomic_info_filler.cpp   : 112)   Mean coverage was calculated as 51.0559
  0:00:41.111   256M / 18G   INFO    General                 (genomic_info_filler.cpp   : 127)   EC coverage threshold value was calculated as 30
  0:00:41.111   256M / 18G   INFO    General                 (genomic_info_filler.cpp   : 128)   Trusted kmer low bound: 0
  0:00:41.111   256M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Raw Simplification
  0:00:41.111   256M / 18G   INFO    General                 (simplification.cpp        : 128)   PROCEDURE == InitialCleaning
  0:00:41.111   256M / 18G   INFO    General                 (graph_simplification.hpp  : 669)   Flanking coverage based disconnection disabled
  0:00:41.111   256M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Self conjugate edge remover
  0:00:41.128   256M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Self conjugate edge remover triggered 0 times
  0:00:41.129   256M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification
  0:00:41.129   256M / 18G   INFO    General                 (simplification.cpp        : 357)   Graph simplification started
  0:00:41.129   256M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:00:41.129   256M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 1
  0:00:41.129   256M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:41.160   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 710 times
  0:00:41.160   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:45.175   156M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 29107 times
  0:00:45.175   156M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:47.399   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 92363 times
  0:00:47.399   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 2
  0:00:47.399   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:47.425   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 531 times
  0:00:47.425   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.075   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 6878 times
  0:00:49.075   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.161   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 4277 times
  0:00:49.161   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 3
  0:00:49.161   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.164   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 3 times
  0:00:49.164   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.519   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1380 times
  0:00:49.519   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.530   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 664 times
  0:00:49.530   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 4
  0:00:49.530   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.531   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.531   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.603   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 273 times
  0:00:49.603   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.606   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 186 times
  0:00:49.606   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 5
  0:00:49.606   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.607   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.607   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.626   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 69 times
  0:00:49.626   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.627   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 70 times
  0:00:49.627   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 6
  0:00:49.627   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.627   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.627   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.635   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 31 times
  0:00:49.635   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.635   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 35 times
  0:00:49.635   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 7
  0:00:49.635   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.635   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.635   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.640   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 16 times
  0:00:49.640   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.641   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 11 times
  0:00:49.641   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 8
  0:00:49.641   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.641   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.641   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.641   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 3 times
  0:00:49.641   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 9 times
  0:00:49.642   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 9
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1 times
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 12 times
  0:00:49.642   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 10
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.642   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.643   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 3 times
  0:00:49.643   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.643   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 24 times
  0:00:49.643   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 11
  0:00:49.643   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.648   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 2 times
  0:00:49.648   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.729   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 406 times
  0:00:49.729   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.733   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:00:49.734   152M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 12
  0:00:49.734   152M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.734   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.734   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.734   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:00:49.734   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:00:49.734   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:00:49.734   148M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification Cleanup
  0:00:49.734   148M / 18G   INFO    General                 (simplification.cpp        : 196)   PROCEDURE == Post simplification
  0:00:49.734   148M / 18G   INFO    General                 (graph_simplification.hpp  : 458)   Disconnection of relatively low covered edges disabled
  0:00:49.734   148M / 18G   INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:00:49.734   148M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:00:49.734   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.737   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.737   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.789   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 35 times
  0:00:49.789   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.791   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.791   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.841   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 7 times
  0:00:49.841   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:49.841   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:00:49.841   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:00:49.841   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:00:49.841   148M / 18G   INFO    General                 (simplification.cpp        : 330)   Disrupting self-conjugate edges
  0:00:49.863   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Removing isolated edges
  0:00:49.865   148M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Removing isolated edges triggered 50 times
  0:00:49.865   148M / 18G   INFO    General                 (simplification.cpp        : 470)   Counting average coverage
  0:00:49.872   148M / 18G   INFO    General                 (simplification.cpp        : 476)   Average coverage = 55.1933
  0:00:49.873   148M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Contig Output
  0:00:49.873   148M / 18G   INFO    General                 (contig_output.hpp         :  22)   Outputting contigs to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K33/simplified_contigs.fasta
  0:00:50.194   148M / 18G   INFO    General                 (launch.hpp                : 151)   SPAdes finished
  0:00:50.309   128M / 18G   INFO    General                 (main.cpp                  : 109)   Assembling time: 0 hours 0 minutes 50 seconds

== Running assembler: K55

  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K55/configs/config.info
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K55/configs/careful_mode.info
  0:00:00.000     4M / 4M    INFO    General                 (memory_limit.cpp          :  49)   Memory limit set to 140 Gb
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  87)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  88)   Maximum k-mer length: 128
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  89)   Assembling dataset (/data/gpfs/assoc/bch709/spiderman/gee/spades_output/dataset.info) with K=55
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  90)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  52)   SPAdes started
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  59)   Starting from stage: construction
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  66)   Two-step RR enabled: 0
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  77)   Will need read mapping, kmer mapper will be attached
  0:00:00.000     4M / 4M    INFO   StageManager             (stage.cpp                 : 132)   STAGE == de Bruijn graph construction
  0:00:00.001     4M / 4M    INFO    General                 (read_converter.hpp        :  59)   Binary reads detected
  0:00:00.002     4M / 4M    INFO    General                 (construction.cpp          : 111)   Max read length 130
  0:00:00.002     4M / 4M    INFO    General                 (construction.cpp          : 117)   Average read length 129.408
  0:00:00.002     4M / 4M    INFO    General                 (stage.cpp                 : 101)   PROCEDURE == k+1-mer counting
  0:00:00.002     4M / 4M    INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 1024 files using 32 threads. This might take a while.
  0:00:00.003     4M / 4M    INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:00.003     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45829 Gb
  0:00:00.003     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 32768
  0:00:07.186    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 289)   Processed 12090422 reads
  0:00:07.186    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 295)   Adding contigs from previous K
  0:00:11.082   132M / 18G   INFO    General                 (kmer_splitters.hpp        : 308)   Used 12090422 reads
  0:00:11.082   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:12.153   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 18215028 kmers in total.
  0:00:12.153   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:13.077   128M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Extension index construction
  0:00:13.077   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:13.077   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:13.077   128M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:13.077   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45703 Gb
  0:00:13.077   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 65536
  0:00:16.511    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 380)   Processed 18215028 kmers
  0:00:16.511    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 385)   Used 18215028 kmers.
  0:00:18.436   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:18.970   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 18150930 kmers in total.
  0:00:18.970   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:19.684   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:20.206   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:00:20.461   132M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 8424816 bytes occupied (3.71323 bits per kmer).
  0:00:20.472   152M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build:  99)   Building k-mer extensions from k+1-mers
  0:00:20.855   152M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 103)   Building k-mer extensions from k+1-mers finished.
  0:00:20.857   148M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Condensing graph
  0:00:20.876   148M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 355)   Extracting unbranching paths
  0:00:32.290   180M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 374)   Extracting unbranching paths finished. 828017 sequences extracted
  0:00:32.695   176M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 310)   Collecting perfect loops
  0:00:32.842   180M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 343)   Collecting perfect loops finished. 0 loops collected
  0:00:33.120   360M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Filling coverage indices (PHM)
  0:00:33.120   360M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:33.120   360M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:33.354   396M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 8463032 bytes occupied (3.71694 bits per kmer).
  0:00:33.374   468M / 18G   INFO    General                 (construction.cpp          : 388)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:38.330   464M / 18G   INFO    General                 (construction.cpp          : 508)   Filling coverage and flanking coverage from PHM
  0:00:38.743   464M / 18G   INFO    General                 (construction.cpp          : 464)   Processed 1656029 edges
  0:00:38.804   336M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == EC Threshold Finding
  0:00:38.804   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 181)   Kmer coverage valley at: 8
  0:00:38.804   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 201)   K-mer histogram maximum: 38
  0:00:38.804   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 237)   Estimated median coverage: 39. Coverage mad: 7.413
  0:00:38.804   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 259)   Fitting coverage model
  0:00:38.860   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 2
  0:00:38.996   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 4
  0:00:39.504   336M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 8
  0:00:40.461   324M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 16
  0:00:40.461   324M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 309)   Fitted mean coverage: 39.1293. Fitted coverage std. dev: 6.42551
  0:00:40.463   324M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 334)   Probability of erroneous kmer at valley: 0.996761
  0:00:40.463   324M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 358)   Preliminary threshold calculated as: 18
  0:00:40.463   324M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 362)   Threshold adjusted to: 18
  0:00:40.463   324M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 375)   Estimated genome size (ignoring repeats): 10980914
  0:00:40.463   324M / 18G   INFO    General                 (genomic_info_filler.cpp   : 112)   Mean coverage was calculated as 39.1293
  0:00:40.463   324M / 18G   INFO    General                 (genomic_info_filler.cpp   : 127)   EC coverage threshold value was calculated as 18
  0:00:40.463   324M / 18G   INFO    General                 (genomic_info_filler.cpp   : 128)   Trusted kmer low bound: 0
  0:00:40.463   324M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Gap Closer
  0:00:40.463   324M / 18G   INFO    General                 (graph_pack.hpp            : 101)   Index refill
  0:00:40.464   324M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:40.464   324M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:40.464   324M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:40.464   324M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45504 Gb
  0:00:40.464   324M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 65536
  0:00:43.671    17G / 18G   INFO    General                 (edge_index_builders.hpp   :  77)   Processed 1656029 edges
  0:00:43.671    17G / 18G   INFO    General                 (edge_index_builders.hpp   :  82)   Used 1656029 sequences.
  0:00:45.372   292M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:45.961   292M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 18215028 kmers in total.
  0:00:45.961   292M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:46.685   292M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:47.235   340M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:00:47.490   340M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 8454640 bytes occupied (3.71326 bits per kmer).
  0:00:47.756   760M / 18G   INFO    General                 (edge_index_builders.hpp   : 107)   Collecting edge information from graph, this takes a while.
  0:00:48.048   760M / 18G   INFO    General                 (edge_index.hpp            :  92)   Index refilled
  0:00:48.049   760M / 18G   INFO    General                 (gap_closer.cpp            : 159)   Preparing shift maps
  0:00:48.751   812M / 18G   INFO    General                 (gap_closer.cpp            : 119)   Processing paired reads (takes a while)
  0:00:50.278   848M / 18G   INFO    General                 (gap_closer.cpp            : 138)   Used 3022500 paired reads
  0:00:50.278   848M / 18G   INFO    General                 (gap_closer.cpp            : 140)   Merging paired indices
  0:00:50.388   784M / 18G   INFO   GapCloser                (gap_closer.cpp            : 346)   Closing short gaps
  0:00:51.816   784M / 18G   INFO   GapCloser                (gap_closer.cpp            : 380)   Closing short gaps complete: filled 0 gaps after checking 892 candidates
  0:00:51.851   784M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Raw Simplification
  0:00:51.854   356M / 18G   INFO    General                 (simplification.cpp        : 128)   PROCEDURE == InitialCleaning
  0:00:51.854   356M / 18G   INFO    General                 (graph_simplification.hpp  : 669)   Flanking coverage based disconnection disabled
  0:00:51.854   356M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Self conjugate edge remover
  0:00:51.901   356M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Self conjugate edge remover triggered 0 times
  0:00:51.902   356M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification
  0:00:51.902   356M / 18G   INFO    General                 (simplification.cpp        : 357)   Graph simplification started
  0:00:51.902   356M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:00:51.902   356M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 1
  0:00:51.902   356M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:00:55.508   360M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 325678 times
  0:00:55.508   360M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:00.215   392M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 16019 times
  0:01:00.215   392M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:01.197   384M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 35150 times
  0:01:01.197   384M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 2
  0:01:01.197   384M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:01.251   380M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 2686 times
  0:01:01.251   380M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.324   396M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 3609 times
  0:01:03.324   396M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.500   396M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 5678 times
  0:01:03.500   396M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 3
  0:01:03.500   396M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.503   396M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 83 times
  0:01:03.503   396M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.901   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1087 times
  0:01:03.901   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.914   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 415 times
  0:01:03.914   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 4
  0:01:03.914   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.915   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.915   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.959   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 155 times
  0:01:03.959   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.960   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 55 times
  0:01:03.960   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 5
  0:01:03.961   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.961   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.961   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.967   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 26 times
  0:01:03.967   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.968   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 16 times
  0:01:03.968   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 6
  0:01:03.968   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.968   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.968   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.969   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 2 times
  0:01:03.969   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.969   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 8 times
  0:01:03.969   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 7
  0:01:03.969   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.969   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.969   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.970   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1 times
  0:01:03.970   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.970   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 7 times
  0:01:03.970   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 8
  0:01:03.970   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.970   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.970   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1 times
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 10 times
  0:01:03.971   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 9
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:03.971   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 8 times
  0:01:03.972   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 10
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 10 times
  0:01:03.972   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 11
  0:01:03.972   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.974   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 3 times
  0:01:03.974   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.028   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 104 times
  0:01:04.028   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.030   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:01:04.030   400M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 12
  0:01:04.030   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.031   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.031   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.031   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:04.031   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.031   400M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:01:04.031   400M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Gap Closer
  0:01:04.031   400M / 18G   INFO    General                 (graph_pack.hpp            : 101)   Index refill
  0:01:04.032   400M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:01:04.032   400M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:01:04.033   400M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:01:04.033   400M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45426 Gb
  0:01:04.033   400M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 65536
  0:01:07.072    18G / 18G   INFO    General                 (edge_index_builders.hpp   :  77)   Processed 37240 edges
  0:01:07.072    18G / 18G   INFO    General                 (edge_index_builders.hpp   :  82)   Used 37240 sequences.
  0:01:08.858   404M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:01:09.320   404M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 11126350 kmers in total.
  0:01:09.320   404M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:01:09.827   404M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:01:10.144   404M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:01:10.296   404M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 5167360 bytes occupied (3.7154 bits per kmer).
  0:01:10.458   660M / 18G   INFO    General                 (edge_index_builders.hpp   : 107)   Collecting edge information from graph, this takes a while.
  0:01:10.663   660M / 18G   INFO    General                 (edge_index.hpp            :  92)   Index refilled
  0:01:10.664   660M / 18G   INFO    General                 (gap_closer.cpp            : 159)   Preparing shift maps
  0:01:10.679   660M / 18G   INFO    General                 (gap_closer.cpp            : 119)   Processing paired reads (takes a while)
  0:01:12.327   660M / 18G   INFO    General                 (gap_closer.cpp            : 138)   Used 3022500 paired reads
  0:01:12.327   660M / 18G   INFO    General                 (gap_closer.cpp            : 140)   Merging paired indices
  0:01:12.329   660M / 18G   INFO   GapCloser                (gap_closer.cpp            : 346)   Closing short gaps
  0:01:12.354   660M / 18G   INFO   GapCloser                (gap_closer.cpp            : 380)   Closing short gaps complete: filled 0 gaps after checking 2 candidates
  0:01:12.355   656M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification Cleanup
  0:01:12.355   656M / 18G   INFO    General                 (simplification.cpp        : 196)   PROCEDURE == Post simplification
  0:01:12.355   656M / 18G   INFO    General                 (graph_simplification.hpp  : 458)   Disconnection of relatively low covered edges disabled
  0:01:12.355   656M / 18G   INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:01:12.355   656M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:01:12.356   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:12.358   660M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:12.358   660M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:12.385   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 12 times
  0:01:12.385   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:12.386   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:12.386   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:12.407   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1 times
  0:01:12.407   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:12.407   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:12.407   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:12.407   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:12.407   656M / 18G   INFO    General                 (simplification.cpp        : 330)   Disrupting self-conjugate edges
  0:01:12.418   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Removing isolated edges
  0:01:12.572   656M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Removing isolated edges triggered 3448 times
  0:01:12.572   656M / 18G   INFO    General                 (simplification.cpp        : 470)   Counting average coverage
  0:01:12.577   656M / 18G   INFO    General                 (simplification.cpp        : 476)   Average coverage = 40.7622
  0:01:12.577   656M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Contig Output
  0:01:12.577   656M / 18G   INFO    General                 (contig_output.hpp         :  22)   Outputting contigs to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K55/simplified_contigs.fasta
  0:01:12.782   656M / 18G   INFO    General                 (launch.hpp                : 151)   SPAdes finished
  0:01:13.385   132M / 18G   INFO    General                 (main.cpp                  : 109)   Assembling time: 0 hours 1 minutes 13 seconds

== Running assembler: K77

  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K77/configs/config.info
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  74)   Loaded config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/K77/configs/careful_mode.info
  0:00:00.000     4M / 4M    INFO    General                 (memory_limit.cpp          :  49)   Memory limit set to 140 Gb
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  87)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  88)   Maximum k-mer length: 128
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  89)   Assembling dataset (/data/gpfs/assoc/bch709/spiderman/gee/spades_output/dataset.info) with K=77
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  90)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  52)   SPAdes started
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  59)   Starting from stage: construction
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  66)   Two-step RR enabled: 0
  0:00:00.000     4M / 4M    INFO    General                 (launch.hpp                :  77)   Will need read mapping, kmer mapper will be attached
  0:00:00.000     4M / 4M    INFO   StageManager             (stage.cpp                 : 132)   STAGE == de Bruijn graph construction
  0:00:00.001     4M / 4M    INFO    General                 (read_converter.hpp        :  59)   Binary reads detected
  0:00:00.002     4M / 4M    INFO    General                 (construction.cpp          : 111)   Max read length 130
  0:00:00.002     4M / 4M    INFO    General                 (construction.cpp          : 117)   Average read length 129.408
  0:00:00.002     4M / 4M    INFO    General                 (stage.cpp                 : 101)   PROCEDURE == k+1-mer counting
  0:00:00.003     4M / 4M    INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 1024 files using 32 threads. This might take a while.
  0:00:00.003     4M / 4M    INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:00.004     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45829 Gb
  0:00:00.004     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 21845
  0:00:06.511    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 289)   Processed 12090422 reads
  0:00:06.511    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 295)   Adding contigs from previous K
  0:00:10.716   132M / 18G   INFO    General                 (kmer_splitters.hpp        : 308)   Used 12090422 reads
  0:00:10.716   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:11.967   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 18597210 kmers in total.
  0:00:11.967   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:13.165   128M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Extension index construction
  0:00:13.166   128M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:13.166   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:13.166   128M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:13.166   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45703 Gb
  0:00:13.166   128M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 43690
  0:00:16.891    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 380)   Processed 18597210 kmers
  0:00:16.891    17G / 18G   INFO    General                 (kmer_splitters.hpp        : 385)   Used 18597210 kmers.
  0:00:18.620   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:19.389   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 18620971 kmers in total.
  0:00:19.389   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:20.408   132M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:21.131   132M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:00:21.522   132M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 8642840 bytes occupied (3.71316 bits per kmer).
  0:00:21.534   152M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build:  99)   Building k-mer extensions from k+1-mers
  0:00:21.961   152M / 18G   INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 103)   Building k-mer extensions from k+1-mers finished.
  0:00:21.964   148M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Condensing graph
  0:00:21.981   148M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 355)   Extracting unbranching paths
  0:00:40.222   176M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 374)   Extracting unbranching paths finished. 737063 sequences extracted
  0:00:40.672   176M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 310)   Collecting perfect loops
  0:00:40.839   176M / 18G   INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 343)   Collecting perfect loops finished. 0 loops collected
  0:00:41.128   316M / 18G   INFO    General                 (stage.cpp                 : 101)   PROCEDURE == Filling coverage indices (PHM)
  0:00:41.128   316M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:41.128   316M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:41.388   368M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 8640160 bytes occupied (3.71676 bits per kmer).
  0:00:41.408   440M / 18G   INFO    General                 (construction.cpp          : 388)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:45.192   428M / 18G   INFO    General                 (construction.cpp          : 508)   Filling coverage and flanking coverage from PHM
  0:00:45.638   432M / 18G   INFO    General                 (construction.cpp          : 464)   Processed 1474126 edges
  0:00:45.713   284M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == EC Threshold Finding
  0:00:45.713   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 181)   Kmer coverage valley at: 4
  0:00:45.713   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 201)   K-mer histogram maximum: 26
  0:00:45.713   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 237)   Estimated median coverage: 27. Coverage mad: 7.413
  0:00:45.713   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 259)   Fitting coverage model
  0:00:45.748   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 2
  0:00:45.841   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 4
  0:00:46.186   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 8
  0:00:46.884   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 295)   ... iteration 16
  0:00:47.040   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 309)   Fitted mean coverage: 27.3626. Fitted coverage std. dev: 5.40493
  0:00:47.041   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 334)   Probability of erroneous kmer at valley: 0.987151
  0:00:47.041   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 358)   Preliminary threshold calculated as: 10
  0:00:47.041   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 362)   Threshold adjusted to: 10
  0:00:47.041   284M / 18G   INFO    General                 (kmer_coverage_model.cpp   : 375)   Estimated genome size (ignoring repeats): 11080100
  0:00:47.041   284M / 18G   INFO    General                 (genomic_info_filler.cpp   : 112)   Mean coverage was calculated as 27.3626
  0:00:47.041   284M / 18G   INFO    General                 (genomic_info_filler.cpp   : 127)   EC coverage threshold value was calculated as 10
  0:00:47.041   284M / 18G   INFO    General                 (genomic_info_filler.cpp   : 128)   Trusted kmer low bound: 0
  0:00:47.041   284M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Gap Closer
  0:00:47.041   284M / 18G   INFO    General                 (graph_pack.hpp            : 101)   Index refill
  0:00:47.042   284M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:47.042   284M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:47.047   284M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:47.047   284M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45544 Gb
  0:00:47.047   284M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 43690
  0:00:50.417    17G / 18G   INFO    General                 (edge_index_builders.hpp   :  77)   Processed 1474126 edges
  0:00:50.417    17G / 18G   INFO    General                 (edge_index_builders.hpp   :  82)   Used 1474126 sequences.
  0:00:52.117   284M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:00:52.860   284M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 18597210 kmers in total.
  0:00:52.860   284M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:00:53.868   284M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:00:54.546   296M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:00:54.930   296M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 8631624 bytes occupied (3.71308 bits per kmer).
  0:00:55.199   724M / 18G   INFO    General                 (edge_index_builders.hpp   : 107)   Collecting edge information from graph, this takes a while.
  0:00:55.492   724M / 18G   INFO    General                 (edge_index.hpp            :  92)   Index refilled
  0:00:55.493   724M / 18G   INFO    General                 (gap_closer.cpp            : 159)   Preparing shift maps
  0:00:56.190   784M / 18G   INFO    General                 (gap_closer.cpp            : 119)   Processing paired reads (takes a while)
  0:00:57.425   816M / 18G   INFO    General                 (gap_closer.cpp            : 138)   Used 3022500 paired reads
  0:00:57.425   816M / 18G   INFO    General                 (gap_closer.cpp            : 140)   Merging paired indices
  0:00:57.574   752M / 18G   INFO   GapCloser                (gap_closer.cpp            : 346)   Closing short gaps
  0:00:58.942   752M / 18G   INFO   GapCloser                (gap_closer.cpp            : 380)   Closing short gaps complete: filled 0 gaps after checking 2972 candidates
  0:00:58.980   748M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Raw Simplification
  0:00:58.984   288M / 18G   INFO    General                 (simplification.cpp        : 128)   PROCEDURE == InitialCleaning
  0:00:58.984   288M / 18G   INFO    General                 (simplification.cpp        :  62)   Most init cleaning disabled on main iteration
  0:00:58.984   288M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Self conjugate edge remover
  0:00:59.030   288M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Self conjugate edge remover triggered 0 times
  0:00:59.030   288M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification
  0:00:59.030   288M / 18G   INFO    General                 (simplification.cpp        : 357)   Graph simplification started
  0:00:59.030   288M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:00:59.030   288M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 1
  0:00:59.030   288M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:03.227   300M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 336002 times
  0:01:03.227   300M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.126   320M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 2409 times
  0:01:04.126   320M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.272   320M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 4437 times
  0:01:04.272   320M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 2
  0:01:04.272   320M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.306   320M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 1817 times
  0:01:04.306   320M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.612   324M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 681 times
  0:01:04.612   324M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.698   324M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 2071 times
  0:01:04.698   324M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 3
  0:01:04.698   324M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.705   324M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 254 times
  0:01:04.705   324M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.888   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 330 times
  0:01:04.888   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.899   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 316 times
  0:01:04.899   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 4
  0:01:04.899   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.900   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 15 times
  0:01:04.900   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.922   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 65 times
  0:01:04.922   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.923   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 54 times
  0:01:04.923   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 5
  0:01:04.923   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.923   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.923   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.929   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 15 times
  0:01:04.929   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.929   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 10 times
  0:01:04.929   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 6
  0:01:04.929   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.929   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.929   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 3 times
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 4 times
  0:01:04.930   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 7
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:04.930   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.931   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 10 times
  0:01:04.931   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 8
  0:01:04.931   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.931   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.931   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.932   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1 times
  0:01:04.932   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.932   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 10 times
  0:01:04.932   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 9
  0:01:04.932   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.932   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.932   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.933   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:04.933   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.933   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 13 times
  0:01:04.933   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 10
  0:01:04.933   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.933   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.933   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.934   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 2 times
  0:01:04.934   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.935   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 20 times
  0:01:04.935   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 11
  0:01:04.935   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.942   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 2 times
  0:01:04.942   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.979   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 45 times
  0:01:04.979   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:01:04.985   332M / 18G   INFO    General                 (simplification.cpp        : 362)   PROCEDURE == Simplification cycle, iteration 12
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Low coverage edge remover
  0:01:04.985   332M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Low coverage edge remover triggered 0 times
  0:01:04.985   332M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Gap Closer
  0:01:04.985   332M / 18G   INFO    General                 (graph_pack.hpp            : 101)   Index refill
  0:01:04.986   332M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:01:04.986   332M / 18G   INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:01:04.987   332M / 18G   INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:01:04.987   332M / 18G   INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45496 Gb
  0:01:04.987   332M / 18G   INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 43690
  0:01:08.217    17G / 18G   INFO    General                 (edge_index_builders.hpp   :  77)   Processed 89330 edges
  0:01:08.217    17G / 18G   INFO    General                 (edge_index_builders.hpp   :  82)   Used 89330 sequences.
  0:01:10.092   332M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
  0:01:10.702   332M / 18G   INFO    General                 (kmer_index_builder.hpp    : 127)   K-mer counting done. There are 12987359 kmers in total.
  0:01:10.702   332M / 18G   INFO    General                 (kmer_index_builder.hpp    : 133)   Merging temporary buckets.
  0:01:11.468   332M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 314)   Building perfect hash indices
  0:01:11.952   332M / 18G   INFO    General                 (kmer_index_builder.hpp    : 150)   Merging final buckets.
  0:01:12.223   332M / 18G   INFO   K-mer Index Building     (kmer_index_builder.hpp    : 336)   Index built. Total 6030448 bytes occupied (3.71466 bits per kmer).
  0:01:12.409   632M / 18G   INFO    General                 (edge_index_builders.hpp   : 107)   Collecting edge information from graph, this takes a while.
  0:01:12.615   632M / 18G   INFO    General                 (edge_index.hpp            :  92)   Index refilled
  0:01:12.617   632M / 18G   INFO    General                 (gap_closer.cpp            : 159)   Preparing shift maps
  0:01:12.681   644M / 18G   INFO    General                 (gap_closer.cpp            : 119)   Processing paired reads (takes a while)
  0:01:13.984   644M / 18G   INFO    General                 (gap_closer.cpp            : 138)   Used 3022500 paired reads
  0:01:13.984   644M / 18G   INFO    General                 (gap_closer.cpp            : 140)   Merging paired indices
  0:01:14.006   632M / 18G   INFO   GapCloser                (gap_closer.cpp            : 346)   Closing short gaps
  0:01:14.073   632M / 18G   INFO   GapCloser                (gap_closer.cpp            : 380)   Closing short gaps complete: filled 0 gaps after checking 8 candidates
  0:01:14.081   632M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Simplification Cleanup
  0:01:14.081   632M / 18G   INFO    General                 (simplification.cpp        : 196)   PROCEDURE == Post simplification
  0:01:14.081   632M / 18G   INFO    General                 (graph_simplification.hpp  : 458)   Disconnection of relatively low covered edges disabled
  0:01:14.081   632M / 18G   INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:01:14.081   632M / 18G   INFO    General                 (graph_simplification.hpp  : 641)   Creating parallel br instance
  0:01:14.082   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:14.087   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:14.087   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:14.109   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 2 times
  0:01:14.109   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:14.114   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:14.114   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:14.135   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 1 times
  0:01:14.135   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Tip clipper
  0:01:14.135   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Tip clipper triggered 0 times
  0:01:14.135   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Bulge remover
  0:01:14.135   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Bulge remover triggered 0 times
  0:01:14.135   632M / 18G   INFO    General                 (simplification.cpp        : 330)   Disrupting self-conjugate edges
  0:01:14.163   632M / 18G   INFO   Simplification           (parallel_processing.hpp   : 165)   Running Removing isolated edges
  0:01:15.504   552M / 18G   INFO   Simplification           (parallel_processing.hpp   : 167)   Removing isolated edges triggered 36890 times
  0:01:15.504   552M / 18G   INFO    General                 (simplification.cpp        : 470)   Counting average coverage
  0:01:15.508   552M / 18G   INFO    General                 (simplification.cpp        : 476)   Average coverage = 28.0327
  0:01:15.508   552M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Mismatch Correction
  0:01:15.508   552M / 18G   INFO    General                 (graph_pack.hpp            : 109)   Normalizing k-mer map. Total 587550 kmers to process
  0:01:15.829   552M / 18G   INFO    General                 (graph_pack.hpp            : 111)   Normalizing done
  0:01:20.571   620M / 18G   INFO    General                 (mismatch_shall_not_pass.hp: 189)   Finished collecting potential mismatches positions
  0:01:20.797   528M / 18G   INFO    General                 (mismatch_shall_not_pass.hp: 290)   All edges processed
  0:01:20.801   528M / 18G   INFO    General                 (mismatch_correction.cpp   :  27)   Corrected 1 nucleotides
  0:01:20.802   528M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Contig Output
  0:01:20.802   528M / 18G   INFO    General                 (contig_output_stage.cpp   :  45)   Writing GFA to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/assembly_graph_with_scaffolds.gfa
  0:01:20.925   528M / 18G   INFO    General                 (contig_output.hpp         :  22)   Outputting contigs to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/before_rr.fasta
  0:01:21.100   528M / 18G   INFO    General                 (contig_output_stage.cpp   :  56)   Outputting FastG graph to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/assembly_graph.fastg
  0:01:21.620   528M / 18G   INFO    General                 (contig_output.hpp         :  22)   Outputting contigs to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/final_contigs.fasta
  0:01:21.789   528M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Paired Information Counting
  0:01:21.811   528M / 18G   INFO    General                 (graph_pack.hpp            : 109)   Normalizing k-mer map. Total 587550 kmers to process
  0:01:22.130   528M / 18G   INFO    General                 (graph_pack.hpp            : 111)   Normalizing done
  0:01:22.134   528M / 18G   INFO    General                 (pair_info_count.cpp       : 323)   Min edge length for estimation: 34607
  0:01:22.134   528M / 18G   INFO    General                 (pair_info_count.cpp       : 334)   Estimating insert size for library #0
  0:01:22.134   528M / 18G   INFO    General                 (pair_info_count.cpp       : 190)   Estimating insert size (takes a while)
  0:01:22.271     1G / 18G   INFO    General                 (pair_info_count.cpp       :  39)   Selecting usual mapper
  0:01:24.405     1G / 18G   INFO    General                 (sequence_mapper_notifier.h:  98)   Total 3022500 reads processed
  0:01:24.657     1G / 18G   INFO    General                 (pair_info_count.cpp       : 209)   Edge pairs: 67108864 (rough upper limit)
  0:01:24.657     1G / 18G   INFO    General                 (pair_info_count.cpp       : 213)   1452051 paired reads (48.0414% of all) aligned to long edges
  0:01:24.660   532M / 18G   INFO    General                 (pair_info_count.cpp       : 357)     Insert size = 499.501, deviation = 10.0018, left quantile = 487, right quantile = 512, read length = 130
  0:01:24.709   724M / 18G   INFO    General                 (pair_info_count.cpp       : 374)   Filtering data for library #0
  0:01:24.710   724M / 18G   INFO    General                 (pair_info_count.cpp       :  39)   Selecting usual mapper
  0:01:26.155   724M / 18G   INFO    General                 (sequence_mapper_notifier.h:  98)   Total 3022500 reads processed
  0:01:26.156   724M / 18G   INFO    General                 (pair_info_count.cpp       : 386)   Mapping library #0
  0:01:26.156   724M / 18G   INFO    General                 (pair_info_count.cpp       : 388)   Mapping paired reads (takes a while)
  0:01:26.156   724M / 18G   INFO    General                 (pair_info_count.cpp       : 289)   Left insert size quantile 487, right insert size quantile 512, filtering threshold 2, rounding threshold 0
  0:01:26.162   736M / 18G   INFO    General                 (pair_info_count.cpp       :  39)   Selecting usual mapper
  0:01:27.697   780M / 18G   INFO    General                 (sequence_mapper_notifier.h:  98)   Total 3022500 reads processed
  0:01:27.704   576M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Distance Estimation
  0:01:27.704   576M / 18G   INFO    General                 (distance_estimation.cpp   : 173)   Processing library #0
  0:01:27.704   576M / 18G   INFO    General                 (distance_estimation.cpp   : 149)   Weight Filter Done
  0:01:27.704   576M / 18G   INFO   DistanceEstimator        (distance_estimation.hpp   : 116)   Using SIMPLE distance estimator
  0:01:27.845   588M / 18G   INFO    General                 (distance_estimation.cpp   :  34)   Filtering info
  0:01:27.845   588M / 18G   INFO    General                 (pair_info_filters.hpp     : 242)   Start filtering; index size: 143368
  0:01:27.890   580M / 18G   INFO    General                 (pair_info_filters.hpp     : 263)   Done filtering
  0:01:27.890   576M / 18G   INFO    General                 (distance_estimation.cpp   : 156)   Refining clustered pair information
  0:01:27.919   576M / 18G   INFO    General                 (distance_estimation.cpp   : 158)   The refining of clustered pair information has been finished
  0:01:27.919   576M / 18G   INFO    General                 (distance_estimation.cpp   : 160)   Improving paired information
  0:01:29.068   624M / 18G   INFO   PairInfoImprover         (pair_info_improver.hpp    : 103)   Paired info stats: missing = 47407; contradictional = 850
  0:01:29.523   624M / 18G   INFO   PairInfoImprover         (pair_info_improver.hpp    : 103)   Paired info stats: missing = 4335; contradictional = 94
  0:01:29.523   624M / 18G   INFO    General                 (distance_estimation.cpp   : 103)   Filling scaffolding index
  0:01:29.523   624M / 18G   INFO   DistanceEstimator        (distance_estimation.hpp   : 116)   Using SMOOTHING distance estimator
  0:01:29.611   624M / 18G   INFO    General                 (distance_estimation.cpp   :  34)   Filtering info
  0:01:29.612   624M / 18G   INFO    General                 (pair_info_filters.hpp     : 242)   Start filtering; index size: 16062
  0:01:29.618   624M / 18G   INFO    General                 (pair_info_filters.hpp     : 263)   Done filtering
  0:01:29.618   624M / 18G   INFO    General                 (distance_estimation.cpp   : 182)   Clearing raw paired index
  0:01:29.635   624M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Repeat Resolving
  0:01:29.635   624M / 18G   INFO    General                 (repeat_resolving.cpp      :  69)   Using Path-Extend repeat resolving
  0:01:29.635   624M / 18G   INFO    General                 (launcher.cpp              : 481)   ExSPAnder repeat resolving tool started
  0:01:29.667   648M / 18G   INFO    General                 (launcher.cpp              : 392)   Creating main extenders, unique edge length = 2000
  0:01:29.667   648M / 18G   INFO    General                 (extenders_logic.cpp       : 278)   Estimated coverage of library #0 is 28.0327
  0:01:29.667   648M / 18G   INFO    General                 (extenders_logic.cpp       : 278)   Estimated coverage of library #0 is 28.0327
  0:01:29.670   648M / 18G   INFO    General                 (extenders_logic.cpp       : 475)   Using 1 paired-end library
  0:01:29.670   648M / 18G   INFO    General                 (extenders_logic.cpp       : 476)   Using 1 paired-end scaffolding library
  0:01:29.670   648M / 18G   INFO    General                 (extenders_logic.cpp       : 477)   Using 0 single read libraries
  0:01:29.670   648M / 18G   INFO    General                 (launcher.cpp              : 420)   Total number of extenders is 3
  0:01:29.670   648M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 0 paths from 7608 (0%)
  0:01:30.713   668M / 18G   INFO    General                 (path_extender.hpp         : 883)   Processed 128 paths from 7608 (1%)
  0:01:32.303   692M / 18G   INFO    General                 (path_extender.hpp         : 883)   Processed 256 paths from 7608 (3%)
  0:01:34.152   716M / 18G   INFO    General                 (path_extender.hpp         : 883)   Processed 512 paths from 7608 (6%)
  0:01:34.653   720M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 761 paths from 7608 (10%)
  0:01:35.561   728M / 18G   INFO    General                 (path_extender.hpp         : 883)   Processed 1024 paths from 7608 (13%)
  0:01:36.197   736M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 1522 paths from 7608 (20%)
  0:01:36.297   736M / 18G   INFO    General                 (path_extender.hpp         : 883)   Processed 2048 paths from 7608 (26%)
  0:01:36.316   740M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 2283 paths from 7608 (30%)
  0:01:36.331   740M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 3044 paths from 7608 (40%)
  0:01:36.366   740M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 3805 paths from 7608 (50%)
  0:01:36.376   744M / 18G   INFO    General                 (path_extender.hpp         : 883)   Processed 4096 paths from 7608 (53%)
  0:01:36.401   744M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 4566 paths from 7608 (60%)
  0:01:36.408   744M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 5327 paths from 7608 (70%)
  0:01:36.414   744M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 6088 paths from 7608 (80%)
  0:01:36.418   748M / 18G   INFO    General                 (path_extender.hpp         : 885)   Processed 6849 paths from 7608 (90%)
  0:01:36.422   748M / 18G   INFO    General                 (launcher.cpp              : 234)   Finalizing paths
  0:01:36.422   748M / 18G   INFO    General                 (launcher.cpp              : 236)   Deduplicating paths
  0:01:36.495   748M / 18G   INFO    General                 (launcher.cpp              : 240)   Paths deduplicated
  0:01:36.495   748M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 295)   Removing overlaps
  0:01:36.495   748M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 298)   Sorting paths
  0:01:36.496   748M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 305)   Marking overlaps
  0:01:36.496   748M / 18G   INFO   OverlapRemover           (pe_resolver.hpp           : 130)   Marking start/end overlaps
  0:01:36.601   748M / 18G   INFO   OverlapRemover           (pe_resolver.hpp           : 133)   Marking remaining overlaps
  0:01:36.707   748M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 308)   Splitting paths
  0:01:36.720   748M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 313)   Deduplicating paths
  0:01:36.729   748M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 315)   Overlaps removed
  0:01:36.737   748M / 18G   INFO    General                 (launcher.cpp              : 257)   Paths finalized
  0:01:36.737   748M / 18G   INFO    General                 (launcher.cpp              : 427)   Closing gaps in paths
  0:01:36.746   752M / 18G   INFO    General                 (launcher.cpp              : 455)   Gap closing completed
  0:01:36.752   756M / 18G   INFO    General                 (launcher.cpp              : 286)   Traversing tandem repeats
  0:01:36.780   756M / 18G   INFO    General                 (launcher.cpp              : 296)   Traversed 0 loops
  0:01:36.780   756M / 18G   INFO    General                 (launcher.cpp              : 234)   Finalizing paths
  0:01:36.780   756M / 18G   INFO    General                 (launcher.cpp              : 236)   Deduplicating paths
  0:01:36.783   756M / 18G   INFO    General                 (launcher.cpp              : 240)   Paths deduplicated
  0:01:36.783   756M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 295)   Removing overlaps
  0:01:36.783   756M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 298)   Sorting paths
  0:01:36.784   756M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 305)   Marking overlaps
  0:01:36.784   756M / 18G   INFO   OverlapRemover           (pe_resolver.hpp           : 130)   Marking start/end overlaps
  0:01:36.795   756M / 18G   INFO   OverlapRemover           (pe_resolver.hpp           : 133)   Marking remaining overlaps
  0:01:36.807   756M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 308)   Splitting paths
  0:01:36.809   756M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 313)   Deduplicating paths
  0:01:36.814   756M / 18G   INFO   PEResolver               (pe_resolver.hpp           : 315)   Overlaps removed
  0:01:36.817   756M / 18G   INFO    General                 (launcher.cpp              : 257)   Paths finalized
  0:01:36.817   756M / 18G   INFO    General                 (launcher.cpp              : 534)   ExSPAnder repeat resolving tool finished
  0:01:37.096   632M / 18G   INFO   StageManager             (stage.cpp                 : 132)   STAGE == Contig Output
  0:01:37.096   632M / 18G   INFO    General                 (contig_output_stage.cpp   :  45)   Writing GFA to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/assembly_graph_with_scaffolds.gfa
  0:01:37.226   628M / 18G   INFO    General                 (contig_output.hpp         :  22)   Outputting contigs to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/before_rr.fasta
  0:01:37.397   628M / 18G   INFO    General                 (contig_output_stage.cpp   :  56)   Outputting FastG graph to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/assembly_graph.fastg
  0:01:37.923   636M / 18G   INFO    General                 (contig_output_stage.cpp   :  20)   Outputting FastG paths to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/final_contigs.paths
  0:01:38.152   628M / 18G   INFO    General                 (contig_output_stage.cpp   :  20)   Outputting FastG paths to /data/gpfs/assoc/bch709/spiderman/gee/spades_output//K77/scaffolds.paths
  0:01:38.380   632M / 18G   INFO    General                 (launch.hpp                : 151)   SPAdes finished
  0:01:38.643   132M / 18G   INFO    General                 (main.cpp                  : 109)   Assembling time: 0 hours 1 minutes 38 seconds

===== Assembling finished. Used k-mer sizes: 21, 33, 55, 77


===== Mismatch correction started.

== Processing of contigs


== Running contig polishing tool: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-corrector-core /data/gpfs/assoc/bch709/spiderman/gee/spades_output/mismatch_corrector/contigs/configs/corrector.info /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_contigs.fasta


== Dataset description file was created: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/mismatch_corrector/contigs/configs/corrector.info

/data/gpfs/assoc/bch709/spiderman/gee/spades_output/mismatch_corrector/contigs/configs/log.properties  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  58)   Starting MismatchCorrector, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  59)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.000     4M / 4M    INFO   DatasetProcessor         (dataset_processor.cpp     : 195)   Splitting assembly...
  0:00:00.000     4M / 4M    INFO   DatasetProcessor         (dataset_processor.cpp     : 196)   Assembly file: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_contigs.fasta
  0:00:00.818     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 203)   Processing paired sublib of number 0
  0:00:00.818     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 206)   /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
  0:00:00.818     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 140)   Running bwa index ...: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa index -a is /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_contigs.fasta
[bwa_index] Pack FASTA... 0.07 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.29 seconds elapse.
[bwa_index] Update BWT... 0.06 sec
[bwa_index] Pack forward-only FASTA... 0.05 sec
[bwa_index] Construct SA from BWT and Occ... 1.14 sec
[main] Version: 0.7.12-r1039
[main] CMD: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa index -a is /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_contigs.fasta
[main] Real time: 3.655 sec; CPU: 3.623 sec
  0:00:04.571     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 149)   Running bwa mem ...:/data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa mem  -v 1 -t 32 /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_contigs.fasta /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz  > /data/gpfs/assoc/bch709/spiderman/gee/spades_output/tmp/corrector_pjylxnuu/lib0_QpuVB2/tmp.sam
[main] Version: 0.7.12-r1039
[main] CMD: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa mem -v 1 -t 32 /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_contigs.fasta /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
[main] Real time: 46.713 sec; CPU: 975.503 sec
  0:00:51.355     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 209)   Adding samfile /data/gpfs/assoc/bch709/spiderman/gee/spades_output/tmp/corrector_pjylxnuu/lib0_QpuVB2/tmp.sam
  0:00:57.140    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 1000000reads, flushing
  0:01:01.054    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 2000000reads, flushing
  0:01:04.946    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 3000000reads, flushing
  0:01:08.710    56M / 56M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 4000000reads, flushing
  0:01:13.415    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 5000000reads, flushing
  0:01:17.977    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 6000000reads, flushing
  0:01:22.580    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 7000000reads, flushing
  0:01:27.102    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 8000000reads, flushing
  0:01:30.809    20M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 235)   Processing contigs
  0:01:35.391  1020M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_32_length_100231_cov_28.025910 processed with 10 changes in thread 15
  0:01:35.394  1020M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_31_length_103791_cov_28.020364 processed with 0 changes in thread 31
  0:01:35.653     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_30_length_106493_cov_28.020758 processed with 0 changes in thread 28
  0:01:35.771  1020M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_29_length_111476_cov_27.660661 processed with 6 changes in thread 18
  0:01:36.177     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_27_length_120177_cov_27.622864 processed with 1 changes in thread 23
  0:01:36.185  1016M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_28_length_115951_cov_28.186867 processed with 8 changes in thread 25
  0:01:36.256  1008M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_26_length_123004_cov_27.572543 processed with 0 changes in thread 21
  0:01:36.295  1012M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_25_length_123051_cov_27.713362 processed with 0 changes in thread 6
  0:01:36.508  1000M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_24_length_124056_cov_28.037369 processed with 0 changes in thread 12
  0:01:36.878  1000M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_22_length_135158_cov_27.729584 processed with 1 changes in thread 17
  0:01:36.878  1000M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_23_length_126795_cov_28.705006 processed with 19 changes in thread 10
  0:01:36.906   988M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_21_length_136933_cov_27.768399 processed with 0 changes in thread 29
  0:01:37.728   980M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_20_length_151736_cov_27.759876 processed with 0 changes in thread 14
  0:01:37.791   968M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_19_length_155476_cov_27.867908 processed with 0 changes in thread 22
  0:01:37.936   956M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_18_length_159924_cov_28.024380 processed with 0 changes in thread 5
  0:01:38.561   944M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_17_length_167555_cov_27.962300 processed with 8 changes in thread 24
  0:01:39.188   932M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_16_length_187382_cov_27.878311 processed with 0 changes in thread 2
  0:01:39.657   916M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_34_length_97565_cov_27.827866 processed with 0 changes in thread 31
  0:01:39.974   912M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_39_length_85253_cov_27.310768 processed with 11 changes in thread 21
  0:01:40.081   912M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_33_length_99028_cov_27.570626 processed with 0 changes in thread 15
  0:01:40.139   908M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_14_length_211535_cov_27.807910 processed with 0 changes in thread 27
  0:01:40.195   888M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_15_length_208551_cov_27.780543 processed with 0 changes in thread 1
  0:01:40.207   868M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_38_length_85579_cov_27.934820 processed with 10 changes in thread 25
  0:01:40.250   868M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_36_length_95026_cov_27.726264 processed with 0 changes in thread 18
  0:01:40.272   868M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_40_length_84832_cov_27.739673 processed with 2 changes in thread 6
  0:01:40.308   868M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_13_length_216511_cov_27.573015 processed with 0 changes in thread 19
  0:01:40.354   848M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_35_length_97052_cov_27.949667 processed with 0 changes in thread 28
  0:01:40.381   844M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_37_length_88431_cov_27.806789 processed with 2 changes in thread 23
  0:01:40.503   844M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_42_length_78756_cov_27.292175 processed with 3 changes in thread 17
  0:01:40.519   844M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_41_length_81614_cov_27.881293 processed with 15 changes in thread 12
  0:01:40.541   844M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_12_length_218160_cov_27.809701 processed with 1 changes in thread 30
  0:01:40.556   824M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_44_length_77052_cov_28.014238 processed with 0 changes in thread 29
  0:01:40.718   820M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_43_length_78741_cov_29.515420 processed with 0 changes in thread 10
  0:01:41.175   816M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_11_length_229454_cov_27.764013 processed with 1 changes in thread 4
  0:01:41.458   792M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_46_length_76363_cov_27.402708 processed with 0 changes in thread 22
  0:01:41.460   792M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_47_length_75643_cov_29.253778 processed with 19 changes in thread 5
  0:01:41.479   784M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_45_length_77024_cov_28.199761 processed with 0 changes in thread 14
  0:01:42.178   780M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_48_length_74517_cov_27.268041 processed with 3 changes in thread 24
  0:01:42.489   776M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_49_length_72177_cov_27.475770 processed with 0 changes in thread 2
  0:01:42.653   772M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_64_length_45765_cov_27.452438 processed with 0 changes in thread 29
  0:01:42.662   768M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_65_length_41581_cov_27.623313 processed with 0 changes in thread 10
  0:01:42.797   764M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_50_length_70916_cov_28.097446 processed with 0 changes in thread 31
  0:01:42.868   756M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_63_length_50390_cov_27.559239 processed with 0 changes in thread 30
  0:01:42.884   748M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_60_length_54589_cov_27.195040 processed with 0 changes in thread 23
  0:01:42.976   744M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_58_length_60029_cov_27.926158 processed with 5 changes in thread 19
  0:01:42.986   736M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_51_length_68486_cov_27.836045 processed with 0 changes in thread 21
  0:01:42.996   728M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_61_length_53502_cov_28.034160 processed with 4 changes in thread 17
  0:01:43.033   720M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_53_length_66218_cov_27.581273 processed with 0 changes in thread 27
  0:01:43.049   712M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_62_length_51899_cov_27.889719 processed with 10 changes in thread 12
  0:01:43.081   704M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_57_length_60433_cov_27.708413 processed with 0 changes in thread 6
  0:01:43.086   696M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_66_length_41395_cov_27.518515 processed with 2 changes in thread 4
  0:01:43.121   692M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_55_length_62466_cov_27.183654 processed with 12 changes in thread 25
  0:01:43.130   684M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_59_length_58518_cov_27.619685 processed with 0 changes in thread 28
  0:01:43.187   676M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_68_length_39716_cov_27.294331 processed with 1 changes in thread 5
  0:01:43.238   672M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_69_length_38001_cov_27.395950 processed with 0 changes in thread 14
  0:01:43.269   668M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_54_length_65767_cov_27.560405 processed with 0 changes in thread 1
  0:01:43.282   660M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_52_length_66518_cov_28.325206 processed with 1 changes in thread 15
  0:01:43.423   652M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_67_length_41278_cov_27.997185 processed with 2 changes in thread 22
  0:01:43.561   648M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_56_length_61036_cov_28.728736 processed with 21 changes in thread 18
  0:01:44.009   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_70_length_37821_cov_28.005749 processed with 0 changes in thread 24
  0:01:44.017   636M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_78_length_22525_cov_29.405381 processed with 7 changes in thread 21
  0:01:44.033   636M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_74_length_28325_cov_27.853122 processed with 0 changes in thread 31
  0:01:44.070   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_73_length_30293_cov_27.808148 processed with 0 changes in thread 10
  0:01:44.108   636M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_77_length_24855_cov_29.877593 processed with 5 changes in thread 19
  0:01:44.132   636M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_79_length_22386_cov_28.218656 processed with 9 changes in thread 17
  0:01:44.178   636M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_72_length_32657_cov_27.883579 processed with 0 changes in thread 29
  0:01:44.179   632M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_76_length_25909_cov_28.677145 processed with 2 changes in thread 23
  0:01:44.188   636M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_75_length_28232_cov_23.780572 processed with 0 changes in thread 30
  0:01:44.281   632M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_10_length_299705_cov_27.603101 processed with 0 changes in thread 26
  0:01:44.295   596M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_71_length_36771_cov_29.120565 processed with 9 changes in thread 2
  0:01:45.225   592M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_9_length_319036_cov_27.763183 processed with 0 changes in thread 8
  0:01:45.555   556M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_6_length_334086_cov_27.613208 processed with 0 changes in thread 16
  0:01:45.607   516M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_7_length_333801_cov_27.746605 processed with 0 changes in thread 3
  0:01:45.660   476M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_8_length_329062_cov_27.843978 processed with 3 changes in thread 9
  0:01:45.927   440M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_5_length_335062_cov_27.756111 processed with 0 changes in thread 20
  0:01:48.425   400M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_4_length_431556_cov_27.699452 processed with 1 changes in thread 7
  0:01:51.587   348M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_3_length_542043_cov_27.751257 processed with 0 changes in thread 11
  0:01:53.382   284M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_2_length_613985_cov_27.733595 processed with 0 changes in thread 13
  0:01:54.703   216M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_1_length_666314_cov_27.644287 processed with 4 changes in thread 0
  0:01:54.714   140M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 255)   Gluing processed contigs
  0:01:56.309   128M / 1G    INFO    General                 (main.cpp                  :  72)   Correcting time: 0 hours 1 minutes 56 seconds

== Processing of scaffolds


== Running contig polishing tool: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-corrector-core /data/gpfs/assoc/bch709/spiderman/gee/spades_output/mismatch_corrector/scaffolds/configs/corrector.info /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_scaffolds.fasta


== Dataset description file was created: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/mismatch_corrector/scaffolds/configs/corrector.info

/data/gpfs/assoc/bch709/spiderman/gee/spades_output/mismatch_corrector/scaffolds/configs/log.properties  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  58)   Starting MismatchCorrector, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  59)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.000     4M / 4M    INFO   DatasetProcessor         (dataset_processor.cpp     : 195)   Splitting assembly...
  0:00:00.000     4M / 4M    INFO   DatasetProcessor         (dataset_processor.cpp     : 196)   Assembly file: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_scaffolds.fasta
  0:00:01.572     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 203)   Processing paired sublib of number 0
  0:00:01.572     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 206)   /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
  0:00:01.572     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 140)   Running bwa index ...: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa index -a is /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_scaffolds.fasta
[bwa_index] Pack FASTA... 0.07 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.35 seconds elapse.
[bwa_index] Update BWT... 0.05 sec
[bwa_index] Pack forward-only FASTA... 0.05 sec
[bwa_index] Construct SA from BWT and Occ... 1.16 sec
[main] Version: 0.7.12-r1039
[main] CMD: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa index -a is /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_scaffolds.fasta
[main] Real time: 3.718 sec; CPU: 3.687 sec
  0:00:05.440     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 149)   Running bwa mem ...:/data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa mem  -v 1 -t 32 /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_scaffolds.fasta /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz  > /data/gpfs/assoc/bch709/spiderman/gee/spades_output/tmp/corrector_govedfg7/lib0_QiPZI2/tmp.sam
[main] Version: 0.7.12-r1039
[main] CMD: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-bwa mem -v 1 -t 32 /data/gpfs/assoc/bch709/spiderman/gee/spades_output/misc/assembled_scaffolds.fasta /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
[main] Real time: 46.552 sec; CPU: 969.363 sec
  0:00:52.072     4M / 8M    INFO   DatasetProcessor         (dataset_processor.cpp     : 209)   Adding samfile /data/gpfs/assoc/bch709/spiderman/gee/spades_output/tmp/corrector_govedfg7/lib0_QiPZI2/tmp.sam
  0:00:56.688    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 1000000reads, flushing
  0:01:00.368    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 2000000reads, flushing
  0:01:04.180    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 3000000reads, flushing
  0:01:07.884    52M / 52M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 4000000reads, flushing
  0:01:12.524    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 5000000reads, flushing
  0:01:17.296    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 6000000reads, flushing
  0:01:22.013    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 7000000reads, flushing
  0:01:26.618    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 8000000reads, flushing
  0:01:31.097    60M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 105)   processed 9000000reads, flushing
  0:01:31.452    20M / 60M   INFO   DatasetProcessor         (dataset_processor.cpp     : 235)   Processing contigs
  0:01:36.094     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_31_length_103791_cov_28.020364 processed with 0 changes in thread 31
  0:01:36.140     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_32_length_100231_cov_28.025910 processed with 9 changes in thread 21
  0:01:36.223     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_30_length_106493_cov_28.020758 processed with 0 changes in thread 27
  0:01:36.498     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_29_length_111476_cov_27.660661 processed with 6 changes in thread 25
  0:01:36.780     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_28_length_115951_cov_28.186867 processed with 8 changes in thread 23
  0:01:36.806     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_25_length_123051_cov_27.713362 processed with 0 changes in thread 1
  0:01:36.818     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_27_length_120177_cov_27.622864 processed with 1 changes in thread 11
  0:01:36.997     1G / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_26_length_123004_cov_27.572543 processed with 0 changes in thread 16
  0:01:37.180  1020M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_24_length_124056_cov_28.037369 processed with 0 changes in thread 2
  0:01:37.438  1020M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_22_length_136933_cov_27.768399 processed with 0 changes in thread 14
  0:01:37.476  1012M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_23_length_126795_cov_28.705006 processed with 19 changes in thread 4
  0:01:38.153  1004M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_21_length_151736_cov_27.759876 processed with 0 changes in thread 5
  0:01:38.371   992M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_20_length_155476_cov_27.867908 processed with 0 changes in thread 7
  0:01:38.543   980M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_18_length_159924_cov_28.024380 processed with 0 changes in thread 8
  0:01:38.626   968M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_19_length_155894_cov_28.740105 processed with 0 changes in thread 29
  0:01:39.131   956M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_17_length_167555_cov_27.962300 processed with 8 changes in thread 28
  0:01:39.872   944M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_16_length_187382_cov_27.878311 processed with 0 changes in thread 24
  0:01:40.597   928M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_38_length_85579_cov_27.934820 processed with 10 changes in thread 1
  0:01:40.678   928M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_33_length_99028_cov_27.570626 processed with 0 changes in thread 31
  0:01:40.742   924M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_34_length_97565_cov_27.827866 processed with 0 changes in thread 21
  0:01:40.814   920M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_39_length_85253_cov_27.310768 processed with 11 changes in thread 11
  0:01:40.853   920M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_35_length_97052_cov_27.949667 processed with 0 changes in thread 27
  0:01:40.860   916M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_42_length_78756_cov_27.292175 processed with 3 changes in thread 14
  0:01:40.917   916M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_15_length_208551_cov_27.780543 processed with 0 changes in thread 0
  0:01:40.931   884M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_14_length_211535_cov_27.807910 processed with 0 changes in thread 17
  0:01:40.951   876M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_36_length_95026_cov_27.726264 processed with 0 changes in thread 25
  0:01:40.957   868M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_37_length_88431_cov_27.806789 processed with 2 changes in thread 23
  0:01:41.040   876M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_40_length_84832_cov_27.739673 processed with 1 changes in thread 16
  0:01:41.169   876M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_41_length_81614_cov_27.881293 processed with 15 changes in thread 2
  0:01:41.173   868M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_43_length_77024_cov_28.199761 processed with 0 changes in thread 4
  0:01:41.201   876M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_12_length_218160_cov_27.809701 processed with 1 changes in thread 13
  0:01:41.246   852M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_13_length_216511_cov_27.573015 processed with 0 changes in thread 22
  0:01:41.724   828M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_11_length_229454_cov_27.764013 processed with 1 changes in thread 15
  0:01:41.764   804M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_45_length_76363_cov_27.402708 processed with 0 changes in thread 7
  0:01:41.874   800M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_44_length_76498_cov_28.139778 processed with 10 changes in thread 5
  0:01:42.050   796M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_46_length_75643_cov_29.253778 processed with 19 changes in thread 8
  0:01:42.158   792M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_47_length_74517_cov_27.268041 processed with 3 changes in thread 29
  0:01:42.525   788M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_48_length_72177_cov_27.475770 processed with 0 changes in thread 28
  0:01:43.163   780M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_49_length_70916_cov_28.097446 processed with 0 changes in thread 24
  0:01:43.246   772M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_64_length_41596_cov_27.812616 processed with 2 changes in thread 22
  0:01:43.297   768M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_63_length_45765_cov_27.452438 processed with 0 changes in thread 13
  0:01:43.465   764M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_59_length_54589_cov_27.195040 processed with 0 changes in thread 23
  0:01:43.477   756M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_62_length_50390_cov_27.559239 processed with 0 changes in thread 4
  0:01:43.553   748M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_66_length_41395_cov_27.518515 processed with 2 changes in thread 7
  0:01:43.584   744M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_60_length_53502_cov_28.034160 processed with 4 changes in thread 16
  0:01:43.595   736M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_50_length_68486_cov_27.836045 processed with 0 changes in thread 1
  0:01:43.619   728M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_57_length_60029_cov_27.926158 processed with 5 changes in thread 17
  0:01:43.668   720M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_61_length_51899_cov_27.889719 processed with 10 changes in thread 2
  0:01:43.674   712M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_65_length_41581_cov_27.623313 processed with 0 changes in thread 15
  0:01:43.675   712M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_69_length_32657_cov_27.883579 processed with 0 changes in thread 29
  0:01:43.678   704M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_68_length_38001_cov_27.395950 processed with 0 changes in thread 8
  0:01:43.693   700M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_58_length_58518_cov_27.619685 processed with 0 changes in thread 25
  0:01:43.769   692M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_56_length_60433_cov_27.708413 processed with 0 changes in thread 0
  0:01:43.780   684M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_54_length_62466_cov_27.183654 processed with 12 changes in thread 27
  0:01:43.781   684M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_51_length_66518_cov_28.325206 processed with 1 changes in thread 31
  0:01:43.810   668M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_67_length_41278_cov_27.997185 processed with 2 changes in thread 5
  0:01:43.848   664M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_52_length_66218_cov_27.581273 processed with 0 changes in thread 21
  0:01:43.933   656M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_53_length_65767_cov_27.560405 processed with 0 changes in thread 11
  0:01:43.936   648M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_55_length_61036_cov_28.728736 processed with 23 changes in thread 14
  0:01:43.967   644M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_70_length_30293_cov_27.808148 processed with 0 changes in thread 28
  0:01:44.472   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_71_length_28325_cov_27.853122 processed with 0 changes in thread 24
  0:01:44.552   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_75_length_22525_cov_29.405381 processed with 7 changes in thread 4
  0:01:44.575   644M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_72_length_28232_cov_23.780572 processed with 0 changes in thread 22
  0:01:44.587   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_76_length_22386_cov_28.218656 processed with 9 changes in thread 7
  0:01:44.596   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_73_length_25909_cov_28.677145 processed with 2 changes in thread 13
  0:01:44.628   640M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_10_length_299705_cov_27.603101 processed with 0 changes in thread 12
  0:01:44.665   604M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_74_length_24855_cov_29.877593 processed with 5 changes in thread 23
  0:01:45.712   600M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_9_length_319036_cov_27.763183 processed with 0 changes in thread 3
  0:01:45.940   564M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_5_length_335062_cov_27.756111 processed with 0 changes in thread 9
  0:01:45.962   524M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_8_length_329062_cov_27.843978 processed with 0 changes in thread 20
  0:01:45.966   516M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_7_length_333801_cov_27.746605 processed with 0 changes in thread 18
  0:01:45.967   508M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_6_length_334086_cov_27.613208 processed with 0 changes in thread 30
  0:01:48.627   408M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_4_length_431556_cov_27.699452 processed with 1 changes in thread 6
  0:01:53.901   356M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_3_length_613985_cov_27.733595 processed with 0 changes in thread 10
  0:01:54.850   288M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_2_length_666314_cov_27.644287 processed with 4 changes in thread 26
  0:01:55.437   216M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 251)   Contig NODE_1_length_677246_cov_27.741934 processed with 3 changes in thread 19
  0:01:55.449   140M / 1G    INFO   DatasetProcessor         (dataset_processor.cpp     : 255)   Gluing processed contigs
  0:01:57.018   128M / 1G    INFO    General                 (main.cpp                  :  72)   Correcting time: 0 hours 1 minutes 57 seconds

===== Mismatch correction finished.

 * Corrected reads are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/
 * Assembled contigs are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/contigs.fasta
 * Assembled scaffolds are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/scaffolds.fasta
 * Paths in the assembly graph corresponding to the contigs are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/scaffolds.paths
 * Assembly graph is in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/assembly_graph.fastg
 * Assembly graph in GFA format is in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/assembly_graph_with_scaffolds.gfa

======= SPAdes pipeline finished.

SPAdes log can be found here: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/spades.log

Thank you for using SPAdes!
```

### Assembly statistics

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


### Download or find your results.
```bash
https://www.dropbox.com/s/fo41zymzyb0222p/scaffolds.fasta
https://www.dropbox.com/s/epjj00zpizs56d0/contigs.fasta
https://www.dropbox.com/s/0dgyfmz4gbkdx3q/assembly_graph.fastg
```

```bash
conda install -c bioconda assembly-stats
assembly-stats scaffolds.fasta
assembly-stats contigs.fasta
```
### Assembly statistics result
```bash
stats for scaffolds.fasta
sum = 11241648, n = 2345, ave = 4793.88, largest = 677246
N50 = 167555, n = 17
N60 = 124056, n = 24
N70 = 97565, n = 34
N80 = 72177, n = 48
N90 = 38001, n = 68
N100 = 78, n = 2345
N_count = 308
Gaps = 4
```

```bash
stats for contigs.fasta
sum = 11241391, n = 2349, ave = 4785.61, largest = 666314
N50 = 167555, n = 17
N60 = 123051, n = 25
N70 = 95026, n = 36
N80 = 70916, n = 50
N90 = 36771, n = 71
N100 = 78, n = 2349
N_count = 0
Gaps = 0
```
```bash
cat scaffolds.paths
```
```bash
NODE_1_length_677246_cov_27.741934
5330914+,39250+,5099246-;
4754344-,4601428-,5180688-,1424894+,5327688+,732820-,5237058-,5052460-,4723018+,4800852+,5331930-,732820-,5019006-,5052460-,4755300-,4800852+,5331932-,5060558+,5185654-,5071338-,5178452+,5178460+,5178468+,5178476+,5254862-,5325448+,88806-,5243982-,5053698+,1425522-,5239940+,5238056-,4867204-,5331654+
NODE_1_length_677246_cov_27.741934'
5331654-,4867204+,5238056+,5239940-,1425522+,5053698-,5243982+,88806+,5325448-,5254862+,5178476-,5178468-,5178460-,5178452-,5071338+,5185654+,5060558-,5331932+,4800852-,4755300+,5052460+,5019006+,732820+,5331930+,4800852-,4723018-,5052460+,5237058+,732820+,5327688-,1424894-,5180688+,4601428+,4754344+;
5099246+,39250-,5330914-
NODE_2_length_666314_cov_27.644287
5327204+,103640+,4836832-,4851626+,5361530-,5361524-,5361528-,5329856-,1126012-,5329854-,1126012-,5236812+,5052228+,5236810+,5052228+,5099424-,4479150+,4812968+,4479150+,5062132+,414588+,5051858+,414588+,5331378+,4760742-,4925978+,5327370-,474724-,5094420-,4653402-,5331664-,4961012-,5018412+,5072166+,5040312+,5030606-,4961012-,5236106+,5072166+,5072168+,4915570-,5050466-,4765966-,5059218-,4915570-,5058386-,5236104+,5095538-,5095540-,5095538-,5149466+,4774822-,5017666+,4774822-,5065186-,876886-,4799048-,876886-,5332178-
NODE_2_length_666314_cov_27.644287'
5332178+,876886+,4799048+,876886+,5065186+,4774822+,5017666-,4774822+,5149466-,5095538+,5095540+,5095538+,5236104-,5058386+,4915570+,5059218+,4765966+,5050466+,4915570+,5072168-,5072166-,5236106-,4961012+,5030606+,5040312-,5072166-,5018412-,4961012+,5331664+,4653402+,5094420+,474724+,5327370+,4925978-,4760742+,5331378-,414588-,5051858-,414588-,5062132-,4479150-,4812968-,4479150-,5099424+,5052228-,5236810-,5052228-,5236812-,1126012+,5329854+,1126012+,5329856+,5361528+,5361524+,5361530+,4851626-,4836832+,103640-,5327204-
NODE_3_length_613985_cov_27.733595
5250014-,5121298+,5057128-,4953418-,5238246+,5264238+,5264242+,4468126+,5331520+,4813546-,4676908-,4813546-,5059540-,4862238+,5032536-,4862238+,5045932+,1122610+,4827200-,928516+,5031788-,4629584+,5007546-,1271448+,4907228+,1271448+,5099418-,5331326-,5030236-,5236282-,5100426-,5100418-,5100430-,139162-,4675324-,5354642-,372-,374+,5044194+,5058512+,5325918+,4544022-,4816684-,427838+,5238146+,269904+,117192-
NODE_3_length_613985_cov_27.733595'
117192+,269904-,5238146-,427838-,4816684+,4544022+,5325918-,5058512-,5044194-,374-,372+,5354642+,4675324+,139162+,5100430+,5100418+,5100426+,5236282+,5030236+,5331326+,5099418+,1271448-,4907228-,1271448-,5007546+,4629584-,5031788+,928516-,4827200+,1122610-,5045932-,4862238-,5032536+,4862238-,5059540+,4813546+,4676908+,4813546+,5331520-,4468126-,5264242-,5264238-,5238246-,4953418+,5057128+,5121298-,5250014+
NODE_4_length_431556_cov_27.699452
```
### FASTG file format

```bash
>EDGE_5360468_length_246_cov_13.568047:EDGE_5284398_length_327_cov_11.636000,EDGE_5354800_length_230_cov_14.470588';
GCTTCTTCTTGCTTCTCAAAGCTTTGTTGGTTTAGCCAAAGTCCAGATGAGTCTTTATCT
TTGTATCTTCTAACAAGGAAACACTACTTAGGCTTTTAGGATAAGCTTGCGGTTTAAGTT
TGTATACTCAATCATACACATGACATCAAGTCATATTCGACTCCAAAACACTAACCAAGC
TTCTTCTTGCACCTCAAAGCTTTGTTGGTTTAGCCAAAGTCCATATAAGTCTTTGTCTTT
GTATCT
>EDGE_5360470_length_161_cov_15.607143:EDGE_5332762_length_98_cov_43.619048';
GCTTCTTCTTGCTTCTCAAAGCTTTGTTGGTTTAGCCAAAGTCCAGATGAGTCTTTATCT
TTGTATCTTCTAACAAGAAAACACTACTTACGCTTTTAGGATAATGTTGCGGTTTAAGTT
CTTATACTCAATCATACACATGACATCAAGTCATATTCGAC
>EDGE_5354230_length_92_cov_267.066667':EDGE_5354222_length_86_cov_252.444444',EDGE_5355724_length_1189_cov_26.724820;
AAGCAAAGACTAAGTTTGGGGGAGTTGATAAGTGTGTATTTTGCATGTTTTGAGCATCCA
TTTGTCATCACTTTAGCATCATATCATCACTG
>EDGE_5344586_length_373_cov_22.574324:EDGE_5360654_length_82_cov_117.400000';
GCTAAAGTGATGACAAATGGATGCTCAAAACATGCAAAATACACACTTATCAACTCCCCC
AAACTTAGTCTTTGCTTAAGAACAAGCTGGAGGTGAGGTTTGAAAGCGGGGACTCAGAGC
CAAAGCAGCAGATAAACCAGATGAAATCAATGTCCAAGTTGATAGTTCTAAGTTGCGATA
TGATCGAATTCTACTCAAAAACGTTAGCCATGCCTTTTTATCAATCAATCCGACTCATAT
GCTCGACCTACACGTGTTTTCAAATCTACCAATCCCTTTAACATTCATTAGCTCTAGAAC
GTGAATCAAGCAATGCATCATCAATGAACTCATTTGGCTAAGGTAAAAGGTCAAGAGACA
AAGATGGTCCCTT
>EDGE_5354236_length_91_cov_242.857143:EDGE_5350728_length_80_cov_275.666667';
GCTAAAGTGATGACAAATGGATGCTCAAAACATGCAAAATACACACTTATCAACTCCCCC
AAACTTAGTCTTTGCTTGCCCTCAAGCAAAC

```
![bandage]({{site.baseurl}}/fig/bandage.png)
![assembly_spades]({{site.baseurl}}/fig/assembly_spades.png)

[bandage](https://rrwick.github.io/Bandage/)


## Homework 11/05/2019
Please calculate N50 of `before_rr.fasta` and check how many reads in "BCH709_0001.fastq.gz" send to `wyim@unr.edu`


## PacBio assembly
```bash
conda deactivate
conda activate preprocessing
conda install -c bioconda nanostat nanoplot 
```
### Reads download
```bash
mkdir PacBio
cd PacBio
```

```bash
https://www.dropbox.com/s/7coua2gedbuykl6/BCH709_Pacbio_1.fastq.gz
https://www.dropbox.com/s/fniub0rxv48hupp/BCH709_Pacbio_2.fastq.gz
```

### Check PacBio reads statistics
```bash
NanoStat --fastq BCH709_Pacbio_1.fastq.gz
NanoPlot -t 2 --fastq  BCH709_Pacbio_1.fastq.gz --maxlength 40000 --plots hex dot pauvre -o pacbio_stat
```

### Transfer your result
```bash
*.png *.html *.txt
```

![HistogramReadlength]({{site.baseurl}}/fig/HistogramReadlength.png)
![LengthvsQualityScatterPlot_hex]({{site.baseurl}}/fig/LengthvsQualityScatterPlot_hex.png)



### PacBio reads statistics
```
General summary:        
Mean read length:               9,698.6
Mean read quality:                  6.6
Median read length:             8,854.0
Median read quality:                6.6
Number of reads:               58,497.0
Read length N50:               10,901.0
Total bases:              567,339,600.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:	58497 (100.0%) 567.3Mb
>Q7:	8260 (14.1%) 78.6Mb
>Q10:	0 (0.0%) 0.0Mb
>Q12:	0 (0.0%) 0.0Mb
>Q15:	0 (0.0%) 0.0Mb
Top 5 highest mean basecall quality scores and their read lengths
1:	8.5 (13544)
2:	8.5 (14158)
3:	8.5 (9590)
4:	8.5 (10741)
5:	8.5 (7529)
Top 5 longest reads and their mean basecall quality score
1:	24999 (6.4)
2:	24992 (6.0)
3:	24983 (7.0)
4:	24980 (6.5)
5:	24977 (6.6)
```

### How can we calculate coverage?



## *De novo* assembly
![De-novo-assembly-High-level-diagram-of-long-read-assembly-pipeline-The-assembly-has](De-novo-assembly-High-level-diagram-of-long-read-assembly-pipeline-The-assembly-has.png)

## Compare Assembly
![dotplot2]({{site.baseurl}}/fig/dotplot2.png)


## Canu
Canu (Koren et al. 2017) is a fork of the celera assembler and improves upon the earlier PBcR pipeline into a single, comprehensive assembler. Highly repetitive k-mers, which are abundant in all the reads, can be non-informative. Hence term frequency, inverse document frequency (tf-idf), a weighting statistic was added to MinHashing, giving weightage to non-repetitive k-mers as minimum values in the MinHash sketches, and sensitivity has been demonstrated to reach up to 89% without any parameter adjustment. By retrospectively inspecting the assembly graphs and also statistically filtering out repeat-induced overlaps, the chances of mis-assemblies are reduced.
![canu]({{site.baseurl}}/fig/canu.png)


### Genome assembly Spades
```bash
mkdir Spades_pacbio
cd Spades_pacbio
conda activate genomeassembly

```
### Submit below job

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
zcat  <LOCATION_BCH709_Pacbio_1.fastq.gz> <LOCATION_BCH709_Pacbio_1.fastq.gz> >> merged_pacbio.fastq
spades.py -k 21,33,55,77 --careful -1 <trim_galore output> -2 <trim_galore output> --pacbio merged_pacbio.fastq -o spades_output --memory 140 --threads 64
```



## Canu assembly
```
conda activate genomeassembly
```
### Submit below job
```bash
#!/bin/bash
#SBATCH --job-name=Canu
#SBATCH --cpus-per-task=64
#SBATCH --time=2:00:00
#SBATCH --mem=140g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o Canu.out # STDOUT
#SBATCH -e Canu.err # STDERR

canu -p canu -d canu_outdir genomeSize=11m corThreads=64 -pacbio-raw <LOCATION_BCH709_Pacbio_1.fastq.gz> <LOCATION_BCH709_Pacbio_1.fastq.gz> <LOCATION_BCH709_0003.fastq.gz> useGrid=false
```
