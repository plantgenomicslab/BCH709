---
layout: page
title: Genome assembly
published: true
---

## Genome assembly


![Genome assembly]({{site.baseurl}}/fig/genome_assembly.png))

After DNA sequencing is complete, the fragments of DNA that come out of the machine are all jumbled up. Like a jigsaw puzzle we need to take the pieces of the genome and put them back together.

### What’s the challenge?
- The technology of DNA sequencing? is not 100 per cent accurate and therefore there are likely to be errors in the DNA sequence that is produced.  

- So, to account for the errors that could potentially occur, each base in the genome? is sequenced a number of times over, this is called coverage. For example, 30 times (30-fold) coverage means each base? is sequenced 30 times.  

- Effectively, the more times you sequence, or “read”, the same section of DNA, the more confidence you have that the final sequence is correct.
30- to 50-fold coverage is currently the standard used when sequencing human genomes to a high level of accuracy.  

- During the Human Genome Project? coverage was only between 5- and 10-fold and used a different sequencing technology to those used today.  

- Coverage has increased because of a few reasons:
    Although most current sequencing techniques are now faster than they were during the Human Genome Project, some sequencing technologies have a higher error rate.
    Some sequencing technologies deal with shorter reads of DNA which means that gaps are more likely to occur when the genome is assembled. Having a higher coverage reduces the likelihood of there being gaps in the final assembled sequence.
    It is also much cheaper to carry out sequencing to a higher coverage than it was at the time of the Human Genome Project.

- High coverage means that after sequencing DNA we have lots and lots of pieces of DNA sequence (reads).

- To put this into perspective, once a human genome has been fully sequenced we have around 100 gigabases (100,000,000,000 bases) of sequence data.

- Like the pieces of a jigsaw puzzle, these DNA reads are jumbled up so we need to piece them together and put them in the correct order to assemble the genome sequence.

### What do we need to do?

- Put the pieces together in the correct order to construct the complete genome sequence and identify any areas of interest.  
- This is done using processes called alignment and assembly:  
    Alignment is when the new DNA sequence is compared to existing DNA sequences to find any similarities or discrepancies between them and then arranged to show these features. Alignment is a vital part of assembly.

    Assembly involves taking a large number of DNA reads, looking for areas in which they overlap with each other and then gradually piecing together the ‘jigsaw’. It is an attempt to reconstruct the original genome. This is primarily carried out for de novo sequences?.

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

Genome assembly refers to the process of taking a large number of short DNA sequences and putting them back together to create a representation of the original chromosomes from which the DNA originated. De novo genome assemblies assume no prior knowledge of the source DNA sequence length, layout or composition. In a genome sequencing project, the DNA of the target organism is broken up into millions of small pieces and read on a sequencing machine. These “reads” vary from 20 to 1000 nucleotide base pairs (bp) in length depending on the sequencing method used. Typically for Illumina type short read sequencing, reads of length 36 - 150 bp are produced. These reads can be either `single ended` as described above or `paired end`. A good summary of other types of DNA sequencing can be found below.

![Genome assembly]({{site.baseurl}}/fig/genome_assembly2.png))

Paired end reads are produced when the fragment size used in the sequencing process is much longer (typically 250 - 500 bp long) and the ends of the fragment are read in towards the middle. This produces two “paired” reads. One from the left hand end of a fragment and one from the right with a known separation distance between them. (The known separation distance is actually a distribution with a mean and standard deviation as not all original fragments are of the same length.) This extra information contained in the paired end reads can be useful for helping to tie pieces of sequence together during the assembly process.

The goal of a sequence assembler is to produce long contiguous pieces of sequence (contigs) from these reads. The contigs are sometimes then ordered and oriented in relation to one another to form scaffolds. The distances between pairs of a set of paired end reads is useful information for this purpose.


The mechanisms used by assembly software are varied but the most common type for short reads is assembly by de Bruijn graph. See this document for an explanation of the de Bruijn graph genome assembler “Velvet.”

Genome assembly is a very difficult computational problem, made more difficult because many genomes contain large numbers of identical sequences, known as repeats. These repeats can be thousands of nucleotides long, and some occur in thousands of different locations, especially in the large genomes of plants and animals.

### Why do we want to assemble an organism’s DNA?
Determining the DNA sequence of an organism is useful in fundamental research into why and how they live, as well as in applied subjects. Because of the importance of DNA to living things, knowledge of a DNA sequence may be useful in practically any biological research. For example, in medicine it can be used to identify, diagnose and potentially develop treatments for genetic diseases. Similarly, research into pathogens may lead to treatments for contagious diseases.


The protocol in a nutshell:
- Obtain sequence read file(s) from sequencing machine(s).
- Look at the reads - get an understanding of what you’ve got and what the quality is like.
Raw data cleanup/quality trimming if necessary.
Choose an appropriate assembly parameter set.
Assemble the data into contigs/scaffolds.
Examine the output of the assembly and assess assembly quality



## Conda
>Conda is a package manager, which helps you find and install packages such as numpy or scipy. It also serves as an environment manager, and allows you to have multiple isolated environments for different projects on a single machine. Each environment has its own installation directories, that doesn’t share packages with other environments.
>
>For example, you need python 2.7 and Biopython 1.60 in project A, while you also work on another project B, which needs python 3.5 and Biopython 1.68. You can use conda to create two separate environments for each project, and you can switch between different versions of packages easily to run your project code.
{: .callout}

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
 ent is installed in your home directory on computer, /home/\<your username\>. The newly created environment <test> will be installed in the directory /home/wyim/miniconda3/envs/test*

### Install from bioconda channel (example: stringtie)
Bioconda is another channel of conda, focusing on bioinformatics software. Instead of adding “-c” to search a channel only one time, “add channels” tells Conda to always search in this channel, so you don’t need to specify the channel every time. Remember to add channel in this order, so that Bioconda channel has the highest priority. Channel orders will be explained in next part.
```bash
 $ conda config --add channels conda-forge
 $ conda config --add channels defaults
 $ conda config --add channels r
 $ conda config --add channels bioconda
```

#### List all package installed
This will show all the packages and versions you’ve installed.
```bash
$ conda list
```


>## RNA-Seq
>Sequence based assays of transcriptomes (RNA-seq) are in wide use because of their favorable properties for quantification, transcript discovery and splice isoform identification, as well as adaptability for numerous more specialized measurements. RNA-Seq studies present some challenges that are shared with prior methods such as microarrays and SAGE tagging, and they also present new ones that are specific to high-throughput sequencing platforms and the data they produce. This document is part of an ongoing effort to provide the community with standards and guidelines that will be updated as RNASeq matures and to highlight unmet challenges. The intent is to revise this document periodically to capture new advances and increasingly consolidate standards and best practices.
>
>RNA-Seq experiments are diverse in their aims and design goals, currently including multiple types of RNA isolated from whole cells or from specific sub-cellular compartments or biochemical classes, such as total polyA+ RNA, polysomal RNA, nuclear ribosome-depleted RNA, various size fractions of RNA and a host of others. The goals of individual experiments range from major transcriptome “discovery” that seeks to define and quantify all RNA species in a starting RNA sample to experiments that simply need to detect significant changes in the more abundant RNA classes across many samples.  
{: .prereq}

### Seven stages to data science
1. Define the question of interest
2. Get the data
3. Clean the data
4. Explore the data
5. Fit statistical models
6. Communicate the results
7. Make your analysis reproducible



>## de Bruijn graph (For transcriptome assembly)
>1. de Bruijn graph construction
> - Draw (by hand) the de Bruijn graph for the following reads using k=3 (assume all reads are from the forward strand, no sequencing errors)  
>AGT  
>ATG  
>GTA  
>TAT  
>TGA  
>GAT  
{: .solution}

### Typical Trinity Command Line

A typical Trinity command for assembling non-strand-specific RNA-seq data would be like so, running the entire process on a single high-memory server (aim for \~1G RAM per \~1M \~76 base Illumina paired reads, but often *much* less memory is required):

Run Trinity like so:

     Trinity --seqType fq --max_memory 50G \
             --left reads_1.fq.gz  --right reads_2.fq.gz --CPU 6

If you have multiple sets of fastq files, such as corresponding to multiple tissue types or conditions, etc., you can indicate them to Trinity like so:

     Trinity --seqType fq --max_memory 50G  \
             --left condA_1.fq.gz,condB_1.fq.gz,condC_1.fq.gz \
             --right condA_2.fq.gz,condB_2.fq.gz,condC_2.fq.gz \
             --CPU 6



### Trinity run
```
pwd  
## should be /data/gpfs/assoc/bch709/<YOUR_NAME>/rnaseq/Trinity  

nano <JOBNAME>.sh  

```

```
#!/bin/bash
#SBATCH --job-name=<JOBNAME>
#SBATCH --cpus-per-task=64
#SBATCH --time=2:00:00
#SBATCH --mem=100g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=<YOUR ID>@nevada.unr.edu
#SBATCH -o <JOBNAME>.out # STDOUT
#SBATCH -e <JOBNAME>.err # STDERR

Trinity --seqType fq  --CPU 64 --max_memory 100G --left /data/gpfs/assoc/bch709/<YOUR_NAME>/rnaseq/rawdata_rnaseq/fastq/pair1.fastq.gz --right /data/gpfs/assoc/bch709/<YOUR_NAME>/rnaseq/rawdata_rnaseq/fastq/pair2.fastq.gz --no_normalize_reads 
```

### Submit job
```bash
sbatch <JOBNAME>.sh 
```


### RNA-Seq reads Count Analysis job script
```
#!/bin/bash
#SBATCH --job-name=<JOBNAME>
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o <JOBNAME>.out # STDOUT
#SBATCH -e <JOBNAME>.err # STDERR

align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left /data/gpfs/assoc/bch709/<YOUR_NAME>/rnaseq/rawdata_rnaseq/fastq/pair1.fastq.gz --right /data/gpfs/assoc/bch709/<YOUR_NAME>/rnaseq/rawdata_rnaseq/fastq/pair2.fastq.gz  --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir rsem_outdir  --thread_count  16
```

### Job submission
```bash
sbatch <JOBNAME>.sh  
```


### RSEM results check
```bash
less rsem_outdir/RSEM.genes.results
```

**Your current location is***
/data/gpfs/assoc/bch709/spiderman/rnaseq/DEG  


#### FPKM
![FPKM]({{site.baseurl}}/fig/FPKM.png)

X = mapped reads count
N = number of reads
L = Length of transcripts

```bash
head -n 2 rsem_outdir/RSEM.genes.results
```


#### Reads count
```bash
samtools flagstat rsem_outdir/bowtie2.bam
```
### Call Python
```bash
python
```


```python
X = 861
Number_Reads_mapped = 1485483
Length = 3475
fpkm= X*(1000/Length)*(1000000/Number_Reads_mapped)
fpkm
```

#### ten to the ninth power = 10\*\*9


### TPM
 Transcripts Per Million

![TPM]({{site.baseurl}}/fig/TPM.png)

![TPM2]({{site.baseurl}}/fig/TPM2.png)


### Sum of FPKM
```bash
cat rsem_outdir/RSEM.genes.results | egrep -v FPKM | awk '{ sum+=$7} END {print sum}'
```
### TPM calculation from FPKM

```python
FPKM = 180.61
SUM_FPKM = 646089
TPM=(FPKM/SUM_FPKM)*10**6
TPM
```

### TPM calculation from reads count
```bash
cat rsem_outdir/RSEM.genes.results | egrep -v FPKM | awk '{ sum+=$5/$3} END {print sum}'
```

```python
sum_count_per_length = 714.037

TPM = (X/Length)*(1/sum_count_per_length )*10**6

```


## File transfer
- There are a number of ways to transfer data to and from HPC clusters. Which you should use depends on several factors, including the ease of use for you personally, connection speed and bandwidth, and the size and number of files which you intend to transfer. Most common options include scp, rsync (command line) and SCP and SFTP clients (GUI). scp (secure copy) is a simple way of transferring files between two machines that use the SSH (Secure SHell) protocol. You may use scp to connect to any system where you have SSH (login) access. scp is available as a protocol choice in some graphical file transfer programs and also as a command line program on most Linux, UNIX, and Mac OS X systems. scp can copy single files, but will also recursively copy directory contents if given a directory name. scp can be used as follows:
```
    scp sourcefile username@pronghorn.rc.unr.edu:somedirectory/
```
    (to a remote system from local)
```
    scp username@pronghorn.rc.unr.edu:somedirectory/sourcefile destinationfile
```
    (from a remote system to local)
```
    scp -r SourceDirectory/ username@pronghorn.rc.unr.edu:somedirectory/
```    
    (***recursive*** directory copy to a remote system from local) 
