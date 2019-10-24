---
layout: page
title: Midterm_Review
published: true
---

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
