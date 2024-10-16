---
layout: page
title: midterm_review
published: true
---

# Midterm 
- Six questions  
- Thursday (10:00AM) to Monday (2PM)


## What is file permission?
```
 -    (rw-)   (rw-)   (r--) 1 john sap
|      |       |       |
type  owner  group   others
```

## Conda  
## Conda  
- Dependencies is one of the main reasons to use Conda.  
Sometimes, install a package is not as straight forward as you think. Imagine a case like this: You want to install package Matplotlib, when installing, it asks you to install Numpy, and Scipy, because Matplotlib need these Numpy and Scipy to work. They are called the dependencies of Matplotlib. For Numpy and Scipy, they may have their own dependencies. These require even more packages.  

## Conda env clean  

```bash
conda clean --all
```
### Conda create enviroment  
```bash
conda create -n review python=3  
```

### Conda activate enviroment  
```bash
conda activate review  
```

### Install software  
```
conda install -c bioconda -c anaconda trinity samtools multiqc fastqc rsem jellyfish bowtie2 salmon trim-galore fastqc bioconductor-ctc bioconductor-deseq2 bioconductor-edger bioconductor-biobase  bioconductor-qvalue r-ape  r-gplots  r-fastcluster
conda install -c anaconda openblas
conda install nano
conda install -c eumetsat tree
conda install -c lmfaber transrate
```


### Check installation  
```
conda list
```
### Conda installation  
```bash
search on web browser [software name & conda] on Google ex: 'hisat2 conda'  
```
```bash
conda install [package name]
```
### Conda export  
export your environment to rnaseq.yaml
```bash
conda env export  > [Name].yaml
```

## Sequencing  
### Illumina sequencing  
![Illumina](https://www.1010genome.com/wp-content/uploads/2018/06/image002.jpg)  

### PacBio sequencing
![PacBio](https://www.pacb.com/wp-content/uploads/Infographic_SMRT-Sequencing-How-it-Works.png)  



## File format  
### Fasta file  
![fasta](https://www.researchgate.net/profile/Morteza-Hosseini-6/publication/309134977/figure/fig1/AS:417452136648705@1476539753111/A-sample-of-the-Multi-FASTA-file.png)
### Fastq file  
![Fastq_file]({{{site.baseurl}}/fig/fastq.png)
![basequality]({{{site.baseurl}}/fig/basequality.png)

### GFF format  
![GFF](https://learn.gencore.bio.nyu.edu/wp-content/uploads/2018/01/Screen-Shot-2018-01-07-at-10.08.12-PM-768x276.png)
1. Sequence ID  
2. Source  
- Describes the algorithm or the procedure that generated this feature. Typically Genescane or Genebank, respectively.  
3. Feature Type  
- Describes what the feature is (mRNA, domain, exon, etc.). These terms are constrained to the [Sequence Ontology terms](http://www.sequenceontology.org/).  
4. Feature Start  
5. Feature End  
6. Score  
- Typically E-values for sequence similarity and P-values for predictions.  
7. Strand  
8. Phase  
- Indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
9 .Atributes   
- A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent . You can see the full list [here](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).   

## RNA Sequencing  
![RNA Sequencing]({{{site.baseurl}}/fig/rnaseq.png)
1. The transcriptome is spatially and temporally dynamic  
2. Data comes from functional units (coding regions)  
3. Only a tiny fraction of the genome  


>## Paper reading  
>Please read this paper  
>https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8  
{: .prereq}  


## Basic Unix/Linux command  

### cd
```
cd /data/gpfs/assoc/bch709-4/<YOUR_ID>
```
### mkdir
```
mkdir BCH709_midterm
cd BCH709_midterm
```

### pwd
```bash
pwd
```

### wget
file download
```bash
wget https://www.dropbox.com/s/yqvfm70yz79jvij/fasta.zip https://www.dropbox.com/s/jjz6aip3euh0d7q/fastq.tar
wget https://www.dropbox.com/s/szzyb3l4243xcsu/bch709.py 
```

### Decompress tar file
```
tar xvf fastq.tar

ls
```
### Decompress zip file
```
unzip fasta.zip

ls
```
### gz file  

### zcat  

### pipe   

### wc  

### rm  

## RNA-Seq  
![RNA Sequencing workflow]({{{site.baseurl}}/fig/rnaseq_workflow.png)

## Advanced bioinformatics tools  
### Seqkit  
https://plantgenomicslab.github.io/BCH709/seqkit_tutorial/index.html

## Trim-Galore  


### References:
- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

