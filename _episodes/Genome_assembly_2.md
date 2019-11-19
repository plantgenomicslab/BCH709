---
layout: page
title: Genome assembly 2
published: true
---
# Meeting schedule
Please feel free to fill up below link to review our course.
Meeting schedule [Link](https://docs.google.com/spreadsheets/d/1c4RzQle8AZPRdayYW5Ov3b16uWKMyUVXl8-iNnuCDSI/edit?usp=sharing)

## HPC 
- S1
- S2
- S3



>## BLAST assignment solution
>
>```bash
>
>#!/bin/bash
>#SBATCH --job-name=BLASTX
>#SBATCH --cpus-per-task=64
>#SBATCH --time=10:00:00
>#SBATCH --mem=140g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=wyim@unr.edu
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>#SBATCH -o Canu.out # STDOUT
>#SBATCH -e Canu.err # STDERR
>```
>conda activate alignment 
>blastx 
>
>
{: .solution}


>## Canu assembly solution
>```
>conda activate genomeassembly
>```
>### Submit below job
>```bash
>#!/bin/bash
>#SBATCH --job-name=Canu
>#SBATCH --cpus-per-task=64
>#SBATCH --time=2:00:00
>#SBATCH --mem=140g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=wyim@unr.edu
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>#SBATCH -o Canu.out # STDOUT
>#SBATCH -e Canu.err # STDERR
>canu -p bch709 -d canu_outdir genomeSize=11m -pacbio-raw <LOCATION_BCH709_Pacbio_02.fastq.gz> <LOCATION_BCH709_Pacbio_01.fastq.gz>executiveMemory=64 gridOptions='--time=12-00:00:00 -p cpu-s2-core-0 -A cpu-s2-bch709-0'
>```
{: .solution}

>## Genome assembly Spades (Illumina + PacBio)
>```bash
>conda activate genomeassembly
>```
>### Submit below job
>
>```bash
>#!/bin/bash
>#SBATCH --job-name=Spades
>#SBATCH --cpus-per-task=64
>#SBATCH --time=2:00:00
>#SBATCH --mem=140g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=wyim@unr.edu
>#SBATCH -o Spades.out # STDOUT
>#SBATCH -e Spades.err # STDERR
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>zcat  <LOCATION_BCH709_Pacbio_1.fastq.gz> <LOCATION_BCH709_Pacbio_1.fastq.gz> >> merged_pacbio.fastq
spades.py -k 21,33,55,77 --careful -1 <trim_galore output> -2 <trim_galore output> --pacbio merged_pacbio.fastq -o spades_output --memory 140 --threads 64
>```
{: .solution}

## Check the Quality of Genome Assembly

```bash
mkdir genomeassembly_results/
```

## Activate Your environment

```bash
conda activate genomeassembly
```

## Copy Your Assembly Results

```bash
assembly-stats canu.contigs.fasta pacbio_illumina_spades.fasta pacbio_spades.fasta
```

## Which one is the best?


>## If you don't have download below link
>```
>https://www.dropbox.com/s/38xjzxptvj2awjv/assembly.tar.gz
>```
{: .solution}


## Install Global Alignmnet Software
```bash
conda install mummer
```
Open source MUMmer 3.0 is described in "Versatile and open software for comparing large genomes." S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot, M. Shumway, C. Antonescu, and S.L. Salzberg, Genome Biology (2004), 5:R12.

MUMmer 2.1, NUCmer, and PROmer are described in "Fast Algorithms for Large-scale Genome Alignment and Comparision." A.L. Delcher, A. Phillippy, J. Carlton, and S.L. Salzberg, Nucleic Acids Research (2002), Vol. 30, No. 11 2478-2483.

MUMmer 1.0 is described in "Alignment of Whole Genomes." A.L. Delcher, S. Kasif, R.D. Fleischmann, J. Peterson, O. White, and S.L. Salzberg, Nucleic Acids Research, 27:11 (1999), 2369-2376.

Space efficent suffix trees are described in "Reducing the Space Requirement of Suffix Trees." S. Kurtz, Software-Practice and Experience, 29(13): 1149-1171, 1999.

## Run Genome Wide Global Alignmnet Software
```bash
#!/bin/bash
#SBATCH --job-name=nucmer
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o nucmer.out # STDOUT
#SBATCH -e nucmer.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0
nucmer  --coords -p canu_spades_pacbio_illumina <canu.contigs.fasta> <pacbio_illumina_spades.fasta>
nucmer  --coords -p canu_spades_illumina <canu.contigs.fasta> <illumina_spades.fasta>
```


### Dot 
Dot is an interactive dot plot viewer for genome-genome alignments.

Dot is publicly available here: https://dnanexus.github.io/dot And can also be used locally by cloning this repository and simply opening the index.html file in a web browser.


After aligning genome assemblies or finished genomes against each other with MUMmer's nucmer, the alignments can be visualized with Dot. Instead of generating static dot plot images on the command-line, Dot lets you interact with the alignments by zooming in and investigating regions in detail.

To prepare a .delta file (nucmer output) for Dot, run this python (3.6) script first: https://dnanexus.github.io/dot/DotPrep.py

The DotPrep.py script will apply a unique anchor filtering algorithm to mark alignments as unique or repetitive. This algorithm analyzes all of the alignments, and it needs to see unfiltered data to determine which alignments are repetitive, so make sure to run nucmer without any filtering options and without running delta-filter on the .delta file before passing it into DotPrep.py.

```bash
wget https://dnanexus.github.io/dot/DotPrep.py
nano DotPrep.py
```


### Improve nanorc
```bash
nano ~/.nanorc
```
```
set nowrap
set softwrap
set const
## Nanorc files
include "/usr/share/nano/nanorc.nanorc"

## C/C++
include "/usr/share/nano/c.nanorc"

## HTML
include "/usr/share/nano/html.nanorc"

## TeX
include "/usr/share/nano/tex.nanorc"

## Quoted emails (under e.g. mutt)
include "/usr/share/nano/mutt.nanorc"

## Patch files
include "/usr/share/nano/patch.nanorc"

## Manpages
include "/usr/share/nano/man.nanorc"

## Groff
include "/usr/share/nano/groff.nanorc"

## Perl
include "/usr/share/nano/perl.nanorc"

## Python
include "/usr/share/nano/python.nanorc"

## Ruby
include "/usr/share/nano/ruby.nanorc"

## Java
include "/usr/share/nano/java.nanorc"

## Assembler
include "/usr/share/nano/asm.nanorc"

## Bourne shell scripts
include "/usr/share/nano/sh.nanorc"

## POV-Ray
include "/usr/share/nano/pov.nanorc"
```

```bash
nano DotPrep.py
```

```bash
#!/bin/bash
#SBATCH --job-name=dot
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o nucmer.out # STDOUT
#SBATCH -e nucmer.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0
python DotPrep.py  --delta canu_spades_pacbio_illumina.delta
python DotPrep.py  --delta canu_spades_illumina.delta
```

The output of DotPrep.py includes the \*.coords and \*.coords.idx that should be used with Dot for visualization.


## Visualization
- Transfer \*.coords.\* files
- Go to  https://dnanexus.github.io/dot/

![dotplot4]({{site.baseurl}}/fig/dotplot4.png)



## Which one is the best?



## BUSCO
BUSCO assessments are implemented in open-source software, with a large selection of lineage-specific sets of Benchmarking Universal Single-Copy Orthologs. These conserved orthologs are ideal candidates for large-scale phylogenomics studies, and the annotated BUSCO gene models built during genome assessments provide a comprehensive gene predictor training set for use as part of genome annotation pipelines.

https://busco.ezlab.org/v2/


```bash
conda install busco
```

```bash
https://busco.ezlab.org/v2/datasets/embryophyta_odb9.tar.gz
```

```bash 
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o busco.out # STDOUT
#SBATCH -e busco.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

export AUGUSTUS_CONFIG_PATH="~/miniconda3/envs/genomeassembly/config/

run_busco -i <canu.contigs.fasta> --cpu 16  -o canu  -l embryophyta_odb9 -m geno
run_busco -i <pacbio_illumina_spades> --cpu 16  -o pacbio_illumina_spades  -l embryophyta_odb9 -m geno
run_busco -i <pacbio_spades> --cpu 16  -o pacbio_spades  -l embryophyta_odb9 -m geno

```



### Investigate taxa
```bash
conda install kraken kraken2 
```
We will be using a tool called Kraken2. This tool uses k-mers to assign a taxonomic labels in form of NCBI Taxonomy to the sequence (if possible). The taxonomic label is assigned based on similar k-mer content of the sequence in question to the k-mer content of reference genome sequence. The result is a classification of the sequence in question to the most likely taxonomic label. If the k-mer content is not similar to any genomic sequence in the database used, it will not assign any taxonomic label.

We can also use another tool by the same group called Centrifuge. This tool uses a novel indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index, optimized specifically for the metagenomic classification problem to assign a taxonomic labels in form of NCBI Taxonomy to the sequence (if possible). The result is a classification of the sequence in question to the most likely taxonomic label. If the search sequence is not similar to any genomic sequence in the database used, it will not assign any taxonomic label.


## Pilon
```bash
conda install -c bioconda pilon
```
Pilon uses read alignment analysis to diagnose, report, and automatically improve de novo genome assemblies as well as call variants.
Pilon then outputs a FASTA file containing an improved representation of the genome from the read data and an optional VCF file detailing variation seen between the read data and the input genome.

To aid manual inspection and improvement by an analyst, Pilon can optionally produce tracks that can be displayed in genome viewers such as IGV and GenomeView, and it reports other events (such as possible large collapsed repeat regions) in its standard output.

```bash 
#!/bin/bash
#SBATCH --job-name=Illumina_mapping
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o mapping.out # STDOUT
#SBATCH -e mapping.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

bwa index canu.contigs.fasta
bwa mem canu.contigs.fasta  <illumina_R1> <illumina_R2>  -o canu_illumina.sam
samtools view -Sb canu_illumina.sam -o canu_illumina.bam
samtools sort canu_illumina.bam -o canu_illumina_sort.bam
samtools index canu_illumina_sort.bam

```

```bash
nano /data/gpfs/home/<YOURID>/miniconda3/envs/genomeassembly/bin/pilon 


default_jvm_mem_opts = ['-Xms8g', '-Xmx40g']
```

```bash 
#!/bin/bash
#SBATCH --job-name=Illumina_mapping
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o mapping.out # STDOUT
#SBATCH -e mapping.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

pilon --genome canu.contigs.fasta --frags canu_illumina_sort.bam --output canu.illumina  --vcf --changes
```






### How can we improve these genome assemblies?

