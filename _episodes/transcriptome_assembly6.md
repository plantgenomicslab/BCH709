---
layout: page
title: 9_Transcriptome Assembly
published: true
---

## Trinity installation

```bash
conda clean --all -y

conda create -n transcriptome_assembly python=3.6

conda activate transcriptome_assembly

conda install -y -c anaconda boost

conda install -y -c bioconda -c conda-forge salmon=0.9.1

conda install -y -c bioconda samtools openssl=1.0 bowtie2 bowtie

conda install -y -c r -c conda-forge -c anaconda -c bioconda  bioconductor-ctc bioconductor-deseq2 bioconductor-edger bioconductor-biobase  bioconductor-qvalue  r-ape  r-gplots  r-fastcluster

conda install -y -c bioconda trinity -y
```

## Processing location and input files

```bash
mkdir -p /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trim

cd /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trim

cp /data/gpfs/assoc/bch709-1/Course_material/2020/RNASeq_trimmed_fastq/*.gz .

cd ../

pwd
## current directory is "/data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/"

ls
```

## Create job submission script

```bash
nano trinity.sh
```
**Please change <SOMETHING> to your input**

```bash
#!/bin/bash
#SBATCH --job-name=<TRINITY>
#SBATCH --time=10:15:00
#SBATCH --account=cpu-s2-bch709-1
#SBATCH --partition=cpu-s2-core-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=80g
#SBATCH --mail-type=fail
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH --output=<TRINITY>.out


Trinity --seqType fq  --CPU 8 --max_memory 80G --left trim/DT1_R1_val_1.fq.gz,trim/DT2_R1_val_1.fq.gz,trim/DT3_R1_val_1.fq.gz,trim/WT1_R1_val_1.fq.gz,trim/WT2_R1_val_1.fq.gz,trim/WT3_R1_val_1.fq.gz --right trim/DT1_R2_val_2.fq.gz,trim/DT2_R2_val_2.fq.gz,trim/DT3_R2_val_2.fq.gz,trim/WT1_R2_val_2.fq.gz,trim/WT2_R2_val_2.fq.gz,trim/WT3_R2_val_2.fq.gz
```

### Job submission

```
sbatch trinity.sh
```


### Please check the result
```bash

cd /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trinity_out_dir/

egrep -c ">" Trinity.fasta

TrinityStats.pl /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trinity_out_dir/Trinity.fasta  > <YOURID>.trinity.stat

cat <YOURID>.trinity.stat
```

