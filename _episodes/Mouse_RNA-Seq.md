---
layout: page
title: 	Mouse RNA-Seq
published: true
---
# Environment creation
```bash
conda create -n BCH709_RNASeq -c bioconda -c conda-forge -c anaconda python=3.7 mamba 

conda activate BCH709_RNASeq

mamba install -c bioconda -c conda-forge -c anaconda trim-galore=0.6.7 sra-tools=2.11.0 STAR htseq=1.99.2 subread=2.0.1 multiqc=1.11 snakemake=7.5.0 parallel-fastq-dump=0.6.7 bioconductor-tximport samtools=1.14 r-ggplot2 trinity=2.13.2 hisat2 bioconductor-qvalue sambamba graphviz gffread tpmcalculator lxml rsem
```

# Mouse RNA-Seq
![](https://i.imgur.com/rqNc6OA.png)

https://www.sciencedirect.com/science/article/pii/S2211124722011111


## Slurm

![](https://i.imgur.com/LynACgh.png)


![](https://i.imgur.com/XEMRbJe.png)


[CheatSheet](https://slurm.schedmd.com/pdfs/summary.pdf)


*Slurm provides resource management for the processors allocated to a job, so that multiple job steps can be simultaneously submitted and queued until there are available resources within the job's allocation.*


## Working directory (Pronghorn)

```bash
echo $USER

cd /data/gpfs/assoc/bch709-3/${USER}

mkdir mouse
mkdir /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq
mkdir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref 
mkdir /data/gpfs/assoc/bch709-3/${USER}/mouse/trim
mkdir /data/gpfs/assoc/bch709-3/${USER}/mouse/bam
mkdir /data/gpfs/assoc/bch709-3/${USER}/mouse/readcount
mkdir /data/gpfs/assoc/bch709-3/${USER}/mouse/DEG
```

## Reference Download
![]({{{site.baseurl}}/fig/mouse_ref.png)


https://hgdownload.soe.ucsc.edu/downloads.html


```bash
## change working directory
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/ref 

## download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz

## decompress
gunzip mm39.fa.gz
gunzip refGene.gtf.gz

## Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-3/${USER}/mouse/ref/ref_build.sh
```

## STAR aligner reference build example
![]({{{site.baseurl}}/fig/star.jpeg)
```bash!
STAR  --runThreadN [CPU] --runMode genomeGenerate --genomeDir . --genomeFastaFiles [GENOMEFASTA] --sjdbGTFfile [GENOME_GTF]  --sjdbOverhang 99   --genomeSAindexNbases 12
```

Bioinformatics. 2013 Jan; 29(1): 15–21 [10.1093/bioinformatics/bts635](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)

## STAR aligner reference build on Pronghorn
```bash
#open text editor
nano ref_build.sh

# Add below command to ref_build.sh

STAR  --runThreadN 4 --runMode genomeGenerate --genomeDir . --genomeFastaFiles mm39.fa --sjdbGTFfile refGene.gtf  --sjdbOverhang 99   --genomeSAindexNbases 12
```

## Submit job to HPC
```bash
#submit job
sbach ref_build.sh

#check job
squeue
```

## FASTQ file 
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/

## Link file (without copy)
ln -s /data/gpfs/assoc/bch709-3/Course_materials/mouse/fastq/* /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq

ls /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq
```

## Create file list
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/

ls -1 /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq/*.gz 

ls -1 /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq/*.gz | sed 's/_R.\.fastq\.gz//g'

ls -1 /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq/*.gz | sed 's/_R.\.fastq\.gz//g' | sort -u

ls -1 /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq/*.gz | sed 's/_R.\.fastq\.gz//g' | sort -u >> ../filelist

cat ../filelist
```
## Regular expression

https://regex101.com/

## Trim reads
```bash
trim_galore --paired  --three_prime_clip_R1 [integer] --three_prime_clip_R2 [integer]  --cores [integer]   --max_n [integer]   --fastqc --gzip -o /data/gpfs/assoc/bch709-3/${USER}/mouse/trim {READ_R1} {READ_R2}
```

### example
```bash
trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-3/${USER}/mouse/trim {READ_R1} {READ_R2}
```

## Prepare templet
```bash
cp /data/gpfs/assoc/bch709-3/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq/trim.sh
```

## Edit templet
```bash
nano /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq/trim.sh
```

## Batch submission
```bash
# Check file list
cat ../filelist


# Loop file list
## Add Forward read to variable
## Add reverse read from forward read name substitution

for i in `cat ../filelist`
    do

    read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        echo $read1 $read2
done


# Loop file list
## add file name from variable to trim-galore

for i in `cat ../filelist`
    do
        read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        echo $read1 $read2
        echo "trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-3/${USER}/mouse/trim $read1 $read2" 
done

## merge trim-galore command and trim.sh
for i in `cat ../filelist`
    do
        read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        echo $read1 $read2
        echo "trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-3/${USER}/mouse/trim $read1 $read2" | cat trim.sh - 
done


## add trim-galore command and trim.sh to new file
for i in `cat ../filelist`
    do
        read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        echo $read1 $read2
        echo "trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-3/${USER}/mouse/trim $read1 $read2" | cat trim.sh - > ${i}_trim.sh
done
```
## Batch submission
```bash
ls *.sh
ls -1 *.sh

## Loop *.sh printing
for i in `ls -1 *.sh`
do
    echo $i
done

## Loop *.sh submission
for i in `ls -1 *.sh`
do
    sbatch $i
done
```

## Check submission
```bash
squeue -u ${USER}
```


## RNA-Seq Alignment
```bash
### Move to trim folder
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/trim

### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/mapping.sh

### Edit templet
nano mapping.sh
```

## Check output
```bash
ls -algh /data/gpfs/assoc/bch709-3/${USER}/mouse/trim
```

## Output example
```bash
[FILENAME]_R1_val_1.fq.gz [FILENAME]_R2_val_2.fq.gz
```

## STAR RNA-Seq read alignment example
```bash
STAR --runMode alignReads --runThreadN [CPU_NUMBER] --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --quantMode TranscriptomeSAM GeneCounts --genomeDir [GENOME_DIR] --readFilesCommand gunzip -c --readFilesIn [FORWARD_READ]  [REVERSE_READ] --outSAMtype BAM SortedByCoordinate --outFileNamePrefix [BAMFILENAME_LOCATION]
```

## STAR RNA-Seq alignment
```bash
STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --quantMode TranscriptomeSAM GeneCounts --genomeDir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/[FILENAME]_R1_val_1.fq.gz  /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/[FILENAME]_R2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/[FILENAME].bam
```


## STAR RNA-Seq alignment batch file test
```bash
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read1} /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/${i}.bam"
    done
```

## STAR RNA-Seq alignment batch file 

```bash
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read1} /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/${i}.bam" | cat mapping.sh - > ${i}_mapping.sh
    done
```


## Job submission dependency

```bash
squeue --noheader --format %i --name [YOUR_PREVIOUS_JOBNAME] 
squeue --noheader --format %i --name Trim | tr '\n'  ':'
```


##  Job submission dependency on Trim
```bash
for i in `ls -1 *.sh`
do
    sbatch --dependency=afterany:$(squeue --noheader --format %i --name Trim | tr '\n'  ':')1 $i
done

```


## FeatureCounts
[Bioinformatics, Volume 30, Issue 7, 1 April 2014, Pages 923–930](https://doi.org/10.1093/bioinformatics/btt656)
![]({{{site.baseurl}}/fig/featurecount.png)

```bash!
featureCounts -o [output] -T [threads] -Q 1 -p -M  -g gene_id -a [GTF] [BAMs]
```
## FeatureCounts location

```bash
### Move to trim folder
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/bam 

### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/count.sh


```


## FeatureCounts read bam file
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/bam 
ls -1 *.bam 
ls -1 *.sortedByCoord.out.bam| tr '\n' ' '
```
## Edit templet
```bash
nano count.sh
```

```bash
#paste this to count.sh
featureCounts -o /data/gpfs/assoc/bch709-3/${USER}//mouse/readcount/featucount -T 4 -Q 1 -p -M  -g gene_id -a /data/gpfs/assoc/bch709-3/${USER}/mouse/ref/refGene.gtf $(ls -1 *.sortedByCoord.out.bam| tr '\n' ' ')
```
## Submit
```
sbatch count.sh
```

# Assignment
Please do samething with human data
https://www.sciencedirect.com/science/article/pii/S2211124720316235#fig1


Fastq location is `/data/gpfs/assoc/bch709-3/Course_materials/human/fastq`
Reference location is `/data/gpfs/assoc/bch709-3/Course_materials/human/ref`

```bash
/data/gpfs/assoc/bch709-3/Course_materials/human/
```







