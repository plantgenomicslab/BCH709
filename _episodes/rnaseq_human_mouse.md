---
layout: page
title: RNA-Seq
published: true
---
# Environment creation
```bash
conda create -n BCH709_RNASeq -c bioconda -c conda-forge -c anaconda python=3.7 mamba 

conda activate BCH709_RNASeq

mamba install -c bioconda -c conda-forge -c anaconda trim-galore=0.6.7 sra-tools=2.11.0 STAR htseq=1.99.2 subread=2.0.1 multiqc=1.11 snakemake=7.5.0 parallel-fastq-dump=0.6.7 bioconductor-tximport samtools=1.14 r-ggplot2 trinity=2.13.2 hisat2 bioconductor-qvalue sambamba graphviz gffread tpmcalculator lxml rsem
```
### Slurm

![](https://i.imgur.com/LynACgh.png)


![](https://i.imgur.com/XEMRbJe.png)


[CheatSheet](https://slurm.schedmd.com/pdfs/summary.pdf)


*Slurm provides resource management for the processors allocated to a job, so that multiple job steps can be simultaneously submitted and queued until there are available resources within the job's allocation.*



# Mouse RNA-Seq
![](https://i.imgur.com/rqNc6OA.png)

https://www.sciencedirect.com/science/article/pii/S2211124722011111


### Working directory (Pronghorn)

```bash
echo $USER

cd /data/gpfs/assoc/bch709-3/${USER}/mouse
```


### Create file list
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/fastq

ls -1 *.gz 

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g'

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u > /data/gpfs/assoc/bch709-3/${USER}/mouse/filelist

cat /data/gpfs/assoc/bch709-3/${USER}/mouse/filelist
```



## RNA-Seq Alignment
```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/trim

#### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/mapping.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Trim/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/mapping.sh

#### Edit templet
nano mapping.sh
```

### Check output
```bash
ls -algh /data/gpfs/assoc/bch709-3/${USER}/mouse/trim
```

### Output example
```bash
[FILENAME]_R1_val_1.fq.gz [FILENAME]_R2_val_2.fq.gz
```

### STAR RNA-Seq read alignment example
```bash
STAR --runMode alignReads --runThreadN [CPU_NUMBER] --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --quantMode TranscriptomeSAM GeneCounts --genomeDir [GENOME_DIR] --readFilesCommand gunzip -c --readFilesIn [FORWARD_READ]  [REVERSE_READ] --outSAMtype BAM SortedByCoordinate --outFileNamePrefix [BAMFILENAME_LOCATION]
```

### STAR RNA-Seq alignment
```bash
STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --quantMode TranscriptomeSAM GeneCounts --genomeDir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/[FILENAME]_R1_val_1.fq.gz  /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/[FILENAME]_R2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/[FILENAME].bam
```


### STAR RNA-Seq alignment batch file test
```bash
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read1} /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/${i}.bam"
    done
```

### STAR RNA-Seq alignment batch file 

```bash
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-3/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read1} /data/gpfs/assoc/bch709-3/${USER}/mouse/trim/${read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/${i}.bam" | cat mapping.sh - > ${i}_mapping.sh
    done
```

###  Job submission dependency on Mapping
```bash
for i in `ls -1 *_mapping.sh`
do
    sbatch $i 
done

```

## FeatureCounts
```bash
featureCounts -o [output] -T [threads] -Q 1 -p -M  -g gene_id -a [GTF] [BAMs]
```

### FeatureCounts execute location

```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/bam 

ls -1 *.bam

#### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/count.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Count/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-3/${USER}/mouse/bam/count.sh

```


### FeatureCounts read bam file
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/mouse/bam 
ls -1 *.sortedByCoord.out.bam
ls -1 *.sortedByCoord.out.bam| tr '\n' ' '
```
### Edit templet
```bash
nano count.sh
```

```bash
#paste this to count.sh
featureCounts -o /data/gpfs/assoc/bch709-3/${USER}//mouse/readcount/featucount -T 4 -Q 1 -p -M  -g gene_id -a /data/gpfs/assoc/bch709-3/${USER}/mouse/ref/refGene.gtf $(for i in `cat /data/gpfs/assoc/bch709-3/wyim/mouse/filelist`; do echo ${i}.bamAligned.sortedByCoord.out.bam| tr '\n' ' ';done)
```


### Job submission dependency

```bash
squeue --noheader --format %i --user ${USER} 
```

### Submit
```bash
jobid=$(squeue --noheader --format %i --user ${USER} | tr '\n'  ':')1

sbatch --dependency=afterany:${jobid} count.sh 
```


# Human RNA-Seq
[***Transcriptome alterations in myotonic dystrophy frontal cortex***](https://www.sciencedirect.com/science/article/pii/S2211124720316235)

![](https://i.imgur.com/BrugOCz.png)

[](https://doi.org/10.1016/j.celrep.2020.108634)


## Environment activation
```bash
conda activate BCH709_RNASeq
```

## Working directory (Pronghorn)

```bash
echo $USER

cd /data/gpfs/assoc/bch709-3/${USER}
```


## Reference Download

https://www.ncbi.nlm.nih.gov/genome/guide/human/


```bash
### change working directory
cd /data/gpfs/assoc/bch709-3/${USER}/human/ref 


### download
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz


### decompress
gunzip GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.gff.gz



gffread  GRCh38_latest_genomic.gff  --keep-exon-attrs -F -T -o GRCh38_latest_genomic.gtf
```

## STAR reference build
### STAR aligner reference build on Pronghorn
```bash
### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/human/run.sh /data/gpfs/assoc/bch709-3/${USER}/human/ref/ref_build.sh

#open text editor
### PLEASE RENAME EMAIL AND JOB NAME
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/ref_build/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-3/${USER}/human/ref/ref_build.sh
nano ref_build.sh

# Add below command to ref_build.sh

STAR  --runThreadN 4 --runMode genomeGenerate --genomeDir . --genomeFastaFiles GRCh38_latest_genomic.fna --sjdbGTFfile GRCh38_latest_genomic.gtf  --sjdbOverhang 99   --genomeSAindexNbases 12
```


### Submit job to HPC
```bash
#submit job
sbach ref_build.sh

#check job
squeue -u ${USER}
```


### FASTQ file 
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/human/

### Link file (without copy)
ln -s /data/gpfs/assoc/bch709-3/Course_materials/human/fastq/* /data/gpfs/assoc/bch709-3/${USER}/human/fastq

ls /data/gpfs/assoc/bch709-3/${USER}/human/fastq
```

### Create file list
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/human/fastq

ls -1 *.gz 

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g'

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u > /data/gpfs/assoc/bch709-3/${USER}/human/filelist

cat /data/gpfs/assoc/bch709-3/${USER}/human/filelist
```
### Regular expression

https://regex101.com/

## Trim reads

### Prepare templet
```bash
cp /data/gpfs/assoc/bch709-3/Course_materials/human/run.sh /data/gpfs/assoc/bch709-3/${USER}/human/fastq/trim.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Trim/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-3/${USER}/human/fastq/trim.sh

```
### Edit templet
```bash
nano /data/gpfs/assoc/bch709-3/${USER}/human/fastq/trim.sh
```

### Batch submission file
```bash
# Check file list
cat ../filelist
nano trim.sh

# Loop file list
### Add Forward read to variable
### Add reverse read from forward read name substitution
### add file name from variable to trim-galore
### merge trim-galore command and trim.sh
### add trim-galore command and trim.sh to new file
for i in `cat ../filelist`
    do
        read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        
        echo "trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-3/${USER}/human/trim $read1 $read2" | cat trim.sh - > ${i}_trim.sh
        echo "$read1 $read2 trim file has been created."
done
```
### Batch submission
```bash
ls -1 *_trim.sh
### Loop *.sh submission
for i in `ls -1 *_trim.sh`
do
    sbatch $i
done
```

### Check submission
```bash
squeue -u ${USER}
```

## RNA-Seq Alignment templet file
```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-3/${USER}/human/trim

#### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/human/run.sh /data/gpfs/assoc/bch709-3/${USER}/human/trim/mapping.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Mapping/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-3/${USER}/human/trim/mapping.sh


#### Edit templet
nano mapping.sh
```


### Check output
```bash
ls -algh /data/gpfs/assoc/bch709-3/${USER}/human/trim
```
### Output example
```bash
[FILENAME]_R1_val_1.fq.gz [FILENAME]_R2_val_2.fq.gz
```


### STAR RNA-Seq alignment
```bash
STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --quantMode TranscriptomeSAM GeneCounts --genomeDir /data/gpfs/assoc/bch709-3/${USER}/human/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/human/trim/[FILENAME]_R1_val_1.fq.gz  /data/gpfs/assoc/bch709-3/${USER}/human/trim/[FILENAME]_R2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/human/bam/[FILENAME].bam
```
### STAR RNA-Seq alignment batch file 

```bash
cd /data/gpfs/assoc/bch709-3/${USER}/human/trim
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-3/${USER}/human/ref --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-3/${USER}/human/trim/${read1} /data/gpfs/assoc/bch709-3/${USER}/human/trim/${read2} --outFileNamePrefix /data/gpfs/assoc/bch709-3/${USER}/human/bam/${i}.bam" | cat mapping.sh - > ${i}_mapping.sh
    done
```


### Job submission dependency

```bash
squeue --noheader --format %i --user ${USER}  
squeue --noheader --format %i --user ${USER} | tr '\n'  ':'
```


###  Job submission dependency on Trim
```bash

jobid=$(squeue --noheader --format %i --user ${USER} | tr '\n'  ':')1
for i in `ls -1 *_mapping.sh`
do
    sbatch  --dependency=afterany:${jobid} $i 
done

```

### FeatureCounts
[Bioinformatics, Volume 30, Issue 7, 1 April 2014, Pages 923–930](https://doi.org/10.1093/bioinformatics/btt656)
![]({{{site.baseurl}}/fig/featurecount.png)

```bash
featureCounts -o [output] -T [threads] -Q 1 -p -M  -g gene_id -a [GTF] [BAMs]
```
### FeatureCounts location

```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-3/${USER}/human/bam 

#### Copy templet
cp /data/gpfs/assoc/bch709-3/Course_materials/human/run.sh /data/gpfs/assoc/bch709-3/${USER}/human/bam/count.sh

sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Count/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g"  /data/gpfs/assoc/bch709-3/${USER}/human/bam/count.sh
```


### FeatureCounts command to count.sh
#### LOOP example
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/human/bam 

ls -1 *.bam 

for i in `cat /data/gpfs/assoc/bch709-3/wyim/human/filelist`
do 
echo ${i}.bamAligned.sortedByCoord.out.bam | tr '\n' ' '
done
```



### FeatureCount 
```bash

echo "featureCounts -o /data/gpfs/assoc/bch709-3/${USER}//mouse/readcount/featucount -T 4 -Q 1 -p -M  -g gene_id -a /data/gpfs/assoc/bch709-3/${USER}/human/ref/GRCh38_latest_genomic.gtf $(for i in `cat /data/gpfs/assoc/bch709-3/wyim/human/filelist`; do echo ${i}.bamAligned.sortedByCoord.out.bam| tr '\n' ' ';done)" >> count.sh
```


### Job submission dependency

```bash
squeue --noheader --format %i --user ${USER} 
squeue --noheader --format %i --user ${USER} | tr '\n'  ':'
```


###  Job submission dependency on Align
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/human/bam 

jobid=$(squeue --noheader --format %i --user ${USER} | tr '\n'  ':')1
sbatch  --dependency=afterany:${jobid} count.sh 
```



## Differential expression
DESeq2
edgeR (Neg-binom > GLM > Test)
Limma-Voom (Neg-binom > Voom-transform > LM > Test)



### DESeq
DESeq: This normalization method is included in the DESeq Bioconductor package (version 1.6.0) and is based on the hypothesis that most genes are not DE. A DESeq scaling factor for a given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. The underlying idea is that non-DE genes should have similar read counts across samples, leading to a ratio of 1. **Assuming most genes are not DE, the median of this ratio for the lane provides an estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis.** By calling the estimateSizeFactors() and sizeFactors() functions in the DESeq Bioconductor package, this factor is computed for each lane, and raw read counts are divided by the factor associated with their sequencing lane.  
[DESeq2](https://www.ncbi.nlm.nih.gov/pubmed/22988256)

### EdgeR
Trimmed Mean of M-values (TMM): This normalization method is implemented in the edgeR Bioconductor package (version 2.4.0). It is also based on the hypothesis that most genes are not DE. The TMM factor is computed for each lane, with one lane being considered as a reference sample and the others as test samples. For each test sample, TMM is computed as the weighted mean of log ratios between this test and the reference, after exclusion of the most expressed genes and the genes with the largest log ratios. **According to the hypothesis of low DE, this TMM should be close to 1. If it is not, its value provides an estimate of the correction factor that must be applied to the library sizes (and not the raw counts) in order to fulfill the hypothesis.** The calcNormFactors() function in the edgeR Bioconductor package provides these scaling factors. To obtain normalized read counts, these normalization factors are re-scaled by the mean of the normalized library sizes. Normalized read counts are obtained by dividing raw read counts by these re-scaled normalization factors.  
[EdgeR](https://www.ncbi.nlm.nih.gov/pubmed/22988256)

## DESeq2 vs EdgeR Statistical tests for differential expression
### DESeq2
DESeq2 uses raw counts, rather than normalized count data, and models the normalization to fit the counts within a Generalized Linear Model (GLM) of the negative binomial family with a logarithmic link. Statistical tests are then performed to assess differential expression, if any.  

### EdgeR
Data are normalized to account for sample size differences and variance among samples. The normalized count data are used to estimate per-gene fold changes and to perform statistical tests of whether each gene is likely to be differentially expressed.  
EdgeR uses an exact test under a negative binomial distribution (Robinson and Smyth, 2008). The statistical test is related to Fisher's exact test, though Fisher uses a different distribution.  

## DESeq2 vs EdgeR Normalization method
DESeq and EdgeR are very similar and both assume that no genes are differentially expressed. DEseq uses a "geometric" normalisation strategy, whereas EdgeR is a weighted mean of log ratios-based method. Both normalise data initially via the calculation of size / normalisation factors.

## Negative binormal
### DESeq2 
ϕ was assumed to be a function of μ (population mean) determined by nonparametric regression. The recent version used in this paper follows a more versatile procedure. Firstly, for each transcript, an estimate of the dispersion is made, presumably using maximum likelihood. Secondly, the estimated dispersions for all transcripts are fitted to the functional form:  
s
### EdgeR
edgeR recommends a “tagwise dispersion” function, which estimates the dispersion on a gene-by-gene basis, and implements an empirical Bayes strategy for squeezing the estimated dispersions towards the common dispersion. Under the default setting, the degree of squeezing is adjusted to suit the number of biological replicates within each condition: more biological replicates will need to borrow less information from the complete set of transcripts and require less squeezing.  

***[DEG software comparison paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-484)***

*** DEseq uses a "geometric" normalisation strategy, whereas EdgeR is a weighted mean of log ratios-based method. Both normalise data initially via the calculation of size / normalisation factors ***


## Functional analysis • GO
Gene enrichment analysis (Hypergeometric test)
Gene set enrichment analysis (GSEA)
Gene ontology / Reactome databases


## FPKM
Fragments per Kilobase of transcript per million mapped reads

![FPKM]({{site.baseurl}}/fig/FPKM.png)

X = mapped reads count
N = number of reads
L = Length of transcripts

### Featurecount output (Read count)
```bash
head -n 2 /data/gpfs/assoc/bch709-3/Course_materials/mouse/featurecount.txt
```

### Sum of FPKM for 12WK_R6-2_Rep_1.bamAligned.sortedByCoord.out.bam
```bash
cat /data/gpfs/assoc/bch709-3/Course_materials/mouse/featurecount.txt | egrep -v Geneid | awk '{ sum+=$7} END {print sum}'
```

### Call Python
```bash
python
```

```python
X = 553

total_umber_Reads_mapped = 47414569

Length = 3634

fpkm = X*(1000/Length)*(1000000/total_umber_Reads_mapped)

fpkm
```

#### ten to the ninth power = 10\*\*9

```python
fpkm=X/(Number_Reads_mapped*Length)*10**9
fpkm
```


### TPM
 Transcripts Per Million

![TPM]({{site.baseurl}}/fig/TPM.png)

![TPM2]({{site.baseurl}}/fig/TPM2.png)


### TPM calculation from reads count
```bash
cat /data/gpfs/assoc/bch709-3/Course_materials/mouse/featurecount.txt | egrep -v Geneid | awk '{ if($6 >= 0) sum+=$7/$6} END {print sum}'
```

```python

sum_count_per_length=17352.8
X = 553
Length = 3634

TPM = (X/Length)*(1/sum_count_per_length )*10**6

```

### Paper read
Fu, Yu, et al. "Elimination of PCR duplicates in RNA-seq and small RNA-seq using unique molecular identifiers." BMC genomics 19.1 (2018): 531
Parekh, Swati, et al. "The impact of amplification on differential expression analyses by RNA-seq." Scientific reports 6 (2016): 25533
Klepikova, Anna V., et al. "Effect of method of deduplication on estimation of differential gene expression using RNA-seq." PeerJ 5 (2017): e3091

