---
layout: page
title: 	Final_review
published: true
---


## Conda environment for RNA-Seq
```bash
conda create -n RNASEQ_bch709 -c bioconda -c conda-forge  -c r  sra-tools minimap2 trinity star trim-galore gffread seqkit kraken2 samtools multiqc subread
conda activate RNASEQ_bch709
```
## Conda Environment for DEG
```bash
conda create -n DEG_bch709 -y

conda activate DEG_bch709
conda install -y -c bioconda -c conda-forge mamba
mamba install -y -c bioconda -c conda-forge r-gplots r-fastcluster=1.1.25  bioconductor-ctc  bioconductor-deseq2 bioconductor-qvalue  bioconductor-limma bioconductor-edger bioconductor-genomeinfodb bioconductor-deseq2 r-rcurl trinity bedtools intervene r-UpSetR r-corrplot r-Cairo
```

## Publication (Arabidopsis)
> 
>A Vitis vinifera basic helix–loop–helix transcription factor enhances plant cell size, vegetative biomass and reproductive yield Sung Don Lim,Won Choel Yim,Degao Liu,Rongbin Hu,Xiaohan Yang,John C. Cushman
>https://doi.org/10.1111/pbi.12898
> 
{: .callout}

## DEG analaysis
The question will provide 12 RNA-Seq reads files associated with four different conditions.
You might need to compare empty vector vs. CEB1 transformation line in leaf and root samples.
The reads and reference file will be provided.

1. Trim the reads by Trim-Galore
2. Align the reads by STAR
3. Count reads per gene by FeatureCount2
4. Quality control by PtR
5. DEG calculation by DESeq2

## MultiQC report
Generate MultiQC report

## Draw Venn diagram
The question will ask you to draw 4-way Venn diagram from DEG analysis.

## Gene expression 
The question will ask you to provide TPM value for one gene.

## Gene Ontology analysis
The question will ask you to use Metascape http://metascape.org/gp/index.html for DEG set.

## Seqkit and BLAST
The question will ask you to find protein sequence from one of DEG gene and ask you to run BLAST analaysis.

## Submitted batch file
Slurm submission files need to be uploaded.



### Trim-galore
```bash
cd  ~/bch709_scratch/RNA-Seq_example/ATH
nano trim.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=trim_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761506 raw_data/SRR1761506_1.fastq.gz raw_data/SRR1761506_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761507 raw_data/SRR1761507_1.fastq.gz raw_data/SRR1761507_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761508 raw_data/SRR1761508_1.fastq.gz raw_data/SRR1761508_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761509 raw_data/SRR1761509_1.fastq.gz raw_data/SRR1761509_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761510 raw_data/SRR1761510_1.fastq.gz raw_data/SRR1761510_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761511 raw_data/SRR1761511_1.fastq.gz raw_data/SRR1761511_2.fastq.gz --fastqc
```

### Create reference index
```bash
cd  ~/bch709_scratch/RNA-Seq_example/ATH/reference
ls -algh
nano index.sh
```
```bash
#!/bin/bash
#SBATCH --job-name=index_ATH
#SBATCH --cpus-per-task=12
#SBATCH --time=2-15:00:00
#SBATCH --mem=48g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o index.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles   phytozome/phyto_mirror/Athaliana_167_10/assembly/Athaliana_167.fa  --sjdbGTFfile TAIR10_GFF3_genes.gtf --sjdbOverhang 99   --genomeSAindexNbases 12
```

## Mapping the reads to genome index
```bash
cd  ~/bch709_scratch/RNA-Seq_example/ATH/
ls -algh
nano align.sh
```
```bash
#!/bin/bash
#SBATCH --job-name=align_ATH
#SBATCH --cpus-per-task=8
#SBATCH --time=2-15:00:00
#SBATCH --mem=32g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o align.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0
#SBATCH --dependency=afterok:<PREVIOUS_JOBID(trim_ATH)>

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761506.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761507_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761507_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761507.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761508_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761508_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761508.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761509_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761509_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761509.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761510_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761510_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761510.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761511_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761511_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761511.bam
```

## Featurecount
```
featureCounts -p  -a <GENOME>.gtf <SAMPLE1>.bam <SAMPLE2>.bam <SAMPLE3>.bam  ...... -o counts.txt
```

```bash
conda activate RNASEQ_bch709
cd ~/bch709_scratch/RNA-Seq_example/ATH/bam
featureCounts -o ATH.featureCount.cnt -p  -a ~/bch709_scratch/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf SRR1761506.bamAligned.sortedByCoord.out.bam  SRR1761509.bamAligned.sortedByCoord.out.bam SRR1761507.bamAligned.sortedByCoord.out.bam  SRR1761510.bamAligned.sortedByCoord.out.bam SRR1761508.bamAligned.sortedByCoord.out.bam  SRR1761511.bamAligned.sortedByCoord.out.bam
```

```bash
conda activate RNASEQ_bch709
cd ~/bch709_scratch/RNA-Seq_example/Mmusculus/bam
featureCounts -o Mmusculus.featureCount.cnt -p  -a ~/bch709_scratch/RNA-Seq_example/Mmusculus/reference/GCF_000001635.27_GRCm39_genomic.gtf -g "gene_name"  <YOUR BAM FILES>
```

### TPM and FPKM calculation

```bash
cut -f1,6-  ATH.featureCount.cnt |  egrep -v "#" | sed 's/\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g' > ATH.featureCount_count_length.cnt

python /data/gpfs/assoc/bch709-2/Course_material/script/tpm_raw_exp_calculator.py -count ATH.featureCount_count_length.cnt

```

## ATH DEG
```bash

cd ~/bch709_scratch/RNA-Seq_example/ATH
mkdir DEG
cd DEG
cp ~/bch709_scratch/RNA-Seq_example/ATH/bam/ATH.featureCount* .

cut -f1,7- ATH.featureCount.cnt | egrep -v "#" | sed 's/\.bamAligned\.sortedByCoord\.out\.bam//g; s/\.TAIR10//g' > ATH.featureCount_count_only.cnt 
```


### PtR (Quality Check Your Samples and Biological Replicates)

Once you've performed transcript quantification for each of your biological replicates, it's good to examine the data to ensure that your biological replicates are well correlated, and also to investigate relationships among your samples. If there are any obvious discrepancies among your sample and replicate relationships such as due to accidental mis-labeling of sample replicates, or strong outliers or batch effects, you'll want to identify them before proceeding to subsequent data analyses (such as differential expression). 
```bash
PtR  --matrix ATH.featureCount_count_only.cnt  --samples samples.txt --CPM  --log2 --min_rowSums 10   --sample_cor_matrix --compare_replicates

```

### DEG calculation
```bash
run_DE_analysis.pl --matrix Drosophila.featureCount_count_only.cnt  --method DESeq2 --samples_file samples.txt --output rnaseq
```


### DEG subset
```bash
cd rnaseq
## 4-fold and p-value 0.01
analyze_diff_expr.pl --samples ~/bch709_scratch/RNA-Seq_example/ATH/DEG/samples.txt  --matrix ~/bch709_scratch/RNA-Seq_example/ATH/DEG/ATH.featureCount_count_length.cnt.tpm.tab -P 0.01 -C 2 --output ATH

## 2-fold and p-value 0.01
analyze_diff_expr.pl --samples  ~/bch709_scratch/RNA-Seq_example/ATH/DEG/samples.txt   --matrix ~/bch709_scratch/RNA-Seq_example/ATH/DEG/ATH.featureCount_count_length.cnt.tpm.tab -P 0.01 -C 1 --output ATH
```


### DEG output
```
ATH.matrix.log2.centered.sample_cor_matrix.pdf
ATH.matrix.log2.centered.genes_vs_samples_heatmap.pdf

ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C2.ABA-UP.subset
ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C2.Control-UP.subset
ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C2.DE.subset

ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C1.ABA-UP.subset
ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C1.Control-UP.subset
ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C1.DE.subset
```

## Venn diagram

### Intervene installation
```bash
mamba install -c bioconda bedtools intervene r-UpSetR=1.4.0 r-corrplot r-Cairo
```

```bash
cd ~/bch709_scratch/RNA-Seq_example/ATH/DEG/rnaseq
cut -f 1 ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C2.ABA-UP.subset |  grep -v sample > DESeq.UP_4fold.subset
cut -f 1 ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C2.Control-UP.subset  |  grep -v sample > DESeq.DOWN_4fold.subset 

cut -f 1 ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C1.ABA-UP.subset |  grep -v sample > DESeq.UP_2fold.subset
cut -f 1 ATH.featureCount_count_only.cnt.ABA_vs_Control.DESeq2.DE_results.P0.01_C1.Control-UP.subset  |  grep -v sample >DESeq.DOWN_2fold.subset
```

```bash
 wc -l DESeq*subset
```
```
  701 DESeq.DOWN_2fold.subset
  227 DESeq.DOWN_4fold.subset
 1218 DESeq.UP_2fold.subset
  463 DESeq.UP_4fold.subset
 2609 total
```

```bash
intervene venn --type list --save-overlaps -i DESeq.DOWN_2fold.subset DESeq.DOWN_4fold.subset DESeq.UP_2fold.subset DESeq.UP_4fold.subset
intervene upset --type list --save-overlaps -i DESeq.DOWN_2fold.subset DESeq.DOWN_4fold.subset DESeq.UP_2fold.subset DESeq.UP_4fold.subset 
```

### Seqkit
```bash
seqkit grep -p {ID}  {Protein sequence}
seqkit grep -p AT4G28110.1  Athaliana_167_TAIR10.cds_primaryTranscriptOnly.fa   -o AT4G28110.1.aa
```


### BLAST
>PR1_CDS
MNFTGYSRFLIVFVALVGALVLPSKAQDSPQDYLRVHNQARGAVGVGPMQWDERVAAYARSYAEQLRGNCRLIHSGGPYGENLAWGSGDLSGVSAVNMWVSEKANYNYAANTCNGVCGHYTQVVWRKSVRLGCAKVRCNNGGTIISCNYDPRGNYVNEKPY

>PR1_CDS
ATGAATTTTACTGGCTATTCTCGATTTTTAATCGTCTTTGTAGCTCTTGTAGGTGCTCTTGTTCTTCCCTCGAAAGCTCAAGATAGCCCACAAGATTATCTAAGGGTTCACAACCAGGCACGAGGAGCGGTAGGCGTAGGTCCCATGCAGTGGGACGAGAGGGTTGCAGCCTATGCTCGGAGCTACGCAGAACAACTAAGAGGCAACTGCAGACTCATACACTCTGGTGGGCCTTACGGGGAAAACTTAGCCTGGGGTAGCGGTGACTTGTCTGGCGTCTCCGCCGTGAACATGTGGGTTAGCGAGAAGGCTAACTACAACTACGCTGCGAACACGTGCAATGGAGTTTGTGGTCACTACACTCAAGTTGTTTGGAGAAAGTCAGTGAGACTCGGATGTGCCAAAGTGAGGTGTAACAATGGTGGAACCATAATCAGTTGCAACTATGATCCTCGTGGGAATTATGTGAACGAGAAGCCATACTAA