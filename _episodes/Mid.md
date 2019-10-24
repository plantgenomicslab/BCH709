---
layout: page
title: Midterm_exam
published: True
---

>## Please Read This First
>
> Please make folder `/data/gpfs/assoc/bch709/<YOUR_ID>/Midterm`
> Process everything under `/data/gpfs/assoc/bch709/<YOUR_ID>/Midterm`
> [Link](https://forms.gle/XxbuCXKzzhFsbPMi7)
{: .callout}




Question 1: Describe how the process of RNA sequencing enables the determination of the transcriptome for an entire genome?  and discuss one example of bioinformatics problems where transcriptome assembly provides a good approach. (max. 200 words)

```
AAT  
AGT   
ATT  
ATA  
TAG  
TAT  
TAC  
TGA  
TTA  
```

Question 2: Given the sequence shown above, if it takes 3 kmer to assemble a pair of reads by de bruijn algorithm, please provide assembled sequence.

Question 3: Which sequencing platform is the best for RNA-Seq? Justify your answer.

Question 4: Describe the fastq and fasta format.

Question 5: Given the following file permission attributes from `ls -algh`, please describe the meaning
```
-rwxr-xr-x
```
Question 6: Please create "midterm" Conda environment, activate and install the following software. After installation, please provide the list of packages in the environment by 'conda list'

```
python=3 trim-galore jellyfish salmon cmake htop rsem nano conda-build trinity r-fastcluster fastqc multiqc bioconductor-edger biopython
```

Question 7: Please download following files, provide the number of reads (WT: Wildtype, TR: Treatment, R1: Left; forward, R2: Right; reverse) (Ex: wt_r2:185,  wt_r2:5023 etc)

```
https://www.dropbox.com/s/p0nxhe78wicybv0/wt_r2.fastq.gz  
https://www.dropbox.com/s/da30aue6f4erlq5/wt_r1.fastq.gz  
https://www.dropbox.com/s/64c2ybx5fpyaqs4/tr_r2.fastq.gz  
https://www.dropbox.com/s/y23c9jjfdup5e79/tr_r1.fastq.gz  
```
Question 8: Please run `multiqc` after processing `fastqc`, check GC contents for each reads.


Question 9: Please merge all fastq.gz files by `zcat`  and run Trinity with the following options. Provide the number of assembled sequences.

```
#!/bin/bash  
#SBATCH --job-name=<YOUR_NAME>  
#SBATCH --cpus-per-task=64  
#SBATCH --time=2:00:00  
#SBATCH --mem=120g  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=<YOUR_EMAIL>  
#SBATCH -o Trinity.out # STDOUT  
#SBATCH -e Trinity.err # STDERR  
  
  
Trinity --seqType fq  --CPU 64 --max_memory 100G --left <MERGED_R1_FILE> --right <MERGED_R2_FILE> --no_normalize_reads  
```


Question 10: Provide the statistics of assembly by `TrinityStats.pl`


Question 11: Please run alignment with the following options and provide the most abundance transcripts ID in WT and TR
```
#!/bin/bash  
#SBATCH --job-name=<YOUR_NAME>  
#SBATCH --cpus-per-task=64  
#SBATCH --time=2:00:00  
#SBATCH --mem=120g  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=<YOUR_EMAIL>  
#SBATCH -o Trinity.out # STDOUT  
#SBATCH -e Trinity.err # STDERR  
  
align_and_estimate_abundance.pl --transcripts <YOUR_Trinity_OUTPUT> --seqType fq --left <wt_r1.fastq._LOCATION> --right <wt_r2.fastq._LOCATION> --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir WT_rsem_outdir  --thread_count  64

align_and_estimate_abundance.pl --transcripts <YOUR_Trinity_OUTPUT> --seqType fq --left <tr_r1.fastq._LOCATION> --right <tr_r2.fastq._LOCATION> --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir TR_rsem_outdir  --thread_count  64

```

Question 12: Please calculate FPKM value of `TRINITY_DN0_c0_g1` before RSEM process from `WT_rsem_outdir/RSEM.isoforms.results` and `TR_rsem_outdir/RSEM.isoforms.results`


Question 13: Please normalize expression value by TMM and provide the expression matrix (`RSEM.isoform.TMM.EXPR.matrix`).
```
#!/bin/bash  
#SBATCH --job-name=<YOUR_NAME>  
#SBATCH --cpus-per-task=64  
#SBATCH --time=2:00:00  
#SBATCH --mem=120g  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=<YOUR_EMAIL>   
#SBATCH -o Trinity.out # STDOUT  
#SBATCH -e Trinity.err # STDERR  

abundance_estimates_to_matrix.pl  --est_method RSEM --gene_trans_map none --name_sample_by_basedir  --cross_sample_norm TMM  <WT RSEM.isoforms.results OUTPUT>  <TR RSEM.isoforms.results OUTPUT> 

```
Question 14: Please normalize expression value by TMM and provide the expression matrix (`RSEM.isoform.TMM.EXPR.matrix`).



