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