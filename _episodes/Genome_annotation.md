---
layout: page
title: Genome annoation
published: true
---

# Genome Annotation

After the sections of DNA sequence have been assembled into a complete genome sequence we need to identify where the genes and key features are. We have our aligned and assembled genome sequence but how do we identify where the genes and other functional regions of the genome are located.

- Annotation involves marking where the genes start and stop in the DNA sequence and also where other relevant and interesting regions are in the sequence.

- Although genome annotation pipelines can differ from one another, for example, some elements can be manual while others have to be automated, they all share a core set of features.

- They are generally divided into two distinct phases: gene prediction and manual annotation.

## Gene prediction
There are two types of gene prediction: 
Ab initio – this technique relies on signals within the DNA sequence. It is an automated process whereby a computer is given instructions for finding genes in the sequence and is then left to find them. The computer looks for common sequences known to be found at the start and end of genes such as promoter sequences (where proteins bind that switch on genes), start codons (where the code for the gene product, RNA ?or protein, starts) and stop codons (where the code for the gene product ends).  

## Evidence-based
This technique relies on evidence beyond the DNA sequence. It involves gathering various pieces of genetic information from the transcript sequence (mRNA), and known protein sequences of the genome. With these pieces of evidence it is then possible to get an idea of the original DNA sequence by working backwards through transcription? and translation? (reverse transcription/translation). For example, if you have the protein sequence it is possible to work out the family of possible DNA sequences it could be derived from by working out which amino acids? make up the protein and then which combination of codons could code for those amino acids and so on, until you get to the DNA sequence. 
The information taken from these two prediction methods is then combined and lined up with the sequenced genome.

## Key Point
- Gene annotation is one of the core mechanisms through which we decipher the information that is contained in genome sequences.

- Gene annotation is complicated by the existence of 'transcriptional complexity', which includes extensive alternative splicing and transcriptional events outside of protein-coding genes.

- The annotation strategy for a given genome will depend on what it is hoped to achieve, as well as the resources available.

- The availability of next-generation data sets has transformed gene annotation pipelines in recent years, although their incorporation is rarely straightforward.

- Even human gene annotation is far from complete: transcripts are missing and existing models are truncated. Most importantly, 'functional annotation' — the description of what transcripts actually do — remains far from comprehensive.

- Efforts are now under way to integrate gene annotation pipelines with projects that seek to describe regulatory sequences, such as promoter and enhancer elements.

- Gene annotation is producing increasingly complex resources. This can present a challenge to usability, most notably in a clinical context, and annotation projects must find ways to resolve such problems.

### The core annotation workflows for different gene types.
![nrg.2016.119-f2]({{site.baseurl}}/fig/nrg.2016.119-f2.jpg)

These workflows illustrate general annotation principles rather than the specific pipelines of any particular genebuild. a | Protein-coding genes within reference genomes were generally annotated on the basis of the computational genomic alignment of Sanger-sequenced transcripts and protein-coding sequences, followed by manual annotation using interface tools such as Zmap, WebApollo, Artemis and the Integrative Genomics Viewer. Transcripts were typically taken from GenBank129 and proteins from Swiss-Prot. b | Protein-coding genes within non-reference genomes are usually annotated based on fewer resources; in this case, RNA sequencing (RNA-seq) data are used in combination with protein homology information that has been extrapolated from a closely related genome. RNA-seq pipelines for read alignment include STAR and TopHat, whereas model creation is commonly carried out by Cufflinks23. c | Long non-coding RNA (lncRNA) structures can be annotated in a similar manner to protein-coding transcripts (parts a and b), although coding potential must be ruled out. This is typically done by examining sequence conservation with PhyloCSF or using experimental data sets, such as mass spectrometry or ribosome profiling. In this example, 5′ Cap Analysis of Gene Expression (CAGE) and polyA-seq data are also incorporated to obtain true transcript end points. Designated lncRNA pipelines include PLAR. d | Small RNAs are typically added to genebuilds by mining repositories such as RFAM or miRBase. However, these entries can be used to search for additional loci based on homology. e | Pseudogene annotation is based on the identification of loci with protein homology to either paralogous or orthologous protein-coding genes. Computational annotation pipelines include PseudoPipe53, although manual annotation is more accurate. Finally, all annotation methods can be thwarted by the existence of sequence gaps in the genome assembly (right-angled arrow). EST, expressed sequence tag.

### High-level strategies for gene annotation projects.
![nrg.2016.119-f3]({{site.baseurl}}/fig/nrg.2016.119-f3.jpg)

This schematic details the annotation pathways for reference and novel genomes. Coding sequences (CDSs) are outlined in green, nonsense-mediated decay (NMD) is shown in purple and untranslated regions (UTRs) are filled in in red. The core evidence sets that are used at each stage are listed, although their availability and incorporation can vary across different projects. The types of evidence used for reference genebuilds have evolved over time: RNA sequencing (RNA-seq) has replaced Sanger sequencing, conservation-based methodologies have become more powerful and proteogenomic data sets are now available. By contrast, novel genebuilds are constructed based on RNA-seq and/or ab initio modelling, in combination with the projection of annotation from other species (which is known as liftover) and the use of other species evidence sets. In fact, certain novel genebuilds such as those of pigs and rats now incorporate a modest amount of manual annotation, and could perhaps be described as 'intermediate' in status between 'novel' and 'reference'. Furthermore, such genebuilds have also been improved by community annotation; this process typically follows the manual annotation workflows for reference genomes, although on a smaller scale. Although all reference genebuilds are 'mature' in our view, progress into the 'extended genebuild' phase is most advanced for humans. A promoter is indicated by the blue circle, an enhancer is indicated by the orange circle, and binding sites for transcription factors (TFs) or RNA-binding proteins (RBPs) are shown as orange triangles. Gene expression can be analysed on any genebuild regardless of quality, although it is more effective when applied to accurate transcript catalogues. Clearly, the results of expression analyses have the potential to reciprocally improve the efficacy of genebuilds, although it remains to be seen how this will be achieved in practice (indicated by the question mark). 3D, three-dimensional; EST, expressed sequence tag.


## General Considerations

### Bacteria use ATG as their main start codon, but GTG and TTG are also fairly common, and a few others are occasionally used.
– Remember that start codons are also used internally: the actual start codon may not be
the first one in the ORF.

Please check [Genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

### The stop codons are the same as in eukaryotes: TGA, TAA, TAG
– stop codons are (almost) absolute: except for a few cases of programmed frameshifts and the use of TGA for selenocysteine, the stop codon at the end of an ORF is the end of protein translation.

### Genes can overlap by a small amount. Not much, but a few codons of overlap is common enough so that you can’t just eliminate overlaps as impossible.
- Cross-species homology works well for many genes. It is very unlikely that noncoding
sequence will be conserved.
– But, a significant minority of genes (say 20%) are unique to a given species.

### Translation start signals (ribosome binding sites; Shine-Dalgarno sequences) are often found just upstream from the start codon
– however, some aren’t recognizable
– genes in operons sometimes don’t always have a separate ribosome binding site for each gene


## Based on hidden Markov model (HMM)
- As you move along the DNA sequence, a given nucleotide can be in an exon or an
intron or in an intergenic region.
- The oversimplified model on this slide doesn’t have the ”non-gene” state
- Use a training set of known genes (from the same or closely related species) to determine transmission and emission probabilities.

Very simple HMM: each base is either in an intron or an exon, and gets emitted with different
frequencies depending on which state it is in.

## A simple HMM for modeling eukaryotic genes.
![CG-10-402_F1]({{site.baseurl}}/fig/CG-10-402_F1.jpg)

The posterior probability P{yn = i | x, Θ} can be computed from

![Annotationhmm]({{site.baseurl}}/fig/Annotationhmm.png)

![HMM2]({{site.baseurl}}/fig/hmm2.png)



## Genome Annotation
![maker]({{site.baseurl}}/fig/maker.png)


```bash
mkdir /data/gpfs/assoc/bch709/<YOURID>/tmp

mkdir /data/gpfs/assoc/bch709/<YOURID>/genomeannotation

cd !$
```

### Conda environment
```bash
conda create -n genomeannotation busco maker kraken2 bracken krona repeatmasker snap rmblast
conda activate genomeannotation
```

### Setting augustus
```bash
export AUGUSTUS_CONFIG_PATH="~/miniconda3/envs/genomeannotation/config/"
```

### Setting RepeatMaker
```bash
cd ~/miniconda3/envs/genomeannotation/share/RepeatMasker
./configure
```

### RMBLAST
```bash
/data/gpfs/home/wyim/miniconda3/envs/genomeannotation/bin/
```



### Copy your draft genome
Copy `bch709_assembly.fasta` (from HiC) to current folder (genomeannotation).

### MAKER software
MAKER is a portable and easily configurable genome annotation pipeline. Its purpose is to allow smaller eukaryotic and prokaryotic genome projects to independently annotate their genomes and to create genome databases. MAKER identifies repeats, aligns ESTs and proteins to a genome, produces ab-initio gene predictions and automatically synthesizes these data into gene annotations having evidence-based quality values. MAKER is also easily trainable: outputs of preliminary runs can be used to automatically retrain its gene prediction algorithm, producing higher quality gene-models on seusequent runs. MAKER's inputs are minimal and its ouputs can be directly loaded into a GMOD database. They can also be viewed in the Apollo genome browser; this feature of MAKER provides an easy means to annotate, view and edit individual contigs and BACs without the overhead of a database. MAKER should prove especially useful for emerging model organism projects with minimal bioinformatics expertise and computer resources.

### Test MAKER
```bash
maker --help
```

### Download Uniprot database
https://www.uniprot.org/

```bash
wget -O uniprot.fasta.gz "https://www.uniprot.org/uniprot/?query=arabidopsis&format=fasta&force=true&sort=score&fil=rev
iewed:yes&compress=yes

gunzip 
```

### Generate Control files
```bash
maker -CTL
```

### Update Control files
```
est=/data/gpfs/assoc/bch709/spiderman/rnaseq/assembly_quality/transrate_results/Trinity/good.Trinity.fasta
rmlib=/data/gpfs/assoc/pgl/bin/maker/data/te_proteins.fasta
augustus_species=arabidopsis
est2genome=1
protein2genome=1
TMP=/data/gpfs/assoc/bch709/spiderman/tmp 
```

### Run Maker
```bash
nano maker.sh
```

```bash 
#!/bin/bash
#SBATCH --job-name=maker
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o maker.out # STDOUT
#SBATCH -e maker.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

maker -cpus 24 -base bch709  -genome bch709_assembly.fasta
```


