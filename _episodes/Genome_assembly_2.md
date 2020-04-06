---
layout: page
title: 11_Genome assembly 2
published: true
---


### Genome assembly Spades
```bash
cd /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/Illumina/
mkdir Spades
cd Spades
conda create -n genomeassembly -y 
conda activate genomeassembly
conda install  -c r -c conda-forge -c anaconda -c bioconda  spades canu pacbio_falcon samtools minimap2 multiqc -y
conda install  -c r -c conda-forge -c anaconda -c bioconda  r-ggplot2 r-stringr r-scales r-argparse -y

```

```bash
#!/bin/bash
#SBATCH --job-name=Spades
#SBATCH --cpus-per-task=32
#SBATCH --time=2:00:00
#SBATCH --mem=64g
#SBATCH --account=cpu-s2-bch709-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o Spades.out # STDOUT
#SBATCH -e Spades.err # STDERR

spades.py -k 21,33,55,77 --careful -1 <trim_galore output> -2 <trim_galore output> -o spades_output --memory 64 --threads 32
```

### Spade
![spades]({{site.baseurl}}/fig/spades.jpg)
![spades2]({{site.baseurl}}/fig/spades2.jpg)
## Log
```bash
Command line: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades.py -k21,33,55,77     --careful       -1      /data/gpfs/assoc/bch709/wyim/gee/trimmed_fastq/WGS_R1_val_1.fq.gz    -2      /data/gpfs/assoc/bch709/wyim/gee/trimmed_fastq/WGS_R2_val_2.fq.gz    -o      /data/gpfs/assoc/bch709/wyim/gee/spades_output       --memory        120     --threads       32

System information:
  SPAdes version: 3.13.1
  Python version: 3.7.3
  OS: Linux-3.10.0-957.27.2.el7.x86_64-x86_64-with-centos-7.6.1810-Core

Output dir: /data/gpfs/assoc/bch709/spiderman/gee/spades_output
Mode: read error correction and assembling
Debug mode is turned OFF

Dataset parameters:
  Multi-cell mode (you should set '--sc' flag if input data was obtained with MDA (single-cell) technology or --meta flag if processing metagenomic dataset)
  Reads:
    Library number: 1, library type: paired-end
      orientation: fr
      left reads: ['/data/gpfs/assoc/bch709/wyim/gee/trimmed_fastq/WGS_R1_val_1.fq.gz']
      right reads: ['/data/gpfs/assoc/bch709/wyim/gee/trimmed_fastq/WGS_R2_val_2.fq.gz']
      interlaced reads: not specified
      single reads: not specified
      merged reads: not specified
Read error correction parameters:
  Iterations: 1
  PHRED offset will be auto-detected
  Corrected reads will be compressed
Assembly parameters:
  k: [21, 33, 55, 77]
  Repeat resolution is enabled
  Mismatch careful mode is turned ON
  MismatchCorrector will be used
  Coverage cutoff is turned OFF
Other parameters:
  Dir for temp files: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/tmp
  Threads: 64
  Memory limit (in Gb): 140


======= SPAdes pipeline started. Log can be found here: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/spades.log


===== Read error correction started.


== Running read error correction tool: /data/gpfs/home/wyim/miniconda3/envs/genomeassembly/bin/spades-hammer /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/configs/config.info

  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  75)   Starting BayesHammer, built from N/A, git revision N/A
  0:00:00.000     4M / 4M    INFO    General                 (main.cpp                  :  76)   Loading config from /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/configs/config.info
  0:00:00.001     4M / 4M    INFO    General                 (main.cpp                  :  78)   Maximum # of threads to use (adjusted due to OMP capabilities): 32
  0:00:00.001     4M / 4M    INFO    General                 (memory_limit.cpp          :  49)   Memory limit set to 140 Gb
  0:00:00.001     4M / 4M    INFO    General                 (main.cpp                  :  86)   Trying to determine PHRED offset
  0:00:00.038     4M / 4M    INFO    General                 (main.cpp                  :  92)   Determined value is 33
  0:00:00.038     4M / 4M    INFO    General                 (hammer_tools.cpp          :  36)   Hamming graph threshold tau=1, k=21, subkmer positions = [ 0 10 ]
  0:00:00.038     4M / 4M    INFO    General                 (main.cpp                  : 113)   Size of aux. kmer data 24 bytes
     === ITERATION 0 begins ===
  0:00:00.042     4M / 4M    INFO   K-mer Index Building     (kmer_index_builder.hpp    : 301)   Building kmer index
  0:00:00.042     4M / 4M    INFO    General                 (kmer_index_builder.hpp    : 117)   Splitting kmer instances into 512 files using 32 threads. This might take a while.
  0:00:00.043     4M / 4M    INFO    General                 (file_limit.hpp            :  32)   Open file limit set to 65536
  0:00:00.043     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  89)   Memory available for splitting buffers: 1.45829 Gb
  0:00:00.043     4M / 4M    INFO    General                 (kmer_splitters.hpp        :  97)   Using cell size of 131072
  0:00:02.300    17G / 17G   INFO   K-mer Splitting          (kmer_data.cpp             :  97)   Processing /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R1_val_1.fq.gz
  0:00:19.373    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             : 107)   Processed 3022711 reads
  0:00:19.373    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             :  97)   Processing /data/gpfs/assoc/bch709/spiderman/gee/trimmed_fastq/WGS_R2_val_2.fq.gz
  0:00:37.483    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             : 107)   Processed 6045422 reads
  0:00:37.483    17G / 18G   INFO   K-mer Splitting          (kmer_data.cpp             : 112)   Total 6045422 reads processed
  0:00:39.173   128M / 18G   INFO    General                 (kmer_index_builder.hpp    : 120)   Starting k-mer counting.
 
===== Mismatch correction finished.

 * Corrected reads are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/corrected/
 * Assembled contigs are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/contigs.fasta
 * Assembled scaffolds are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/scaffolds.fasta
 * Paths in the assembly graph corresponding to the contigs are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/scaffolds.paths
 * Assembly graph is in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/assembly_graph.fastg
 * Assembly graph in GFA format is in /data/gpfs/assoc/bch709/spiderman/gee/spades_output/assembly_graph_with_scaffolds.gfa

======= SPAdes pipeline finished.

SPAdes log can be found here: /data/gpfs/assoc/bch709/spiderman/gee/spades_output/spades.log

Thank you for using SPAdes!
```

### Assembly statistics

## N50 example
N50 is a measure to describe the quality of assembled genomes that are fragmented in contigs of different length. The N50 is defined as the minimum contig length needed to cover 50% of the genome.


|Contig Length|Cumulative Sum|  
|------|------|
| 1280 | 1280 |
| 1278 | 2558 |
| 1020 | 3578 |
| 990  | 4568 |
| 950  | 5518 |
| 852  | 6370 |
| 750  | 7120 |
| 400  | 7520 |
| 230  | 7750 |
| 200  | 7950 |
| 100  | 8050 | 


```bash
conda install -c bioconda assembly-stats
cd /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/Illumina/Spades
cd spades_output
assembly-stats scaffolds.fasta
assembly-stats contigs.fasta
```

![sccaffold]({{site.baseurl}}/fig/sccaffold.png)



### Assembly statistics result
```bash
stats for scaffolds.fasta
sum = 11241648, n = 2345, ave = 4793.88, largest = 677246
N50 = 167555, n = 17
N60 = 124056, n = 24
N70 = 97565, n = 34
N80 = 72177, n = 48
N90 = 38001, n = 68
N100 = 78, n = 2345
N_count = 308
Gaps = 4
```

```bash
stats for contigs.fasta
sum = 11241391, n = 2349, ave = 4785.61, largest = 666314
N50 = 167555, n = 17
N60 = 123051, n = 25
N70 = 95026, n = 36
N80 = 70916, n = 50
N90 = 36771, n = 71
N100 = 78, n = 2349
N_count = 0
Gaps = 0
```
```bash
cat scaffolds.paths
```
```bash
NODE_1_length_677246_cov_27.741934
5330914+,39250+,5099246-;
4754344-,4601428-,5180688-,1424894+,5327688+,732820-,5237058-,5052460-,4723018+,4800852+,5331930-,732820-,5019006-,5052460-,4755300-,4800852+,5331932-,5060558+,5185654-,5071338-,5178452+,5178460+,5178468+,5178476+,5254862-,5325448+,88806-,5243982-,5053698+,1425522-,5239940+,5238056-,4867204-,5331654+
NODE_1_length_677246_cov_27.741934'
5331654-,4867204+,5238056+,5239940-,1425522+,5053698-,5243982+,88806+,5325448-,5254862+,5178476-,5178468-,5178460-,5178452-,5071338+,5185654+,5060558-,5331932+,4800852-,4755300+,5052460+,5019006+,732820+,5331930+,4800852-,4723018-,5052460+,5237058+,732820+,5327688-,1424894-,5180688+,4601428+,4754344+;
5099246+,39250-,5330914-
NODE_2_length_666314_cov_27.644287
5327204+,103640+,4836832-,4851626+,5361530-,5361524-,5361528-,5329856-,1126012-,5329854-,1126012-,5236812+,5052228+,5236810+,5052228+,5099424-,4479150+,4812968+,4479150+,5062132+,414588+,5051858+,414588+,5331378+,4760742-,4925978+,5327370-,474724-,5094420-,4653402-,5331664-,4961012-,5018412+,5072166+,5040312+,5030606-,4961012-,5236106+,5072166+,5072168+,4915570-,5050466-,4765966-,5059218-,4915570-,5058386-,5236104+,5095538-,5095540-,5095538-,5149466+,4774822-,5017666+,4774822-,5065186-,876886-,4799048-,876886-,5332178-
NODE_2_length_666314_cov_27.644287'
5332178+,876886+,4799048+,876886+,5065186+,4774822+,5017666-,4774822+,5149466-,5095538+,5095540+,5095538+,5236104-,5058386+,4915570+,5059218+,4765966+,5050466+,4915570+,5072168-,5072166-,5236106-,4961012+,5030606+,5040312-,5072166-,5018412-,4961012+,5331664+,4653402+,5094420+,474724+,5327370+,4925978-,4760742+,5331378-,414588-,5051858-,414588-,5062132-,4479150-,4812968-,4479150-,5099424+,5052228-,5236810-,5052228-,5236812-,1126012+,5329854+,1126012+,5329856+,5361528+,5361524+,5361530+,4851626-,4836832+,103640-,5327204-
NODE_3_length_613985_cov_27.733595
5250014-,5121298+,5057128-,4953418-,5238246+,5264238+,5264242+,4468126+,5331520+,4813546-,4676908-,4813546-,5059540-,4862238+,5032536-,4862238+,5045932+,1122610+,4827200-,928516+,5031788-,4629584+,5007546-,1271448+,4907228+,1271448+,5099418-,5331326-,5030236-,5236282-,5100426-,5100418-,5100430-,139162-,4675324-,5354642-,372-,374+,5044194+,5058512+,5325918+,4544022-,4816684-,427838+,5238146+,269904+,117192-
NODE_3_length_613985_cov_27.733595'
117192+,269904-,5238146-,427838-,4816684+,4544022+,5325918-,5058512-,5044194-,374-,372+,5354642+,4675324+,139162+,5100430+,5100418+,5100426+,5236282+,5030236+,5331326+,5099418+,1271448-,4907228-,1271448-,5007546+,4629584-,5031788+,928516-,4827200+,1122610-,5045932-,4862238-,5032536+,4862238-,5059540+,4813546+,4676908+,4813546+,5331520-,4468126-,5264242-,5264238-,5238246-,4953418+,5057128+,5121298-,5250014+
NODE_4_length_431556_cov_27.699452
```
### FASTG file format

```bash
>EDGE_5360468_length_246_cov_13.568047:EDGE_5284398_length_327_cov_11.636000,EDGE_5354800_length_230_cov_14.470588';
GCTTCTTCTTGCTTCTCAAAGCTTTGTTGGTTTAGCCAAAGTCCAGATGAGTCTTTATCT
TTGTATCTTCTAACAAGGAAACACTACTTAGGCTTTTAGGATAAGCTTGCGGTTTAAGTT
TGTATACTCAATCATACACATGACATCAAGTCATATTCGACTCCAAAACACTAACCAAGC
TTCTTCTTGCACCTCAAAGCTTTGTTGGTTTAGCCAAAGTCCATATAAGTCTTTGTCTTT
GTATCT
>EDGE_5360470_length_161_cov_15.607143:EDGE_5332762_length_98_cov_43.619048';
GCTTCTTCTTGCTTCTCAAAGCTTTGTTGGTTTAGCCAAAGTCCAGATGAGTCTTTATCT
TTGTATCTTCTAACAAGAAAACACTACTTACGCTTTTAGGATAATGTTGCGGTTTAAGTT
CTTATACTCAATCATACACATGACATCAAGTCATATTCGAC
>EDGE_5354230_length_92_cov_267.066667':EDGE_5354222_length_86_cov_252.444444',EDGE_5355724_length_1189_cov_26.724820;
AAGCAAAGACTAAGTTTGGGGGAGTTGATAAGTGTGTATTTTGCATGTTTTGAGCATCCA
TTTGTCATCACTTTAGCATCATATCATCACTG
>EDGE_5344586_length_373_cov_22.574324:EDGE_5360654_length_82_cov_117.400000';
GCTAAAGTGATGACAAATGGATGCTCAAAACATGCAAAATACACACTTATCAACTCCCCC
AAACTTAGTCTTTGCTTAAGAACAAGCTGGAGGTGAGGTTTGAAAGCGGGGACTCAGAGC
CAAAGCAGCAGATAAACCAGATGAAATCAATGTCCAAGTTGATAGTTCTAAGTTGCGATA
TGATCGAATTCTACTCAAAAACGTTAGCCATGCCTTTTTATCAATCAATCCGACTCATAT
GCTCGACCTACACGTGTTTTCAAATCTACCAATCCCTTTAACATTCATTAGCTCTAGAAC
GTGAATCAAGCAATGCATCATCAATGAACTCATTTGGCTAAGGTAAAAGGTCAAGAGACA
AAGATGGTCCCTT
>EDGE_5354236_length_91_cov_242.857143:EDGE_5350728_length_80_cov_275.666667';
GCTAAAGTGATGACAAATGGATGCTCAAAACATGCAAAATACACACTTATCAACTCCCCC
AAACTTAGTCTTTGCTTGCCCTCAAGCAAAC

```
![bandage]({{site.baseurl}}/fig/bandage.png)
![assembly_spades]({{site.baseurl}}/fig/assembly_spades.png)

[bandage](https://rrwick.github.io/Bandage/)

## Assignment
Please upload Bandage output 



## PacBio assembly
```bash
conda deactivate
conda activate preprocessing
conda install -c bioconda nanostat nanoplot 
```
### Reads download
```bash
cd /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/
mkdir PacBio
cd PacBio
```

```bash
https://www.dropbox.com/s/7coua2gedbuykl6/BCH709_Pacbio_1.fastq.gz
https://www.dropbox.com/s/fniub0rxv48hupp/BCH709_Pacbio_2.fastq.gz
```

### Check PacBio reads statistics
```bash
NanoStat --fastq BCH709_Pacbio_1.fastq.gz
NanoPlot -t 2 --fastq  BCH709_Pacbio_1.fastq.gz --maxlength 40000 --plots hex dot pauvre -o pacbio_stat
```

### Transfer your result
```bash
*.png *.html *.txt
```

![HistogramReadlength]({{site.baseurl}}/fig/HistogramReadlength.png)
![LengthvsQualityScatterPlot_hex]({{site.baseurl}}/fig/LengthvsQualityScatterPlot_hex.png)



### PacBio reads statistics
```
General summary:        
Mean read length:               9,698.6
Mean read quality:                  6.6
Median read length:             8,854.0
Median read quality:                6.6
Number of reads:               58,497.0
Read length N50:               10,901.0
Total bases:              567,339,600.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:	58497 (100.0%) 567.3Mb
>Q7:	8260 (14.1%) 78.6Mb
>Q10:	0 (0.0%) 0.0Mb
>Q12:	0 (0.0%) 0.0Mb
>Q15:	0 (0.0%) 0.0Mb
Top 5 highest mean basecall quality scores and their read lengths
1:	8.5 (13544)
2:	8.5 (14158)
3:	8.5 (9590)
4:	8.5 (10741)
5:	8.5 (7529)
Top 5 longest reads and their mean basecall quality score
1:	24999 (6.4)
2:	24992 (6.0)
3:	24983 (7.0)
4:	24980 (6.5)
5:	24977 (6.6)
```




### Hybrid Genome assembly Spades
```bash
cd /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/

mkdir Spades_Illumina_Pacbio
cd Spades_Illumina_Pacbio
conda activate genomeassembly

```
### Submit below job

```bash
#!/bin/bash
#SBATCH --job-name=Spades
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --mem=140g
#SBATCH --account=cpu-s2-bch709-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o Spades.out # STDOUT
#SBATCH -e Spades.err # STDERR
zcat  <LOCATION_BCH709_Pacbio_1.fastq.gz> <LOCATION_BCH709_Pacbio_1.fastq.gz> >> merged_pacbio.fastq
spades.py -k 21,33,55,77 --careful -1 <trim_galore output> -2 <trim_galore output> --pacbio merged_pacbio.fastq -o spades_output --memory 140 --threads 64
```

## *De novo* assembly
![De-novo-assembly-High-level-diagram-of-long-read-assembly-pipeline-The-assembly-has](De-novo-assembly-High-level-diagram-of-long-read-assembly-pipeline-The-assembly-has.png)


## Canu
Canu (Koren et al. 2017) is a fork of the celera assembler and improves upon the earlier PBcR pipeline into a single, comprehensive assembler. Highly repetitive k-mers, which are abundant in all the reads, can be non-informative. Hence term frequency, inverse document frequency (tf-idf), a weighting statistic was added to MinHashing, giving weightage to non-repetitive k-mers as minimum values in the MinHash sketches, and sensitivity has been demonstrated to reach up to 89% without any parameter adjustment. By retrospectively inspecting the assembly graphs and also statistically filtering out repeat-induced overlaps, the chances of mis-assemblies are reduced.
![canu]({{site.baseurl}}/fig/canu.png)


## Canu assembly
```bash
cd /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/PacBio
conda activate genomeassembly
```
### Submit below job
```bash
#!/bin/bash
#SBATCH --job-name=Canu
#SBATCH --cpus-per-task=64
#SBATCH --time=2:00:00
#SBATCH --mem=140g
#SBATCH --account=cpu-s2-bch709-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o Canu.out # STDOUT
#SBATCH -e Canu.err # STDERR

canu -p canu -d canu_outdir genomeSize=11m corThreads=64 -pacbio-raw <LOCATION_BCH709_Pacbio_1.fastq.gz> <LOCATION_BCH709_Pacbio_1.fastq.gz>
```

## Check the Quality of Genome Assembly

```bash
cd /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/
mkdir genomeassembly_results/
```

### Activate Your environment

```bash
conda activate genomeassembly
```

### Copy Your Assembly Results
```bash
canu.contigs.fasta
spades_pacbio_illumina.fasta
spades_illumina.fasta
```

### Check Your Assembly Results

```bash
assembly-stats canu.contigs.fasta spades_pacbio_illumina.fasta spades_illumina.fasta
```

## Compare Assemblies
![dotplot2]({{site.baseurl}}/fig/dotplot2.png)


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
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOURID>@unr.edu
#SBATCH -o nucmer.out # STDOUT
#SBATCH -e nucmer.err # STDERR
#SBATCH --account=cpu-s6-test-0 
#SBATCH --partition=cpu-s6-core-0
nucmer  --coords -p canu_pacbio_Spades_illumina <canu.contigs> <illumina_spades_scaffold_file>
nucmer  --coords -p canu_pacbio_Spades_illumina_pacbio <canu.contigs> <Spades_illumina_pacbio_scaffold_file>
```

### Dot 
Dot is an interactive dot plot viewer for genome-genome alignments.

Dot is publicly available here: https://dnanexus.github.io/dot And can also be used locally by cloning this repository and simply opening the index.html file in a web browser.


After aligning genome assemblies or finished genomes against each other with MUMmer's nucmer, the alignments can be visualized with Dot. Instead of generating static dot plot images on the command-line, Dot lets you interact with the alignments by zooming in and investigating regions in detail.

To prepare a .delta file (nucmer output) for Dot, run this python (3.6) script first: https://dnanexus.github.io/dot/DotPrep.py

The DotPrep.py script will apply a unique anchor filtering algorithm to mark alignments as unique or repetitive. This algorithm analyzes all of the alignments, and it needs to see unfiltered data to determine which alignments are repetitive, so make sure to run nucmer without any filtering options and without running delta-filter on the .delta file before passing it into DotPrep.py.


```bash
wget https://dnanexus.github.io/dot/DotPrep.py
```

```bash
nano DotPrep.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=dot
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o nucmer.out # STDOUT
#SBATCH -e nucmer.err # STDERR
#SBATCH --account=cpu-s6-test-0 
#SBATCH --partition=cpu-s6-core-0
python DotPrep.py  --delta canu_pacbio_Spades_illumina.delta --out  canu_pacbio_Spades_illumina
python DotPrep.py  --delta canu_pacbio_Spades_illumina_pacbio.delta  --out canu_pacbio_Spades_illumina_pacbio
```

The output of DotPrep.py includes the \*.coords and \*.coords.idx that should be used with Dot for visualization.


## Visualization
- Transfer \*.coords.\* files
- Go to  https://dnanexus.github.io/dot/

![dotplot4]({{site.baseurl}}/fig/dotplot4.png)


The output of DotPrep.py includes the \*.coords and \*.coords.idx that should be used with Dot for visualization.


## Which one is the best?

![alignment_reference]({{site.baseurl}}/fig/alignment_reference.png)
![structure]({{site.baseurl}}/fig/structure.png)


## Assignment
```bash
wget https://www.dropbox.com/s/xpnhj6j99dr1kum/Athaliana_subset_BCH709.fa
```



## Pilon
```bash
conda create -n postprocess python=3 -y
conda activate postprocess
conda install -y -c bioconda pilon bwa samtools

mkdir pilon ## at your genome assembly folder
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
#SBATCH --mail-user=<YOUREMAIL>
#SBATCH -o mapping.out # STDOUT
#SBATCH -e mapping.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

bwa index canu.contigs.fasta
bwa mem canu.contigs.fasta  <Trimmed_illumina_R1 val something> <Trimmed_illumina_R2 val something>  -o canu_illumina.sam
samtools view -Sb canu_illumina.sam -o canu_illumina.bam
samtools sort canu_illumina.bam -o canu_illumina_sort.bam
samtools index canu_illumina_sort.bam

```

```bash
nano /data/gpfs/home/<YOURID>/miniconda3/envs/postprocess/bin/pilon 


default_jvm_mem_opts = ['-Xms8g', '-Xmx40g']
```
```bash
nano pilon.sh
```

```bash 
#!/bin/bash
#SBATCH --job-name=pilon
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUREMAIL>
#SBATCH -o pilon.out # STDOUT
#SBATCH -e pilon.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

pilon --genome canu.contigs.fasta --frags canu_illumina_sort.bam --output canu.illumina  --vcf --changes
```


```bash
sbatch --dependency=afterok:<123456> pilon.sh
```

### How can we improve these genome assemblies?

![illumina]({{site.baseurl}}/fig/mate.png)
![optical mapping]({{site.baseurl}}/fig/bionano.jpg)
![pacbio_scaff]({{site.baseurl}}/fig/pacbio_scaff.png)




## BUSCO
BUSCO assessments are implemented in open-source software, with a large selection of lineage-specific sets of Benchmarking Universal Single-Copy Orthologs. These conserved orthologs are ideal candidates for large-scale phylogenomics studies, and the annotated BUSCO gene models built during genome assessments provide a comprehensive gene predictor training set for use as part of genome annotation pipelines.

https://busco.ezlab.org/v2/


```bash
conda create -n busco python=3
conda activate busco
conda install busco multiqc
```

```bash
https://busco.ezlab.org/v2/datasets/embryophyta_odb9.tar.gz
```

```bash 
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o busco.out # STDOUT
#SBATCH -e busco.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

export AUGUSTUS_CONFIG_PATH="~/miniconda3/envs/busco/config/"

run_busco -i <canu.contigs.fasta> --cpu 24  -o canu  -l embryophyta_odb9 -m geno -s arabidopsis
run_busco -i <pacbio_illumina_spades.fasta> --cpu 24  -o pacbio_illumina_spades  -l embryophyta_odb9 -m geno -s arabidopsis
multiqc . -n assembly
```

## BUSCO results
```
 run_busco -i canu.contigs.fasta  --cpu 16  -o canu  -l embryophyta_odb9 -m geno -r
INFO    ****************** Start a BUSCO 3.0.2 analysis, current time: 11/21/2019 00:14:45 ******************
INFO    Configuration loaded from /data/gpfs/home/wyim/miniconda3/envs/busco/bin/../config/config.ini
INFO    Init tools...
INFO    Check dependencies...
INFO    Check input file...
INFO    To reproduce this run: python /data/gpfs/home/wyim/miniconda3/envs/busco/bin/run_busco -i canu.contigs.fasta -o canu -l embryophyta_odb9/ -m genome -c 16 -sp arabidopsis
INFO    Mode is: genome
INFO    The lineage dataset is: embryophyta_odb9 (eukaryota)
INFO    Temp directory is ./tmp/
WARNING Restarting an uncompleted run
INFO    ****** Phase 1 of 2, initial predictions ******
INFO    Phase 1 was already completed.
INFO    Results:
INFO    C:9.2%[S:8.7%,D:0.5%],F:0.8%,M:90.0%,n:1440
INFO    132 Complete BUSCOs (C)
INFO    125 Complete and single-copy BUSCOs (S)
INFO    7 Complete and duplicated BUSCOs (D)
INFO    11 Fragmented BUSCOs (F)
INFO    1297 Missing BUSCOs (M)
INFO    1440 Total BUSCO groups searched

```




## Chromosome assembly
```bash
cd  /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/
mkdir /data/gpfs/assoc/bch709/<YOURID>/Genome_assembly/hic   ## will go to genome assembly folder adn make hic folder
cd !$
```



### HiC

```bash
conda create -n hic

conda activate hic

conda install samtools bedtools matplotlib numpy scipy bwa
```


### Download
```bash
https://www.dropbox.com/s/0waw9b2uy4iarq2/hic_r1.fastq.gz
https://www.dropbox.com/s/tq0iy4815hw473z/hic_r2.fastq.gz
https://www.dropbox.com/s/2vku066402h5una/allhic.zip
```
### Decompress
```
unzip allhic.zip
chmod -R 775 allhic 
```
## canu.illumina.fasta
**The input of HiC is the output of Pilon. If you don't have it please do Pilon first.**


### Run AllHIC ->  hic.sh
```bash
#!/bin/bash
#SBATCH --job-name=AllHIC
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mem=48g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUREMAIL>
#SBATCH -o hic.out # STDOUT
#SBATCH -e hic.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0
unzip allhic.zip
cp /data/gpfs/assoc/bch709/<YOURID>/<YOUR GENOMEASSEMBLY FOLDER>/pilon/canu.illumina.fasta*  .
samtools faidx canu.illumina.fasta
bwa index canu.illumina.fasta
bwa mem -t 24 -SPM canu.illumina.fasta hic_r1.fastq.gz hic_r2.fastq.gz  > hic.sam
samtools view -Sb hic.sam -o hic.bam -@ 24
chmod 775 ALL* all*
./allhic extract hic.bam canu.illumina.fasta
./allhic partition hic.counts_GATC.txt hic.pairs.txt 2
./allhic optimize hic.counts_GATC.2g1.txt  hic.clm
./allhic optimize hic.counts_GATC.2g2.txt  hic.clm
./allhic  build hic.counts_GATC.2g1.tour hic.counts_GATC.2g2.tour canu.illumina.fasta bch709_assembly
cut -f 1,2 canu.illumina.fasta.fai >> chrn.list
ALLHiC_plot  hic.bam bch709_assembly.agp chrn.list 10k pdf
```


```bash
sbatch --dependency=afterok:<123456> hic.sh
```



### Investigate taxa
```bash
conda install kraken kraken2 
```
We will be using a tool called Kraken2. This tool uses k-mers to assign a taxonomic labels in form of NCBI Taxonomy to the sequence (if possible). The taxonomic label is assigned based on similar k-mer content of the sequence in question to the k-mer content of reference genome sequence. The result is a classification of the sequence in question to the most likely taxonomic label. If the k-mer content is not similar to any genomic sequence in the database used, it will not assign any taxonomic label.

We can also use another tool by the same group called Centrifuge. This tool uses a novel indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index, optimized specifically for the metagenomic classification problem to assign a taxonomic labels in form of NCBI Taxonomy to the sequence (if possible). The result is a classification of the sequence in question to the most likely taxonomic label. If the search sequence is not similar to any genomic sequence in the database used, it will not assign any taxonomic label.





###########################

### Assignment by 11/26
## Chromosome assembly
```bash
mkdir /data/gpfs/assoc/bch709/<YOURID>/<g TAB>/hic  
cd !$
```
### HiC

```bash
conda create -n hic

conda activate hic

conda install samtools bedtools matplotlib numpy scipy bwa
```


### Download
```bash
https://www.dropbox.com/s/0waw9b2uy4iarq2/hic_r1.fastq.gz
https://www.dropbox.com/s/tq0iy4815hw473z/hic_r2.fastq.gz
https://www.dropbox.com/s/2vku066402h5una/allhic.zip
```
### Decompress
```
unzip allhic.zip
```

### Run AllHIC ->  hic.sh
```bash
#!/bin/bash
#SBATCH --job-name=AllHIC
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mem=48g
#SBATCH --mail-type=all
#SBATCH --mail-user=<YOUREMAIL>
#SBATCH -o hic.out # STDOUT
#SBATCH -e hic.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0
unzip allhic.zip
cp /data/gpfs/assoc/bch709/<YOURID>/<YOUR GENOMEASSEMBLY FOLDER>/genomeassembly_results/canu.illumina.fasta*  .
samtools faidx canu.illumina.fasta
bwa index canu.illumina.fasta
bwa mem -t 24 -SPM canu.illumina.fasta hic_r1.fastq.gz hic_r2.fastq.gz  > hic.sam
samtools view -Sb hic.sam -o hic.bam -@ 24
chmod 775 ALL* all*
./allhic extract hic.bam canu.illumina.fasta
./allhic partition hic.counts_GATC.txt hic.pairs.txt 2
./allhic optimize hic.counts_GATC.2g1.txt  hic.clm
./allhic optimize hic.counts_GATC.2g2.txt  hic.clm
./allhic  build hic.counts_GATC.2g1.tour hic.counts_GATC.2g2.tour canu.illumina.fasta bch709_assembly
cut -f 1,2 canu.illumina.fasta.fai >> chrn.list
ALLHiC_plot  hic.bam groups.agp chrn.list 10k pdf
```


```bash
sbatch --dependency=afterok:<123456> hic.sh
```

![hic1]({{site.baseurl}}/fig/hic1.png)
![hic1]({{site.baseurl}}/fig/hic2.png)
![hic1]({{site.baseurl}}/fig/hic3.png)
![hic1]({{site.baseurl}}/fig/hic4.png)
![hic1]({{site.baseurl}}/fig/hic5.png)
![hic1]({{site.baseurl}}/fig/hic6.png)
![hic1]({{site.baseurl}}/fig/hic7.png)
![hic1]({{site.baseurl}}/fig/hic8.png)
![hic1]({{site.baseurl}}/fig/hic9.png)
![hic1]({{site.baseurl}}/fig/hic10.png)

[!][(http://img.youtube.com/vi/-MxEw3IXUWU/0.jpg)](http://www.youtube.com/watch?v=-MxEw3IXUWU " ")


![![hic1]({{site.baseurl}}/fig/hic10.png)]({{site.baseurl}}/fig/hicmovie.gif)





