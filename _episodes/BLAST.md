---
layout: page
title: 14_BLAST
published: true
---


## location
```
 mkdir /data/gpfs/assoc/bch709-3/<YOURID>/BLAST
 cd $!
```

## ENV
```bash
 conda create -n blast -c bioconda  -c conda-forge blast seqkit -y
```

## BLAST
Basic Local Alignment Search Tool (Altschul et al., 1990 & 1997) is a sequence comparison algorithm optimized for speed used to search sequence databases for optimal local alignments to a query. The initial search is done for a word of length "W" that scores at least "T" when compared to the query using a substitution matrix. Word hits are then extended in either direction in an attempt to generate an alignment with a score exceeding the threshold of "S". The "T" parameter dictates the speed and sensitivity of the search.

### Rapidly compare a sequence Q to a database to find all sequences in the database with an score above some cutoff S.
- Which protein is most similar to a newly sequenced one?
- Where does this sequence of DNA originate?
- Speed achieved by using a procedure that typically finds *most* matches with scores > S.
- Tradeoff between sensitivity and specificity/speed
- Sensitivity – ability to find all related sequences
- Specificity – ability to reject unrelated sequences 


### Homologous sequence are likely to contain a short high scoring word pair, a seed.
– Unlike Baeza-Yates, BLAST *doesn't* make explicit guarantees 

### BLAST then tries to extend high scoring word pairs to compute maximal high scoring segment pairs (HSPs).
– Heuristic algorithm but evaluates the result statistically.

![seed]({{site.baseurl}}/fig/seed.png)
![seed]({{site.baseurl}}/fig/I.7_1_blast_illustration.png)


### E-value
E-value = the number of HSPs having score S (or higher) expected to occur by chance.

Smaller E-value, more significant in statistics
Bigger E-value , by chance

**E[# occurrences of a string of length m in reference of length L] ~ L/4m**



## PAM and BLOSUM Matrices
Two different kinds of amino acid scoring matrices, *PAM (Percent Accepted Mutation) and BLOSUM (BLOcks SUbstitution Matrix)*, are in wide use. The PAM matrices were created by Margaret Dayhoff and coworkers and are thus sometimes referred to as the Dayhoff matrices. These scoring matrices have a strong theoretical component and make a few evolutionary assumptions. The BLOSUM matrices, on the other hand, are more empirical and derive from a larger data set. Most researchers today prefer to use BLOSUM matrices because in silico experiments indicate that searches employing BLOSUM matrices have higher sensitivity.

There are several PAM matrices, each one with a numeric suffix. The PAM1 matrix was constructed with a set of proteins that were all 85 percent or more identical to one another. The other matrices in the PAM set were then constructed by multiplying the PAM1 matrix by itself: 100 times for the PAM100; 160 times for the PAM160; and so on, in an *attempt to model the course of sequence evolution*. Though highly theoretical (and somewhat suspect), it is certainly a reasonable approach. There was little protein sequence data in the 1970s when these matrices were created, so this approach was a good way to extrapolate to larger distances.

Protein databases contained many more sequences by the 1990s so a more empirical approach was possible. The BLOSUM matrices were constructed by extracting ungapped segments, or blocks, from a set of multiply aligned protein families, and then further clustering these blocks on the basis of their percent identity. The blocks used to derive the BLOSUM62 matrix, for example, all have at least 62 percent identity to some other member of the block.

![PAM-250-and-Blosum-62-matrices]({{site.baseurl}}/fig/PAM-250-and-Blosum-62-matrices.png)

![codon]({{site.baseurl}}/fig/codon.jpg)

### BLAST has a number of possible programs to run depending on whether you have nucleotide or protein sequences:

nucleotide query and nucleotide db - blastn
nucleotide query and nucleotide db - tblastx (includes six frame translation of query and db sequences)
nucleotide query and protein db - blastx (includes six frame translation of query sequences)
protein query and nucleotide db - tblastn (includes six frame translation of db sequences)
protein query and protein db - blastp

![blasttype]({{site.baseurl}}/fig/blasttype.png)
### BLAST Process


![step1]({{site.baseurl}}/fig/step1.png)
![step2]({{site.baseurl}}/fig/step2.png)
![step3]({{site.baseurl}}/fig/step3.png)
![step4]({{site.baseurl}}/fig/step4.png)


![blast]({{site.baseurl}}/fig/blast.gif)



### NCBI BLAST
https://blast.ncbi.nlm.nih.gov/Blast.cgi

### Uniprot

https://www.uniprot.org/


### BLASTN example
Run blastn against the nt database.


```

ATGAAAGCGAAGGTTAGCCGTGGTGGCGGTTTTCGCGGTGCGCTGAACTA
CGTTTTTGACGTTGGCAAGGAAGCCACGCACACGAAAAACGCGGAGCGAG
TCGGCGGCAACATGGCCGGGAATGACCCCCGCGAACTGTCGCGGGAGTTC
TCAGCCGTGCGCCAGTTGCGCCCGGACATCGGCAAGCCCGTCTGGCATTG
CTCGCTGTCACTGCCTCCCGGCGAGCGCCTGAGCGCCGAGAAGTGGGAAG
CCGTCGCGGCTGACTTCATGCAGCGCATGGGCTTTGACCAGACCAATACG
CCGTGGGTGGCCGTGCGCCACCAGGACACGGACAAGGATCACATCCACAT
CGTGGCCAGCCGGGTAGGGCTGGACGGGAAAGTGTGGCTGGGCCAGTGGG
AAGCCCGCCGCGCCATCGAGGCGACCCAAGAGCTTGAGCATACCCACGGC
CTGACCCTGACGCCGGGGCTGGGCGATGCGCGGGCCGAGCGCCGGAAGCT
GACCGACAAGGAGATCAACATGGCCGTGAGAACGGGCGATGAACCGCCGC
GCCAGCGTCTGCAACGGCTGCTGGATGAGGCGGTGAAGGACAAGCCGACC
GCGCTAGAACTGGCCGAGCGGCTACAGGCCGCAGGCGTAGGCGTCCGGGC
AAACCTCGCCAGCACCGGGCGCATGAACGGCTTTTCCTTCGAGGTGGCCG
GAGTGCCGTTCAAAGGCAGCGACTTGGGCAAGGGCTACACATGGGCGGGG
CTACAGAAAGCAGGGGTGACTTATGACGAAGCTAGAGACCGTGCGGGCCT
TGAACGATTCAGGCCCACAGTTGCAGATCGTGGAGAGCGTCAGGACGTTG
CAGCAGTCCGTGAGCCTGATGCACGAGGACTTGAAGCGCCTACCGGGCGC
AGTCTCGACCGAGACGGCGCAGACCTTGGAACCGCTGGCCCGACTCCGGC
AGGACGTGACGCAGGTTCTGGAAGCCTACGACAAGGTGACGGCCATTCAG
CGCAAGACGCTGGACGAGCTGACGCAGCAGATGAGCGCGAGCGCGGCGCA
GGCCTTCGAGCAGAAGGCCGGGAAGCTGGACGCGACCATCTCCGACCTGT
CGCGCAGCCTGTCAGGGCTGAAAACGAGCCTCAGCAGCATGGAGCAGACC
GCGCAGCAGGTGGCGACCTTGCCGGGCAAGCTGGCGAGCGCACAGCAGGG
CATGACGAAAGCCGCCGACCAACTGACCGAGGCAGCGAACGAGACGCGCC
CGCGCCTTTGGCGGCAGGCGCTGGGGCTGATTCTGGCCGGGGCCGTGGGC
GCGATGCTGGTAGCGACTGGGCAAGTCGCTTTAAACAGGCTAGTGCCGCC
AAGCGACGTGCAGCAGACGGCAGACTGGGCCAACGCGATTTGGAACAAGG
CCACGCCCACGGAGCGCGAGTTGCTGAAACAGATCGCCAATCGGCCCGCG
AACTAGACCCGACCGCCTACCTTGAGGCCAGCGGCTACACCGTGAAGCGA
GAAGGGCGGCACCTGTCCGTCAGGGCGGGCGGTGATGAGGCGTACCGCGT
GACCCGGCAGCAGGACGGGCGCTGGCTCTGGTGCGACCGCTACGGCAACG
ACGGCGGGGACAATATCGACCTGGTGCGCGAGATCGAACCCGGCACCGGC
TACGCCGAGGCCGTCTATCGGCTTTCAGGTGCGCCGACAGTCCGGCAGCA
ACCGCGCCCGAGCGAGCCGAAGCGCCAACCGCCGCAGCTACCGGCGCAAG
GGCTGGCAGCCCGCGAGCATGGCCGCGACTACCTCAAGGGCCGGGGCATC
AGCCAGGACACCATCGAGCACGCCGAGAAGGCGGGCATGGTGCGCTATGC
AGACGGTGGAGTGCTGTTCGTCGGCTACGACCGTGCAGGCACCGCGCAGA
ACGCCACACGCCGCGCCATTGCCCCCGCTGACCCGGTGCAGAAGCGCGAC
CTACGCGGCAGCGACAAGAGCTATCCGCCGATCCTGCCGGGCGACCCGGC
AAAGGTCTGGATCGTGGAAGGTGGCCCGGATGCGCTGGCCCTGCACGACA
TCGCCAAGCGCAGCGGCCAGCAGCCGCCCACCGTCATCGTGTCAGGCGGG
GCGAACGTGCGCAGCTTCTTGGAGCGGGCCGACGTGCAAGCGATCCTGAA
GCGGGCCGAGCGCGTCACCGTGGCCGGGGAAAACGAGAAGAACCCCGAGG
CGCAGGCAAAGGCCGACGCCGGGCACCAGAAGCAGGCGCAGCGGGTGGCC
AAAATCACCGGGCGCGAGGTGCGCCAATGGACGCCGAAGCCCGAGCACGG
CAAGGACTTGGCCGACATGAACGCCCGGCAGGTGGCAGAGATCGAGCGCA
AGCGACAGGCCGAGATCGAGGCCGAAAGAGCACGAAACCGCGAGCTTTCA
CGCAAGAGCCGGAGGTATGATGGCCCCAGCTTCGGCAGATAA
```

### BLASTP Query
Do a BLASTP on NCBI website with the following protein against nr, but limit the organism to cetartiodactyla using default parameters:

```
MASGPGGWLGPAFALRLLLAAVLQPVSAFRAEFSSESCRELGFSSNLLCSSCDLLGQFSL
LQLDPDCRGCCQEEAQFETKKYVRGSDPVLKLLDDNGNIAEELSILKWNTDSVEEFLSEK
LERI
```


Have a look at the multiple sequence alignment, can you explain the results?  

Do a similar blastp vs UniProtKB (UniProt) without post filtering.


## Running a standalone BLAST program
### location
```
cd /data/gpfs/assoc/bch709-3/<YOURID>/BLAST
```

### ENV
```bash
conda activate blast
```

### Running a standalone BLAST program
Create the index for the target database using makeblastdb;
Choose the task program: blastn, blastp, blastx, tblatx, psiblast or deltablast;
Set the configuration for match, mismatch, gap-open penalty, gap-extension penalty or scoring matrix;
Set the word size;
Set the E-value threshold;
Set the output format and the number of output results

### Standalone BLAST 
In addition to providing BLAST sequence alignment services on the web, NCBI also makes these sequence alignment utilities available for download through FTP. This allows BLAST searches to be performed on local platforms against databases downloaded from NCBI or created locally. These utilities run through DOS-like command windows and accept input through text-based command line switches. There is no graphic user interface

https://www.ncbi.nlm.nih.gov/books/NBK52640/

ftp://ftp.ncbi.nlm.nih.gov/blast/db/

### NR vs NT

At NCBI they are two different things as well. 'nr' is a database of protein sequences and 'nt' is nucleotide. At one time 'nr' meant non-redundant but it stopped being non-redundant a while ago. nt is a nucleotide database, while nr is a protein database (in amino acids)



### Standalone BLAST
1. Download the database.
2. Use makeblastdb to build the index.
3. Change the scoring matrix, record the changes in the alignment results and interpret the results.

### Download Database
```
wget ftp://ftp.ncbi.nih.gov/refseq/release/plant/plant.1.protein.faa.gz
```
### How many sequences in `plant.1.protein.faa.gz`


### Input file
```bash
/data/gpfs/assoc/bch709-3/Course_materials/BLAST/Athaliana_167_TAIR10.cds.fa.gz
```


### Example Input sequence

```bash
seqkit stats /data/gpfs/assoc/bch709-3/Course_materials/BLAST/Athaliana_167_TAIR10.cds.fa.gz -T
file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
/data/gpfs/assoc/bch709-3/Course_materials/BLAST/Athaliana_167_TAIR10.cds.fa.gz	FASTA	DNA	35386	43546761	22	1230.6	16182
```


### Subsampling by SeqKit

FASTA and FASTQ are basic and ubiquitous formats for storing nucleotide and protein sequences. Common manipulations of FASTA/Q file include converting, searching, filtering, deduplication, splitting, shuffling, and sampling. Existing tools only implement some of these manipulations, and not particularly efficiently, and some are only available for certain operating systems. Furthermore, the complicated installation process of required packages and running environments can render these programs less user friendly.

This project describes a cross-platform ultrafast comprehensive toolkit for FASTA/Q processing. SeqKit provides executable binary files for all major operating systems, including Windows, Linux, and Mac OS X, and can be directly used without any dependencies or pre-configurations. SeqKit demonstrates competitive performance in execution time and memory usage compared to similar tools. The efficiency and usability of SeqKit enable researchers to rapidly accomplish common FASTA/Q file manipulations.

https://bioinf.shenwei.me/seqkit/

https://bioinf.shenwei.me/seqkit/tutorial/



## Run BLAST
### Make BLAST DB

```bash
makeblastdb -in your-nucleotide-db.fa -dbtype nucl ###for nucleotide sequence
```

```bash
makeblastdb -in your-protein-db.fas -dbtype prot ###for protein sequence
```

### Run BLASTX
```bash
cd /data/gpfs/assoc/bch709-3/<YOURID>/BLAST
gunzip plant.1.protein.faa.gz
makeblastdb -in plant.1.protein.faa -dbtype prot
seqkit sample -n 100 /data/gpfs/assoc/bch709-3/Course_materials/BLAST/Athaliana_167_TAIR10.cds.fa.gz > ATH_100.fasta
blastx -query ATH_100.fasta  -db plant.1.protein.faa -outfmt 8
```


### Tab output

	qseqid 		Query sequence ID
	sseqid		Subject (ie DB) sequence ID
	pident		Percent Identity across the alignment
	length 		Alignment length
	mismatch 	# of mismatches
	gapopen 	Number of gap openings
	qstart 		Start of alignment in query
	qend 		End of alignment in query 
	sstart 		Start of alignment in subject
	send		End of alignment in subject
	evalue 		E-value
	bitscore	Bit score


># Question
- find the option below within BLASTX
1. Set output to file
2. Set tabular output format
3. Set maximum target sequence to one
4. Set threads (CPU) to 32
5. Set evalue threshold to 1e-30
{: .prereq}



# DCBLAST

The Basic Local Alignment Search Tool (BLAST) is by far best the most widely used tool in for sequence analysis for rapid sequence similarity searching among nucleic acid or amino acid sequences. Recently, cluster, HPC, grid, and cloud environmentshave been are increasing more widely used and more accessible as high-performance computing systems. Divide and Conquer BLAST (DCBLAST) has been designed to perform run on grid system with query splicing which can run National Center for Biotechnology Information (NCBI) BLASTBLAST search comparisons  over withinthe cluster, grid, and cloud computing grid environment by using a query sequence distribution approach NCBI BLAST. This is a promising tool to accelerate BLAST job dramatically accelerates the execution of BLAST query searches using a simple, accessible, robust, and practical approach. 

- DCBLAST can run BLAST job across HPC.
- DCBLAST suppport all NCBI-BLAST+ suite.
- DCBLAST generate exact same NCBI-BLAST+ result.
- DCBLAST can use all options in NCBI-BLAST+ suite.


![blast](https://raw.githubusercontent.com/wyim-pgl/DCBLAST/master/fig/fig-1-2x.jpg)


## Requirement
Following basic softwares are needed to run

- Perl (Any version 5+)

```bash
which perl
perl --version
```

- NCBI-BLAST+ (Any version)
for easy approach, you can download binary version of blast from below link.
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST

For using recent version, please update BLAST path in config.ini

```bash
which blastn
```

- Sun Grid Engine (Any version)
```bash
which qsub
```

- Slurm
```bash
which sbatch
```

- Grid cloud or distributed computing system.


## Prerequisites

The following Perl modules are required:
```bash
- Path::Tiny
- Data::Dumper
- Config::Tiny
```
Install prerequisites with the following command:

```bash
cpan `cat requirement`
```
or
```bash
cpanm `cat requirement`
```
or 
```bash
cpanm Path::Tiny Data::Dumper Config::Tiny
```
We strongly recommend to use Perlbrew http://perlbrew.pl/ to avoid having to type sudo

We also recommend to use 'cpanm' https://github.com/miyagawa/cpanminus

## Prerequisites by Conda

```bash
conda activate blast
conda install -c bioconda perl-path-tiny blast perl-data-dumper perl-config-tiny -y
```

## Installation

The program is a single file Perl scripts. Copy it into executive directories.

We recommend to copy it on scratch disk.


```bash
cp /data/gpfs/assoc/bch709-3/Course_materials/BLAST/Athaliana_167_TAIR10.cds.fa.gz /data/gpfs/assoc/bch709-3/<YOURID>/BLAST

gunzip /data/gpfs/assoc/bch709-3/<YOURID>/BLAST/Athaliana_167_TAIR10.cds.fa.gz 


cd /data/gpfs/assoc/bch709-3/<YOURID>/BLAST

mkdir /data/gpfs/assoc/bch709-3/<YOURID>/BLAST/

git clone git@github.com:wyim-pgl/DCBLAST.git

cd DCBLAST/DCBLAST-SLURM

pwd

chmod 775 dcblast.pl

perl dcblast.pl
```

### Help

```bash
Usage : dcblast.pl --ini config.ini --input input-fasta --size size-of-group --output output-filename-prefix  --blast blast-program-name

  --ini <ini filename> ##config file ex)config.ini

  --input <input filename> ##query fasta file

  --size <output size> ## size of chunks usually all core x 2, if you have 160 core all nodes, you can use 320. please check it to your admin.

  --output <output filename> ##output folder name

  --blast <blast name> ##blastp, blastx, blastn and etcs.

  --dryrun Option will only split fasta file into chunks
```


### Configuration

**Please edit config.ini with `nano` before you run!!**

```bash
[dcblast]
##Name of job (will use for SGE job submission name)
job_name_prefix=dcblast

[blast]
##BLAST options

##BLAST path (your blast+ path); $ which blastn; then remove "blastn"
path=~/miniconda3/envs/blast/bin/

##DB path (build your own BLAST DB)
##example
##makeblastdb -in example/test_db.fas -dbtype nucl (for nucleotide sequence)
##makeblastdb -in example/your-protein-db.fas -dbtype prot (for protein sequence)
db=/data/gpfs/assoc/bch709-3/<YOURID>/BLAST/plant.1.protein.faa 

##Evalue cut-off (See BLAST manual)
evalue=1e-05

##number of threads in each job. If your CPU is AMD it needs to be set 1.
num_threads=2

##Max target sequence output (See BLAST manual)
max_target_seqs=10

##Output format (See BLAST manual)
outfmt=6

##any other option can be add it this area
#matrix=BLOSUM62
#gapopen=11
#gapextend=1

[oldsge]
##Grid job submission commands
##please check your job submission scripts
##Especially Queue name (q) and Threads option (pe) will be different depends on your system

pe=SharedMem 1
M=your@email
q=common.q
j=yes
o=log
cwd=

[slurm]
time=04:00:00
cpus-per-task=1
mem-per-cpu=800M
ntasks=1
output=log
error=error
partition=cpu-core-0
account=cpu-s5-bch709-3
mail-type=all
mail-user=<YOURID>@unr.edu

```
If you need any other options for your enviroment please contant us or admin

PBS & LSF need simple code hack. If you need it please request through issue.

## Run DCBLAST

### Run (--dryrun option will only split fasta file into chunks)
```bash
perl dcblast.pl --ini config.ini --input /data/gpfs/assoc/bch709-3/<YOURID>/BLAST/Athaliana_167_TAIR10.cds.fa --output test --size 100 --blast blastx 
```

```bash
squeue
```

**It usually finish within up to 20min depends on HPC status and CPU speed.**



## Citation
Won C. Yim and John C. Cushman (2017) Divide and Conquer BLAST: using grid engines to accelerate BLAST and other sequence analysis tools. PeerJ 10.7717/peerj.3486 https://peerj.com/articles/3486/

