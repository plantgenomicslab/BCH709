---
layout: page
title: Genome description and enrichment test
published: true
---


## Course Evaluation

Students will have access to them until 11:59 PM on Wed, Dec 11, 2019 PST.

**You can log in with your NetID to [http://www.unr.edu/evaluate](http://www.unr.edu/evaluate) and check live updating response rates for your course evaluations. Our institutional goal is to achieve an 85% response rate for all evaluations, and to help us achieve that, we rely on you as well as the students.**


# Assignment due by next Tuesday.

>## Environment 
>```bash
>conda activate genomeannotation
>```
>## Make assignment directory
>```bash
>/data/gpfs/home/wyim/bch709/<YOURID>/genome_transcriptome_comparison
>```
>
>
>
>
>## BLAST comparison
>```bash
>
>cp /data/gpfs/home/wyim/bch709/Course_material/genome_transcriptome_comparison/* .
>```
>
>## Edit bch709.yml first
>*protein fasta name is <YOURFULLPATH>/Trinity.aa*
>
>![ahrd_yml]({{site.baseurl}}/fig/ahrd_yml.png)
>
>
>
>### homework.sh
>```bash
>#!/bin/bash
>#SBATCH --job-name=assignment
>#SBATCH --cpus-per-task=8
>#SBATCH --time=12:00:00
>#SBATCH --mem=20g
>#SBATCH --mail-type=all
>#SBATCH --mail-user=<EMAIL>@unr.edu
>#SBATCH -o assignment.out # STDOUT
>#SBATCH -e assignment.err # STDERR
>#SBATCH -p cpu-s2-core-0 
>#SBATCH -A cpu-s2-bch709-0
>
>TransDecoder.LongOrfs -t Trinity.fasta
>cp Trinity.fasta.transdecoder_dir/longest_orfs.pep Trinity.aa
>
>makeblastdb -in uniprot.faa -dbtype prot
>
>blastp -outfmt 6 -query Trinity.aa -db  uniprot.faa  -out transcriptome_vs_uniprot_blastp.txt -num_threads 8
>
>makeblastdb -in arabidopsis.faa -dbtype prot
>
>blastp -query Trinity.aa -db arabidopsis.faa -out transcriptome_vs_arabidopsis_blastp.txt -num_threads 8 -outfmt 6
>
>
>
>perl blast_parse.pl  -q Trinity.aa  -i transcriptome_vs_arabidopsis_blastp.txt -o transcriptome_vs_arabidopsis_blastp_filtered
>
>java -Xmx20g -jar ahrd.jar bch709.yml
>
>cut -f 1,2 transcriptome_vs_arabidopsis_blastp_filtered | sed 's/\.p.*//g' > id_map
>map_fasta_ids id_map  Trinity.fasta
>
>map_data_ids  id_map RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.WT-UP.subset
>map_data_ids  id_map RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.DT-UP.subset
>map_data_ids  id_map RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.DE.subset
>
>cut -f 1 RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.WT-UP.subset > WT-UP.subset
>cut -f 1 RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.DT-UP.subset > DT-UP.subset
>cut -f 1 RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.DE.subset > all.subset
>```
>## Results
>
>1. compress and send me bch709.transcription.go.tab file
>
>2. Please analyze three files independently ( all.subset, WT-UP.subset and DT-UP.subset ) with metascape and send me.
>
>
>
{: .callout}

### Environment
```bash
conda activate genomeannotation

```

### BLAST comparison
```bash

cp /data/gpfs/assoc/bch709/Course_material/Genome_description_and_enrichment_test/* .
```


### blast.sh
```bash
#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o blast.out # STDOUT
#SBATCH -e blast.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0

makeblastdb -in uniprot.faa -dbtype prot

blastp -outfmt 6 -query bch709.all.maker.proteins.fasta -db  uniprot.faa  -out genome_vs_uniprot_blastp.txt -num_threads 8

makeblastdb -in arabidopsis.faa -dbtype prot

blastp -query bch709.all.maker.proteins.fasta -db arabidopsis.faa -out genome_vs_arabidopsis_blastp.txt -num_threads 8 -outfmt 6

```



## BLAST (Basic Local Alignment Search Tool) 
BLAST is a popular program for searching biosequences against databases. BLAST was developed and is maintained by a group at the National Center for Biotechnology Information (NCBI). Salient characteristics of BLAST are:

### Local alignments
BLAST tries to find patches of regional similarity, rather than trying to find the best alignment between your entire query and an entire database sequence.
### Ungapped alignments
Alignments generated with BLAST do not contain gaps. BLAST's speed and statistical model depend on this, but in theory it reduces sensitivity. However, BLAST will report multiple local alignments between your query and a database sequence.

### Explicit statistical theory
BLAST is based on an explicit statistical theory developed by Samuel Karlin and Steven Altschul (PNAS 87:2284-2268. 1990) The original theory was later extended to cover multiple weak matches between query and database entry PNAS 90:5873. 1993).

CAUTION: the repetitive nature of many biological sequences (particularly naive translations of DNA/RNA) violates assumptions made in the Karlin & Altschul theory. While the P values provided by BLAST are a good rule-of-thumb for initial identification of promising matches, care should be taken to ensure that matches are not due simply to biased amino acid composition.

CAUTION: The databases are contaminated with numerous artifacts. The intelligent use of filters can reduce problems from these sources. Remember that the statistical theory only covers the likelihood of finding a match by chance under particular assumptions; it does not guarantee biological importance.

### Heuristic
BLAST is not guaranteed to find the best alignment between your query and the database; it may miss matches. This is because it uses a strategy which is expected to find most matches, but sacrifices complete sensitivity in order to gain speed. However, in practice few biologically significant matches are missed by BLAST which can be found with other sequence search programs. BLAST searches the database in two phases. First it looks for short subsequences which are likely to produce significant matches, and then it tries to extend these subsequences.
A substitution matrix is used during all phases of protein searches (BLASTP, BLASTX, TBLASTN)
Both phases of the alignment process (scanning & extension) use a substitution matrix to score matches. This is in contrast to FASTA, which uses a substitution matrix only for the extension phase. Substitution matrices greatly improve sensitivity.

## Popular BLAST software
### BLASTP
search a Protein Sequence against a Protein Database.
### BLASTN
search a Nucleotide Sequence against a Nucleotide Database.
### TBLASTN
search a Protein Sequence against a Nucleotide Database, by translating each database Nucleotide sequence in all 6 reading frames.
### BLASTX
search a Nucleotide Sequence against a Protein Database, by first translating the query Nucleotide sequence in all 6 reading frames.

## AHRD

*edit bch709.yml*

### ahrd.sh
```bash
#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH -o blast.out # STDOUT
#SBATCH -e blast.err # STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-0
 java -Xmx8g -jar ahrd.jar bch709.yml
```

## AHRD
High throughput protein function annotation with Human Readable Description (HRDs) and Gene Ontology (GO) Terms. https://www.cropbio.uni-bonn.de/
![AHRD]({{site.baseurl}}/fig/AHRD.jpg)

1. Choose best scoring blast results, 200 from each database searched
2. Filter description lines of above blast-results using regular expressions:
	- Reject those matched by any regex given in e.g. ./test/resources/blacklist_descline.txt,
	- Delete those parts of each description line, matching any regex in e.g. ./test/resources/filter_descline_sprot.txt.
3. Divide each description line into tokens (characters of collective meaning)
	-	In terms of score ignore any tokens matching regexs given e.g. in ./test/resources/blacklist_token.txt.
4. Token score (calculated from: bitscore, database weight, overlap score)
5. Lexical score (calculated from: Token score, High score factor, Pattern factor, Correction factor)
6. Description score (calculated from: Lexical score and Blast score)
7. Choose best scoring description line


### Find best BLAST match

```bash
 perl blast_parse.pl  -q bch709.all.maker.proteins.fasta -i genome_vs_arabidopsis_blastp.txt -o genome_vs_arabidopsis_blastp_filtered
```

### Replace Name of Data

```bash
cut -f 1,2 genome_vs_arabidopsis_blastp_filtered | sed 's/-mRNA-1//g' > id_map
map_fasta_ids id_map  bch709.all.maker.transcripts.fasta
map_fasta_ids id_map  bch709.all.maker.proteins.fasta
map_gff_ids id_map bch709_all.gff

map_data_ids  id_map BCH709.featureCount_count_only.cnt.DT_vs_WT.DESeq2.DE_results.P0.001_C2.DT-UP.subset
map_data_ids  id_map BCH709.featureCount_count_only.cnt.DT_vs_WT.DESeq2.DE_results.P0.001_C2.WT-UP.subset
map_data_ids  id_map BCH709.featureCount_count_only.cnt.DT_vs_WT.DESeq2.DE_results.P0.001_C2.DE.subset
```

![GO]({{site.baseurl}}/fig/GO.png)

### Gene Ontology
Gene Ontology project is a major bioinformatics initiative Gene ontology is an annotation system The project provides the controlled and consistent vocabulary of terms and gene product annotations, i.e. terms occur only once, and there is a dictionary of allowed words
GO describes how gene products behave in a cellular context A consistent description of gene products attributes in terms of their associated biological processes, cellular components and molecular functions in a species-independent manner Each GO term consists of a unique alphanumerical identifier, a common name, synonyms (if applicable), and a definition Each term is assigned to one of the three ontologies Terms have a textual definition When a term has multiple meanings depending on species, the GO uses a "sensu" tag to differentiate among them (trichome differentiation (sensu Magnoliophyta) 


![kegg]({{site.baseurl}}/fig/kegg.png)

### hypergeometric test
The hypergeometric distribution is the lesser-known cousin of the binomial distribution, which describes the probability of k successes in n draws with replacement. The hypergeometric distribution describes probabilities of drawing marbles from the jar without putting them back in the jar after each draw.
The hypergeometric probability mass function is given by (using the original variable convention


![hyper_geo]({{site.baseurl}}/fig/hyper_geo.png)
![combination]({{site.baseurl}}/fig/combination.png)
![FWER]({{site.baseurl}}/fig/FWER.png)

#### FWER
The FWER for the other tests is computed in the same way: the gene-associated variables (scores or counts) are permuted while the annotations of genes to GO-categories stay fixed. Then the statistical tests are evaluated again for every GO-category.


## Hypergeometric Test Example 1
Suppose we randomly select 2 cards without replacement from an ordinary deck of playing cards. What is the probability of getting exactly 2 cards you want (i.e., Ace or 10)?

Solution: This is a hypergeometric experiment in which we know the following:

N = 52; since there are 52 cards in a deck.
k = 16; since there are 16 Ace or 10 cards in a deck.
n = 2; since we randomly select  cards from the deck.
x = 2; since 2 of the cards we select are red.
We plug these values into the hypergeometric formula as follows:

h(x; N, n, k) = [ kCx ] [ N-kCn-x ] / [ NCn ]

h(2; 52, 2, 16) = [ 16C2 ] [ 48C1 ] / [ 52C2 ]

h(2; 52, 2, 16) = [ 325 ] [  1 ] / [ 1,326 ]

h(2; 52, 2, 16) = 0.0904977

Thus, the probability of randomly selecting 2 Ace or 10 cards is 9%

|category|probability|
| -- | -- |
|probability mass f|	0.09049773755656108597285|
|lower cumulative P|	1|
|upper cumulative Q|	0.09049773755656108597285|
|Expectation|	0.6153846153846153846154|


## Hypergeometric Test Example 2
Suppose we have 30 DEGs in human genome (200). What is the probability of getting 10 oncogene?

*An oncogene is a gene that has the potential to cause cancer.*

Solution: This is a hypergeometric experiment in which we know the following:

N = 200; since there are 200 genes in human genome 
k = 10; since there are 10 oncogenes in human
n = 30; since 30 DEGs
x = 5; since 5 of the oncogenes in DEGs.

We plug these values into the hypergeometric formula as follows:

h(x; N, n, k) = [ kCx ] [ N-kCn-x ] / [ NCn ]

h(5; 200, 30, 10) = [ 10C5 ] [ 190C25 ] / [ 200C30 ]

h(5; 200, 30, 10) = [ 252 ] [  11506192278177947613740456466942 ] / [ 409681705022127773530866523638950880 ]

h(5; 200, 30, 10) = 0.007078

Thus, the probability of oncogene 0.7%.

hypergeometry.png

### hypergeometric distribution	value
|category|probability|
| -- | -- |
|probability mass f|	0.0070775932109153651831923063371216961166297|
|lower cumulative P|	0.99903494867072865323201131115533112651846|
|upper cumulative Q|	0.0080426445401867119511809951817905695981658|
|Expectation|	1.5|


### False Discovery Rate (FDR) q-value
The false discovery rate (FDR) is a method of conceptualizing the rate of type I errors in null hypothesis testing when conducting multiple comparisons. FDR-controlling procedures are designed to control the expected proportion of "discoveries" (rejected null hypotheses) that are false (incorrect rejections).

- Benjamini–Yekutieli 
- Benjamini–Hochberg 
- Bonferroni-Selected–Bonferroni
- Bonferroni and Sidak 

### REViGO 
http://revigo.irb.hr/revigo.jsp

### cleverGO 
http://www.tartaglialab.com/GO_analyser/tutorial

### MetaScape
http://metascape.org/gp/index.html#/main/step1

### DAVID
https://david.ncifcrf.gov/







### AHRD References
High quality genome projects that used AHRD:

1. Young, Nevin D., Frédéric Debellé, Giles E. D. Oldroyd, Rene Geurts, Steven B. Cannon, Michael K. Udvardi, Vagner A. Benedito, et al. “The Medicago Genome Provides Insight into the Evolution of Rhizobial Symbioses.” Nature 480, no. 7378 (December 22, 2011): 520–24. doi:10.1038/nature10625.

2. The Tomato Genome Consortium. “The Tomato Genome Sequence Provides Insights into Fleshy Fruit Evolution.” Nature 485, no. 7400 (May 31, 2012): 635–41. doi:10.1038/nature11119.

3. International Wheat Genome Sequencing Consortium (IWGSC). “A Chromosome-Based Draft Sequence of the Hexaploid Bread Wheat (Triticum Aestivum) Genome.” Science (New York, N.Y.) 345, no. 6194 (July 18, 2014): 1251788. doi:10.1126/science.1251788.

4. International Barley Genome Sequencing Consortium, Klaus F. X. Mayer, Robbie Waugh, John W. S. Brown, Alan Schulman, Peter Langridge, Matthias Platzer, et al. “A Physical, Genetic and Functional Sequence Assembly of the Barley Genome.” Nature 491, no. 7426 (November 29, 2012): 711–16. doi:10.1038/nature11543.

5. Wang, W., G. Haberer, H. Gundlach, C. Gläßer, T. Nussbaumer, M. C. Luo, A. Lomsadze, et al. “The Spirodela Polyrhiza Genome Reveals Insights into Its Neotenous Reduction Fast Growth and Aquatic Lifestyle.” Nature Communications 5 (2014): 3311. doi:10.1038/ncomms4311.

6. Guo, Shaogui, Jianguo Zhang, Honghe Sun, Jerome Salse, William J. Lucas, Haiying Zhang, Yi Zheng, et al. “The Draft Genome of Watermelon (Citrullus Lanatus) and Resequencing of 20 Diverse Accessions.” Nature Genetics 45, no. 1 (January 2013): 51–58. doi:10.1038/ng.2470.

AHRD was applied on all plant genomes present in the public database PlantsDB:

7. Spannagl, Manuel, Thomas Nussbaumer, Kai C. Bader, Mihaela M. Martis, Michael Seidel, Karl G. Kugler, Heidrun Gundlach, and Klaus F. X. Mayer. “PGSB PlantsDB: Updates to the Database Framework for Comparative Plant Genome Research.” Nucleic Acids Research 44, no. D1 (January 4, 2016): D1141–47. doi:10.1093/nar/gkv1130.
