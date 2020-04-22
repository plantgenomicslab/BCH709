---
layout: page
title: 13_Sequence_Alignment
published: true
---

## Pairwise alignment
The process or result of matching up the nucleotide or amino acid residues of two or more biological sequences to achieve maximal levels of identity and, in the case of amino acid sequences, conservation, for the purpose of assessing the degree of similarity and the possibility of homology.
 
One of the most fundamental problems in bioinformatics is determining how "similar" a pair of biological sequences are. There are many applications for this, including inferring the function or source organism of an unknown gene sequence, developing hypotheses about the relatedness of organisms, or grouping sequences from closely related organisms. On the surface this seems like a pretty straight-forward problem, not one that would have been at the center of decades of research and the subject of one of the most cited papers in modern biology. In this chapter we'll explore why determining sequence similarity is harder than it might initially seem, and learn about pairwise sequence alignment, the standard approach for determining sequence similarity.

## Protein vs. DNA sequence alignment
### Protein amino acid sequences are preferred over DNA sequences for a list of reasons.

- Protein residues are more informative - a change in DNA (especially the 3rd position) does not necessarily change the AA.
- The larger number of amino acids than nucleic acids makes it easier to find significance.
- Some amino acids share related biochemical properties, which can be accounted for when scoring multiple pairwise alignments.
- Protein sequence comparisons can link back to over a billion years ago, whereas DNA sequence comparisons can only go back up to 600 mya. Thus, protein sequences are far better for evolutionary studies.

### However, there are some obvious instances when DNA alignments are needed.

- When confirming the identity of cDNA (forensic sequencing).
- When studying noncoding regions of DNA. These regions evolve at a faster rate than coding DNA, while mitochondrial noncoding DNA evolves even faster.
- When studying DNA mutations.
- When researching on very similar organisms such as Neanderthals and modern humans.



![insertion]({{site.baseurl}}/fig/insertion.png)
![mismatch]({{site.baseurl}}/fig/mismatch.png)
![deletion]({{site.baseurl}}/fig/deletion.png)

## Hamming Distance 

### Metric to compare sequences (DNA, AA, ASCII, binary, etc…)
- Non-negative, identity, symmetry, triangle equality
- How many characters are different between the 2 strings?
- Minimum number of substitutions required to change transform A into B
- Traditionally defined for end-to-end comparisons

![hamming]({{site.baseurl}}/fig/hamming.png)


### Here end-to-end (global) for query, partial (local) for reference
- Find all occurrences of GATTACA with Hamming Distance ≤ 1
 ![editdistance]({{site.baseurl}}/fig/editdistance.png)

## Local vs Global alignment
![ALIGNMENT]({{site.baseurl}}/fig/ALIGNMENT.png)

![local_example]({{site.baseurl}}/fig/local_example.png)


The *Smith–Waterman* algorithm finds the segments in two sequences that have similarities while the *Needleman–Wunsch* algorithm aligns two complete sequences. Therefore, they serve different purposes. Both algorithms use the concepts of a substitution matrix, a gap penalty function, a scoring matrix, and a traceback process. 

![alignment_compare]({{site.baseurl}}/fig/alignment_compare.png)
![gap]({{site.baseurl}}/fig/gap.png)
![transitions]({{site.baseurl}}/fig/transitions.png)


Take the alignment of DNA sequences TGTTACGG and GGTTGACTA as an example. Use the following scheme:


![penalty]({{site.baseurl}}/fig/penalty.png)

Initialize and fill the scoring matrix, shown as below. This figure shows the scoring process of the first three elements. The yellow color indicates the bases that are being considered. The red color indicates the highest possible score for the cell being scored.


![Smith-Waterman-Algorithm-Example-Step1]({{site.baseurl}}/fig/Smith-Waterman-Algorithm-Example-Step1.png)
![matrix_finish]({{site.baseurl}}/fig/matrix_finish.png)
![backtrace]({{site.baseurl}}/fig/backtrace.png)

## location
```
 mkdir /data/gpfs/assoc/bch709/<YOURID>/alignment
 cd $!
```

## Install software 
create `alignment` env, activate and install `blast muscle clustalo` from bioconda


>## Actual Alignment
> Match score:	3   
> Mismatch score: -3  
> Gap penalty:	-2               
>TTAGTAT  
>TAGCTAT
{: .callout}


| |T|T|A|G|T|A|T|
|---|---|---|---|---|---|---|---|
|T| | | | | | | |
|A| | | | | | | |
|G| | | | | | | |
|C| | | | | | | |
|T| | | | | | | |
|A| | | | | | | |
|T| | | | | | | |


[.](https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/data/Needleman-Wunsch%20Algorithm_BCH709.xls)


## Example Protein sequences
- save it as example.aa
```
>sp|P10242|MYB_HUMAN Transcriptional activator Myb OS=Homo sapiens OX=9606 GN=MYB PE=1 SV=2
MARRPRHSIYSSDEDDEDFEMCDHDYDGLLPKSGKRHLGKTRWTREEDEKLKKLVEQNGT
DDWKVIANYLPNRTDVQCQHRWQKVLNPELIKGPWTKEEDQRVIELVQKYGPKRWSVIAK
HLKGRIGKQCRERWHNHLNPEVKKTSWTEEEDRIIYQAHKRLGNRWAEIAKLLPGRTDNA
IKNHWNSTMRRKVEQEGYLQESSKASQPAVATSFQKNSHLMGFAQAPPTAQLPATGQPTV
NNDYSYYHISEAQNVSSHVPYPVALHVNIVNVPQPAAAAIQRHYNDEDPEKEKRIKELEL
LLMSTENELKGQQVLPTQNHTCSYPGWHSTTIADHTRPHGDSAPVSCLGEHHSTPSLPAD
PGSLPEESASPARCMIVHQGTILDNVKNLLEFAETLQFIDSFLNTSSNHENSDLEMPSLT
STPLIGHKLTVTTPFHRDQTVKTQKENTVFRTPAIKRSILESSPRTPTPFKHALAAQEIK
YGPLKMLPQTPSHLVEDLQDVIKQESDESGIVAEFQENGPPLLKKIKQEVESPTDKSGNF
FCSHHWEGDSLNTQLFTQTSPVADAPNILTSSVLMAPASEDEDNVLKAFTVPKNRSLASP
LQPCSSTWEPASCGKMEEQMTSSSQARKYVNAFSARTLVM

>sp|P46200|MYB_BOVIN Transcriptional activator Myb OS=Bos taurus OX=9913 GN=MYB PE=2 SV=1
MARRPRHSIYSSDEDDEDIEMCDHDYDGLLPKSGKRHLGKTRWTREEDEKLKKLVEQNGT
DDWKVIANYLPNRTDVQCQHRWQKVLNPELIKGPWTKEEDQRVIELVQKYGPKRWSVIAK
HLKGRIGKQCRERWHNHLNPEVKKTSWTEEEDRIIYQAHKRLGNRWAEIAKLLPGRTDNA
IKNHWNSTMRRKVEQEGYLQESSKASQPAVTTSFQKNSHLMGFTHAPPSAQLPPAGQPSV
NSDYPYYHISEAQNVSSHVPYPVALHVNIVNVPQPAAAAIQRHYNDEDPEKEKRIKELEL
LLMSTENELKGQQALPTQNHTCSYPGWHSTTIADHTRPHGDSAPVSCLEEHHSTPSLPAD
PGSLPEESASPARCMIFHQSTILDNVKNLLEFAETLQFIDSFLNTSNNHENLDLEMPSLT
STPLNGHKLTVTTPFHRDQTVKIQKENTIFRTPAIKRSILEGSPRTPTPFKHALTAQEIK
YGPLKMLPQTPSHLVEDLQDEIKQESDESGIVAEFQENGQPLLKKIKQEVESPTDKAGNF
FCSNHWEGDSLNTQLFTQASPVADMPNILTSSVLMTPVSEDEDNVLKAFTVPKSRSLASP
LQPCNGAWESASCGKTDDQMTASGQSRKYVNAFSTRTLVM

>tr|A0A2I2ZIL5|A0A2I2ZIL5_GORGO MYB proto-oncogene, transcription factor OS=Gorilla gorilla gorilla OX=9595 PE=4 SV=1
MARRPRHSIYSSDEDDEDFEMCDHDYDGLLPKSGKRHLGKTRWTREEDEKLKKLVEQNGT
DDWKVIANYLPNRTDVQCQHRWQKVLNPELIKGPWTKEEDQRVIELVQKYGPKRWSVIAK
HLKGRIGKQCRERWHNHLNPEVKKTSWTEEEDRIIYQAHKRLGNRWAEIAKLLPGRTDNA
IKNHWNSTMRRKVEQEGYLQESSKASQPAVATSFQKNSHLMGFAQAPPTAQLPATGQPTV
NNDYSYYHISEAQNVSSHVPYPVALHVNIVNVPQPAAAAIQRHYNDEDPEKEKRIKELEL
LLMSTENELKGQQVLPTQNHTCSYPGWHSTTIADHTRPHGDSAPVSCLGEHHSTPSLPAD
PGSLPEESASPARCMIVHQGTILDNVKNLLEFAETLQFIDSFLNTSSNHENSDLEMPSLT
STPLIGHKLTVTTPFHRDQTVKTQKENTVFRTPAIKRSILESSPRTPTPFKHALAAQEIK
YGPLKMLPQTPSHLVEDLQDVIKQESDESGIVAEFQENGPPLLKKIKQEVESPTDKSGNF
FCSHHWEGDGLNTQLFTQTSPVADAPNILTSSVLMAPASEDEDNVLKAFTVPKNRSLASP
LQPCSSTWEPASCGKLEEQMTASSQARKYVNAFSARTLVM

>XP_012418097.1 PREDICTED: transcriptional activator Myb isoform X5 [Odobenus rosmarus divergens]
MARRPRHSIYSSDEDDEDIEMCDHDYDGLLPKSGKRHLGKTRWTREEDEKLKKLVEQNGTDDWKVIANYL
PNRTDVQCQHRWQKVLNPELIKGPWTKEEDQRVIELVQKYGPKRWSVIAKHLKGRIGKQCRERWHNHLNP
EVKKTSWTEEEDRIIYQAHKRLGNRWAEIAKLLPGRTDNAIKNHWNSTMRRKVEQEGYLQESSKASPPTV
ATSFQKNSHLMGFAHAPPSAHLPPAGQPSVNNDYSYYHISEAQNVSSHVPYPVALHVNIVNVPQPAAAAI
QRHYNDEDPEKEKRIKELELLLMSTENELKGQQTQNHTCSYPGWHSTTIADHTRPHGDSAPVSCLEEHHS
TPSLPVDPGSLPEESASPARCMIVHQGTILDNVKNLLEFAETLQFIDSFLNTSGNHENLDLEMPSLTSTP
LNGHKLTVTTPFHRDQTVKTQKENTIFRTPAIKRSILESSPRTPTPFKHGLAAQEIKYGPLKMLPQTPSH
LVEDLQDVIKQESDEPGIVAEFQENGPPLLKKIKQEVESPTDKAGNFFCSSHWEGESLNTQLFPQALPVT
DVPNILTSSVLMTPVSEDEDNVLKAFTVPKNRSLASPLQPCGGAWEAASCGKTEDQMTASGQARKYVNAF
STRTLVM


```

### Example Nucleotide Sequence
- save it as example.fasta
```
>NM_001130173.1
ATGGCCCGAAGACCCCGGCACAGCATATATAGCAGTGACGAGGATGATGAGGACTTTGAGATGTGTGACC
ATGACTATGATGGGCTGCTTCCCAAGTCTGGAAAGCGTCACTTGGGGAAAACAAGGTGGACCCGGGAAGA
GGATGAAAAACTGAAGAAGCTGGTGGAACAGAATGGAACAGATGACTGGAAAGTTATTGCCAATTATCTC
CCGAATCGAACAGATGTGCAGTGCCAGCACCGATGGCAGAAAGTACTAAACCCTGAGCTCATCAAGGGTC
CTTGGACCAAAGAAGAAGATCAGAGAGTGATAGAGCTTGTACAGAAATACGGTCCGAAACGTTGGTCTGT
TATTGCCAAGCACTTAAAGGGGAGAATTGGAAAACAATGTAGGGAGAGGTGGCATAACCACTTGAATCCA
GAAGTTAAGAAAACCTCCTGGACAGAAGAGGAAGACAGAATTATTTACCAGGCACACAAGAGACTGGGGA
ACAGATGGGCAGAAATCGCAAAGCTACTGCCTGGACGAACTGATAATGCTATCAAGAACCACTGGAATTC
TACAATGCGTCGGAAGGTCGAACAGGAAGGTTATCTGCAGGAGTCTTCAAAAGCCAGCCAGCCAGCAGTG
GCCACAAGCTTCCAGAAGAACAGTCATTTGATGGGTTTTGCTCAGGCTCCGCCTACAGCTCAACTCCCTG
CCACTGGCCAGCCCACTGTTAACAACGACTATTCCTATTACCACATTTCTGAAGCACAAAATGTCTCCAG
TCATGTTCCATACCCTGTAGCGTTACATGTAAATATAGTCAATGTCCCTCAGCCAGCTGCCGCAGCCATT
CAGAGACACTATAATGATGAAGACCCTGAGAAGGAAAAGCGAATAAAGGAATTAGAATTGCTCCTAATGT
CAACCGAGAATGAGCTAAAAGGACAGCAGACACAGAACCACACATGCAGCTACCCCGGGTGGCACAGCAC
CACCATTGCCGACCACACCAGACCTCATGGAGACAGTGCACCTGTTTCCTGTTTGGGAGAACACCACTCC
ACTCCATCTCTGCCAGCGGATCCTGGCTCCCTACCTGAAGAAAGCGCCTCGCCAGCAAGGTGCATGATCG
TCCACCAGGGCACCATTCTGGATAATGTTAAGAACCTCTTAGAATTTGCAGAAACACTCCAATTTATAGA
TTCTTTCTTAAACACTTCCAGTAACCATGAAAACTCAGACTTGGAAATGCCTTCTTTAACTTCCACCCCC
CTCATTGGTCACAAATTGACTGTTACAACACCATTTCATAGAGACCAGACTGTGAAAACTCAAAAGGAAA
ATACTGTTTTTAGAACCCCAGCTATCAAAAGGTCAATCTTAGAAAGCTCTCCAAGAACTCCTACACCATT
CAAACATGCACTTGCAGCTCAAGAAATTAAATACGGTCCCCTGAAGATGCTACCTCAGACACCCTCTCAT
CTAGTAGAAGATCTGCAGGATGTGATCAAACAGGAATCTGATGAATCTGGAATTGTTGCTGAGTTTCAAG
AAAATGGACCACCCTTACTGAAGAAAATCAAACAAGAGGTGGAATCTCCAACTGATAAATCAGGAAACTT
CTTCTGCTCACACCACTGGGAAGGGGACAGTCTGAATACCCAACTGTTCACGCAGACCTCGCCTGTGGCA
GATGCACCGAATATTCTTACAAGCTCCGTTTTAATGGCACCAGCATCAGAAGATGAAGACAATGTTCTCA
AAGCATTTACAGTACCTAAAAACAGGTCCCTGGCGAGCCCCTTGCAGCCTTGTAGCAGTACCTGGGAACC
TGCATCCTGTGGAAAGATGGAGGAGCAGATGACATCTTCCAGTCAAGCTCGTAAATACGTGAATGCATTC
TCAGCCCGGACGCTGGTCATGTGA

>NM_005375.3
ATGGCCCGAAGACCCCGGCACAGCATATATAGCAGTGACGAGGATGATGAGGACTTTGAGATGTGTGACC
ATGACTATGATGGGCTGCTTCCCAAGTCTGGAAAGCGTCACTTGGGGAAAACAAGGTGGACCCGGGAAGA
GGATGAAAAACTGAAGAAGCTGGTGGAACAGAATGGAACAGATGACTGGAAAGTTATTGCCAATTATCTC
CCGAATCGAACAGATGTGCAGTGCCAGCACCGATGGCAGAAAGTACTAAACCCTGAGCTCATCAAGGGTC
CTTGGACCAAAGAAGAAGATCAGAGAGTGATAGAGCTTGTACAGAAATACGGTCCGAAACGTTGGTCTGT
TATTGCCAAGCACTTAAAGGGGAGAATTGGAAAACAATGTAGGGAGAGGTGGCATAACCACTTGAATCCA
GAAGTTAAGAAAACCTCCTGGACAGAAGAGGAAGACAGAATTATTTACCAGGCACACAAGAGACTGGGGA
ACAGATGGGCAGAAATCGCAAAGCTACTGCCTGGACGAACTGATAATGCTATCAAGAACCACTGGAATTC
TACAATGCGTCGGAAGGTCGAACAGGAAGGTTATCTGCAGGAGTCTTCAAAAGCCAGCCAGCCAGCAGTG
GCCACAAGCTTCCAGAAGAACAGTCATTTGATGGGTTTTGCTCAGGCTCCGCCTACAGCTCAACTCCCTG
CCACTGGCCAGCCCACTGTTAACAACGACTATTCCTATTACCACATTTCTGAAGCACAAAATGTCTCCAG
TCATGTTCCATACCCTGTAGCGTTACATGTAAATATAGTCAATGTCCCTCAGCCAGCTGCCGCAGCCATT
CAGAGACACTATAATGATGAAGACCCTGAGAAGGAAAAGCGAATAAAGGAATTAGAATTGCTCCTAATGT
CAACCGAGAATGAGCTAAAAGGACAGCAGGTGCTACCAACACAGAACCACACATGCAGCTACCCCGGGTG
GCACAGCACCACCATTGCCGACCACACCAGACCTCATGGAGACAGTGCACCTGTTTCCTGTTTGGGAGAA
CACCACTCCACTCCATCTCTGCCAGCGGATCCTGGCTCCCTACCTGAAGAAAGCGCCTCGCCAGCAAGGT
GCATGATCGTCCACCAGGGCACCATTCTGGATAATGTTAAGAACCTCTTAGAATTTGCAGAAACACTCCA
ATTTATAGATTCTGATTCTTCATCATGGTGTGATCTCAGCAGTTTTGAATTCTTTGAAGAAGCAGATTTT
TCACCTAGCCAACATCACACAGGCAAAGCCCTACAGCTTCAGCAAAGAGAGGGCAATGGGACTAAACCTG
CAGGAGAACCTAGCCCAAGGGTGAACAAACGTATGTTGAGTGAGAGTTCACTTGACCCACCCAAGGTCTT
ACCTCCTGCAAGGCACAGCACAATTCCACTGGTCATCCTTCGAAAAAAACGGGGCCAGGCCAGCCCCTTA
GCCACTGGAGACTGTAGCTCCTTCATATTTGCTGACGTCAGCAGTTCAACTCCCAAGCGTTCCCCTGTCA
AAAGCCTACCCTTCTCTCCCTCGCAGTTCTTAAACACTTCCAGTAACCATGAAAACTCAGACTTGGAAAT
GCCTTCTTTAACTTCCACCCCCCTCATTGGTCACAAATTGACTGTTACAACACCATTTCATAGAGACCAG
ACTGTGAAAACTCAAAAGGAAAATACTGTTTTTAGAACCCCAGCTATCAAAAGGTCAATCTTAGAAAGCT
CTCCAAGAACTCCTACACCATTCAAACATGCACTTGCAGCTCAAGAAATTAAATACGGTCCCCTGAAGAT
GCTACCTCAGACACCCTCTCATCTAGTAGAAGATCTGCAGGATGTGATCAAACAGGAATCTGATGAATCT
GGAATTGTTGCTGAGTTTCAAGAAAATGGACCACCCTTACTGAAGAAAATCAAACAAGAGGTGGAATCTC
CAACTGATAAATCAGGAAACTTCTTCTGCTCACACCACTGGGAAGGGGACAGTCTGAATACCCAACTGTT
CACGCAGACCTCGCCTGTGGCAGATGCACCGAATATTCTTACAAGCTCCGTTTTAATGGCACCAGCATCA
GAAGATGAAGACAATGTTCTCAAAGCATTTACAGTACCTAAAAACAGGTCCCTGGCGAGCCCCTTGCAGC
CTTGTAGCAGTACCTGGGAACCTGCATCCTGTGGAAAGATGGAGGAGCAGATGACATCTTCCAGTCAAGC
TCGTAAATACGTGAATGCATTCTCAGCCCGGACGCTGGTCATGTGA

>NM_001130172.1
ATGGCCCGAAGACCCCGGCACAGCATATATAGCAGTGACGAGGATGATGAGGACTTTGAGATGTGTGACC
ATGACTATGATGGGCTGCTTCCCAAGTCTGGAAAGCGTCACTTGGGGAAAACAAGGTGGACCCGGGAAGA
GGATGAAAAACTGAAGAAGCTGGTGGAACAGAATGGAACAGATGACTGGAAAGTTATTGCCAATTATCTC
CCGAATCGAACAGATGTGCAGTGCCAGCACCGATGGCAGAAAGTACTAAACCCTGAGCTCATCAAGGGTC
CTTGGACCAAAGAAGAAGATCAGAGAGTGATAGAGCTTGTACAGAAATACGGTCCGAAACGTTGGTCTGT
TATTGCCAAGCACTTAAAGGGGAGAATTGGAAAACAATGTAGGGAGAGGTGGCATAACCACTTGAATCCA
GAAGTTAAGAAAACCTCCTGGACAGAAGAGGAAGACAGAATTATTTACCAGGCACACAAGAGACTGGGGA
ACAGATGGGCAGAAATCGCAAAGCTACTGCCTGGACGAACTGATAATGCTATCAAGAACCACTGGAATTC
TACAATGCGTCGGAAGGTCGAACAGGAAGGTTATCTGCAGGAGTCTTCAAAAGCCAGCCAGCCAGCAGTG
GCCACAAGCTTCCAGAAGAACAGTCATTTGATGGGTTTTGCTCAGGCTCCGCCTACAGCTCAACTCCCTG
CCACTGGCCAGCCCACTGTTAACAACGACTATTCCTATTACCACATTTCTGAAGCACAAAATGTCTCCAG
TCATGTTCCATACCCTGTAGCGTTACATGTAAATATAGTCAATGTCCCTCAGCCAGCTGCCGCAGCCATT
CAGAGACACTATAATGATGAAGACCCTGAGAAGGAAAAGCGAATAAAGGAATTAGAATTGCTCCTAATGT
CAACCGAGAATGAGCTAAAAGGACAGCAGGTGCTACCAACACAGAACCACACATGCAGCTACCCCGGGTG
GCACAGCACCACCATTGCCGACCACACCAGACCTCATGGAGACAGTGCACCTGTTTCCTGTTTGGGAGAA
CACCACTCCACTCCATCTCTGCCAGCGGATCCTGGCTCCCTACCTGAAGAAAGCGCCTCGCCAGCAAGGT
GCATGATCGTCCACCAGGGCACCATTCTGGATAATGTTAAGAACCTCTTAGAATTTGCAGAAACACTCCA
ATTTATAGATTCTTTCTTAAACACTTCCAGTAACCATGAAAACTCAGACTTGGAAATGCCTTCTTTAACT
TCCACCCCCCTCATTGGTCACAAATTGACTGTTACAACACCATTTCATAGAGACCAGACTGTGAAAACTC
AAAAGGAAAATACTGTTTTTAGAACCCCAGCTATCAAAAGGTCAATCTTAGAAAGCTCTCCAAGAACTCC
TACACCATTCAAACATGCACTTGCAGCTCAAGAAATTAAATACGGTCCCCTGAAGATGCTACCTCAGACA
CCCTCTCATCTAGTAGAAGATCTGCAGGATGTGATCAAACAGGAATCTGATGAATCTGGAATTGTTGCTG
AGTTTCAAGAAAATGGACCACCCTTACTGAAGAAAATCAAACAAGAGGTGGAATCTCCAACTGATAAATC
AGGAAACTTCTTCTGCTCACACCACTGGGAAGGGGACAGTCTGAATACCCAACTGTTCACGCAGACCTCG
CCTGTGGCAGATGCACCGAATATTCTTACAAGCTCCGTTTTAATGGCACCAGCATCAGAAGATGAAGACA
ATGTTCTCAAAGCATTTACAGTACCTAAAAACAGGTCCCTGGCGAGCCCCTTGCAGCCTTGTAGCAGTAC
CTGGGAACCTGCATCCTGTGGAAAGATGGAGGAGCAGATGACATCTTCCAGTCAAGCTCGTAAATACGTG
AATGCATTCTCAGCCCGGACGCTGGTCATGTGA
```


## Muscle

![muscle]({{site.baseurl}}/fig/muscle.jpg)
Iterative Progressive Alignment

```bash
muscle --help
muscle -in example.aa -out example.aa.muscle
muscle -in example.aa -out example.aa.muscle.clw -clw
muscle -in example.fasta -out example.fasta.muscle
muscle -in example.fasta -out example.fasta.muscle.clw -clw
```


## ClustalO
![HMM]({{site.baseurl}}/fig/HMM.jpeg)
Hidden Markov Model

![Global_vs_Local]({{site.baseurl}}/fig/Global_vs_Local.png)

![local]({{site.baseurl}}/fig/local.png)


https://www.ebi.ac.uk/Tools/msa/clustalo/

```bash
clustalo --help
clustalo -i example.aa -o example.aa.clustao
clustalo -i example.aa -o example.aa.clustao.clw --outfmt=clu
clustalo -i example.fasta -o example.fasta.clustalo
clustalo -i example.fasta -o example.fasta.clustalo.clw --outfmt=clu

```



### Alignmnet Visualization
https://www.ebi.ac.uk/Tools/msa/mview/



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

** E[# occurrences of a string of length m in reference of length L] ~ L/4m **



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

Select and align all proteins. Can you explain the differences?

What do you think of the Bos taurus sequence (A8YXY3) and the pig sequence (A1Z623)?


### Running a standalone BLAST program
Create the index for the target database using makeblastdb;
Choose the task program: blastn, blastp, blastx, tblatx, psiblast or deltablast;
Set the configuration for match, mismatch, gap-open penalty, gap-extension penalty or scoring matrix;
Set the word size;
Set the E-value threshold;
Set the output format and the number of output results


### Standalone BLAST
1. Download the database.
2. Use makeblastdb to build the index.
3. Change the scoring matrix, record the changes in the alignment results and interpret the results.

### Download Database
```
ftp://ftp.ncbi.nih.gov/refseq/release/plant/plant.1.protein.faa.gz
```
### How many sequences in `plant.1.protein.faa.gz`


### Run BLAST
```
makeblastdb -in plant.1.protein.faa -dbtype prot
blastx -query /data/gpfs/assoc/bch709/<YOURID>/rnaseq_slurm/trinity_out_dir/Trinity.fasta -db plant.1.protein.faa 
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


### Assignment Run standalone BLASTX
- find output to file, tabular output format, set maximu target sequence to one, and threads (CPU) to four 
