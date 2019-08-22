---
layout: page
title: CRISPR/Cas9 database build
published: true
---
# *I. scapularis* genome for CRISPR/Cas9
Convert your genome file to ChopChop format
Requested by Dr. Monika Gulia-Nuss and Dr. Arvind Sharma

[Gulia-Nuss Lab.](https://naes.unr.edu/gulia)


---


![Dr. Monika Gulia-Nuss](https://naes.unr.edu/gulia/wp-content/uploads/graduation.jpg)
Figure 1. Jeremiah Reyes (Grad special), Preston, **Dr. Arvind Sharma**, **Dr. Monika Gulia-Nuss** and Dr. Andrew Nuss.

![Ticks](https://media.springernature.com/m685/nature-static/assets/v1/image-assets/ncomms10507-f1.jpg)
Figure 2. Genes associated with the unique parasitic lifestyle of *Ixodes scapularis*.
[Genomic insights into the Ixodes scapularis tick vector of Lyme disease](https://www.nature.com/articles/ncomms10507) 
 
## ChopChop
http://chopchop.cbu.uib.no/

## Genome submission instructions
https://chopchop.cbu.uib.no/submissions

## Citations for ChopChop
Kornel Labun; Tessa G. Montague; James A. Gagnon; Summer B. Thyme; Eivind Valen. (2016). CHOPCHOP v2: a web tool for the next generation of CRISPR genome engineering. Nucleic Acids Research; doi:10.1093/nar/gkw398

Tessa G. Montague; Jose M. Cruz; James A. Gagnon; George M. Church; Eivind Valen. (2014). CHOPCHOP: a CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic Acids Res. 42. W401-W407


## Requirement 
* Any linux system or Any terminal
* Web browser with internet
* Your passion


## Processing
### 1. Download Tick genome from Vectorbase
[https://www.vectorbase.org/downloads?field_organism_taxonomy_tid%5B%5D=340&field_download_file_type_tid%5B%5D=412&field_download_file_format_tid=All&field_status_value=Current](https://https://www.vectorbase.org/downloads?field_organism_taxonomy_tid%5B%5D=340&field_download_file_type_tid%5B%5D=412&field_download_file_format_tid=All&field_status_value=Current)

#### Prepare following files.

* Gene feature file (GFF)

**Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.6.gff3.gz**

* Genome sequence (scaffold)

**Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fa.gz**

[Temporary download links](https://www.vectorbase.org/sites/default/files/ftp/vbo_archive_20180731_0.zip)

* Download
>
wget https://www.vectorbase.org/sites/default/files/ftp/vbo_archive_20180731_0.zip
{: .bash}

* Unzip data
>
unzip vbo_archive_20180731_0.zip
{: .bash}

* Check the files
use "space bar" and quit with "q" for "less"
``` 
zcat  Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fa.gz | less 
```
```
zcat Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.6.gff3.gz | less
```
* Check the number of scaffold in GFF and scaffold(fa)
```
zcat Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.6.gff3.gz | egrep -c  "sequence-region"
```
```

zcat  Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fa.gz | egrep -c ">"
```
* Result (number of scaffolds in GFF and fasta files)

**369492**


---
### 2. Prepare gff3ToGenePred.
#### gff3ToGenePred is validating and GFF file format to Gp (GenePred) format

*Download gff3ToGenePred
>
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
{: .bash}
  
*Check the excute
>
chmod 775  ./gff3ToGenePred ##permission change for excute
./gff3ToGenePred ##excute 
{: .bash}
---

### 3. Convert gff3 to Gp

*Convert with gff3ToGenePred
>
./gff3ToGenePred Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.6.gff3.gz Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.pred -bad=Ixodes-scapularis.bad
{: .bash}

*Check the BAD file (Should be nothing)
>less Ixodes-scapularis.bad
{: .bash}
```
* Check the Gp file (Should include exon locations, gene name and scaffold ID)
```
$ less Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.pred 
```
---

### 4. Prepare the file for ChopChop
```
mkdir Ixodes-scapularis_ChopChop ## make directory
```

```
mv Ixodes-scapularis-Wikel_BASEFEATURES_IscaW1.pred Ixodes-scapularis_ChopChop ## move files
```

```
mv Ixodes-scapularis-Wikel_SCAFFOLDS_IscaW1.fa.gz Ixodes-scapularis_ChopChop ## move files
```
```
tar -cvzf Ixodes-scapularis_ChopChop.tar.gz Ixodes-scapularis_ChopChop ## compress
```

```
ls -lh Ixodes-scapularis_ChopChop.tar.gz
```
---

### 5. Contact to ChopChop
**The file size is around 453Mb after compress.**
**We need to use BOX to share this files.**
**I prefer to share it by link.**
[How to share it by link](https://community.box.com/t5/Using-Shared-Links/Creating-Shared-Links/ta-p/19523)

Send email to Tessa G. Montague

tessa.chopchop@gmail.com

---

## HOMEWORK for Arvind

1. what is "'|'" ?
2. What is "zcat" ?
3. How to compress and decompress files? explain option (ex cvf, xzvf)
4. Why do we need to "Check the number of scaffold in GFF and scaffold" ?
5. Explain GFF
6. Can you do same thing in *Anopheles stephensi*?
7. Explain about CHOPCHOP statistical approach
