---
layout: page
title: 15_RNA-Seq tutorial
published: true
---

## Importing Data from the NCBI Sequence Read Archive (SRA) using the DE



## Course Evaluation

Students will have access to course evaluation
You can log in with your NetID to http://www.unr.edu/evaluate and check live updating response rates for your course evaluations. Our institutional goal is to achieve an 85% response rate for all evaluations, and to help us achieve that, we rely on you as well as the students.

**If we can achieve 100% response rate for evaluation, I will give you additional points for all of you.**


## Discussion is open
https://unr.canvaslms.com/courses/56453/discussion_topics/514795

The due date for the question is November 23rd  11:59pm
The due date for discussion is December 4th  11:59pm


1. Define the biological hypotheses or bottleneck you wish to address which is related to your research, state the approach of your experiment, also state your system, study organism, or study site, and provide justification for what is the goal of your biological hypotheses. Please provide enough background information that the other students can understand your biological hypotheses or bottleneck. If your experiments are complicated, consider briefly explaining the experimental design with reason. If you get more like will get points. (30 Points)

2. Please provide the bioinformatics suggestion that you want to suggest for other people's research hypotheses or bottleneck. It should be scientifically valid methods even if it does not exist. Provide enough information to create an experiment and if you want to create software, please provide reasons and explain what kind of software we need, which part of the hypotheses or bottleneck can be solved. If the software doesn't exist, please provide the design or roadmap of your software. Citation is optional but recommended. Please provide an obstacle to other people's suggestions. In addition, insights and addition will also get points. If you get more like will get points. (10 points per valid answer with reference or concept or hypothesis, a total of 70 Points, seven replies are needed )

3.  The suggestion needs to reply as threaded format.


**Examples are below**

**Example**
 
Tef (or Teff) is a warm season, C4-photosynthesis grass that is gaining popularity in the U.S. as a high-quality summer forage, fodder, and gluten-free grain. However, Tef has relatively tiny seeds compare to other C4 grass.  Currently, the primary goal of my research is to determine the loci of seed color and size. We are trying to use Genome-wide Association Study (GWAS)  https://en.wikipedia.org/wiki/Genome-wide_association_study) (Links to an external site.)
to identify the locus of seed color and size. In our lab, we have 386 teff accession and all of them have different seed colors.  We extracted all of 386 teff accession DNA and sequencing was done. But I don't know how to check the size and colors. The phenotyping is the most important but the main bottleneck of our experiment.  How can we facilitate this task?
 
 
**Student A (This example answer will get 3 points)**
I cannot find any solution but you can use a similar approach such as colony counting.  The accurate counting of plates with high numbers of CFUs is error-prone since it requires a high level of attention by the counter. In the microbiome and general biology field use colony count software to analyze whole plate count. The examples are below.

https://www.nature.com/articles/s41598-018-24916-9 (Links to an external site.)

http://opencfu.sourceforge.net/ (Links to an external site.)

 
**Student B (This example answer will get 5 points)**
I found one software, especially for seed size and color. GrainScan software was designed for seed size and color estimation.  GrainScan uses a grayscale image is derived from the scanned color image by converting Red and Green color channel averaging. Based on the grayscaled image, the dimension measurements will be provided which include area, perimeter, and surrogates for length and width the major and minor axes of the best fit ellipse. Another great point of this software will provide color measurements for each seed in CIELAB values based on user provide color calibration options.

https://plantmethods.biomedcentral.com/articles/10.1186/1746-4811-10-23 (Links to an external site.)

**Student C (This example answer will get 5 points)**
I don't know how to code python, but there are several image analysis packages such as scikit-image (Links to an external site.)

Base on scikit software, you can calculate circularity with "4 * pi * props.area / props.perimeter ** 2" 

The props area can be calculated with number of pixels from the centroid approach. The axis location and length can be converted by using orientation value from props and axis value can be estimated cos(orientation) * length/ 2 and sin (orientation) * length/ 2.



### WORKTING PATH
```bash
mkdir cd /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/
cd /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/
```



### Conda environment
```bash
CONDA_INSTRUMENTATION_ENABLED=1 conda create -n BCH709 python=3.7

conda activate BCH709

CONDA_INSTRUMENTATION_ENABLED=1 conda install -y -c bioconda -c conda-forge  sra-tools bowtie2 minimap2 trinity star multiqc=1.9 samtools=1.9 trim-galore gffread seqkit genometools biopython==1.76 kraken2


CONDA_INSTRUMENTATION_ENABLED=1 conda install -y -c bioconda -c conda-forge -c r openssl=1.0 r-base icu=58.2 bioconductor-ctc  bioconductor-deseq2=1.20.0 bioconductor-biobase=2.40.0  bioconductor-qvalue=2.16.0 r-ape  r-gplots r-fastcluster=1.1.25 libiconv
```

> ## SRA
> Sequence Read Archive (SRA) data, available through multiple cloud providers and NCBI servers, is the largest publicly available repository of high throughput sequencing data. The archive accepts data from all branches of life as well as metagenomic and environmental surveys.
> 
> Searching the SRA: Searching the SRA can be complicated. Often a paper or reference will specify the accession number(s) connected to a dataset. You can search flexibly using a number of terms (such as the organism name) or the filters (e.g. DNA vs. RNA). The SRA Help Manual provides several useful explanations. It is important to know is that projects are organized and related at several levels, and some important terms include:
> 
> Bioproject: A BioProject is a collection of biological data related to a single initiative, originating from a single organization or from a consortium of coordinating organizations; see for example Bio Project 272719
> Bio Sample: A description of the source materials for a project
> Run: These are the actual sequencing runs (usually starting with SRR); see for example SRR1761506

{: .prereq}


> ### Publication (Arabidopsis)
> 
> [Kim JS et al., "ROS1-Dependent DNA Demethylation Is Required for ABA-Inducible NIC3 Expression.", Plant Physiol, 2019 Apr;179(4):1810-1821](http://www.plantphysiol.org/content/179/4/1810)
> 
{: .callout}


### SRA Bioproject site 

```bash
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA272719
```


### Runinfo



| Run        | ReleaseDate     | LoadDate        | spots    | bases      | spots_with_mates | avgLength | size_MB | AssemblyName | download_path                                                                           | Experiment | LibraryName | LibraryStrategy | LibrarySelection | LibrarySource  | LibraryLayout | InsertSize | InsertDev | Platform | Model               | SRAStudy  | BioProject  | Study_Pubmed_id | ProjectID | Sample    | BioSample    | SampleType | TaxID | ScientificName       | SampleName | g1k_pop_code | source | g1k_analysis_group | Subject_ID | Sex | Disease | Tumor | Affection_Status | Analyte_Type | Histological_Type | Body_Site | CenterName | Submission | dbgap_study_accession | Consent | RunHash                          | ReadHash                         |
|------------|-----------------|-----------------|----------|------------|------------------|-----------|---------|--------------|-----------------------------------------------------------------------------------------|------------|-------------|-----------------|------------------|----------------|---------------|------------|-----------|----------|---------------------|-----------|-------------|-----------------|-----------|-----------|--------------|------------|-------|----------------------|------------|--------------|--------|--------------------|------------|-----|---------|-------|------------------|--------------|-------------------|-----------|------------|------------|-----------------------|---------|----------------------------------|----------------------------------|
| SRR1761506 | 1/15/2016 15:51 | 1/15/2015 12:43 | 7379945  | 1490748890 | 7379945          | 202       | 899     |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761506/SRR1761506.1 | SRX844600  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820503 | SAMN03285048 | simple     | 3702  | Arabidopsis thaliana | GSM1585887 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | F335FB96DDD730AC6D3AE4F6683BF234 | 12818EB5275BCB7BCB815E147BFD0619 |
| SRR1761507 | 1/15/2016 15:51 | 1/15/2015 12:43 | 9182965  | 1854958930 | 9182965          | 202       | 1123    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761507/SRR1761507.1 | SRX844601  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820504 | SAMN03285045 | simple     | 3702  | Arabidopsis thaliana | GSM1585888 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 00FD62759BF7BBAEF123BF5960B2A616 | A61DCD3B96AB0796AB5E969F24F81B76 |
| SRR1761508 | 1/15/2016 15:51 | 1/15/2015 12:47 | 19060611 | 3850243422 | 19060611         | 202       | 2324    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761508/SRR1761508.1 | SRX844602  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820505 | SAMN03285046 | simple     | 3702  | Arabidopsis thaliana | GSM1585889 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | B75A3E64E88B1900102264522D2281CB | 657987ABC8043768E99BD82947608CAC |
| SRR1761509 | 1/15/2016 15:51 | 1/15/2015 12:51 | 16555739 | 3344259278 | 16555739         | 202       | 2016    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761509/SRR1761509.1 | SRX844603  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820506 | SAMN03285049 | simple     | 3702  | Arabidopsis thaliana | GSM1585890 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 27CA2B82B69EEF56EAF53D3F464EEB7B | 2B56CA09F3655F4BBB412FD2EE8D956C |
| SRR1761510 | 1/15/2016 15:51 | 1/15/2015 12:46 | 12700942 | 2565590284 | 12700942         | 202       | 1552    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761510/SRR1761510.1 | SRX844604  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820508 | SAMN03285050 | simple     | 3702  | Arabidopsis thaliana | GSM1585891 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | D3901795C7ED74B8850480132F4688DA | 476A9484DCFCF9FFFDAADAAF4CE5D0EA |
| SRR1761511 | 1/15/2016 15:51 | 1/15/2015 12:44 | 13353992 | 2697506384 | 13353992         | 202       | 1639    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761511/SRR1761511.1 | SRX844605  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820507 | SAMN03285047 | simple     | 3702  | Arabidopsis thaliana | GSM1585892 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 5078379601081319FCBF67C7465C404A | E3B4195AFEA115ACDA6DEF6E4AA7D8DF |
| SRR1761512 | 1/15/2016 15:51 | 1/15/2015 12:44 | 8134575  | 1643184150 | 8134575          | 202       | 1067    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761512/SRR1761512.1 | SRX844606  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820509 | SAMN03285051 | simple     | 3702  | Arabidopsis thaliana | GSM1585893 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | DDB8F763B71B1E29CC9C1F4C53D88D07 | 8F31604D3A4120A50B2E49329A786FA6 |
| SRR1761513 | 1/15/2016 15:51 | 1/15/2015 12:43 | 7333641  | 1481395482 | 7333641          | 202       | 960     |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761513/SRR1761513.1 | SRX844607  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820510 | SAMN03285053 | simple     | 3702  | Arabidopsis thaliana | GSM1585894 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 4068AE245EB0A81DFF02889D35864AF2 | 8E05C4BC316FBDFEBAA3099C54E7517B |
| SRR1761514 | 1/15/2016 15:51 | 1/15/2015 12:44 | 6160111  | 1244342422 | 6160111          | 202       | 807     |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761514/SRR1761514.1 | SRX844608  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820511 | SAMN03285059 | simple     | 3702  | Arabidopsis thaliana | GSM1585895 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 0A1F3E9192E7F9F4B3758B1CE514D264 | 81BFDB94C797624B34AFFEB554CE4D98 |
| SRR1761515 | 1/15/2016 15:51 | 1/15/2015 12:44 | 7988876  | 1613752952 | 7988876          | 202       | 1048    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761515/SRR1761515.1 | SRX844609  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820512 | SAMN03285054 | simple     | 3702  | Arabidopsis thaliana | GSM1585896 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 39B37A0BD484C736616C5B0A45194525 | 85B031D74DF90AD1815AA1BBBF1F12BD |
| SRR1761516 | 1/15/2016 15:51 | 1/15/2015 12:44 | 8770090  | 1771558180 | 8770090          | 202       | 1152    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761516/SRR1761516.1 | SRX844610  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820514 | SAMN03285055 | simple     | 3702  | Arabidopsis thaliana | GSM1585897 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | E4728DFBF0F9F04B89A5B041FA570EB3 | B96545CB9C4C3EE1C9F1E8B3D4CE9D24 |
| SRR1761517 | 1/15/2016 15:51 | 1/15/2015 12:44 | 8229157  | 1662289714 | 8229157          | 202       | 1075    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761517/SRR1761517.1 | SRX844611  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820513 | SAMN03285058 | simple     | 3702  | Arabidopsis thaliana | GSM1585898 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | C05BC519960B075038834458514473EB | 4EF7877FC59FF5214DBF2E2FE36D67C5 |
| SRR1761518 | 1/15/2016 15:51 | 1/15/2015 12:44 | 8760931  | 1769708062 | 8760931          | 202       | 1072    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761518/SRR1761518.1 | SRX844612  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820515 | SAMN03285052 | simple     | 3702  | Arabidopsis thaliana | GSM1585899 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 7D8333182062545CECD5308A222FF506 | 382F586C4BF74E474D8F9282E36BE4EC |
| SRR1761519 | 1/15/2016 15:51 | 1/15/2015 12:44 | 6643107  | 1341907614 | 6643107          | 202       | 811     |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761519/SRR1761519.1 | SRX844613  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820516 | SAMN03285056 | simple     | 3702  | Arabidopsis thaliana | GSM1585900 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 163BD8073D7E128D8AD1B253A722DD08 | DFBCC891EB5FA97490E32935E54C9E14 |
| SRR1761520 | 1/15/2016 15:51 | 1/15/2015 12:44 | 8506472  | 1718307344 | 8506472          | 202       | 1040    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761520/SRR1761520.1 | SRX844614  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820517 | SAMN03285062 | simple     | 3702  | Arabidopsis thaliana | GSM1585901 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 791BD0D8840AA5F1D74E396668638DA1 | AF4694425D34F84095F6CFD6F4A09936 |
| SRR1761521 | 1/15/2016 15:51 | 1/15/2015 12:46 | 13166085 | 2659549170 | 13166085         | 202       | 1609    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761521/SRR1761521.1 | SRX844615  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820518 | SAMN03285057 | simple     | 3702  | Arabidopsis thaliana | GSM1585902 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 47C40480E9B7DB62B4BEE0F2193D16B3 | 1443C58A943C07D3275AB12DC31644A9 |
| SRR1761522 | 1/15/2016 15:51 | 1/15/2015 12:49 | 9496483  | 1918289566 | 9496483          | 202       | 1162    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761522/SRR1761522.1 | SRX844616  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820519 | SAMN03285061 | simple     | 3702  | Arabidopsis thaliana | GSM1585903 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | BB05DF11E1F95427530D69DB5E0FA667 | 7706862FB2DF957E4041D2064A691CF6 |
| SRR1761523 | 1/15/2016 15:51 | 1/15/2015 12:46 | 14999315 | 3029861630 | 14999315         | 202       | 1832    |              | https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1761523/SRR1761523.1 | SRX844617  |             | RNA-Seq         | cDNA             | TRANSCRIPTOMIC | PAIRED        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP052302 | PRJNA272719 | 3               | 272719    | SRS820520 | SAMN03285060 | simple     | 3702  | Arabidopsis thaliana | GSM1585904 |              |        |                    |            |     |         | no    |                  |              |                   |           | GEO        | SRA232612  |                       | public  | 101D3A151E632224C09A702BD2F59CF5 | 0AC99FAA6B8941F89FFCBB8B1910696E |


### Subset of data

| Sample information | Run        |
|--------------------|------------|
| WT_rep1            | SRR1761506 |
| WT_rep2            | SRR1761507 |
| WT_rep3            | SRR1761508 |
| ABA_rep1           | SRR1761509 |
| ABA_rep2           | SRR1761510 |
| ABA_rep3           | SRR1761511 |


```bash
mkdir ATH && cd ATH

mkdir raw_data

mkdir trim
```

### fastq-dump submission
```bash
#!/bin/bash
#SBATCH --job-name=fastqdump_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o fastq-dump.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1

fastq-dump SRR1761506 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761507 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761508 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761509 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761510 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761511 --split-3 --outdir ./raw_data  --gzip
```

## Trim-galore
```bash
#!/bin/bash
#SBATCH --job-name=trim_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1
#SBATCH --Dependency=afterok:<PREVIOUS_JOBID>

trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim raw_data/SRR1761506_1.fastq.gz raw_data/SRR1761506_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim raw_data/SRR1761507_1.fastq.gz raw_data/SRR1761507_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim raw_data/SRR1761508_1.fastq.gz raw_data/SRR1761508_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim raw_data/SRR1761509_1.fastq.gz raw_data/SRR1761509_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim raw_data/SRR1761510_1.fastq.gz raw_data/SRR1761510_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim raw_data/SRR1761511_1.fastq.gz raw_data/SRR1761511_2.fastq.gz --fastqc
```

### Reference downloads
https://www.arabidopsis.org/

https://phytozome.jgi.doe.gov/pz/portal.html

```bash

cd /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH

mkdir bam

mkdir reference && cd reference

wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff

wget ftp://ftp.arabidopsis.org/home/tair/Sequences/ATH_cDNA_EST_sequences_FASTA/ATH_cDNA_sequences_20101108.fas

seqkit stats TAIR10_chr_all.fas 

seqkit stats ATH_cDNA_sequences_20101108.fas

```
## Reference index
```bash
#!/bin/bash
#SBATCH --job-name=reference_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1

gffread  TAIR10_GFF3_genes.gff -T -o TAIR10_GFF3_genes.gtf

STAR  --runThreadN <YOUR THREAD> --runMode genomeGenerate --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --genomeFastaFiles  /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/TAIR10_chr_all.fas --sjdbGTFfile  /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf  --sjdbOverhang 99   --genomeSAindexNbases 12
```

> ### ERROR
> Fatal INPUT FILE error, no valid exon lines in the GTF file: /data/gpfs/assoc/bch709-1/wyim/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf
> Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.
{: .prereq} 

```bash
head  /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf
head  /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/TAIR10_chr_all.fas
sed -i 's/Chr//g' TAIR10_GFF3_genes.gtf
```

## Mapping
```bash
#!/bin/bash
#SBATCH --job-name=mapping_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o mapping.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1
#SBATCH --Dependency=afterok:<PREVIOUS_JOBID(trim_ATH)>

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761506_1_val_1.fq.gz /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761506_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/bam/SRR1761506.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761507_1_val_1.fq.gz /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761507_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/bam/SRR1761507.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761508_1_val_1.fq.gz /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761508_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/bam/SRR1761508.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761509_1_val_1.fq.gz /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761509_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/bam/SRR1761509.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761510_1_val_1.fq.gz /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761510_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/bam/SRR1761510.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761511_1_val_1.fq.gz /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/trim/SRR1761511_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/ATH/bam/SRR1761511.bam

```






> ### Publication (Drosophila)
> 
> [Ramond E et al., "Comparative RNA-Seq analyses of Drosophila plasmatocytes reveal gene specific signatures in response to clean injury and septic injury", Plos one, 2019 June 29, 2020](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0235294#sec008)
> 
{: .callout}


### SRA Bioproject site 

```bash
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA638422
```

| Run         | ReleaseDate    | LoadDate       | spots    | bases      | spots_with_mates | avgLength | size_MB | AssemblyName | download_path                                                             | Experiment | LibraryName | LibraryStrategy | LibrarySelection | LibrarySource  | LibraryLayout | InsertSize | InsertDev | Platform | Model               | SRAStudy  | BioProject  | Study_Pubmed_id | ProjectID | Sample     | BioSample    | SampleType | TaxID | ScientificName          | SampleName    | g1k_pop_code | source | g1k_analysis_group | Subject_ID | Sex     | Disease | Tumor | Affection_Status | Analyte_Type | Histological_Type | Body_Site | CenterName                                     | Submission | dbgap_study_accession | Consent | RunHash                          | ReadHash                         |
|-------------|----------------|----------------|----------|------------|------------------|-----------|---------|--------------|---------------------------------------------------------------------------|------------|-------------|-----------------|------------------|----------------|---------------|------------|-----------|----------|---------------------|-----------|-------------|-----------------|-----------|------------|--------------|------------|-------|-------------------------|---------------|--------------|--------|--------------------|------------|---------|---------|-------|------------------|--------------|-------------------|-----------|------------------------------------------------|------------|-----------------------|---------|----------------------------------|----------------------------------|
| SRR11968960 | 6/9/2020 17:10 | 6/9/2020 17:09 | 12256307 | 1237887007 | 0                | 101       | 378     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra23/SRR/011688/SRR11968960 | SRX8512716 | 4w1118-ci   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811237 | SAMN15192434 | simple     | 7227  | Drosophila melanogaster | w1118-ci-rep1 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 3D3A8EBF0A13F90F9305C5DD917E9AE2 | A111523A7FB7106EE54D2D8337D2E8F2 |
| SRR11968959 | 6/9/2020 17:09 | 6/9/2020 17:07 | 14144827 | 1428627527 | 0                | 101       | 432     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra1/SRR/011688/SRR11968959  | SRX8512717 | 5w1118-ci   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811238 | SAMN15192435 | simple     | 7227  | Drosophila melanogaster | w1118-ci-rep2 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 5515CADB5697C29CDC396F942C24F387 | 6D312D3B5BF5001309FF93CB968E584B |
| SRR11968958 | 6/9/2020 17:11 | 6/9/2020 17:09 | 16118803 | 1627999103 | 0                | 101       | 495     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra60/SRR/011688/SRR11968958 | SRX8512718 | 6w1118-ci   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811239 | SAMN15192436 | simple     | 7227  | Drosophila melanogaster | w1118-ci-rep3 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | FCC81714EB524E34632C58BDC1E4C162 | 9F494AF29716E7175EA1E4652B08F0B7 |
| SRR11968957 | 6/9/2020 17:07 | 6/9/2020 17:05 | 6215784  | 627794184  | 0                | 101       | 188     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra47/SRR/011688/SRR11968957 | SRX8512719 | 7w1118-ec   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811240 | SAMN15192437 | simple     | 7227  | Drosophila melanogaster | w1118-ec-rep1 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 9FA82BA9A828F9BDDE839810689EFA4F | CC558BFAAE5EE65BDD1CC1C690575F9D |
| SRR11968956 | 6/9/2020 19:58 | 6/9/2020 19:56 | 46628659 | 4709494559 | 0                | 101       | 1573    |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRR/011688/SRR11968956 | SRX8512720 | 8w1118-ec   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811241 | SAMN15192438 | simple     | 7227  | Drosophila melanogaster | w1118-ec-rep2 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 16E3AE29FAC6BDFDB2B60F5300A02302 | 0F7770A244784C635FEC2DC814A1040C |
| SRR11968955 | 6/9/2020 17:13 | 6/9/2020 17:11 | 16299093 | 1646208393 | 0                | 101       | 496     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra62/SRR/011688/SRR11968955 | SRX8512721 | 9w1118-ec   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811242 | SAMN15192439 | simple     | 7227  | Drosophila melanogaster | w1118-ec-rep3 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | CFA33A602A41E07AC4EFBEED3D2A0FE3 | 4F5983B317885D4E8FFC4B3D312B7674 |
| SRR11968964 | 6/9/2020 17:15 | 6/9/2020 17:12 | 22436848 | 2266121648 | 0                | 101       | 843     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra49/SRR/011688/SRR11968964 | SRX8512712 | 22w1118-l3  | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811233 | SAMN15192443 | simple     | 7227  | Drosophila melanogaster | w1118-l3-rep1 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 2D2BB637C1817EC80B369D1EF0B39615 | 2136B5CFE75B7833A2A6927CF26E701E |
| SRR11968963 | 6/9/2020 19:33 | 6/9/2020 17:14 | 19826612 | 2002487812 | 0                | 101       | 740     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra45/SRR/011688/SRR11968963 | SRX8512713 | 23w1118-l3  | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811234 | SAMN15192444 | simple     | 7227  | Drosophila melanogaster | w1118-l3-rep2 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | A832FE389916D06C16C6F21DB93AD77A | CE81D58F649BBBCECE551891816145AF |
| SRR11968962 | 6/9/2020 17:15 | 6/9/2020 17:12 | 20056763 | 2025733063 | 0                | 101       | 750     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra11/SRR/011688/SRR11968962 | SRX8512714 | 24w1118-l3  | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811235 | SAMN15192445 | simple     | 7227  | Drosophila melanogaster | w1118-l3-rep3 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 3F651F739352EAC0B28096237F2254EC | 98BD1486600E785F9B5F8AC7DBCD4EA6 |
| SRR11968954 | 6/9/2020 17:13 | 6/9/2020 17:10 | 16301608 | 1646462408 | 0                | 101       | 499     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra20/SRR/011688/SRR11968954 | SRX8512722 | 10w1118-sa  | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811243 | SAMN15192440 | simple     | 7227  | Drosophila melanogaster | w1118-sa-rep1 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 55D22FA303406FAE40145D8A1E62598B | 3E5BEB8C6FF03B853BA64D10989542E1 |
| SRR11968966 | 6/9/2020 17:10 | 6/9/2020 17:08 | 16076977 | 1623774677 | 0                | 101       | 485     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra50/SRR/011688/SRR11968966 | SRX8512710 | 11w1118-sa  | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811231 | SAMN15192441 | simple     | 7227  | Drosophila melanogaster | w1118-sa-rep2 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 9789FA28D07EBFD979E3DAE45E9D8CDF | 54D9D5C5343EB9C9A7817434F1D4BB8B |
| SRR11968965 | 6/9/2020 17:10 | 6/9/2020 17:08 | 10379871 | 1048366971 | 0                | 101       | 316     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra76/SRR/011688/SRR11968965 | SRX8512711 | 12w1118-sa  | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811232 | SAMN15192442 | simple     | 7227  | Drosophila melanogaster | w1118-sa-rep3 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 908A23B3924A405F6BE6D5362130E7B3 | 8BD0419DF27093A94546D02492F6661C |
| SRR11968968 | 6/9/2020 17:11 | 6/9/2020 17:09 | 16112703 | 1627383003 | 0                | 101       | 494     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra51/SRR/011688/SRR11968968 | SRX8512708 | 1w1118-uc   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811229 | SAMN15192431 | simple     | 7227  | Drosophila melanogaster | w1118-uc-rep1 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 39E07B4F04BC5A14AE664312E4DD5E67 | 27280649B15E4B766D86363C23679BE1 |
| SRR11968967 | 6/9/2020 17:08 | 6/9/2020 17:06 | 9828233  | 992651533  | 0                | 101       | 302     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra46/SRR/011688/SRR11968967 | SRX8512709 | 2w1118-uc   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811230 | SAMN15192432 | simple     | 7227  | Drosophila melanogaster | w1118-uc-rep2 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | C9868075DE213901336D4DD9D22A9B72 | 558FD3E03B55FA60A231537D5E8EE198 |
| SRR11968961 | 6/9/2020 17:15 | 6/9/2020 17:11 | 16343251 | 1650668351 | 0                | 101       | 498     |              | https://sra-download.ncbi.nlm.nih.gov/traces/sra70/SRR/011688/SRR11968961 | SRX8512715 | 3w1118-uc   | RNA-Seq         | Oligo-dT         | TRANSCRIPTOMIC | SINGLE        | 0          | 0         | ILLUMINA | Illumina HiSeq 2500 | SRP266662 | PRJNA638422 |                 | 638422    | SRS6811236 | SAMN15192433 | simple     | 7227  | Drosophila melanogaster | w1118-uc-rep3 |              |        |                    |            | unknown |         | no    |                  |              |                   |           | SWISS FEDERAL INSTITUTE OF TECHNOLOGY LAUSANNE | SRA1085163 |                       | public  | 1725A91FA94755464378D8FF0F18A197 | 870D1C8B738F5A31C179B44124757B27 |



IMAGE


Fig 2. Transcriptome summaries from unchallenged whole larvae and hemocytes from unchallenged and infected larvae.
(A) Transcriptome summary showing the number of reads for each triplicate in all experimental conditions with their corresponding number of mapped reads and the average percentage of alignment to the D. melanogaster genome. (B) Venn diagram representing the quantity of shared genes between all experimental treatments: Unchallenged wandering L3 larvae, hemocytes from unchallenged larvae, hemocytes from clean-pricked larvae (CI), hemocytes from larvae pricked with Escherichia coli (Ec), hemocytes from larvae pricked with Staphylococcus aureus (Sa). 



### Subset of data


| Sample information  | Run         |
|---------------------|-------------|
| 22w1118-l3          | SRR11968964 |
| 23w1118-l3          | SRR11968963 |
| 24w1118-l3          | SRR11968962 |
| 10w1118-sa          | SRR11968954 |
| 11w1118-sa          | SRR11968966 |
| 12w1118-sa          | SRR11968965 |



```bash
cd /data/gpfs/assoc/bch709-1/wyim/RNA-Seq_example/
mkdir Drosophila && cd Drosophila

mkdir raw_data

mkdir trim
```

### fastq-dump submission
```bash
#!/bin/bash
#SBATCH --job-name=fastqdump_Droso
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o fastq-dump.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1

fastq-dump  SRR11968964 --outdir ./raw_data  --gzip
fastq-dump  SRR11968963 --outdir ./raw_data  --gzip
fastq-dump  SRR11968962 --outdir ./raw_data  --gzip
fastq-dump  SRR11968954 --outdir ./raw_data  --gzip
fastq-dump  SRR11968966 --outdir ./raw_data  --gzip
fastq-dump  SRR11968965 --outdir ./raw_data  --gzip
```

## Trim-galore
```bash
#!/bin/bash
#SBATCH --job-name=trim_Droso
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1
#SBATCH --Dependency=afterok:<PREVIOUS_JOBID>

trim_galore --cores 2  --max_n 40  --gzip -o trim raw_data/SRR11968964.fastq.gz --fastqc
trim_galore --cores 2  --max_n 40  --gzip -o trim raw_data/SRR11968963.fastq.gz --fastqc
trim_galore --cores 2  --max_n 40  --gzip -o trim raw_data/SRR11968962.fastq.gz --fastqc
trim_galore --cores 2  --max_n 40  --gzip -o trim raw_data/SRR11968954.fastq.gz --fastqc
trim_galore --cores 2  --max_n 40  --gzip -o trim raw_data/SRR11968966.fastq.gz --fastqc
trim_galore --cores 2  --max_n 40  --gzip -o trim raw_data/SRR11968965.fastq.gz --fastqc
```



### Reference downloads
https://flybase.org/


```bash

cd /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila

mkdir bam

mkdir reference && cd reference

wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/fasta/dmel-all-chromosome-r6.36.fasta.gz
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/gtf/dmel-all-r6.36.gtf.gz
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/fasta/dmel-all-CDS-r6.36.fasta.gz

gunzip dmel-all-chromosome-r6.36.fasta.gz
gunzip dmel-all-r6.36.gtf.gz
gunzip dmel-all-CDS-r6.36.fasta.gz

seqkit stats dmel-all-chromosome-r6.36.fasta
seqkit stats dmel-all-CDS-r6.36.fasta  
```
## Reference index
```bash
#!/bin/bash
#SBATCH --job-name=reference_Droso
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1

 STAR  --runThreadN 24 --runMode genomeGenerate --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --genomeFastaFiles  /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/dmel-all-chromosome-r6.36.fasta --sjdbGTFfile /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/dmel-all-r6.36.gtf  --sjdbOverhang 99 -genomeSAindexNbases 12

```
## Mapping reads
```bash
#!/bin/bash
#SBATCH --job-name=mapping_Droso
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH -p cpu-s2-core-0 
#SBATCH -A cpu-s2-bch709-1
#SBATCH --Dependency=afterok:<trim_Droso>

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/trim/SRR11968964_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/bam/SRR11968964.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/trim/SRR11968963_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/bam/SRR11968963.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/trim/SRR11968962_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/bam/SRR11968962.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/trim/SRR11968954_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/bam/SRR11968954.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/trim/SRR11968966_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/bam/SRR11968966.bam
STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/reference/ --readFilesIn /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/trim/SRR11968965_trimmed.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/bam/SRR11968965.bam
```


## Assignment
```bash
cd /data/gpfs/assoc/bch709-1/wyim/RNA-Seq_example
multiqc . -n rnaseq1
```
Please upload rnaseq1.html to Webcampus.



### Investigate taxa

Here we introduce a software called Kraken2. This tool uses k-mers to assign a taxonomic labels in form of NCBI Taxonomy to the sequence (if possible). The taxonomic label is assigned based on similar k-mer content of the sequence in question to the k-mer content of reference genome sequence. The result is a classification of the sequence in question to the most likely taxonomic label. If the k-mer content is not similar to any genomic sequence in the database used, it will not assign any taxonomic label.


### Donwload most recent database
```bash
https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20200919.tar.gz
```

```bash
cd /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example

mkdir Kraken2 && cd Kraken2
kraken2-inspect --db EXAMPLE_DB | head -5

kraken2  --threads 24 --report SRR11968954 --db /data/gpfs/assoc/bch709-1/Course_material/database/ /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/raw_data/SRR11968954.fastq.gz
```
```
"C"/"U": a one letter code indicating that the sequence was either classified or unclassified.
```
```bash
https://fbreitwieser.shinyapps.io/pavian/
```
