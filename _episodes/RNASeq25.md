---
layout: page
title: RNA-SEq
published: true
---

## HPC Clusters

This exercise focuses on using HPC clusters for large-scale data analysis (e.g., Next Generation Sequencing, genome annotation, evolutionary studies). These clusters contain multiple processors with large amounts of RAM, making them ideal for computationally intensive tasks. The operating system is primarily UNIX, accessed via the command line. All the commands you've learned in previous exercises can be used here.

Pronghorn High Performance Computing offers shared infrastructure for researchers and students at UNR. You can find available resources [here](https://www.unr.edu/research-computing/hpc). Request access through your department or advisor. All attendees of this workshop will have their accounts set up on the HPC class education cluster using their UNR NetID and password. You should have received a confirmation email with connection instructions. This exercise covers connecting to a remote HPC server, transferring files, and running programs by requesting resources.

To log into the HPC front-end/job-submission system (pronghorn.rc.unr.edu), use your UNR NetID and password. Windows users will need an SSH client, while Mac/Linux users have SSH built-in.

```bash
ssh <YOUR_NET_ID>@pronghorn.rc.unr.edu
## First login will prompt for key confirmation. Choose 'yes.'
```

### Pronghorn HPC Cluster Overview

Pronghorn is UNR's new GPU-accelerated HPC cluster, supporting general research across NSHE. Comprising CPU, GPU, and storage subsystems, Pronghorn's main features include:

- **CPU Partition**: 93 nodes, 2,976 CPU cores, 21TiB of memory.
- **GPU Partition**: 44 NVIDIA Tesla P100 GPUs, 352 CPU cores, 2.75TiB of memory.
- **Storage**: 1PB high-performance storage using IBM SpectrumScale.

Pronghorn is located at Switch Citadel Campus, 25 miles east of UNR. Switch is renowned for sustainable data center operations.

![Pronghorn system map](../fig/pronghorn.png){: width="70%" height="70%"}


## Customizing the Bash Prompt for Pronghorn
```bash
echo '###BCH709' >> ~/.bashrc
echo 'tty -s && export PS1="\[\033[38;5;164m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;231m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;172m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bashrc
echo "alias ls='ls --color=auto'" >> ~/.bashrc
source ~/.bashrc
```


## File Transfer Methods

Several methods exist for transferring files to/from HPC clusters, including command-line tools (scp, rsync) and graphical clients (SCP, SFTP). For secure copying:

### Transferring from Local to Remote System
```bash
scp <source_file> <username>@pronghorn.rc.unr.edu:<target_location>
```
Example:
```bash
mkdir ~/bch709
cd ~/bch709
echo "hello world" > test_uploading_file.txt
scp test_uploading_file.txt <username>@pronghorn.rc.unr.edu:~/
```

### Transferring from Remote to Local System
```bash
scp <username>@pronghorn.rc.unr.edu:<source_file> <destination_file>
```
Example:
```bash
scp <username>@pronghorn.rc.unr.edu:~/test_downloading_file.txt ~/
```

### Recursive Directory Transfer (Local to Remote)
```bash
scp -r <source_directory> <username>@pronghorn.rc.unr.edu:<target_directory>
```
Example:
```bash
scp -r ../bch709 <username>@pronghorn.rc.unr.edu:~/
```

### Opening Location

- **Windows (WSL)**
```bash
cd ~/bch709
explorer.exe .
```

- **Mac**
```bash
cd ~/bch709
open .
```

### Rsync Usage

- Sync local directory to remote server:
```bash
rsync -avhP <source_directory> <username>@pronghorn.rc.unr.edu:<target_directory>
```

- Sync remote directory to local:
```bash
rsync -avhP <username>@pronghorn.rc.unr.edu:<source_directory> <target_directory>
```

## Conda Installation on Pronghorn

Install Miniconda3 using the following commands:
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
~/miniconda3/bin/conda init
```
Follow the on-screen instructions to complete the installation.



## Using Conda Environment

```bash
conda install mamba

mamba create -y -n RNASEQ_bch709 -c bioconda -c conda-forge sra-tools=3.1.1 minimap2 star trim-galore gffread seqkit samtools multiqc subread tree
conda activate RNASEQ_bch709
```


### Environment create and installation in Pronghorn

```bash
# Create a new conda environment named "rnaseq".
# Add two channels to fetch the required packages:
# - bioconda: A channel specializing in bioinformatics software
# - conda-forge: A community-maintained collection of conda packages
-c bioconda -c conda-forge 

# List of packages/software to be installed in the "rnaseq" environment:
# - fastqc: A tool for quality control checks on raw sequence data
# - trim-galore: A wrapper tool around Cutadapt and FastQC to consistently apply adapter and quality trimming
# - hisat2: A fast and sensitive alignment program for mapping next-generation sequencing reads to a population of genomes
# - samtools: A suite of programs for interacting with high-throughput sequencing data
# - subread: A toolkit for processing next-gen sequencing read data, including feature counting
# - bioconductor-deseq2: A package for differential expression analysis based on the negative binomial distribution
# - bc: An arbitrary precision calculator language
# The "-y" flag allows the command to proceed without asking for user confirmation.



## Connecting Scratch Disk

```bash
mkdir /data/gpfs/assoc/bch709-5/students/$(whoami)
cd /data/gpfs/assoc/bch709-5/students/$(whoami)
ln -s /data/gpfs/assoc/bch709-5/students/$(whoami) ~/scratch
```

## Job Submission with SBATCH

Create a job submission script named `submit.sh`:

```bash
nano submit.sh
```

Add the following content:

```bash
#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --ntasks=1
#SBATCH --mem=1g
#SBATCH --time=8:10:00
#SBATCH --account=cpu-s5-bch709-5
#SBATCH --partition=cpu-core-0

for i in {1..1000}; do
  echo $i;
  sleep 1;
done
```

Submit the job:
```bash
chmod 775 submit.sh
sbatch submit.sh
```

To cancel the job:
```bash
scancel <JOB_ID>
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

### SRA Data Access

SRA (Sequence Read Archive) is a repository of high-throughput sequencing data. To download sequencing data from SRA, use `fastq-dump`.

### preparing Fastq-Dump Job

Create and submit a script `fastq-dump.sh` to download RNA-Seq data:

```bash
mkdir ~/scratch/raw_data
cd ~/scratch/raw_data
nano fastq-dump.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=fastqdump_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o fastq-dump.out
#SBATCH --account=cpu-s5-bch709-5
#SBATCH --partition=cpu-core-0

fastq-dump SRR1761506 --split-3 --outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data --gzip
fastq-dump SRR1761507 --split-3 --outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data --gzip
fastq-dump SRR1761508 --split-3 --outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data --gzip
fastq-dump SRR1761509 --split-3 --outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data --gzip
fastq-dump SRR1761510 --split-3 --outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data --gzip
fastq-dump SRR1761511 --split-3 --outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data --gzip
```

### submit Fastq-Dump Job

```bash
sbatch fastq-dump.sh
squeue
```
### Explain the command
```
This command uses the `fastq-dump` tool from the SRA (Sequence Read Archive) Toolkit to download and process sequencing data associated with the accession number **SRR1761510**. Let’s break it down:

1. **`fastq-dump SRR1761510`**: This part tells `fastq-dump` to retrieve the sequencing data associated with the accession **SRR1761510**. Accession numbers like these correspond to specific datasets available in the NCBI SRA database.

2. **`--split-3`**: This option splits paired-end reads into separate files (for example, `_1` for the first read in a pair and `_2` for the second read). If there are also unpaired reads, they’ll be stored in a separate file as well.

3. **`--outdir /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data`**: The `--outdir` option specifies the output directory for the downloaded FASTQ files. The path:
   
   ```
   /data/gpfs/assoc/bch709-5/students/$(whoami)/raw_data
   ```
   
   includes `$(whoami)`, which dynamically inserts the current username. This makes sure each user has their own dedicated output directory under `raw_data`.

4. **`--gzip`**: This option compresses the output FASTQ files in **GZIP** format, saving storage space and making downstream data handling more efficient.

In summary, this command will download the sequencing data for **SRR1761510** from SRA, split the paired-end reads, and save them in a specific directory with GZIP compression. 



```

### Trimming Reads with Trim-Galore

Submit a trimming job:

```bash
mkdir  -p ~/scratch/RNA-Seq_example/ATH/trim
cd  ~/scratch/RNA-Seq_example/ATH
nano trim.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=trim_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o trim.out
#SBATCH --account=cpu-s5-bch709-5
#SBATCH --partition=cpu-core-0

trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o  ~/scratch/RNA-Seq_example/ATH/trim --basename SRR1761506 raw_data/SRR1761506_1.fastq.gz raw_data/SRR1761506_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o  ~/scratch/RNA-Seq_example/ATH/trim --basename SRR1761507  raw_data/SRR1761507_1.fastq.gz raw_data/SRR1761507_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o  ~/scratch/RNA-Seq_example/ATH/trim --basename SRR1761508 raw_data/SRR1761508_1.fastq.gz raw_data/SRR1761508_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o  ~/scratch/RNA-Seq_example/ATH/trim --basename SRR1761509 raw_data/SRR1761509_1.fastq.gz raw_data/SRR1761509_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o  ~/scratch/RNA-Seq_example/ATH/trim --basename SRR1761510 raw_data/SRR1761510_1.fastq.gz raw_data/SRR1761510_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o  ~/scratch/RNA-Seq_example/ATH/trim --basename SRR1761511 raw_data/SRR1761511_1.fastq.gz raw_data/SRR1761511_2.fastq.gz --fastqc
```
### Submit trimming
```bash
squeue -u $(whoami)
sbatch --dependency=afterok:######## trim.sh 
```

### Explain the command
```
This command runs **Trim Galore**, a tool used for trimming adapters and low-quality sequences from high-throughput sequencing data, with specific parameters for paired-end RNA-Seq reads. Here’s a detailed breakdown:

1. **`trim_galore --paired`**: This specifies that the input data consists of paired-end reads. Trim Galore will process both forward and reverse reads together, ensuring paired-end compatibility after trimming.

2. **`--three_prime_clip_R1 5` and `--three_prime_clip_R2 5`**: These options clip (remove) 5 bases from the 3' end of both reads in the pair. **`_R1`** applies to the first read and **`_R2`** applies to the second read, which can help remove low-quality or unwanted bases at the end of each read.

3. **`--cores 2`**: This specifies the number of cores (CPUs) to use, allowing Trim Galore to run with 2 parallel threads for faster processing.

4. **`--max_n 40`**: This sets the maximum number of ambiguous bases (`N`) allowed per read. Reads with more than 40 `N` bases will be discarded, helping improve the quality of the data.

5. **`--gzip`**: This option compresses the output files in **GZIP** format, saving space and making downstream handling more efficient.

6. **`-o ~/scratch/RNA-Seq_example/ATH/trim`**: This specifies the output directory for the processed files. Here, the results will be saved in `~/scratch/RNA-Seq_example/ATH/trim`.

7. **`--basename SRR1761511`**: This sets the base name for the output files, which will start with **SRR1761511**. This is useful for organizing results by sample name.

8. **`raw_data/SRR1761511_1.fastq.gz raw_data/SRR1761511_2.fastq.gz`**: These are the input files: the paired-end FASTQ files (forward and reverse) that will be trimmed.

9. **`--fastqc`**: This tells Trim Galore to run **FastQC** on the trimmed reads, generating a quality control report. FastQC provides information on read quality, adapter content, and other metrics.

### In summary:

This command will trim adapters and low-quality bases from paired-end RNA-Seq reads in **SRR1761511**, remove excess ambiguous bases, clip 5 bases from the 3' ends of both reads, and save the results (compressed) to the specified output directory with **FastQC** quality reports. This helps prepare high-quality reads for further analysis.
```
### Trim the reads
- Trim IF necessary
   - Synthetic bases can be an issue for SNP calling
   - Insert size distribution may be more important for assemblers
- Trim/Clip/Filter reads
- Remove adapter sequences
- Trim reads by quality
- Sliding window trimming
- Filter by min/max read length
- Remove reads less than ~18nt
- Demultiplexing/Splitting

![Trimming]({{{site.baseurl}}/fig/trim.png)  

[Cutadapt](https://github.com/marcelm/cutadapt/)  
[fastp](https://github.com/OpenGene/fastp)  
[Skewer](https://github.com/relipmoc/skewer)  
[Prinseq](http://prinseq.sourceforge.net/)  
[Trimmomatics](http://www.usadellab.org/cms/?page=trimmomatic)  
[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)  


## Download reference file
```
cd ~/scratch/RNA-Seq_example/ATH

mkdir bam
mkdir reference
cd reference

wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz -O TAIR10_chr_all.fas.gz --no-check-certificate

wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff  -O TAIR10_GFF3_genes.gff --no-check-certificate

seqkit stats TAIR10_chr_all.fas.gz
```

### Explain the command
```
This command uses **wget** to download a GFF3 file from the Arabidopsis Information Resource (TAIR) website. Here’s what each part does:

1. **`wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff`**: This is the download command with the full URL of the TAIR10 GFF3 file. **GFF3** (General Feature Format, version 3) files are used for genome annotations and contain information about gene locations, exons, coding regions, and other genomic features.

2. **`-O TAIR10_GFF3_genes.gff`**: This option specifies the output filename. Without this option, `wget` would use the default name provided by the URL, which might be complex or include additional characters. **TAIR10_GFF3_genes.gff** is set as the filename to make it easier to reference.

3. **`--no-check-certificate`**: This option tells `wget` to ignore SSL certificate verification errors. It can be useful when a website’s SSL certificate is expired or not properly recognized. Here, it ensures the download proceeds without interruptions.

### Summary:

This command downloads the **TAIR10_GFF3_genes.gff** file from the TAIR database without SSL verification, saving it as **TAIR10_GFF3_genes.gff**. This GFF3 file will contain essential genome annotation data for Arabidopsis, commonly used in genomic analysis.
```

## Convert GFF to GTF
```bash
cd ~/scratch/RNA-Seq_example/ATH/reference
gffread TAIR10_GFF3_genes.gff -T -F --keep-exon-attrs -o TAIR10_GFF3_genes.gtf
```
### Explain the command
```
This command uses **gffread** to convert a **GFF3** file into a **GTF** file format. **GTF** (Gene Transfer Format) is similar to **GFF** but often preferred in certain bioinformatics tools, especially for transcript-based analyses. Here’s a breakdown of each component:

1. **`gffread TAIR10_GFF3_genes.gff`**: This specifies the input file in **GFF3** format, here **TAIR10_GFF3_genes.gff**.

2. **`-T`**: This option tells `gffread` to output the file in **GTF** format instead of **GFF3**.

3. **`-F`**: This forces the inclusion of features that might be incomplete or missing certain attributes (e.g., lacking start or stop codons). It ensures that all exons are included in the output.

4. **`--keep-exon-attrs`**: This option retains additional attributes associated with exon features, which are often discarded in standard conversions. Keeping exon attributes can be valuable for certain analyses where more detailed annotation is needed.

5. **`-o TAIR10_GFF3_genes.gtf`**: This specifies the output filename. Here, the GTF-formatted file will be saved as **TAIR10_GFF3_genes.gtf**.

### Summary:

This command converts the **TAIR10** genome annotation file from **GFF3** to **GTF** format, ensuring that all exons are retained along with their attributes, even if they lack some standard annotations. This output file, **TAIR10_GFF3_genes.gtf**, will be useful for downstream transcript-based analysis in pipelines or tools that prefer GTF format.
```


## Create reference index
```bash
cd  ~/scratch/RNA-Seq_example/ATH/reference
ls -algh
gunzip TAIR10_chr_all.fas.gz

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
#SBATCH --account=cpu-s5-bch709-5
#SBATCH --partition=cpu-core-0

STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles  ~/scratch/RNA-Seq_example/ATH/reference/TAIR10_chr_all.fas --sjdbGTFfile ~/scratch/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf --sjdbOverhang 99   --genomeSAindexNbases 12
```

### Explain the command
```
This command uses **STAR** (Spliced Transcripts Alignment to a Reference) to generate a genome index, which is a critical step for efficiently aligning RNA-Seq reads to a reference genome. Here’s a detailed breakdown of each parameter:

1. **`STAR`**: This calls the STAR program, which is an RNA-Seq read aligner optimized for high accuracy and speed.

2. **`--runThreadN 48g`**: This option specifies the number of threads (CPUs) STAR should use. However, the argument should be an integer (like `48`) rather than `48g`. Assuming you intended to use 48 threads, the correct syntax would be `--runThreadN 48`.

3. **`--runMode genomeGenerate`**: This tells STAR to run in genome generation mode, which creates an index for the reference genome. This index is needed for subsequent alignment steps.

4. **`--genomeDir .`**: This sets the directory where the generated genome index files will be stored. Using `.` specifies the current directory.

5. **`--genomeFastaFiles ~/scratch/RNA-Seq_example/ATH/reference/TAIR10_chr_all.fas`**: This provides the path to the reference genome FASTA file. Here, **TAIR10_chr_all.fas** contains the reference sequences for **Arabidopsis thaliana**.

6. **`--sjdbGTFfile ~/scratch/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf`**: This specifies the path to the annotation file in **GTF** format (TAIR10_GFF3_genes.gtf). STAR uses this file to incorporate known splice junctions into the index, which improves alignment accuracy, especially for spliced reads in RNA-Seq data.

7. **`--sjdbOverhang 99`**: This defines the length of the sequence to be used for junctions. Ideally, it should be set to the length of the read minus 1. For instance, if the RNA-Seq reads are 100 bp, `sjdbOverhang` should be 99. This value helps STAR optimize the alignment of reads that span splice junctions.

8. **`--genomeSAindexNbases 12`**: This parameter controls the size of the suffix array index used by STAR. **12** is typical for a smaller genome (such as Arabidopsis) to balance between memory usage and indexing speed. Larger values reduce memory requirements but can slow down the alignment step slightly.

### Summary:

This command sets STAR to generate a genome index for **Arabidopsis thaliana** using 48 threads, based on the reference FASTA and GTF files. The generated index will be saved in the current directory and includes splice junction information, which is crucial for accurately mapping RNA-Seq reads that contain introns. This index will be used in future alignment steps to align RNA-Seq reads quickly and accurately to the genome.
```



## Mapping the reads to genome index
```bash
cd  ~/scratch/RNA-Seq_example/ATH/
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
#SBATCH --account=cpu-s5-bch709-5
#SBATCH --partition=cpu-core-0

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_1.fq.gz ~/scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/scratch/RNA-Seq_example/ATH/bam/SRR1761506.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/scratch/RNA-Seq_example/ATH/trim/SRR1761507_val_1.fq.gz ~/scratch/RNA-Seq_example/ATH/trim/SRR1761507_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/scratch/RNA-Seq_example/ATH/bam/SRR1761507.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/scratch/RNA-Seq_example/ATH/trim/SRR1761508_val_1.fq.gz ~/scratch/RNA-Seq_example/ATH/trim/SRR1761508_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/scratch/RNA-Seq_example/ATH/bam/SRR1761508.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/scratch/RNA-Seq_example/ATH/trim/SRR1761509_val_1.fq.gz ~/scratch/RNA-Seq_example/ATH/trim/SRR1761509_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/scratch/RNA-Seq_example/ATH/bam/SRR1761509.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/scratch/RNA-Seq_example/ATH/trim/SRR1761510_val_1.fq.gz ~/scratch/RNA-Seq_example/ATH/trim/SRR1761510_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/scratch/RNA-Seq_example/ATH/bam/SRR1761510.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/scratch/RNA-Seq
```


### Submit mapping
```bash
squeue -u $(whoami)
sbatch --dependency=afterok:########:####### align.sh
```

### Explain the command
```
This command uses **STAR** in **alignment mode** to align paired-end RNA-Seq reads to a pre-built genome index. Here’s a detailed breakdown:

1. **`STAR --runMode alignReads`**: Specifies that STAR should run in **alignment mode**, which aligns RNA-Seq reads to the reference genome.

2. **`--runThreadN 8`**: Sets the number of threads to use (in this case, 8), which will speed up the alignment process by utilizing multiple CPU cores.

3. **`--readFilesCommand zcat`**: Instructs STAR to use `zcat` to decompress the input files since they are **GZIP** compressed (`.fq.gz` format).

4. **`--outFilterMultimapNmax 10`**: Sets the maximum number of loci a read can map to. If a read maps to more than 10 locations, it will be discarded. This helps control the level of ambiguity in alignment, especially in repetitive regions.

5. **`--alignIntronMin 25`**: Specifies the minimum allowed length for introns. STAR will ignore introns shorter than 25 bp, which reduces false alignments in smaller repetitive regions.

6. **`--alignIntronMax 10000`**: Sets the maximum allowed intron length to 10,000 bp, accommodating typical intron lengths found in **Arabidopsis**. This helps STAR avoid spurious alignments that would involve unusually long gaps.

7. **`--genomeDir ~/scratch/RNA-Seq_example/ATH/reference/`**: Specifies the directory containing the STAR genome index, created in the previous genome generation step.

8. **`--readFilesIn ~/scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_1.fq.gz ~/scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_2.fq.gz`**: These are the paths to the input FASTQ files for paired-end reads (forward and reverse), which have been trimmed and compressed.

9. **`--outSAMtype BAM SortedByCoordinate`**: Specifies that the output should be in **BAM** format and sorted by genomic coordinates. BAM is a binary, compressed format for alignment data, commonly used for downstream analysis.

10. **`--outFileNamePrefix ~/scratch/RNA-Seq_example/ATH/bam/SRR1761506.bam`**: Sets the output filename prefix, so the results will be saved with the prefix `SRR1761506.bam` in the specified BAM directory. STAR will automatically append additional information as needed.

### Summary:

This command aligns the trimmed, paired-end RNA-Seq reads for **SRR1761506** to the **Arabidopsis** reference genome, utilizing 8 threads. It decompresses the input files, filters out highly multimapping reads, and limits intron size for optimized mapping. The output is saved in sorted **BAM** format, ready for downstream analysis.
```

## BW algorithm
![algorithm]({{{site.baseurl}}/fig/algorithm.png)
![banana]({{{site.baseurl}}/fig/banana.png)
![banana]({{{site.baseurl}}/fig/BackwardMatching.png)

####################################




# Mouse RNA-Seq
![](https://i.imgur.com/rqNc6OA.png)

https://www.sciencedirect.com/science/article/pii/S2211124722011111

Benraiss A et al., "A TCF7L2-responsive suppression of both homeostatic and compensatory remyelination in Huntington disease mice.", Cell Rep, 2022 Aug 30;40(9):111291


### Working directory (Pronghorn)

```bash
echo $USER

cd /data/gpfs/assoc/bch709-5/students/${USER}

mkdir mouse
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/fastq
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref 
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/bam
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/readcount
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/DEG
```


## Reference Download
![]({{{site.baseurl}}/fig/mouse_ref.png)


https://hgdownload.soe.ucsc.edu/downloads.html


```bash
### change working directory
cd /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref 

### download
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz

### decompress
gunzip mm39.fa.gz
gunzip refGene.gtf.gz


### index
nano index.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=index_mouse
#SBATCH --cpus-per-task=12
#SBATCH --time=2-15:00:00
#SBATCH --mem=48g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o index.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-5
#SBATCH --partition=cpu-core-0

STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles  /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref/mm39.fa  --sjdbGTFfile /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref/refGene.gtf --sjdbOverhang 99   --genomeSAindexNbases 12
```



## FASTQ file 
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/mouse/

### Link file (without copy)
ln -s /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/fastq/* /data/gpfs/assoc/bch709-5/students/${USER}/mouse/fastq

ls /data/gpfs/assoc/bch709-5/students/${USER}/mouse/fastq
```

### Create file list
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/mouse/fastq

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u > /data/gpfs/assoc/bch709-5/students/${USER}/mouse/filelist

cat /data/gpfs/assoc/bch709-5/students/${USER}/mouse/filelist
```
### Regular expression

https://regex101.com/

## Trim reads
```bash
trim_galore --paired  --three_prime_clip_R1 [integer] --three_prime_clip_R2 [integer]  --cores [integer]   --max_n [integer]   --fastqc --gzip -o /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim {READ_R1} {READ_R2}
```

####  Trim reads example
```bash
trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim {READ_R1} {READ_R2}
```

### Prepare templet
```bash
cat /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/run.sh | sed "s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Trim/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g"  > /data/gpfs/assoc/bch709-5/students/${USER}/mouse/fastq/trim.sh
```

### Edit templet
```bash
nano /data/gpfs/assoc/bch709-5/students/${USER}/mouse/fastq/trim.sh
```

### Batch submission
```bash
# Check file list
cat ../filelist


# Loop file list
### Add Forward read to variable
### Add reverse read from forward read name substitution

for i in `cat ../filelist`
    do

    read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        echo $read1 $read2
done


# Loop file list
### add file name from variable to trim-galore

### merge trim-galore command and trim.sh
for i in `cat ../filelist`
    do
        read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        echo $read1 $read2
        echo "trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim $read1 $read2" | cat trim.sh - > ${i}_trim.sh
done
```

```bash
### Batch submission

ls *.sh
ls -1 *.sh

### Loop *.sh printing
for i in `ls -1 *.sh`
do
    echo $i
done

### Loop *.sh submission
for i in `ls -1 *.sh`
do
    sbatch $i
done
```

### Check submission
```bash
squeue -u ${USER}
```


## Environment activation
```bash
conda activate BCH709_RNASeq
```

## Copy files
cp /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/ref/*  /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref/

cp /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/trim/*  /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim/


## RNA-Seq Alignment
```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim

#### Copy templet
cat /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/run.sh | sed "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Trim/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" > /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim/mapping.sh

#### Edit templet
nano mapping.sh
```

### Check output
```bash
ls -algh /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim
```

### Output example
```bash
[FILENAME]_R1_val_1.fq.gz [FILENAME]_R2_val_2.fq.gz
```


### STAR RNA-Seq alignment batch file test

```bash
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim/${read1} /data/gpfs/assoc/bch709-5/students/${USER}/mouse/trim/${read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-5/students/${USER}/mouse/bam/${i}.bam" | cat mapping.sh - > ${i}_mapping.sh
    done
```



###  Job submission dependency on Mapping
```bash

for i in `ls -1 *_mapping.sh`
do
    sbatch  $i 
done

```


## Reads count
In the case of RNA-Seq, the features are typically genes, where each gene is considered here as the union of all its exons. Counting RNA-seq reads is complex because of the need to accommodate exon splicing. The common approach is to summarize counts at the gene level, by counting all reads that overlap any exon for each gene. In this method, gene annotation file from RefSeq or Ensembl is often used for this purpose. So far there are two major feature counting tools: featureCounts (Liao et al.) and htseq-count (Anders et al.)

![featurecount]({{site.baseurl}}/fig/featurecount.png)


#### Copy templet
mkdir  /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/readcount
cd /data/gpfs/assoc/bch709-5/students/wyim/RNA-Seq_example/ATH

cat  /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/run.sh | sed "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Count/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" > /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/count.sh


### FeatureCounts read bam file

```bash
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/readcount/
cd /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/raw_data
ls -1 *.gz | sed 's/_.*//g' | sort -u > /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/filelist
cd  /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/bam
ls -1 *.sortedByCoord.out.bam
ls -1 *.sortedByCoord.out.bam| tr '\n' ' '
cp /data/gpfs/assoc/bch709-5/students/Course_materials/ATH/reference/*  /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf
```

###########

```bash
featureCounts -o [output] -T [threads] -Q 1 -p -M  -g gene_id -a [GTF] [BAMs]
```

### Edit file

```bash
#paste this to count.sh
featureCounts -o /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/readcount/ATH_featucount -T 4 -Q 1 -p -M  -g gene_id -a /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf $(for i in `cat  /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH/filelist`; do echo ${i}.bamAligned.sortedByCoord.out.bam| tr '\n' ' ';done)
```

```bash
nano count.sh
```

## RNA-Seq report

## How to make a report?
![MultiQC]({{{site.baseurl}}/fig/multiqc.png)
[MultiQC](https://multiqc.info/)

### MultiQC
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/RNA-Seq_example/ATH
 multiqc --pdf -n test .
```
*Slurm provides resource management for the processors allocated to a job, so that multiple job steps can be simultaneously submitted and queued until there are available resources within the job's allocation.*
############


## RNA Sequencing
![RNA Sequencing]({{{site.baseurl}}/fig/rnaseq.png)

1. The transcriptome is spatially and temporally dynamic
2. Data comes from functional units (coding regions)
3. Only a tiny fraction of the genome

>## Introduction
>Sequence based assays of transcriptomes (RNA-seq) are in wide use because of their favorable properties for quantification, transcript discovery and splice isoform identification, as well as adaptability for numerous more specialized measurements. RNA-Seq studies present some challenges that are shared with prior methods such as microarrays and SAGE tagging, and they also present new ones that are specific to high-throughput sequencing platforms and the data they produce. This document is part of an ongoing effort to provide the community with standards and guidelines that will be updated as RNASeq matures and to highlight unmet challenges. The intent is to revise this document periodically to capture new advances and increasingly consolidate standards and best practices.
>
>RNA-Seq experiments are diverse in their aims and design goals, currently including multiple types of RNA isolated from whole cells or from specific sub-cellular compartments or biochemical classes, such as total polyA+ RNA, polysomal RNA, nuclear ribosome-depleted RNA, various size fractions of RNA and a host of others. The goals of individual experiments range from major transcriptome “discovery” that seeks to define and quantify all RNA species in a starting RNA sample to experiments that simply need to detect significant changes in the more abundant RNA classes across many samples.  
{: .prereq}



![RNA Sequencing workflow]({{{site.baseurl}}/fig/rnaseq_workflow.png)

### Seven stages to data science
1. Define the question of interest
2. Get the data
3. Clean the data
4. Explore the data
5. Fit statistical models
6. Communicate the results
7. Make your analysis reproducible

### What do we need to prepare ?
### Sample Information
a. What kind of material it is should be noted: Tissue, cell line, primary cell type, etc…
b. It’s ontology term (a DCC wrangler will work with you to obtain this)
c. If any treatments or genetic modifications (TALENs, CRISPR, etc…) were done to the sample
prior to RNA isolation.
d. If it’s a subcellular fraction or derived from another sample. If derived from another sample,
that relationship should be noted.
e. Some sense of sample abundance: RNA-Seq data from “bulk” vs. 10,000 cell equivalents can
give very different results, with lower input samples typically being less reproducible. Having
a sense of the amount of starting material here is useful.
f. If you received a batch of primary or immortalized cells, the lot #, cat # and supplier should be
noted.
g. If cells were cultured out, the protocol and methods used to propagate the cells should be
noted.
h. If any cell phenotyping or other characterizations were done to confirm it’s identify, purity,
etc.. those methods should be noted.


### RNA Information: 
RNAs come in all shapes and sizes. Some of the key properties to report are:
a. Total RNA, Poly-A(+) RNA, Poly-A(-) RNA
b. Size of the RNA fraction: we typically have a + 200 and – 200 cutoff, but there is a wide
range, i.e. microRNA-sized, etc…
c. If the RNA was treated with Ribosomal RNA depletion kits (RiboMinus, RiboZero): please
note the kit used.

### Protocols: 
There are several methods used to isolate RNAs with that work fine for the purposes of RNA-Seq. For all the ENCODE libraries that we make, we provide a document that lists in detail:
a. The RNA isolation methods,
b. Methods of size selections
c. Methods of rRNA removal
d. Methods of oligo-dT selections
e. Methods of DNAse I treatments

### Experimental Design
- Balanced design
- Technical replicates not necessary (Marioni et al., 2008)
- Biological replicates: 6 - 12 (Schurch et al., 2016)
- Power analysis

>## Reading materials
>[Paul L. Auer and R. W. Doerge "Statistical Design and Analysis of RNA Sequencing Data" Genetics June 1, 2010 vol. 185 no.2 405-416](https://www.genetics.org/content/185/2/405)  
>[Busby, Michele A., et al. "Scotty: a web tool for designing RNA-Seq experiments to measure differential gene expression." Bioinformatics 29.5 (2013): 656-657](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582267/)  
>[Marioni, John C., et al. "RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays." Genome research (2008)](https://genome.cshlp.org/content/18/9/1509.full.html)  
>[Schurch, Nicholas J., et al. "How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?." Rna (2016)](https://rnajournal.cshlp.org/content/22/6/839.long)  
>[Zhao, Shilin, et al. "RnaSeqSampleSize: real data based sample size estimation for RNA sequencing." BMC bioinformatics 19.1 (2018): 191](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5)  
{: .prereq}

### Replicate number
In all cases, experiments should be performed with two or more biological replicates, unless there is a
compelling reason why this is impractical or wasteful (e.g. overlapping time points with high temporal
resolution). A biological replicate is defined as an independent growth of cells/tissue and subsequent
analysis. Technical replicates from the same RNA library are not required, except to evaluate cases
where biological variability is abnormally high. In such instances, separating technical and biological
variation is critical. In general, detecting and quantifying low prevalence RNAs is inherently more variable
than high abundance RNAs. As part of the ENCODE pipeline, annotated transcript and genes are
quantified using RSEM and the values are made available for downstream correlation analysis. Replicate
concordance: the gene level quantification should have a Spearman correlation of >0.9 between
isogenic replicates and >0.8 between anisogenic replicates.

### RNA extraction
- Sample processing and storage
- Total RNA/mRNA/small RNA
- DNAse treatment
- Quantity & quality
- RIN values (Strong effect)
- Batch effect
- Extraction method bias (GC bias)

>## Reading materials
>Romero, Irene Gallego, et al. "RNA-seq: impact of RNA degradation on transcript quantification." BMC biology 12.1 (2014): 42  
>Kim, Young-Kook, et al. "Short structured RNAs with low GC content are selectively lost during extraction from a small number of cells." Molecular cell 46.6 (2012): 893-89500481-9).  
{: .prereq}

**RNA Quantification and Quality Control: When working with bulk samples, throughout the various steps we periodically assess the quality and quantity of the RNA. This is typically done on a BioAnalyzer. Points to check are:**
a. Total RNA
b. After oligo-dT size selections
c. After rRNA-depletions
d. After library construction

### Library prep
- PolyA selection
- rRNA depletion
- Size selection
- PCR amplification (See section PCR duplicates)
- Stranded (directional) libraries
   - Accurately identify sense/antisense transcript
   - Resolve overlapping genes
- Exome capture
- Library normalisation
- Batch effect

![RNA library]({{{site.baseurl}}/fig/library.png)

![RNA Sequencing tool]({{{site.baseurl}}/fig/frag.png)  


### Sequencing: 
There are several sequencing platforms and technologies out there being used. It is important to provide the following pieces of information:
a. Platform: Illumina, PacBio, Oxford Nanopore, etc…
b. Format: Single-end, Pair-end,
c. Read Length: 101 bases, 125 bases, etc…
d. Unusual barcode placement and sequence: Some protocols introduce barcodes in noncustomary places. If you are going to deliver a FASTQ file that will contain the barcode
sequences in it or other molecular markers – you will need to report both the position in the
read(s) where they are and their sequence(s).
e. Please provide the sequence of any custom primers that were used to sequence the library

![RNA library]({{{site.baseurl}}/fig/sequencing.png)



### Sequencing depth.
The amount of sequencing needed for a given sample is determined by the goals of the experiment and
the nature of the RNA sample. Experiments whose purpose is to evaluate the similarity between the
transcriptional profiles of two polyA+ samples may require only modest depths of sequencing.
Experiments whose purpose is discovery of novel transcribed elements and strong quantification of
known transcript isoforms requires more extensive sequencing.  
• Each Long RNA-Seq library must have a minimum of 30 million aligned reads/mate-pairs.  
• Each RAMPAGE library must have a minimum of 20 million aligned reads/mate-pairs.  
• Each small RNA-Seq library must have a minimum of 30 million aligned reads/mate-pairs.  


### Quantitative Standards (spike-ins).
It is highly desirable to include a ladder of RNA spike-ins to calibrate quantification, sensitivity, coverage
and linearity. Information about the spikes should include the stage of sample preparation that the spiked
controls were added, as the point of entry affects use of spike data in the output. In general, introducing
spike-ins as early in the process as possible is the goal, with more elaborate uses of different spikes at
different steps being optional (e.g. before poly A+ selection, at the time of cDNA synthesis, or just prior to
sequencing). Different spike-in controls are needed for each of the RNA types being analyzed (e.g. long
RNAs require different quantitative controls from short RNAs). Such standards are not yet available for all
RNA types. Information about quantified standards should also include:
a) A FASTA (or other standard format) file containing the sequences of each spike in.
b) Source of the spike-ins (home-made, Ambion, etc..)
c) The concentration of each of the spike-ins in the pool used. 

[Hong et al., 2016, Principles of metadata organization at the ENCODE data coordination center.](https://academic.oup.com/database/article-lookup/doi/10.1093/database/baw001)



![RNA Sequencing tool]({{{site.baseurl}}/fig/rnasoftware.png)


### QC FAIL?
https://sequencing.qcfail.com/


### Fastq format
FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores.


The format is similar to fasta though there are differences in syntax as well as integration of quality scores. Each sequence requires at least 4 lines:

1. The first line is the sequence header which starts with an ‘@’ (not a ‘>’!).
Everything from the leading ‘@’ to the first whitespace character is considered the sequence identifier.
Everything after the first space is considered the sequence description
2. The second line is the sequence.
3. The third line starts with ‘+’ and can have the same sequence identifier appended (but usually doesn’t anymore).
4. The fourth line are the quality scores

The FastQ sequence identifier generally adheres to a particular format, all of which is information related to the sequencer and its position on the flowcell. The sequence description also follows a particular format and holds information regarding sample information.

```
@A00261:180:HL7GCDSXX:2:1101:30572:1047/2
AAAATACATTGATGACCATCTAAAGTCTACGGCGTATGCGACTGATGAAGTATATTGCACCACCTGAGGGTGATGCTAATACTACTGTTGACGATAATGCTGATCTTCTTGCTAAGCTTAATATTGTTGGTGTTGAACCTAATGTTGGTG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF
@A00261:180:HL7GCDSXX:2:1101:21088:1094/2
ATCTCACATCGTTCCCTCAAGATTCTGAATTTTGGCAGCTCATTGCATTCTGTGCCGGCACTGGTGGTTCGATGCTTGTCATTGGTTCTGCTGCTGGTGTAGCCTTCATGGGGATGGAGAAAGTCGATTTCTTTTGGTATTTCCGAAAGG
+
FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00261:180:HL7GCDSXX:2:1101:21251:1125/2
CGGTGGAAAAGGAAACAGCTTTGGAAGGTTGATTCCATTACAGATTCGATTCGAAACTATGGTTCAGATTTCCGATCTTCCACGGGATTTGACAGAGGAGGTGCTCTCTAGGATTCCGGTGACATCTATGAGAGCAGTGAGATTTACTTG
```

- A00261  : Instrument name
- 180 : run ID
- HL7GCDSXX : Flowcell ID
- 2 : Flowcell lane
- 1101 : tile number within the flowcell lane
- 30572 : X-coordinate of the cluster within tile
- 1047 : Y-coordinate of the cluster within tile
- /2 : member of a pair 1 or 2 (Paired end reads only)


### Quality Scores
Quality scores are a way to assign confidence to a particular base within a read. Some sequencers have their own proprietary quality encoding but most have adopted Phred-33 encoding. Each quality score represents the probability of an incorrect basecall at that position.

### Phred Quality Score Encoding
Quality scores started as numbers (0-40) but have since changed to an ASCII encoding to reduce filesize and make working with this format a bit easier, however they still hold the same information. ASCII codes are assigned based on the formula found below. This table can serve as a lookup as you progress through your analysis.

### Quality Score Interpretation
Once you know what each quality score represents you can then use this chart to understand the confidence in a particular base.


![FASTQ quality]({{{site.baseurl}}/fig/quality.png)


## BW algorithm
![algorithm]({{{site.baseurl}}/fig/algorithm.png)
![banana]({{{site.baseurl}}/fig/banana.png)
![banana]({{{site.baseurl}}/fig/BackwardMatching.png)


### Reads QC
- Number of reads
- Per base sequence quality
- Per sequence quality score
- Per base sequence content
- Per sequence GC content
- Per base N content
- Sequence length distribution
- Sequence duplication levels
- Overrepresented sequences
- Adapter content
- Kmer content

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)


## FeatureCounts execute location

```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-5/students/${USER}/mouse/bam 

ls -1 *.bam

#### Copy templet
cp /data/gpfs/assoc/bch709-5/students/Course_materials/mouse/run.sh /data/gpfs/assoc/bch709-5/students/${USER}/mouse/bam/count.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Count/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-5/students/${USER}/mouse/bam/count.sh

```


### FeatureCounts read bam file
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/mouse/bam 
ls -1 *.sortedByCoord.out.bam
ls -1 *.sortedByCoord.out.bam| tr '\n' ' '
```
### Edit templet
```bash
nano count.sh
```

```bash
#paste this to count.sh
featureCounts -o /data/gpfs/assoc/bch709-5/students/${USER}//mouse/readcount/featucount -T 4 -Q 1 -p -M  -g gene_id -a /data/gpfs/assoc/bch709-5/students/${USER}/mouse/ref/refGene.gtf $(for i in `cat /data/gpfs/assoc/bch709-5/students/${USER}/mouse/filelist`; do echo ${i}.bamAligned.sortedByCoord.out.bam| tr '\n' ' ';done)
```


### Job submission dependency

```bash
squeue --noheader --format %i --user ${USER} 
```

### Submit
```bash
jobid=$(squeue --noheader --format %i --user ${USER} | tr '\n'  ':')1

sbatch --dependency=afterany:${jobid} count.sh 
```


# Human RNA-Seq
[***Transcriptome alterations in myotonic dystrophy frontal cortex***](https://www.sciencedirect.com/science/article/pii/S2211124720316235)

![](https://i.imgur.com/BrugOCz.png)

[](https://doi.org/10.1016/j.celrep.2020.108634)


## Environment activation
```bash
conda activate BCH709_RNASeq
```

## Working directory (Pronghorn)

```bash
echo $USER

cd /data/gpfs/assoc/bch709-5/students/${USER}

mkdir human
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/human/ref 
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/human/trim
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/human/bam
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/human/readcount
mkdir /data/gpfs/assoc/bch709-5/students/${USER}/human/DEG
```

## Reference Download

https://www.ncbi.nlm.nih.gov/genome/guide/human/


```bash
### change working directory
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/ref 


### download
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz


### decompress
gunzip GRCh38_latest_genomic.fna.gz
gunzip hg38.refGene.gtf.gz
```

## STAR reference build
### STAR aligner reference build on Pronghorn
```bash
### Copy templet
cp /data/gpfs/assoc/bch709-5/students/Course_materials/human/run.sh /data/gpfs/assoc/bch709-5/students/${USER}/human/ref/ref_build.sh

#open text editor
### PLEASE RENAME EMAIL AND JOB NAME
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/ref_build/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-5/students/${USER}/human/ref/ref_build.sh
nano ref_build.sh

# Add below command to ref_build.sh

STAR  --runThreadN 4 --runMode genomeGenerate --genomeDir . --genomeFastaFiles GRCh38_latest_genomic.fna --sjdbGTFfile GRCh38_latest_genomic.gtf  --sjdbOverhang 99   --genomeSAindexNbases 12
```


### Submit job to HPC
```bash
#submit job
sbach ref_build.sh

#check job
squeue -u ${USER}
```


### FASTQ file 
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/

### Link file (without copy)
ln -s /data/gpfs/assoc/bch709-5/students/Course_materials/human/fastq/* /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq

ls /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq
```

### Create file list
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq

ls -1 *.gz 

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g'

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u

ls -1 *.gz | sed 's/_R.\.fastq\.gz//g' | sort -u > /data/gpfs/assoc/bch709-5/students/${USER}/human/filelist

cat /data/gpfs/assoc/bch709-5/students/${USER}/human/filelist
```
### Regular expression

https://regex101.com/

## Trim reads

### Prepare templet
```bash
cp /data/gpfs/assoc/bch709-5/students/Course_materials/human/run.sh /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq/trim.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Trim/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq/trim.sh

```
### Edit templet
```bash
nano /data/gpfs/assoc/bch709-5/students/${USER}/human/fastq/trim.sh
```

### Batch submission
```bash
# Check file list
cat ../filelist
nano trim.sh

# Loop file list
### Add Forward read to variable
### Add reverse read from forward read name substitution
### add file name from variable to trim-galore
### merge trim-galore command and trim.sh
### add trim-galore command and trim.sh to new file
for i in `cat ../filelist`
    do
        read1=${i}_R1.fastq.gz
        read2=${read1//_R1.fastq.gz/_R2.fastq.gz}
        
        echo "trim_galore --paired  --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --fastqc --gzip -o /data/gpfs/assoc/bch709-5/students/${USER}/human/trim $read1 $read2" | cat trim.sh - > ${i}_trim.sh
        echo "$read1 $read2 trim file has been created."
done
```
### Batch submission
```bash
ls -1 *_trim.sh
### Loop *.sh submission
for i in `ls -1 *_trim.sh`
do
    sbatch $i
done
```

### Check submission
```bash
squeue -u ${USER}
```

## RNA-Seq Alignment
```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/trim

#### Copy templet
cp /data/gpfs/assoc/bch709-5/students/Course_materials/human/run.sh /data/gpfs/assoc/bch709-5/students/${USER}/human/trim/mapping.sh
sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Mapping/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-5/students/${USER}/human/trim/mapping.sh


#### Edit templet
nano mapping.sh
```


### Check output
```bash
ls -algh /data/gpfs/assoc/bch709-5/students/${USER}/human/trim
```
### Output example
```bash
[FILENAME]_R1_val_1.fq.gz [FILENAME]_R2_val_2.fq.gz
```


### STAR RNA-Seq alignment
```bash
STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --quantMode TranscriptomeSAM GeneCounts --genomeDir /data/gpfs/assoc/bch709-5/students/${USER}/human/ref  --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-5/students/${USER}/human/trim/[FILENAME]_R1_val_1.fq.gz  /data/gpfs/assoc/bch709-5/students/${USER}/human/trim/[FILENAME]_R2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/gpfs/assoc/bch709-5/students/${USER}/human/bam/[FILENAME].bam
```
### STAR RNA-Seq alignment batch file 

```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/trim
for i in `cat ../filelist`
    do
        read1=${i}_R1_val_1.fq.gz
        read2=${read1//_R1_val_1.fq.gz/_R2_val_2.fq.gz}
        echo $read1 $read2
        echo "STAR --runMode alignReads --runThreadN 4 --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 --genomeDir /data/gpfs/assoc/bch709-5/students/${USER}/human/ref --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --readFilesIn /data/gpfs/assoc/bch709-5/students/${USER}/human/trim/${read1} /data/gpfs/assoc/bch709-5/students/${USER}/human/trim/${read2} --outFileNamePrefix /data/gpfs/assoc/bch709-5/students/${USER}/human/bam/${i}.bam" | cat mapping.sh - > ${i}_mapping.sh
    done
```


### Job submission dependency

```bash
squeue --noheader --format %i --user ${USER}  
squeue --noheader --format %i --user ${USER} | tr '\n'  ':'
```


###  Job submission dependency on Trim
```bash

jobid=$(squeue --noheader --format %i --user ${USER} | tr '\n'  ':')1
for i in `ls -1 *_mapping.sh`
do
    sbatch  --dependency=afterany:${jobid} $i 
done

```

### FeatureCounts
[Bioinformatics, Volume 30, Issue 7, 1 April 2014, Pages 923–930](https://doi.org/10.1093/bioinformatics/btt656)
![]({{{site.baseurl}}/fig/featurecount.png)

```bash
featureCounts -o [output] -T [threads] -Q 1 -p -M  -g gene_id -a [GTF] [BAMs]
```
### FeatureCounts location

```bash
#### Move to trim folder
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/bam 

#### Copy templet
cp /data/gpfs/assoc/bch709-5/students/Course_materials/human/run.sh /data/gpfs/assoc/bch709-5/students/${USER}/human/bam/count.sh

sed -i "s/16g/64g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=4/g; s/\[NAME\]/Count/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g"  /data/gpfs/assoc/bch709-5/students/${USER}/human/bam/count.sh
```


### FeatureCounts command to count.sh
#### LOOP example
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/bam 

ls -1 *.bam 

for i in `cat /data/gpfs/assoc/bch709-5/students/${USER}/human/filelist`
do 
echo ${i}.bamAligned.sortedByCoord.out.bam | tr '\n' ' '
done
```



### FeatureCount 
```bash

echo "featureCounts -o /data/gpfs/assoc/bch709-5/students/${USER}//mouse/readcount/featucount -T 4 -Q 1 -p -M  -g gene_id -a /data/gpfs/assoc/bch709-5/students/${USER}/human/ref/GRCh38_latest_genomic.gtf $(for i in `cat /data/gpfs/assoc/bch709-5/students/${USER}/human/filelist`; do echo ${i}.bamAligned.sortedByCoord.out.bam| tr '\n' ' ';done)" >> count.sh
```


### Job submission dependency

```bash
squeue --noheader --format %i --user ${USER} 
squeue --noheader --format %i --user ${USER} | tr '\n'  ':'
```


###  Job submission dependency on Align
```bash
cd /data/gpfs/assoc/bch709-5/students/${USER}/human/bam 

jobid=$(squeue --noheader --format %i --user ${USER} | tr '\n'  ':')1
sbatch  --dependency=afterany:${jobid} count.sh 
```

####################################

### Slurm

![](https://i.imgur.com/LynACgh.png)


![](https://i.imgur.com/XEMRbJe.png)


[CheatSheet](https://slurm.schedmd.com/pdfs/summary.pdf)





