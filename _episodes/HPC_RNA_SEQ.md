---
layout: page
title: 	HPC_RNA_Seq
published: true
---


## Using Pronghorn (High-Performance Computing)

**Pronghorn** is the University of Nevada, Reno's new High-Performance Computing (HPC) cluster. The GPU-accelerated system is designed, built and maintained by the Office of Information Technology's HPC Team. Pronghorn and the HPC Team supports general research across the Nevada System of Higher Education (NSHE).

Pronghorn is composed of CPU, GPU, and Storage subsystems interconnected by a 100Gb/s non-blocking Intel Omni-Path fabric. The CPU partition features 93 nodes, 2,976 CPU cores, and 21TiB of memory. The GPU partition features 44 NVIDIA Tesla P100 GPUs, 352 CPU cores, and 2.75TiB of memory. The storage system uses the IBM SpectrumScale file system to provide 1PB of high-performance storage. The computational and storage capabilities of Pronghorn will regularly expand to meet NSHE computing demands.

Pronghorn is collocated at the Switch Citadel Campus located 25 miles East of the University of Nevada, Reno. Switch is the definitive leader of sustainable data center design and operation. The Switch Citadel is rated Tier 5 Platinum, and will be the largest, most advanced data center campus on the planet.

![Pronghorn system map](../fig/pronghorn.png){: width="70%" height="70%"}


## Slurm Start Tutorial
Resource sharing on a supercomputer dedicated to technical and/or scientific computing is often organized by a piece of software called a resource manager or job scheduler. Users submit jobs, which are scheduled and allocated resources (CPU time, memory, etc.) by the resource manager.

Slurm is a resource manager and job scheduler designed to do just that, and much more. It was originally created by people at the Livermore Computing Center, and has grown into a full-fledge open-source software backed up by a large community, commercially supported by the original developers, and installed in many of the Top500 supercomputers.

Gathering information
Slurm offers many commands you can use to interact with the system. For instance, the sinfo command gives an overview of the resources offered by the cluster, while the squeue command shows to which jobs those resources are currently allocated.

By default, sinfo lists the partitions that are available. A partition is a set of compute nodes (computers dedicated to... computing) grouped logically. Typical examples include partitions dedicated to batch processing, debugging, post processing, or visualization.

### sinfo
```bash
sinfo
```

```
PARTITION      AVAIL  TIMELIMIT  NODES  STATE NODELIST
cpu-s2-core-0     up 14-00:00:0      2    mix cpu-[8-9]
cpu-s2-core-0     up 14-00:00:0      7  alloc cpu-[1-2,4-6,78-79]
cpu-s2-core-0     up 14-00:00:0     44   idle cpu-[0,3,7,10-47,64,76-77]
cpu-s3-core-0*    up    2:00:00      2    mix cpu-[8-9]
cpu-s3-core-0*    up    2:00:00      7  alloc cpu-[1-2,4-6,78-79]
cpu-s3-core-0*    up    2:00:00     44   idle cpu-[0,3,7,10-47,64,76-77]
gpu-s2-core-0     up 14-00:00:0     11   idle gpu-[0-10]
cpu-s6-core-0     up      15:00      2   idle cpu-[65-66]
cpu-s1-pgl-0      up 14-00:00:0      1    mix cpu-49
cpu-s1-pgl-0      up 14-00:00:0      1  alloc cpu-48
cpu-s1-pgl-0      up 14-00:00:0      2   idle cpu-[50-51]

```
In the above example, we see two partitions, named batch and debug. The latter is the default partition as it is marked with an asterisk. All nodes of the debug partition are idle, while two of the batch partition are being used.

The sinfo command also lists the time limit (column TIMELIMIT) to which jobs are subject. On every cluster, jobs are limited to a maximum run time, to allow job rotation and let every user a chance to see their job being started. Generally, the larger the cluster, the smaller the maximum allowed time. You can find the details on the cluster page.

You can actually specify precisely what information you would like sinfo to output by using its --format argument. For more details, have a look at the command manpage with man sinfo.

## squeue
The squeue command shows the list of jobs which are currently running (they are in the RUNNING state, noted as ‘R’) or waiting for resources (noted as ‘PD’, short for PENDING).
```bash
squeue
```
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            983204 cpu-s2-co    neb_K jzhang23  R 6-09:05:47      1 cpu-6
            983660 cpu-s2-co   RT3.sl yinghanc  R   12:56:17      1 cpu-9
            983659 cpu-s2-co   RT4.sl yinghanc  R   12:56:21      1 cpu-8
            983068 cpu-s2-co Gd-bound   dcantu  R 7-06:16:01      2 cpu-[78-79]
            983067 cpu-s2-co Gd-unbou   dcantu  R 1-17:41:56      2 cpu-[1-2]
            983472 cpu-s2-co   ub-all   dcantu  R 3-10:05:01      2 cpu-[4-5]
            982604 cpu-s1-pg     wrap     wyim  R 12-14:35:23      1 cpu-49
            983585 cpu-s1-pg     wrap     wyim  R 1-06:28:29      1 cpu-48
            983628 cpu-s1-pg     wrap     wyim  R   13:44:46      1 cpu-49
```


## Text editor



## SBATCH
Now the question is: How do you create a job?

A job consists in two parts: resource requests and job steps. Resource requests consist in a number of CPUs, computing expected duration, amounts of RAM or disk space, etc. Job steps describe tasks that must be done, software which must be run.

The typical way of creating a job is to write a submission script. A submission script is a shell script, e.g. a Bash script, whose comments, if they are prefixed with SBATCH, are understood by Slurm as parameters describing resource requests and other submissions options. You can get the complete list of parameters from the sbatch manpage man sbatch.

>## Important
>
>The SBATCH directives must appear at the top of the submission file, before any other line except for the very first line which should be the shebang (e.g. #!/bin/bash).
>The script itself is a job step. Other job steps are created with the srun command.
>For instance, the following script, hypothetically named submit.sh,
{: checklist}

```bash
nano submit.sh
```

```
#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mail-type=all
#SBATCH --mail-user=wyim@unr.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1g
#SBATCH --time=8:10:00
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

for i in {1..1000}; 
do 
    echo $i;
    sleep 1; 
done

```
would request one CPU for 10 minutes, along with 1g of RAM, in the default queue. When started, the job would run a first job step srun hostname, which will launch the UNIX command hostname on the node on which the requested CPU was allocated. Then, a second job step will start the sleep command. Note that the --job-name parameter allows giving a meaningful name to the job and the --output parameter defines the file to which the output of the job must be sent.

Once the submission script is written properly, you need to submit it to slurm through the sbatch command, which, upon success, responds with the jobid attributed to the job. (The dollar sign below is the shell prompt)
```bash
chmod 775 submit.sh
sbatch submit.sh
sbatch: Submitted batch job 99999999
```



## How to cancel the job?

```bash
scancel <JOB ID>
```


## Create scratch disk space
![Pronghorn system map](../fig/pronghorn.png){: width="70%" height="70%"}

```bash
cd /data/gpfs/assoc/bch709-2/
mkdir $(whoami) 
cd ~/
ln -s /data/gpfs/assoc/bch709-2/$(whoami) bch709_scratch
cd bch709_scratch
```

## Importing Data from the NCBI Sequence Read Archive (SRA) using the DE

### WORKTING PATH
```bash
cd bch709_scratch
mkdir  RNA-Seq_example/
cd ~/bch709_scratch/RNA-Seq_example/
pwd
```



### Conda environment
```bash
conda create -n RNASEQ_bch709 -c bioconda -c conda-forge  -c r  sra-tools minimap2 trinity star trim-galore gffread seqkit kraken2 samtools multiqc subread
conda activate RNASEQ_bch709
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
mkdir ~/bch709_scratch/RNA-Seq_example/
cd ~/bch709_scratch/RNA-Seq_example/
mkdir ATH && cd ATH
mkdir raw_data
mkdir trim
pwd
```

### fastq-dump submission
```bash
cd ~/bch709_scratch/RNA-Seq_example/ATH
nano fastq-dump.sh
```
```bash
#!/bin/bash
#SBATCH --job-name=fastqdump_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o fastq-dump.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

fastq-dump SRR1761506 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761507 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761508 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761509 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761510 --split-3 --outdir ./raw_data  --gzip
fastq-dump SRR1761511 --split-3 --outdir ./raw_data  --gzip
```


## Trim-galore
```bash
cd  ~/bch709_scratch/RNA-Seq_example/ATH
nano trim.sh
```
```bash
#!/bin/bash
#SBATCH --job-name=trim_ATH
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761506 raw_data/SRR1761506_1.fastq.gz raw_data/SRR1761506_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761507  raw_data/SRR1761507_1.fastq.gz raw_data/SRR1761507_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761508 raw_data/SRR1761508_1.fastq.gz raw_data/SRR1761508_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761509 raw_data/SRR1761509_1.fastq.gz raw_data/SRR1761509_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761510 raw_data/SRR1761510_1.fastq.gz raw_data/SRR1761510_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR1761511 raw_data/SRR1761511_1.fastq.gz raw_data/SRR1761511_2.fastq.gz --fastqc
```

## Reference downloads

### Please create account in JGI

https://contacts.jgi.doe.gov/registration/new


```bash
cd ~/bch709_scratch/RNA-Seq_example/ATH
mkdir bam
mkdir reference && cd reference
pwd
```

## Download Arabidopsis thaliana TAIR10
https://phytozome-next.jgi.doe.gov/info/Athaliana_TAIR10

```
Athaliana_167_gene.gff3.gz
Athaliana_167.fa.gz 
```


## Unzip file
```bash
cd ~/bch709_scratch/RNA-Seq_example/ATH/reference
unzip download.#######.zip
zcat phytozome/phyto_mirror/Athaliana_167_10/assembly/Athaliana_167.fa.gz | head
zcat phytozome/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.gene.gff3.gz | head
gunzip phytozome/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.gene.gff3.gz 
gunzip phytozome/phyto_mirror/Athaliana_167_10/assembly/Athaliana_167.fa.gz
```
**Your location might be different**

## Convert GFF to GTF
```bash

gffread phytozome/Athaliana/TAIR10/annotation/Athaliana_167_TAIR10.gene.gff3 -T -F --keep-exon-attrs -o TAIR10_GFF3_genes.gtf
```
## Create reference index
```bash
cd  ~/bch709_scratch/RNA-Seq_example/ATH/reference
ls -algh
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
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles   phytozome/phyto_mirror/Athaliana_167_10/assembly/Athaliana_167.fa  --sjdbGTFfile TAIR10_GFF3_genes.gtf --sjdbOverhang 99   --genomeSAindexNbases 12
```

## Mapping the reads to genome index
```bash
cd  ~/bch709_scratch/RNA-Seq_example/ATH/
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
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0
#SBATCH --dependency=afterok:<PREVIOUS_JOBID(trim_ATH)>

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761506_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761506.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761507_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761507_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761507.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761508_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761508_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761508.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761509_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761509_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761509.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761510_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761510_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761510.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 --genomeDir ~/bch709_scratch/RNA-Seq_example/ATH/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761511_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/ATH/trim/SRR1761511_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/ATH/bam/SRR1761511.bam
```
```bash
conda install -c conda-forge tree
```
# Drosophila
> ## Publication (Drosophila)
> 
> Not published.
> 
{: .callout}


### SRA Bioproject site 

```bash
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770108
```
**Gene expression profiling of Drosophila melanogaster larval brains after chronich alcohol exposure (fruit fly)**

We sequenced mRNA extracted from brains of (1) D. melanogaster larvae exposed to food containing 5% ethanol (v/v) for 6 conscutive days, and (2) an age-matched untreated control larvae, that grew in regular food. Differential gene expression between the two groups was calculated and reported. Each group consisted of 3 biological replicates of 30 brains each. Overall design: Examination of mRNA levels in brains of D. melanogaster larvae after chronich ethanol exposure was performed using next generation sequencing (NGS) technology (RNA-seq)


## Subset of data


| Sample information  | Run         |
|---------------------|-------------|
| Control         | SRR16287545 |
| Control          | SRR16287546 |
| Control          | SRR16287547 |
| Ethanol treatment         | SRR16287549 |
| Ethanol treatment          | SRR16287548 |
| Ethanol treatment          | SRR16287550 |


```bash
cd  ~/bch709_scratch/RNA-Seq_example/
mkdir Drosophila && cd Drosophila
mkdir raw_data trim bam reference
pwd
```



## fastq donwload

```bash
cd ~/bch709_scratch/RNA-Seq_example/Drosophila

nano fastq-dump.sh
```
```bash
#!/bin/bash
#SBATCH --job-name=fastqdump_Drosophila
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<youremail>
#SBATCH -o fastq-dump.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

fastq-dump SRR16287545 --split-3 --outdir ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data --gzip
fastq-dump SRR16287546 --split-3 --outdir ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data --gzip
fastq-dump SRR16287547 --split-3 --outdir ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data --gzip
fastq-dump SRR16287549 --split-3 --outdir ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data --gzip
fastq-dump SRR16287548 --split-3 --outdir ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data --gzip
fastq-dump SRR16287550 --split-3 --outdir ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data --gzip
```


## fastq trim
```bash
cd ~/bch709_scratch/RNA-Seq_example/Drosophila
mkdir trim
nano trim.sh

```

```bash
#!/bin/bash
#SBATCH --job-name=trim_Drosophila
#SBATCH --cpus-per-task=2
#SBATCH --time=2-15:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o trim.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR16287545 ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287545_1.fastq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287545_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR16287546  ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287546_1.fastq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287546_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR16287547 ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287547_1.fastq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287547_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR16287549 ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287549_1.fastq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287549_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR16287548 ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287548_1.fastq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287548_2.fastq.gz --fastqc
trim_galore --paired   --three_prime_clip_R1 5 --three_prime_clip_R2 5 --cores 2  --max_n 40  --gzip -o trim --basename SRR16287550 ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287550_1.fastq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/raw_data/SRR16287550_2.fastq.gz --fastqc
```
## Reference donwload

```bash
cd  ~/bch709_scratch/RNA-Seq_example/Drosophila/reference
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.42_FB2021_05/fasta/dmel-all-chromosome-r6.42.fasta.gz 
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.42_FB2021_05/gtf/dmel-all-r6.42.gtf.gz
gunzip dmel-all-chromosome-r6.42.fasta.gz
gunzip dmel-all-r6.42.gtf.gz
ls -algh
```

## Reference index

```
nano index.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=index_Drosophila
#SBATCH --cpus-per-task=12
#SBATCH --time=2-15:00:00
#SBATCH --mem=48g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o index.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0

STAR  --runThreadN 48g --runMode genomeGenerate --genomeDir . --genomeFastaFiles  dmel-all-chromosome-r6.42.fasta --sjdbGTFfile dmel-all-r6.42.gtf --sjdbOverhang 99   --genomeSAindexNbases 12
```


## Mapping
```
nano mapping.sh
```
```bash
#!/bin/bash
#SBATCH --job-name=align_Drosophila
#SBATCH --cpus-per-task=8
#SBATCH --time=2-15:00:00
#SBATCH --mem=32g
#SBATCH --mail-type=all
#SBATCH --mail-user=<PLEASE CHANGE THIS TO YOUR EMAIL>
#SBATCH -o align.out # STDOUT & STDERR
#SBATCH --account=cpu-s5-bch709-2
#SBATCH --partition=cpu-core-0
#SBATCH --dependency=afterok:<PREVIOUS_JOBID(trim_Drosophila)>

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir ~/bch709_scratch/RNA-Seq_example/Drosophila/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287547_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287547_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/Drosophila/bam/SRR16287547.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir ~/bch709_scratch/RNA-Seq_example/Drosophila/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287548_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287548_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/Drosophila/bam/SRR16287548.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir ~/bch709_scratch/RNA-Seq_example/Drosophila/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287549_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287549_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/Drosophila/bam/SRR16287549.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir ~/bch709_scratch/RNA-Seq_example/Drosophila/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287550_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287550_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/Drosophila/bam/SRR16287550.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir ~/bch709_scratch/RNA-Seq_example/Drosophila/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287545_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287545_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/Drosophila/bam/SRR16287545.bam

STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 100000 --genomeDir ~/bch709_scratch/RNA-Seq_example/Drosophila/reference/ --readFilesIn ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287546_val_1.fq.gz ~/bch709_scratch/RNA-Seq_example/Drosophila/trim/SRR16287546_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/bch709_scratch/RNA-Seq_example/Drosophila/bam/SRR16287546.bam
```



# Mus Musculus
## Data Download
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA773499 

**CCR2-dependent monocyte-derived cells restrict SARS-CoV-2 infection (house mouse)**  

SARS-CoV-2 has caused a historic pandemic of respiratory disease (COVID-19) and current evidence suggests severe disease is associated with dysregulated immunity within the respiratory tract1,2. However, the innate immune mechanisms that mediate protection during COVID-19 are not well defined. Here we characterize a mouse model of SARS-CoV-2 infection and find that early CCR2-dependent infiltration of monocytes restricts viral burden in the lung. We find that a recently developed mouse-adapted MA-SARS-CoV-2 strain, as well as the emerging B.1.351 variant, trigger an inflammatory response in the lung characterized by expression of pro-inflammatory cytokines and interferon-stimulated genes. Using intravital antibody labeling, we demonstrate that MA-SARS-CoV-2 infection leads to increases in circulating monocytes and an influx of CD45+ cells into the lung parenchyma that is dominated by monocyte-derived cells. scRNA-seq analysis of lung homogenates identified a hyper-inflammatory monocyte profile. We utilize this model to demonstrate that mechanistically, CCR2 signaling promotes infiltration of classical monocytes into the lung and expansion of monocyte-derived cells. Parenchymal monocyte-derived cells appear to play a protective role against MA-SARS-CoV-2, as mice lacking CCR2 showed higher viral loads in the lungs, increased lung viral dissemination, and elevated inflammatory cytokine responses. These studies have identified a CCR2-monocyte axis that is critical for promoting viral control and restricting inflammation within the respiratory tract during SARS-CoV-2 infection. Overall design: 8 samples in total corresponding to different mice. 4 samples are from mock, control mice. 4 samples are from SARS-CoV-2 infected mice.

```bash
cd  ~/bch709_scratch/RNA-Seq_example/
mkdir Mmusculus && cd Mmusculus
mkdir raw_data trim bam reference
pwd
```

| Run ID      | LibraryName                   |
|-------------|-------------------------------|
| SRR16526489 | Mock 1; Mus musculus; RNA-Seq |
| SRR16526488 | Mock 2; Mus musculus; RNA-Seq |
| SRR16526486 | Mock 3; Mus musculus; RNA-Seq |
| SRR16526483 | Mock 4; Mus musculus; RNA-Seq |
| SRR16526477 | CoV2 3; Mus musculus; RNA-Seq |
| SRR16526479 | CoV2 2; Mus musculus; RNA-Seq |
| SRR16526481 | CoV2 1; Mus musculus; RNA-Seq |
| SRR16526475 | CoV2 4; Mus musculus; RNA-Seq |



## Reference download
https://www.ncbi.nlm.nih.gov/genome/?term=Mus+musculus

### Download files
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz



# Solanum lycopersicum

## Project site
Whole genome sequencing and transcriptome sequencing of Solanum lycopersicum, M82
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA753098

```bash
cd  ~/bch709_scratch/RNA-Seq_example/
mkdir Slycopersium && cd Slycopersium
mkdir raw_data trim bam reference
pwd
```


| Run ID      | LibraryName              |
|-------------|--------------------------|
| SRR15607542 | Root control Rep1        |
| SRR15607543 | Root control Rep1        |
| SRR15607544 | Root control Rep1        |
| SRR15607552 | Root Salt treatment Rep1 |
| SRR15607553 | Root Salt treatment Rep2 |
| SRR15607554 | Root Salt treatment Rep3 |



## Reference Download
https://phytozome-next.jgi.doe.gov/info/Slycopersicum_ITAG4_0


### Download files
Slycopersicum_691_ITAG4.0.gene.gff3.gz  
Slycopersicum_691_SL4.0.fa.gz   





# Expression values and Normalization

CPM, RPKM, FPKM, TPM, RLE, MRN, Q, UQ, TMM, VST, RLOG, VOOM ... Too many...  

CPM: Controls for sequencing depth when dividing by total count. Not for within-sample comparison or DE.  

Counts per million (CPM) mapped reads are counts scaled by the number of fragments you sequenced (N) times one million. This unit is related to the FPKM without length normalization and a factor of 10^6:  
![CPM]({{site.baseurl}}/fig/CPM.png)

RPKM/FPKM: Controls for sequencing depth and gene length. Good for technical replicates, not good for sample-sample due to compositional bias. Assumes total RNA output is same in all samples. Not for DE.  

TPM: Similar to RPKM/FPKM. Corrects for sequencing depth and gene length. Also comparable between samples but no correction for compositional bias.  

TMM/RLE/MRN: Improved assumption: The output between samples for a core set only of genes is similar. Corrects for compositional bias. Used for DE. RLE and MRN are very similar and correlates well with sequencing depth. edgeR::calcNormFactors() implements TMM, TMMwzp, RLE & UQ.   DESeq2::estimateSizeFactors implements median ratio method (RLE). Does not correct for gene length.  

VST/RLOG/VOOM: Variance is stabilised across the range of mean values. For use in exploratory analyses. Not for DE. vst() and rlog() functions from DESeq2. voom() function from Limma converts data to normal distribution.  

geTMM: Gene length corrected TMM.  

**For DEG using DEG R packages (DESeq2, edgeR, Limma etc), use raw counts  
For visualisation (PCA, clustering, heatmaps etc), use TPM or TMM  
For own analysis with gene length correction, use TPM (maybe geTMM?)  
Other solutions: spike-ins/house-keeping genes**


## FPKM
![FPKM]({{site.baseurl}}/fig/FPKM.png)

X = mapped reads count
N = number of reads
L = Length of transcripts


'length' is this transcript's sequence length (poly(A) tail is not counted). 'effective_length' counts only the positions that can generate a valid fragment.



### FPKM 
Fragments per Kilobase of transcript per million mapped reads


```python
X = 3752
Number_Reads_mapped = 559192
Length = 651.04
fpkm= X*(1000/Length)*(1000000/Number_Reads_mapped)
fpkm
```

#### ten to the ninth power = 10\*\*9


### TPM
 Transcripts Per Million

![TPM]({{site.baseurl}}/fig/TPM.png)

![TPM2]({{site.baseurl}}/fig/TPM2.png)



### Paper read
[Li et al., 2010, RSEM](http://bioinformatics.oxfordjournals.org/content/26/4/493.long)  

[Dillies et al., 2013](http://bib.oxfordjournals.org/content/14/6/671.full)



# Featurecount
```
featureCounts -p  -a <GENOME>.gtf <SAMPLE1>.bam <SAMPLE2>.bam <SAMPLE3>.bam  ...... -o counts.txt
```
```bash
conda activate RNASEQ_bch709
cd ~/bch709_scratch/RNA-Seq_example/ATH/bam
featureCounts -o ATH.featureCount.cnt -p  -a ~/bch709/wyim/RNA-Seq_example/ATH/reference/TAIR10_GFF3_genes.gtf SRR1761506.bamAligned.sortedByCoord.out.bam  SRR1761509.bamAligned.sortedByCoord.out.bam SRR1761507.bamAligned.sortedByCoord.out.bam  SRR1761510.bamAligned.sortedByCoord.out.bam SRR1761508.bamAligned.sortedByCoord.out.bam  SRR1761511.bamAligned.sortedByCoord.out.bam
```




### FPKM 
Fragments per Kilobase of transcript per million mapped reads

```
awk 'FNR > 2 { sum+=$7 } END {print sum}' ATH.featureCount.cnt
```
### Example AT1G01060.TAIR10 
```python
Length = 976 
X = 500
Number_Reads_mapped = 5949384
fpkm= X*(1000/Length)*(1000000/Number_Reads_mapped)
fpkm
```

#### ten to the ninth power = 10\*\*9

```python
fpkm=X/(Number_Reads_mapped*Length)*10**9
fpkm
```
## TPM
### sum_count_per_length
```
awk 'FNR > 2 { sum+=$7/$6 } END {print sum}' ATH.featureCount.cnt
```
### TPM calculation from reads count
```python

sum_count_per_length =  4747.27
X = 500
Length = 976
TPM = (X/Length)*(1/sum_count_per_length )*10**6
TPM
```


### TPM and FPKM calculation

```bash
cut -f1,6-  ATH.featureCount.cnt |  egrep -v "#" | sed 's/\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g' > ATH.featureCount_count_length.cnt

python /data/gpfs/assoc/bch709-1/Course_material/script/tpm_raw_exp_calculator.py -count ATH.featureCount_count_length.cnt

```


### TPM calculation from FPKM

```python
FPKM = 86.10892858272605
SUM_FPKM = 797942
TPM=(FPKM/SUM_FPKM)*10**6
TPM
```



# Mosquito (Anopheles stephensi)
RNAseq from adult male and female Anopheles stephensi
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA277477


## Folder preparation
```bash
cd  ~/bch709_scratch/RNA-Seq_example/  
mkdir Astephensi && cd Astephensi  
mkdir raw_data trim bam reference  
pwd 
```

## SRA read download

| Run ID     | LibraryName                                   |
|------------|-----------------------------------------------|
| SRR1851022 | Anopheles stephensi male RNAseq replicate 1   |
| SRR1851024 | Anopheles stephensi male RNAseq replicate 2   |
| SRR1851026 | Anopheles stephensi male RNAseq replicate 3   |
| SRR1851027 | Anopheles stephensi female RNAseq replicate 1 |
| SRR1851028 | Anopheles stephensi female RNAseq replicate 2 |
| SRR1851030 | Anopheles stephensi female RNAseq replicate 3 |

## Reference genome (VectorBase)
https://vectorbase.org/vectorbase/app/record/dataset/TMPTX_asteIndian

### Reference Download Link
https://vectorbase.org/common/downloads/Current_Release/AstephensiSDA-500/fasta/data/VectorBase-54_AstephensiSDA-500_Genome.fasta
https://vectorbase.org/common/downloads/Current_Release/AstephensiSDA-500/gff/data/VectorBase-54_AstephensiSDA-500.gff



