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
conda create -n RNASEQ_bch709 -c bioconda -c conda-forge  -c r  sra-tools minimap2 trinity star trim-galore gffread seqkit kraken2 samtools multiqc
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
#SBATCH -o index.out # STDOUT & STDERR
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
cd  ~/bch709_scratch/RNA-Seq_example/
mkdir Drosophila && cd Drosophila
mkdir raw_data trim bam reference
pwd
```

## Download reference
```
http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.42_FB2021_05/fasta/dmel-all-chromosome-r6.42.fasta.gz

http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.42_FB2021_05/gtf/dmel-all-r6.42.gtf.gz
```


