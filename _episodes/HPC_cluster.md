---
layout: page
title: 	HPC
published: true
---

## HPC clusters
>This exercise mainly deals with using HPC clusters for large scale data (Next Generation Sequencing analysis, Genome annotation, evolutionary studies etc.). These clusters have several processors with large amounts of RAM (compared to typical desktop/laptop), which makes it ideal for running programs that are computationally intensive. The operating system of these clusters are primarily UNIX and are mainly operated via command line. All the commands that you have learned in the previous exercises can be used on HPC.
>
>Pronghorn High Performance Computing offers shared cluster computing infrastructure for researchers and students at UNR. Brief descriptions for the available resources can be found here: https://www.unr.edu/research-computing/hpc. To begin with, you need to request permission for accessing these resources either through your department or through your advisor. All workshop attendees will have their account setup on HPC class education cluster and they can use their UNR NetID and the password for logging-in. You should have already received a confirmation email about your account creation with instructions on how to connect to the cluster. In this exercise we will specifically teach you how to connect to a remote server (HPC), transfer files in and out of the server, and running programs by requesting resources.
>
>You can log onto its front-end/job-submission system (pronghorn.rc.unr.edu) using your UNR NetID and password. Logging into HPC class requires an SSH client if you are >using Windows but Mac/Linux have these built into their OS. There are several available for download for the Windows platform.
>
>```bash
>ssh <YOURID>@pronghorn.rc.unr.edu
>```
{: .prereq}



**Pronghorn** is the University of Nevada, Reno's new High-Performance Computing (HPC) cluster. The GPU-accelerated system is designed, built and maintained by the Office of Information Technology's HPC Team. Pronghorn and the HPC Team supports general research across the Nevada System of Higher Education (NSHE).

Pronghorn is composed of CPU, GPU, and Storage subsystems interconnected by a 100Gb/s non-blocking Intel Omni-Path fabric. The CPU partition features 93 nodes, 2,976 CPU cores, and 21TiB of memory. The GPU partition features 44 NVIDIA Tesla P100 GPUs, 352 CPU cores, and 2.75TiB of memory. The storage system uses the IBM SpectrumScale file system to provide 1PB of high-performance storage. The computational and storage capabilities of Pronghorn will regularly expand to meet NSHE computing demands.

Pronghorn is collocated at the Switch Citadel Campus located 25 miles East of the University of Nevada, Reno. Switch is the definitive leader of sustainable data center design and operation. The Switch Citadel is rated Tier 5 Platinum, and will be the largest, most advanced data center campus on the planet.

![Pronghorn system map](../fig/pronghorn.png){: width="70%" height="70%"}



## File transfer

There are a number of ways to transfer data to and from HPC clusters. Which you should use depends on several factors, including the ease of use for you personally, connection speed and bandwidth, and the size and number of files which you intend to transfer. Most common options include scp, rsync (command line) and SCP and SFTP clients (GUI). scp (secure copy) is a simple way of transferring files between two machines that use the SSH (Secure SHell) protocol. You may use scp to connect to any system where you have SSH (login) access. scp is available as a protocol choice in some graphical file transfer programs and also as a command line program on most Linux, UNIX, and Mac OS X systems. scp can copy single files, but will also recursively copy directory contents if given a directory name. scp can be used as follows:

### to a remote system from local
```bash
scp sourcefile username@pronghorn.rc.unr.edu:somedirectory/
```
- example
```bash
echo "hello world" >> test_uploading_file.txt
scp test_uploading_file.txt username@pronghorn.rc.unr.edu:~/
```
  
### from a remote system to local
```bash
scp username@pronghorn.rc.unr.edu:somedirectory/sourcefile destinationfile
```
- example
```bash
scp username@pronghorn.rc.unr.edu:~/test_uploading_file.txt ~/
```

### recursive directory copy to a remote system from local
```bash
scp -r SourceDirectory/ username@pronghorn.rc.unr.edu:somedirectory/
```

## Open location
- Windows (WSL)
```bash
explorer.exe .
```
- Mac
```bash
open .
```
  
### Rsync
  
***rsync*** is a fast and extraordinarily versatile file copying tool. It can synchronize file trees across local disks, directories or across a network

- Synchronize a local directory with the remote server directory
```bash
rsync -avhP path/to/SourceDirectory username@pronghorn.rc.unr.edu:somedirectory/
```

- Synchronize a remote directory with the local directory
```bash
rsync -avhP username@hpronghorn.rc.unr.edu:SourceDirectory/ path/to/Destination/
```

> ## Pronghorn conda installation.
> 
>To install miniconda3, type:
>
>~~~
>curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
>bash Miniconda3-latest-Linux-x86_64.sh
>~~~
>{: .bash}
> Then, follow the instructions that you are prompted with on the screen to install Miniconda3. 
{: .solution}

![conda1]({{site.baseurl}}/fig/conda_excute.png)
![conda2]({{site.baseurl}}/fig/conda_excute2.png)

  
## Prompt Customization for Pronghorn

```bash
echo '###BCH709 ' >> ~/.bashrc

echo 'tty -s && export PS1="\[\033[38;5;164m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;231m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;172m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bashrc
echo "alias ls='ls --color=auto'" >> ~/.bashrc

```

```bash
source ~/.bashrc
```

## Conda cheatsheet
[cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)
  
  
>## Assignment due next Monday
>### Create your environment name `rnaseq` on Pronghorn, install following tools through Conda
> 1. Login to Pronghorn
> 2. Create environment `rnaseq`
> 3. Activate your environment
> 4. Install below software
>- star  
>- fastqc  
>- rsem  
>- subread  
>- hisat2  
>- samtools  
>- bowtie2  
>- trim-galore   
>- multiqc  
###  export your environment to rnaseq.yaml
>```bash
> conda env export  > rnaseq.yaml
>```
>### copy contents of rnaseq.yaml and paste to Webcanvas 
>```bash
> cat rnaseq.yaml
>```
{: .solution}



### References:
- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

