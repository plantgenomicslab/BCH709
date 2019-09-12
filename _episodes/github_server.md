---
layout: page
title: Compile and Software installation
published: true
---

{% include gh_variables.html %}


## Last class question
1. I can't paste (Ubuntu)
![paste]({{{site.baseurl}}/fig/paste.png)
![paste2]({{{site.baseurl}}/fig/paste2.png)

2. ChIP-Seq
Maybe we need a special session.


## Check your CPUs and Memory

```
$ lscpu
$ free
$ htop
```

If your have error ***command not found***
Please check below site.
https://command-not-found.com/

```Mac
brew install htop util-linux
```

```Ubuntu
sudo apt install htop util-linux
```

### Prompt Customization

```
echo '###BCH709 ' >> ~/.bashrc

echo 'export PS1="\[\033[38;5;2m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;33m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;166m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;4m\]\\$\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bashrc

echo "alias ls='ls --color=auto'" >> ~/.bashrc

source ~/.bashrc
```
More information is [here](https://plantgenomicslab.github.io/BCH709/bash/index.html)

```
echo '###BCH709 ' >> ~/.bas_profile

echo 'export PS1="\[\033[38;5;2m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;33m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;166m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;4m\]\\$\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bash_profile

echo "alias ls='ls --color=auto'" >> ~/.bash_profile

source ~/.bash_profile
```


## Software

| Software | Version | Manual | Available for | Description |
| -------- | ------------ | ------ | ------------- | ----------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.7 | [Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)| Linux, MacOS, Windows | Quality control tool for high throughput sequence data. |
| [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) | 2.1.0| [Link](https://ccb.jhu.edu/software/hisat2/index.shtml) | Linux, MacOS, Windows | Mapping RNA sequences against genome |
| [BWA](http://bio-bwa.sourceforge.net/) | 0.7.17 | [Link](http://bio-bwa.sourceforge.net/bwa.shtml) | Linux, MacOS | Mapping DNA sequences against reference genome. |

## QuickStart Software Installation Instructions
​
First, I would recommend to put all the software under `~/bch709/bin` folder.
Please use `mkdir -p ~/bch709/bin` 


### Download example
Second please make example folder below `~/bch709/`
Such as `~/bch709/example`  and download below files.

#### Download reads example 
```
https://www.dropbox.com/s/p11a5cw5lb5y0un/test.fq
```

#### Download genome example
```
https://www.dropbox.com/s/q2srdymrp76jevm/test_genome.fasta
```

## Current location
Please check your current location and move to   `~/bch709/example`


## FastQC Source Code Installation
Please download and unzip
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
> ## MacOS
> 
>To install unzip: 
>
>~~~
>$ brew install unzip
>~~~
>{: .bash}
> . 
{: .solution}

> ## On Ubuntu systems
> 
>To install unzip: 
>
>~~~
>$ sudo apt install unzip
>~~~
>{: .bash}
> . 
{: .solution}

### Please check  INSTALL.txt & README.md
```
less <FILENAME>
```

### Run FastQC

**Test your installation by running:**

```
./fastqc 
```

### Which error did you get?



### JAVA Installation
```bash
sudo apt install default-jre
```

```bash
brew cask install java
```


### Run FastQC
```
./fastqc <your input>
```

### Result
- HTML file

### Link Your Directory (Short Cut in WINDOWS)
```bash
$ mkdir /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop 
$ ln -s /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop ~/bch709/results
```

```bash
$ mkdir ~/Desktop/BCH709_Desktop
$ ln -s ~/Desktop/BCH709_Desktop ~/bch709/results
```

### Check your results
unzip your results




## HISAT2
https://ccb.jhu.edu/software/hisat2/index.shtml

![hisat2]({{{site.baseurl}}/fig/hisat2.png)


### Check HISAT2 manual
```
less MANUAL
```

### compile your HISAT2
use make

### How to run HISAT2?
- need to build index first `hisat2-build`
```bash
./hisat2-build <YOUR_GENOME_SEQUENCE> <YOUR_GENOME_INDEX>
```
- map your reads by `hisat2`
```bash
./hisat2 -x <YOUR_GENOME_INDEX> -U <SEQUENCING_READS> -S output.sam
```

## Different way to install software?
## Github
### What is Git and why use it?
Git is an open source distributed version control system that facilitates GitHub activities on your laptop or desktop. Version control using Git is the most reasonable way to keep track of changes in code, manuscripts, presentations, and data analysis projects.  

### Why Version Control?
Version control is essential in creating any project that takes longer than 5 minutes to complete. Even if your memory is longer than 5 minutes, next month you are not likely to be able to retrace your steps.  
![github-workflow]({{{site.baseurl}}/fig/git_overview.png)

![github]({{{site.baseurl}}/fig/github_dri.png)

## HISAT2 installation from Github
```
https://github.com/DaehwanKimLab/hisat2

```
![github-workflow]({{{site.baseurl}}/fig/hisat2_github.png)
![github-workflow]({{{site.baseurl}}/fig/hisat2_github2.png)

### Download
```bash

$ cd ~/bch709/bin

$ git clone https://github.com/DaehwanKimLab/hisat2.git

$ cd hisat2

$ ls -algh

$ less MANUAL
```
![github-workflow]({{{site.baseurl}}/fig/hisat2_github3.png)

### Compile
```bash
$ make -j <YOUR CPU>
```

## HISAT2 installation from binary
![github-workflow]({{{site.baseurl}}/fig/hisat2_binary.png)


### BWA
​​
> ## BWA Source Code Installation
>
> If you prefer to install from source, follow the instructions below:
>
> ~~~
> $ cd ~/bch709/bin
> $ curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
> $ tar xvf bwa-0.7.17.tar.bz2
> $ cd bwa-0.7.17
> $ make
> ~~~
> {: .bash}
{: .solution}
​
### Please check README.md
```
less <FILENAME>
```
#### Tell me your error


>## How to solve it?
>~~~
>sudo apt install zlib1g-dev
>~~~
> {: .bash}
>~~~
>brew install zlib
>~~~
> {: .bash}
{: .solution}

**Test your installation by running:**
​
~~~
$ ./bwa
~~~
{: .bash}

### BWA Github Installation
***Seach BWA on Google***


### How to run BWA?
```
bwa index <YOUR_GENOME_SEQUENCE>
bwa mem  <YOUR_GENOME_SEQUENCE> <SEQUENCING_READS>
```


## Conda?
- Dependencies is one of the main reasons to use Conda.
Sometimes, install a package is not as straight forward as you think. Imagine a case like this: You want to install package Matplotlib, when installing, it asks you to install Numpy, and Scipy, because Matplotlib need these Numpy and Scipy to work. They are called the dependencies of Matplotlib. For Numpy and Scipy, they may have their own dependencies. These require even more packages.
 
- Conda provide a solution for this situation: when you install package Matplotlib, it will automatically install all the dependencies like Numpy and Scipy. So you don’t have to install them one by one, manually. This can save you great amount of time.
 
- The other advantage of conda, is that conda can have multiple environments for different projects. As mentioned at the very beginning, it can have two separate environments of different versions of software.
Using conda environment on BioHPC

### Installing Packages Using Conda
>
>Conda is a package manager, which helps you find and install packages such as numpy or scipy. It also serves as an environment manager, and allows you to have multiple isolated environments for different projects on a single machine. Each environment has its own installation directories, that doesn’t share packages with other environments.
>
>For example, you need python 2.7 and Biopython 1.60 in project A, while you also work on another project B, which needs python 3.5 and Biopython 1.68. You can use conda to create two separate environments for each project, and you can switch between different versions of packages easily to run your project code.
{: .callout}

### Anaconda or Miniconda?  
- Anaconda includes both Python and conda, and additionally bundles a suite of other pre-installed packages geared toward scientific computing. Because of the size of this bundle, expect the installation to consume several gigabytes of disk space.

- Miniconda gives you the Python interpreter itself, along with a command-line tool called conda which operates as a cross-platform package manager geared toward Python packages, similar in spirit to the apt or yum tools that Linux users might be familiar with.

### Miniconda3

Miniconda is a package manager that simplifies the installation process. Please first install miniconda3 (installation instructions below), and then proceed to the installation of individual tools. 
​​
### Install Miniconda
Visit the [miniconda](https://docs.conda.io/en/latest/miniconda.html) page and get ready to download the installer of your choice/system.

There are several different env in this world.
***https://repo.anaconda.com/miniconda/***
![package]({{site.baseurl}}/fig/slide_package.png)

> ## MacOS
> 
>To install miniconda3, type:
>
>~~~
>$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
>$ bash Miniconda3-latest-MacOSX-x86_64.sh
>~~~
>{: .bash}
> Then, follow the instructions that you are prompted with on the screen to install Miniconda3. 
{: .solution}

> ## Ubuntu
> 
>To install miniconda3, type:
>
>~~~
>$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
>$ bash Miniconda3-latest-Linux-x86_64.sh
>~~~
>{: .bash}
> Then, follow the instructions that you are prompted with on the screen to install Miniconda3. 
{: .solution}

![conda1]({{site.baseurl}}/fig/conda_excute.png)
![conda2]({{site.baseurl}}/fig/conda_excute2.png)


### Reload your enviroment
#### Ubuntu
```bash
$ source ~/.bashrc
```

#### Ubuntu
```bash
$ source ~/.bash_profile
```
​
### Initialize Miniconda3

```bash
$ conda init
```



### Create new conda environment

#### Create a conda environment named test with latest anaconda package.
```bash 
$ conda create -n bch709 python=3
```
#### Alternatively you can specify python version
```bash
$ conda create -n snowdeer_env python=2.7.16
```

*Usually, the conda environment is installed in your home directory on computer, /home/\<your username\>. The newly created environment <test> will be installed in the directory /home/wyim/miniconda3/envs/test*
 


#### Use the environment you just created
Activate your environment:
```bash  
$ conda activate bch709
```
It will show your environment name at the beginning of the prompt.

![conda3]({{site.baseurl}}/fig/conda.png)

### Install packages in the conda environment

Install from default conda channel
You can search if your package is in the default source from Anaconda collection. Besides the 200 pre-built Anaconda packages, it contains over 600 extra scientific and analytic packages. All the dependencies will also be installed automatically.
``` 
$ conda search <package>
$ conda install <package>
```
### Install from conda-forge channel (example: hisat2)
Conda channels are the remote repository that conda takes to search or download the packages. If you want to install a package that is not in the default Anaconda channel, you can tell conda which channel containing the package, so that conda can find and install.
Conda-forge is a GitHub community-led conda channel, containing general packages which are not in the default Anaconda source. All the packages from conda-forge is listed at https://bioconda.github.io/conda-recipe_index.html

```bash
$ conda search hisat2
$ conda search -c bioconda hisat2
```

### Install from bioconda channel (example: hisat2)
Bioconda is another channel of conda, focusing on bioinformatics software. Instead of adding “-c” to search a channel only one time, “add channels” tells Conda to always search in this channel, so you don’t need to specify the channel every time. Remember to add channel in this order, so that Bioconda channel has the highest priority. Channel orders will be explained in next part.
```bash
 $ conda config --add channels conda-forge
 $ conda config --add channels defaults
 $ conda config --add channels r
 $ conda config --add channels bioconda
```
Adding channels will not generate any command line output.
Then, you can install Stringtie from the Bioconda channel
```bash   
 $ conda install hisat2
```
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/hisat2/README.html)


All the bioconda packages can be found here: https://bioconda.github.io/conda-recipe_index.html


### Install R and R packages
The Conda package manager is not limited to Python. R and R packages are well supported by a conda channel maintained by the developers of Conda. The R-essentials include all the popular R packages with all of their dependencies. The command below opens R channel by “-c r”, and then install the r-essentials using R channel.
```bash
$ conda install -c r r-essentials
```

#### Update R packages
```bash 
$ conda update -c r r-essentials
$ conda update -c r r-<package name>
```

### More conda commands:
#### See all available environments
You can check the list of all separate environments, and it will show * at your current environment. In the figure below, it shows root, since I’m not in any conda environment.
```bash   
$ conda env list
```

#### List all package installed
This will show all the packages and versions you’ve installed.
```bash
$ conda list
```
#### Update packages or conda itself
This will update to the newest version of one package, or conda itself.
update package
```bash
$ conda update <package>
```
update package in env
```bash
conda update  --name <ENV_name> <package>
```
update conda itself
```bash
$ conda update -n test --all
$ conda update -n bch709 --all
```

#### Uninstall package from the environment
```   
$ conda uninstall <package name>
```

#### Exit current environment:
You can exit, when you finish your work in the current environment.
```bash   
$ conda deactivate
```

#### Remove environment
```bash   
$ conda env remove --name bch709
```
When you finish your project, you might want to remove the environment. However, it is not recommended because you might want to update some work in this project in the future.

#### Enviroment export
```bash
conda env export  --name <ENVIRONMENT> --file <outputfilename>.yaml
```
#### Envrioment import 
```bash
conda env create --file <outputfilename>.yaml  
``` 


### HPC clusters
>This exercise mainly deals with using HPC clusters for large scale data (Next Generation Sequencing analysis, Genome annotation, evolutionary studies etc.). These clusters have several processors with large amounts of RAM (compared to typical desktop/laptop), which makes it ideal for running programs that are computationally intensive. The operating system of these clusters are primarily UNIX and are mainly operated via command line. All the commands that you have learned in the previous exercises can be used on HPC.
>
>Pronghorn High Performance Computing offers shared cluster computing infrastructure for researchers and students at UNR. Brief descriptions for the available resources can be found here: https://www.unr.edu/research-computing/hpc. To begin with, you need to request permission for accessing these resources either through your department or through your advisor. All workshop attendees will have their account setup on HPC class education cluster and they can use their UNR NetID and the password for logging-in. You should have already received a confirmation email about your account creation with instructions on how to connect to the cluster. In this exercise we will specifically teach you how to connect to a remote server (HPC), transfer files in and out of the server, and running programs by requesting resources.
>
>
You can log onto its front-end/job-submission system (pronghorn.rc.unr.edu) using your UNR NetID and password. Logging into HPC class requires an SSH client if you are >using Windows but Mac/Linux have these built into their OS. There are several available for download for the Windows platform.
>
>```
>ssh <YOURID>@pronghorn.rc.unr.edu
>```
{: .prereq}


There are a number of ways to transfer data to and from HPC clusters. Which you should use depends on several factors, including the ease of use for you personally, connection speed and bandwidth, and the size and number of files which you intend to transfer. Most common options include scp, rsync (command line) and SCP and SFTP clients (GUI). scp (secure copy) is a simple way of transferring files between two machines that use the SSH (Secure SHell) protocol. You may use scp to connect to any system where you have SSH (login) access. scp is available as a protocol choice in some graphical file transfer programs and also as a command line program on most Linux, UNIX, and Mac OS X systems. scp can copy single files, but will also recursively copy directory contents if given a directory name. scp can be used as follows:

- to a remote system from local
```
scp sourcefile username@pronghorn.rc.unr.edu:somedirectory/
```
- from a remote system to local
```
scp username@pronghorn.rc.unr.edu:somedirectory/sourcefile destinationfile
```
- recursive directory copy to a remote system from local
```
scp -r SourceDirectory/ username@pronghorn.rc.unr.edu:somedirectory/
```

***rsync*** is a fast and extraordinarily versatile file copying tool. It can synchronize file trees across local disks, directories or across a network

- Synchronize a local directory with the remote server directory
```
rsync -avhP path/to/SourceDirectory username@pronghorn.rc.unr.eduu:somedirectory/
```

- Synchronize a remote directory with the local directory
```
rsync -avhP username@hpronghorn.rc.unr.edu:SourceDirectory/ path/to/Destination/
```



>## HOME WORK
>1. Create your enviroment name bch709, install following tools
>- bwa
>- trinity
>- hisat2
>- samtools
>- bowtie2
>- canu
>- trim-galore 
>- multiqc 
>2. Then export your enviroment to bch.yaml (conda env export ......)
>3. copy contents of bch.yaml to [this site](https://docs.google.com/forms/d/e/1FAIpQLSe0faV4UHKEQ8CnqJ-iXOTAmKM2IxN6g7ZNSSmtcw5FuPwmWA/viewform?usp=sf_link).
{: .solution}

### Reference:

- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

