---
layout: page
title: Compile and Software installation
published: true
---

{% include gh_variables.html %}


![software_compile]({{{site.baseurl}}/fig/software-compiler.png)
![software_compile](https://pbs.twimg.com/media/CYIT_SJWQAIExU8.png)

## macOS

> ## Install Xcode
>Inside the Terminal window, copy and paste (or type) the following command, and press the return key on your keyboard:
>```
>$ xcode-select --install
>```
{: .prereq}

> ## Install Homebrew
>```bash
>$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> ```
{: .prereq}

> ## Brew update
>```bash
>$ brew update
>```
{: .prereq}

### Install Prerequire Software
```bash
$ brew install openssl readline sqlite3 xz wget
```

## Ubuntu on Windows
```bash
$ sudo apt update
$ sudo apt-get install -y make build-essential libssl-dev  libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl vim debootstrap
```

## Let's install test package!
### On Ubuntu systems:  
```bash
$ apt install screenfetch
```

### On macOS systems:  
```bash
$ brew install screenfetch
```

```bash
screenfetch
```

### Programming languages

|-----|------|
|Compiled|FORTRAN, C, C++, Java|
|Interpretive|Unix-Shell, awk, Basic, Perl, Tcl, Scheme, Ruby, Python, R|



#### **Perl**
Flexible, by a global repository (CPAN), thus it is small install new modules. It has Bio per, one of the first biological unit repositories that upsurge the usability from, for instance, change setups to do the phylogenetic investigation. There are several biological software those usages Perl such as GBrowse thus might be an exciting language if you need toward interact with it. Good test units Perl is still what a lot of persons use, but it is declining out of use since Python accomplishes the similar tasks and is easier toward write code for, particularly for beginners.


#### **R & Python**
R is great for all the reasons, but if you like coding more than statistics, you may enjoy Python’s style a lot more. That sounds like a contradiction: How could you possibly know you enjoy coding more than statistics when you are choosing your first programming language? I would suggest trying them both and seeing what you like best. I personally enjoy coding in Python more than in R because its rules make more sense and it feels more like a programming language. In my experience, it is also much easier to make a command-line tool in Python than in R, and Python also has some packages for bioinformatics that are quite useful.


#### **Bash**
It is also very important for bioinformaticians to learn Bash, which for all of our intents and purposes is interchangeable with shell, the command-line, or the terminal. Bash is the primary way to access your data on your institution’s cluster and to run most genomics and bioinformatics software. It is also very powerful for manipulating your data like sorting, filtering, or doing calculations between columns, which is available through various utilities.

In my experience, and everyone I have talked to about it, bash was confusing and scary at first, but when you get the hang out it you start to feel this power surging through you, and you can do things in second that would take you hours to do by hand. Even two years into it I would still learn something new in bash that would blow my mind and I would kick myself for wasting time having programmed it from scratch in Python.

### Python, Perl, R and bash
In summary, for wet-lab people who want to add bioinformatics to their toolbox, focus on learning R first and applying it to your own work. For people who want to focus on bioinformatics as a career and make their own tools too, I would actually recommend learning the trifecta of R, Python, and Bash, though you could get away with choosing between R and Python as long as you still learn Bash too. I can go into more depth on any of these topics or give an introduction to any of these languages if you let me know in the comments.


### Other programming languages
There are many other languages out there, so before I end here I’m going to give a brief reason why these are not recommended for bioinformatics, beginners, or anyone at all in some cases.

#### C and C++
C or C++ are great for making super optimized command-line tools like aligners and variant-callers, but you will have a much easier time learning Python first and then going to these high-performance languages for a particular problem in the future, since they are harder to learn, more finicky, and take a lot more code to do the same thing.

#### Ruby
Ruby is one of those hot languages right now, for good reason largely because of the power of Ruby on Rails for making database-driven web applications like blogs or twitter. Ruby however is not great for bioinformatics because it lacks the community support in terms of packages that R and Python have, so you would be better off learning Python instead of Ruby.

#### JavaScript or PHP
JavaScript and PHP are great languages for web applications, but bioinformatics web applications should never be your first project. You could make a computational method in Python or R and then later make it into a web application, but that is not a project for a beginner. HTML and CSS by the way are not programming languages, but actually markup and styling languages that you will use along with JavaScript and PHP for that web application someday.

#### Java
Java is a popular language that most people have heard of. In bioinformatics, a notable example is the genome browser IGV. However, I would not recommend for beginners to learn Java due to many issues including memory management and that Python and R have many more bioinformaticians who build packages and answer questions online.


![language]({{{site.baseurl}}/fig/language.png)


### Package Library Module

#### Library 
Most often will refer to the general library or another collection created with a similar format and use. The General Library is the sum of 'standard', popular and widely used Modules, witch can be thought of as single file tools, for now or short cuts making things possible or faster. The general library is an option most people enable when installing Python. Because it has this name "Python General Library" it is used often with similar structure, and ideas. Witch is simply to have a bunch of Modules, maybe even packages grouped together, usually in a list. The list is usually to download them. Generally it is just related files, with similar interests. That is the easiest way to describe it.

#### Module 
A Module refers to a file. The file has script 'in it' and the name of the file is the name of the module, Python files end with .py. All the file contains is code that ran together makes something happen, by using functions, strings ect. Main modules you probably see most often are popular because they are special modules that can get info from other files/modules. It is confusing because the name of the file and module are equal and just drop the .py. Really it's just code you can use as a shortcut written by somebody to make something easier or possible.

#### Package
This is a termis used to generally sometimes, although context makes a difference. The most common use from my experience is multiple modules (or files) that are grouped together. Why they are grouped together can be for a few reasons, that is when context matters. These are ways I have noticed the term package(s) used. They are a group of Downloaded, created and/or stored modules. Which can all be true, or only 1, but really it is just a file that references other files, that need to be in the correct structure or format, and that entire sum is the package itself, installed or may have been included in the python general library. A package can contain modules(.py files) because they depend on each other and sometimes may not work correctly, or at all. There is always a common goal of every part (module/file) of a package, and the total sum of all of the parts is the package itself.

### Programming languages module or library manager
Python - pip  
Perl - cpan  
R - native manager  

## File Permission

![language]({{{site.baseurl}}/fig/linux_file_permissions.png)

### Understanding of attribute which can be out put by `ls -l`:
![language]({{{site.baseurl}}/fig/file_permission.png)
![language]({{{site.baseurl}}/fig/file_permission2.png)

### Chmod (Change mode)
chmod is the command and system call which is used to change the access permissions of file system objects (files and directories). It is also used to change special mode flags. The request is filtered by the umask. The name is an abbreviation of change mode.

### Applying Permission:
![language]({{{site.baseurl}}/fig/file_permission3.png)
![language]({{{site.baseurl}}/fig/file_permission4.png)

### Using Octal number:
![language]({{{site.baseurl}}/fig/file_permission5.png)


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

### Link Your Directory 
```bash
$ mkdir /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop 
$ mkdir /mnt/c/Users/<YOURID_WINDOWSID>/Desktop/BCH709_Desktop .
```

```bash
$ mkdir ~/Desktop/BCH709
$ ln -s ~/Desktop/BCH709_Desktop .
```

### Check your results
unzip your results


## HISAT2
https://ccb.jhu.edu/software/hisat2/index.shtml


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


### BWA
​​
> ## BWA Source Code Installation
>
> If you prefer to install from source, follow the instructions below:
>
> ~~~
> $ cd ~/src
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
$ bwa
~~~
{: .bash}


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

### Install from bioconda channel (example: stringtie)
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
>2. Then export your enviroment to bch.yaml
>3. copy contents of bch.yaml to [this site](https://docs.google.com/forms/d/e/1FAIpQLSe0faV4UHKEQ8CnqJ-iXOTAmKM2IxN6g7ZNSSmtcw5FuPwmWA/viewform?usp=sf_link).
{: .solution}

### Reference:

- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

