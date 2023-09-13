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


## Check your CPUs and Memory

```
$ lscpu
$ free
$ htop
```
### RC file such as .bashrc, .zshrc

In Unix-like operating systems such as Linux, rc is a file name extension which stands for "run commands." An rc file contains several statements or commands, one per line, to be executed or evaluated by the shell. Its purpose is to configure the environment and prepare the system to run specific software.

### Connect to Pronghorn
You can log onto its front-end/job-submission system (pronghorn.rc.unr.edu) using your UNR NetID and password. Logging into HPC class requires an SSH client if you are >using Windows but Mac/Linux have these built into their OS. There are several available for download for the Windows platform.

```bash
ssh <YOURID>@pronghorn.rc.unr.edu
```


### Prompt Customization for Pronghorn and Ubuntu

```bash
echo '###BCH709 ' >> ~/.bashrc

echo 'tty -s && export PS1="\[\033[38;5;164m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;231m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;172m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"' >> ~/.bashrc
echo "alias ls='ls --color=auto'" >> ~/.bashrc

```


```bash
source ~/.bashrc
```

### Prompt Customization for Mac
```bash
sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"
```

```bash
source ~/.zshrc
```

#### *More information is* [here](https://plantgenomicslab.github.io/BCH709/bash/index.html)


## Linux family tree
[Linux family tree](https://en.wikipedia.org/wiki/List_of_Linux_distributions)

> ## Package Management Concepts
>![package](https://community-cdn-digitalocean-com.global.ssl.fastly.net/assets/tutorials/images/large/Package_Management_tw_mostov.png)
> Contemporary distributions of Linux-based operating systems install software in pre-compiled packages, which are archives that contain binaries of software, configuration files, and information about dependencies. Furthermore, package management tools keep track of updates and upgrades so that the user doesn’t have to hunt down information about bug and security fixes.
>
> Without package management, users must ensure that all of the required dependencies for a piece of software are installed and up-to-date, compile the software from the source code (which takes time and introduces compiler-based variations from system to system), and manage configuration for each piece of software. Without package management, application files are located in the standard locations for the system to which the developers are accustomed, regardless of which system they’re using.
>
> Package management systems attempt to solve these problems and are the tools through which developers attempt to increase the overall quality and coherence of a Linux-based operating system. The features that most package management applications provide are:
>
> - Package downloading: Operating-system projects provide package repositories which allow users to download their packages from a single, trusted provider. When you download from a package manager, the software can be authenticated and will remain in the repository even if the original source becomes unreliable.
> - Dependency resolution: Packages contain metadata which provides information about what other files are required by each respective package. This allows applications and their dependencies to be installed with one command, and for programs to rely on common, shared libraries, reducing bulk and allowing the operating system to manage updates to the packages.
> - A standard binary package format: Packages are uniformly prepared across the system to make installation easier. While some distributions share formats, compatibility issues between similarly formatted packages for different operating systems can occur.
> - Common installation and configuration locations: Linux distribution developers often have conventions for how applications are configured and the layout of files in the /etc/ and /etc/init.d/ directories; by using packages, distributions are able to enforce a single standard.
> - Additional system-related configuration and functionality: Occasionally, operating system developers will develop patches and helper scripts for their software which get distributed within the packages. These modifications can have a significant impact on user experience.
> - Quality control: Operating-system developers use the packaging process to test and ensure that the software is stable and free of bugs that might affect product quality and that the software doesn’t cause the system to become unstable. The subjective judgments and community standards that guide packaging and package management also guide the “feel” and “stability” of a given system.
> In general, we recommend that you install the versions of software available in your distribution’s repository and packaged for your operating system. If packages for the application or software that you need to install aren’t available, we recommend that you find packages for your operating system, when available, before installing from source code.
>
> The remainder of this guide will cover how to use specific package management systems and how to compile and package software yourself.
{: .prereq}  

> ## Advanced Packaging Tool (APT)
> You may already be familiar with apt-get, a command which uses the advanced packaging tool to interact with the operating system’s package system. The most relevant and useful commands are (to be run with root privileges):
>
> - 'apt-get install package-name(s)' - Installs the package(s) specified, along with any dependencies.
> - 'apt-get remove package-name(s)' - Removes the package(s) specified, but does not remove dependencies.
> - 'apt-get autoremove' - Removes any orphaned dependencies, meaning those that remain installed but are no longer required.
> - 'apt-get clean' - Removes downloaded package files (.deb) for software that is already installed.
> - 'apt-get purge package-name(s)' - Combines the functions of remove and clean for a specific package, as well as configuration files.
> - 'apt-get update' - Reads the /etc/apt/sources.list file and updates the system’s database of packages available for installation. Run this after changing sources.list.
> - 'apt-get upgrade' - Upgrades all packages if there are updates available. Run this after running apt-get update.
> While apt-get provides the most often-used functionality, APT provides additional information in the apt-cache command.
> 
> - 'apt-cache search package-name(s)' - If you know the name of a piece of software but apt-get install fails or points to the wrong software, this looks for other possible names.
> - 'apt-cache show package-name(s)' - Shows dependency information, version numbers and a basic description of the package.
> - 'apt-cache depends package-name(s)' - Lists the packages that the specified packages depends upon in a tree. These are the packages that will be installed with the apt-get install command.
> - 'apt-cache rdepends package-name(s)' - Outputs a list of packages that depend upon the specified package. This list can often be rather long, so it is best to pipe its output through a command, like less.
> - 'apt-cache pkgnames' - Generates a list of the currently installed packages on your system. This list is often rather long, so it is best to pipe its output through a program, like less, or direct the output to a text file.
> Combining most of these commands with apt-cache show can provide you with a lot of useful information about your system, the software that you might want to install, and the software that you have already installed.
{: .callout}

> ## Aptitude
> Aptitude is another front-end interface for APT. In addition to a graphical interface, Aptitude provides a combined command-line interface for most APT functionality. Some notable commands are:
{: .callout}

> ## Using dpkg
> Apt-get and apt-cache are merely frontend programs that provide a more usable interface and connections to repositories for the underlying package management tools called dpkg and debconf. These tools are quite powerful, and fully explaining their functionality is beyond the scope of this document. However, a basic understanding of how to use these tools is useful. Some important commands are:
> 
> - 'dpkg -i package-file-name.deb' - Installs a .deb file.
> - 'dpkg --list search-pattern' - Lists packages currently installed on the system.
> - 'dpkg --configure package-name(s)' - Runs a configuration interface to set up a package.
> - 'dpkg-reconfigure package-name(s)' - Runs a configuration interface on an already installed package
{: .callout}

> ## Fedora and CentOS Package Management
> Fedora and CentOS are closely related distributions, being upstream and downstream (respectively) from Red Hat Enterprise Linux (RHEL). Their main differences stem from how packages are chosen for inclusion in their repositories.
>
> CentOS uses yum, Yellowdog Updater, Modified, as a front end to interact with system repositories and install dependencies, and also includes a lower-level tool called rpm, which allows you to interact with individual packages.
>
> Starting with version 22, Fedora uses the dnf package manager instead of YUM to interact with rpm. DNF supports many of the same commands as YUM, with some slight changes.
>
> Note: Many operating systems aside from RedHat use rpm packages. These include OpenSuSE, AIX, and Mandriva. While it may be possible to install an RPM packaged for one operating system on another, this is not supported or recommended, and the results of this action can vary greatly.
{: .callout}

> ## How about macOS?
> Homebrew is package manager for Macs which makes installing lots of different software like Git, Ruby, and Node simpler. Homebrew lets you avoid possible security problems associated with using the sudo command to install software like Node.
> Homebrew has made extensive use of GitHub to expand the support of several packages through user contributions. In 2010, Homebrew was the third-most-forked repository on GitHub. In 2012, Homebrew had the largest number of new contributors on GitHub. In 2013, Homebrew had both the largest number of contributors and issues closed of any project on GitHub.
> Homebrew has spawned several sub-projects such as Linuxbrew, a Linux port now officially merged into Homebrew; Homebrew Cask, which builds upon Homebrew and focuses on the installation of GUI applications and "taps" dedicated to specific areas or programming languages like PHP.
{: .prereq}

> ## macOS Requirements
> A 64-bit Intel CPU 
> macOS 10.12 (or higher)
> [Command Line Tools (CLT) for Xcode](https://plantgenomicslab.github.io/BCH709/CLT/index.html)
> A Bourne-compatible shell for installation (e.g. bash or zsh)
{: .prereq}  

[![homebrew](https://brew.sh/assets/img/homebrew-256x256.png)](https://brew.sh/)


> ## Install Homebrew
> ```bash
> /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> ```
{: .prereq}


> ## How to Update Homebrew
> New versions of Homebrew come out frequently, so make sure you update it before updating any of the other software components that you’ve installed using Homebrew. * In Terminal type ```brew update```
{: .callout}

> ## How to Uninstall Homebrew
> Open the Terminal app
> Type ```ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)```
> This downloads and runs the uninstaller script. Follow the instructions and Homebrew will be removed from your computer.
{: .callout}

## List installed packages
---
We can get a list of all the installed packages on a Debian / Ubuntu server by issuing:

```bash
$ sudo dpkg --get-selections
```
Ubuntu server by issuing:
```bash
$ apt list --installed
```
on macOS
```bash
$ brew list
```

On RPM systems:
```bash
$ yum list installed
```
On BSD systems:
```bash
$ pkg_version
```
It is good practice to save this file as it can be useful when migrating, so we pipe it into a file:
```bash
$ dpkg --get-selections > ~/package_list
 #yum list installed
 #pkg_version
```
To search for a specific package run:
```bash
$ dpkg --get-selections | grep <package>
$ yum list installed "package_name"
```

## Search packages
On Ubuntu systems: 
```
apt search <package-name>
```
```bash
$ apt search firefox
$ apt search ^firefox 
```

On macOS systems:
```
brew search <package-name>
```
```bash
$ brew search firefox
$ brew search /^firefox/
```
\^ means regular expressions start of the line.

## Install packages
### Install single packages:
On Ubuntu systems:
```
apt install <package-name> 
```
On macOS systems:
```
brew install <package-name> 
```
### Install multiple packages:
On Ubuntu systems:
```
apt install <package-name> <package-name> ...
```
On macOS systems:
```
brew install <package-name> <package-name> ...
```
## Install specific version
### Search version
On Ubuntu systems:
```bash
$ apt-cache policy <package-name> 
```
On macOS systems:
```bash
$ brew search <package-name>
```
### Install specific version
On Ubuntu systems:
```bash
$ apt install firefox=68.0.1+build1-0ubuntu0.18.04.1
```
On macOS systems:
```bash
$ brew install firefox@68.0.2
```


## Software

| Software | Version | Manual | Available for | Description |
| -------- | ------------ | ------ | ------------- | ----------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.7 | [Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)| Linux, MacOS, Windows | Quality control tool for high throughput sequence data. |
| [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) | 2.1.0| [Link](https://ccb.jhu.edu/software/hisat2/index.shtml) | Linux, MacOS, Windows | Mapping RNA sequences against genome |
| [BWA](http://bio-bwa.sourceforge.net/) | 0.7.17 | [Link](http://bio-bwa.sourceforge.net/bwa.shtml) | Linux, MacOS | Mapping DNA sequences against reference genome. |


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

### Install Miniconda
Visit the [miniconda](https://docs.conda.io/en/latest/miniconda.html) page and get ready to download the installer of your choice/system.

There are several different env in this world.
***https://repo.anaconda.com/miniconda/***
![package]({{site.baseurl}}/fig/slide_package.png)

> ## Linux
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


## Reload your Conda enviroment
### Linux
```bash
$ source ~/.bashrc
```

### Mac OS
```bash
$ source ~/.bash_profile ## OLD MAC
$ source ~/.zshrc  ## NEW MAC
```

### Initialize Miniconda3

```bash
$ conda init
```

### Create new conda environment

#### Create a conda environment named test with latest anaconda package.
```bash 
$ conda create -n bch709 
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
```bash
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


Adding channels will not generate any command line output.
Then, you can install Stringtie from the Bioconda channel
```bash   
 $ conda install -c bioconda hisat2
```
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/hisat2/README.html)


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


### More conda commands:.
#### search packages
This will search the packages in conda
```bash
$ conda search hisat2
```

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
```bash
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

### Reference:

- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

