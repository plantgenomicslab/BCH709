---
layout: page
title: Compile and Software installation
published: true
---

{% include gh_variables.html %}


![software_compile]({{{site.baseurl}}/fig/software-compiler.png)
![software_compile](https://pbs.twimg.com/media/CYIT_SJWQAIExU8.png)

## macOS

> ### Install Homebrew
> ```bash
> $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> ```
{: .prereq}


> ## Install Homebrew
>```bash
>$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> ```
{: .prereq}

> ### Brew update
> ```bash
> $ brew update
> ```
{: .prereq}

### Install Prerequisite Software
```bash
$ brew install openssl readline sqlite3 xz wget
```

## Ubuntu on Windows
```bash
$ sudo apt update
$ sudo apt-get install -y make build-essential libssl-dev  libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl vim debootstrap
```

## Install test package!
### On Ubuntu systems:  
```bash
$ apt install screenfetch
$ screenfetch

```

### On macOS systems:  
```bash
$ brew install screenfetch
$ screenfetch

```


### Programming languages

|-----|------|
|Compiled|FORTRAN, C, C++, Java|
|Interpretive|Unix-Shell, awk, Basic, Perl, Tcl, Scheme, Ruby, Python, R|



#### **Perl**
Perl is flexible and has a global repository (CPAN), which makes it easy to install new modules. It also has BioPerl, one of the first biological unit repositories, enhancing usability for tasks such as phylogenetic analysis. Although Perl was widely used, Python has become more popular due to its ease of use, especially for beginners.



#### **R & Python**
R is excellent for statistical analysis, but if you prefer coding, Python might suit you more. Python's rules are easier to follow, making it more beginner-friendly. It's also easier to develop command-line tools in Python, and there are useful bioinformatics packages available in Python.



#### **Bash**
Bash (or shell scripting) is essential for bioinformaticians. It’s a powerful tool for data manipulation (sorting, filtering, etc.) and is often used on institutional clusters. It may seem intimidating at first, but with time, you will find it very efficient for repetitive tasks and system administration.


### Python, Perl, R and bash
For wet-lab researchers starting with bioinformatics, R is a good choice to learn first. If you aim for a bioinformatics career, knowing R, Python, and Bash is recommended. For beginners, focusing on either R or Python, while learning Bash, can still be effective.


### Other programming languages

#### C and C++
C and C++ are great for high-performance tools like aligners, but they are harder to learn and take more code to accomplish tasks that can be done more simply in Python.


#### Ruby
Ruby is popular for web applications but lacks the package support for bioinformatics that Python and R have.


#### JavaScript or PHP
These languages are better suited for web applications. Bioinformatics should start with Python or R before considering web development languages.


#### Java
Java has some uses in bioinformatics (e.g., IGV genome browser), but it's not beginner-friendly, especially when compared to Python or R.



![language]({{{site.baseurl}}/fig/language.png)


### Package Library Module

#### Library 
Refers to a collection of related packages or modules. It’s often used to describe the Python Standard Library, which contains many modules that provide additional functionality.


#### Module 
A module is a single file containing Python code. When you import a module, Python executes the code inside that file. Modules make code more reusable and manageable.


#### Package
A package is a collection of related modules that work together. It usually includes several files and folders organized in a specific structure.


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

### Using Octal number for Permissions:
![language]({{{site.baseurl}}/fig/file_permission5.png)


## Check your CPUs and Memory

```
$ lscpu
$ free
$ htop
```
### RC file such as .bashrc, .zshrc

RC files configure the environment and prepare the system to run specific software. These are commonly used in Unix-like systems to automate shell configurations.


### Connect to Pronghorn
You can log onto its front-end/job-submission system (pronghorn.rc.unr.edu) using your UNR NetID and password. Logging into HPC class requires an SSH client if you are >using Windows but Mac/Linux have these built into their OS. There are several available for download for the Windows platform.

```bash
ssh <YOURID>@pronghorn.rc.unr.edu
```


### Prompt Customization for Windows
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
source ~/.zshrc

```

#### *More information is here* [here](https://plantgenomicslab.github.io/BCH709/bash/index.html)


## Linux family tree
[Linux family tree](https://en.wikipedia.org/wiki/List_of_Linux_distributions)

> ## Package Management Concepts
>![package](https://community-cdn-digitalocean-com.global.ssl.fastly.net/assets/tutorials/images/large/Package_Management_tw_mostov.png)
> Package management in Linux allows for easier installation and updating of software. It handles dependencies and ensures proper installation across systems. Popular tools include APT (Debian/Ubuntu), YUM (CentOS/Fedora), and Homebrew (macOS).
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
dpkg --get-selections | grep <package>
yum list installed "package_name"
```

## Search packages
On Ubuntu systems: 
```
apt search <package-name>
```
```bash
apt search firefox
apt search ^firefox 
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
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.7 | [Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)| Linux, MacOS, Windows | Quality control tool for high throughput sequence data. |
| [HISAT2](https://daehwankimlab.github.io/hisat2/) | 2.1.0| [Link](https://daehwankimlab.github.io/hisat2/) | Linux, MacOS, Windows | Mapping RNA sequences against genome |
| [BWA](http://bio-bwa.sourceforge.net/) | 0.7.17 | [Link](http://bio-bwa.sourceforge.net/bwa.shtml) | Linux, MacOS | Mapping DNA sequences against reference genome. |


## Conda?
Conda helps manage package dependencies and environments, making it easier to install packages and maintain reproducibility.

- Dependencies is one of the main reasons to use Conda.
Sometimes, install a package is not as straight forward as you think. Imagine a case like this: You want to install package Matplotlib, when installing, it asks you to install Numpy, and Scipy, because Matplotlib need these Numpy and Scipy to work. They are called the dependencies of Matplotlib. For Numpy and Scipy, they may have their own dependencies. These require even more packages.
 
- Conda provide a solution for this situation: when you install package Matplotlib, it will automatically install all the dependencies like Numpy and Scipy. So you don’t have to install them one by one, manually. This can save you great amount of time.
 
- The other advantage of conda, is that conda can have multiple environments for different projects. As mentioned at the very beginning, it can have two separate environments of different versions of software.
Using conda environment on BioHPC

### Installing Packages Using Miniconda
>
>Conda is a package manager, which helps you find and install packages such as numpy or scipy. It also serves as an environment manager, and allows you to have multiple isolated environments for different projects on a single machine. Each environment has its own installation directories, that doesn’t share packages with other environments.
>
>For example, you need python 2.7 and Biopython 1.60 in project A, while you also work on another project B, which needs python 3.5 and Biopython 1.68. You can use conda to create two separate environments for each project, and you can switch between different versions of packages easily to run your project code.
{: .callout}

### Anaconda or Miniconda?  
- Anaconda includes both Python and conda, and additionally bundles a suite of other pre-installed packages geared toward scientific computing. Because of the size of this bundle, expect the installation to consume several gigabytes of disk space.

- Miniconda gives you the Python interpreter itself, along with a command-line tool called conda which operates as a cross-platform package manager geared toward Python packages, similar in spirit to the apt or yum tools that Linux users might be familiar with.

Here is the rewritten content in markdown source code format:

```markdown
### Miniconda3

Miniconda is a lightweight package manager that simplifies the installation process for environments and packages. First, install Miniconda3 by following the instructions below, then proceed to install individual tools.

### Install Miniconda

Visit the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) page to download the installer for your system.

For direct access, visit:  
***https://repo.anaconda.com/miniconda/***

![package]({{site.baseurl}}/fig/slide_package.png)

> ## Linux
> 
> To install Miniconda3 on Linux, run:
> 
> ~~~bash
> mkdir -p ~/miniconda3
> wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
> bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
> rm ~/miniconda3/miniconda.sh
> ~~~
> 
> Follow the on-screen instructions to complete the installation.
{: .solution}

![conda1]({{site.baseurl}}/fig/conda_excute.png)
![conda2]({{site.baseurl}}/fig/conda_excute2.png)

> ## macOS Intel
> 
> To install Miniconda3 on macOS with Intel architecture, run:
> 
> ~~~bash
> mkdir -p ~/miniconda3
> curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda3/miniconda.sh
> bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
> rm ~/miniconda3/miniconda.sh
> ~~~
> 
> Follow the on-screen instructions to complete the installation.
{: .solution}

> ## macOS M1 (Apple Silicon)
> 
> To install Miniconda3 on macOS with M1 (Apple Silicon) architecture, run:
> 
> ~~~bash
> mkdir -p ~/miniconda3
> curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
> bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
> rm ~/miniconda3/miniconda.sh
> ~~~
> 
> Follow the on-screen instructions to complete the installation.
{: .solution}

## Reload Your Conda Environment

### Linux
```bash
source ~/.bashrc
```

### macOS
```bash
source ~/.bash_profile ## For older macOS versions
source ~/.zshrc         ## For newer macOS versions
```

### Initialize Miniconda3

To initialize Conda after installation, run:

```bash
$ conda init
```


### Creating and Using Conda Environments

To create a new Conda environment with Python 3.8 and activate it, use the following commands:

```bash 
$ conda create -n bch709 python=3.8
$ conda activate bch709
```

### Installing Packages in Conda

You can install packages in your active environment using:

```bash
$ conda install <package-name>
```

*By default, Conda environments are installed in your home directory, typically under `/home/<your_username>/miniconda3/envs/<environment_name>`. For example, the environment `bch709` would be installed in `/home/<your_username>/miniconda3/envs/bch709`.*

### Deactivating and Removing Environments

To deactivate the current environment, run:

```bash
$ conda deactivate
```

To remove the environment `bch709`, use:

```bash
$ conda env remove --name bch709
```

### Activating the Environment

Once created, activate your environment using:

```bash  
$ conda activate bch709
```

You will see the environment name in the prompt.

![conda3]({{site.baseurl}}/fig/conda.png)

### Installing Packages from Conda Channels

#### Installing from Default Conda Channel

Search for a package in the default Anaconda repository:

```bash
$ conda search <package>
```

Install the package:

```bash
$ conda install <package>
```

#### Installing from Conda-Forge Channel (Example: HISAT2)

Conda channels are remote repositories that contain packages. To install a package from Conda-Forge:

```bash
$ conda search hisat2
$ conda search -c conda-forge hisat2
```

#### Installing from Bioconda Channel (Example: HISAT2)

Bioconda is a channel dedicated to bioinformatics software. To ensure Bioconda has the highest priority, add the channel order:

```bash
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

Now, install HISAT2 from Bioconda:

```bash   
$ conda install -c bioconda hisat2
```

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/hisat2/README.html)

### Installing R and R Packages

Conda also supports R and R packages. You can install the R-essentials package (a collection of popular R packages) using:

```bash
$ conda install -c r r-essentials
```

#### Updating R Packages

To update R packages, run:

```bash 
$ conda update -c r r-essentials
$ conda update -c r r-<package-name>
```

### Additional Conda Commands

#### Search Packages

To search for a package:

```bash
$ conda search hisat2
```

#### List Available Environments

To view all available environments:

```bash   
$ conda env list
```

#### List Installed Packages

To list all installed packages in the current environment:

```bash
$ conda list
```

#### Updating Packages or Conda Itself

To update a specific package:

```bash
$ conda update <package>
```

To update a package in a specific environment:

```bash
$ conda update --name <ENV_name> <package>
```

To update Conda itself:

```bash
$ conda update -n test --all
$ conda update -n bch709 --all
```

#### Uninstalling a Package

To uninstall a package from the environment:

```bash
$ conda uninstall <package-name>
```

#### Exiting the Current Environment

To exit the active environment:

```bash   
$ conda deactivate
```

#### Removing an Environment

To remove an environment:

```bash   
$ conda env remove --name bch709
```

#### Exporting and Importing Environments

To export an environment to a YAML file:

```bash
$ conda env export --name <ENVIRONMENT> --file <outputfilename>.yaml
```

To import an environment from a YAML file:

```bash
$ conda env create --file <outputfilename>.yaml  
```

### References

- Conda documentation: https://docs.conda.io/en/latest/
- Conda-forge: https://conda-forge.github.io/
- BioConda: https://bioconda.github.io/
