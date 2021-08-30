---
layout: page
title: 2_Linux Enviroment and command line
published: true
---

{% include gh_variables.html %}
>## Reading and Slack
>
>- Reading 1: [Luscombe et al., 2001](http://archive.gersteinlab.org/papers/e-print/whatis-mim/text.pdf)  
>- Reading 2: [Attwood 2000](https://science.sciencemag.org/content/290/5491/471)  
>- Reading 3: [Smith 2018](https://www.embopress.org/doi/full/10.15252/embr.201846262)  
{: .callout}


>## Assignments
> Please finish below assignment
> - [DataCamp Introduction to Shell](https://learn.datacamp.com/courses/introduction-to-shell) 

![bioinformatics_DNA]({{site.baseurl}}/fig/DNA.jpg)

## History of UNIX 
The Unix operating system was conceived and implemented in 1969 at **AT&T's Bell Laboratories** in the United States by Ken Thompson, Dennis Ritchie, Douglas McIlroy, and Joe Ossanna. It was first released in 1971 and was initially entirely written in assembly language, a common practice at the time. Later, in a key pioneering approach in 1973, Unix was re-written in the programming language C by Dennis Ritchie (with exceptions to the kernel and I/O). The availability of an operating system written in a high-level language allowed easier portability to different computer platforms. With a legal glitch forcing AT&T to license the operating system's source code to anyone who asked. Unix quickly grew and became widely adopted by academic institutions and businesses. In 1984, AT&T divested itself of Bell Labs. Free of the legal glitch requiring free licensing, Bell Labs began selling Unix as a proprietary product.


## Why UNIX
- Unix is an important OS in the history of computing
- Two major OS variants, Unix-based and Windows-based
- Used in a lot of back-end systems and personal computing
- Unix derivatives like Linux are open source, and well known to the community and developed in the open where we can study and understand them.
- The skills you learn on Unix will easily translate to many other OS platforms because all Unix-based systems share standard characteristic

![Unix_family tree]({{site.baseurl}}/fig/unix-simple.png)

## The kernel
This is called as the hub of the operating system, serving as allocator of time and memory to programs and handling the filestore and communications in response to system calls.

## BSD
BSD (Berkeley Software Distribution) is a version of UNIX developed by folks at UC Berkeley starting from the original Bell Labs code. Several other derivatives have come from BSD, notably FreeBSD, OpenBSD, and NetBSD. OS X (macOS) was brough to Apple from Steve Jobs' previous company NeXT, and was built using pieces from BSD around a kernel called Mach, which incidentally is also the basis of GNU Hurd. Jobs wanted to get a new computer (with a new operating system) to market quickly. To save time, the NeXT team used the Mach kernel from Carnegie Mellon and parts of the BSD code base to created the NeXTSTEP operating system. Also, PS4 is BSD.



## LINUX  
Linux is a Unix-like computer operating system assembled under the model of free and open source software development and distribution. The defining component of Linux is the Linux kernel, an operating system kernel first released 5 October 1991 by Linus Torvalds. Linux was originally developed as a free operating system for Intel x86-based personal computers.

 - He had been built Linux kernel **Just for Fun**, when he originally played with MINIX. Linux means *Linus's MINIX*. He is very famous with **flamming** like flamming in Reddit.

 - "You people on the East coast think you have it bad, with snow-storms and whatever.
That's nothing. My coffee maker broke, and calling the service hotline says "we're not open today due to inclement weather".
You guys get a little snow, and suddenly civilization breaks down.
My coffee maker is broken and nobody is answering the phone.
And CNN just keeps talking about snow. What about my coffee? Priorities, people, priorities.
What am I going to do without my coffee maker? I'm going to sit here in a corner, crying, that's what." --Linus Torvalds
 
 - "Talk is cheap. Show me the code." --Linus Torvalds

 - Git was developed by Linus Torvalds

It has since been ported to more computer hardware platforms than any other operating system. It is a leading operating system on servers and other big iron systems such as mainframe computers and supercomputers:**more than 90% of today's 500 fastest supercomputers** run some variant of Linux,including the 10 fastest. Linux also runs on embedded systems (devices wherthe operating system is typically built into the firmware and highly tailored to the system) such as mobile phones, tablet computers, network routers, televisions and video game consoles; the **Android** system in wide use on mobile devices is built on the Linux kernel. Also, Nintendo Switch (Linux) is based on Linux.

![Lnix_family tree](https://aerojsoft.files.wordpress.com/2016/02/linus-distribution-family-tree.jpg)


[Lnix_family tree](https://en.wikipedia.org/wiki/List_of_Linux_distributions)  


## Unix/Linux Main Components
### The Unix/Linux computer ecosystem can be divided into three main parts:
- **User Space**: Defines the applications, libraries, and standard utilities that are user accessible. When we write a program, it is from this perspective that we operate, without concern for the underlying components. For example, writing a "Hello World" program on any computer is the same from the user-perspective, but might be different when it comes to actually executing the program and writing "Hello World" to the terminal.
- **Kernel Space**: This refers to the operations of OS that manage the interface between user actions and the hardware. It is the central part of the OS, and its primary job is to pair user applications with the underlying hardware and allow multiple programs to share singular hardware components. For example, how does a user input event, such as typing 'a' on the keyboard, get translated into 'a' appearing on the screen? Or, how does two programs both read from disc at the same time or run on the CPU at the same time?  
- **Hardware**: The underlying physical components of the computer. These include Input/Output devices, like keyboards and monitors, the CPU which does calculations, the memory components, and the network interface.

-----------------------

![anatomy]({{site.baseurl}}/fig/anatomy.jpg)

## What is GNU?
GNU is an operating system that is free software—that is, it respects users' freedom. The GNU operating system consists of GNU packages (programs specifically released by the GNU Project) as well as free software released by third parties. The development of GNU made it possible to use a computer without software that would trample your freedom. GNU is a recursive acronym for "GNU's Not Unix!", chosen because GNU's design is Unix-like, but differs from Unix by being free software and containing no Unix code. Richard Matthew Stallman often known by his initials, rms, he started this project. The GNU General Public License (GNU GPL or GPL) is a widely-used free software license, which guarantees end users the freedom to run, study, share and modify the software.


## How about macOS? (XNU)
Some people might think that there are similarities between the macOS and the Linux kernel because they can handle similar commands and similar software. Some people even think that Apple’s macOS is based on Linux. The truth is that both kernels have very different histories and features. 

The macOS kernel is officially known as XNU. The acronym stands for “XNU is Not Unix.” According to Apple’s Github page, XNU is “a hybrid kernel combining the Mach kernel developed at Carnegie Mellon University with components from FreeBSD and C++ API for writing drivers”. The BSD subsystem part of the code is “typically implemented as user-space servers in microkernel systems”. The Mach part is responsible for low-level work, such as multitasking, protected memory, virtual memory management, kernel debugging support, and console I/O.

![GNU](http://www.linuxandubuntu.com/wp-content/uploads/2019/07/What-is-GNU-in-GNULinux.jpg)
GNU and Tux

![MacOS](https://upload.wikimedia.org/wikipedia/commons/thumb/f/f2/Diagram_of_Mac_OS_X_architecture.svg/1280px-Diagram_of_Mac_OS_X_architecture.svg.png){: width="50%" height="50%"}



### Operating Systems Tasks
The operating system's primary task is to manage services as an interface between the user and the hardware. Examples include:
- File System: managing files on the user
- Device I/O: managing input from devices
- Processes: Starting, running, and stopping programs, and allowing multiple programs to run at once, i.e., program multiprogramming.
- Memory Management: Allocating runtime memory for process and separating memory between process and between user-space and the kernel-space.

![OS]({{site.baseurl}}/fig/OS.png)



### I/O
In computing, input/output or I/O (or, informally, io or IO) is the communication between an information processing system, such as a computer, and the outside world, possibly a human or another information processing system. Inputs are the signals or data received by the system and outputs are the signals or data sent from it. The term can also be used as part of an action; to "perform I/O" is to perform an input or output operation.

### As a user of the OS, you will see these interactions from two perspectives:
- Shell: You will use the shell to interact with the OS
- System Call API: You will program in C to interact with the OS
The big part of this interaction comes from the System Call API, which you will use the C programming language. Why C?
  - C is a low level language
  - The OS is written in C
  - Understanding the OS and C together is a natural process and will make you a better programmer


-----------------------

## The shell
This serves as an interface between the user and the kernel. When a user logs in, the login program checks the username and password, then start another program called the shell. The shell is a command line interpreter. It interprets the commands the user run and arranges for them to be executed. The commands are the programs to be run. The shell itself has different types of shell with its own set of commands and functions.

### Shell Prompt?
The prompt “$” is called the command prompt. While the prompt is being displayed, you can run a command.

### Shell Types?
- $ character is the default prompt.
- Bourne Shell(sh)
- Korn shell(ksh)
- Bourne Again shell(bash)
- POSIX shell(sh)
- C Shell — % character is the default prompt.
- C shell(csh)
- TENEX/TOPS C shell(tcsh)

![shell types](https://d1jnx9ba8s6j9r.cloudfront.net/blog/wp-content/uploads/2019/05/Evolution-of-Linux-Shells-Types-of-Shells-in-Linux-Edureka.png)


**The original Unix shell was made in the mid 19’s by Stephen R. Bourne. Bourne shell was the first shell to show up in Unix world. Bourne shell is commonly installed as /bin/sh on most versions of Unix.**

### Shell structutre
![shell structure](https://d1jnx9ba8s6j9r.cloudfront.net/blog/wp-content/uploads/2019/05/Shell-Architecture-Types-of-Shells-in-Linux-Edureka-528x205.png)

### BASH
![BASH](https://miro.medium.com/proxy/0*L0nhgi_19dlQJtzb.png)
Bash stands for Bourne Again Shell and it is the default shell on many Linux distributions today. It is also a sh-compatible shell and offers practical improvements over sh for programming and interactive use which includes:

- Command line editing
- Job Control
- Unlimited size command history
- Shell Functions and Aliases
- Unlimited size Indexed arrays
- Integer arithmetic in any base from two to sixty-four


## Text Editor
![text editor]({{site.baseurl}}/fig/texteditor.png){: width="150%" height="150%"}

## Text Editing on Unix and Linux 
- nano
- emacs
- vim

### nano
nano was first created in 1999 with the name TIP (This isn't Pico), by Chris Allegretta. His motivation was to create a free software replacement for Pico, which was not distributed under a free software license. The name was changed to nano on January 10, 2000 to avoid a naming conflict with the existing Unix utility tip. The name comes from the system of SI prefixes, in which nano is 1000 times larger than pico, In February 2001, nano became a part of the GNU Project. nano implements some features that Pico lacks, including colored text, regular expression search and replace, smooth scrolling, multiple buffers, rebindable key support, and (experimental) undoing and redoing of edit changes. On August 11, 2003, Chris Allegretta officially handed the source code maintenance for nano to David Lawrence Ramsey. On December 20, 2007, Ramsey stepped down as nano's maintainer. 

### emacs
GNU Emacs is an extensible, customizable text editor—and more. At its core is an interpreter for Emacs Lisp, a dialect of the Lisp programming language with extensions to support text editing. The features of GNU Emacs include:
- Content-sensitive editing modes, including syntax coloring, for a variety of file types including plain text, source code, and HTML.
- Complete built-in documentation, including a tutorial for new users.
- Full Unicode support for nearly all human languages and their scripts.
- Highly customizable, using Emacs Lisp code or a graphical interface. 
- A large number of extensions that add other functionality, including a project planner, mail and news reader, debugger interface, calendar, and more. Many of these extensions are distributed with GNU Emacs; others are available separately. 

### VIM

vi is a screen-oriented text editor originally created for the Unix operating system. The portable subset of the behavior of vi and programs based on it, and the ex editor language supported within these programs, is described by (and thus standardized by) the Single Unix Specification and POSIX. 

Over the years since its creation, vi became the de facto standard Unix editor and a nearly undisputed number one editor until the rise of Emacs after about 1984. The Single UNIX Specification specifies vi, so every conforming system must have it. 

A 2009 survey of Linux Journal readers found that vi was the most widely used text editor among respondents, beating gedit, the second most widely used editor by nearly a factor of two (36% to 19%). 


### Text editor Cheat sheet

- [VIM](https://preview.redd.it/ve1jv3m3qqj21.png?width=960&crop=smart&auto=webp&s=deb6dc83a462dc54523d703574e953638598af19)
- [nano](https://www.cheatography.com/bipinthite/cheat-sheets/nano-editor/)
- [Emacs](https://sachachua.com/blog/wp-content/uploads/2013/05/How-to-Learn-Emacs-v2-Large.png)


```
#!/bin/bash
echo "hello world"
```

### Shebang line

In order to make it possible to execute scripts as though they were first class executables, UNIX systems will looks for what we refer to as a shebang line at the top of the file. The origin of the name is murky. Some think it came from sharp-bang or hash-bang – contractions of # (“sharp”) and ! (“bang”). Others think the “SH” is in reference to the first UNIX shell, named “sh”.

In any case, if a UNIX system sees that the first line of an executable begins with #!, then it will execute the file using whatever command is specified in the rest of the line. For example, if there’s a file named /path/to/bar that looks like:
```
#!/bin/bash
```

```
#!/bin/python
```

```
#!/bin/perl
```

#### Advanced shebang line

For perl
```
#!/usr/bin/env perl
```
For python
```
#!/usr/bin/env python
```


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

### Let's install test package!
On Ubuntu systems:
```bash
$ sudo apt install screenfetch
```

On macOS systems: 
```bash
$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
$ brew update

$ brew install screenfetch
```

```bash
screenfetch
```

### Let's install prerequisite packages
On Ubuntu systems:
```bash
$ sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl vim
```

On macOS systems:
```bash
## Install Homebrew


$ brew install openssl readline sqlite3 xz zlib
```

## Explainshell
[Explainshell](https://explainshell.com/)



## Your first Unix command.
```bash
echo "Hello World!"
```

## Connect to Pronghorn
```bash
ssh <YOURID>@pronghorn.rc.unr.edu
```

## Basic commands

|Category|comnmand|
|---|---|
|Navigation| cd, ls, pwd|
|File creation|touch,nano,mkdir,cp,mv,rm,rmdir|
|Reading|more,less,head,tail,cat|
|Compression|zip,gzip,bzip2,tar,compress|
|Uncompression|unzip,gunzip,bunzip2,uncompress|
|Permissions|chmod|
|Help|man|

>## pwd
>Returns you the `p`resent `w`orking `d`irectory (print working directory)
>```bash
>$ pwd
>```
>```output
>/home/{username}
>```
> This means, you are now working in the home directory, which is located under your ID. The directory that will be the first directory when you log-in. You can also avoid writing the full path by using `~` in front of your username or simply. `~`
{: .keypoints}

>## mkdir
>To create a directory, mkdir `make directory` can be used.
>$ mkdir <DIRECTORY>
>```bash
>$ mkdir bch709_test
>```
> Unlike PC/Mac folders, here you can have space in different way (some special characters are okay). You can also specify the path where you want to create your new folder. 
>```bash
>$ mkdir bch709\ test
>```
> Like most Unix commands, mkdir supports command-line options which let you alter its behavior and functionality. Command-like options are – as the name suggests – optional arguments that are placed after the command name. They often take the form of single letters (following a dash). If we had used the -p option of the mkdir command we could have done this in one step.
>```bash
>$ mkdir -p bch709_test/help
>```
{: .keypoints}

>## cd
>This is `c`hanges `d`irectory. 
>$ cd <DIRECTORY>
>```bash
>$ cd bch709_test
>```
>```bash
>$ pwd
>```
>```output
>/home/{username}/bch709_test
>```
>Changes your present location to the parent directory
>```bash
>$ cd ..
>```
>Please check your current location with `pwd`.
>
>Changes your present location to the Grand parent (two levels upward) directory.
>```bash
>$ cd ../../
>```
>Please check your current location with `pwd`.
>**If this is confusing, please use our [command line server](http://35.199.189.11:8081/) to understand the location and folders.**
>Go to the root directory
>```bash
>$ cd /
>```
>Go to your home directory
>```bash
>$ cd ~
>```
>Go to specific directory
>```bash
>$ cd ~/bch709_test/help 
>```
>Go to upward directory
>$ cd ../
>
{: .keypoints}

>## Absolute and relative paths
>`cd` .. allows us to change directory relative to where we are now. You can also always change to a directory based on its absolute location. E.g. if you are working in the `~/bch709_test` directory and you want to change to the `~/bch709_test/help` directory, then you could do either of the following:
> **absolute path**
>```bash
>$ cd /home/<username>/bch709_test/help
>```
> **relative path**
>```bash
>$ cd ../
>```
{: .checklist}

>## The way back home
>Again, your home directory is `/home/<username>/`. . they take you back to your home directory (from wherever you were). You will frequently want to jump straight back to your home directory, and typing cd is a very quick way to get there.
>```bash
> cd ~
>```
>```bash
>cd
>```
>
>```bash
>cd ~/
>```
{: .checklist} 


>## Pressing < TAB > key
>You can type in first few letters of the directory name and then press tab to auto complete rest of the name (especially useful when the file/directory name is long).
>If there is more than one matching files/directories, pressing tab twice will list all the matching names.
>Saving keystrokes may not seem important now, but the longer that you spend typing in a terminal window, the happier you will be if you can reduce the time you spend at the keyboard. Especially as prolonged typing is not good for your body. So the best Unix tip to learn early on is that you can [tab complete][] the names of files and programs on mostUnix systems. Type enough letters to uniquely identify the name of a file, directory, or program and press tab – Unix will do the rest. E.g. if you type ‘tou’ and then press tab, Unix should autocomplete the word to ‘touch’ (this is a command which we will learn more
about in a minute). In this case, tab completion will occur because there are no other Unix commands that start with ‘tou’. If pressing tab doesn’t do anything, then you have not havetyped enough unique characters. In this case pressing tab twice will show you all possible completions. This trick can save you a LOT of typing!
{: .checklist}

>## Pressing <Arrow up/down> key
>You can also recall your previous commands by pressing up/down arrow or browse all your previously used commands by typing history on your terminal.
{: .checklist} 

>## Copy & Paste
>If you drag it will copy the text in terminal. If you use right click, then it will paste.
>In macOS, it depends on setting. You can still use command + c and command + v.
{: .checklist} 


>## ls
>
>The contents of a directory can be viewed using ls (list) command. 
>
>![ls]({{site.baseurl}}/fig/ls.png)   
>
>Like any other command, you can use absolute path or abbreviated path. There are also various options available for ls command.
>Some very useful options include:
>Lists files in lengthy or detailed view
>```bash
>$ ls –l
>```
>Lists files, sorted based on creation time
>```bash
>$ ls –t
>```
>Lists files, sorted based on size
>```bash
>$ ls –S
>```
>Lists all the files, do not ignore entries starting with .
>```bash
>$ ls -a
>```
>For each file or directory we now see more information (including file ownership and modification times). The ‘d’ at the start of each line indicates that these are directories. There
>are many, many different options for the ls command. Try out the following (against any
>directory of your choice) to see how the output changes
>
>```bash
>$ ls -l
>$ ls -R
>$ ls -l -t -r
>$ ls -ltr
>$ ls -lh
>```
>*Because your directory is almost empty. They will not show anything a lot.*
>Let's check other folder. Somthing like below.
>```bash
>$ ls -l /usr/bin/
>```
>![ls]({{site.baseurl}}/fig/ls2.png)  
{: .keypoints}


>## man & help
>If every Unix command has so many options, you might be wondering how you find out
>what they are and what they do. Well, thankfully every Unix command has an associated
>‘manual’ that you can access by using the man command. E.g.
>
>```bash
>$ man ls
>$ man cd
>$ man man
>```
>or
>
>you can call help directrly from command.
>```bash
>$ ls --help
>$ cd --help
>```
{: .keypoints}

>## rmdir
>To remove a directory, rmdir `remove directory` can be used.
>$ rmdir <DIRECTORY>
>```bash
>$ cd ~/
>$ ls
>```
>![ls3]({{site.baseurl}}/fig/ls3.png)
>```bash
>$ rmdir help
>```
>![ls4]({{site.baseurl}}/fig/ls4.png)
>```bash
>$ ls bch709_test
>$ rmdir bch709_test/help
>$ ls bch709_test
>```
>*You have to be outside a directory before you can remove it with rmdir*
{: .keypoints}

>## touch
>The following sections will deal with Unix commands that help us to work with files, i.e. copy files to/from places, move files, rename files, remove files, and most importantly, look at files. First, we need to have some files to play with. The Unix command touch will let us create a new, empty file.
>```bash
>$ touch --help
>$ touch test.txt
>$ touch exam.txt
>$ touch ETA.txt
>$ ls
>```
>![ls5]({{site.baseurl}}/fig/ls5.png)
>```bash
>$ ls -a
>```
{: .keypoints}


>## mv 
>We want to move these files to a new directory (‘temp’). We will do this using the Unix `mv`(move) command. Remember to use tab completion:
>```bash
>$ cd ~/
>$ mkdir Hello
>$ mv test.txt Hello
>$ mv exam.txt Hello
>$ mv ETA.txt Hello
>$ ls Hello
>```
>For the mv command, we always have to specify a source file (or directory) that we want to
>move, and then specify a target location. If we had wanted to, we could have moved both
>files in one go by typing any of the following commands:
>```bash
>$ mv *.txt Hello
>$ mv *t Hello
>$ mv *e* Hello
>```
>
>The asterisk * acts as a wild-card character, essentially meaning ‘match anything’. The second example works because there are no other files or directories in the directory that end with the letter ‘t’ (if there were, then they would be moved too). So all the 'bch709 test', bch709_test will be move. Likewise, the third example works because only those two files contain the letters ‘e’ in their names and 'bch709 test', bch709_test and Hello folder has e so it will move everything except for ETA.txt and will give you error `cannot move 'Hello' to a subdirectory of itself`.  Using wild-card characters can save you a lot of typing but need to be careful The ? character is also a wild-card but with a slightly different meaning. See if you can work out what it does.
{: .keypoints}

>## Renaming files by mv
In the earlier example, the destination for the mv command was a directory name (Hello). So we moved a file from its source location to a target location, but note that the target could have also been a (different) file name, rather than a directory. E.g. let’s make a new file and move it whilst renaming it at the same time
>```bash
>$ cd ~/
>$ touch rags
>$ ls
>$ mv rags Hello/riches
>$ ls Hello/
>```
>![ls6]({{site.baseurl}}/fig/ls6.png)
>*Of course you can do without moving*
>
>In this example we create a new file (‘rags’) and move it to a new location and in the process change the name (to ‘riches’). So mv can rename a file as well as move it. The logical extension of this is using `mv` to rename a file without moving it. You may also have access to a tool called `rename`, type man rename for more information. 
>``` bash
>$ mv Hello/riches Hello/rags
>```
{: .keypoints}

>## Moving directory by mv
>It is important to understand that as long as you have specified a ‘source’ and a ‘target’ location when you are moving a file, then it doesn’t matter what your current directory is. You can move or copy things within the same directory or between different directories regardless of whether you are in any of those directories. Moving directories is just like moving files:
>![ls6]({{site.baseurl}}/fig/ls6.png)
>```bash
>$ mv Hello/bch709_test .
>```
>`.` means current directory
>
>```bash
>$ ls Hello
>$ ls .
>```
>![ls7]({{site.baseurl}}/fig/ls7.png)
{: .keypoints}

>## rm
>`R`e`m`oving file
>
>You’ve seen how to remove a directory with the rmdir command, but rmdir won’t remove directories if they contain any files. So how can we remove the files we have created (inside
temp)? In order to do this, we will have to use the rm (remove) command. Please read the next section VERY carefully. Misuse of the rm command can lead to needless death & destruction. Potentially, rm is a very >dangerous command; if you delete something with rm, you will not get it back! It is possible to delete everything in your home directory (all directories and subdirectories) with rm.
>
>That is why it is such a dangerous command.  Let me repeat that last part again. It is possible to delete EVERY file you have ever created with the rm command. 
>Luckily there is a way of making  rm a little bit safer. We can use it with the -i (interactive) command-line option which will ask for confirmation before deleting anything (remember to use tab-completion):
>```bash
>$ cd ~/
>$ cd Hello 
>$ ls
>$ rm -i  ETA.txt exam.txt rags
>```
>*will ask permission for each step:*
>![ls8]({{site.baseurl}}/fig/ls8.png)
>```bash
>ls
>```
{: .keypoints}

>## cp 
>Copying files with the cp (copy) command has a similar syntax as mv, but the file will remain at the source and be copied to the target location. Remember to always specify a source and a target location. Let’s create a new file and make a copy of it:
>```bash
>$ cd ~/
>$ touch file1
>$ cp file1 file2
>$ ls
>```
>
>What if we wanted to copy files from a different directory to our current directory? Let’s put a file in our home directory (specified by \~, remember) and copy it to the lecture directory.
>```bash
>$ cd ~/
>$ touch ~/file3
>$ cp ~/file3 ~/Hello
>```
> The cp command also allows us (with the use of a command-line option) to copy entire directories. Use man cp to see how the -R or -r options let you copy a directory recursively.
> Please check help `cp --help`
{: .keypoints}



>## using . < Dot> ?
>In Unix, the current directory can be represented by a . (dot) character. You will often use for copying files to the directory that you are in. Compare the following:
>```bash
>ls
>ls .
>ls ./
>```
>In this case, using the dot is somewhat pointless because ls will already list the contents of the current directory by default. Also note how the trailing slash is optional.
{: .checklist} 

>## Clean up and start new
>Before we start other session, let's clean up folders
>```bash
>$ cd ~/
>$ ls 
>```
>Please do it first and check the solution below.
{: .challenge}

>## Clean the folder with contents
>How can we clean up `bch709_test`  `Hello` ?
>```bash
>rm -R  <folder name>
>```
{: .solution}


>## Downloading file
>Let's downloading one file from website. There are several command to download file. Such as 'wget','curl','rsync' etcs. but this time we will use `curl`.
>File location
>```
>https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/bch709_student.txt
>```
>How to use curl?
>![curl]({{site.baseurl}}/fig/curl.png)
>
>How to download from web file?
>
>```
>curl -L -o <output_name> <link>
>```
>```bash
>curl -L -o bch709_student.txt https://raw.githubusercontent.com/plantgenomicslab/BCH709/6055e8c5faf400119dc78b4d36f6b93814d9f76d/bch709_student.txt
> ls bch709_student.txt
>```
{: .keypoints}

>## Viewing file contents
>There are various commands to print the contents of the file in bash. Most of these commands are often used in specific contexts. All these commands when executed with filenames displays the contents on the screen. Most common ones are less, more, cat, head and tail.
>`less` FILENAME try this: `less bch709_student.txt` Displays file contents on the screen with line scrolling (to scroll you can use arrow keys, PgUp/PgDn keys, space bar or Enter key). When you are done press q to exit.  
>`more` FILENAME try this: `more bch709_student.txt` Like less command, also, displays file contents on the screen with line scrolling but uses only space bar or Enter key to scroll. When you are done press q to exit.  
>`cat` FILENAME try this: `cat bch709_student.txt` Simplest form of displaying contents. It catalogs the entire contents of the file on the screen. In case of large files, entire file will scroll on the screen without pausing  
>`head` FILENAME try this: `head bch709_student.txt` Displays only the starting lines of a file. The default is first ten lines. But, any number of lines can be displayed using –n option (followed by required number of lines).  
>`tail` FILENAME try this: `tail bch709_student.txtSimilar` to head, but displays the last 10 lines. Again –n option can be used to change this. More information about any of these commands can be found in man pages (man command)  
{: .keypoints}

>## Let's Download bigger data
>Please go to below wesite.
>```
>ftp://ftp.arabidopsis.org/home/tair/Sequences/ATH_cDNA_EST_sequences_FASTA/
>```
>Please download yellow marked file.
>![download]({{site.baseurl}}/fig/download.png)
{: .keypoints}

>## How to download
>![download2]({{site.baseurl}}/fig/download2.png)
>```bash
>curl -L -O ftp://ftp.arabidopsis.org/home/tair/Sequences/ATH_cDNA_EST_sequences_FASTA/ATH_cDNA_sequences_20101108.fas
>```
{: .solution}

>## Viewing file contents
>What is the difference in `less`, `tail`, `more` and `head`?
{: .discussion}


>## What are flags (parameters / options)?
>A “flag” in Unix terminology is a parameter added to the command. See for example
>```
>$ ls
>```
>versus
>```
>$ ls -l
>```
>Typically flags make programs behave differently or report their information in different ways.
{: .checklist} 

>## How to find out what flags are available?
>You can use the manual to learn more about a command:
>```bash
>$ man ls
>```
>will produce
>![man]({{site.baseurl}}/fig/man.png)
{: .checklist} 

>## What if the tools does not have manual page?
>Only tools that come with Unix will have a manual. For other software the -h, -help or
>--help command will typically print instructions for their use.
>```bash
>$ curl --help
>```
>If the help flag is not recognized you will need to find the manual for the software.
{: .checklist} 

>## What are flag formats?
>Traditionally Unix tools use two flag forms:
>- short form: single minus - then one letter, like -o, ’-L‘
>
>- long form: double minus -- then a word, like --output, --Location
>In each case the parameter may stand alone as a toggle (on/off) or may take additional values after the flag. -o <filename> or -L
>Now some bioinformatics tools do not follow this tradition and use a single - character for both short and long options. -g and -genome.
>**Using flags is a essential to using Unix. Bioinformatics tools typically take a large number of parameters specified *via* flags**. The correctness of the results critically depends on the proper use of these parameters and flags.
{: .checklist} 


>## How many lines does the file have?
>You can pipe the output of you stream into another program rather than the screen. Use the | (Vertical Bar) character to connect the programs. For example, the wc program is the word counter.
>```bash
>cat bch709_student.txt | wc
```
>prints the number of lines, words, and characters in the stream:
>```
>     41     141     900
>```
>How many lines?
>```bash
>cat bch709_student.txt | wc -l
>```
>```output
>41
>```
> *of course we can use `wc` directly*
>```bash
>wc -l bch709_student.txt 
>```
> That is equivalent to:
>```bash
>cat bch709_student.txt  | wc -l
>```
>In general, it is a better option to open a stream with cat then pipe the flow into the next program. Later you will see that it is easier to design, build and understand more complex pipelines when the data stream is opened at the beginning as a separate step.
>
>Again, let's do head.
>```bash
>cat bch709_student.txt  | head
>```
>![head]({{site.baseurl}}/fig/head.png)
>>## is this equivalent to?
>>```bash
>>head bch709_student.txt
>>```
>{: .solution}
{: checklist}

>## grep 
>Globally search a Regular Expression and Print is one of the most useful commands in UNIX and it is commonly used to filter a file/input, line by line, against a pattern eg., to print each line of a file which contains a match for pattern.
>Please check option with :
>```bash
>grep --help
>```
>With options, syntax is
>```
>grep [OPTIONS] PATTERN FILENAME
>```
>Let's find how many people use macOS?
>First we need to check file.
>```bash
>$ wc -l bch709_student.txt
>```
>```bash
>$ less bch709_student.txt
>```
>How can we find the macOS people?
>Now we can use `grep`
>```bash
>$ cat bch709_student.txt | grep MacOS
>```
>How can we count?
>```bash
>$ cat bch709_student.txt | grep MacOS | wc -l
>```
>How about this?
>```bash
>$ grep MacOS bch709_student.txt | wc -l
>```
>How about this?
>```bash
>$ grep -c MacOS bch709_student.txt
>```
>With other options (flags)
>```bash
>$ grep -v Windows  bch709_student.txt
>```
>With multiple options (flags)
>```bash
>$ grep -c -v Windows  bch709_student.txt
>```
>
>```bash
>$ grep --color  -i macos  bch709_student.txt
>```
{: checklist}

>## How do I store the results in a new file?  
>The > character is the redirection.
>```bash 
>$ grep  -i macos  bch709_student.txt > mac_student
>```
>```bash
>$ cat bch709_student.txt | grep -i windows > windows_student
>```
>Please check with `cat` or `less`
{: checklist}

>## Do you want to check difference?
>```bash
>$ diff mac_student windows_student
>```
>```bash
>$ diff -y  mac_student windows_student
>```
{: .keypoints}

>## How can I select name only? (cut)
>```bash
>$ cat bch709_student.txt  | cut -f 1
>```
>```bash
>$ cut -f 1 bch709_student.txt
>```
{: .keypoints}

>## How can I sort it ? (sort)
>```bash
>cut -f 1 bch709_student.txt | sort
>```
>```bash
>$ sort bch709_student.txt
>$ sort -k 2 bch709_student.txt
>$ sort -k 1 bch709_student.txt
>```
>```bash
>$ sort -k 1 bch709_student.txt | cut -f 1
>```
>Save?
>```bash
>sort -k 1 bch709_student.txt | cut -f 1 > name_sort
>```
{: .keypoints}


>## uniq 
>uniq command removes duplicate lines from a **sorted file**, retaining only one instance of the running matching lines. Optionally, it can show only lines that appear exactly once, or lines that appear more than once. uniq requires sorted input since it compares only consecutive lines.
>```bash
>$ cut -f 2 bch709_student.txt > os.txt
>```
>```bash
>$ uniq -c os.txt
>```
>![uniq]({{site.baseurl}}/fig/uniq.png)
>```bash
>$ sort os.txt | uniq -c
>```
>![uniq2]({{site.baseurl}}/fig/uniq2.png)
Of course you can use `sort` independently.
>```bash
>$ sort os.txt > os_sort.txt
>$ uniq -c os_sort.txt
>```
>```bash
>$ uniq --help
>```
>![uniq3]({{site.baseurl}}/fig/uniq3.png)
{: checklist}

>## diff 
>diff (difference) reports differences between files. A simple example for diff usage would be
>`
>$ diff FILEA FILEB
>`
>When you try to use new command, please check with `--help`.
>![diff]({{site.baseurl}}/fig/diff.png)
>```bash
>$ diff  os.txt  os_sort.txt
>```
>```bash
>$ diff -y os.txt  os_sort.txt
>```
>```bash
>$ sort os.txt | diff -y - os_sort.txt
>```
{: checklist}


### There are still a lot of command that you can use. Such as `paste`, `comm`, `join`, `split` etc.

### Oneliners
Oneliner, textual input to the command-line of an operating system shell that performs some function in just one line of input. This need to be done with "|".
[For advanced usage, please check this](https://plantgenomicslab.github.io/BCH709/onliner/index.html) 

### FASTA format
The original FASTA/Pearson format is described in the documentation for the FASTA suite of programs. It can be downloaded with any free distribution of FASTA (see fasta20.doc, fastaVN.doc or fastaVN.me—where VN is the Version Number).

The first line in a FASTA file started either with a ">" (greater-than; Right angle braket) symbol or, less frequently, a ";" (semicolon) was taken as a comment. Subsequent lines starting with a semicolon would be ignored by software. Since the only comment used was the first, it quickly became used to hold a summary description of the sequence, often starting with a unique library accession number, and with time it has become commonplace to always use ">" for the first line and to not use ";" comments (which would otherwise be ignored).

Following the initial line (used for a unique description of the sequence) is the actual sequence itself in standard one-letter character string. Anything other than a valid character would be ignored (including spaces, tabulators, asterisks, etc...). Originally it was also common to end the sequence with an "\*" (asterisk) character (in analogy with use in PIR formatted sequences) and, for the same reason, to leave a blank line between the description and the sequence.

![fasta]({{site.baseurl}}/fig/fasta.png)

#### Description line
The description line (defline) or header/identifier line, which begins with '>', gives a name and/or a unique identifier for the sequence, and may also contain additional information. In a deprecated practice, the header line sometimes contained more than one header, separated by a ^A (Control-A) character. In the original Pearson FASTA format, one or more comments, distinguished by a semi-colon at the beginning of the line, may occur after the header. Some databases and bioinformatics applications do not recognize these comments and follow the [NCBI FASTA specification](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp). 

### FASTA file handling with command line.
Please check one fasta file
```bash
$ ls ATH_cDNA_sequences_20101108.fas
```

>## previous example
>![download2]({{site.baseurl}}/fig/download2.png)
>```bash
>curl -L -O ftp://ftp.arabidopsis.org/home/tair/Sequences/ATH_cDNA_EST_sequences_FASTA/ATH_cDNA_sequences_20101108.fas
>```
{: .solution}

### count cDNA
How many cDNA in this fasta file?
Please use `grep`  `wc` to find number.
>## fasta count
>![fasta_count]({{site.baseurl}}/fig/fasta_count.png)
{: .solution}

>## count DNA letter
>How many sequences (DNA letter) in this fasta file?
>Please use `grep`  `wc` to find number.
{: .discussion}


### what is GI and GB?
[NCBI identifiers](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/)
![NCBI identifiers]({{site.baseurl}}/fig/NCBI_identifiers.png)

### Collect GI
How can I collect GI from FASTA description line?
Please use `grep`  `cut` to find number.


>## Sequence redundancy
>Does GI have any redundancy?
>Please use `grep`  `wc` `diff` to solve.
{: .discussion}

### GFF file
The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines. The following documentation is based on the Version 3 (http://gmod.org/wiki/GFF3) specifications.

Please download below file
`
http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz
`
What is `.gz` ?
`
file MGI.gff3.gz
`


### Compression
There are several options for archiving and compressing groups of files or directories. Compressed files are not only easier to handle (copy/move) but also occupy less size on the disk (less than 1/3 of the original size). In Linux systems you can use zip, tar or gz for archiving and compressing files/directories.

>## ZIP compression/extraction
>zip OUTFILE.zip INFILE.txt Compress INFILE.txt
>zip -r OUTDIR.zip DIRECTORY Compress all files in a DIRECTORY into one archive file (OUTDIR.zip)
>zip -r OUTFILE.zip . -i \*.txt Compress all txt files in a DIRECTORY into one archive file (OUTFILE.zip)
>unzip SOMEFILE.zip
{: checklist}

>## TAR compression/extraction
>tar (tape archive) utility saves many files together into a single archive file, and restores individual files from the archive. It also includes automatic archive compression/decompression options and special features for incremental and full backups.
tar -cvf OUTFILE.tar INFILE
archive INFILE
tar -czvf OUTFILE.tar.gz INFILE
archive and compress file INFILE
tar -tvf SOMEFILE.tar
list contents of archive SOMEFILE.tar
tar -xvf SOMEFILE.tar
>
>
>extract contents of SOMEFILE.tar
>tar -xzvf SOMEFILE.tar.gz
>extract contents of gzipped archive SOMEFILE.tar.gz
>tar -czvf OUTFILE.tar.gz DIRECTORY
>archive and compress all files in a directory into one archive file
>tar -czvf OUTFILE.tar.gz \*.txt
>archive and compress all ".txt" files in current directory into one archive file
>tar -czvf backup.tar.gz BACKUP_WORKSHOP
{: checklist}

## Gzip compression/extraction
>gzip (gnu zip) compression utility designed as a replacement for compress, with much better compression >and no patented algorithms. The standard compression system for all GNU software.
>gzip SOMEFILE compress SOMEFILE (also removes uncompressed file)
>gunzip SOMEFILE.gz uncompress SOMEFILE.gz (also removes compressed file)
>
>gzip the file MGI.gff3.gz and examine the size. gunzip it back so that you can use this file for thelater exercises.
>```bash
>$ gunzip MGI.gff3.gz
>$ ls –lh
>$ gzip MGI.gff3
>$ ls -lh
>$ gunzip MGI.gff3.gz
>```
{: checklist}


## GFF3 Annotations


Print all sequences annotated in a GFF3 file.

    cut -s -f 1,9 MGI.gff3 | grep $'\t' | cut -f 1 | sort | uniq


Determine all feature types annotated in a GFF3 file.

    grep -v '^#' MGI.gff3 | cut -s -f 3 | sort | uniq


Determine the number of genes annotated in a GFF3 file.

    grep -c $'\tgene\t' MGI.gff3


Extract all gene IDs from a GFF3 file.

    grep $'\tgene\t' MGI.gff3 | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'

Print all CDS.
    
    cat MGI.gff3 | cut -f 3 | grep CDS | 

Print CDS and ID
```bash
    cat MGI.gff3 | cut -f 1,3,4,5,7,9 | head  
    cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | head  
    cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | sed 's/;.*//g' | head  
    cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | sed 's/;.*//g' | sed 's/ID=//g' | head  
    cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep $'\tCDS\t' | sed 's/;.*//g' | sed 's/ID=//g' | head  
```
Print length of each gene in a GFF3 file.

    grep $'\tgene\t' MGI.gff3 | cut -s -f 4,5 | perl -ne '@v = split(/\t/); printf("%d\n", $v[1] - $v[0] + 1)'

Extract all gene IDs from a GFF3 file.

    grep $'\tgene\t' MGI.gff3 | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'

Time and again we are surprised by just how many applications it has, and how frequently
problems can be solved by sorting, collapsing identical values, then resorting by the collapsed
counts.
The skill of using Unix is not just that of understanding the commands themselves. It is
more about recognizing when a pattern, such as the one that we show above, is the solution
to the problem that you wish to solve. The easiest way to learn to apply these patterns is
by looking at how others solve problems, then adapting it to your needs.

<!--
FASTA header lines to GFF format (assuming the length is in the header as an appended "\_length" as in [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/) assembled transcripts):

    grep '>' file.fasta | awk -F "_" 'BEGIN{i=1; print "##gff-version 3"}{ print $0"\t BLAT\tEXON\t1\t"$10"\t95\t+\t.\tgene_id="$0";transcript_id=Transcript_"i;i++ }' > file.gff


###

16 centromere
16 centromere_DNA_Element_I
16 centromere_DNA_Element_II
16 centromere_DNA_Element_III
8 external_transcribed_spacer_region
24 five_prime_UTR_intron
To find out the most frequent unique elements, we need to sort the output of uniq.
cat types.txt | sort | uniq -c | sort -rn | head
it now prints:
7074 CDS
6604 ORF
484 noncoding_exon
383 long_terminal_repeat
377 intron
352 ARS
299 tRNA_gene
196 ARS_consensus_sequence
91 transposable_element_gene
77 snoRNA_gene




cat MGI.gff3 | cut -f 2 | head
Build your commands one step at a time, always checking that you are on the right track:
cat MGI.gff3 | head
cat MGI.gff3 | cut -f 2 | head
cat MGI.gff3 | cut -f 2 | grep ORF | head
12.3.16 How many genes are there?
cat MGI.gff3 | cut -f 2 | grep ORF | wc -l
12.3.17 Can I select multiple columns?
cat MGI.gff3 | cut -f 2,3,4 | grep ORF | head
What does this do?
cat MGI.gff3 | cut -f 2,3,4 | grep ORF | grep -v Dubious | wc -l
12.4 How many feature types are in this data?
We are going to use this data a lot, so place it into a separate file for now.
cat MGI.gff3 | cut -f 2 > types.txt
Sorting places identical consecutive entries next to one another.
cat types.txt | sort | head
Find unique words. The uniq command collapses consecutive identical words into one.
cat types.txt | sort | uniq | head
Using -c flag to uniq will not only collapse consecutive entries it will print their counts.
cat types.txt | sort | uniq -c | head
prints:
352 ARS
196 ARS_consensus_sequence
6 blocked_reading


12.5 The single most useful Unix pattern
The pattern sort | uniq -c | sort -rn is perhaps the most useful simple unix command
Time and again we are surprised by just how many applications it has, and how frequently
problems can be solved by sorting, collapsing identical values, then resorting by the collapsed
counts.
The skill of using Unix is not just that of understanding the commands themselves. It is
more about recognizing when a pattern, such as the one that we show above, is the solution
to the problem that you wish to solve. The easiest way to learn to apply these patterns is
by looking at how others solve problems, then adapting it to your needs.pattern that you will ever learn.
Time and again we are surprised by just how many applications it has, and how frequently
problems can be solved by sorting, collapsing identical values, then resorting by the collapsed
counts.
The skill of using Unix is not just that of understanding the commands themselves. It is
more about recognizing when a pattern, such as the one that we show above, is the solution
to the problem that you wish to solve. The easiest way to learn to apply these patterns is
by looking at how others solve problems, then adapting it to your needs.


-->


This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.

, a tutorial that teaches you how to work at the command-line. You'll learn all the basic skills needed to start being productive in the UNIX terminal.


This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.
