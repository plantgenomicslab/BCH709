---
layout: page
title: Linux Enviroment and command line
published: true
---

{% include gh_variables.html %}


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

### The kernel
This is called as the hub of the operating system, serving as allocator of time and memory to programs and handling the filestore and communications in response to system calls.

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
![text editor](https://w.namu.la/s/faf332270510cbe6d8359870f39e006df2f67b7da39fa6bb994a8b43f86ca2c82c35c165e084aa611c18ed754995f1ac3c2db553e2a64cff5945b582bba3ed13dc812a4340d702808d45887eb414f905cab17ba7a7ec73026b41bb16dfe764f0){: width="300%" height="300%"}

### Text Editing on Unix and Linux 
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
$ apt install screenfetch
```

On macOS systems: 
```bash
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
$ brew install openssl readline sqlite3 xz zlib
```

This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.


## Command Line Bootcamp
[Command-line bootcamp](http://35.199.189.11:8081/), a tutorial that teaches you how to work at the command-line. You'll learn all the basic skills needed to start being productive in the UNIX terminal.


This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.
