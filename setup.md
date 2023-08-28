---
layout: page
title: Setup
published: true
---

# Overview

This lecture is designed to be run on a terminal program. With the exception of a GUI program, all of the software and data used in the terminal. Please follow the instructions below to prepare your computer for the lecture.

**Please do not buy a new laptop. Your level is not that high, and our lecture will not need a new fancy computer.**
**Please read below and think about the best way.**


## Required Operating System + Laptop

> ## Windows 10 or 11
> - Install MSOffice by going to [the installation page](https://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/) or [connect this site](https://oit.unr.edu/services-and-support/data-storage/office-365/). The version for Windows should automatically be selected. Once the installer is downloaded, double click on it and MSOffice should install.
> ### Windows Subsystem for Linux (WSL2)
> Developers can access the power of both Windows and Linux at the same time on a Windows machine. The Windows Subsystem for Linux (WSL) lets developers install a Linux distribution (such as Ubuntu, OpenSUSE, Kali, Debian, Arch Linux, etc) and use Linux applications, utilities, and Bash command-line tools directly on Windows, unmodified, without the overhead of a traditional virtual machine or dualboot setup.
> ### Prerequisites
> -------------
> You must be running Windows 10 version 2004 and higher (Build 19041 and higher) or Windows 11 to use the commands below. If you are on earlier versions please see [the manual install page](https://learn.microsoft.com/en-us/windows/wsl/install-manual).
> ### Install WSL command
> -------------------
> You can now install everything you need to run WSL with a single command. Open PowerShell or Windows Command Prompt in **administrator** mode by right-clicking and selecting "Run as administrator", enter the wsl --install command, then restart your machine.
>    ```wsl --install --distribution ubuntu```  
> This command will enable the features necessary to run WSL and install the Ubuntu distribution of Linux. ([This default distribution can be changed](https://learn.microsoft.com/en-us/windows/wsl/basic-commands#install)).
> The first time you launch a newly installed Linux distribution, a console window will open and you'll be asked to wait for files to de-compress and be stored on your machine. All future launches should take less than a second.
> 
> ### Set up your Linux user info
> ---------------------------
> Once you have installed WSL, you will need to create a user account (Please use your NETID @ UNR) and password for your newly installed Linux distribution. See the [Best practices for setting up a WSL development environment](https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password) guide to learn more.
> ### Set up and best practices
> -------------------------
> We recommend following our [Best practices for setting up a WSL development environment](https://learn.microsoft.com/en-us/windows/wsl/setup/environment) guide for a step-by-step walk-through of how to set up a user name and password for your installed Linux distribution(s), using basic WSL commands, installing and customizing Windows Terminal, set up for Git version control, code editing and debugging using the VS Code remote server, good practices for file storage, setting up a database, mounting an external drive, setting up GPU acceleration, and more.
{: .solution}


> ## Windows 7 and under
> - Install MSOffice by going to [the installation page](hhttps://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version for Windows should automatically be selected. Once the installer is downloaded, double click on it and MSOffice should install.
> - Install Putty by going to [the installation page](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html). For most newer computers, click on putty-64bit-X.XX-installer.msi to download the 64-bit version. If you have an older laptop, you may need to get the 32-bit version putty-X.XX-installer.msi. If you aren't sure whether you need the 64 or 32-bit version, you can check your laptop version by following [the instructions here](https://support.microsoft.com/en-us/help/15056/windows-32-64-bit-faq)
> - Once the installer is downloaded, double click on it, and PuTTY should install.
{: .solution}


> ## Mac OS X
> - Install MSOffice by going to [the installation page](hhttps://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version cannot be selected. Once the installer is downloaded, double click on it and MSOffice should install.
> - Once the installer is downloaded, double click on it and MSOffice should install.
>
> - Mac has native Terminal (command + space) and search terminal
> - However, I personally recommend to use [iTerm2](https://www.iterm2.com/) for a better experience.
> Additionally, I would recommend to install [MacOS command-line tools with the following instruction](https://plantgenomicslab.github.io/BCH709/CLT/index.html).
{: .solution}

> ## Linux
>  - Install LibreOffice by going to [the installation page](https://www.libreoffice.org/download/libreoffice-fresh/). The version for Linux should automatically be selected. Click Download Version X.X.X (whichever is the most recent version). You will go to a page that asks about a donation, but you don't need to make one. Your download should begin automatically.  
> - Once the installer is downloaded, double click on it and LibreOffice should install.
>
>
> - Please find terminal  Ctrl+Alt+T or terminal or xterm or uxterm
{: .solution}

> ## Old laptop
>If you're laptop is old; We still have an option. Please choose one option in High-Performance Computing and Cloud service and connect to [Google Cloud SSh](https://ssh.cloud.google.com/)
{: .solution}

> ## Chromebook
>Of course you can do it with [Secure Shell App](https://chrome.google.com/webstore/detail/secure-shell-app/pnhechapfaindjhompbnflcldabbghjo)
{: .solution}

> ## Phone (Android/IOS)
>Of course you can do it with [Termius](https://www.termius.com/)
>Not recommended.
{: .solution}


> ## Pencil 
>Of course you can do it with [Termius](https://www.termius.com/)
>Not recommended.
{: .solution}

## Required Additional Software

This lesson requires a working web browser, terminal, spreadsheet program. If you don't have a spreadsheet program already, you can use MSOffice for free with your affiliated email. 
The following guide will provide the specific information about the required software. Please select your [OS system](https://en.wikipedia.org/wiki/Operating_system) to install the software.

## High-Performance Computing and Cloud service

### Option A: Using the lessons on your local machine

It is possible to work through the lessons on your local machine (i.e., without using cloud or Pronghorn). However, I would recommend experiencing in High-Performance Computing and Cloud service.

### Option B: Using Pronghorn (High-Performance Computing)

**Pronghorn** is the University of Nevada, Reno's new High-Performance Computing (HPC) cluster. The GPU-accelerated system is designed, built and maintained by the Office of Information Technology's HPC Team. Pronghorn and the HPC Team supports general research across the Nevada System of Higher Education (NSHE).

Pronghorn is composed of CPU, GPU, and Storage subsystems interconnected by a 100Gb/s non-blocking Intel Omni-Path fabric. The CPU partition features 93 nodes, 2,976 CPU cores, and 21TiB of memory. The GPU partition features 44 NVIDIA Tesla P100 GPUs, 352 CPU cores, and 2.75TiB of memory. The storage system uses the IBM SpectrumScale file system to provide 1PB of high-performance storage. The computational and storage capabilities of Pronghorn will regularly expand to meet NSHE computing demands.

Pronghorn is collocated at the Switch Citadel Campus located 25 miles East of the University of Nevada, Reno. Switch is the definitive leader of sustainable data center design and operation. The Switch Citadel is rated Tier 5 Platinum, and will be the largest, most advanced data center campus on the planet.

Pronghorn is available to all University of Nevada, Reno faculty, staff, students, and sponsored affiliates. Priority access to the system is available for purchase.
[Please apply your account here](https://www.unr.edu/research-computing/hpc-accounts)

![Pronghorn system map](./fig/pronghorn.png){: width="70%" height="70%"}


### Option C: Using Google Cloud Platform (GCP)
Currently, this lecture is supported by the Google Cloud Platform. If you want to try, we will provide free credit for your usage (about USD $1.50 per user, per day). If you choose this option, you need to have a Google account. We will mainly use the compute engine to run our bioinformatics project. Please check below.
[![Run on Google Cloud](https://storage.googleapis.com/cloudrun/button.svg)](https://cloud.google.com/compute/)
Won Yim has *no* control over GCP pricing structure and provides this
cost estimate with no guarantees. Please read the documentation on pricing for up-to-date information.

### Option D: Using Amazon Web Service
If you would like to work through Amazon Web Service (AWS). I will not stop you.
The cost of using this AWS for a few days, with the t2.medium instance type is very low (about USD $1.50 per user, per day). Won Yim has *no* control over AWS pricing structure and provides this
cost estimate with no guarantees. Please read GCP documentation on pricing for up-to-date information.

### Data
The data used in this lecture will be available by direct link.
