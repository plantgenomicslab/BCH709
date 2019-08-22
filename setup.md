---
layout: page
title: Setup
published: true
---

# Overview

This lecture is designed to be run on terminal program. With the exception of a GUI program, all of the software and data used in terminal. Please follow the instructions below to prepare your computer for the lecture.

## Required additional software

This lesson requires a working web browser, terminal, spreadsheet program. If you don't have a spreadsheet program already, you can use MSOffice for free with your affilated email. 
Following guide will provide the specific information about required software. Please select your [OS system](https://en.wikipedia.org/wiki/Operating_system) to install software.

> ## Windows 10
> - Install MSOffice by going to [the installation page](hhttps://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version for Windows should automatically be selected. Once the installer is downloaded, double click on it and MSOffice should install.
>
>
> - Install [Bash on Windows](https://www.linux.com/news/bash-windows-what-does-it-mean/) by Enable “Windows Subsystem for Linux” feature.
> - Open powershell with right click and "choose Run it as administrator"
>![powershell](https://i1.wp.com/itsfoss.com/wp-content/uploads/2016/08/Powershell-Ubuntu-install.jpeg?w=800&ssl=1){: width="50%" height="50%"}
>Once you have the PowerShell running, use the command below to enable Bash in Windows 10.
>~~~
>Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
>~~~
>- You’ll be asked to confirm your choice. Type Y or press enter:
>![powershell](https://i1.wp.com/itsfoss.com/wp-content/uploads/2016/08/Powershell-Ubuntu-install-2.jpeg?w=799&ssl=1){: width="70%" height="70%"}
> - Now you should be asked to reboot. Even if you are not asked to, you must restart your system. Once your system has rebooted, go to the Windows Store and search for “Linux.”
>![Windows Store](https://i2.wp.com/itsfoss.com/wp-content/uploads/2016/08/install-ubuntu-windows-10-linux-subsystem-3-1.jpeg?w=800&ssl=1){: width="70%" height="70%"}
> - **Please install Ubuntu 18.04 LTS.**
>![Ubuntu](https://i0.wp.com/itsfoss.com/wp-content/uploads/2016/08/install-ubuntu-windows-10-linux-subsystem-7.jpeg?w=800&ssl=1){: width="70%" height="70%"}
> - Once you choose the distribution of your choice, you’ll see the option to install it. Do note that it will download files of around 1Gb in size. So you should have a good internet connection here.
> - Press "windows key + s" and search "ubuntu". *If you have Contana, you can use it at home*.
> - Please be patience then it will ask id and password. 
>![Linux](https://i2.wp.com/itsfoss.com/wp-content/uploads/2016/08/install-ubuntu-windows-10-linux-subsystem-10.jpeg?w=800&ssl=1){: width="70%" height="70%"}
>- If you find following message, please go back to top of this article and do it again.
>~~~
>The WSL optional component is not enabled. Please enable it and try again.
>See https://aka.ms/wslinstall for details.
>Error: 0x8007007e
>Press any key to continue...
>~~~
>- If you cannot get the Fall Creator’s update on Windows 10 for some reason, please update your Windows 10 first. If you still have a trouble, please follow "Windows 7 and under section"
{: .solution}


> ## Windows 7 and under
> - Install MSOffice by going to [the installation page](hhttps://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version for Windows should automatically be selected. Once the installer is downloaded, double click on it and MSOffice should install.
> - Install Putty by going to [the installation page](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html). For most newer computers, click on putty-64bit-X.XX-installer.msi to download the 64-bit version. If you have an older laptop, you may need to get the 32-bit version putty-X.XX-installer.msi. If you aren't sure whether you need the 64 or 32 bit version, you can check your laptop version by following [the instructions here](https://support.microsoft.com/en-us/help/15056/windows-32-64-bit-faq)
> - Once the installer is downloaded, double click on it, and PuTTY should install.
{: .solution}


> ## Mac OS X
> - Install MSOffice by going to [the installation page](hhttps://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version cannot be selected. Once the installer is downloaded, double click on it and MSOffice should install.
> - Once the installer is downloaded, double click on it and MSOffice should install.
>
> - Mac has native Terminal (command + space) and search terminal
> - However I personally recommend to use [iTerm2](https://www.iterm2.com/) for better experience.
{: .solution}

> ## Linux
>  - Install LibreOffice by going to [the installation page](https://www.libreoffice.org/download/libreoffice-fresh/). The version for Linux should automatically be selected. Click Download Version X.X.X (whichever is the most recent version). You will go to a page that asks about a donation, but you don't need to make one. Your download should begin automatically.  
> - Once the installer is downloaded, double click on it and LibreOffice should install.
>
>
> - Please find terminal  Ctrl+Alt+T or terminal or xterm or uxterm
{: .solution}


## High Performance Computing and Cloud service

### Option A: Using the lessons on your local machine

It is possible to work through the lessons on your local machine (i.e. without using cloud or Pronghorn). However I would recomment to experience in High Performance Computing and Cloud service.

### Option B: Using Pronghorn (High-Performance Computing)

**Pronghorn** is the University of Nevada, Reno's new High-Performance Computing (HPC) cluster. The GPU-accelerated system is designed, built and maintained by the Office of Information Technology's HPC Team. Pronghorn and the HPC Team supports general research across the Nevada System of Higher Education (NSHE).

Pronghorn is composed of CPU, GPU, and Storage subsystems interconnected by a 100Gb/s non-blocking Intel Omni-Path fabric. The CPU partition features 93 nodes, 2,976 CPU cores, and 21TiB of memory. The GPU partition features 44 NVIDIA Tesla P100 GPUs, 352 CPU cores, and 2.75TiB of memory. The storage system uses the IBM SpectrumScale file system to provide 1PB of high-performance storage. The computational and storage capabilities of Pronghorn will regularly expand to meet NSHE computing demands.

Pronghorn is collocated at the Switch Citadel Campus located 25 miles East of the University of Nevada, Reno. Switch is the definitive leader of sustainable data center design and operation. The Switch Citadel is rated Tier 5 Platinum, and will be the largest, most advanced data center campus on the planet.

Pronghorn is available to all University of Nevada, Reno faculty, staff, students, and sponsored affiliates. Priority access to the system is available for purchase.
[Please apply your account here](https://www.unr.edu/research-computing/hpc-accounts)

### Option C: Using Google Cloud Platform
![Run on Google Cloud](https://lh5.googleusercontent.com/Isv8FIUtRG4H-gJL5kVtzptP7Db3UcjuAOohBWD_6VVwmDsp6B_g0ZsjBWFpuVstt-AT_dtu9O4FEykGfQ52vA5-iaZ5x4lA4KEEvUQ6iONTwuJ_3Q=w472)
Currently this lecture is supported by Google Cloud Platform. If you want to try, we will provide free credit for your usage.
[![Run on Google Cloud](https://storage.googleapis.com/cloudrun/button.svg)](https://cloud.google.com/)



If you are signed up to take a Genomics Data Carpentry workshop, you do *not* need to worry about setting up an AMI instance. The Carpentries
staff will create an instance for you and this will be provided to you at no cost. This is true for both self-organized and centrally-organized workshops. Your Instructor will provide instructions for connecting to the AMI instance at the workshop.

If you would like to work through these lessons independently, outside of a workshop, you will need to start your own AMI instance. 
Follow these [instructions on creating an Amazon instance](https://datacarpentry.org/genomics-workshop/AMI-setup/). Use the AMI `ami-0985860a69ae4cb3d` (Data Carpentry Genomics Beta 2.0 (April 2019)) listed on the Community AMIs page. Please note that you must set your location as `N. Virginia` in order to access this community AMI. You can change your location in the upper right corner of the main AWS menu bar. The cost of using this AMI for a few days, with the t2.medium instance type is very low (about USD $1.50 per user, per day). Data Carpentry has *no* control over AWS pricing structure and provides this
cost estimate with no guarantees. Please read AWS documentation on pricing for up-to-date information.

If you're an Instructor or Maintainer or want to contribute to these lessons, please get in touch with us [team@carpentries.org](mailto:team@carpentries.org) and we will start instances for you. 


### Data
The data used in this lecture will be available by direct link.

[FigShare Data Carpentry Genomics Beta 2.0](https://figshare.com/articles/Data_Carpentry_Genomics_beta_2_0/7726454)

More information about these data will be presented in the [first lesson of the workshop](http://www.datacarpentry.org/organization-genomics/data/).
