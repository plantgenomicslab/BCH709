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

> ## Windows 10 and Windows 11
> - Install MSOffice by going to [the installation page](https://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/) or [connect this site](https://oit.unr.edu/services-and-support/data-storage/office-365/). The version for Windows should automatically be selected. Once the installer is downloaded, double click on it and MSOffice should install.
>
> **Modern WSL Installation (Windows 10 version 2004+ or Windows 11)**
>
> The easiest way to install WSL is using a single command:
> - Open PowerShell as Administrator (right click and choose "Run as administrator")
>![powershell](https://i1.wp.com/itsfoss.com/wp-content/uploads/2016/08/Powershell-Ubuntu-install.jpeg?w=800&ssl=1){: width="50%" height="50%"}
> - Run the following command:
>~~~
>wsl --install
>~~~
> - This will install WSL and Ubuntu (latest LTS version) by default
> - Restart your computer when prompted
>
> **Install Ubuntu 24.04 LTS specifically:**
>
> If you want to install Ubuntu 24.04 LTS specifically, use:
>~~~
>wsl --install -d Ubuntu-24.04
>~~~
>
> You can also install Ubuntu 24.04 LTS from the [Microsoft Store](https://learn.microsoft.com/windows/wsl/install) or search for "Ubuntu 24.04" in the Microsoft Store app.
> ![Ubuntu](https://i0.wp.com/itsfoss.com/wp-content/uploads/2016/08/install-ubuntu-windows-10-linux-subsystem-7.jpeg?w=800&ssl=1){: width="70%" height="70%"}
>
> **After Installation:**
> - Once installed, launch Ubuntu from the Start menu
> - The first time you launch Ubuntu, it will ask you to create a username and password
> - Please be patient during the initial setup
>![Linux](https://i2.wp.com/itsfoss.com/wp-content/uploads/2016/08/install-ubuntu-windows-10-linux-subsystem-10.jpeg?w=800&ssl=1){: width="70%" height="70%"}
>
> **Troubleshooting:**
> - If you encounter issues, ensure your Windows is up to date
> - For older versions of Windows 10 (before version 2004), please update Windows first
> - Visit [Microsoft's official WSL documentation](https://learn.microsoft.com/windows/wsl/install) for detailed troubleshooting
>
> **Visual Studio Code (Required):**
> - **Download:** Visit [https://code.visualstudio.com/](https://code.visualstudio.com/) and download for Windows
> - Run the installer (.exe file) and follow the installation wizard
> - During installation, check these options:
>   - "Add to PATH" (recommended)
>   - "Create a desktop icon" (optional)
>   - "Add 'Open with Code' action to context menu" (recommended)
> - **WSL Integration:** After installing VS Code, install the "WSL" extension
>   - Open VS Code and press `Ctrl+Shift+X` to open Extensions
>   - Search for "WSL" and install the official extension by Microsoft
>   - This allows you to open WSL directories directly in VS Code
> - **Why VS Code for Bioinformatics?**
>   - Built-in terminal for running commands
>   - Syntax highlighting for Python, R, Bash, and bioinformatics formats
>   - Git integration for version control
>   - Remote development support (WSL, SSH, containers)
>   - Free and open source
>{: .solution}


> ## Windows 7 and under
> - Install MSOffice by going to [the installation page](https://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version for Windows should automatically be selected. Once the installer is downloaded, double click on it and MSOffice should install.
> - Install Putty by going to [the installation page](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html). For most newer computers, click on putty-64bit-X.XX-installer.msi to download the 64-bit version. If you have an older laptop, you may need to get the 32-bit version putty-X.XX-installer.msi. If you aren't sure whether you need the 64 or 32-bit version, you can check your laptop version by following [the instructions here](https://support.microsoft.com/en-us/help/15056/windows-32-64-bit-faq)
> - Once the installer is downloaded, double click on it, and PuTTY should install.
{: .solution}


> ## Mac OS X
> - Install MSOffice by going to [the installation page](https://oit.unr.edu/services-and-support/software-and-online-applications/software-purchasing-and-installation/microsoft-office-365-for-personal-computers/install-microsoft-office-for-home-student/). The version cannot be selected. Once the installer is downloaded, double click on it and MSOffice should install.
>
> **Terminal Application (Required):**
> - Mac has native Terminal (press command + space and search for "terminal")
> - **Recommended:** I strongly recommend using [iTerm2](https://www.iterm2.com/) for a better terminal experience
>   - iTerm2 offers split panes, search, autocomplete, and better customization
>   - Download from [https://iterm2.com/](https://www.iterm2.com/)
>   - Free and widely used by developers
>
> **Code Editor / IDE (Required):**
> - **Visual Studio Code (Required):**
>   - **Download:** Visit [https://code.visualstudio.com/](https://code.visualstudio.com/) and download for macOS
>   - Open the downloaded .zip file and drag "Visual Studio Code.app" to Applications folder
>   - Launch VS Code from Applications or via Spotlight (`Cmd+Space`, type "Visual Studio Code")
>   - **First-time setup:** VS Code may ask to install command line tools - click "Install" to enable `code` command in Terminal
>   - **Recommended Extensions for Bioinformatics:**
>     - Python (Microsoft)
>     - R (REditorSupport)
>     - Remote - SSH (for connecting to HPC clusters)
>     - Jupyter (for notebook support)
>   - **Why VS Code?**
>     - Free and open source
>     - Integrated terminal
>     - Excellent for Python, R, Bash scripting
>     - Git integration
>     - Remote development (SSH to Pronghorn HPC)
>     - Large extension ecosystem
>
> - **Alternative - Sublime Text (Optional):**
>   - If you prefer a lighter text editor instead of a full IDE
>   - Fast, lightweight, and powerful
>   - Download from [https://www.sublimetext.com/](https://www.sublimetext.com/)
>   - Free to evaluate (unlimited trial), license available for purchase
>   - Note: VS Code is required for the course, but you may use Sublime Text for simple text editing
>
> **Command Line Tools (Required):**
> - Install Xcode Command Line Tools. Open Terminal (or iTerm2) and run:
>~~~
>xcode-select --install
>~~~
> - For detailed instructions, see [Command Line Tools Installation Guide](https://plantgenomicslab.github.io/BCH709/CLT/index.html)
>
> **Homebrew Package Manager (Highly Recommended):**
> - Homebrew is the most popular package manager for macOS, essential for installing bioinformatics software
> - Install Homebrew by running this command in Terminal:
>~~~
>/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
>~~~
> - After installation, follow any instructions to add Homebrew to your PATH
> - Verify installation:
>~~~
>brew --version
>~~~
> - For more information, visit [https://brew.sh/](https://brew.sh/)
>
> **Why Homebrew?**
> - Easy installation of bioinformatics tools (samtools, bwa, BLAST, etc.)
> - Automatic dependency management
> - Simple updates with `brew update` and `brew upgrade`
>
> **Essential Tools via Homebrew (Recommended):**
> - After installing Homebrew, install essential development tools:
>~~~
>brew install git curl wget vim
>~~~
> - **git**: Version control system
> - **curl/wget**: Download tools for fetching files from the internet
> - **vim**: Powerful text editor for terminal
> - These tools may already be installed via Command Line Tools, but Homebrew versions are often more up-to-date
{: .solution}

> ## Linux
> **Office Suite:**
> - Install LibreOffice by going to [the installation page](https://www.libreoffice.org/download/libreoffice-fresh/). The version for Linux should automatically be selected. Click Download Version X.X.X (whichever is the most recent version). You will go to a page that asks about a donation, but you don't need to make one. Your download should begin automatically.
> - Once the installer is downloaded, double click on it and LibreOffice should install.
>
> **Terminal:**
> - Open terminal using Ctrl+Alt+T or search for "terminal", "xterm", or "uxterm" in your application menu
>
> **Code Editor / IDE (Required):**
> - **Visual Studio Code (Required):**
>   - **Installation via Snap (Recommended for Ubuntu 24.04):**
>~~~
>sudo snap install code --classic
>~~~
>   - **Alternative: Installation via .deb package:**
>     - Download from [https://code.visualstudio.com/](https://code.visualstudio.com/)
>     - Or install via terminal:
>~~~
>sudo apt update
>sudo apt install wget gpg
>wget -qO- https://packages.microsoft.com/keys/microsoft.asc | gpg --dearmor > packages.microsoft.gpg
>sudo install -D -o root -g root -m 644 packages.microsoft.gpg /etc/apt/keyrings/packages.microsoft.gpg
>sudo sh -c 'echo "deb [arch=amd64,arm64,armhf signed-by=/etc/apt/keyrings/packages.microsoft.gpg] https://packages.microsoft.com/repos/code stable main" > /etc/apt/sources.list.d/vscode.list'
>rm -f packages.microsoft.gpg
>sudo apt update
>sudo apt install code
>~~~
>   - Launch VS Code by typing `code` in terminal or searching for "Visual Studio Code" in your application menu
>   - **Recommended Extensions for Bioinformatics:**
>     - Python (Microsoft)
>     - R (REditorSupport)
>     - Remote - SSH (for HPC connections)
>     - Jupyter (for notebook support)
>   - **Why VS Code for Linux?**
>     - Free and open source
>     - Integrated terminal
>     - Perfect for Python, R, and Bash scripting
>     - Git integration
>     - Remote development capabilities
>
> **System Update (Required):**
> - Before installing any software, update your system. Open Terminal and run:
>~~~
>sudo apt update
>sudo apt upgrade -y
>~~~
> - The first command updates package lists
> - The second command upgrades installed packages
> - You may be prompted for your password
>
> **Build Essential Tools (Required):**
> - Install build-essential package for compiling software:
>~~~
>sudo apt install build-essential -y
>~~~
> - This package includes:
>   - GCC/G++ compilers
>   - Make and other build tools
>   - Standard C/C++ libraries
> - Verify installation:
>~~~
>gcc --version
>make --version
>~~~
>
> **Additional Recommended Packages:**
>~~~
>sudo apt install git curl wget vim -y
>~~~
> - **git**: Version control system
> - **curl/wget**: Download tools
> - **vim**: Text editor for terminal
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
