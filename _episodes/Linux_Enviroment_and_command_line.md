---
layout: page
title: 2_Linux Environment and Command Line
published: true
---

{% include gh_variables.html %}

> ## Reading and Watching
>
> - Reading 1: [Luscombe et al., 2001](http://archive.gersteinlab.org/papers/e-print/whatis-mim/text.pdf)  
> - Reading 2: [Attwood 2000](https://science.sciencemag.org/content/290/5491/471)  
> - Reading 3: [Smith 2018](https://www.embopress.org/doi/full/10.15252/embr.201846262) 
> - [File Permission and chmod](https://www.youtube.com/watch?v=3gcSeDoQ_rU)
> - [Nano Text Editor Basics](https://www.youtube.com/watch?v=Jf0ZJZJ8jlI)
{: .callout}

> ## Assignments
>
> Please complete the following assignments before the due date:  
> - [DataCamp: Introduction to Shell](https://app.datacamp.com/learn/courses/introduction-to-shell)
> - [Data Processing in Shell](https://app.datacamp.com/learn/courses/data-processing-in-shell)
{: .prereq}

![bioinformatics_DNA]({{site.baseurl}}/fig/DNA.jpg)

## History of UNIX 
Unix was conceived and implemented in 1969 at **AT&T's Bell Laboratories** by Ken Thompson, Dennis Ritchie, Douglas McIlroy, and Joe Ossanna. Initially released in 1971 and written in assembly language, Unix was re-written in C in 1973 by Dennis Ritchie, making it more portable across platforms. A legal issue forced AT&T to license the source code, which led to its widespread adoption in academia and industry. In 1984, AT&T began selling Unix as a proprietary product after divesting Bell Labs.

## Why UNIX
- Unix is historically significant in computing.
- Two dominant OS families: Unix-based and Windows-based.
- Widely used in back-end systems and personal computing.
- Unix derivatives like Linux are open source and community-developed.
- Skills learned in Unix transfer easily to other platforms.

![Unix_family tree]({{site.baseurl}}/fig/unix-simple.png)

## The Kernel
The kernel is the core of the operating system, managing memory, time, file storage, and communications in response to system calls.

## BSD
BSD (Berkeley Software Distribution) is a Unix variant developed at UC Berkeley. Derivatives like FreeBSD, OpenBSD, and NetBSD have emerged from BSD. OS X (macOS) and PS4 also have roots in BSD.

## LINUX  
Linux, released by Linus Torvalds in 1991, is a Unix-like, open-source operating system. Initially built for Intel x86 PCs, Linux has since been ported to more platforms than any other OS. It is now widely used on servers, supercomputers, mobile phones (Android), and gaming consoles like the Nintendo Switch.

- Linux means *Linus's MINIX*. Linus Torvalds is known for his frank communication style.
- Famous quotes from Linus:
  - "Talk is cheap. Show me the code."
  - "What am I going to do without my coffee maker? I'm going to sit here in a corner, crying, that's what."
  
![Linux_family tree](https://aerojsoft.files.wordpress.com/2016/02/linus-distribution-family-tree.jpg)

[Linux Family Tree](https://en.wikipedia.org/wiki/List_of_Linux_distributions)

## Unix/Linux Main Components
### Unix/Linux systems consist of three main parts:
- **User Space**: User-accessible programs, libraries, and utilities.
- **Kernel Space**: Manages interactions between user actions and hardware.
- **Hardware**: Physical components like CPU, memory, and I/O devices.

![anatomy]({{site.baseurl}}/fig/anatomy.jpg)


## What is GNU?
GNU is a free operating system that respects users' freedom. It is Unix-like but contains no Unix code. The GNU project was started by Richard Stallman, and the GNU General Public License ensures software freedom to use, modify, and share.

## How about macOS (XNU)?
macOS's kernel is XNU (XNU is Not Unix), a hybrid of the Mach kernel and BSD components. While macOS and Linux may seem similar, they have distinct histories and features.

![GNU](http://www.linuxandubuntu.com/wp-content/uploads/2019/07/What-is-GNU-in-GNULinux.jpg)
GNU and Tux

![MacOS](https://upload.wikimedia.org/wikipedia/commons/thumb/f/f2/Diagram_of_Mac_OS_X_architecture.svg/1280px-Diagram_of_Mac_OS_X_architecture.svg.png){: width="50%" height="50%"}

### Operating Systems Tasks
OS tasks include managing file systems, device I/O, processes, memory management, and more.

![OS]({{site.baseurl}}/fig/OS.png)

### I/O
I/O (Input/Output) refers to the communication between a system and the outside world, such as with a human or another processing system.

### The Shell
The shell is an interface between the user and the kernel, interpreting commands and arranging their execution.

### Shell Types
Different shells have unique features:
- **Bourne Shell (sh)**
- **Korn Shell (ksh)**
- **Bourne Again Shell (bash)** - Default in most Linux distributions.
- **C Shell (csh)**
- **TENEX/TOPS C Shell (tcsh)**

![shell types](https://d1jnx9ba8s6j9r.cloudfront.net/blog/wp-content/uploads/2019/05/Evolution-of-Linux-Shells-Types-of-Shells-in-Linux-Edureka.png)

**Bourne Shell** was created in the mid-1970s by Stephen R. Bourne.

### BASH
Bash (Bourne Again Shell) offers command-line editing, job control, and more, making it a powerful interactive shell.

![BASH](https://miro.medium.com/proxy/0*L0nhgi_19dlQJtzb.png)

## Text Editor Options
Common text editors in Unix/Linux:
- **nano**
- **emacs**
- **vim**

![text editor]({{site.baseurl}}/fig/texteditor.png){: width="150%" height="150%"}

### Using nano
Nano, created in 1999, is a free replacement for Pico and includes features like colored text and multiple buffers.

### Using emacs
GNU Emacs is a highly customizable text editor with an extensive feature set, including syntax highlighting and a built-in tutorial.

### Using vim
Vi is the standard Unix text editor and a powerful tool for text manipulation. Vim (Vi Improved) adds more features.

### Text Editor Cheat Sheets
- [VIM](https://preview.redd.it/ve1jv3m3qqj21.png?width=960&crop=smart&auto=webp&s=deb6dc83a462dc54523d703574e953638598af19)
- [nano](https://www.cheatography.com/bipinthite/cheat-sheets/nano-editor/)
- [Emacs](https://sachachua.com/blog/wp-content/uploads/2013/05/How-to-Learn-Emacs-v2-Large.png)


#### go to nano
```bash
#!/bin/bash
echo "hello world"
```
#### Save it "first.sh"

####
```bash
sh first.sh
```

### Shebang line
A shebang line (e.g., `#!/bin/bash`) at the top of a script tells the OS which interpreter to use for executing the file.

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

> ## pwd
> Returns the `p`resent `w`orking `d`irectory (print working directory).
> ```bash
> $ pwd
> ```
> ```output
> /home/{username}
> ```
> This means you are now working in the home directory, located under your username. You can also avoid writing the full path by using `~` in front of your username, or simply `~`.
{: .keypoints}

> ## mkdir
> To create a directory, the `mkdir` (make directory) command can be used.
> ```bash
> $ mkdir <DIRECTORY>
> ```
> Example:
> ```bash
> $ mkdir bch709_test
> ```
> Unlike PC/Mac folders, Unix allows spaces in directory names, using special characters. You can also specify the path where you want to create the new folder:
> ```bash
> $ mkdir bch709\ test
> ```
> Like most Unix commands, `mkdir` supports command-line options. For example, the `-p` option allows you to create parent directories in one step:
> ```bash
> $ mkdir -p bch709_test/help
> ```
{: .keypoints}

> ## cd
> The `cd` command stands for `c`hange `d`irectory.
> ```bash
> $ cd <DIRECTORY>
> ```
> Example:
> ```bash
> $ cd bch709_test
> ```
> To check your current directory, use:
> ```bash
> $ pwd
> ```
> ```output
> /home/{username}/bch709_test
> ```
> To move up one level to the parent directory:
> ```bash
> $ cd ..
> ```
> To move up two levels to the grandparent directory:
> ```bash
> $ cd ../../
> ```
> To navigate to the root directory:
> ```bash
> $ cd /
> ```
> To return to your home directory:
> ```bash
> $ cd ~
> ```
> To navigate to a specific directory:
> ```bash
> $ cd ~/bch709_test/help
> ```
> To move up one level again:
> ```bash
> $ cd ../
> ```
{: .keypoints}

> ## Absolute and relative paths
> The `cd` command allows you to change directories relative to your current location. However, you can also specify an absolute path to navigate directly to a folder.
> - **Absolute path**:
> ```bash
> $ cd /home/<username>/bch709_test/help
> ```
> - **Relative path**:
> ```bash
> $ cd ../
> ```
{: .checklist}

> ## The way back home
> Your home directory is `/home/<username>`. To quickly return to your home directory from any location, you can use:
> ```bash
> $ cd ~
> ```
> Or simply:
> ```bash
> $ cd
> ```
{: .checklist}

> ## Pressing `<TAB>` key
> The `<TAB>` key helps auto-complete file or directory names when typing in the terminal. If multiple matches exist, pressing `<TAB>` twice will show all possible options. This feature saves time and keystrokes.
{: .checklist}

> ## Pressing `<Arrow up/down>` key
> You can recall previously used commands by pressing the `up/down` arrow keys, or view your command history using the `history` command in the terminal.
{: .checklist}

> ## Copy & Paste
> In most terminal environments, dragging text automatically copies it, and right-clicking pastes it. On macOS, you can use `Command + C` and `Command + V` depending on your settings.
{: .checklist}

> ## ls
> The `ls` command lists the contents of a directory.
> ![ls]({{site.baseurl}}/fig/ls.png)
> The `ls` command has various useful options:
> - List files in long format:
> ```bash
> $ ls -l
> ```
> - List files sorted by creation time:
> ```bash
> $ ls -t
> ```
> - List files sorted by size:
> ```bash
> $ ls -S
> ```
> - List all files, including hidden files:
> ```bash
> $ ls -a
> ```
> You can also try combining options:
> ```bash
> $ ls -ltr
> ```
> To list files in another directory:
> ```bash
> $ ls -l /usr/bin/
> ```
> ![ls]({{site.baseurl}}/fig/ls2.png)
{: .keypoints}

> ## man & help
> Every Unix command comes with a manual, accessible via the `man` command. To view the manual for `ls`:
> ```bash
> $ man ls
> ```
> Alternatively, many commands also support a `--help` option:
> ```bash
> $ ls --help
> ```
{: .keypoints}

> ## rmdir
> The `rmdir` (remove directory) command deletes empty directories.
> ```bash
> $ rmdir <DIRECTORY>
> ```
> Example:
> ```bash
> $ rmdir help
> ```
> ![ls4]({{site.baseurl}}/fig/ls4.png)
> Note: You must be outside the directory to remove it with `rmdir`.
{: .keypoints}

> ## touch
> The `touch` command creates an empty file. To use it:
> ```bash
> $ touch test.txt
> $ touch exam.txt
> $ touch ETA.txt
> ```
> To list files, use:
> ```bash
> $ ls
> ```
> ![ls5]({{site.baseurl}}/fig/ls5.png)
{: .keypoints}

> ## mv
> The `mv` (move) command moves files or directories from one location to another.
> ```bash
> $ mv test.txt Hello
> ```
> You can also use wildcards like `*` to move multiple files:
> ```bash
> $ mv *.txt Hello
> ```
> The `*` wildcard matches any sequence of characters. This allows you to move files that follow a particular pattern.
{: .keypoints}

> ## Renaming files with mv
> The `mv` command can also be used to rename files:
> ```bash
> $ mv Hello/riches Hello/rags
> ```
{: .keypoints}

> ## Moving directories with mv
> You can move directories in the same way you move files:
> ```bash
> $ mv Hello/bch709_test .
> ```
> Here, `.` represents the current directory.
{: .keypoints}

> ## rm
> The `rm` (remove) command deletes files. Be cautious when using it, as deleted files cannot be recovered.
> To make `rm` safer, use the `-i` option for interactive deletion:
> ```bash
> $ rm -i ETA.txt exam.txt rags
> ```
> This will ask for confirmation before deleting each file.
> ![ls8]({{site.baseurl}}/fig/ls8.png)
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

> ## Clean up and start new
> Before we start the next session, let’s clean up the existing folders.
> ```bash
> $ cd ~/
> $ ls 
> ```
> Please do this first, then check the solution below.
{: .challenge}

> ## Clean the folder with contents
> How can we clean up the `bch709_test` and `Hello` folders?
> ```bash
> $ rm -R <folder_name>
> ```
> This command will recursively remove the specified folder and its contents.
{: .solution}

> ## Downloading a file
> Let's download a file from a website. There are several commands to download files, such as `wget`, `curl`, and `rsync`. In this case, we will use `curl`.
> File location:
> ```
> https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/bch709_student.txt
> ```
> How to use `curl`:
> ![curl]({{site.baseurl}}/fig/curl.png)
>
> To download the file using `curl`, use the following syntax:
> ```bash
> curl -L -o <output_name> <link>
> ```
> Example:
> ```bash
> curl -L -o bch709_student.txt https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/bch709_student.txt
> ls bch709_student.txt
> ```
{: .keypoints}


> ## Viewing file contents
> There are various commands to print the contents of a file in bash. Each command is often used in specific contexts, and when executed with filenames, they display the contents on the screen. Common commands include `less`, `more`, `cat`, `head`, and `tail`.
>
> - **`less` FILENAME**: Try this: `less bch709_student.txt`. It displays the file contents with line scrolling (use arrow keys, PgUp/PgDn, space bar, or Enter to scroll). Press `q` to exit.
> - **`more` FILENAME**: Try this: `more bch709_student.txt`. Similar to `less`, but you scroll using only the space bar or Enter. Press `q` to exit.
> - **`cat` FILENAME**: Try this: `cat bch709_student.txt`. This command displays the entire file content at once, which may result in the file scrolling off the screen for large files.
> - **`head` FILENAME**: Try this: `head bch709_student.txt`. It shows only the first 10 lines by default, but you can specify a different number using the `-n` option (e.g., `head -n 20`).
> - **`tail` FILENAME**: Try this: `tail bch709_student.txt`. It displays the last 10 lines by default, and similar to `head`, you can modify the number of lines with `-n`. 
{: .keypoints}

> ## How many lines does the file have?
> You can pipe the output of your stream into another program rather than displaying it on the screen. Use the `|` (Vertical Bar) character to connect the programs. For example, the `wc` (word count) program:
> ```bash
> cat bch709_student.txt | wc
> ```
> prints the number of lines, words, and characters in the stream:
> ```output
>     106     140    956
> ```
> To count just the lines:
> ```bash
> cat bch709_student.txt | wc -l
> ```
> ```output
> 106
> ```
> *Of course, we can use `wc` directly:*
> ```bash
> wc -l bch709_student.txt 
> ```
> This is equivalent to:
> ```bash
> cat bch709_student.txt | wc -l
> ```
> In general, it is better to open a stream with `cat` and then pipe it into the next program. This method simplifies building and understanding more complex pipelines.
>
> Let's also check the first few lines of the file using `head`:
> ```bash
> cat bch709_student.txt | head
> ```
> ![head]({{site.baseurl}}/fig/head.png)
> 
> Is this equivalent to running?
> ```bash
> head bch709_student.txt
> ```
> {: .solution}
{: .checklist}

> ## grep 
> `grep` (Global Regular Expression Print) is one of the most useful commands in Unix. It is commonly used to filter a file/input, line by line, against a pattern. It prints each line of the file containing a match for the pattern.
> Check available options with:
> ```bash
> grep --help
> ```
> Syntax:
> ```
> grep [OPTIONS] PATTERN FILENAME
> ```
> Let's find how many people use macOS. First, check the file:
> ```bash
> wc -l bch709_student.txt
> ```
> ```bash
> less bch709_student.txt
> ```
> To find the macOS users:
> ```bash
> cat bch709_student.txt | grep MacOS
> ```
> To count the number of macOS users:
> ```bash
> cat bch709_student.txt | grep MacOS | wc -l
> ```
> Alternatively:
> ```bash
> grep MacOS bch709_student.txt | wc -l
> ```
> Or simply:
> ```bash
> grep -c MacOS bch709_student.txt
> ```
> Using flags to filter lines that don’t contain "Windows":
> ```bash
> grep -v Windows bch709_student.txt
> ```
> Combining multiple flags:
> ```bash
> grep -c -v Windows bch709_student.txt
> ```
> With case-insensitive search and colored output:
> ```bash
> grep --color -i macos bch709_student.txt
> ```
{: .checklist}

> ## How do I store the results in a new file?
> Use the `>` character for redirection:
> ```bash 
> grep -i macos bch709_student.txt > mac_student
> ```
> ```bash
> cat bch709_student.txt | grep -i windows > windows_student
> ```
> You can check the new files with `cat` or `less`.
{: .checklist}

> ## Do you want to check the differences between two files?
> ```bash
> diff mac_student windows_student
> ```
> ```bash
> diff -y mac_student windows_student
> ```
{: .keypoints}

> ## How can I select the name only? (cut)
> To extract specific columns, use `cut`:
> ```bash
> cat bch709_student.txt | cut -f 1
> ```
> Or:
> ```bash
> cut -f 1 bch709_student.txt
> ```
{: .keypoints}

> ## How can I sort it? (sort)
> Sorting names:
> ```bash
> cut -f 1 bch709_student.txt | sort
> ```
> Sorting by the second field:
> ```bash
> sort -k 2 bch709_student.txt
> ```
> Sorting by the first field:
> ```bash
> sort -k 1 bch709_student.txt
> ```
> Sorting and extracting names:
> ```bash
> sort -k 1 bch709_student.txt | cut -f 1
> ```
> Save the sorted names:
> ```bash
> sort -k 1 bch709_student.txt | cut -f 1 > name_sort
> ```
{: .keypoints}

> ## uniq
> The `uniq` command removes duplicate lines from a **sorted file**, retaining only one instance of matching lines. Optionally, it can show lines that appear exactly once or more than once. Note that `uniq` requires sorted input.
> ```bash
> cut -f 2 bch709_student.txt > os.txt
> ```
> To count duplicates:
> ```bash
> uniq -c os.txt
> ```
> ![uniq]({{site.baseurl}}/fig/uniq.png)
> To sort and then count:
> ```bash
> sort os.txt | uniq -c
> ```
> ![uniq2]({{site.baseurl}}/fig/uniq2.png)
> Of course, you can sort independently:
> ```bash
> sort os.txt > os_sort.txt
> uniq -c os_sort.txt
> ```
> To learn more about `uniq`:
> ```bash
> uniq --help
> ```
> ![uniq3]({{site.baseurl}}/fig/uniq3.png)
{: .checklist}

> ## diff
> The `diff` command compares the differences between two files.
> Example usage:
> ```bash
> diff FILEA FILEB
> ```
> When trying new commands, always check with `--help`:
> ![diff]({{site.baseurl}}/fig/diff.png)
> Compare the contents of `os.txt` and `os_sort.txt`:
> ```bash
> diff os.txt os_sort.txt
> ```
> Side-by-side comparison:
> ```bash
> diff -y os.txt os_sort.txt
> ```
> You can also pipe sorted data into `diff`:
> ```bash
> sort os.txt | diff -y - os_sort.txt
> ```
{: .checklist}


### There are still a lot of command that you can use. Such as `paste`, `comm`, `join`, `split` etc.
> ## Let's Download Bigger Data
> Please visit the following website to download larger data:
> ```
> https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
> ```
> Download the file marked in yellow.
> ![download]({{site.baseurl}}/fig/download.png)
{: .keypoints}

> ## How to Download
> ![download2]({{site.baseurl}}/fig/download2.png)
> To download a file from a URL, use the following command:
> ```bash
> curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
> gunzip mrna.fa.gz
> ```
{: .solution}

> ## Viewing File Contents
> What is the difference between `less`, `tail`, `more`, and `head`?
{: .discussion}

> ## What Are Flags (Parameters/Options)?
> A "flag" in Unix terminology is a parameter that is added to a command to modify its behavior. For example:
> ```bash
> $ ls
> ```
> versus
> ```bash
> $ ls -l
> ```
> The `-l` flag changes the output of `ls` to display detailed information about each file. Flags help commands behave differently or report information in various formats.
{: .checklist}

> ## How to Find Available Flags
> You can use the manual (man) to learn more about a command and its available flags:
> ```bash
> $ man ls
> ```
> This will display the manual page for `ls`, detailing all available flags and their purposes.
> ![man]({{site.baseurl}}/fig/man.png)
{: .checklist}

> ## What if the Tool Doesn’t Have a Manual Page?
> Not all tools include a manual page, especially third-party software. In those cases, use the `-h`, `-help`, or `--help` options:
> ```bash
> $ curl --help
> ```
> If these flags do not work, you may need to refer to external documentation or online resources.
{: .checklist}

> ## What Are Flag Formats?
> Unix tools typically follow two flag formats:
> - **Short form**: A single minus `-` followed by a single letter, like `-o`, `-L`.
> - **Long form**: Double minus `--` followed by a word, like `--output`, `--Location`.
> 
> Flags may act as toggles (on/off) or accept additional values (e.g., `-o <filename>` or `--output <filename>`).
> Some bioinformatics tools diverge from this format and use a single `-` for both short and long options (e.g., `-g`, `-genome`).
> 
> **Using flags is essential in Unix, especially in bioinformatics, where tools rely on a large number of parameters. Proper flag usage ensures the accuracy of results.**
{: .checklist}

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
$ ls mrna.fa
```

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
```bash
http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz
```
What is `.gz` ?
```bash
file MGI.gff3.gz
```

### Compression
There are several options for archiving and compressing groups of files or directories. Compressed files are not only easier to handle (copy/move) but also occupy less size on the disk (less than 1/3 of the original size). In Linux systems you can use zip, tar or gz for archiving and compressing files/directories.

>## ZIP compression/extraction
>```bash
>zip OUTFILE.zip INFILE.txt # Compress INFILE.txt  
>zip -r OUTDIR.zip DIRECTORY # Compress all files in a DIRECTORY into one archive file (OUTDIR.zip)  
>zip -r OUTFILE.zip . -i \*.txt # Compress all txt files in a DIRECTORY into one archive file (OUTFILE.zip)   
>unzip SOMEFILE.zip
>```
{: .checklist}

## TAR Compression and Extraction

The `tar` (tape archive) utility is used to bundle multiple files into a single archive file and to extract individual files from that archive. It offers options for automatic compression and decompression, along with special features for incremental and full backups.

### Common Commands
- To extract the contents of a gzipped TAR file:
  ```bash
  tar -xzvf SOMEFILE.tar.gz
  ```
- To create a gzipped TAR archive from a directory:
  ```bash
  tar -czvf OUTFILE.tar.gz DIRECTORY
  ```
- To archive and compress all `.txt` files in the current directory:
  ```bash
  tar -czvf OUTFILE.tar.gz *.txt
  ```
- To create a backup archive of a specific directory:
  ```bash
  tar -czvf backup.tar.gz BACKUP_WORKSHOP
  ```

## Gzip Compression and Extraction

The `gzip` (GNU zip) compression utility is designed as a replacement for the `compress` program, offering much better compression without using patented algorithms. It is the standard compression system for all GNU software.

### Commands

- To compress a file:
  ```bash
  gzip SOMEFILE  # This also removes the uncompressed file
  ```

- To uncompress a file:
  ```bash
  gunzip SOMEFILE.gz  # This also removes the compressed file
  ```

### Example

Compress the file `MGI.gff3.gz` and examine the size. Then, uncompress it so that you can use this file for later exercises.

```bash
$ gunzip MGI.gff3.gz
$ ls -lh
$ gzip MGI.gff3
$ ls -lh
$ gunzip MGI.gff3.gz
```
{: .checklist}


## GFF3 Annotations


Print all sequences annotated in a GFF3 file.
```bash
cut -s -f 1,9 MGI.gff3 | grep $'\t' | cut -f 1 | sort | uniq
```

Determine all feature types annotated in a GFF3 file.
```bash
grep -v '^#' MGI.gff3 | cut -s -f 3 | sort | uniq
```

Determine the number of genes annotated in a GFF3 file.
```bash
grep -c $'\tgene\t' MGI.gff3
```

Extract all gene IDs from a GFF3 file.
```bash
grep $'\tgene\t' MGI.gff3 | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'
```

Print all CDS.
```bash  
cat MGI.gff3 | cut -f 3 | grep CDS | 
```
Print CDS and ID
```bash
cat MGI.gff3 | cut -f 1,3,4,5,7,9 | head
cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | head
cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | sed 's/;.*//g' | head
cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | sed 's/;.*//g' | sed 's/ID=//g' | head
cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep $'\tCDS\t' | sed 's/;.*//g' | sed 's/ID=//g' | head

```
Print length of each gene in a GFF3 file.
```bash
grep $'\tgene\t' MGI.gff3 | cut -s -f 4,5 | perl -ne '@v = split(/\t/); printf("%d\n", $v[1] - $v[0] + 1)'
```

Extract all gene IDs from a GFF3 file.
```bash
grep $'\tgene\t' MGI.gff3 | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'
```

Time and again we are surprised by just how many applications it has, and how frequently problems can be solved by sorting, collapsing identical values, then resorting by the collapsed
counts. 
The skill of using Unix is not just that of understanding the commands themselves. It is more about recognizing when a pattern, such as the one that we show above, is the solution to the problem that you wish to solve. The easiest way to learn to apply these patterns is by looking at how others solve problems, then adapting it to your needs.

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


## Cloud Platform Choices

There are several cloud providers to choose from. Some scientific clouds may either be free or allocate resources competitively. Commercial clouds can be very powerful, but choice can be overwhelming. Availability of public and commercial cloud resources also vary by country and region.

The major tradeoff between platforms is between flexibility and cost. Generally speaking, services that allow you more flexibility and autonomy will be more expensive than highly managed services.

Below, we have highlighted three types of computing resources: Clusters, Open Science Clouds, and Commercial Clouds, which are widely available to researchers around the world. However, the availability of any specific cluster or cloud may be region-specific, and this is not meant to be an exhaustive list. We encourage researchers to use this list as a starting point for learning about cloud resources and suggest checking with your local or regional government to see what other Open Science Clouds might be available to you. The cloud resources listed here should be available to any scientist based in the US but may be unavailable or have different pricing in other countries.

### University/Corporate Computing Clusters

Many universities and businesses operate their own computing clusters that are available to students and staff at low or no cost. If your employer maintains a computing cluster, this will almost always be the least expensive option.

However, most HPCCs (High Performance Computing Clusters) put limits on:
- The number of processors a user can utilize at once
- The amount of disk storage per user
- The amount of time a single process can run
- What programs can be installed, and by whom
- Who can have accounts and access data

HPCCs are also a shared resource, so even when you have access, your programs are unlikely to run immediately. Most HPCCs run some kind of scheduler to which you submit your processing jobs, and it runs them as resources become available. In order to submit a job, you generally will need to know not only what program you want to run, but with how many processors and for how long. While interacting with the scheduler is no more difficult than interacting with the shell, it will have its own set of commands and syntax that you'll need to learn; these vary widely among HPCCs.

There are also many upsides to using an HPCC. As previously mentioned, they're generally the least expensive option, but they often come with more perks. For instance, many HPCCs offer free or low-cost training, storage space, backup options, and technical support. If your HPCC has a scheduler, you can also queue up many sequential jobs all at once, and you don't have to worry about racking up fees on instances that are sitting idle. It's often also much easier to pay for HPCC use than to pay for Amazon using grant money; however, universities are getting better about AWS payments.

### Open Science Clouds

#### [XSEDE](https://www.xsede.org/)

The Extreme Science and Engineering Discovery Environment (XSEDE) is an NSF-funded HPCC, so it is open to any US-based researcher and shares most of the same benefits and drawbacks of a university or corporate HPCC. If your university or corporation doesn't have its own HPCC resources, XSEDE will likely be your cheapest option.

Although any US-based researcher can use XSEDE, first [they'll need an account](https://portal.xsede.org/#/guest). Like the HPCC options described above, XSEDE uses a scheduler to start jobs and puts limits on how many resources any one user can utilize at once.

XSEDE can also be a bit intimidating at first because you will need to know what resources you need and for how long before you get started. XSEDE runs like a mini version of the NSF grant system. In order to qualify to submit large jobs, you'll have to submit a [allocation request](https://portal.xsede.org/allocations/research) in the form of a short proposal. Also like an NSF grant, if your proposal is accepted, that means you have access to whatever resources you were approved for, for the time frame you requested.

Don't let that paragraph scare you off, though. XSEDE has two different allocation tracks. If you aren't sure exactly what you'll need for your big project, you can request a [startup allocation](https://portal.xsede.org/allocations/startup), which only requires an abstract rather than a proposal, and grants you a year to try out your new pipeline or analysis. These are usually granted in a week or so and are intended for you to test your pipeline so you know what to ask for in your allocation proposal.

If that still sounds a little too daunting, XSEDE also has [trial allocations](https://iujetstream.atlassian.net/wiki/spaces/JWT/pages/76149919/Jetstream+Trial+Access+Allocation) which give you access to only a tiny fraction of XSEDE's power but are plenty large enough to test your code and see if a larger allocation is worth pursuing. These allocations are granted more or less immediately by simply filling in a form and agreeing to the usage rules.

If you're interested in using XSEDE, check to see if your workplace has a [Campus Champion](https://www.xsede.org/community-engagement/campus-champions). These are people who have had extensive training on both the XSEDE system and the allocation program and can help you figure out how to apply and what you need.

#### [Open Science Grid](https://opensciencegrid.org)

The Open Science Grid (OSG) is an NSF-funded national network of computing centers that have pooled their resources together and made them available to various research groups. The OSG is usable by any researcher based at a US institution and is accessible for free without an allocation. It can provide millions of computing hours for researchers who have problems that fit well on its setup.

Certain projects and universities have direct access to the Open Science Grid, but any researcher can access it through the [OSG Connect](https://osgconnect.net/) entry point. If you apply for OSG access through that website, you will have a consultation with someone who can help you determine if your analysis is a good fit and how to get started.

The OSG is a great fit for problems that can be broken into lots of independent pieces. One good example is read alignment: the starting read data can be broken into several pieces, each of them aligned, and then the results combined. Another good problem type for the OSG are multi-start simulations or statistical analyses where you need to run the same model or simulation many, many times. The payoff of using this approach is being able to run on many hundreds (sometimes thousands!) of computers at once, accelerating your analysis.

Note that you don't access a specific computing center through OSG -- unlike XSEDE, where you apply for time and then run on a specific HPCC resource, the OSG sits on top of many resources, and when you submit your work, it could run almost anywhere in the overall system.

#### [Open Science Data Cloud (OSDC)](https://www.opensciencedatacloud.org/)

The Open Science Data Cloud provides the scientific community with resources for storing, sharing, and analyzing terabyte and petabyte-scale scientific datasets. OSDC's Bionimbus Protected Data Cloud (PDC) is a platform designed with the sole purpose of analyzing and sharing protected genomics data.

#### [Atmosphere](https://pods.iplantcollaborative.org/wiki/display/atmman/Getting+Started)

#### [CyVerse (iPlant Collaborative) Atmosphere](http://www.cyverse.org/atmosphere)

#### [JetStream](http://jetstream-cloud.org/)

### Commercial Clouds

Computing architecture is moving (albeit at a slow pace) to the Model-to-Data paradigm. This means that scientists should be encouraged to bring their compute to where the data is stored, instead of the other way around. The following outlines the general differences between the three major commercial cloud providers: Amazon Web Services (AWS), Google Cloud Platform (GCP), and Microsoft Azure.

Essentially all cloud providers provide extremely similar computing and storage options; you can "rent" or provision computing infrastructure with very similar specifications across all three cloud vendors. Even the costs are highly comparable. What governs how to choose the right cloud computing vendor is highly opportunistic: (1) funding options, (2) solidarity with collaborating/similar scientific groups, (3) location of datasets that a particular research group works with, and (4) familiarity with cloud vendor services.

1. **Funding options**: Does your grant stipulate where you should build your computing pipeline? For example, the NIH often partners with specific cloud vendors to provide cloud credits that allow researchers to compute for free. Some cloud vendors also provide research credits.
   
2. **Solidarity with collaborating/similar scientific groups**: Are other research groups in your field drawn to a specific cloud vendor? It might make sense to utilize the same cloud service to minimize transfer (egress) costs, especially if you are sharing large datasets. You may also be able to make use of existing pipelines without reinventing the wheel if you choose the same cloud provider that your collaborators are using.
   
3. **Location of datasets that a particular research group works with**: Again, thinking of bringing your models to where the data is stored helps minimize costs and saves you time in having to download and store data separately.
   
4. **Services**: Here, services refer to cloud vendor add-ons that take away the need for a user to set up their own computing infrastructure. A fully managed database (e.g., AWS RDS, GCP CloudSQL, Azure SQL DB) is an example of a service. If you are used to SQL Server, you may want to look into options provided by Azure. Are you more familiar with Postgres SQL? Then AWS and GCP might provide cheaper options for you.

#### [Amazon EC2](http://aws.amazon.com/ec2/)

The Amazon Web Service (AWS) that you've been using is the Elastic Compute (EC2) cloud. There are actually lots of other cloud and storage solutions under the AWS umbrella, but when most data scientists say AWS, they mean [EC2](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html). With EC2, you can rent access to a cloud computing resource as small as your laptop or as large as a 64-processor machine with 488GB of memory, and with a number of different operating systems. These instances can be optimized for jobs that are memory intensive, require a lot of bandwidth, or [almost any other specific need](https://aws.amazon.com/ec2/instance-types/). There are so many options that we can't cover them all here, but these are a few popular ones:

##### On-Demand

All this variety and optimization makes EC2 much more expensive than an average HPCC; however, depending on your needs, it can be [quite affordable](https://aws.amazon.com/ec2/pricing/). If you want to start an EC2 instance whenever you want and have instant access, you can rent a quite large on-demand instance with 8 processors and 32 GB of memory for ~40 cents an hour, and tiny instances are only about half a cent per hour.

##### Spot Instances

If your program can tolerate pauses and you don't need the analysis done as fast as possible, you can request a spot instance. Essentially, whenever Amazon has computing capacity that no one is paying them for, they lower the prices for renting some systems. If you request a spot-instance, that means you specify the size and parameters of your machine rental and set a limit on how much you're willing to spend per hour. Then, whenever the rental rate dips below your maximum, your instance turns on and runs until the price goes back up. If it's not an important shopping season, and you aren't in a hurry, you can often run spot-instances for less than half their normal cost.

##### Free Tier

There are also [free options](https://aws.amazon.com/free/), which allow you to test out the interface and make sure it will meet your needs before you start renting.

Just remember that with EC2 and all other commercial services, you're paying for renting the computer, whether you're using it or not. If you leave an instance on idly for months after your pipeline has finished, you'll still have to pay for that time.

#### [Google Cloud](https://cloud.google.com/): [Getting Started](https://cloud.google.com/compute/docs/quickstart)

GCP offers very competitive prices for compute and storage (as of July 2019, their compute pricing is lower than that of AWS and Azure for instances of comparable specifications). If you are looking to dabble in cloud computing but do not need a vast catalog of services, GCP would be a good place to start looking.

Their version of "Spot Instances" are known as pre-emptible instances and offer very competitive pricing. GCP also has TPUs.

#### [Microsoft Azure](https://azure.microsoft.com/en-us/)

If your software requires Microsoft Windows, it may be cheaper to use MS Azure due to licensing issues. Azure's computing instances are known as Azure Virtual Machines and often come at a slightly higher cost than other cloud computing vendors' offerings. If a lot of your computing pipeline is Windows dependent, it may make sense to build everything on MS Azure from the get-go.

#### [IBM Cloud](https://www.ibm.com/cloud)

IBM Cloud offers more than 11 million bare metal configurations in virtual mode which are customizable RAM and SSDs on bare metal. They also have on-demand provisioning for all servers, with management and monitoring included along with direct and cost-free tech support.

## How to Choose

As you can see, highly managed systems (HPCCs, XSEDE, etc.) usually are free or cheap, but relatively inflexible. There may be certain programs you can't install, or there may be long wait times. Commercial systems are generally more flexible because you can make them look however you want, but they can be quite expensive, especially if you run for a long time or have a lot of data. However, there are other things to consider.

Another way to think about this is not *whether* you want to spend your time and money, but *where* you want to spend them. AWS will let you install anything you want, but that also means that you'll have to spend some time up-front installing programs and testing configurations. Your HPCC jobs might take a week to start, but you don't have to do any of the systems administration, and you can work on other projects while your job sits in the queue.

Your familiarity with the pipeline can also be a factor. Let's say you want to run a program you've never run before. If you have no way to estimate how long your job will take, a service like AWS might save you money because you can run your program until it's done, no matter how long that is. Running it on an HPCC can be really frustrating because you'll need to submit your job with the amount of time and resources it will use. An HPCC will only let your job run for the amount of time you requested, so if you underestimate the amount of time it will take, even by one minute, the system kills your program, and you need to resubmit it. On the other hand, if you want to run that same program, but you can easily estimate the runtime, an HPCC is going to be a much cheaper choice.

In my work, I often use both commercial and non-commercial services. I tend to use AWS for testing, with small amounts of data, until I know how the program behaves. Then I port the pipeline to my university HPCC for running the full dataset.

> ## Discussion
>
> In small groups or on your own, plot out your next bioinformatics project. With guidance from your instructors and the above references, try to determine not only what types of resources you'll need but what platform will best suit your project.
>
> Some things to consider:
>
> - How much data do you have?
> - What computational steps will it need?
>   - What is the *largest* computational step?
>   - Can any steps be done in parallel?
> - What is your timeframe?
> - Who will be doing most of the computational work?
>   - What computational skills do they have?
>   - Do you need to share the data across many labs?
> - How many times will you need to run this pipeline?
{: .challenge}

> ## Human Genomic Data & Security
>
> Note that if you are working with human genomics data, there might be ethical and legal considerations that affect your choice of cloud resources to use. The terms of use, and/or the legislation under which you are handling the genomic data, might impose heightened information security measures for the computing environment in which you intend to process it. This is too broad a topic to discuss in detail here, but in general terms you should think through the technical and procedural measures needed to ensure that the confidentiality and integrity of the human data you work with is not breached. If there are laws that govern these issues in the jurisdiction in which you work, be sure that the cloud service provider you use can certify that they support the necessary measures. Also note that there might exist restrictions for the use of cloud service providers that operate in other jurisdictions than your own, either by how the data was consented by the research subjects or by the jurisdiction under which you operate. Do consult the legal office of your institution for guidance when processing human genomic data.
{: .callout}

### Other Resources

Learn more about cloud computing in bioinformatics:

Fusaro VA, Patil P, Gafni E, Wall DP, Tonellato PJ (2011) **Biomedical Cloud Computing With Amazon Web Services**. PLoS Comput Biol 7(8): e1002147. doi: [10.1371/journal.pcbi.1002147](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002147)

Charlebois K, Palmour N, Knoppers BM (2016) **The Adoption of Cloud Computing in the Field of Genomics Research: The Influence of Ethical and Legal Issues**. PLoS ONE 11(10): e0164347. https://doi.org/10.1371/journal.pone.0164347

Langmead B, Nellore A (2018) **Cloud Computing for Genomic Data Analysis and Collaboration**. Nature Reviews Genetics 19 (208). doi: 10.1038/nrg.2017.113 (https://www.nature.com/articles/nrg.2017.113)

This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.
