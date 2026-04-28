---
layout: page
title: 02_Linux Environment and Command Line
published: true
---

{% include gh_variables.html %}

> ## Reading and Watching
>
> - Reading 1: [Luscombe et al., 2001](http://archive.gersteinlab.org/papers/e-print/whatis-mim/text.pdf)
> - Reading 2: [Attwood 2000](https://doi.org/10.1126/science.290.5491.471)
> - Reading 3: [Smith 2018](https://doi.org/10.15252/embr.201846262)
> - [Nano Text Editor Basics](https://www.youtube.com/watch?v=Jf0ZJZJ8jlI)
> - [Edit files on command line](https://www.youtube.com/watch?v=y4SzAr33st0)
{: .callout}

> ## Before you start
> - Access to a Unix-like terminal (Linux, macOS, or WSL) with internet connectivity.
> - Permission to install software (either `sudo` on Linux or Homebrew on macOS); on managed systems, check with your admin first.
> - Run install commands on your own machine or assigned training server—not on shared clusters unless instructed.
> - Comfort with copy/paste in your terminal and a basic text editor (nano/vim/emacs) available.
{: .prereq}

> ## Assignments
>
> Please complete the following assignments before the due date:
> - [DataCamp: Introduction to Shell](https://app.datacamp.com/learn/courses/introduction-to-shell)
> - [Data Processing in Shell](https://app.datacamp.com/learn/courses/data-processing-in-shell)
> - [Introduction to Git](https://campus.datacamp.com/courses/introduction-to-git)
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

> ## Every Linux Concept Explained in 8 Minutes
> [![Every Linux Concept Explained](../fig/youtube-linux-concepts.jpg)](https://www.youtube.com/watch?v=gBuTIrEG87s)
>
> [Watch: Every LINUX Concept Explained in 8 Minutes](https://www.youtube.com/watch?v=gBuTIrEG87s)
{: .callout}

- Linux means *Linus's MINIX*. Linus Torvalds is known for his frank communication style.
- Famous quotes from Linus:
  - "Talk is cheap. Show me the code."
  - "What am I going to do without my coffee maker? I'm going to sit here in a corner, crying, that's what."

![Linux_family tree](../fig/linux-family-tree.jpg)

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

![GNU](../fig/gnu-linux.png)
GNU and Tux

![MacOS](../fig/macos-architecture.png){: width="50%" height="50%"}

### Operating Systems Tasks
OS tasks include managing file systems, device I/O, processes, memory management, and more.

![OS]({{site.baseurl}}/fig/OS.png)

### I/O
I/O (Input/Output) refers to the communication between a system and the outside world, such as with a human or another processing system.

### The Shell
The shell is an interface between the user and the kernel, interpreting commands and arranging their execution.

### Shell Types

> ## Bash in 100 Seconds
> [![Bash in 100 Seconds](../fig/youtube-bash-100sec.jpg)](https://www.youtube.com/watch?v=I4EWvMFj37g)
>
> [Watch: Bash in 100 Seconds](https://www.youtube.com/watch?v=I4EWvMFj37g)
{: .callout}

Different shells have unique features:
- **Bourne Shell (sh)**
- **Korn Shell (ksh)**
- **Bourne Again Shell (bash)** - Default in most Linux distributions.
- **C Shell (csh)**
- **TENEX/TOPS C Shell (tcsh)**

![shell types](../fig/shell-types-edureka.png)

**Bourne Shell** was created in the mid-1970s by Stephen R. Bourne.

### BASH
Bash (Bourne Again Shell) offers command-line editing, job control, and more, making it a powerful interactive shell.

![BASH](../fig/bash-medium.png)

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
- [VIM Cheatsheet](../fig/vim-cheatsheet.png)
- [nano](https://www.cheatography.com/bipinthite/cheat-sheets/nano-editor/)
- [Emacs Cheatsheet](../fig/emacs-cheatsheet.png)


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

In order to make it possible to execute scripts as though they were first class executables, UNIX systems will looks for what we refer to as a shebang line at the top of the file. The origin of the name is murky. Some think it came from sharp-bang or hash-bang – contractions of # ("sharp") and ! ("bang"). Others think the "SH" is in reference to the first UNIX shell, named "sh".

In any case, if a UNIX system sees that the first line of an executable begins with #!, then it will execute the file using whatever command is specified in the rest of the line. For example, if there's a file named /path/to/bar that looks like:
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

## Basic commands

|Category|command|
|---|---|
|Navigation| cd, ls, pwd|
|File creation|touch,nano,mkdir,cp,mv,rm,rmdir|
|Reading|more,less,head,tail,cat|
|Compression|zip,gzip,bzip2,tar,compress|
|Uncompression|unzip,gunzip,bunzip2,uncompress|
|Permissions|chmod|
|Help|man|

> ## Open Current Folder in File Manager
> You can open your current terminal directory in a graphical file manager:
>
> **Mac (iTerm2/Terminal):**
> ```bash
> open .
> ```
> This opens the current directory in Finder.
>
> **Windows WSL:**
> ```bash
> explorer.exe .
> ```
> This opens the current directory in Windows Explorer.
>
> **Accessing WSL Files from Windows Explorer:**
> 1. Open Windows Explorer
> 2. Type in the address bar:
> ```
> \\wsl$
> ```
> 3. You'll see your Linux distributions listed (e.g., `Ubuntu`)
> 4. Navigate to your home directory: `\\wsl$\Ubuntu\home\{username}`
>
> **Tip:** Pin `\\wsl$\Ubuntu\home\{username}` to Quick Access for easy access!
>
> You can also type this directly in Explorer's address bar:
> ```
> \\wsl.localhost\Ubuntu\home\{username}
> ```
{: .callout}

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
> $ cd ~
> $ mkdir bch709_test
> ```
> Like most Unix commands, `mkdir` supports command-line options. For example, the `-p` option allows you to create parent directories in one step:
> ```bash
> $ mkdir -p bch709_test/help
> ```
> Note: Unix allows spaces in directory names using `\` (e.g., `mkdir bch709\ test`), but this is not recommended.
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
> Example (assuming you're in `bch709_test` directory):
> ```bash
> $ cd ~/bch709_test
> $ rmdir help
> ```
> ![ls4]({{site.baseurl}}/fig/ls4.png)
> Note: The directory must be empty, and you must be outside the directory to remove it.
{: .keypoints}

> ## touch
> The `touch` command creates an empty file. Make sure you're in `bch709_test`:
> ```bash
> $ cd ~/bch709_test
> $ touch test.txt
> $ touch exam.txt
> $ touch ETA.txt
> ```
> To list files, use:
> ```bash
> $ ls
> ```
> ```output
> ETA.txt  exam.txt  test.txt
> ```
{: .keypoints}

> ## mv
> The `mv` (move) command moves files or directories from one location to another.
> ```bash
> $ mkdir Hello
> $ mv test.txt Hello
> $ ls Hello
> ```
> ```output
> test.txt
> ```
> You can also use wildcards like `*` to move multiple files:
> ```bash
> $ mv *.txt Hello
> $ ls Hello
> ```
> ```output
> ETA.txt  exam.txt  test.txt
> ```
> The `*` wildcard matches any sequence of characters. This allows you to move files that follow a particular pattern.
{: .keypoints}

> ## Renaming files with mv
> The `mv` command can also be used to rename files:
> ```bash
> $ mv Hello/test.txt Hello/renamed.txt
> $ ls Hello
> ```
> ```output
> ETA.txt  exam.txt  renamed.txt
> ```
{: .keypoints}

> ## Moving files back
> You can move files back to the current directory using `.`:
> ```bash
> $ mv Hello/renamed.txt .
> $ ls
> ```
> Here, `.` represents the current directory.
{: .keypoints}

> ## rm
> The `rm` (remove) command deletes files. Be cautious when using it, as deleted files cannot be recovered.
> To make `rm` safer, use the `-i` option for interactive deletion:
> ```bash
> $ rm -i renamed.txt
> ```
> ```output
> rm: remove regular empty file 'renamed.txt'? y
> ```
> This will ask for confirmation before deleting each file.
>
> To remove multiple files:
> ```bash
> $ rm Hello/ETA.txt Hello/exam.txt
> ```
{: .keypoints}

> ## cp
> The `cp` (copy) command copies files. Unlike `mv`, the original file remains at the source.
> ```bash
> $ touch file1
> $ cp file1 file2
> $ ls
> ```
> ```output
> file1  file2  Hello
> ```
>
> Copy a file into a directory:
> ```bash
> $ cp file1 Hello/
> $ ls Hello
> ```
> ```output
> file1
> ```
>
> Copy an entire directory using `-r` (recursive):
> ```bash
> $ cp -r Hello Hello_backup
> $ ls
> ```
> ```output
> file1  file2  Hello  Hello_backup
> ```
> Use `cp --help` or `man cp` for more options.
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
> Before we start the next session, let's clean up the files and folders we created.
> ```bash
> $ cd ~/bch709_test
> $ ls
> ```
> You should see: `Hello`, `Hello_backup`, `file1`, `file2`
{: .challenge}

> ## Clean the folder with contents
> To remove folders and their contents, use `rm -r` (recursive):
> ```bash
> $ cd ~
> $ rm -r bch709_test
> $ ls
> ```
> This removes `bch709_test` and everything inside it (Hello, Hello_backup, file1, file2).
>
> **Warning:** `rm -r` permanently deletes folders and all contents. Use with caution!
{: .solution}

> ## Downloading a file
> Let's download a file from a website. There are several commands to download files, such as `wget`, `curl`, and `rsync`. In this case, we will use `curl`.
>
> First, create a working directory:
> ```bash
> $ cd ~
> $ mkdir -p bch709_data
> $ cd bch709_data
> ```
>
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
> $ curl -L -o bch709_student.txt https://raw.githubusercontent.com/plantgenomicslab/BCH709/gh-pages/bch709_student.txt
> $ ls
> ```
> ```output
> bch709_student.txt
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
> Let's find how many people use Mac. First, check the file:
> ```bash
> wc -l bch709_student.txt
> ```
> ```bash
> less bch709_student.txt
> ```
> To find the Mac users:
> ```bash
> cat bch709_student.txt | grep Mac
> ```
> To count the number of Mac users:
> ```bash
> cat bch709_student.txt | grep Mac | wc -l
> ```
> Alternatively:
> ```bash
> grep Mac bch709_student.txt | wc -l
> ```
> Or simply:
> ```bash
> grep -c Mac bch709_student.txt
> ```
> Using flags to filter lines that don't contain "Windows":
> ```bash
> grep -v Windows bch709_student.txt
> ```
> Combining multiple flags:
> ```bash
> grep -c -v Windows bch709_student.txt
> ```
> With case-insensitive search and colored output:
> ```bash
> grep --color -i mac bch709_student.txt
> ```
{: .checklist}

> ## How do I store the results in a new file?
> Use the `>` character for redirection:
> ```bash
> grep -i mac bch709_student.txt > mac_student
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
> Make sure you're in the working directory, then download:
> ```bash
> $ cd ~/bch709_data
> $ curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
> $ gunzip mrna.fa.gz
> $ ls -lh mrna.fa
> ```
> ```output
> -rw-r--r-- 1 user group 416M Jan 20 10:00 mrna.fa
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

> ## What if the Tool Doesn't Have a Manual Page?
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

### Split fasta
Split a multi-sequence FASTA file into individual files (one sequence per file):
```bash
# Using mrna.fa downloaded earlier
awk '/^>/{f=++d".fasta"} {print > f}' mrna.fa
ls *.fasta | head
```

### Merge fasta
Combine multiple FASTA files into one:
```bash
cat 1.fasta 2.fasta 3.fasta >> myfasta.fasta
# Or use wildcards:
cat ?.fasta > single_digit.fasta   # matches 1.fasta, 2.fasta, etc.
cat ??.fasta > double_digit.fasta  # matches 10.fasta, 11.fasta, etc.
```

### Search fasta
Search for specific sequences or patterns:
```bash
# Find exact sequence
grep -n --color "GAATTC" mrna.fa | head
# Find pattern with regex (? = 0 or 1 of preceding char)
grep -n --color -E 'GAA?TTC' mrna.fa | head
```

### Regular Expression
A regular expression is a pattern that the regular expression engine attempts to match in input text. A pattern consists of one or more character literals, operators, or constructs.
Please play [this](https://regexone.com/lesson)


### GFF file
The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines. The following documentation is based on the Version 3 [(http://gmod.org/wiki/GFF3)](http://gmod.org/wiki/GFF3) specifications.

Download the GFF file:
```bash
$ cd ~/bch709_data
$ curl -L -o MGI.gff3.gz http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz
$ ls MGI.gff3.gz
```
```output
MGI.gff3.gz
```
What is `.gz` ?
```bash
$ file MGI.gff3.gz
```
```output
MGI.gff3.gz: gzip compressed data
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

Uncompress the file, examine the size difference, then keep it uncompressed for later exercises:

```bash
$ cd ~/bch709_data
$ ls -lh MGI.gff3.gz
```
```output
-rw-r--r-- 1 user group 12M Jan 20 10:00 MGI.gff3.gz
```
```bash
$ gunzip MGI.gff3.gz
$ ls -lh MGI.gff3
```
```output
-rw-r--r-- 1 user group 95M Jan 20 10:00 MGI.gff3
```
Notice the uncompressed file is ~8x larger! You can recompress with `gzip MGI.gff3` if needed.
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

Print all CDS lines:
```bash
$ cat MGI.gff3 | cut -f 3 | grep CDS | head
```

Print CDS and ID (step by step):
```bash
# Step 1: Select relevant columns
$ cat MGI.gff3 | cut -f 1,3,4,5,7,9 | head

# Step 2: Filter for CDS only
$ cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | head

# Step 3: Remove everything after semicolon
$ cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep CDS | sed 's/;.*//g' | head

# Step 4: Remove ID= prefix
$ cat MGI.gff3 | cut -f 1,3,4,5,7,9 | grep $'\tCDS\t' | sed 's/;.*//g' | sed 's/ID=//g' | head
```
Print length of each gene in a GFF3 file.
```bash
grep $'\tgene\t' MGI.gff3 | cut -s -f 4,5 | perl -ne '@v = split(/\t/); printf("%d\n", $v[1] - $v[0] + 1)'
```

Extract all gene IDs from a GFF3 file.
```bash
grep $'\tgene\t' MGI.gff3 | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'
```

## GFF3 file format
- Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
- seqid - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an
- source - name of the program that generated this feature, or the data source (database or project name)
- type - type of feature. Must be a term or accession from the SOFA sequence ontology
- start - Start position of the feature, with sequence numbering starting at 1.
- end - End position of the feature, with sequence numbering starting at 1.
- score - A floating point value.
- strand - defined as + (forward) or - (reverse).
- phase - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
attributes - A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation for [more details](http://gmod.org/wiki/GFF3).

Returns all lines on Chr 10 between 5.2MB and 5.45MB in  MGI.gff3. (assumes) chromosome in column 1 and position in column 4:
```bash
cat MGI.gff3 | awk '$1=="10"' | awk '$4>=5200000' | awk '$4<=5450000'

cat MGI.gff3 | awk '$1=="10"' | awk '$4>=5200000' | awk '$4<=5450000' |  grep mRNA

cat MGI.gff3 | awk '$1=="10"' | awk '$4>=5200000' | awk '$4<=5450000' |  grep mRNA | awk '{print $0,$5-$4}'
```
Returns specific lines
```bash
sed -n '1,10p' MGI.gff3

sed -n '52p' MGI.gff3
```

Time and again we are surprised by just how many applications it has, and how frequently problems can be solved by sorting, collapsing identical values, then resorting by the collapsed counts. The skill of using Unix is not just that of understanding the commands themselves. It is more about recognizing when a pattern, such as the one that we show above, is the solution to the problem that you wish to solve. The easiest way to learn to apply these patterns is by looking at how others solve problems, then adapting it to your needs.


## Quick reminder
- Learn basic Bash. Actually, type `man bash` and at least skim the whole thing; it's pretty easy to follow and not that long. Alternate shells can be nice, but Bash is powerful and always available (learning *only* zsh, fish, etc., while tempting on your own laptop, restricts you in many situations, such as using existing servers).

- Learn at least one text-based editor well. The `nano` editor is one of the simplest for basic editing (opening, editing, saving, searching). However, for the power user in a text terminal, there is no substitute for Vim (`vi`), the hard-to-learn but venerable, fast, and full-featured editor. Many people also use the classic Emacs, particularly for larger editing tasks. (Of course, any modern software developer working on an extensive project is unlikely to use only a pure text-based editor and should also be familiar with modern graphical IDEs and tools.)

- Finding documentation:
  - Know how to read official documentation with `man` (for the inquisitive, `man man` lists the section numbers, e.g. 1 is "regular" commands, 5 is files/conventions, and 8 are for administration). Find man pages with `apropos`.
  - Know that some commands are not executables, but Bash builtins, and that you can get help on them with `help` and `help -d`. You can find out whether a command is an executable, shell builtin or an alias by using `type command`.
  - `curl cheat.sh/command` will give a brief "cheat sheet" with common examples of how to use a shell command.
- Learn about redirection of output and input using `>` and `<` and pipes using `|`. Know `>` overwrites the output file and `>>` appends. Learn about stdout and stderr.

- Basic file management: `ls` and `ls -l` (in particular, learn what every column in `ls -l` means), `less`, `head`, `tail` and `tail -f` (or even better, `less +F`), `ln` and `ln -s` (learn the differences and advantages of hard versus soft links), `chown`, `chmod`, `du` (for a quick summary of disk usage: `du -hs *`). For filesystem management, `df`, `mount`, `fdisk`, `mkfs`, `lsblk`. Learn what an inode is (`ls -i` or `df -i`).

- Know regular expressions well, and the various flags to `grep`/`egrep`. The `-i`, `-o`, `-v`, `-A`, `-B`, and `-C` options are worth knowing.

## Everyday use
- In Bash, use **Tab** to complete arguments or list all available commands and **ctrl-r** to search through command history (after pressing, type to search, press **ctrl-r** repeatedly to cycle through more matches, press **Enter** to execute the found command, or hit the right arrow to put the result in the current line to allow editing).

- Use `alias` to create shortcuts for commonly used commands. For example, `alias ll='ls -latr'` creates a new alias `ll`.

- Save aliases, shell settings, and functions you commonly use in `~/.bashrc`, and [arrange for login shells to source it](http://superuser.com/a/183980/7106). This will make your setup available in all your shell sessions.

- To see recent commands, use `history`. Follow with `!n` (where `n` is the command number) to execute again. There are also many abbreviations you can use, the most useful probably being `!$` for last argument and `!!` for last command (see "HISTORY EXPANSION" in the man page). However, these are often easily replaced with **ctrl-r** and **alt-.**.


## Obscure but useful

- `expr`: perform arithmetic or boolean operations or evaluate regular expressions

- `m4`: simple macro processor

- `yes`: print a string a lot

- `cal`: nice calendar

- `env`: run a command (useful in scripts)

- `printenv`: print out environment variables (useful in debugging and scripts)

- `look`: find English words (or lines in a file) beginning with a string

- `cut`, `paste` and `join`: data manipulation

- `fmt`: format text paragraphs

- `pr`: format text into pages/columns

- `fold`: wrap lines of text

- `column`: format text fields into aligned, fixed-width columns or tables

- `expand` and `unexpand`: convert between tabs and spaces

- `nl`: add line numbers

- `seq`: print numbers

- `bc`: calculator

- `factor`: factor integers

- [`gpg`](https://gnupg.org/): encrypt and sign files

- `toe`: table of terminfo entries

- `nc`: network debugging and data transfer

- `socat`: socket relay and tcp port forwarder (similar to `netcat`)

- `dd`: moving data between files or devices

- `file`: identify type of a file

- `tree`: display directories and subdirectories as a nesting tree; like `ls` but recursive

- `stat`: file info

- `time`: execute and time a command

- `timeout`: execute a command for specified amount of time and stop the process when the specified amount of time completes.

- `lockfile`: create semaphore file that can only be removed by `rm -f`

- `logrotate`: rotate, compress and mail logs.

- `watch`: run a command repeatedly, showing results and/or highlighting changes

- [`when-changed`](https://github.com/joh/when-changed): runs any command you specify whenever it sees file changed. See `inotifywait` and `entr` as well.

- `tac`: print files in reverse

- `comm`: compare sorted files line by line

- `strings`: extract text from binary files

- `tr`: character translation or manipulation

- `iconv` or `uconv`: conversion for text encodings

- `split` and `csplit`: splitting files

- `sponge`: read all input before writing it, useful for reading from then writing to the same file, e.g., `grep -v something some-file | sponge some-file`

- `units`: unit conversions and calculations; converts furlongs per fortnight to twips per blink (see also `/usr/share/units/definitions.units`)

- `apg`: generates random passwords

- `xz`: high-ratio file compression

- `ldd`: dynamic library info

- `nm`: symbols from object files

- `ab` or [`wrk`](https://github.com/wg/wrk): benchmarking web servers

- `strace`: system call debugging

- [`mtr`](http://www.bitwizard.nl/mtr/): better traceroute for network debugging

- `cssh`: visual concurrent shell

- `rsync`: sync files and folders over SSH or in local file system

- [`wireshark`](https://wireshark.org/) and [`tshark`](https://www.wireshark.org/docs/wsug_html_chunked/AppToolstshark.html): packet capture and network debugging

- [`ngrep`](http://ngrep.sourceforge.net/): grep for the network layer

- `host` and `dig`: DNS lookups

- `lsof`: process file descriptor and socket info

- `dstat`: useful system stats

- [`glances`](https://github.com/nicolargo/glances): high level, multi-subsystem overview

- `iostat`: Disk usage stats

- `mpstat`: CPU usage stats

- `vmstat`: Memory usage stats

- `htop`: improved version of top

- `last`: login history

- `w`: who's logged on

- `id`: user/group identity info

- [`sar`](http://sebastien.godard.pagesperso-orange.fr/): historic system stats

- [`iftop`](http://www.ex-parrot.com/~pdw/iftop/) or [`nethogs`](https://github.com/raboof/nethogs): network utilization by socket or process

- `ss`: socket statistics

- `dmesg`: boot and system error messages

- `sysctl`: view and configure Linux kernel parameters at run time

- `hdparm`: SATA/ATA disk manipulation/performance

- `lsblk`: list block devices: a tree view of your disks and disk partitions

- `lshw`, `lscpu`, `lspci`, `lsusb`, `dmidecode`: hardware information, including CPU, BIOS, RAID, graphics, devices, etc.

- `lsmod` and `modinfo`: List and show details of kernel modules.

- `fortune`, `ddate`, and `sl`: um, well, it depends on whether you consider steam locomotives and Zippy quotations "useful"


## Recent unix command

{{< rawhtml >}}
<h1 align="center">Modern Unix</h1>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/sharkdp/bat"><code>bat</code></a>
  </h1>
  <p align="center">A <code>cat</code> clone with syntax highlighting and Git integration.</p>
  <p align="center">
    <img src="../fig/2lSW4RE.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/ogham/exa"><code>exa</code></a>
  </h1>
  <p align="center">A modern replacement for <code>ls</code>.</p>
  <p align="center">
    <img src="../fig/exa-screenshots.png" width="700" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/Peltoche/lsd"><code>lsd</code></a>
  </h1>
  <p align="center">The next gen file listing command. Backwards compatible with <code>ls</code>.</p>
  <p align="center">
    <img src="../fig/lsd-screen.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/dandavison/delta"><code>delta</code></a>
  </h1>
  <p align="center">A viewer for <code>git</code> and <code>diff</code> output</p>
  <p align="center">
    <img src="../fig/delta-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/bootandy/dust"><code>dust</code></a>
  </h1>
  <p align="center">A more intuitive version of <code>du</code> written in rust.</p>
  <p align="center">
    <img src="../fig/dust-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/muesli/duf"><code>duf</code></a>
  </h1>
  <p align="center">A better <code>df</code> alternative </p>
  <p align="center">
    <img src="../fig/duf-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/Canop/broot"><code>broot</code></a>
  </h1>
  <p align="center">A new way to see and navigate directory <code>tree</code>s</p>
  <p align="center">
    <img src="../fig/broot-overview.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/sharkdp/fd"><code>fd</code></a>
  </h1>
  <p align="center">A simple, fast and user-friendly alternative to <code>find</code>.</p>
  <p align="center">
    <img src="../fig/fd-screencast.svg" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/BurntSushi/ripgrep"><code>ripgrep</code></a>
  </h1>
  <p align="center">An extremely fast alternative to <code>grep</code> that respects your gitignore</p>
  <p align="center">
    <img src="../fig/ripgrep-demo.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/ggreer/the_silver_searcher"><code>ag</code></a>
  </h1>
  <p align="center">A code searching tool similar to <code>ack</code>, but faster.</p>
  <p align="center">
    <img src="../fig/ag-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/junegunn/fzf"><code>fzf</code></a>
  </h1>
  <p align="center">A general purpose command-line fuzzy finder.</p>
  <p align="center">
    <img src="../fig/fzf-preview.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/cantino/mcfly"><code>mcfly</code></a>
  </h1>
  <p align="center">Fly through your shell <code>history</code>. Great Scott! </p>
  <p align="center">
    <img src="../fig/mcfly-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/theryangeary/choose"><code>choose</code></a>
  </h1>
  <p align="center"> A human-friendly and fast alternative to <code>cut</code> and (sometimes) <code>awk</code> </p>
  <p align="center">
    <img src="../fig/choose-asciinema.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/stedolan/jq"><code>jq</code></a>
  </h1>
  <p align="center">
    <code>sed</code> for JSON data.
  </p>
  <p align="center">
    <img src="../fig/jq-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/chmln/sd"><code>sd</code></a>
  </h1>
  <p align="center">An intuitive find & replace CLI (<code>sed</code> alternative).</p>
  <p align="center">
    <img src="../fig/sd-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/cheat/cheat"><code>cheat</code></a>
  </h1>
  <p align="center">Create and view interactive cheatsheets on the command-line.</p>
  <p align="center">
    <img src="../fig/cheat-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/tldr-pages/tldr"><code>tldr</code></a>
  </h1>
  <p align="center">A community effort to simplify <code>man</code> pages with practical examples.</p>
  <p align="center">
    <img src="../fig/tldr-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/ClementTsang/bottom"><code>bottom</code></a>
  </h1>
  <p align="center">Yet another cross-platform graphical process/system monitor.</p>
  <p align="center">
    <img src="../fig/bottom-demo.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/nicolargo/glances"><code>glances</code></a>
  </h1>
  <p align="center">Glances an Eye on your system. A <code>top</code>/<code>htop</code> alternative for GNU/Linux, BSD, Mac OS and Windows operating systems.</p>
  <p align="center">
    <img src="../fig/glances-summary.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/aksakalli/gtop"><code>gtop</code></a>
  </h1>
  <p align="center">System monitoring dashboard for terminal.</p>
  <p align="center">
    <img src="../fig/gtop-demo.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/sharkdp/hyperfine"><code>hyperfine</code></a>
  </h1>
  <p align="center">A command-line benchmarking tool.</p>
  <p align="center">
    <img src="../fig/z19OYxE.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/orf/gping"><code>gping</code></a>
  </h1>
  <p align="center"><code>ping</code>, but with a graph.</p>
  <p align="center">
    <img src="../fig/gping-demo.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/dalance/procs"><code>procs</code></a>
  </h1>
  <p align="center">A modern replacement for <code>ps</code> written in Rust.</p>
  <p align="center">
    <img src="../fig/procs-screenshot.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/httpie/httpie"><code>httpie</code></a>
  </h1>
  <p align="center">A modern, user-friendly command-line HTTP client for the API era.</p>
  <p align="center">
    <img src="../fig/httpie-demo.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/rs/curlie"><code>curlie</code></a>
  </h1>
  <p align="center">The power of <code>curl</code>, the ease of use of <code>httpie</code>.</p>
  <p align="center">
    <img src="../fig/curlie-get.png" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/ducaale/xh"><code>xh</code></a>
  </h1>
  <p align="center">A friendly and fast tool for sending HTTP requests. It reimplements as much as possible of HTTPie's excellent design, with a focus on improved performance.</p>
  <p align="center">
    <img src="../fig/xh-demo.gif" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/ajeetdsouza/zoxide"><code>zoxide</code></a>
  </h1>
  <p align="center">A smarter <code>cd</code> command inspired by <code>z</code>.</p>
  <p align="center">
    <img src="../fig/zoxide-tutorial.webp" width="600" />
  </p>
</p>

<p align="center">
  <h1 align="center">
    <a href="https://github.com/ogham/dog"><code>dog</code></a>
  </h1>
  <p align="center">A user-friendly command-line DNS client. <code>dig</code> on steroids</p>
  <p align="center">
    <img src="../fig/dog-screenshot.png" width="700" />
  </p>
</p>
{{< /rawhtml >}}


## Advanced Linux Topics

The following topics are covered in separate pages for more detailed learning:

| Topic | Description |
|-------|-------------|
| [File Permissions](Linux_File_Permissions) | Understanding and managing file permissions with chmod |
| [Environment Variables](Linux_Environment_Variables) | Working with system and user environment variables |
| [Process Management](Linux_Process_Management) | Managing processes, background jobs, and signals |
| [Symbolic Links](Linux_Symbolic_Links) | Creating and understanding symbolic and hard links |
| [Find Command](Linux_Find_Command) | Searching for files by name, type, size, and more |
| [sed - Stream Editor](Linux_sed_Stream_Editor) | Text transformations and find-and-replace operations |
| [awk - Pattern Processing](Linux_awk_Pattern_Processing) | Processing structured text data and columns |
| [Shell Scripting Basics](Linux_Shell_Scripting_Basics) | Automating tasks with bash scripts |
| [Standard Streams and Redirection](Linux_Standard_Streams_Redirection) | Working with stdin, stdout, stderr, and pipes |
| [Bioinformatics File Formats](Bioinformatics_File_Formats) | SAM/BAM, BED, VCF formats and tools |


## Screen and tmux - Terminal Multiplexers

Run long jobs that continue after you disconnect from the server.

> ## tmux Quick Start
> ```bash
> # Start a new named session
> tmux new -s analysis
>
> # Run your long command
> # ... your command here ...
>
> # Detach: Press Ctrl+B, then D
>
> # List sessions
> tmux ls
> ```
> ```output
> analysis: 1 windows (created Mon Jan 20 10:00:00 2026)
> ```
> ```bash
> # Reattach later
> tmux attach -t analysis
>
> # Kill session when done
> tmux kill-session -t analysis
> ```
{: .keypoints}

### tmux Key Bindings (Ctrl+B, then)

| Key | Action |
|-----|--------|
| `d` | Detach from session |
| `c` | Create new window |
| `n` | Next window |
| `p` | Previous window |
| `%` | Split vertically |
| `"` | Split horizontally |


## Git Basics - Version Control

Track changes to your code and collaborate with others.

> ## Step 1: Configure Git
> ```bash
> # Set your identity (one-time setup)
> git config --global user.name "Your Name"
> git config --global user.email "your.email@example.com"
>
> # Verify
> git config --list | grep user
> ```
> ```output
> user.name=Your Name
> user.email=your.email@example.com
> ```
{: .keypoints}

> ## Step 2: Create a Repository
> ```bash
> # Create and initialize a project
> mkdir my_analysis
> cd my_analysis
> git init
> ```
> ```output
> Initialized empty Git repository in /home/user/my_analysis/.git/
> ```
> ```bash
> # Create some files
> echo "# My Analysis" > README.md
> echo "data/" > .gitignore
>
> # Check status
> git status
> ```
> ```output
> On branch main
> Untracked files:
>   README.md
>   .gitignore
> ```
{: .keypoints}

> ## Step 3: Make Your First Commit
> ```bash
> # Add files to staging
> git add README.md .gitignore
>
> # Commit with a message
> git commit -m "Initial commit: add README and gitignore"
> ```
> ```output
> [main (root-commit) abc1234] Initial commit: add README and gitignore
>  2 files changed, 2 insertions(+)
> ```
> ```bash
> # View history
> git log --oneline
> ```
> ```output
> abc1234 Initial commit: add README and gitignore
> ```
{: .keypoints}

### Git Commands Reference

| Command | Description |
|---------|-------------|
| `git init` | Initialize repository |
| `git status` | Check file status |
| `git add FILE` | Stage changes |
| `git commit -m "msg"` | Commit changes |
| `git log` | View history |
| `git diff` | See changes |
| `git clone URL` | Clone repository |
