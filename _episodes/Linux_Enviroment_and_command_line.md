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
> Please finish below assignment before due date  
> - [DataCamp Introduction to Shell](https://learn.datacamp.com/courses/introduction-to-shell) Due by 9/21/22  
> - [DataCamp Downloading Data on the Command Line](https://campus.datacamp.com/courses/data-processing-in-shell/downloading-data-on-the-command-line?ex=1) Due by 9/19/22 
> - [File permission and chmod](https://www.youtube.com/watch?v=3gcSeDoQ_rU) Due by 9/19/22 
> - [Nano Text Editor Basics](https://www.youtube.com/watch?v=Jf0ZJZJ8jlI) Due by 9/19/22 
{: .prereq}

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
{: .checklist}

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
{: .checklist}

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
{: .checklist}


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
>```bash
>zip OUTFILE.zip INFILE.txt Compress INFILE.txt  
>zip -r OUTDIR.zip DIRECTORY Compress all files in a DIRECTORY into one archive file (OUTDIR.zip)  
>zip -r OUTFILE.zip . -i \*.txt Compress all txt files in a DIRECTORY into one archive file (OUTFILE.zip)   
>unzip SOMEFILE.zip
>```
{: .checklist}

>## TAR compression/extraction
>```bash
>tar (tape archive) utility saves many files together into a single archive file, and restores individual files from the archive. It also includes automatic archive compression/decompression options and special features for incremental and full backups.
>tar -xzvf SOMEFILE.tar.gz # extract contents of SOMEFILE.tar
>tar -czvf OUTFILE.tar.gz DIRECTORY #extract contents of gzipped archive SOMEFILE.tar.gz
>tar -czvf OUTFILE.tar.gz \*.txt #archive and compress all files in a directory into one archive file
>tar -czvf backup.tar.gz BACKUP_WORKSHOP #archive and compress all ".txt" files in current directory into one archive file
>```
{: .checklist}

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


## Cloud platform choices

There are several cloud providers to choose from. Some scientific clouds may either be free or allocate resources competitively. Commercial clouds can be very powerful, but choice can be overwhelming. Availablity of public and commercial cloud resources also vary by country and region.

The major tradeoff between platforms is between flexibility and cost. Generally speaking,
services that allow you more flexibility and autonomy, will be more expensive than
highly managed services.

Below, we have highlighted the three types of computing resources: Clusters, Open Science Clouds, and Commercial Clouds, which are widely available to researchers around the world. However, the availability of any specific cluster or cloud may be region specific, and this is not meant to be an exhaustive list. We encourage researchers to use this list as a starting point for learning about cloud resources and suggest checking with your local or regional government to see what other Open Science Clouds might be available to you. The cloud resources listed here should be available to any scientist based in the US, but may be unavailable, or have different pricing, in other countries. 

### University/Corporate Computing Clusters

Many universities and businesses operate their own computing clusters that are available to students and staff
at low or no cost. If your employer maintains a computing cluster, this will almost always be the least
expensive option.

However, most HPCCs (High Performance Computing Clusters) put limits on:
 - The number of processors a user can utilize at once
 - The amount of disk storage per user
 - The amount of time a single process can run
 - What programs can be installed, and by whom
 - Who can have accounts and access data

 HPCCs are also a shared resource, so even when you have access, your programs are unlikely to run
 immediately. Most HPCCs run some kind of scheduler, that you submit your processing jobs to, and
 it runs them as resources become available. In order to submit a job, you generally will need to
 know not only what program you want to run, but with how many processors and for how long.
 While interacting with the scheduler is no more difficult than interacting with the shell,
 it will have its own set of commands and syntax that you'll need to learn; and these vary
 widely among HPCCs.

 There are also many upsides to using an HPCC. As previously mentioned, they're generally the least
 expensive option, but they often come with more perks. For instance, many HPCCs offer free or low-cost
 training, storage space, back-up options, and technical support. If your HPCC has a scheduler, you
 can also queue up many sequential jobs all at once, and you don't have to worry about racking up
 fees on instances that are sitting idle. It's often also much easier to pay for HPCC use than to
 pay for Amazon using grant money, however universities are getting better about AWS payments.

### Open Science Clouds

#### [XSEDE](https://www.xsede.org/)

The Extreme Science and Engineering Discovery Environment (XSEDE) is an NSF funded HPCC, so
it is open to any US-based researcher, and shares most of the same benefits and drawbacks
of a university or corporate HPCC. If your university or corporation doesn't have it's
own HPCC resources, XSEDE will likely be your cheapest option.

Although any US-based researcher can use XSEDE, first [they'll need an account](https://portal.xsede.org/#/guest).
Like the HPCC options described above, XSEDE uses a scheduler to start jobs, and puts limits on
how many resources any one user can utilize at once.

XSEDE can also be a bit intimidating at first because you will need to know what resources
you need, and for how long, before you get started. XSEDE runs like a mini version of the
NSF grant system. In order to qualify to submit large jobs, you'll have to submit a [allocation request](https://portal.xsede.org/allocations/research), in the form of a short proposal.
Also like an NSF grant, if your proposal is accepted, that means you have access to whatever
resources you were approved for, for the time frame you requested.

Don't let that paragraph scare you off though. XSEDE has two different allocation tracks. If
you aren't sure exactly what you'll need for your big project, you can request a [startup allocation](https://portal.xsede.org/allocations/startup) which only requires an abstract
rather than a proposal, and grants you a year to try out your new pipeline or analysis. These
are usually granted in a week or so, and are intended for you to test your pipeline so you
know what to ask for in your allocation proposal.

If that still sounds a little too daunting, XSEDE also has [trial allocations](https://iujetstream.atlassian.net/wiki/spaces/JWT/pages/76149919/Jetstream+Trial+Access+Allocation)
which give you access to only a tiny fraction of XSEDES power, but are plenty large enough to
test your code and see if a larger allocation is worth pursuing. These allocations are granted
more or less immediately by simply filling in a form and agreeing to the usage rules.

If you're interested in using XSEDE, check to see if your workplace has a [Campus Champion](https://www.xsede.org/community-engagement/campus-champions). These are people who
have had extensive training on both the XSEDE system and the allocation program, and can
help you figure out how to apply and what you need.

#### [Open Science Grid](https://opensciencegrid.org)

The Open Science Grid (OSG) is an NSF-funded national network of computing centers that
have pooled their resources together and made them available to various research groups.
The OSG is usable by any researcher based at a US institution, and is accessible
for free without an allocation. It can provide millions of computing hours for researchers
who have problems that fit well on its setup.

Certain projects and universities have direct access to the Open Science Grid, but any
researcher can access it through the [OSG Connect](https://osgconnect.net/) entry point.
If you apply for OSG access through that website, you will have a consultation with
someone who can help you determine if your analysis is a good fit for and how to get started.

The OSG is a great fit for problems that can be broken into lots of independent pieces.
One good example is read alignment: the starting read data can be broken into several
pieces, each of them aligned, and then the results combined. Another good problem type
for the OSG are multi-start simulations or statistical analyses where you need to run the
same model or simulation many, many times. The payoff of using this
approach is being able to run on many hundreds (sometimes thousands!) of computers at
once, accelerating your analysis.

Note that you don't access a specific computing center through OSG -- unlike XSEDE, where
you apply for time and then run on a specific HPCC resource, the OSG sits on top of many
resources and when you submit your work, it could run almost anywhere in the overall system.

#### [Open Science Data Cloud (OSDC)](https://www.opensciencedatacloud.org/)

The Open Science Data Cloud provides the scientific community with resources for storing, sharing, and analyzing terabyte and petabyte-scale scientific datasets. OSDC's Bionimbus Protected Data Cloud (PDC) is a platform designed with the sole purpose of analysing and sharing protected genomics data.

#### [Atmosphere](https://pods.iplantcollaborative.org/wiki/display/atmman/Getting+Started)


#### [CyVerse (iPlant Collaborative) Atmosphere](http://www.cyverse.org/atmosphere)


#### [JetStream](http://jetstream-cloud.org/)

### Commercial Clouds
Computing architecture is moving (albeit at a slow pace) to the Model-to-Data paradigm. This means that scientists should be encouraged to bring their compute to where the data is stored, instead of the the other way around. The following outlines the general differences between the three major commercial cloud providers: Amazon Web Services (AWS), Google Cloud Platform (GCP) and Microsoft Azure.

Essentially all cloud providers provide extremely similar computing and storage options; you can "rent" or provision computing infrastructure with very similar specifications across all three cloud vendors. Even the costs are highly comparable. What governs how to choose the right cloud computing vendor is highly opportunistic: (1)funding options, (2)solidarity with collaborating/similar scientific groups, (3)location of datasets that a particular research group works with and (4)familiarity with cloud vendor services.

1. Funding options: Does your grant stipulate where you should build your computing pipeline? For example, the NIH often partners with specific cloud vendors to provide cloud credits that allow researchers to compute for free. Some cloud vendors also provide research credits.
2. Solidarity with collaborating/similar scientific groups: Are other research groups in your field drawn to a specific cloud vendor? It might make sense to utilize the same cloud service to minimize transfer (egress) costs especially if you are sharing large datasets. You may also be able to make use of existing pipelines without reinventing the wheel if you choose the same cloud provider that your collaborators are using.
3. Location of datasets that a particular research group works with: Again, thinking of bringing your models to the where the data is stored helps minimize costs and saves you time in having to download and store data separately.
4. Services here refer to cloud vendor add-ons that take away the need for a user to set up their own computing infrastructure.  A fully managed database (e.g. AWS RDS, GCP CloudSQL, Azure SQL DB) is an example of a service. If you are used to SQL Server, you may want to look into options provided by Azure. Are you more familiar with Postgres SQL? Then AWS and GCP might provide cheaper options for you.

#### [Amazon EC2](http://aws.amazon.com/ec2/)

The Amazon Web Service (AWS) that you've been using is the Elastic Compute (EC2) cloud. There
are actually lots of other cloud and storage solutions under the AWS umbrella, but when most
data scientists say AWS, they mean [EC2](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html).
With EC2, you can rent access to a cloud computing resource as small as your laptop, or as large as a 64 processor
machine with 488GB of memory, and with a number of different operating systems. These instances can
be optimized for jobs that are memory intensive, or require a lot of bandwidth, or [almost any other
specific need](https://aws.amazon.com/ec2/instance-types/). There are so many options that we can't
cover them all here, but these are a few popular ones:

##### On-Demand
All this variety and optimization makes EC2 much more expensive than an average HPCC, however,
depending on your needs it can be [quite affordable](https://aws.amazon.com/ec2/pricing/). If
you want to start an EC2 instance whenever you want and have instant access, you can
rent a quite large on-demand instance with 8 processors and 32 GB of memory for ~40 cents
an hour, and tiny instances are only about half a cent per hour.


##### Spot-Instances
If your program can tolerate pauses, and you don't need the analysis done as fast as
possible, you can request a spot instance. Essentially, whenever Amazon has computing
capacity that no one is paying them for, they lower the prices for renting some systems.
If you request a spot-instance, that means you specify the size and parameters of
your machine rental, and set a limit on how much you're willing to spend per hour. Then,
whenever the rental rate dips below your maximum, your instance turns on and runs until
the price goes back up. If it's not an important shopping season, and you aren't in a hurry,
you can often run spot-instances for less than half their normal cost.

##### Free Tier
There are also [free options](https://aws.amazon.com/free/), which allow you to test out
the interface and make sure it will meet your needs before you start renting.

Just remember that with EC2 and all other commercial services, you're paying for renting the computer,
whether you're using it or not. If you leave an instance on idly for months after your pipeline has finished,
you'll still have to pay for that time.

#### [Google Cloud](https://cloud.google.com/): [getting started](https://cloud.google.com/compute/docs/quickstart)
GCP offers very competitive prices for compute and storage (as of July 2019, their compute pricing is lower than that of AWS and Azure for instances of comparable specifications). If you are looking to dabble in cloud computing but do not need a vast catalog of services, GCP would be a good place to start looking.

Their version of "Spot Intances" are known as pre-emptible instances and offer very competitive pricing. GCP also has TPUs.

#### [Microsoft Azure](https://azure.microsoft.com/en-us/)
If your software requires Microsoft Windows, it may be cheaper to use MS Azure due to licensing issues. Azure's computing instances are known as Azure Virtual Machines and often come at a slightly higher cost than other cloud computing vendors' offerings. If a lot of your computing pipeling is Windows dependent, it may make sense to build everything on MS Azure from the get go.

#### [IBM Cloud](https://www.ibm.com/cloud)
IBM Cloud offers more than 11 million bare metal configurations in virtual mode which are customizable RAM and SSDs on bare metal. They also have an on-demand provisioning for all servers whose management and monitoring included along with the direct and cost-free tech support

## How to Choose

As you can see, highly managed systems (HPCCs, XSEDE, etc) usually are free or cheap, but
relatively inflexible. There may be certain programs you can't install, or there may be long
wait times. Commercial systems are generally more flexible because you can make them look
however you want, but they can be quite expensive, especially if you run for a long time, or have a
lot of data. However, there are other things to consider.

Another way to think about this is not *whether* you want to spend your time and money, but *where* you
want to spend them. AWS will let you install anything you want, but that also means that you'll
have to spend some time up-front installing programs and testing configurations. Your HPCC jobs
might take a week to start, but you don't have to do any of the systems administration, and you
can work on other projects while your job sits in the queue.

Your familiarity with the pipeline can also be a factor. Let's say you want to run a program
you've never run before. If you have no way to estimate how
long your job will take, a service like AWS might save you money, because you can run your
program until its done, no matter how long that is. Running it on an HPCC can be really
frustrating, because you'll need to submit your job with the amount of time and resources it will use.
An HPCC will only let your job run the amount of time you requested, so if you underestimate the
amount of time it will take, even by one minute, the system kills your program and you need to resubmit it.
On the other hand, if you want to run that same program, but you can easily estimate the runtime, an
HPCC is going to be a much cheaper choice.

In my work, I often use both commercial and non-commercial services. I tend to use AWS for testing,
with small amounts of data, until I know how the program behaves. Then I port the pipeline to
my university HPCC for running the full dataset.

> ## Discussion
>
> In small groups or on your own, plot out your next bioinformatics project. With guidance
> from your instructors and the above references, try to determine not only what types of
> resources you'll need, but what platform will best suit your project.
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

> ## Human genomic data & Security
>
> Note that if you are working with human genomics data there might be ethical and legal
> considerations that affect your choice of cloud resources to use. The terms of use, and/or
> the legislation under which you are handling the genomic data, might impose heightened information
> security measures for the computing environment in which you intend to process it. This is a too broad
> topic to discuss in detail here, but in general terms you should think through the technical and
> procedural measures needed to ensure that the confidentiality and integrity of the human data you work
> with is not breached. If there are laws that govern these issues in the jurisdiction in which you work,
> be sure that the cloud service provider you use can certify that they support the necessary measures.
> Also note that there might exist restrictions for use of cloud service providers that operate in other
> jurisdictions than your own, either by how the data was consented by the research subjects or by the
> jurisdiction under which you operate. Do consult the legal office of your institution for guidance
> when processing human genomic data.
{: .callout}

### Other Resources:


Learn more about cloud computing in bioinformatics:

Fusaro VA, Patil P, Gafni E, Wall DP, Tonellato PJ (2011) **Biomedical Cloud Computing With Amazon Web Services**. PLoS Comput Biol 7(8): e1002147. doi: [10.1371/journal.pcbi.1002147](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002147)

Charlebois K, Palmour N, Knoppers BM (2016) **The Adoption of Cloud Computing in the Field of Genomics Research: The Influence of Ethical and Legal Issues**. PLoS ONE 11(10): e0164347. https://doi.org/10.1371/journal.pone.0164347

Langmead B, Nellore A (2018) **Cloud computing for genomic data analysis and collaboration** Nature Reviews Genetics 19 (208). doi: 10.1038/nrg.2017.113 (https://www.nature.com/articles/nrg.2017.113)





This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.

, a tutorial that teaches you how to work at the command-line. You'll learn all the basic skills needed to start being productive in the UNIX terminal.


This manual was adapted from [Linode](http://www.linode.com). Linode is the BEST knowledge site ever.
