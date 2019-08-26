---
layout: page
title: Linux Enviroment and command line
published: true
---

{% include gh_variables.html %}


![bioinformatics_DNA]({{site.baseurl}}/fig/DNA.jpg)

## History of UNIX 
The Unix operating system was conceived and implemented in 1969 at AT&T's Bell Laboratories in the United States by Ken Thompson, Dennis Ritchie, Douglas McIlroy, and Joe Ossanna. It was first released in 1971 and was initially entirely written in assembly language, a common practice at the time. Later, in a key pioneering approach in 1973, Unix was re-written in the programming language C by Dennis Ritchie (with exceptions to the kernel and I/O). The availability of an operating system written in a high-level language allowed easier portability to different computer platforms. With a legal glitch forcing AT&T to license the operating system's source code to anyone who asked. Unix quickly grew and became widely adopted by academic institutions and businesses. In 1984, AT&T divested itself of Bell Labs. Free of the legal glitch requiring free licensing, Bell Labs began selling Unix as a proprietary product.


## Why UNIX
- Unix is an important OS in the history of computing
- Two major OS variants, Unix-based and Windows-based
- Used in a lot of back-end systems and personal computing
- Unix derivatives like Linux are open source, and well known to the community and developed in the open where we can study and understand them.
- The skills you learn on Unix will easily translate to many other OS platforms because all Unix-based systems share standard characteristic

![Unix_family tree]({{site.baseurl}}/fig/unix-simple.png)

### The kernel
This is called as the hub of the operating system, serving as allocator of time and memory to programs and handling the filestore and communications in response to system calls.

### I/O
In computing, input/output or I/O (or, informally, io or IO) is the communication between an information processing system, such as a computer, and the outside world, possibly a human or another information processing system. Inputs are the signals or data received by the system and outputs are the signals or data sent from it. The term can also be used as part of an action; to "perform I/O" is to perform an input or output operation.

## Unix/Linux Main Components
### The Unix/Linux computer ecosystem can be divided into three main parts:
- **User Space**: Defines the applications, libraries, and standard utilities that are user accessible. When we write a program, it is from this perspective that we operate, without concern for the underlying components. For example, writing a "Hello World" program on any computer is the same from the user-perspective, but might be different when it comes to actually executing the program and writing "Hello World" to the terminal.
- **Kernel Space**: This refers to the operations of OS that manage the interface between user actions and the hardware. It is the central part of the OS, and its primary job is to pair user applications with the underlying hardware and allow multiple programs to share singular hardware components. For example, how does a user input event, such as typing 'a' on the keyboard, get translated into 'a' appearing on the screen? Or, how does two programs both read from disc at the same time or run on the CPU at the same time?
- **Hardware**: The underlying physical components of the computer. These include Input/Output devices, like keyboards and monitors, the CPU which does calculations, the memory components, and the network interface.



### Operating Systems Tasks
The operating system's primary task is to manage services as an interface between the user and the hardware. Examples include:
- File System: managing files on the user
- Device I/O: managing input from devices
- Processes: Starting, running, and stopping programs, and allowing multiple programs to run at once, i.e., program multiprogramming.
- Memory Management: Allocating runtime memory for process and separating memory between process and between user-space and the kernel-space

![OS]({{site.baseurl}}/fig/OS.png)


### As a user of the OS, you will see these interactions from two perspectives:
- Shell: You will use the shell to interact with the OS
- System Call API: You will program in C to interact with the OS
The big part of this interaction comes from the System Call API, which you will use the C programming language. Why C?
- C is a low level language
- The OS is written in C
- Understanding the OS and C together is a natural process and will make you a better programmer


**The most common versions are Solaris, Linux, and MacOS X.**

## LINUX  
Linux is a Unix-like computer operating system assembled under the model of free and open source software development and distribution. The defining component of Linux is the Linux kernel, an operating system kernel first released 5 October 1991 by Linus Torvalds. Linux was originally developed as a free operating system for Intel x86-based personal computers. It has since been ported to more computer hardware platforms than any other operating system. It is a leading operating system on servers and other big iron systems such as mainframe computers and supercomputers:more than 90% of today's 500 fastest supercomputers run some variant of Linux,including the 10 fastest. Linux also runs on embedded systems (devices where the operating system is typically built into the firmware and highly tailored to the system) such as mobile phones, tablet computers, network routers, televisions and video game consoles; the Android system in wide use on mobile devices is built on the Linux kernel.

![Lnix_family tree](https://aerojsoft.files.wordpress.com/2016/02/linus-distribution-family-tree.jpg)



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
![text editor](https://w.namu.la/s/faf332270510cbe6d8359870f39e006df2f67b7da39fa6bb994a8b43f86ca2c82c35c165e084aa611c18ed754995f1ac3c2db553e2a64cff5945b582bba3ed13dc812a4340d702808d45887eb414f905cab17ba7a7ec73026b41bb16dfe764f0)

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

#### Cheat sheet

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


## Command Line Bootcamp
[Command-line bootcamp](http://35.199.189.11:8081/), a tutorial that teaches you how to work at the command-line. You'll learn all the basic skills needed to start being productive in the UNIX terminal.







