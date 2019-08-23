---
layout: lesson
title: Linux Enviroment and command line
published: true
---

{% include gh_variables.html %}


![bioinformatics_DNA]({{site.baseurl}}/fig/DNA.jpg)

The Unix operating system was conceived and implemented in 1969 at AT&T's Bell Laboratories in the United States by Ken Thompson, Dennis Ritchie, Douglas McIlroy, and Joe Ossanna. It was first released in 1971 and was initially entirely written in assembly language, a common practice at the time. Later, in a key pioneering approach in 1973, Unix was re-written in the programming language C by Dennis Ritchie (with exceptions to the kernel and I/O). The availability of an operating system written in a high-level language allowed easier portability to different computer platforms. With a legal glitch forcing AT&T to license the operating system's source code to anyone who asked,[22] Unix quickly grew and became widely adopted by academic institutions and businesses. In 1984, AT&T divested itself of Bell Labs. Free of the legal glitch requiring free licensing, Bell Labs began selling Unix as a proprietary product.

Linux is a Unix-like computer operating system assembled under the model of free and open source software development and distribution. The defining component of Linux is the Linux kernel, an operating system kernel first released 5 October 1991 by Linus Torvalds. Linux was originally developed as a free operating system for Intel x86-based personal computers. It has since been ported to more computer hardware platforms than any other operating system. It is a leading operating system on servers and other big iron systems such as mainframe computers and supercomputers:more than 90% of today's 500 fastest supercomputers run some variant of Linux,including the 10 fastest. Linux also runs on embedded systems (devices where the operating system is typically built into the firmware and highly tailored to the system) such as mobile phones, tablet computers, network routers, televisions and video game consoles; the Android system in wide use on mobile devices is built on the Linux kernel.



https://www.usna.edu/Users/cs/bilzor/ic221/calendar.php?type=unit&event=1

https://www.usna.edu/Users/cs/bilzor/ic221/calendar.php?type=unit&event=1

![Unix_family tree]({{site.baseurl}}./fig/unix-simple.png)
![OS]({{site.baseurl}}./fig/OS.png)

https://medium.com/@youngstone89/unix-introduction-shell-980212852897

![Unix_family tree]({{site.baseurl}}/https://aerojsoft.files.wordpress.com/2016/02/linus-distribution-family-tree.jpg)
So, what is UNIX?
UNIX is an OS that was first developed in the 1960s, and has been developed constantly ever since then. It is a collections of programs to make the computer work properly with stable, multi-user, multi-tasking system for servers, desktops, and laptops.
Types of UNIX
The most common versions are Solaris, Linux, and MacOS X.
The UNIX operating system analysis
The kernel
This is called as the hub of the operating system, serving as allocator of time and memory to programs and handling the filestore and communications in response to system calls.
For an instance, in a case that a user runs command “rm myfile”. The shell searches for the file containing the program “rm”, and then requests the kernel via system calls to execute the program “rm” on “myfile”. When the process “rm myfile” gets done, the shell returns the UNIX prompt % to the user, indicating that it is waiting for another command.
The Shell
This serves as an interface between the user and the kernel. When a user logs in, the login program checks the username and password, then start another program called the shell. The shell is a command line interpreter. It interprets the commands the user run and arranges for them to be executed. The commands are the programs to be run. The shell itself has different types of shell with its own set of commands and functions.
Shell Prompt?
The prompt “$” is called the command prompt. While the prompt is being displayed, you can run a command.
Shell Types?
Bourne Shell — $ character is the default prompt.
Bourne Shell(sh)
Korn shell(ksh)
Bourne Again shell(bash)
POSIX shell(sh)
C Shell — % character is the default prompt.
C shell(csh)
TENEX/TOPS C shell(tcsh)
The original Unix shell was made in the mid 19’s by Stephen R. Bourne. Bourne shell was the first shell to show up in Unix world. Bourne shell is commonly installed as /bin/sh on most versions of Unix.
