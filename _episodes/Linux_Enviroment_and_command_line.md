---
layout: lesson
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

![Lnix_family tree](https://aerojsoft.files.wordpress.com/2016/02/linus-distribution-family-tree.jpg)


## Operating Systems Tasks
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


**The most common versions are Solaris, Linux, and MacOS X. **





Linux is a Unix-like computer operating system assembled under the model of free and open source software development and distribution. The defining component of Linux is the Linux kernel, an operating system kernel first released 5 October 1991 by Linus Torvalds. Linux was originally developed as a free operating system for Intel x86-based personal computers. It has since been ported to more computer hardware platforms than any other operating system. It is a leading operating system on servers and other big iron systems such as mainframe computers and supercomputers:more than 90% of today's 500 fastest supercomputers run some variant of Linux,including the 10 fastest. Linux also runs on embedded systems (devices where the operating system is typically built into the firmware and highly tailored to the system) such as mobile phones, tablet computers, network routers, televisions and video game consoles; the Android system in wide use on mobile devices is built on the Linux kernel.


https://www.usna.edu/Users/cs/bilzor/ic221/calendar.php?type=unit&event=1

https://www.usna.edu/Users/cs/bilzor/ic221/calendar.php?type=unit&event=1

![Unix_family tree]({{site.baseurl}}/fig/unix-simple.png)


https://medium.com/@youngstone89/unix-introduction-shell-980212852897





![OS]({{site.baseurl}}/fig/OS.png)
So, what is UNIX?


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



>UNIX is an OS that was first developed in the 1960s, and has been developed constantly ever since then. It is a collections of programs to make the computer work properly with stable, multi-user, multi-tasking system for servers, desktops, and laptops.
{: .callout}


