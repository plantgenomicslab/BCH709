---
layout: page
title: Linux Enviroment and command line (II)
published: true
---

{% include gh_variables.html %}

## Command Line Bootcamp
[Command-line bootcamp](http://35.199.189.11:8081/), a tutorial that teaches you how to work at the command-line. You'll learn all the basic skills needed to start being productive in the UNIX terminal.

>## Basic commands
>|Category|comnmand|
>|---|---|
>|Navigation| cd, ls, pwd|
>|File creation|touch,nano,mkdir,cp,mv,rm,rmdir|
>|Reading|more,less,head,tail,cat|
>|Compression|zip,gzip,bzip2,tar,compress|
>|Uncompression|unzip,gunzip,bunzip2,uncompress|
>|Permissions|chmod|
>|Help|man|
{: .prereq}

>## pwd
>Returns you the `p`resent `w`orking `d`irectory (print working directory)
>```bash
>$ pwd
>```
>```output
>/home/username
>```
> This means, you are now working in the home directory, which is located under your ID. The directory that will be the first directory when you log-in. You can also avoid writing the full path by using ~ in front of your username or simply ~.
{: .keypoints}

>## mkdir
>To create a directory, mkdir `make directory` can be used.
>$ mkdir <DIRECTORY>
>```bash
>$ mkdir bch709_test
>```
>```output
>
>```
> Unlike PC/Mac folders, here you can have space in different way (some special characters are okay). You can also specify the path where you want to create your new folder. 
{: .keypoints}

>## cd
>This is `c`hanges `d`irectory. 
>$ cd <DIRECTORY>
>```bash
>$ cd bch709_test
>```
>```output
>
>```
>```bash
>$ pwd
>```
>```output
>/home/username/bch709_test
>```
>Changes your present location to the parent directory
>```bash
>cd ..
>```
>Please check your current location with `pwd`.
>
>Changes your present location to the Grand parent directory
>```bash
>cd ../../
>```
>```bash
>cd ../../
>```
>Please check your current location with `pwd`.
{: .keypoints}

