---
layout: page
title: Linux Enviroment and command line (II)
published: true
---

{% include gh_variables.html %}

## Command Line Bootcamp
[Command-line bootcamp](http://35.199.189.11:8081/), a tutorial that teaches you how to work at the command-line. You'll learn all the basic skills needed to start being productive in the UNIX terminal.


## Basic commands
|Navigation| cd, ls, pwd|
|File creation|touch,nano,mkdir,cp,mv,rm,rmdir|
|Reading|more,less,head,tail,cat|
|Compression|zip,gzip,bzip2,tar,compress|
|Uncompression|unzip,gunzip,bunzip2,uncompress|
|Permissions|chmod|
|Help|man|

## pwd

```bash
pwd
```
```output
/home/username
```
:::info
Returns you the `p`resent `w`orking `d`irectory (print working directory)

This means, you are now working in the home directory, which is located under your ID. The directory that will be the first directory when you log-in. You can also avoid writing the full path by using ~ in front of your username or simply ~.
:::
