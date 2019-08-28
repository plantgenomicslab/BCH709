---
layout: page
title: Linux Enviroment and command line (II)
published: true
---

{% include gh_variables.html %}

## Your first Unix command.
```bash
echo "Hello World!"
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
>/home/username
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
>/home/username/bch709_test
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
>$ cd /home/<userid>/bch709_test/help
>```
> **relative path**
>```bash
>$ cd ../
>```
{: .checklist}

>## The way back home
>Again, your home directory is `/home/<userid>/`. . they take you back to your home directory (from wherever you were). You will frequently want to jump straight back to your home directory, and typing cd is a very quick way to get there.
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
>$ ls
>$ rmdir help
>```
>*You have to be outside a directory before you can remove it with rmdir*
{: .keypoints}

>## touch
>The following sections will deal with Unix commands that help us to work with files, i.e. copy files to/from places, move files, rename files, remove files, and most importantly, look at files. First, we need to have some files to play with. The Unix command touch will let us create a new, empty file.
```bash
$ touch --help
$ touch test.txt
$ touch exam.txt
$ touch ETA.txt
$ ls
```




