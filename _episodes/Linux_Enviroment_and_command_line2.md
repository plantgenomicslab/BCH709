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
>$ touch file1
>$ cp file1 file2
>$ ls
>```
>
>What if we wanted to copy files from a different directory to our current directory? Let’s put a file in our home directory (specified by \~, remember) and copy it to the lecture directory.
>```bash
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
>https://pastebin.com/raw/9bCszuvD
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
>curl -L -o bch709_student.txt https://pastebin.com/raw/9bCszuvD
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
>          12            33            205
>```
>How many lines?
>```bash
>cat bch709_student.txt | wc -l
>```
>```output
>12
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
{: checklist}

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
>diff mac_student windows_student
>```

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
>cut -f 2 bch709_student.txt > os.txt
>```
{: checklist}