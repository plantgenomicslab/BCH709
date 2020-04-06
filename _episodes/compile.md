---
layout: page
title: Compile and Software installation
published: true
---

{% include gh_variables.html %}


![software_compile]({{{site.baseurl}}/fig/software-compiler.png)
![software_compile](https://pbs.twimg.com/media/CYIT_SJWQAIExU8.png)

## macOS

> ## Install Xcode
>Inside the Terminal window, copy and paste (or type) the following command, and press the return key on your keyboard:
>```
>$ xcode-select --install
>```
{: .prereq}

> ## Install Homebrew
>```bash
>$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> ```
{: .prereq}

> ## Brew update
>```bash
>$ brew update
>```
{: .prereq}

### Install Prerequire Software
```bash
$ brew install openssl readline sqlite3 xz wget
```

## Ubuntu on Windows
```bash
$ sudo apt update
$ sudo apt-get install -y make build-essential libssl-dev  libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl vim debootstrap
```

## Let's install test package!
### On Ubuntu systems:  
```bash
$ apt install screenfetch
```

### On macOS systems:  
```bash
$ brew install screenfetch
```

```bash
screenfetch
```

### Programming languages

|-----|------|
|Compiled|FORTRAN, C, C++, Java|
|Interpretive|Unix-Shell, awk, Basic, Perl, Tcl, Scheme, Ruby, Python, R|



#### **Perl**
Flexible, by a global repository (CPAN), thus it is small install new modules. It has Bio per, one of the first biological unit repositories that upsurge the usability from, for instance, change setups to do the phylogenetic investigation. There are several biological software those usages Perl such as GBrowse thus might be an exciting language if you need toward interact with it. Good test units Perl is still what a lot of persons use, but it is declining out of use since Python accomplishes the similar tasks and is easier toward write code for, particularly for beginners.


#### **R & Python**
R is great for all the reasons, but if you like coding more than statistics, you may enjoy Python’s style a lot more. That sounds like a contradiction: How could you possibly know you enjoy coding more than statistics when you are choosing your first programming language? I would suggest trying them both and seeing what you like best. I personally enjoy coding in Python more than in R because its rules make more sense and it feels more like a programming language. In my experience, it is also much easier to make a command-line tool in Python than in R, and Python also has some packages for bioinformatics that are quite useful.


#### **Bash**
It is also very important for bioinformaticians to learn Bash, which for all of our intents and purposes is interchangeable with shell, the command-line, or the terminal. Bash is the primary way to access your data on your institution’s cluster and to run most genomics and bioinformatics software. It is also very powerful for manipulating your data like sorting, filtering, or doing calculations between columns, which is available through various utilities.

In my experience, and everyone I have talked to about it, bash was confusing and scary at first, but when you get the hang out it you start to feel this power surging through you, and you can do things in second that would take you hours to do by hand. Even two years into it I would still learn something new in bash that would blow my mind and I would kick myself for wasting time having programmed it from scratch in Python.

### Python, Perl, R and bash
In summary, for wet-lab people who want to add bioinformatics to their toolbox, focus on learning R first and applying it to your own work. For people who want to focus on bioinformatics as a career and make their own tools too, I would actually recommend learning the trifecta of R, Python, and Bash, though you could get away with choosing between R and Python as long as you still learn Bash too. I can go into more depth on any of these topics or give an introduction to any of these languages if you let me know in the comments.


### Other programming languages
There are many other languages out there, so before I end here I’m going to give a brief reason why these are not recommended for bioinformatics, beginners, or anyone at all in some cases.

#### C and C++
C or C++ are great for making super optimized command-line tools like aligners and variant-callers, but you will have a much easier time learning Python first and then going to these high-performance languages for a particular problem in the future, since they are harder to learn, more finicky, and take a lot more code to do the same thing.

#### Ruby
Ruby is one of those hot languages right now, for good reason largely because of the power of Ruby on Rails for making database-driven web applications like blogs or twitter. Ruby however is not great for bioinformatics because it lacks the community support in terms of packages that R and Python have, so you would be better off learning Python instead of Ruby.

#### JavaScript or PHP
JavaScript and PHP are great languages for web applications, but bioinformatics web applications should never be your first project. You could make a computational method in Python or R and then later make it into a web application, but that is not a project for a beginner. HTML and CSS by the way are not programming languages, but actually markup and styling languages that you will use along with JavaScript and PHP for that web application someday.

#### Java
Java is a popular language that most people have heard of. In bioinformatics, a notable example is the genome browser IGV. However, I would not recommend for beginners to learn Java due to many issues including memory management and that Python and R have many more bioinformaticians who build packages and answer questions online.


![language]({{{site.baseurl}}/fig/language.png)


### Package Library Module

#### Library 
Most often will refer to the general library or another collection created with a similar format and use. The General Library is the sum of 'standard', popular and widely used Modules, witch can be thought of as single file tools, for now or short cuts making things possible or faster. The general library is an option most people enable when installing Python. Because it has this name "Python General Library" it is used often with similar structure, and ideas. Witch is simply to have a bunch of Modules, maybe even packages grouped together, usually in a list. The list is usually to download them. Generally it is just related files, with similar interests. That is the easiest way to describe it.

#### Module 
A Module refers to a file. The file has script 'in it' and the name of the file is the name of the module, Python files end with .py. All the file contains is code that ran together makes something happen, by using functions, strings ect. Main modules you probably see most often are popular because they are special modules that can get info from other files/modules. It is confusing because the name of the file and module are equal and just drop the .py. Really it's just code you can use as a shortcut written by somebody to make something easier or possible.

#### Package
This is a termis used to generally sometimes, although context makes a difference. The most common use from my experience is multiple modules (or files) that are grouped together. Why they are grouped together can be for a few reasons, that is when context matters. These are ways I have noticed the term package(s) used. They are a group of Downloaded, created and/or stored modules. Which can all be true, or only 1, but really it is just a file that references other files, that need to be in the correct structure or format, and that entire sum is the package itself, installed or may have been included in the python general library. A package can contain modules(.py files) because they depend on each other and sometimes may not work correctly, or at all. There is always a common goal of every part (module/file) of a package, and the total sum of all of the parts is the package itself.

### Programming languages module or library manager
Python - pip  
Perl - cpan  
R - native manager  

## File Permission

![language]({{{site.baseurl}}/fig/linux_file_permissions.png)

### Understanding of attribute which can be out put by `ls -l`:
![language]({{{site.baseurl}}/fig/file_permission.png)
![language]({{{site.baseurl}}/fig/file_permission2.png)

### Chmod (Change mode)
chmod is the command and system call which is used to change the access permissions of file system objects (files and directories). It is also used to change special mode flags. The request is filtered by the umask. The name is an abbreviation of change mode.

### Applying Permission:
![language]({{{site.baseurl}}/fig/file_permission3.png)
![language]({{{site.baseurl}}/fig/file_permission4.png)

### Using Octal number:
![language]({{{site.baseurl}}/fig/file_permission5.png)



### Reference:

- Conda documentation https://docs.conda.io/en/latest/
- Conda-forge https://conda-forge.github.io/
- BioConda https://bioconda.github.io/

