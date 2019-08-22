---
layout: lesson
title: Compile and Software installation
published: true
---

{% include gh_variables.html %}


![bioinformatics_DNA](https://www.google.com/url?sa=i&source=images&cd=&ved=2ahUKEwjI1Jqs0JfkAhVoFjQIHR6gBJQQjRx6BAgBEAQ&url=https%3A%2F%2Fwww.sciencemag.org%2Ffeatures%2F2014%2F06%2Fexplosion-bioinformatics-careers&psig=AOvVaw0JkZ3wkSOWtozvbuAykcHz&ust=1566602597354488)

When we talk about learning bioinformatics, it is useful to divide the students up into two groups: the ones who don’t want to make their own software and the ones who do. Both of these groups will do data analysis, run statistical tests, make plots, and use bioinformatics software made by other scientists. But the second group will also make their own bioinformatics software for the community to use. If you need to make some specialized scripts for your own research but you are not releasing anything for other researchers in your field to use, then you are in the first group.

Bioinformaticians who don’t build tools
For the first group, you are likely going to get the most use out of R. Some people are a little stuck up about R, saying it is not a “real” programming language, but it definitely is, and it has a lot of cool things built into it that also makes it ideal for bioinformatics.

It has a built-in data type called a data frame that has the same column and row setup as an Excel spreadsheet, where your genes, cells, people, time points, etc. will be rows while your variables are columns. This makes a lot of sense as a way to think about most kinds of data, so the Python people have made a package called Pandas to copy some of this functionality into Python, though it doesn’t work as smoothly as data frames do natively in R. The packages available for R to do bioinformatics are great, ranging from RNAseq to phylogenetic trees, and these are super easy to install from CRAN or the BioConductor.

If you use the free Rstudio software as your programming environment then it is even easier to manage what you are doing, and I would highly recommend Rstudio. Another major advantage of R is ggplot2, an awesome package for making plots that gives you results really quickly with even minimal coding skills. I made a video course about ggplot on my personal youtube channel, just search for Plotting in R for Biologists, which includes a good getting started guide for R in general.

Bioinformaticians who build tools
For bioinformaticians who make their own software, I would recommend either R or Python, plus bash.

R is great for all the reasons I just described, but if you like coding more than statistics, you may enjoy Python’s style a lot more. That sounds like a contradiction: How could you possibly know you enjoy coding more than statistics when you are choosing your first programming language? I would suggest trying them both and seeing what you like best. I personally enjoy coding in Python more than in R because its rules make more sense and it feels more like a programming language. In my experience, it is also much easier to make a command-line tool in Python than in R, and Python also has some packages for bioinformatics that are quite useful.

As you can probably tell, I have used both R and Python a lot in my work, where I use R for plotting and statistics, while I use Python for basically everything else, ranging from merging variant call sets to providing back-end algorithms for my web applications.

It is also very important for bioinformaticians to learn Bash, which for all of our intents and purposes is interchangeable with shell, the command-line, or the terminal. Bash is the primary way to access your data on your institution’s cluster and to run most genomics and bioinformatics software. It is also very powerful for manipulating your data like sorting, filtering, or doing calculations between columns, which is available through various utilities.

In my experience, and everyone I have talked to about it, bash was confusing and scary at first, but when you get the hang out it you start to feel this power surging through you, and you can do things in second that would take you hours to do by hand. Even two years into it I would still learn something new in bash that would blow my mind and I would kick myself for wasting time having programmed it from scratch in Python.

R, Python, and bash
In summary, for wet-lab people who want to add bioinformatics to their toolbox, focus on learning R first and applying it to your own work. For people who want to focus on bioinformatics as a career and make their own tools too, I would actually recommend learning the trifecta of R, Python, and Bash, though you could get away with choosing between R and Python as long as you still learn Bash too. I can go into more depth on any of these topics or give an introduction to any of these languages if you let me know in the comments.

Other programming languages
There are many other languages out there, so before I end here I’m going to give a brief reason why these are not recommended for bioinformatics, beginners, or anyone at all in some cases.

C and C++
C or C++ are great for making super optimized command-line tools like aligners and variant-callers, but you will have a much easier time learning Python first and then going to these high-performance languages for a particular problem in the future, since they are harder to learn, more finicky, and take a lot more code to do the same thing.

Perl
Perl is still what a lot of people use, but it is fading out of use because Python accomplishes the same tasks and is easier to write code for, especially for beginners.

Ruby
Ruby is one of those hot languages right now, for good reason largely because of the power of Ruby on Rails for making database-driven web applications like blogs or twitter. Ruby however is not great for bioinformatics because it lacks the community support in terms of packages that R and Python have, so you would be better off learning Python instead of Ruby.

JavaScript or PHP
JavaScript and PHP are great languages for web applications, but bioinformatics web applications should never be your first project. You could make a computational method in Python or R and then later make it into a web application, but that is not a project for a beginner. HTML and CSS by the way are not programming languages, but actually markup and styling languages that you will use along with JavaScript and PHP for that web application someday.

Java
Java is a popular language that most people have heard of. In bioinformatics, a notable example is the genome browser IGV. However, I would not recommend for beginners to learn Java due to many issues including memory management and that Python and R have many more bioinformaticians who build packages and answer questions online.

That’s all I have to say about bioinformatics programming languages for now. If you want to see more videos like this about bioinformatics, then make sure to subscribe on YouTube and sign up for updates below to get new videos, guides, and scripts about bioinformatics delivered to your email inbox every week.

And if you have a question you would like me to answer on the show, you can send it to me by going to omgenomics.com/tv and typing in your question there.

### Software

| Software | Version | Manual | Available for | Description |
| -------- | ------------ | ------ | ------------- | ----------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.7 | [Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)| Linux, MacOS, Windows | Quality control tool for high throughput sequence data. |
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.38 | [Link](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) | Linux, MacOS, Windows | A flexible read trimming tool for Illumina NGS data. |
| [BWA](http://bio-bwa.sourceforge.net/) | 0.7.17 | [Link](http://bio-bwa.sourceforge.net/bwa.shtml) | Linux, MacOS | Mapping DNA sequences against reference genome. |
| [SAMtools](http://samtools.sourceforge.net/) | 1.9 | [Link](http://www.htslib.org/doc/samtools.html) | Linux, MacOS | Utilities for manipulating alignments in the SAM format. |
| [BCFtools](https://samtools.github.io/bcftools/) | 1.8 | [Link](https://samtools.github.io/bcftools/bcftools.html) | Linux, MacOS | Utilities for variant calling and manipulating VCFs and BCFs. |
| [IGV](http://software.broadinstitute.org/software/igv/home) | [Link](https://software.broadinstitute.org/software/igv/download) | [Link](https://software.broadinstitute.org/software/igv/UserGuide) | Linux, MacOS, Windows | Visualization and interactive exploration of large genomics datasets. |

### QuickStart Software Installation Instructions
​
These are the QuickStart installation instructions. They assume familiarity with the command line and with installation in general. As there are different operating systems and many different versions of operating systems and environments, these may not work on your computer. If an installation doesn't work for you, please refer to the user guide for the tool, listed in the table above.
​
We have installed software using [miniconda](https://docs.conda.io/en/latest/miniconda.html). Miniconda is a package manager that simplifies the installation process. Please first install miniconda3 (installation instructions below), and then proceed to the installation of individual tools. 
​
### Miniconda3
​
> ## MacOS
> 
>To install miniconda3, type:
>
>~~~
>$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
>$ bash Miniconda3-latest-MacOSX-x86_64.sh
>~~~
>{: .bash}
> Then, follow the instructions that you are prompted with on the screen to install Miniconda3. 
{: .solution}
​
​
### FastQC
​
> ## MacOS
>
>To install FastQC, type:
>
> ~~~
> $ conda install -c bioconda fastqc=0.11.7=5
> ~~~
>{: .bash}
{: .solution}
​
> ## FastQC Source Code Installation
>
> If you prefer to install from source, follow the directions below:
>
> ~~~
> $ cd ~/src
> $ curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
> $ unzip fastqc_v0.11.7.zip
> ~~~
> {: .bash}
>
> Link the fastqc executable to the ~/bin folder that
> you have already added to the path.
>
> ~~~
> $ ln -sf ~/src/FastQC/fastqc ~/bin/fastqc
> ~~~
> {: .bash}
>
> Due to what seems a packaging error
> the executable flag on the fastqc program is not set.
> We need to set it ourselves.
>
> ~~~
> $ chmod +x ~/bin/fastqc
> ~~~
> {: .bash}
{: .solution}
​
**Test your installation by running:**
​
~~~
$ fastqc -h
~~~
{: .bash}
​
### Trimmomatic
​
> ## MacOS
>
> ~~~
> conda install -c bioconda trimmomatic=0.38=0
> ~~~
>{: .bash}
{: .solution}
​
> ## Trimmomatic Source Code Installation
>
> If you prefer to install from source, follow the directions below:
>
> ~~~
> $ cd ~/src
> $ curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
> $ unzip Trimmomatic-0.38.zip
> ~~~
> {: .bash}
>
> The program can be invoked via:
>
> ~~~
> $ java -jar ~/src/Trimmomatic-0.38/trimmomatic-0.38.jar
> ~~~
>
> The ~/src/Trimmomatic-0.38/adapters/ directory contains
> Illumina specific adapter sequences.
>
> ~~~
> $ ls ~/src/Trimmomatic-0.38/adapters/
> ~~~
> {: .bash}
{: .solution}
​
**Test your installation by running:** (assuming things are installed in ~/src)
​
~~~
$ java -jar ~/src/Trimmomatic-0.38/trimmomatic-0.38.jar
~~~
{: .bash}
​
​
> ## Simplify the Invocation, or to Test your installation if you installed with miniconda3:
>
> To simplify the invocation you could also create a script in the ~/bin folder:
>
> ~~~
> $ echo '#!/bin/bash' > ~/bin/trimmomatic
> $ echo 'java -jar ~/src/Trimmomatic-0.36/trimmomatic-0.36.jar $@' >> ~/bin/trimmomatic
> $ chmod +x ~/bin/trimmomatic
> ~~~
> {: .bash}
>
> Test your script by running:
>
> ~~~
> $ trimmomatic
> ~~~
> {: .bash}
{: .solution}
​
### BWA
​
> ## MacOS
>
>~~~
>conda install -c bioconda bwa=0.7.17=ha92aebf_3
>~~~
>{: .bash}
{: .solution}
​
> ## BWA Source Code Installation
>
> If you prefer to install from source, follow the instructions below:
>
> ~~~
> $ cd ~/src
> $ curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
> $ tar jxvf bwa-0.7.17.tar.bz2
> $ cd bwa-0.7.17
> $ make
> $ export PATH=~/src/bwa-0.7.17:$PATH
> ~~~
> {: .bash}
{: .solution}
​
**Test your installation by running:**
​
~~~
$ bwa
~~~
{: .bash}
​
### SAMtools
​
> ## MacOS
>
>~~~
>$ conda install -c bioconda samtools=1.9=h8ee4bcc_1
>~~~
>{: .bash}
{: .solution}
​
> ## SAMtools Versions
> SAMtools has changed the command line invocation (for the better). But this means that most of the tutorials
> on the web indicate an older and obsolete usage.
>
> Using SAMtools version 1.9 is important to work with the commands we present in these lessons.
{: .callout}
​
> ## SAMtools Source Code Installation
>
> If you prefer to install from source, follow the instructions below:
>
> ~~~
> $ cd ~/src
> $ curl -OkL https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
> $ tar jxvf samtools-1.9.tar.bz2
> $ cd samtools-1.9
> $ make
> ~~~
> {: .bash}
>
> Add directory to the path if necessary:
>
> ~~~
> $ echo export `PATH=~/src/samtools-1.9:$PATH` >> ~/.bashrc
> $ source ~/.bashrc
> ~~~
> {: .bash}
{: .solution}
​
**Test your installation by running:**
​
~~~
$ samtools
~~~
{: .bash}
​
​
### BCFtools
​
> ## MacOS
>
>~~~
>$ conda install -c bioconda bcftools=1.8=h4da6232_3 
>~~~
>{: .bash}
{: .solution}
​
> ## BCF tools Source Code Installation
>
> If you prefer to install from source, follow the instructions below:
>
> ~~~
> $ cd ~/src
> $ curl -OkL https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
> $ tar jxvf bcftools-1.8.tar.bz2
> $ cd bcftools-1.8
> $ make
> ~~~
> {: .bash}
>
> Add directory to the path if necessary:
>
> ~~~
> $ echo export `PATH=~/src/bcftools-1.8:$PATH` >> ~/.bashrc
> $ source ~/.bashrc
> ~~~
> {: .bash}
{: .solution}
​
**Test your installation by running:**
​
~~~
$ bcftools
~~~
{: .bash}
​
​
### IGV
​
- [Download the IGV installation files](https://software.broadinstitute.org/software/igv/download)
- Install and run IGV using the [instructions for your operating system](https://software.broadinstitute.org/software/igv/download).
​
