---
---

The philosophy of 'BCH709 Introduction to Bioinformatics' is to teach researchers basic concepts, skills, and tools for working with data so that student can get their own bioinformatics work done without pain. This teaches concept of bioinformatics analysis and Next Generation Sequencing for genomics including: hands-on bioinformatics, understanding genomics through bioinformatics, use of command-line in multiple platform, use of command-line tools to analyze sequencing data, and connecting to and using cloud computing.


> ## Course description
>As contemporary biologists we have entered an age where the use of computers in our daily work has become all but essential. The manipulation and analysis of DNA, RNA and protein data by electronic means has become a routine task. Further, the amount of DNA, RNA and protein sequence data we are putting into databases every day is expanding at a geometric rate, and with coming advances in sequencing technology this rate is only expected to increase. With all this new data, analysis by individual humans is simply not possible. Thus, in the past 15 years, computational biology has emerged as a field concerned with storage, manipulation, and extraction of valuable information from all this new data. However, because computational biology is an emerging field, organized courses are generally saved for higher-level study, and often are not required parts of an undergraduate curriculum. We seek to fill this void in education and create a course that will introduce students to bioinformatics at an earlier point in their education. This knowledge will prove to be not simply useful, but essential, for any student considering a degree in any area of biology and medical science.
{: .prereq}

> ## Class Room location
> ![Classroom location](./fig/classroom_location.png){: width="50%" height="50%"}
> 
> [PMB 104D](https://goo.gl/maps/GGx8NTwfyi8GZFpb8) (inside of Savitt Medical Library)  
{: .prereq}

> ## Class Schedule
> TuTh 9:00AM - 10:15AM  
> Aug 26, 2019 - Dec 10, 2019
{: .prereq}

> ## Prerequisites
> Online introduction to Linux. Students must complete one of the following online tutorials (or both) before class begins.
> - [Code academy's Intro to Unix](https://www.codecademy.com/learn/learn-the-command-line "Code academy")
> - Computer with ethernet port or wifi (If in case you bring your **desktop**, please do not bring your monitor. we have a monitor in our classroom)
>
> - Terminal software
> - UNR affilated email **\<ID\>@unr.edu or \<ID\>@nevada.unr.edu**
{: .callout}

> ## Getting Started
>
> This lesson assumes that learners have no prior experience with the tools covered in the workshop. 
> However, learners are expected to have some familiarity with biological concepts,
> including the 
> concept of genomic variation within a population. Participants should bring their own laptops and plan to participate actively. 
> 
> To get started, follow the directions in the [Setup](setup.html) tab to 
> get access to the required software and data for this workshop.
> 
{: .prereq}
<!-- 
> ## Data
> 
> This workshop uses data from a long term evolution experiment published in 2016: [Tempo and mode of genome evolution in a 50,000-generation experiment](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4988878/) by Tenaillon O, Barrick JE, Ribeck N, Deatherage DE, Blanchard JL, Dasgupta A, Wu GC, Wielgoss S, Cruveiller S, Médigue C, Schneider D, and Lenski RE. (doi: 10.1038/nature18959)
>
> All of the data used in this workshop can be [downloaded from Figshare](https://figshare.com/articles/Data_Carpentry_Genomics_beta_2_0/7726454). 
> More information about this data is available on the [Data page](https://datacarpentry.org/organization-genomics/data/).
{: .prereq} 
-->

## Class Overview 

| Week    | Tuesday | Thursday|
| ------- | ---------- |---------- |
|Week1|Intro to Bioinformatics|Linux Enviroment and command line|
|Week2|Linux Enviroment and bioawk|Linux Enviroment and Cloud|
|Week3|Github|Compile|
|Week4|Sequencing methods and strategies|RNA-Seq|
|Week5|Transcriptome assembly|Sequence similarity searching|
|Week6|Functional annotation|Database|
|Week7|Introduction of R &  R plotting ([Tong Zhou PhD](http://tongzhoulab.org/))|DEG analaysis
|Week8|DESeq2 / EdgeR|Special topics ([Tong Zhou PhD](http://tongzhoulab.org/))|
|Week9|Quiz: Transcriptome & Database|Gene family analysis and phylogenetics ([David Alvarez-Ponce, PhD](https://genomeevol.wordpress.com/))|
|Week10|Genome assembly|Genome assembly|
|Week11|Genome assembly|How to annotate genomes|
|Week12|Genome assembly and annotation|Genome structure|
|Week13|DEG analaysis|Transcriptome analaysis
|Week14|Variant analaysis|[Thanksgiving](https://giphy.com/search/dance)|
|Week15|Enrichment analysis|Presentation|
|Week16|Final exam|

## Optional Additional Meeting

### _Research Computing Hackathon_  ([Hosted by HPC team](https://www.unr.edu/research-computing/hpc))  
[Mathewson-IGT Knowledge Center, MIKC 414](https://events.unr.edu/mathewson-igt_knowledge_center_508#.XVyb3OhKiiM)  
Every Friday at 2:00pm to 4:00pm

[Hackathons](https://en.wikipedia.org/wiki/Hackathon) provide a space for hands-on training and solution development within a Research Computing environment at the University. This is also a place to get clarification on questions/concerns regarding the HPC environment. Please bring problems to challenge the HPC team, the Office of Information Technology, and research colleagues. If you don't need help, we still encourage you to attend and share your time and expertise with those in need of assistance. You don’t need to be an expert to attend a hackathon. Individuals at all computing skill levels are welcome! Won Yim will attend this hackathon.

### _Meeting with Won Yim_  
[Howard Molecular Science 216](https://goo.gl/maps/o41BMmcawsTPoES57)  
After Thursday class, we can have a meeting. We are not providing financial aid, but we can disccuss your research bottleneck in terms of bioinformatics perspective. **Please schedule after Tuesday class.**   
  
## Optional Reading materials
- [Introduction to Bioinformatics (3rd Edition)](http://app.knovel.com/web/toc.v/cid:kpIBE00007/viewerType:toc/ "Introduction to Bioinformatics (3rd Edition)")

- [Learn Linux Shell Scripting - Fundamentals of Bash 4.4](https://learning.oreilly.com/library/view/learn-linux-shell/9781788995597/ "Learn Linux Shell Scripting - Fundamentals of Bash 4.4")

- [Effective awk Programming, 3rd Edition](https://learning.oreilly.com/library/view/effective-awk-programming/0596000707/ "Effective awk Programming, 3rd Edition")

- [Bioinformatics with Python Cookbook - Second Edition](https://learning.oreilly.com/library/view/bioinformatics-with-python/9781789344691/ "Bioinformatics with Python Cookbook - Second Edition")

- [Basic Applied Bioinformatics](https://learning.oreilly.com/library/view/basic-applied-bioinformatics/9781119244332/ "Basic Applied Bioinformatics")

- [Ubuntu Unleashed 2019 Edition: Covering 18.04, 18.10, 19.04, 13/e](https://learning.oreilly.com/library/view/ubuntu-unleashed-2019/9780134985497/ "Ubuntu Unleashed 2019 Edition: Covering 18.04, 18.10, 19.04, 13/e")

- [Unix and Perl: Keith Bradnam & Ian Korf](https://j.p.gogarten.uconn.edu/mcb5472_2018/current.pdf)

- [Other bioinformatics books from knovel](http://app.knovel.com/web/search.v?q=bioinformatics&search_type=tech-reference&rows=10&offset=0&group_by=true&my_subscription=true&sort_on=default&content_type=all_references&include_synonyms=no "Other bioinformatics books from knovel")

- [Other bioinformatics books from Oreilly](https://learning.oreilly.com/search/?query=bioinformatics&extended_publisher_data=true&highlight=true&include_assessments=false&include_case_studies=true&include_courses=true&include_orioles=true&include_playlists=true&include_collections=false&include_notebooks=false&is_academic_institution_account=false&sort=relevance&facet_json=true "Other bioinformatics books from Oreilly")

**Note: all reading material can be freely accessed and downloaded from the UNR internet.**

## Other Website
[Command line website](http://35.199.189.11:8081/  "Command line website")  
[Plant Genomics Lab](https://www.plantbioinformatics.org/ "Plant Genomics Lab")

  
> ## Frequently Asked Questions  
Read our [FAQ](./_episodes/FAQ/FAQ.md). Currently this page is empty but we will build it through the class.
{: .prereq}

  
>## Teaching Platform
This workshop is designed to be run on [Unix-base system](https://en.wikipedia.org/wiki/Unix) such as 
Ubuntu, mac etc. All the software and data used in the class will be open source All example data will be hosted on an Google Cloud Service. If you want to know how to use Unix-base system on your computer, please follow the directions in the [Setup](setup.html) tab.
{: .prereq}

  

The website theme was adapted from the original by [Data Carpentry](https://datacarpentry.org/). The infrastructure, including adventure-time and docker-browser-server, was built by @maxogden and @mafintosh. The setup of this app was based on the get-dat adventure. This adventure was made by Richard Smith-Unna. The lecture materials were crafted by Won Yim. This work is licensed under a Creative Commons 4.0 International License.
