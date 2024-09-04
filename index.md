---
---

The philosophy of BCH709: Introduction to Bioinformatics is to equip researchers with essential concepts, skills, and tools for handling data, enabling them to carry out their own bioinformatics tasks efficiently. This course covers key aspects of bioinformatics analysis, with a focus on Next Generation Sequencing (NGS) in genomics. It includes hands-on bioinformatics practice, understanding genomics through bioinformatics tools, using the command-line interface across different platforms, applying command-line tools for sequencing data analysis, and accessing cloud computing resources for bioinformatics workflows.


> ## Class Schedule
> MoWe 9:00AM - 10:15AM 
> Aug 26, 2024 - Dec 18 2023
{: .prereq}



> ## Modes of Instruction
> This course will be accommodated in-person. This semester we will not offer Online or Alternative HyFlex course model
{: .callout}



> ## Class Room location
> William N. Pennington Medical Education
> Building Code: PMB
> Building Number: 121
> Classroom: PMB 16
> ![Classroom location](./fig/classroom_location2.png){: width="50%" height="50%"}
{: .prereq}


> ## Course description
> As contemporary biologists we have entered an age where the use of computers in our daily work has become all but essential. The manipulation and analysis of DNA, RNA, and protein data by electronic means has become a routine task. Further, the amount of DNA, RNA and protein sequence data we are putting into databases every day is expanding at a geometric rate, and with coming advances in sequencing technology, this rate is only expected to increase. With all this new data, analysis by individual humans is simply not possible. Thus, in the past 15 years, computational biology has emerged as a field concerned with storage, manipulation, and extraction of valuable information from all this new data. However, because computational biology is an emerging field, organized courses are generally saved for higher-level study, and often are not required parts of an undergraduate curriculum. We seek to fill this void in education and create a course that will introduce students to bioinformatics at an earlier point in their education. This knowledge will prove to be not simply useful, but essential, for any student considering a degree in any area of biology and medical science.
{: .prereq}


> ## Syllabus  
Please read our [Syllabus](./syllabus/BCH709_2024fall.pdf).
{: .prereq}

> ## Rubric  
Please read our [Rubric](./syllabus/BCH709_rubric.pdf).
{: .prereq}


## Tentative Course Schedules

| Week   | Date       | Days      | Subject                                                           |
| ------ | ---------- | --------- | ----------------------------------------------------------------- |
| Week1  | 8/26/2024  | Monday    | Introduction                                                      |
| Week1  | 8/28/2024  | Wednesday | Introduction to Bioinformatics                                    |
| Week2  | 9/2/2024   | Monday    | Labor Day                                                         |
| Week2  | 9/4/2024   | Wednesday | Linux Environment and command line                                |
| Week3  | 9/9/2024   | Monday    | Linux Environment and command line                                |
| Week3  | 9/11/2024  | Wednesday | Gene family analysis and phylogenetics (David Alvarez-Ponce, PhD) |
| Week4  | 9/16/2024  | Monday    | Conda, Compile & Software Installations                           |
| Week4  | 9/18/2024  | Wednesday | Conda, Compile & Software Installations                           |
| Week5  | 9/23/2024  | Monday    | GitHub and server                                                 |
| Week5  | 9/25/2024  | Wednesday | Sequencing methods and strategies                                 |
| Week6  | 9/30/2024  | Monday    | Sequencing methods and strategies                                 |
| Week6  | 10/2/2024  | Wednesday | Sequence manipulation                                             |
| Week7  | 10/7/2024  | Monday    | Sequence manipulation                                             |
| Week7  | 10/9/2024  | Wednesday | Introduction of R & R plotting (Tong Zhou PhD)                    |
| Week8  | 10/14/2024 | Monday    | Transcriptome assembly                                            |
| Week8  | 10/16/2024 | Wednesday | RNA-Seq                                                           |
| Week9  | 10/21/2024 | Monday    | Midterm Exam                                                      |
| Week10 | 10/23/2024 | Wednesday | R in RNA-Seq / DESeq2 / EdgeR                                     |
| Week10 | 10/28/2024 | Monday    | R in RNA-Seq / DESeq2 / EdgeR                                     |
| Week11 | 10/30/2024 | Wednesday | Viral variant identification in NGS data (Richard Tillet, Ph. D)  |
| Week11 | 11/4/2024  | Monday    | BLAST search and gene alignment                                   |
| Week12 | 11/6/2024  | Wednesday | Genome assembly & annotation & structure                          |
| Week12 | 11/11/2024 | Monday    | Veterans Day                                                      |
| Week13 | 11/13/2024 | Wednesday | Variant analysis                                                  |
| Week13 | 11/18/2024 | Monday    | Transcriptome analysis (Genome based)                             |
| Week14 | 11/20/2024 | Wednesday | Nextday is Thanksgiving                                           |
| Week14 | 11/25/2024 | Monday    | Transcriptome analysis (Genome based)                             |
| Week15 | 11/27/2024 | Wednesday | Enrichment analysis                                               |
| Week15 | 12/2/2024  | Monday    | Presentation & Discussions                                        |
| Week16 | 12/4/2024  | Wednesday | Presentation & Discussions                                        |
| Week16 | 12/9/2024  | Monday    | Class Review                                                      |
| Week17 | 12/11/2024 | Wednesday | Prepday                                                           |
| Week17 | 12/16/2024 | Monday    | Final Exam                                                        |


> ## Prerequisites
> - Computer with ethernet port or wifi (If in case you bring your **desktop**, please do not bring your monitor. we have a monitor in our classroom)
> Online introduction to Linux. Students must complete one of the following online tutorials (or both) before class begins.
> - UNR affilated email **\<ID\>@unr.edu or \<ID\>@nevada.unr.edu** - [How to Activate](https://oit.unr.edu/services-and-support/login-ids-and-passwords/netid/netid-activation/)
> - [Setup your computer](https://plantgenomicslab.github.io/BCH709/setup.html)
> - [Setup Slack ID](https://unrrc.slack.com/)
> - [Setup Github ID](https://github.com/)
> - Please register by using UNR email [datacamp](https://www.datacamp.com/)
> - Please fill this [form](https://forms.gle/2sho6Nbh2PQ8akC1A)
>
{: .callout}



> ## Getting Started
>
> This lesson assumes that learners have no prior experience with the tools covered in the workshop. 
> However, learners are expected to have some familiarity with biological concepts,
> including the concept of genomic variation within a population. Participants should bring their own laptops and plan to participate actively. 
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



## Optional Additional Meeting

### _Research Computing Hackathon_  ([Hosted by HPC team](https://www.unr.edu/research-computing/hpc))  
Every Friday at 2:00pm to 4:00pm through SLACK

[Hackathons](https://en.wikipedia.org/wiki/Hackathon) provide a space for hands-on training and solution development within a Research Computing environment at the University. This is also a place to get clarification on questions/concerns regarding the HPC environment. Please bring problems to challenge the HPC team, the Office of Information Technology, and research colleagues. If you don't need help, we still encourage you to attend and share your time and expertise with those in need of assistance. You don’t need to be an expert to attend a hackathon. Individuals at all computing skill levels are welcome! Won Yim will attend this hackathon.

### _Meeting_  
Office [Howard Medical Science 216](https://goo.gl/maps/o41BMmcawsTPoES57)  
I prefer to have online meeting through SLACK.
  
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

[Plant Genomics Lab](https://www.plantgenomicslab.org/ "Plant Genomics Lab")  
[Lecture website](https://bch709.plantgenomicslab.org/)  
[Lecture Github](https://github.com/plantgenomicslab/BCH709)  
  
> ## Frequently Asked Questions  
Read our [FAQ](./_episodes/FAQ/FAQ.md). Currently, this page is empty, but we will build it through the class.
{: .prereq}

  
>## Teaching Platform
This lecture was designed to be run on [Unix-base system](https://en.wikipedia.org/wiki/Unix) such as 
Ubuntu, mac, etc. All the software and data used in the class will be open source. All example data will be hosted on a Google Cloud Service. If you want to know how to use Unix-base system on your computer, please follow the directions in the [Setup](setup.html) tab.
{: .prereq}

  

The website theme was adapted from the original by [Data Carpentry](https://datacarpentry.org/). The infrastructure, including adventure-time and docker-browser-server, was built by @maxogden and @mafintosh. The setup of this app was based on the get-data adventure. This adventure app was made by Richard Smith-Unna. The lecture materials were crafted by Won Yim. This work is licensed under a Creative Commons 4.0 International License.
