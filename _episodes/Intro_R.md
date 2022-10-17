---
layout: page
title: 6_Introduction of R & R plotting 
published: true
---

{% include gh_variables.html %}



## PhD Tong Zhou
Tong Zhou, PhD  
L-207B, Center for Molecular Medicine, MS575  
Department of Physiology and Cell Biology  
University of Nevada, Reno School of Medicine  
1664 North Virginia Street, Reno , NV 89557  


The Zhou Lab carries out translational and theoretical research in bioinformatics and computational biology. Much of our research addresses questions of computational molecular medicine and molecular evolution, in particular about the use of genomic data to understand the pathobiology and develop biomarkers for human diseases and to understand the mechanism of exnoic sequence evolution.


## Introduction to R
R (www.r-project.org) is a commonly used free Statistics software. R allows you to carry out statistical analyses in an interactive mode, as well as allowing simple programming.
![R](https://cran.r-project.org/Rlogo.svg)

## Prepare your laptop
#### Open two terminal
##### Connect one to pronghorn
```ssh <yourID>@pronghorn.rc.unr.edu```

#### Open web browser
##### Connect to https://pastebin.com/6rba4Q9x


## R Installation

### Prepare working folder in Pronghorn
```bash 
cd ~/
mkdir r_plot
cd r_plot
wget -O r_plot.yaml https://pastebin.com/raw/kSAC1AsK
wget -O dataset1.txt https://pastebin.com/raw/N5g8bXg6
conda env create -n r_plot -f r_plot.yaml
conda activate r_plot
```

## Start R
```bash
R
```


## R Plot with Gene Expression Data Read
```bash
R
```

## Getting help with functions and features
### Getting information on any specific named function
```R
help(plot)
?plot
```
### Getting information on a feature specified by special characters
```R
help("%%")
? "%%"
```
### Launching a web browser for help
```R
help.start()
```

## Vectors and assignment
### Setting up a vector named a, consisting of five numbers
```R
a <- c(1, 4, 6, 7, 20)
a = c(1, 4, 6, 7, 20)
```
### Setting up a vector named b, consisting of two sequences of characters (strings)
```R
b <- c("Hello", "world")
b = c("Hello", "world")
```

### Setting up a logical vector named x
```R
x = c(TRUE, TRUE, FALSE, FALSE)
x = c(T, T, F, F)
```

## Arithmetic operators
```R
x = c(8, 4, 2, 1)
y = c(2, 2, 2, 2)
x + y
x - y
x * y
x / y
x ^ y
x %% y
```

### Functions mean(), sd(), and sum()
```R
mean(x)
sd(x)
sum(x
```

## Logical operators
```R
x = c(8, 4, 2, 1)
y = c(2, 2, 2, 2)
x > y; x >= y
x < y; x <= y
x == y
x != y
```

```R
a = 5
b = 0
a > 0 & b > 0 # and
a < 0 | b > 0 # or
! a == 5 # not
```
```bash
TRUE & TRUE == TRUE
TRUE & FALSE == FALSE
FALSE & FALSE == FALSE
TRUE | TRUE == TRUE
TRUE | FALSE == TRUE 
FALSE | FALSE == FALSE
! TRUE == FALSE
! FALSE == TRUE
```

## Generating regular sequences
### Generating a vector c(1, 2, …, 29, 30)
```R
1:30
1:30*2
# What is the difference between 1:30-1 and 1:(30-1) ?
```

### Function seq()
```R
seq(from=3, to=20); seq(3, 20); 3:20
seq(3, 20, by=0.5)
```

### Function rep()
```R
rep(-2, 20)
rep(c(1, 4, 5, 6), 10)
```


## Selecting subsets of a vector
```R
age = c(20, 19, 25, 26, 33, 24)
id= c("Tom", "Susan", "David", "Lauren", "Joe", "Nancy")
names(age) = id
```
### Pick up the first three elements of the “age” vector
```R
age[1:3]
age[c("Tom", "Susan", "David")]
```

### Pick up the last three elements of the “age” vector
```R
age[4:6]
age[-(1:3)]
age[c("Lauren", "Joe", "Nancy")]
```
## Exercise

### Generate a vector containing the odd numbers between 1 and 100
```R
1:50*2-1
seq(1, 100, by=2)
```
### Pick up the integers from 1 to 100 that can be divided by both 3 and 4
```R
x = 1:100
x[x %% 3 == 0 & x %% 4 == 0]
```
### Pick up the integers from 300 to 1000 that can be divided by either 12 or 16
```R
x = 300:1000
x[x %% 12 == 0 | x %% 16 == 0]
```

## Basic function: plot() 
### Scatter plot
```R
pdf("plot.pdf") 
plot(1:10, 1:10)
dev.off()
```
### Line plot

```R
pdf("lineplot.pdf") 
plot(1:10, 1:10, type="l")
dev.off()
```

### Plot a curve

```R
pdf("curveplot.pdf") 
x = seq(0, 2*pi, by=0.01)
plot(x, sin(x), type="l", col="red")
dev.off()
```


### Download results to ***your laptop***
### ***Try this on your laptop**

```bash
###Please do it on your desktop
###Please replace <YOURID> to your id.
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/plot.pdf .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/lineplot.pdf .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/curveplot.pdf .
```
## In your local terminal
### Windows
```explorer.exe .```
### MacOS
```open .```
### Git Bash
```explorer.exe .```



## Dataset
### Dataset (dataset1.txt)

```bash
cat dataset1.txt
```

```bash
Gene expression dataset of lung tissue samples (72×78)
72 genes
78 samples
30 samples from normal tissues
48 samples from tumor tissues
```

## Loading data frame

```R
#1 load the expression data sheet - 78 samples and 72 genes  
expr = read.delim("dataset1.txt", row.names="gene")  

#2 the first 30 samples are normal lung tissue while the last 48 samples are from lung tumors  
cl = c(rep("Normal", 30), rep("Tumor", 48))  
```

## Box plot of gene expression 

```R
#3 boxplot showing the difference in expression of gene TP63 between normal and tumor tissues  
pdf("T63_boxplot.pdf")  
boxplot(t(expr["TP63",])~cl, xlab="", ylab="Expression", main="TP63", col=c("grey80", "grey40"))  
dev.off()  

#4 boxplot showing the difference in expression of gene SLC2A1 between normal and tumor tissues
pdf("SLC2A1_boxplot2.pdf")  
boxplot(t(expr["SLC2A1",])~cl, xlab="", ylab="Expression", main="SLC2A1", col=c("grey80", "grey40"))  
dev.off()  
```

### Download results to ***your laptop***
### ***Try this on your laptop**

```bash
###Please do it on your desktop
###Please replace <YOURID> to your id.

scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/T63_boxplot.pdf .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/SLC2A1_boxplot2.pdf .
```
## In your local terminal
### Windows
```explorer.exe .```
### MacOS
```open .```
### Git Bash
```explorer.exe .```



## Co-expression pattern 1

```R
#5 scatter plot showing the co-expression pattern between gene TP63 and gene SLC2A1
pdf("TP63_SLC2A1_scatter_plot.pdf")  
plot(t(expr["TP63",]), t(expr["SLC2A1",]), xlab="Expression of TP63", ylab="Expression of SLC2A1", pch=19, col="blue", cex=1.2)  
dev.off()  

#6 correlation in expression between between gene TP63 and gene SLC2A1
cor.test(t(expr["TP63",]), t(expr["SLC2A1",]))  

```


### Download results to ***your laptop***
***Try this on your laptop***

```bash
###Please replace <YOURID> to your id.
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/TP63_SLC2A1_scatter_plot.pdf .
```

## Co-expression pattern 2
  
```R

#7 scatter plot showing the co-expression pattern between gene TP63 and gene TSHZ2
pdf("co-expression_pattern_between.pdf")
plot(t(expr["TP63",]), t(expr["TSHZ2",]), xlab="Expression of TP63", ylab="Expression of SLC2A1", pch=21, bg="orange", cex=1.5)
dev.off()

#8 correlation in expression between between gene TP63 and gene TSHZ2
cor.test(t(expr["TP63",]), t(expr["TSHZ2",]))
```


### Download results to ***your laptop***
***Try this on your laptop***

```bash
###Please replace <YOURID> to your id.
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/co-expression_pattern_between.pdf .
```


## Heatmap of the gene expression

```R

#9 gene expression heatmap
pdf("heatmap1.pdf")
heatmap(as.matrix(expr), col=cm.colors(256), scale='row')
dev.off()

pdf("heatmap2.pdf")
heatmap(as.matrix(expr), col=cm.colors(256), scale='row', labRow="", labCol="", ColSideColors=c(rep("slateblue", 30), rep("red", 48)))
dev.off()

pdf("heatmap3.pdf")
heatmap(as.matrix(expr), col=rainbow(256), scale='row', labRow="", labCol="", ColSideColors=c(rep("slateblue", 30), rep("red", 48)))
dev.off()

pdf("heatmap4.pdf")
heatmap(as.matrix(expr), col=heat.colors(256), scale='row', labRow="", labCol="", ColSideColors=c(rep("slateblue", 30), rep("red", 48)))
dev.off()

pdf("heatmap5.pdf")
heatmap(as.matrix(expr), col=c("darkblue", "blue", "lightblue", "white", "lightgreen", "green", "darkgreen"), scale='row', labRow="", labCol="", ColSideColors=c(rep("slateblue", 30), rep("red", 48)))
dev.off()
```

### Download results to ***your laptop***
***Try this on your laptop***

```bash
###Please replace <YOURID> to your id.
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/heatmap1.pdf .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/heatmap2.pdf  .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/heatmap3.pdf  .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/heatmap4.pdf .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/heatmap5.pdf .
```


## PCA of the gene expression

```R
#10 principal component analysis on the gene expression matrix
library(ade4) # if not installed yet, type install.packages(“ade4”)
pca_output = dudi.pca(t(expr), scannf=FALSE, nf=2)
pdf("PCA_plot.pdf")
plot(pca_output$l1$RS1[cl=="Normal"], pca_output$l1$RS2[cl=="Normal"], bg="slateblue", xlim=c(-2, 2), ylim=c(-2, 4.5), xlab='PC1', ylab='PC2', pch=21, cex=1.2)
points(pca_output$l1$RS1[cl=="Tumor"], pca_output$l1$RS2[cl=="Tumor"], bg="red", pch=21, cex=1.2)
legend('topleft', c("Normal tissue", "Tumor tissue"), pch=21, pt.bg=c("slateblue", "red"), pt.cex=1.2, bty='n', ncol=1)
dev.off()

```

### Download results to ***your laptop***
***Try this on your laptop***

```bash
###Please replace <YOURID> to your id.
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/PCA_plot.pdf .
```

## Histogram of the gene expression
```R
#11 histogram of the TP63 gene expression
pdf("histogram1.pdf")
hist(t(expr["TP63",]), xlab="Expression", main="")
dev.off()
```

```R
pdf("histogram2.pdf")
hist(t(expr["TP63",]), xlab="Expression", main="", breaks= -2:12*0.5, col="grey")
dev.off()
```

### Download results to ***your laptop***
***Try this on your laptop***

```bash
###Please replace <YOURID> to your id.
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/histogram1.pdf .
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/histogram2.pdf .
```


## Distribution of the gene expression

```R
#12 distribution of the TP63 gene expression
pdf("TP63_expression.pdf")
plot(density(t(expr["TP63",])), xlab="Expression", main="")
dev.off()
```


### Download results to ***your laptop***
***Try this on your laptop***

```bash
###Please replace <YOURID> to your id.
###Please replace <your_BCH709_Desktop> to your linked BCH709_desktop
scp <YOURID>@pronghorn.rc.unr.edu:~/r_plot/TP63_expression.pdf .
```



>## Reading materials
>The R Project for Statistical Computing https://www.r-project.org/
>R Graphic cookbook https://learning.oreilly.com/library/view/r-graphics-cookbook/9781491978597/
>R workshop https://bioinformatics.ca/workshops/2018-introduction-to-R/
>RStudio https://www.rstudio.com/
>Wikipedia R (programming language) https://en.wikipedia.org/wiki/R_(programming_language)
>Wikipedia RStudio https://en.wikipedia.org/wiki/RStudio
{: .prereq}

