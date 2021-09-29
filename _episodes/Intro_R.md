---
layout: page
title: 6_Introduction of R & R plotting 
published: true
---

{% include gh_variables.html %}


## Introduction to R
R (www.r-project.org) is a commonly used free Statistics software. R allows you to carry out statistical analyses in an interactive mode, as well as allowing simple programming.
![R](https://cran.r-project.org/Rlogo.svg)

## Prepare your laptop
#### Open two terminal
##### Connect one to pronghorn
```ssh <yourID>@pronghorn.rc.unr.edu```

#### Open web browser
##### Connect to https://pastebin.com/6rba4Q9x
##### Connect to https://plantgenomicslab.github.io/BCH709/Intro_R/index.html



## R Installation
```bash
conda env create -n r_plot -f r_plot.yaml
conda activate r_plot
```

## Prepare working folder in Pronghorn
```bash 
mkdir r_plot
cd r_plot
wget -O r_plot.yaml https://pastebin.com/raw/kSAC1AsK
wget -O dataset1.txt https://pastebin.com/raw/N5g8bXg6
conda env create -n r_plot -f r_plot.yaml
conda activate r_plot
```
## In your local terminal
### Windows
```explorer.exe .```
### MacOS
```open .```
### Git Bash
```explorer.exe .```

## Start R
```bash
R
```

## PhD Tong Zhou
Tong Zhou, PhD  
L-207B, Center for Molecular Medicine, MS575  
Department of Physiology and Cell Biology  
University of Nevada, Reno School of Medicine  
1664 North Virginia Street, Reno , NV 89557  


The Zhou Lab carries out translational and theoretical research in bioinformatics and computational biology. Much of our research addresses questions of computational molecular medicine and molecular evolution, in particular about the use of genomic data to understand the pathobiology and develop biomarkers for human diseases and to understand the mechanism of exnoic sequence evolution.


## R Plot with Gene Expression Data Read
```bash
R
```


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

