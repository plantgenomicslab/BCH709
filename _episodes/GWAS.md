---
layout: page
title: GWAS
published: true
---
# Additional point quiz 2

Please read and follow below tutorials and answer the quiz. It will provide additional points.

## GWAS
https://www.genome.gov/17516714/2006-release-about-whole-genome-association-studies

## EasyGWAS
https://academic.oup.com/plcell/article/29/1/5/6099036

https://easygwas.ethz.ch/
1. Please go to EasyGWAS site

2. Click Perform a GWAS
![](https://i.imgur.com/9bCZIkt.png)

3. Create your account
![](https://i.imgur.com/joRmwTC.png)

4. Click Perform a GWAS then login
![](https://i.imgur.com/9bCZIkt.png)

5. Select following options (No human data)
![](https://i.imgur.com/68eG4sQ.png)

6. In the second step we have to select a Phenotype. Here we have the choice to select publicly available Phenotypes or private ones (if available) by searching the name of the phenotype using the autocompletion field. We can select a total of 5 Phenotypes. For this tutorial we select the two Phenotypes LD and LDV by typing their names. Then click Continue. Select phenotype, if you type `L` it will show multiple phenotypes.
![](https://i.imgur.com/Am8cwLr.png)

7. Select Two phenotype `LD` and `LDV`
![](https://i.imgur.com/oZzgSsy.png)

8. In the following view we can see histograms of the phenotypic distributions and a p-value from the Shapiro-Wilk test. The Shapiro-Wilk-Test tests the null hypothesis that the data was drawn from a normal distribution. Here we can transform the phenotype by applying a transformation to the phenotype. The number of available transformations is automatically-determined and depends on the distribution of the selected Phenotypes. For our selected Phenotypes we choose a Log10 transformation and click Continue.

9. To account for hidden confounding we can add Covariates or Principal Components to our experiments. However, this step is optional and in this tutorial we skip it by clicking Continue.
![](https://i.imgur.com/tjwDyCJ.png)

10. Here we can either select all available SNPs for our analysis or we select a set of chromosomes. For our purpose, we select all SNPs and click Continue.
![](https://i.imgur.com/JUiHrX1.png)

11. Next we have to choose the algorithms and filters we would like to use for our analysis. First we select a Minor Allele Frequency Filter of 10%. Next we select, for each Phenotype, the algorithm EMMAX and click Continue.

![](https://i.imgur.com/ehBgx0w.png)

12. In the last step we are asked to confirm if all the parameters we chose for our analyses are correct. We can always adjust them, if necessary. Finally, we can submit our GWAS to the computation servers by clicking Submit GWAS.
![](https://i.imgur.com/VYT7Dm8.png)

13. Now all GWAS are submitted. We can monitor the current status of our GWAS at My temporary history. This view is updated automatically and we will get a notification via e-mail when the computations have finalized. Will take 10 min. Go and drink water.
![](https://i.imgur.com/Pq83ZYm.png)

14. After drinking water, please click GWAS center then click my temporary history
![](https://i.imgur.com/Ze73rdr.png)

15. It will show your GWAS results
![](https://i.imgur.com/5wz5MD3.png)

16. Click LD experiment
![](https://i.imgur.com/5wz5MD3.png)

17. It will show Manhattan plot
![](https://i.imgur.com/Vl2bfOA.png)
A Manhattan plot is a type of plot, usually used to display data with a large number of data-points, many of non-zero amplitude, and with a distribution of higher-magnitude values. The plot is commonly used in genome-wide association studies (GWAS) to display significant SNPs.

Each chromosome will show individually and it will show p-value for each SNPs as dot and Bonferroni thresdhold as green line.

18. By clicking QQ-plots, it will show QQ-plots
 ![](https://i.imgur.com/6EuTZ3p.png)
 The Quantile-Quantile plot is a graphical representation of the deviation of the observed P values from the null hypothesis: the observed P values for each SNP are sorted from largest to smallest and plotted against expected values from a theoretical χ2-distribution.

Quantile-Quantile plot of the GWAS results. In this plot, each dot corresponds to a SNP tested for association where the observed –log10 p values, shown by vertical axis, were plotted by the expected –log10 p values under the null hypothesis. Upper right dots with higher observed significance than expected represent candidate variants for association with the phenotype tested. The genomic control ratio (λ) was 1.033, indicating the lack of strong effect of systematic error such as population stratification.


19. Back to Manhattan plot, click small dot above  Bonferroni thresdhold in Chromosome 2
![](https://i.imgur.com/Xjj0DMy.png)

20. It will show Phenotypic Values for LD and related SNPs that we clicked.

![](https://i.imgur.com/t5sweHR.png)

The allele G to C  might related to LD phenotype.

21. If you look at the table, it will show the effect.
It looks like it is on AT2G22530 and AT2G22560 upstream (promoter region). Please remember gene in genome is not always forwad stranded. It is also located in the intron of three different splicing transcripts (AT2G22540.1, AT2G22540.2 and AT2G22540.3).

HGVS is indicate the reference and sample. All variants are described in relation to a reference, the so called reference sequence, in the examples NM_004006.3 (from the GenBank database) NC_000023.11 (from the GenBank database). After the reference a description of the variant is given, in the examples c.4375C>T and g.32389644G>A.

### HGVS 
Changes in DNA, RNA and protein sequences, also called variants (and sometimes mutations or polymorphisms), are described using a specific language. To prevent confusion regarding its meaning a standard has been developed for this language, the so called HGVS nomenclature standard. The standard is used world-wide, especially in human health and clinical diagnostics. This page will try to explain the standard, briefly and in simple terms. After reading you should be able to understand the basics of the HGVS nomenclature and be able to use the internet to find more information about specific variants. In addition, while searching, you should be able to prevent making mistakes leading to misinterpretation of the variant and its possible consequences. More details, on all subjects, are availble elsewhere on the HGVS nomenclature pages.
https://3billion.io/blog/how-to-properly-read-variant-information-in-the-results-of-genetic-testing-for-rare-diseases-hgvs-nomenclature/

### HGVS c.261+21G>C means
nucleotides at the 5’ end of an intron are numbered relative to the last nucleotide of the directly upstream exon, followed by a “+” (plus) and their position in to the intron, like c.87+1, c.87+2, c.87+3, …
So c indicates coding and 261 bp from 5’ end 
21 bp from upstream exon 
then it changed to G to C.

21. If you click SNP annotation it will show P-value sorted table.
![](https://i.imgur.com/C6dhg69.png)

22. Only chromosome 2 has significant SNPs, this GWAS seems failed. Please pack your bag and run away from lab. Life is difficult.

23. However, we have one more phenotype 'LDV'
![](https://i.imgur.com/tdrSUAo.png)


24. This time we have two SNPs seems interesting.
![](https://i.imgur.com/VNRKZFm.png)
![](https://i.imgur.com/pQZ6hJm.png)

## Quiz
Please find the SNPs that related to LDV.
Describe which one will make sense to do focal CrisPR-Cas9