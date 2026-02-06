---
layout: page
title: awk - Pattern Scanning and Processing
published: true
---

{% include gh_variables.html %}

## awk - Pattern Scanning and Processing

`awk` processes structured text data (columns). Essential for bioinformatics files.

> ## Setup: Create Sample Data
> ```bash
> cat > genes.txt << 'EOF'
> chr1	100	500	geneA	45.2
> chr1	600	900	geneB	78.9
> chr2	200	400	geneC	23.1
> chr2	800	1200	geneD	92.5
> chr3	150	350	geneE	15.8
> EOF
> cat genes.txt
> ```
> ```output
> chr1	100	500	geneA	45.2
> chr1	600	900	geneB	78.9
> chr2	200	400	geneC	23.1
> chr2	800	1200	geneD	92.5
> chr3	150	350	geneE	15.8
> ```
{: .prereq}

> ## Step 1: Print Columns
> ```bash
> # Print first column (chromosome)
> awk '{print $1}' genes.txt
> ```
> ```output
> chr1
> chr1
> chr2
> chr2
> chr3
> ```
> ```bash
> # Print columns 1 and 4 (chromosome and gene name)
> awk '{print $1, $4}' genes.txt
> ```
> ```output
> chr1 geneA
> chr1 geneB
> chr2 geneC
> chr2 geneD
> chr3 geneE
> ```
> ```bash
> # Print last column
> awk '{print $NF}' genes.txt
> ```
> ```output
> 45.2
> 78.9
> 23.1
> 92.5
> 15.8
> ```
{: .keypoints}

> ## Step 2: Filter with Conditions
> ```bash
> # Print only chr1 genes
> awk '$1 == "chr1"' genes.txt
> ```
> ```output
> chr1	100	500	geneA	45.2
> chr1	600	900	geneB	78.9
> ```
> ```bash
> # Print genes with score > 50
> awk '$5 > 50' genes.txt
> ```
> ```output
> chr1	600	900	geneB	78.9
> chr2	800	1200	geneD	92.5
> ```
> ```bash
> # Combine conditions (chr2 AND score > 50)
> awk '$1 == "chr2" && $5 > 50' genes.txt
> ```
> ```output
> chr2	800	1200	geneD	92.5
> ```
{: .keypoints}

> ## Step 3: Calculations
> ```bash
> # Sum of scores (column 5)
> awk '{sum += $5} END {print "Total:", sum}' genes.txt
> ```
> ```output
> Total: 255.5
> ```
> ```bash
> # Average score
> awk '{sum += $5; count++} END {print "Average:", sum/count}' genes.txt
> ```
> ```output
> Average: 51.1
> ```
> ```bash
> # Calculate gene length (end - start)
> awk '{print $4, $3 - $2}' genes.txt
> ```
> ```output
> geneA 400
> geneB 300
> geneC 200
> geneD 400
> geneE 200
> ```
{: .keypoints}

### awk Built-in Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `$0` | Entire line | `awk '{print $0}'` |
| `$1, $2...` | Column 1, 2, etc. | `awk '{print $1}'` |
| `NF` | Number of columns | `awk '{print NF}'` |
| `NR` | Line number | `awk '{print NR, $0}'` |
| `-F` | Set delimiter | `awk -F',' '{print $1}'` |

> ## Challenge: Analyze Gene Data
> Using genes.txt, find:
> 1. All genes on chr2
> 2. The gene with highest score
> 3. Total length of all genes
>
> > ## Solutions
> > ```bash
> > # 1. Genes on chr2
> > awk '$1 == "chr2" {print $4}' genes.txt
> > ```
> > ```output
> > geneC
> > geneD
> > ```
> > ```bash
> > # 2. Gene with highest score
> > awk 'NR==1 || $5 > max {max=$5; gene=$4} END {print gene, max}' genes.txt
> > ```
> > ```output
> > geneD 92.5
> > ```
> > ```bash
> > # 3. Total gene length
> > awk '{len += $3 - $2} END {print "Total length:", len}' genes.txt
> > ```
> > ```output
> > Total length: 1500
> > ```
> {: .solution}
{: .challenge}
