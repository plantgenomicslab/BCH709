---
layout: page
title: Standard Streams and Redirection
published: true
---

{% include gh_variables.html %}

## Standard Streams and Redirection

Every command has three data streams: input (stdin), output (stdout), and errors (stderr).

### The Three Streams

| Stream | Number | Description |
|--------|--------|-------------|
| stdin | 0 | Input (from keyboard/file) |
| stdout | 1 | Normal output |
| stderr | 2 | Error messages |

> ## Step 1: Output Redirection
> ```bash
> # Redirect output to file (overwrites)
> echo "Hello" > output.txt
> cat output.txt
> ```
> ```output
> Hello
> ```
> ```bash
> # Append to file (>>)
> echo "World" >> output.txt
> cat output.txt
> ```
> ```output
> Hello
> World
> ```
> ```bash
> # Redirect errors separately
> ls nonexistent 2> errors.txt
> cat errors.txt
> ```
> ```output
> ls: cannot access 'nonexistent': No such file or directory
> ```
> ```bash
> # Redirect both stdout and stderr
> ls genes.txt nonexistent > all.txt 2>&1
> cat all.txt
> ```
> ```output
> ls: cannot access 'nonexistent': No such file or directory
> genes.txt
> ```
{: .keypoints}

> ## Step 2: Input Redirection
> ```bash
> # Read input from file
> wc -l < genes.txt
> ```
> ```output
> 5
> ```
{: .keypoints}

> ## Step 3: Pipes
> Connect commands: output of one becomes input of next.
> ```bash
> # Count chr1 genes in our data
> cat genes.txt | grep "chr1" | wc -l
> ```
> ```output
> 2
> ```
> ```bash
> # Sort by score (column 5), show top 3
> sort -k5 -rn genes.txt | head -3
> ```
> ```output
> chr2	800	1200	geneD	92.5
> chr1	600	900	geneB	78.9
> chr1	100	500	geneA	45.2
> ```
> ```bash
> # Save intermediate result with tee
> cat genes.txt | grep "chr1" | tee chr1_genes.txt | wc -l
> cat chr1_genes.txt
> ```
> ```output
> 2
> chr1	100	500	geneA	45.2
> chr1	600	900	geneB	78.9
> ```
{: .keypoints}

### Redirection Quick Reference

| Symbol | Description | Example |
|--------|-------------|---------|
| `>` | Redirect stdout (overwrite) | `cmd > file` |
| `>>` | Redirect stdout (append) | `cmd >> file` |
| `2>` | Redirect stderr | `cmd 2> errors` |
| `&>` | Redirect both | `cmd &> all` |
| `<` | Read from file | `cmd < file` |
| `\|` | Pipe to next command | `cmd1 \| cmd2` |
