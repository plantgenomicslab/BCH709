---
layout: page
title: Shell Scripting Basics
published: true
---

{% include gh_variables.html %}

## Shell Scripting Basics

Automate repetitive tasks by combining commands into scripts.

> ## Step 1: Create Your First Script
> ```bash
> # Create a simple script
> cat > hello.sh << 'EOF'
> #!/bin/bash
> echo "Hello, World!"
> echo "Today is $(date +%Y-%m-%d)"
> EOF
>
> # View it
> cat hello.sh
> ```
> ```output
> #!/bin/bash
> echo "Hello, World!"
> echo "Today is $(date +%Y-%m-%d)"
> ```
> ```bash
> # Make it executable and run
> chmod +x hello.sh
> ./hello.sh
> ```
> ```output
> Hello, World!
> Today is 2026-01-20
> ```
{: .keypoints}

> ## Step 2: Using Variables
> ```bash
> cat > variables.sh << 'EOF'
> #!/bin/bash
> # Variables (no spaces around =)
> NAME="Student"
> COUNT=5
>
> echo "Hello, $NAME"
> echo "Count is $COUNT"
>
> # Command substitution
> TODAY=$(date +%Y-%m-%d)
> NUM_FILES=$(ls | wc -l)
> echo "Today: $TODAY, Files in directory: $NUM_FILES"
> EOF
>
> chmod +x variables.sh
> ./variables.sh
> ```
> ```output
> Hello, Student
> Count is 5
> Today: 2026-01-20, Files in directory: 12
> ```
{: .keypoints}

> ## Step 3: Command Line Arguments
> ```bash
> cat > args.sh << 'EOF'
> #!/bin/bash
> echo "Script: $0"
> echo "First arg: $1"
> echo "Second arg: $2"
> echo "All args: $@"
> echo "Num args: $#"
> EOF
>
> chmod +x args.sh
> ./args.sh apple banana cherry
> ```
> ```output
> Script: ./args.sh
> First arg: apple
> Second arg: banana
> All args: apple banana cherry
> Num args: 3
> ```
{: .keypoints}

> ## Step 4: If Statements
> ```bash
> cat > checker.sh << 'EOF'
> #!/bin/bash
> if [ -z "$1" ]; then
>     echo "Usage: $0 <filename>"
>     exit 1
> fi
>
> if [ -f "$1" ]; then
>     echo "$1 is a file"
> elif [ -d "$1" ]; then
>     echo "$1 is a directory"
> else
>     echo "$1 does not exist"
> fi
> EOF
>
> chmod +x checker.sh
> ./checker.sh genes.txt
> ```
> ```output
> genes.txt is a file
> ```
> ```bash
> ./checker.sh find_test
> ```
> ```output
> find_test is a directory
> ```
>
> **Test operators:** `-f` (file exists), `-d` (directory), `-z` (empty string), `-eq` (equal), `-gt` (greater than)
{: .keypoints}

> ## Step 5: For Loops
> ```bash
> cat > loop.sh << 'EOF'
> #!/bin/bash
> # Loop over files
> for file in *.txt; do
>     echo "Found: $file"
> done
>
> # Loop with range
> for i in {1..3}; do
>     echo "Count: $i"
> done
> EOF
>
> chmod +x loop.sh
> ./loop.sh
> ```
> ```output
> Found: genes.txt
> Found: sample.txt
> Found: testfile.txt
> Count: 1
> Count: 2
> Count: 3
> ```
{: .keypoints}

> ## Challenge: Bioinformatics Script
> Create a script that counts lines in all .txt files:
> ```bash
> cat > count_lines.sh << 'EOF'
> #!/bin/bash
> for file in *.txt; do
>     lines=$(wc -l < "$file")
>     echo "$file: $lines lines"
> done
> EOF
>
> chmod +x count_lines.sh
> ./count_lines.sh
> ```
> ```output
> genes.txt: 5 lines
> sample.txt: 5 lines
> testfile.txt: 1 lines
> ```
{: .challenge}
