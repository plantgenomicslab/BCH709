---
layout: page
title: Find Command
published: true
---

{% include gh_variables.html %}

## Find Command

The `find` command searches for files based on name, size, time, and more.

> ## Setup: Create Test Files
> ```bash
> cd ~/permission_practice
> mkdir -p find_test/{dir1,dir2}
> touch find_test/file1.txt find_test/file2.txt find_test/script.sh
> touch find_test/dir1/data.csv find_test/dir2/notes.txt
> echo "sample content" > find_test/file1.txt
> ```
{: .prereq}

> ## Step 1: Find by Name
> ```bash
> # Find all .txt files
> find find_test -name "*.txt"
> ```
> ```output
> find_test/file1.txt
> find_test/file2.txt
> find_test/dir2/notes.txt
> ```
> ```bash
> # Case-insensitive search
> find find_test -iname "*.TXT"
> ```
{: .keypoints}

> ## Step 2: Find by Type
> ```bash
> # Find only files
> find find_test -type f
> ```
> ```output
> find_test/file1.txt
> find_test/file2.txt
> find_test/script.sh
> find_test/dir1/data.csv
> find_test/dir2/notes.txt
> ```
> ```bash
> # Find only directories
> find find_test -type d
> ```
> ```output
> find_test
> find_test/dir1
> find_test/dir2
> ```
{: .keypoints}

> ## Step 3: Find with Actions
> ```bash
> # Make all .sh files executable
> find find_test -name "*.sh" -exec chmod +x {} \;
> ls -l find_test/script.sh
> ```
> ```output
> -rwxr-xr-x 1 user group 0 Jan 20 10:00 find_test/script.sh
> ```
> ```bash
> # Delete all .txt files (careful!)
> # Use -print first to preview
> find find_test -name "*.txt" -print
> ```
{: .keypoints}

### Quick Reference

| Option | Description | Example |
|--------|-------------|---------|
| `-name` | Match filename | `find . -name "*.fa"` |
| `-type f` | Files only | `find . -type f` |
| `-type d` | Directories only | `find . -type d` |
| `-size +10M` | Larger than 10MB | `find . -size +10M` |
| `-mtime -7` | Modified last 7 days | `find . -mtime -7` |
