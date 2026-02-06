---
layout: page
title: sed - Stream Editor
published: true
---

{% include gh_variables.html %}

## sed - Stream Editor

`sed` performs text transformations on files. Great for find-and-replace.

> ## Setup: Create Sample File
> ```bash
> cat > sample.txt << 'EOF'
> Hello World
> hello world
> The quick brown fox
> Line with spaces
> Another line
> EOF
> cat sample.txt
> ```
> ```output
> Hello World
> hello world
> The quick brown fox
> Line with spaces
> Another line
> ```
{: .prereq}

> ## Step 1: Basic Substitution
> ```bash
> # Replace first 'o' with 'O' on each line
> sed 's/o/O/' sample.txt
> ```
> ```output
> HellO World
> hellO world
> The quick brOwn fox
> Line with spaces
> AnOther line
> ```
> ```bash
> # Replace ALL 'o' with 'O' (global flag)
> sed 's/o/O/g' sample.txt
> ```
> ```output
> HellO WOrld
> hellO wOrld
> The quick brOwn fOx
> Line with spaces
> AnOther line
> ```
> ```bash
> # Case-insensitive replace
> sed 's/hello/Hi/i' sample.txt
> ```
> ```output
> Hi World
> Hi world
> The quick brown fox
> Line with spaces
> Another line
> ```
{: .keypoints}

> ## Step 2: Line Operations
> ```bash
> # Print only line 3
> sed -n '3p' sample.txt
> ```
> ```output
> The quick brown fox
> ```
> ```bash
> # Print lines 2-4
> sed -n '2,4p' sample.txt
> ```
> ```output
> hello world
> The quick brown fox
> Line with spaces
> ```
> ```bash
> # Delete lines containing 'world'
> sed '/world/d' sample.txt
> ```
> ```output
> The quick brown fox
> Line with spaces
> Another line
> ```
{: .keypoints}

> ## Step 3: Edit File In-Place
> ```bash
> # Create a backup and edit
> cp sample.txt sample_backup.txt
> sed -i 's/World/Universe/g' sample.txt
> cat sample.txt
> ```
> ```output
> Hello Universe
> hello world
> The quick brown fox
> Line with spaces
> Another line
> ```
{: .checklist}

### sed Quick Reference

| Command | Description |
|---------|-------------|
| `s/old/new/` | Replace first match |
| `s/old/new/g` | Replace all matches |
| `s/old/new/i` | Case-insensitive |
| `-n '5p'` | Print line 5 only |
| `/pattern/d` | Delete matching lines |
| `-i` | Edit file in-place |
