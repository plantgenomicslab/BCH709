---
layout: page
title: Symbolic Links
published: true
---

{% include gh_variables.html %}

## Symbolic Links

Symbolic links (symlinks) are shortcuts that point to other files or directories.

> ## Step 1: Create a Symbolic Link
> ```bash
> # Create a test file
> echo "Original content" > original_file.txt
> cat original_file.txt
> ```
> ```output
> Original content
> ```
> ```bash
> # Create a symbolic link
> ln -s original_file.txt my_shortcut.txt
>
> # View the link (note the -> arrow)
> ls -l my_shortcut.txt
> ```
> ```output
> lrwxrwxrwx 1 user group 17 Jan 20 10:00 my_shortcut.txt -> original_file.txt
> ```
> ```bash
> # Access through the link
> cat my_shortcut.txt
> ```
> ```output
> Original content
> ```
{: .keypoints}

> ## Step 2: Understand Link Behavior
> ```bash
> # Modify through the link
> echo "Added via link" >> my_shortcut.txt
>
> # Original file is modified!
> cat original_file.txt
> ```
> ```output
> Original content
> Added via link
> ```
> ```bash
> # Delete original - link becomes broken
> rm original_file.txt
> cat my_shortcut.txt
> ```
> ```output
> cat: my_shortcut.txt: No such file or directory
> ```
> ```bash
> # Link still exists but is broken (red in colored ls)
> ls -l my_shortcut.txt
> ```
> ```output
> lrwxrwxrwx 1 user group 17 Jan 20 10:00 my_shortcut.txt -> original_file.txt
> ```
{: .keypoints}

### Symbolic vs Hard Links

| Feature | Symbolic Link | Hard Link |
|---------|---------------|-----------|
| Command | `ln -s target link` | `ln target link` |
| Link to directories | Yes | No |
| Original deleted | Link breaks | Link works |
| Cross filesystems | Yes | No |

> ```bash
> # Create files to compare
> echo "test" > testfile.txt
> ln -s testfile.txt symlink.txt    # Symbolic
> ln testfile.txt hardlink.txt      # Hard
>
> # Check inode numbers (-i flag)
> ls -li testfile.txt symlink.txt hardlink.txt
> ```
> ```output
> 123456 -rw-r--r-- 2 user group 5 Jan 20 testfile.txt
> 123457 lrwxrwxrwx 1 user group 12 Jan 20 symlink.txt -> testfile.txt
> 123456 -rw-r--r-- 2 user group 5 Jan 20 hardlink.txt
> ```
> Note: hardlink.txt has the SAME inode (123456) as testfile.txt
{: .checklist}
