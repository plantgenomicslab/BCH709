---
layout: page
title: Environment Variables
published: true
---

{% include gh_variables.html %}

## Environment Variables

Environment variables store system settings and user preferences that programs can access.

> ## Step 1: View Your Environment
> ```bash
> # See common variables
> echo "Home directory: $HOME"
> echo "Current user: $USER"
> echo "Current shell: $SHELL"
> echo "Current directory: $PWD"
> ```
> ```output
> Home directory: /home/username
> Current user: username
> Current shell: /bin/bash
> Current directory: /home/username/permission_practice
> ```
> ```bash
> # View PATH (where system looks for commands)
> echo $PATH
> ```
> ```output
> /usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin
> ```
{: .keypoints}

### Common Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `HOME` | User's home directory | `/home/username` |
| `PATH` | Search path for commands | `/usr/bin:/bin` |
| `USER` | Current username | `username` |
| `SHELL` | Current shell | `/bin/bash` |
| `PWD` | Current directory | `/home/username` |

> ## Step 2: Create and Use Variables
> ```bash
> # Create a variable (no spaces around =)
> MYNAME="Student"
> echo "Hello, $MYNAME"
> ```
> ```output
> Hello, Student
> ```
> ```bash
> # Export makes it available to child processes
> export PROJECT_DIR="$HOME/myproject"
> echo $PROJECT_DIR
> ```
> ```output
> /home/username/myproject
> ```
> ```bash
> # Unset removes a variable
> unset MYNAME
> echo "Name is: $MYNAME"
> ```
> ```output
> Name is:
> ```
{: .keypoints}

> ## Step 3: Make Variables Permanent
> Add to `~/.bashrc` for permanent variables:
> ```bash
> # View current bashrc (last 5 lines)
> tail -5 ~/.bashrc
>
> # Add a custom variable (be careful with >>)
> echo 'export BIOINF_DATA="$HOME/biodata"' >> ~/.bashrc
>
> # Reload bashrc
> source ~/.bashrc
>
> # Verify
> echo $BIOINF_DATA
> ```
> ```output
> /home/username/biodata
> ```
{: .checklist}
