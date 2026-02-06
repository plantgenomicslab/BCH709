---
layout: page
title: Process Management
published: true
---

{% include gh_variables.html %}

## Process Management

Managing processes is essential for running long bioinformatics analyses.

> ## Step 1: View Running Processes
> ```bash
> # Show your current processes
> ps
> ```
> ```output
>   PID TTY          TIME CMD
> 12345 pts/0    00:00:00 bash
> 12400 pts/0    00:00:00 ps
> ```
> ```bash
> # Show all processes (abbreviated output)
> ps aux | head -5
> ```
> ```output
> USER       PID %CPU %MEM    VSZ   RSS TTY   STAT START   TIME COMMAND
> root         1  0.0  0.1  16894  1340 ?     Ss   Jan19   0:02 /sbin/init
> root         2  0.0  0.0      0     0 ?     S    Jan19   0:00 [kthreadd]
> ...
> ```
> ```bash
> # Interactive view (press 'q' to quit)
> top
> ```
{: .keypoints}

> ## Step 2: Run Commands in Background
> ```bash
> # Start a long-running process (sleep simulates a long job)
> sleep 60 &
> ```
> ```output
> [1] 12456
> ```
> ```bash
> # List background jobs
> jobs
> ```
> ```output
> [1]+  Running                 sleep 60 &
> ```
> ```bash
> # Bring to foreground
> fg %1
> # Press Ctrl+C to stop, or Ctrl+Z to suspend
> ```
{: .keypoints}

> ## Step 3: Kill Processes
> ```bash
> # Start a background process
> sleep 300 &
> ```
> ```output
> [1] 12500
> ```
> ```bash
> # Kill by job number
> kill %1
> jobs
> ```
> ```output
> [1]+  Terminated              sleep 300
> ```
> ```bash
> # Or kill by PID
> sleep 300 &
> ps | grep sleep
> ```
> ```output
> 12510 pts/0    00:00:00 sleep
> ```
> ```bash
> kill 12510
> ```
{: .keypoints}

> ## Step 4: Keep Processes Running After Logout
> ```bash
> # nohup keeps process running after you log out
> nohup sleep 120 > mysleep.log 2>&1 &
> ```
> ```output
> [1] 12600
> ```
> ```bash
> # Check it's running
> jobs
> ps aux | grep sleep
>
> # You can now log out and the process continues!
> ```
{: .checklist}

### Process Signals Reference

| Signal | Number | Shortcut | Description |
|--------|--------|----------|-------------|
| SIGINT | 2 | Ctrl+C | Interrupt (stop) |
| SIGTSTP | 20 | Ctrl+Z | Suspend (pause) |
| SIGTERM | 15 | `kill PID` | Terminate gracefully |
| SIGKILL | 9 | `kill -9 PID` | Force kill |
