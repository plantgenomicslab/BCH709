---
layout: page
title: Fix your warning
published: true
---
## `$TERM` warning
```
"tput: No value for $TERM and no -T specified:
```
### Solution
Please add `tty -s &&` in front of your `export PS1`


```bash
tty -s && export PS1="\[\033[38;5;2m\]\u\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;33m\]@\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;166m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;2m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;4m\]\\$\[$(tput sgr0)\]\[\033[38;5;15m\]\n \[$(tput sgr0)\]"
```


## Locale Warning
```bash
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
    LANGUAGE = (unset),
    LC_ALL = (unset),
    LANG = "en_US.UTF-8"
are supported and installed on your system.
perl: warning: Falling back to the standard locale ("C").
```
### Please add this to your `.bashrc` or `bash_profile`  
```bash
   export LANG=en_US.UTF-8
```
