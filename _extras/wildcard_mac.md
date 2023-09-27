---
layout: page
title: Wildcard on Mac
---
## MacOSX Wildcard use on zsh
The shell (both bash & zsh) tries to interpret `scp abc@123:/home/se/exports/201405091107/*` as a glob to match files on your local system. The shell doesn't know what scp is, or that you're trying to match remote files.

The difference between bash and zsh is their default behavior when it comes to failed globbing. In bash, if a glob doesn't match anything, it passes the original glob pattern as an argument. In zsh it throws an error instead.

To address the issue, you need to quote it so the shell doesn't try to interpret it as a local glob.
```bash
scp 'abc@123:/home/se/exports/201405091107/*' .
(other things like ...1107/'*' or ...1107/\* work too)
```

If you want to change it so the zsh no-match behavior is the same as bash, you can do the following
```bash
setopt nonomatch
```

If you need to use is everytime, please follow this
```bash
echo "setopt nonomatch" >> ~/.zshrc
source ~/.zshrc
```

