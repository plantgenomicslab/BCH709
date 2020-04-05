---
layout: page
title: Secure Shell 
published: true
---
** Introduction

The tool [[https://en.wikipedia.org/wiki/Secure_Shell][SSH]] (Secure Shell) can be used to connect to another computer
and execute programs there. The remote machine is often a server that
for example offer more computational power than you laptop. The most
widely used implemenation of SSH is [[http://www.openssh.com/][OpenSSH]].

> ## Basic usage
> In order to connect to the remote machine you need to have an account
on it. In the following examples we assume that you have an account
called =alice= and a remote machine called =bobs_server= (Your name is
not Alice and you don't know any Bob? [[https://en.wikipedia.org/wiki/Alice_and_Bob][Read this]] to understand the
choice of names).
> You can connect to a remote machine and have an interactive session by
calling ssh like this (don't copy the =$= sign in any of the following
commands!):

```bash
$ ssh <your_username>@<server_name_or_ip>
```

In our example this would be

```
$ ssh alice@bobs_server
```

SSH will presend you the host finger print

```
Host key fingerprint is SHA256:ZwDBE5pV7SCqTXrHYzBChkhn3Nkx6cvvF9StX+W25No
+---[ECDSA 256]---+
|oo.o..B*+.       |
|o =. *++o .      |
| o  o..o.o  . .  |
|  . =  . ... . ..|
|   * +. S +   ...|
|  o o =o o . . .+|
|   . o ..   . +.o|
|         . .  .+ |
|        ...  ..E |
+----[SHA256]-----+
```

The first time your are connting to a remove machine you will be asked
to continue. Please type =yes= and press Enter.

```
Are you sure you want to continue connecting (yes/no)?
Warning: Permanently added 'hawaii' (ECDSA) to the list of known hosts
```

Then you will be asked for your password.

```
alice@bobs_server's password:
```

Please enter your password and press Enter. Now you should have an
interactive shell session on the remote machine. You prompt might
indicate this but you can use =hostname= to make that you are on the
remote machine.

```
$ hostname
bobs_server
```

You can disconnect by typing 
```
$ exit
```

You will be informed that hte connection is closed
```
Connection to bobs_server closed.
```

Now you can run commandline tools like on your local machine. If the
server is adequately configured you can even run tools with a
graphical interface (GUI) using so called X11-forwarding. For this
connect with ssh using the -X flag (disconnect first if you are
connected).

```
$ ssh -X alice@bobs_server
```

You can also run SSH non-interactively by calling it with a command to
execute:

```
$ ssh <your_username>@<server_name_or_ip> <command>
```

Example:

```
$ ssh alice@bobs_server ls
```

This will list the file in the home directory of the remote
machine. You will not have an interactive session on that machine but
will stay in the session of you local machine after the file list is
returned.

** Copying files with scp

Basic

```
$ scp <source> <target>
```

Folder

```
$ scp -r <source> <target>
```

Using compression

```
$ scp -C <source> <target>
```

Using compression when sending

```
$ scp -rC <source> <target>
```

Examples. Be aware that =~/= represent your home directory on
GNU/Linux and the different BSD this is usually
=/home/your_username=. On OS X this is =/Users/your_username=. 

This will copy the file =my_local_big_file.txt= from your local
machine into your home directory of the remote machine.

```
$ scp my_local_big_file.txt alice@bobs_server:~/
```

The following command will send the file to the server and rename it
to =my_remote_big_file.txt=.

```
$ scp my_local_big_file.txt alice@bobs_server:~/my_remote_big_file.txt
```

Sure you can use scp to copy from the server to you local machine:

```
$ scp alice@bobs_server:/data/my_local_big_file.txt .
```

** Configure SSH

While you are now capable of working on the remote machine you have to
type a lot to connect to it. /Laziness is virtue/ - in the following
we will make working with SSH more convinient and reduce the need of
typing.

Create in the home directory of you local machine a folder called
=.ssh=.

```
$ mkdir -p ~/.ssh
```

Make it unaccessible (read/write/executable) by other users:

```
$ chmod -R go-rwx ~/.ssh
```

Create a configuration file

```
$ touch ~/.ssh/config
```

Open this configuration file in you favorite text editor (e.g. Emacs,
vi, nano, gedit) and add

```
Host bobs
  User alice
  HostName bobs_server
```

Now you can connect to the server simply by calling

```
$ ssh bobs
```

or copy files without by using the short alias

```
$ scp my_file.txt bobs:~/
```

** Using SSH keys

Currently you have to type your password everytime you connect to the
remote machine. In order to avoid that we will generate once a so
called /SSH key pair/. A private key will be stored at your local
machine and should be *keep private*! The public key has to be put on
the server and can be liberately given away (e.g. a system
administrators might ask you for the public key in order to give you
access to a remote machine without any password). More details
regarding this public/private key approach can be found [[https://en.wikipedia.org/wiki/Public-key_cryptography][here]].

Generate a key pair on your local machine. We use =-b 4096 in order
to generate a stronger key than with the default settings and =-N ""=
to set no passphrase for the key. Generating the key can take a while.

```
$ ssh-keygen -b 4096 -N ""
```

During the procedure you will be asked for the name of the key
file. Use the default by pressing Enter:

```
Generating public/private rsa key pair.
Enter file in which to save the key (/home/your_local_user_name/.ssh/id_rsa):
Your identification has been saved in /home/your_local_user_name/.ssh/id_rsa.
Your public key has been saved in /home/your_local_user_name/.ssh/id_rsa.pub.
The key fingerprint is:
0b:5c:b3:e8:51:d2:0e:3b:b8:cc:6d:97:4d:1d:e7:f0 your_local_user_name@your_local_machine
The key's randomart image is:
+--[ RSA 4096]----+
|                 |
|       .         |
|      o =   o .  |
|     o O o . *   |
|    . B S . . E  |
|   o + + =       |
|    + + + .      |
|     . .         |
|                 |
+-----------------+
```

As =ssh-keygen= reported we have now two files:

- Your private key file =/home/your_local_user_name/.ssh/id_rsa=
- Your public key file =/home/your_local_user_name/.ssh/id_rsa.pub=

Those are just plain test files and you can have a look at them using
=cat= or =less=:

```
$ cat ~/.ssh/id_rsa
$ cat ~/.ssh/id_rsa.pub
```

We now want to copy the public key to the server. To set up a =.ssh=
folder on the server run the following command on your local machine:

```
$ ssh alice@bobs_server "mkdir -p ~/.ssh"
```

Add the public key to the file (as you can see you can even pipe
commands into SSH) by running the following command on your local
machine:

```
$ cat ~/.ssh/*.pub | ssh alice@bobs_server 'cat >> ~/.ssh/authorized_keys'
$ cat ~/.ssh/*.pub | ssh alice@bobs_server 'chmod 700 ~/.ssh'
$ cat ~/.ssh/*.pub | ssh alice@bobs_server 'chmod 600 ~/.ssh/authorized_keys'
```

Add the private key to local SSH key-chain:

```
$ ssh-add ~/.ssh/id_rsa
```

** Persistent session with tmux

Alternatives:

-  nohub
-  screen

The terminal multiplexer [[https://tmux.github.io/][tmux]]

Start tmux for the first time

```
$ tmux
```

List sessions:

```
$ tmux list-sessions
```

Reconnect to running tmux sessions:

```
$ tmux attach
```

or simply

```
$ tmux a
```

Important keys:

- =Ctr-b c= - new window
- =Ctr-b ,= - rename window
- =Ctr-b x= - close window
- =Ctr-b <NUMBER>= - go to window with the given number 
- =Ctr-b [= - switch scrolling mode; use =q= to quit

List keys

```
$ tmux list-keys
```

List commands

```
$ tmux list-commands
```

Configure =~/.tmux.conf=:

```
# Change prefix to Control + o
set -g prefix C-o

# Increase history length
set -g history-limit 1000
```

** Further topics

- [[https://mosh.mit.edu/][Mosh]] (mobile shell) - more robust and responsive alternative to SSH.
- [[https://en.wikipedia.org/wiki/SSHFS][SSHFS]]
  - Under OS X:
    - https://osxfuse.github.io/
    - https://github.com/osxfuse/osxfuse/wiki/SSHFS

```
mkdir -p server_folder sshf alice@bobs_server:/home/alice/
umount server_folder
```

This work was adapted from [Konrad Förstner](https://github.com/konrad)
