---
layout: page
title: Python Version Control
published: true
---
## Python Version Control
 
If you have ever watch my youtube which I would encourage you to watch and subscribe to, you must have noticed that most often that not I have two python environment that keep on conflicting and interfering with my work. The fact that python 2.7 and below will no longer be supported from next year means alot. Means that if you still want to collaborate with your team or on some old projects, then you need to learn on how to manage these two different development environments and.

In this tutorial I am gonna write about how to manage this challenge that often occurs if you have two different versions of python installed or how to install and manage them up. I will be using the pyenv package.

## Why Use pyenv?
pyenv is a wonderful tool for managing multiple Python versions. Even if you already have Python installed on your system, it is worth having pyenv installed so that you can easily try out new language features or help contribute to a project that is on a different version of Python.

### Installing pyenv
Before you install pyenv itself, you’re going to need some OS-specific dependencies. These dependencies are mostly development utilities written in C and are required because pyenvinstalls Python by building from source. For a more detailed breakdown and explanation of the build dependencies, you can check out the official docs. In this tutorial, you’ll see the most common ways to install these dependencies.

### Build Dependencies
pyenv builds Python from source, which means you’ll need build dependencies to actually use pyenv. The build dependencies vary by platform. If you are on Ubuntu/Debian and want to install the build dependencies, you could use the following:

```
 
$ sudo apt-get install -y make build-essential libssl-dev zlib1g-dev \
libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl
```
{: .bash}

 This uses Apt to install all the build dependencies. Let this run, and you’ll be ready to go for Debian systems.

 If you use Fedora/CentOS/RHEL, you could use yum to install your build dependencies:
```
$ sudo yum install zlib-devel bzip2 bzip2-devel readline-devel sqlite \
sqlite-devel openssl-devel xz xz-devel libffi-devel
```
{: .bash}
 This command will install all the build dependencies for Python using yum.

 macOS users can use the following command:

```
$ brew install readline xz
```
{: .bash}  


 This command relies on Homebrew and installs the few dependencies for macOS users.

 If you’re instead using openSUSE then you would run the following:
```
$ zypper in zlib-devel bzip2 libbz2-devel libffi-devel \
libopenssl-devel readline-devel sqlite3 sqlite3-devel xz xz-devel
Once again, this command installs all the Python build dependencies for your system.
```
{: .bash}  

#### Finally, for Alpine users, you can use this:
```
$ apk add libffi-dev ncurses-dev openssl-dev readline-dev \
tk-dev xz-dev zlib-dev
```
{: .bash}  

 This command uses apk as the package manager and will install all build dependencies for Python on Alpine.

Using the pyenv-installer
After you’ve installed the build dependencies, you’re ready to install pyenv itself. I recommend using the pyenv-installer project:
```
$ curl https://pyenv.run | bash
```
{: .bash}  
 This will install pyenv along with a few plugins that are useful:
```
pyenv: The actual pyenv application
pyenv-virtualenv: Plugin for pyenv and virtual environments
pyenv-update: Plugin for updating pyenv
pyenv-doctor: Plugin to verify that pyenv and build dependencies are installed
pyenv-which-ext: Plugin to automatically lookup system commands
Note: The above command is the same as downloading the pyenv-installer script and running it locally. So if you’d like to see exactly what you’re running, you can view the file yourself. Alternatively, if you really don’t want to run a script, you can checkout the manual installation instructions.
```

 At the end of the run, you should see something like this:

*WARNING: seems you still have not added 'pyenv' to the load path.*

#### Load pyenv automatically by adding 
#### the following to ~/.bashrc: 
```
export PATH="$HOME/.pyenv/bin:$PATH" 
eval "$(pyenv init -)" 
eval "$(pyenv virtualenv-init -)" 
```

#### The output will be based on your shell. But you should follow the instructions to add pyenvto your path and to initialize pyenv/pyenv-virtualenv auto completion. Once you’ve done this, you need to reload your shell:

```
$ exec "$SHELL" # Or just restart your terminal
```
{: .bash}  


 That’s it. You now have pyenv and four useful plugins installed.

Using pyenv to Install Python
Now that you have pyenv installed, installing Python is the next step. You have many versions of Python to choose from. If you wanted to see all the available CPython 3.6 through 3.8, you can do this:
```
$ pyenv install --list | grep " 3\.[678]"   
3.6.0   
3.6-dev   
3.6.1   
3.6.2   
3.6.3   
3.6.4   
3.6.5   
3.6.6   
3.6.7   
3.6.8   
3.7.0   
3.7-dev   
3.7.1   
3.7.2   
3.8-dev 
```
{: .bash}  
#### The above shows all the Python versions that pyenv knows about that match the regular expression. In this case, that is all available CPython versions 3.6 through 3.8. Likewise, if you wanted to see all the Jython versions, you could do this:
```
$ pyenv install --list | grep "jython"  
```
{: .bash}  
#### Again, you can see all the Jython versions that pyenv has to offer. If you want to see all the versions, you can do the following:
```
$ pyenv install --list
```
{: .bash}  
...
 There are a lot
 Once you find the version you want, you can install it with a single command:

```
$ pyenv install -v 3.7.2 
```
{: .bash}  

```
/tmp/python-build.20190208022403.30568 ~ 
Downloading Python-3.7.2.tar.xz... 
-> https://www.python.org/ftp/python/3.7.2/Python-3.7.2.tar.xz Installing Python-3.7.2... 
/tmp/python-build.20190208022403.30568/Python-3.7.2 /tmp/python-build.20190208022403.30568 ~ 
[...] 
Installing collected packages: setuptools, pip Successfully installed pip-18.1 setuptools-40.6.2 Installed Python-3.7.2 to /home/realpython/.pyenv/versions/3.7.2 
```
{: .output}
  
 This will take a while because pyenv is building Python from source, but once it’s done, you’ll have Python 3.7.2 available on your local machine. If you don’t want to see all the output, just remove the -v flag. Even development versions of CPython can be installed:
```
$ pyenv install 3.8-dev
```
{: .bash}  

 Pro Tip: If you’ve been using pyenv for a while and don’t see the version you’re looking for, you may need to run pyenv update to update the tool and make sure you have access to the latest versions.

 For the rest of the tutorial, the examples assume you’ve installed 3.6.8 and 2.7.15, but you’re free to substitute these values for the Python versions you actually installed. Also note that the system Python version in the examples is 2.7.12.

#### Installation Location
 As mentioned before, pyenv works by building Python from source. Each version that you have installed is located nicely in your pyenv root directory:
```
$ ls ~/.pyenv/versions/
```
{: .bash}  

```
2.7.15  3.6.8  3.8-dev
```
{: .output}
  
 All of your versions will be located here. This is handy because removing these versions is trivial:
```
$ rm -rf ~/.pyenv/versions/2.7.15
```
{: .bash}  

 Of course pyenv also provides a command to uninstall a particular Python version:
```
$ pyenv uninstall 2.7.15
```

 Using Your New Python
 Now that you’ve installed a couple of different Python versions, let’s see some basics on how to use them. First, check what versions of Python you have available:

```
$ pyenv versions 
```
{: .bash}  
```
* system (set by /home/realpython/.pyenv/version)   
2.7.15   
3.6.8   
3.8-dev 
```
{: .output}
  
 The * indicates that the system Python version is active currently. You’ll also notice that this is set by a file in your root pyenv directory. This means that, by default, you are still using your system Python:

```
$ python -V
```
{: .bash}  

```
Python 2.7.12
```
{: .output}
  
 If you try to confirm this using which, you’ll see this:

```
$ which python
```
{: .bash}  
```
/home/realpython/.pyenv/wyim/python
```
{: .output}  

 This might be surprising, but this is how pyenv works. pyenv inserts itself into your PATH and from your OS’s perspective is the executable that is getting called. If you want to see the actual path, you can run the following:

```
$ pyenv which python
```
{: .bash}  
```
/usr/bin/python
``` 
{: .output}  
 If, for example, you wanted to use version 2.7.15, then you can use the global command:

```
$ pyenv global 2.7.15 
$ python -V 
```
{: .bash}  
```
Python 2.7.15 
```
{: .output}  

```
$ pyenv versions   
```
{: .bash}  
```
system 
* 2.7.15 (set by /home/realpython/.pyenv/version)   
3.6.8   
3.8-dev 
```
{: .output}  
 Pro Tip: A great way to get peace of mind that the version of Python you just installed is working properly is to run the built-in test suite:

$ pyenv global 3.8-dev
$ python -m test
 This will kick off lots of internal Python tests that will verify your installation. You can just kick back and watch the tests pass.

 If you ever want to go back to the system version of Python as the default, you can run this:

```
$ pyenv global system
$ python -V
```
{: .bash}

```
Python 2.7.12
```
{: .output}


 You can now switch between different versions of Python with ease. This is just the beginning. If you have many versions that you want to switch between, typing these commands consistently is tedious.

### Virtual Environments and pyenv
#### Virtual environments are a big part of managing Python installations and applications.

Virtual environments and pyenv are a match made in heaven. pyenv has a wonderful plugin called pyenv-virtualenv that makes working with multiple Python version and multiple virtual environments a breeze. If you’re wondering what the difference is between pyenv, pyenv-virtualenv, and tools like virtualenv or venv, then don’t worry. You’re not alone.

*Here’s what you need to know:*

pyenv manages multiple versions of Python itself.
virtualenv/venv manages virtual environments for a specific Python version.
pyenv-virtualenv manages virtual environments for across varying versions of Python.
If you’re a die-hard virtualenv or venv user, don’t worry: pyenv plays nicely with either. In fact, you can keep the same workflow you’ve had if you’d prefer, though I think pyenv-virtualenv makes for a nicer experience when you’re switching between multiple environments that require different Python versions.

 The good news is that since you used the pyenv-installer script to install pyenv, you already have pyenv-virtualenv installed and ready to go.

Creating Virtual Environments
Creating a virtual environment is a single command:
```
$ pyenv virtualenv <python_version> <environment_name>
```
{: .bash}
Technically, the <python_version> is optional, but you should consider always specifying it so that you’re certain of what Python version you’re using.

The <environment_name> is just a name for you to help keep your environments separate. A good practice is to name your environments the same name as your project. For example, if you were working on myproject and wanted to develop against Python 3.6.8, you would run this:

```
$ pyenv virtualenv 3.6.8 myproject
```
{: .bash}
The output includes messages that show a couple of extra Python packages getting installed, namely wheel, pip, and setuptools. This is strictly for convenience and just sets up a more full featured environment for each of your virtual environments.

Activating Your Versions
Now that you’ve created your virtual environment, using it is the next step. Normally, you should activate your environments by running the following:

```
$ pyenv local myproject
```
{: .bash}

You’ve seen the pyenv local command before, but this time, instead of specifying a Python version, you specify an environment. This creates a .python-version file in your current working directory and because you ran eval "$(pyenv virtualenv-init -)" in your environment, the environment will automatically be activated.

*You can verify this by running the following:*

```
$ pyenv which python
```
{: .bash}

```
/home/realpython/.pyenv/versions/myproject/bin/python
```
{: .output}
You can see a new version has been created called myproject and the python executable is pointing to that version. If you look at any executable this environment provides, you will see the same thing. Take, for example, pip:

```
$ pyenv which pip 
```
{: .bash}

```
/home/realpython/.pyenv/versions/myproject/bin/pip 
```
{: .output}
If you did not configure eval "$(pyenv virtualenv-init -)" to run in your shell, you can manually activate/deactivate your Python versions with this:

```
$ pyenv activate <environment_name> 
$ pyenv deactivate 
```
{: .bash}

The above is what pyenv-virtualenv is doing when it enters or exits a directory with a .python-version file in it.

Working With Multiple Environments
Putting everything you’ve learned together, you can work effectively with multiple environments. Let’s assume you have the following versions installed:
```
$ pyenv versions * 
```
{: .bash}
```
system (set by /home/realpython/.pyenv/version)   
2.7.15   
3.6.8   
3.8-dev 
```
{: .output}

Now you want to work on two different, aptly named, projects:

project1 supports Python 2.7 and 3.6.
project2 supports Python 3.6 and experiments with 3.8-dev.
You can see that, by default, you are using the system Python, which is indicated by the * in the pyenv versions output. First, create a virtual environment for the first project:
```
$ cd project1/ 
$ pyenv which python /usr/bin/python 
$ pyenv virtualenv 3.6.8 project1 ... 
$ pyenv local project1 
$ python -V 
```
{: .bash}
```
/home/realpython/.pyenv/versions/project1/bin/python 
```
{: .output}
Finally, notice that when you cd out of the directory, you default back to the system Python:
```
$ cd $HOME 
$ pyenv which python /usr/bin/python 
```
{: .output}
You can follow the above steps and create a virtual environment for project2:
```
$ cd project2/ 
$ pyenv which python /usr/bin/python 
$ pyenv virtualenv 3.8-dev project2 ... 
$ pyenv local 3.8-dev 
$ pyenv which python 
```
{: .output}

/home/realpython/.pyenv/versions/3.8-dev/bin/python 
These are one time steps for your projects. Now, as you cd between the projects, your environments will automatically activate:
```
$ cd project2/ 
$ python -V Python 3.8.0a0 
$ cd ../project1 
$ python -V 
```
{: .bash}
```
Python 3.6.8 
```
{: .output}
No more remembering to activate environments: you can switch between all your projects, and pyenv will take care of automatically activating the correct Python versions and the correct virtual environments.

Activating Multiple Versions Simultaneously
As described in the example above, project2 uses experimental features in 3.8. Suppose you wanted to ensure that your code still works on Python 3.6. If you try running python3.6, you’ll get this:
```
$ cd project2/ 
$ python3.6 -V 
```
{: .bash}
```
pyenv: python3.6: command not found 

The `python3.6' command exists in these Python versions:   
3.6.8   
3.6.8/envs/project1   
project1 
```
{: .output}
pyenv informs you that, while Python 3.6 is not available in the current active environment, it is available in other environments. pyenv gives you a way to activate multiple environments at once using a familiar command:
```
$ pyenv local project2 3.6.8
```
{: .bash}

This indicates to pyenv that you would like to use the virtual environment project2 as the first option. So if a command, for example python, can be resolved in both environments, it will pick project2 before 3.6.8. Let’s see what happens if you run this:
```
$ python3.6 -V
```
{: .bash}
```
Python 3.6.8
```
{: .output}

Here, pyenv attempts to find the python3.6 command, and because it finds it in an environment that is active, it allows the command to execute. This is extremely useful for tools like tox that require multiple versions of Python to be available on your PATH in order to execute.P

Suppose that in the above example, you’ve found a compatibility problem with your library and would like to do some local testing. The testing requires that you install all the dependencies. You should follow the steps to create a new environment:
```
$ pyenv virtualenv 3.6.8 project2-tmp 
$ pyenv local project2-tmp 
```
{: .bash}
Once you’re satisfied with your local testing, you can easily switch back to your default environment:
```
$ pyenv local project2 3.6.8  
```
{: .bash}

Thank you for reading. let me know on the comment section how helpful it was.

