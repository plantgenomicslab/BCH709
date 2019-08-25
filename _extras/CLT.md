---
layout: page
title: Mac Command Line Tool Installation
published: true
---
## Mac Command Line Tool Installation

Whether you use bioinformatics script that does everything for you, or set everything up manually, it's best that you start with a clean installation of macOS. If you've already tried to install a development environment, I can't guarantee that you won't run into any issues. 

Check your macOS version below to get started:

- Mavericks and above, including the latest Mojave
- Mountain Lion
- Lion
- Snow Leopard

## Historical Background
Up until February 2012, the only way you could get the Command Line Tools required for web development was via the full Xcode package, which is almost 2 GB in size. Since then, Apple started offering the Command Line Tools (CLT) as a separate, much smaller download (~118MB), which benefits those who don't plan on writing Mac or iOS apps.

There is also a third-party option, the osx-gcc-installer by Kenneth Reitz, that supports both Snow Leopard and Lion. However, it is not updated as often as the official Apple tools, and I personally ran into issues with it on Lion.

###The Easy Way for Mavericks and above
I've written an open source script that can set everything up for you, including configuring your Mac to work with GitHub.

If you prefer to do everything manually, start with Step 1 below.

####Step 1: Download and Install the Command Line Tools
Installing the standalone Command Line Tools on Mavericks and above
Most of the work you'll be doing in this tutorial will be in the "Terminal" application. The easiest way to open an application in OS X is to search for it via Spotlight. The default keyboard shortcut for invoking Spotlight is command-Space. Once Spotlight is up, just start typing the first few letters of the app you are looking for, and once it appears, select it, and press return to launch it. See the animated GIF below for an example:

Launch Terminal via Spotlight

Inside the Terminal window, copy and paste (or type) the following command, and press the return key on your keyboard:
```
xcode-select --install
```
{: .bash}

You should see the pop up below on your screen. Click Install when it appears.

- install xcode on mavericks step 1

Click Agree when the License Agreement appears:

- install xcode on mavericks step 2

Your computer will then attempt to find the software, and then will start downloading it. The following popup will appear:

- install xcode on mavericks step 4

Once the software is installed, click Done. That's it! You're now ready to go to Step 2.

- install xcode on mavericks step 5

####Installing the standalone Command Line Tools on Mountain Lion
Go to [https://developer.apple.com/downloads](https://developer.apple.com/downloads) and sign in with your Apple ID (the same one you use for iTunes and app purchases).

sign in to developer.apple.com

Search for "command line tools" (in the search field on the left), then click on the latest version of "Command Line Tools (OS X Mountain Lion) for Xcode," and click on the the .dmg link to download it.

command line tools for mountain lion

Once the .dmg has finished downloading, double-click on it (if it didn't already open automatically). This will mount the disk image and open a window in your Finder that looks like this:

command line tools package installer for mountain lion

Double-click on the "Command Line Tools (Mountain Lion).mpkg" installer and go through the installation. Once the CLT are installed, launch the "Terminal" application via Spotlight (as explained in Step 1), then go to Step 2.

####Installing Xcode on Lion
Click on this link to Xcode on the Mac App Store, then click on "View in Mac App Store."

view in mac app store

It should automatically launch the "App Store" app on your Mac and take you the Xcode page. Click on the "Free" button, then click on "Install App."

Once the installation is complete, go to your Applications folder and double-click on Xcode, then install any required components if asked to.

install mobile component

Go to Xcode's Preferences via the menu bar, or by pressing the command and comma keys.

Go to Xcode Preferences

Click on the "Downloads" icon, then click on the "Install" button next to "Command Line Tools."

Install Command Line Tools

When prompted to log in, you should be able to use the same email and password you use for iTunes and app purchases. Once the Command Line Tools are installed, quit Xcode, launch the "Terminal" application via Spotlight (as explained in Step 1), then go to Step 2.

IMPORTANT NOTE: If you upgraded to Mountain Lion from Lion, and you already had Xcode installed on Lion, and you updated to Xcode 4.4 and updated the Command Line Tools while still on Lion, you will have to go back to Xcode and download the Command Line Tools again after upgrading to Mountain Lion.

####Snow Leopard Instructions
UPDATE: A kind reader (P. Martin) pointed out that the Xcode 4.2 download for Snow Leopard is only available to those registered in the $99/year developer program. I confirmed that the latest version of Xcode for Snow Leopard available to me while signed in with a free account is 3.2.6. I have not tested this setup with Xcode 3.2.6, but I would love to hear from you if you have. Otherwise, I recommend that you upgrade to a newer version of OS X.

Go to [https://developer.apple.com/downloads](https://developer.apple.com/downloads) and sign in with your Apple ID (the same one you use for iTunes and app purchases).

If you are part of the $99/year Apple developer program, search for "xcode 4.2" (in the search field on the left), then click on "Xcode 4.2 for Snow Leopard," and click on the .dmg link to download it.

Otherwise, search for "xcode 3.2", then click on "Xcode 3.2.6 and iOS SDK 4.3 for Snow Leopard," and click on the .dmg link to download it. As mentioned at the beginning of this section, I have not tested this tutorial with Xcode 3.2.6, so I would recommend that you upgrade to a newer version of OS X.

Download Xcode 4.2 for Snow Leopard

Once the .dmg has finished downloading, it should automatically mount the disk image and open a window in your Finder that looks like this:

Xcode package installer

Double-click on the "Xcode" package installer. Once the installer launches, make sure all the checkboxes are checked, as shown in the screenshot below:

Install Xcode

Click "Continue," and go through the rest of the installation. If the installation fails, quit the installer, then run Software Update and install any updates that it finds.

run software update

If no new updates are available, restart your computer and try installing Xcode again. Once Xcode is successfully installed, you can move on to Step 2.
