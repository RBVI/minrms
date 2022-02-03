Installation instructions for MINRMS
==========================

This file explains how to compile and install **minrms**
and it's companion software tools.


## Outline

1) Compile MINRMS
2) Windows and MacOS specific instructions


## Prerequisites:

Before you begin, you will need:

- A modern C++ compiler (C++11 or later).  This software has been tested using
  GCC (v9.3.0) and CLANG (v10.0.0).  Earlier compiler versions may not work.
- The "make" tool.
  (Type "make" into the terminal.  If it doesn't complain "command not found"
  then make is installed.)
- git is recommended.
  (Usually "git" can be installed using your OS' package management
  installer tool, such as "apt", "yum", or "brew", etc...)


## Download the MINRMS code

```
git clone https://github.com/jewettaij/minrms ~/minrms
```

## Compile the MINRMS code

```
cd ~/minrms
```

The next step depends on which compiler you are using.

### If you are using the GCC compiler, enter:
```
gcc --version
```
This will print the version of your compiler to the terminal.
If the version is earlier than 4.9, then stop.
Either update your compiler, or install CLANG.
If your compiler is 4.9 or later, then enter this:
```
source setup_gcc.sh
```
...and skip th the next section.

### If you are using the CLANG compiler, enter this instead:
```
source setup_clang.sh
```

## Now compile the MINRMS software using:

```
make clean
make
```

Copy the MINRMS executable programs to somewhere in your PATH.
In these instructions, I assume you have a ~/bin directory in your path.
If you do not have a writeable directory in your PATH,
(or if you have no idea what I'm talking about), then enter this beforehand:

```
mkdir ~/bin
echo "export PATH=\"$PATH:$HOME/bin\" >> ~/.bashrc
```

Now copy these files to that ~/bin directory:
```
cp -f bin/minrms/minrms ~/bin/
cp -f bin/msf2stat3d/msf2stat3d ~/bin/
cp -f bin/msf_compare/msf_compare ~/bin/
cp -f bin/msf2sequence/msf2sequence ~/bin/
cp -f bin/pdb2sequence/pdb2sequence ~/bin/
cp -f bin/pdb_select/pdb_select ~/bin/
cp -f bin/rotate_pdb/rotate_pdb ~/bin/
cp -f bin/align2msf/align2msf.py ~/bin/
```


## Apple MacOS instructions

If you are using an apple computer, you will need to install a C++11
compatible compiler.  As of 2022-2-02, this is included with Xcode.
It can also be installed using the Homebrew (brew) package management system.


## Windows instructions

If you have a preferred way of compiling C++ code and entering commands
into a terminal, you can probably ignore this section.  Otherwise, read on...

Although it is optional, I recommend installing 
the BASH shell environment on your computer,
along with *make* and either *gcc* or *clang*.
There are several ways to to create a BASH environment in Windows.
If you are using Windows 10 or 11, perhaps the easiest method is to install
[Windows Subsystem for Linux (WSL/WSL2)](https://devblogs.microsoft.com/commandline/a-preview-of-wsl-in-the-microsoft-store-is-now-available/)
*If you have an older version of windows, you can try
[virtualbox](https://www.virtualbox.org) or [CYGWIN](https://www.cygwin.com/).
(In the former case, you will also need to install a linux distribution,
preferably with a lightweight
desktop such as [xubuntu](https://xubuntu.org).)*

### Background
WSL and virtualbox are virtual machines that allow you to run an
alternate operating system from within windows.
In this case that operating system is linux.  The BASH shell and the
compiler tools that you need are usually included with linux.
If not, they can be easily installed from within in linux.
Both WSL and virtualbox also create an alternate filesystem inside windows
where the linux operating system is stored.  Software (like *minrms*)
that you download and install there can access the files in that filesystem.
So you may need to copy your PDB files and other files to this fileystem
beforehand.


*These instructions were written on 2022-2-02.*
