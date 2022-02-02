Installation instructions for MINRMS
==========================

This file explains how to install **minrms**.

## Outline

1) Install MINRMS
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

## STEP 1: Install the MINRMS software tools

### Download the MINRMS code

```
git clone https://github.com/jewettaij/minrms ~/minrms
```

### Compile the MINRMS code

```
cd ~/minrms
```

The next step depends on which compiler you are using.

#### If you are using the GCC compiler, enter:
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

#### If you are using the CLANG compiler, enter this instead:
```
source setup_clang.sh
```

### Now compile the MINRMS software using:

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
cp -f bin/filter_mrc/filter_mrc ~/bin/
cp -f bin/sum_voxels/sum_voxels ~/bin/
cp -r bin/combine_mrc/combine_mrc ~/bin/
cp -r bin/voxelize_mesh/voxelize_mesh.py ~/bin/
```

Then log out and log in again for the changes to take effect.



## STEP 2: Install the python modules that MINRMS needs:

#### Optional, but recommended:

First setup a virtual environment where you can run the MINRMS tools.
```
python3 -m venv ~/bin/venv_minrms    #(use "python" if "python3" is not found)
deactivate # Don't worry if this step prints an error message.
source ~/bin/venv_minrms/bin/activate
```

Then run these commands:
*(Note: If it complains "pip3: command not found",
then try using "pip" instead.)*

```
pip3 install numpy       # (Use "pip" if "pip3" fails.)
pip3 install mrcfile
pip3 install pyvista
pip3 install matplotlib  # optional but recommended
pip3 install skimage     # optional but recommended
pip3 install scipy       # optional but recommended
```




Installation is complete.



### *If you created a virtual environment*

If you followed the suggestion above and created a virtual environment
*(using "python -m venv ~/bin/venv_minrms")*, then you must do this
before using the MINRMS tools:  Start a new terminal and enter:
```
source ~/bin/venv_minrms/bin/activate
```
You must do this every time before you run the MINRMS software.
(Otherwise, the "*voxelize_mesh.py*" program will not work.)
Alternatively, you can add that command to your ~/.bashrc file this way:
```
echo "source ~/bin/venv_minrms/bin/activate" >> ~/.bashrc
```


## Apple MacOS instructions
If you are using an apple computer, you will need to install a C++11
compatible compiler.

## Windows instructions

It is recommended that you install the BASH shell environment on your computer,
along with *make* and either *gcc* or *clang*.  Once you have done that,
you can follow the instructions above for linux users.
There are several ways to to create a BASH environment,
but perhaps the easiest method is to install
[Windows Subsystem for Linux (WSL/WSL2)](https://devblogs.microsoft.com/commandline/a-preview-of-wsl-in-the-microsoft-store-is-now-available/)
***or***
[virtualbox](https://www.virtualbox.org)
(In the later case, you will also need to install a linux distribution,
preferably with a lightweight
desktop such as [xubuntu](https://xubuntu.org).)
Alternatively, you can try 
[Hyper-V](https://www.nakivo.com/blog/run-linux-hyper-v/)
or (if you have an older version of windows)
[CYGWIN](https://www.cygwin.com/).

WSL and virtualbox are virtual machines that allow you to run an
alternate operating system from within windows.
In this case that operating system is linux.  The BASH shell and the
compiler tools that you need can be easily installed from within in linux.
Both WSL and virtualbox also create an alternate filesystem inside windows
where the linux operating system is stored.  Software (like *minrms*)
that you download and install there can access the files in that filesystem.
So you may need to copy your PDB files and other files to this fileystem
beforehand.
