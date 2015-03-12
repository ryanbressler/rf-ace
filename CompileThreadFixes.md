Here is a list of the stable compilations that have been achieved with v1.2.5, branch thread-fixes:

# Mac OS X #

On Mac OS X 10.7 Lion, install [Homebrew](http://mxcl.github.io/homebrew/) as described in the instructions.

From the Mac OS X terminal, install gcc 4.7:
```
brew tap homebrew/versions
brew install gcc-4.7
```

or gcc 4.8 instead, the second command should be:
```
brew install gcc-4.8
```

Follow the instructions in the section [Checkout and Compile RF-ACE](CompileThreadFixes#Checkout_and_Compile_RF-ACE.md)

# Linux (Debian) #

To install gcc-4.7, you can search for and install the package:

```
sudo aptitude search gcc-4.7
sudo apt-get install gcc-4.7
```

at this point, gcc-4.7 should be available.  Test whether the default gcc is pointing to this version
```
gcc --version
```

If it is not pointing at this version, call g++-4.7 directly by making the same change to the Makefile as described in the [Checkout and Compile RF-ACE](CompileThreadFixes#Checkout_and_Compile_RF-ACE.md) section.

# Linux (Red Hat/CentOS) #

To install gcc-4.7, you can search for and install the package:

```
sudo yum search gcc-4.7
sudo yum install gcc-4.7
```

at this point, gcc-4.7 should be available.  Test whether the default gcc is pointing to this version
```
gcc --version
```

If it is not pointing at this version, call g++-4.7 directly by making the same change to the Makefile as described in the [Checkout and Compile RF-ACE](CompileThreadFixes#Checkout_and_Compile_RF-ACE.md) section.

# Checkout and Compile RF-ACE #

Grab the `thread-fixes` branch of rf-ace by cloning to a local directory:
```
svn checkout https://rf-ace.googlecode.com/svn/branches/thread-fixes rf-ace
```

To compile rf-ace successfully, edit the Makefile in the root directory of the project.  Change the first line to:
```
COMPILER = g++-4.7
```

or g++-4.8 is you've installed that.

no you may run
```
make
make test
```

from the terminal.  Tests on Mac OS X 10.7 and Linux Mint 14 using these instructions have been all successful, with the output:
```
ALL DONE! 4887 tests run: 4887 successes and 0 failures
```

If the project fails to compile, try to compile using:
```
make no-threads
```

For questions about these specific instructions, you can email `codefor@systemsbiology.org`