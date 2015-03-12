## Q:Where is Windows support? ##
Since the addition of multithreading, I have been having difficulties compiling the source under Windows. This is because Visual Studio compiler doesn't recognize the `-pthread` option. The lack of support to multithreading might sometimes also arise under Linux, for which there is already a simple solution (check the next question!)

## Q:How do I compile without multithreading? ##
Call
```
make no-threads
```
This will make `-e / --nThreads` obsolete.

## Q:Is it possible to use RF-ACE from R? ##
Yes, check [this](http://code.google.com/p/rf-ace/wiki/RPackageInstallation) out.

## Q: How about other scripting environments? ##
Not yet... stay tuned!