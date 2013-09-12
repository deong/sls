Overview
========

This code formed the experimental basis for my dissertation work on multiobjective
landscape analysis. It is a fairly general platform for implementing multiobjective
metaheuristics for numeric or combinatorial optimization. The focus of the system is
on making it fairly straightforward to describe an algorithm, and this comes
somewhat at the expense of the clarity of the internal code. I'm not too happy with
the system as it stands today, as I think the complexity has outpaced the benefits.
I'm currently working on a ground-up rewrite of the package that will remove much of
the incidental complexity, but as of now, the newer package is not yet ready to use
(and as of mid-2011, I'm no longer devoting much time to this line of research, so
the future of that code is unknown).

The code is all C++, and relies rather heavily on features of modern C++. Most of
the code gets pulled into a single giant template instantiation, which has a number
of bad consequences. The most obvious is that compile times are quite long, and very
little can be done in the way of separate or incremental compilation.

If you want to use this code, just run make from the src directory (you may need to
edit the makefile for your platform). Note that this will take several minutes to
complete. If it builds successfully (it's been tested on Mac, Linux, and Windows,
but not really recently and the Windows port in particular will probably need some
work to set up a Visual Studio project or something), then you should have a
command-line program named "sls". Run it by passing it one or more configuration
files. It will produce output in a roughly human-readable form on standard output.
There are some scripts in the src/scripts directory to do interesting things with
this output, but the scripts will require a functional Unix toolkit.


Getting the code
================
Note that it uses a submodule for the kvparse library, so the checkout procedure is
a little more involved than normal.

    $ git clone https://github.com/deong/sls
    $ cd sls
    $ git submodule init
    $ git submodule update
	
If you use the Github client for Windows (and probably other OSs as well) it
should fetch the submodule for you.


Building on Linux or OS X
=========================
I assume you're capable of installing boost here, but if you don't know where to
start, it should be in your package manager for Linux. On the Mac, it can be
installed via Homebrew (http://brew.sh). If you have put the boost headers and
libraries into a place where the default compiler can find them, then you can
just do the following.

    $ cd src
    $ make

This should, if you have a working compiler with libboost_regex somewhere on the
default library path, give you an executable named "sls". You can test it on some of
the included configuration files.

    $ ./sls ../cfg/nug12.cfg
    $ ./sls ../cfg/rana.cfg


Building on Windows with Visual Studio
======================================
There are project files included for both Visual Studio 2010 and 2012. With
luck, you can just open those projects and click build. I have no idea if that
actually works though, so you may end up needing to set up your own project.
Doing that is a little bit tedious due to the way templates work in C++. Most of
the classes in the system are template classes, and as such you can't really do
normal separate compilation. The way the code is set up is that the .h file for
one of these classes #includes the corresponding .cpp file. This means that if
you add those .cpp files to the project and compile them along with everything
else, you get duplicate symbols. To get it working properly, you need to get
boost_regex installed, and then create a project that excludes the right .cpp
files from the build.

1. Install boost
	- Download the current .7z file and unzip it somewhere

    - Move the boost_1_54_0 directory to somewhere convenient (optional)

2. Start a new empty project
	I named mine sls_vs2k12 and added the project directory underneath the
	top level sls directory, like
    
	    sls  
		    cfg  
    		prob  
    		sls_vs2k12  
    		src  
    		.gitignore  
    		.gitmodules  
    		README.md  
    
3. Add all the header files from the src directory to the project.

4. Add all of the cpp files from the src directory to the project EXCEPT for

    	enumerate.cpp  
    	gapgen.cpp  
    	ksgen.cpp  
    	mergepf.cpp  
    	pearson.cpp  

5. Add the following header files from the src/kvparse directory to the project

	    kvparse.h  
    	kvparse_except.h  

6. Add kvparse.cpp from the src/kvparse directory to the project

7. Select all of the cpp files under "Source Files" EXCEPT for the following

    	factory.cpp  
    	keywords.cpp  
    	kvparse.cpp  
    	mtrandom.cpp  
    	problems.cpp  
    	slsmain.cpp  
    	strategy.cpp  
    	utilities.cpp  

   and then right click on the selected files, choose "Properties", and under
   "Excluded from Build", select "Yes".

8. Right click onthe project in the Solution explorer (not the solution) and 
   select Properties. Make the following changes.
    - Under C/C++ -> All Options
	    * add your boost_1_54_0 directory to Additional Include Directories

	    * add "/bigobj" (without the quotes) to Additional Options
    - Under Linker -> All Options
	    * add your boost_1_54_0/stage/lib directory to Additional Library Directories

9. By default, this all happened for the Debug target. Change the target to Release
   and repeat steps 7 and 8.




Dependencies
============
* boost (specifically boost_regex)

* gtest (https://code.google.com/p/googletest/) -- needed only if you
  want to run the tests for the kvparse library. Note that if you do run the
  test, a few of them currently fail. The core functionality works, but there
  are some enhancements I haven't made yet so the tests for those fail.
