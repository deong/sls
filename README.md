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


Usage
============
Note that it uses a submodule for the kvparse library, so the checkout procedure is
a little more involved than normal.

    $ git clone https://github.com/deong/sls
    $ cd sls
    $ git submodule init
    $ git submodule update
    $ cd src
    $ make

This should, if you have a working compiler with libboost_regex somewhere on the
default library path, give you an executable named "sls". You can test it on some of
the included configuration files.

    $ ./sls ../cfg/nug12.cfg
    $ ./sls ../cfg/rana.cfg

Most of the keywords and values used in the config files are listed in keywords.cpp,
but unfortunately they aren't really documented well at this time.


Dependencies
============
* boost (specifically boost_regex)

* gtest (https://code.google.com/p/googletest/) -- needed only if you
  want to run the tests for the kvparse library