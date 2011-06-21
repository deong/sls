#!/usr/bin/perl -w

use strict;
use warnings;

my $filehnd;

while(<>) {
  if(/trial (\d*)/) 
  {
    open($filehnd, ">", "trial$1.dat");
  }
  else {
    print $filehnd $_;
    #print $filehnd "\n";
  }
}
