#!/usr/bin/gnuplot

set terminal postscript eps color
set output 'progress.eps'
plot 'spea2ts_landscape/spea2ts_landscape.dat' using 1:2 title 'SPEA2+TS Initial' with points,\
     'spea2ts_landscape/spea2ts_landscape.dat' using 3:4 title 'SPEA2+TS Final' with points,\
     'mots_landscape/mots_landscape.dat' using 1:2 title 'MoTS Initial' with points,\
     'mots_landscape/mots_landscape.dat' using 3:4 title 'MoTS Final' with points

