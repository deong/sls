#!/usr/bin/gnuplot

# generate the actual plot
set terminal postscript eps color dashed 14
set output "OUTPUTFILE"

# draw the plot
plot [0:XMAX] "DATAFILE" using 1:3 title "PLOTTITLE" with points

