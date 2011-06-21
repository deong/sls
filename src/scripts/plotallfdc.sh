#!/bin/bash

# move to the appropriate root directory
cd /home/deong/research/sls/prob/mgap/toy/

for size in {3x17,5x12}
do
    # compute the maximum for the x-axis of the plot
    MAXX=${size##*x}

    for r in {-0.4,0.0,0.4}
    do
        # change to the right subdirectory
        cd $size/r$r
        
        for num in {1..30}
        do
            for type in {C,D,E}
            do
                # create the gnuplot commands to build the plot
                sed "s/DATAFILE/toy_$num$type.ifdc/g" /home/deong/research/sls/src/scripts/plotallfdc.gnuplot | \
                    sed "s/OUTPUTFILE/toy_$num$type.eps/g" | \
                    sed "s/XMAX/$MAXX/g" | \
                    sed "s/PLOTTITLE/FDC ${size}x2-$num$type (r=$r)/g" | gnuplot
            done
        done
        
        # go back up to the root
        cd ../..
    done
done