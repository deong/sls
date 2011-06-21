#!/bin/bash

for h in {30,}
do
    for i in {uni,rl}
    do
        for j in {1,2,3}
        do
            for k in {spea2ts,emoeats,emoeatpts}
            do
                cd Gar${h}-2fl-${j}${i}/$k
                for t in trial*.dat
                do
                    cat $t | hypervol.awk > ${t%%.dat}.hv
                done
                cd ../..
            done
        done
    done
done
