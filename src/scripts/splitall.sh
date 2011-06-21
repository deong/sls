#!/bin/bash

for h in {30,100}
do
    for i in {uni,rl}
    do
        for j in {1,2,3}
        do
            for k in {spea2ts,emoeats,emoeatpts}
            do
                cd Gar${h}-2fl-${j}${i}/$k
                split.pl ${k}.out
                cd ../..
            done
        done
    done
done
