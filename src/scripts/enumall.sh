#!/bin/bash

#for r in {-0.4,0.4}
for r in {0.4,}
do
    for num in {1..30}
    do
        for type in {C,D,E}
        do
            sed "s/1C/$num$type/g" ../cfg/enumerate_gap.cfg | sed "s/r-0.4/r$r/g" > current.cfg
            ./enumerate gap current.cfg > ../prob/mgap/toy/5x12/r$r/toy_$num$type.front
        done
    done
done
