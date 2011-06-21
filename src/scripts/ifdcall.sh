#!/bin/bash

for size in {3x17,5x12}
do
    for r in {-0.4,0.4}
    do
        for num in {1..30}
        do
            for type in {C,D,E}
            do
                sed "s/1C/$num$type/g" ../cfg/lndscp.cfg | sed "s/3x17/$size/g" | sed "s/r-0.4/r$r/g" > current.cfg
                ./sls current.cfg > current.out
                grep -v 'trial' current.out | sort -n | uniq | grep -v '0 0 0' > ../prob/mgap/toy/$size/r$r/toy_$num$type.ifdc
                rm current.out
            done
        done
    done
done
