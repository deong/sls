#!/bin/bash

ROOT=/home/deong/research/sls

for size in {3x17,5x12}
do
    for r in {-0.4,0.0,0.4}
    do
        for num in {1..30}
        do
            for type in {C,D,E}
            do
                $ROOT/src/pearson $ROOT/prob/mgap/toy/$size/r$r/toy_$num$type.ifdc > $ROOT/prob/mgap/toy/$size/r$r/toy_$num$type.r
            done
        done
    done
done
