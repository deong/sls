#!/bin/bash

for h in {30,}
do
    for i in {uni,rl}
    do
        for j in {1,2,3}
        do
            cd Gar${h}-2fl-${j}${i}
            echo "plothv" | matlab -nojvm -nosplash
            cd ..
        done
    done
done
