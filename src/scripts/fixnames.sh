#!/bin/bash

for h in {30,}
do
    for i in {uni,rl}
    do
        for j in {1,2,3}
        do
            cd Gar${h}-2fl-${j}${i}
            mv emoeats/emoeats.cfg emoeats/emoeatpts.cfg
            mv emoeats/emoeats.out emoeats/emoeatpts.out
            mv emoeatpts/emoeatpts.cfg emoeatpts/emoeats.cfg
            mv emoeatpts/emoeatpts.out emoeatpts/emoeats.out
            mv emoeats tmp
            mv emoeatpts emoeats
            mv tmp emoeatpts
            cd ..
        done
    done
done
