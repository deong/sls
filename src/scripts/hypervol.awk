#!/bin/awk -f

/evaluations:/ {
    printf "%d ", $2;
}

/generations:/ {
    printf "%d ", $2;
}

/hypervolume:/ {
    print $2;
}

