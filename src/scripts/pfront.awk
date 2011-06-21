#!/bin/awk -f

/fitness:/ {
    print $2,$3;
}
/^$/ {
    print $0;
}

