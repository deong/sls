#!/bin/awk -f

/fitness:/ {
    printf("%d %d ",$2,$3);
}
/^$/ {
    printf("\n");
}

