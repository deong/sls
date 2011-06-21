#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv)
{
    FILE* in;
    char  line[256];
    int   i=0;
    
    if(argc != 2)
    {
        fprintf(stderr,"usage: splitlndscp.exe datafile\n");
        exit(1);
    }
    
    if((in=fopen(argv[1],"r"))==NULL)
    {
        fprintf(stderr,"could not open datafile: %s\n", argv[1]);
        exit(1);
    }
    
    while(!feof(in))
    {
        fgets(line,256,in);
        if(i++==50)
        {
            printf("\n\n");
            i=0;
        }
        else
        {
            fputs(line,stdout);
        }
    }
    fclose(in);
    return 0;
}

    
