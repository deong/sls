/*
makeQAPuni.cc (C) Joshua Knowles 12/3/2002

Compile: gcc makeQAPuni.cc -o makeQAPuni -lm

Usage: ./makeQAPuni -h 

 ****and the command line parameters will be explained****



INTRODUCTION

This program generates symmetric multi-objective Quadratic Assignment Problems 
with one distance matrix and multiple flow matrices.

Details of the generator can be found in:

@TechReport{Knowles2002EMO,
  author =       {Joshua Knowles and David Corne},
  title =        {Instance Generators and Test Suites for the Multiobjective Quadratic Assignment Problem},
  institution =  {IRIDIA},
  number = {TR/IRIDIA/2002-25},
  year =         {2002},
  note =         {(Submitted to the {\em 2003 Evolutionary Multi-criterion Optimization Conference (EMO-2003))}
}

and we have done some simple landscape analysis in:

@InProceedings{Knowles2002HIS,
  author =       {Joshua Knowles and David Corne},
  title =        {Towards {L}andscape {A}nalyses to {I}nform the {D}esign of {H}ybrid {L}ocal
{S}earch for the {M}ultiobjective {Q}uadratic {A}ssignment {P}roblem},
  booktitle =    {Proceedings of the Second International Conference on Hybrid Intelligent Systems (HIS 2002)},
  year =         {2002},
  publisher = {IOS Press},
  note =         {(To appear, but available from http://iridia.ulb.ac.be/~jknowles/publications.html )}
}


** Please contact me - Joshua Knowles - if you have any comments, suggestions
or questions regarding this program or multicriteria QAP problems. My email
address is jknowles@ulb.ac.be

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version. 

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details. 

   The GNU General Public License is available at:
      http://www.gnu.org/copyleft/gpl.html
   or by writing to: 
        The Free Software Foundation, Inc., 
        675 Mass Ave, Cambridge, MA 02139, USA.  

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#define MAX_N 500  // max facilities - change as required
#define MAX_K 5    // max number of objectives i.e. flow matrices - change as required
#define VERBOSE 0

#define RN ran0(&seed)
  
/* Random number generator */
/* Copyright Numerical Recipes in C */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
  
double ran0(long *idum); 
/* End copyright Numerical Recipes in C */

void print_dominance();
void print_output();
double correl_val(double v, double c);

long seed;
long init_seed;
int n_fac = 30;  // number of facilities/locations in the QAP
int n_k = 2;    // number of objectives
double corr = 0;   // correlation between 1st and all other objectives
double amp; // relates to the correlation
double offset; // relates to the correlation
int max_flow = 100;
int max_dist = 100;

int d_matrix[MAX_N][MAX_N];
int f_matrix[MAX_K][MAX_N][MAX_N];

int main(int argc, char *argv[])
{
  int i,j,k;
  double r1,r2;
  int f1,fk;
  int d1;
  seed = 23453464;

  if(argc==2) 
    {
      if(strcmp("-h", argv[1])==0)
        {
          printf("You have requested help with the -h argument.\n\
\nOther optional parameters are:\n\
-n positive integer : number of facilities/locations\n\
-k positive integer : number of objectives\n\
-c real in the interval [-1,1]: correlation between objectives\n\
-f positive integer : maximum flow between facilities\n\
-d positive integer : maximum distance between locations\n\
-s positive long : random seed\n\
default values are:\n\
n_fac = 30\n\
n_k = 2\n\
corr = 0\n\
max_flow = 100\n\
max_dist = 100\n\
seed = 23453464\n");
          exit(1);
        }
    }

  if(argc>2)
    {
      // printf("argc = %d\n", argc);
      for(i=1;i<argc-1;i++)
        {
          if(strcmp("-n", argv[i])==0)
            n_fac = atoi(argv[i+1]);
          else if (strcmp("-k", argv[i])==0)
            n_k = atoi(argv[i+1]);
          else if (strcmp("-c", argv[i])==0)
            corr = atof(argv[i+1]);
          else if (strcmp("-f", argv[i])==0)
            max_flow = atoi(argv[i+1]);
          else if (strcmp("-d", argv[i])==0)
            max_dist = atoi(argv[i+1]);
          else if (strcmp("-s", argv[i])==0)
            seed = atol(argv[i+1]);
        }
    }
  if(VERBOSE)
    {
      printf("n_fac = %d\n", n_fac);
      printf("n_k = %d\n", n_k);
      printf("corr = %lf\n", corr);
      printf("max_flow = %d\n", max_flow);
      printf("max_dist = %d\n", max_dist); 
      printf("seed = %ld\n", seed);
    }

  init_seed=seed;

  if(n_fac>MAX_N)
    {
      fprintf(stderr,"Number of facilities too high. Maximum is currently set at %d\n", MAX_N);
      exit(1);
    }
  if((corr>1.0)||(corr<-1.0))
    {
      fprintf(stderr,"Correlation must be in the interval [-1,1]\n");
      exit(1);
    }
  if(n_k>MAX_K)
    {
      fprintf(stderr,"Number of objectives too high. Maximum is currently set at %d\n", MAX_K);
      exit(1);
    }

  for(i=0;i<n_fac;i++)
    for(j=0;j<n_fac;j++)
      {   
        if(i==j)
          d_matrix[i][i]=0;
        else
          {
            d1=1+(int)(max_dist*RN);
            d_matrix[i][j] = d_matrix[j][i] = d1;
          }
      }

  for(i=0;i<n_fac;i++)
    for(j=0;j<n_fac;j++)
      {   
        if(i==j)
          f_matrix[0][i][i]=0;
        else
          {
            r1 = RN;
            f1 = 1+(int)(max_flow*r1);
            f_matrix[0][i][j] = f_matrix[0][j][i] = f1;
            for(k=1;k<n_k;k++)
              {
                if(i==j)
                  f_matrix[k][i][i]=0;
                else
                  {
                    r2 = RN;
                    if(corr>=0)
                      {                
                        fk = 1+(int)(max_flow*correl_val(r1,corr));                    
                      }
                    else
                      {
                        fk = 1+(int)(max_flow*(1.0-correl_val(r1,-corr)));                     
                      }
                    f_matrix[k][i][j]=f_matrix[k][j][i]=fk;
                  }
              }
          }
      }
  print_output();
  return 0;
}

double correl_val(double v, double c)
{
  double q,w,diff,r,p=0.5;

 if(c==1)
    return(v);

  if(c==0)
    return(RN);

  do
    {
      q= RN;
      diff =q-v;
      if(diff<0)
        diff*=-1;
      w = exp(-(diff*diff)/(2.0*(1.0-pow(c,p))*(1.0-pow(c,p))))/(1.0-pow(c,p))*2.506;

      r=RN*1.0/(1.0-pow(c,p))*2.506;
    }while(w<r);
  return(q);
}

void print_output()
{
  int i,j,k;
  printf("facilities = %d objectives = %d max_distances = %d max flows = %d correlation = %4lf seed = %ld\n", n_fac, n_k, max_dist, max_flow, corr, init_seed);
  for(i=0;i<n_fac;i++)
    {
      for(j=0;j<n_fac;j++)
        {
          printf("%2d ",d_matrix[i][j]);
        }
      printf("\n");
    }
 

  for(k=0;k<n_k;k++)
    {
      printf("\n");
      for(i=0;i<n_fac;i++)
        {
          for(j=0;j<n_fac;j++)
            {
              printf("%2d ",f_matrix[k][i][j]);
            }
          printf("\n");
        }
    }
  print_dominance();
  
}

void print_dominance()
{
  int i,j,k;
  double fd;
  double dd;
  double a, b;
 
  // print flow dominances first
  for(k=0;k<n_k;k++)
    {
      b=0;
      for(i=0;i<n_fac;i++)
        {
          for(j=0;j<n_fac;j++)
            {
              b+=f_matrix[k][i][j];
            }
        }
      b /= (n_fac*n_fac);
      a=0;
      for(i=0;i<n_fac;i++)
        {
          for(j=0;j<n_fac;j++)
            {
              a+=pow((f_matrix[k][i][j]-b),2);
            }
        }
      a /=(n_fac*n_fac);
      a = pow(a,0.5);
      fd = 100 * (a/b);
      
      fprintf(stderr, "fd%d = %lf\n", k+1, fd);
    }
  b=0;
  for(i=0;i<n_fac;i++)
    {
      for(j=0;j<n_fac;j++)
        {
          b+=d_matrix[i][j];
        }
    }
  b /= (n_fac*n_fac);
  a=0;
  for(i=0;i<n_fac;i++)
    {
      for(j=0;j<n_fac;j++)
        {
          a+=pow((d_matrix[i][j]-b),2);
        }
    }
  a /=(n_fac*n_fac);
  a = pow(a,0.5);
  dd = 100 * (a/b);
      
  fprintf(stderr, "dd = %lf\n", dd);

}

/* Copyright Numerical Recipes in C */

double ran0(long *idum)
{
        long k;
        double ans;

        *idum ^= MASK;
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        ans=AM*(*idum);
        *idum ^= MASK;
        return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

/* End copyright Numerical Recipes in C */

