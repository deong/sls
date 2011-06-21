/*
makeQAPrl.cc (C) Joshua Knowles 27/3/2002

Compile: gcc makeQAPrl.cc -o makeQAPrl -lm

Usage: ./makeQAPrl -h 

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


STRUCTURED INSTANCES

The distance and flow matrices generated have structured, "real-world-like" entries 
and correlations can be set between corresponding entries in the two
or more flow matrices. The amount of "overlap" between the flow matrices can
also be controlled. The generator follows procedures for making non-uniform
QAP problems outlined in:

E. D. Taillard. Comparison of iterative searches for the quadratic assignment 
problem. Location Science, 3:87--105, 1995.

DISTANCE MATRIX

The distance matrix entries are the Euclidean distances between points 
in the plane. The points are randomly distributed in small circular regions,
with these regions distributed in a larger circle. The size and number of
the small and larger circles can be controlled by the command line parameters
R_max, r_max, and N_max.

The points are generated using the procedure:

Choose Theta randomly, uniformly between 0 and 2 * PI
Choose R randomly, uniformly between 0 and R_max
Choose N randomly, uniformly between 1 and N_max
Repeat N times:
   Choose theta randomly, uniformly between 0 and 2 * PI
   Choose r randomly, uniformly between 0 and r_max.
   The Euclidean coordinates of the next generated point are:
   (R cos Theta + r cos theta, R sin Theta + r sin theta)


FIRST FLOW MATRIX

The flow entries are non-uniform random values, controlled by two parameters,
A and B, with A < B, and B > 0. Let X be a random variable uniformly distributed
between 0 and 1. Then a flow entry is given by (int)(10^((B-A)*X + A)) -- Equation 1. 
With negative vaues of A the flow matrix is sparse, ie it contains a number of
zero entries. The non-zero entries are non-uniformly distributed. 

OTHER FLOW MATRICES

The number of flow matrices is controlled by the command line parameter, k.
With k >= 2, a multiobjective QAP problem is generated with k flow matrices.
The entries in the nth matrix (k>=n>=2) are generated using the same equation as
above but the random variable X is correlated with the value of X
used in the corresponding entry in the first flow matrix. Correlations between
-1 and 1 can be set
between the first and each of the additional flow matrices using the command line 
parameters.

A degree of overlap between the matrices can also be specified. The overlap
parameter is set between 0 and 1. It controls the fraction of entries in
the nth flow matrix that are correlated with the corresponding entries in the
1st flow matrix. With the overlap parameter set to zero, a random un-correlated
value, calculated using Equation 1, will be placed in each entry of the 
nth flow matrix that corresponds to a zero entry in the first flow matrix.
Similarly, a zero will be placed in each entry of the nth flow matrix that
corresponds to a non-zero value in the first flow matrix. Thus there is no
overlap between the flows of the first and nth matrix. With the overlap
parameter set to 1 all the flows are correlated. With the overlap set to
intermediate values some of the flows will overlap and others will not.

Here is a typical command line for the generator:

./makeQAPrl -n 30 -k 3 -M 1000 -m 10 -K 4 -A -3 -B 4 -c1 0.5 -c2 -0.8 -ov 0.7  

--

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
#define VERBOSE 1
#define PI 3.141592654

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


void print_output();
void point_in_plane(const int M, const int K, const int m, const int max_points);
void print_dominance();
double correl_val(double v, double c);


long seed;
long init_seed;
int n_fac = 30;  // number of facilities/locations in the QAP
int n_k = 2;    // number of objectives
double corr[MAX_K];   // correlations between 1st and each other objective
int max_flow;
int max_dist;

int M=0;
int K=1;
int m=100;
int A=-5;
int B=5;
double overlap=1.0;

int d_matrix[MAX_N][MAX_N];
int f_matrix[MAX_K][MAX_N][MAX_N];
int XYpos[MAX_N][2];

int main(int argc, char *argv[])
{
  int i,j,k;
  double r1,r2,val;
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
-c1 real in the interval [-1,1]: correlation between objectives 1 and 2\n\
-c2 real in the interval [-1,1]: correlation between objectives 1 and 3\n\
-c3 real in the interval [-1,1]: correlation between objectives 1 and 4\n\
-c4 real in the interval [-1,1]: correlation between objectives 1 and 5\n\
-ov real in the interval [0:1]: sets the fraction of flows that are correlated\n\
-A integer : flow control parameter (must be set less than B; negative values cause a sparse matrix)\n\
-B positive integer : flow control parameter\n\
-K positive integer : maximum number of points in a small cluster \n\
-M non-negative integer : radius of large clusters\n\
-m non-negative integer : radius of small clusters\n\
-s positive long : random seed\n\
\n\
Default values are:\n\
n = 30, k = 2, c1 = 0, c2 = 0, c3 = 0, c4 = 0,\n\
A = -5, B = 5, m = 100, M = 0, K = 1, ov = 0.7,\n\
s = 23453464\n");
          exit(1);
        }
    }

  if((argc>2)&&(argc%2==1))
    {
      // printf("argc = %d\n", argc);
      for(i=1;i<argc-1;i+=2)
        {
          if(strcmp("-n", argv[i])==0)
            n_fac = atoi(argv[i+1]);
          else if(strcmp("-A", argv[i])==0)
            A = atoi(argv[i+1]);
          else if(strcmp("-B", argv[i])==0)
            B = atoi(argv[i+1]);
          else if(strcmp("-M", argv[i])==0)
            M = atoi(argv[i+1]);
          else if(strcmp("-K", argv[i])==0)
            K = atoi(argv[i+1]);
          else if(strcmp("-m", argv[i])==0)
            m = atoi(argv[i+1]);
          else if(strcmp("-ov", argv[i])==0)
            overlap = atof(argv[i+1]);
          else if (strcmp("-k", argv[i])==0)
            n_k = atoi(argv[i+1]);
          else if (strcmp("-c1", argv[i])==0)
            corr[1] = atof(argv[i+1]);
          else if (strcmp("-c2", argv[i])==0)
            corr[2] = atof(argv[i+1]);
          else if (strcmp("-c3", argv[i])==0)
            corr[3] = atof(argv[i+1]);
          else if (strcmp("-c4", argv[i])==0)
            corr[4] = atof(argv[i+1]);
          else if (strcmp("-s", argv[i])==0)
            seed = atol(argv[i+1]);
          else
            {
              fprintf(stderr, "Undefined command line parameter entered. Do \"./makeQAPrl -h\" for help with parameters.\n");
              exit(1);
            }
        }
    }
  else
    {
      fprintf(stderr, "No params set / odd number of params set. Using defaults. OK? Press a key to continue\n");
      getchar();
    }

  if(VERBOSE)
    {
      fprintf(stderr,"n_fac = %d\n", n_fac);
      fprintf(stderr,"n_k = %d\n", n_k);
      for(k=1;k<n_k;k++)
        fprintf(stderr,"corr[%d] = %lf\n", k, corr[k]);
      fprintf(stderr,"overlap = %lf\n", overlap);
      fprintf(stderr,"A = %d\n", A);
      fprintf(stderr,"B = %d\n", B);
      fprintf(stderr,"M = %d\n", M);
      fprintf(stderr,"m = %d\n", m);
      fprintf(stderr,"K = %d\n", K);
      fprintf(stderr,"seed = %ld\n", seed);
    }

  init_seed=seed;

  if(n_fac>MAX_N)
    {
      printf("Number of facilities too high. Maximum is currently set at %d\n", MAX_N);
      exit(1);
    }
  for(k=1;k<n_k;k++)
    {
      if((corr[k]>1.0)||(corr[k]<-1.0))
        {
          fprintf(stderr, "Correlations must be in the interval [-1,1]\n");
          exit(1);
        }
    }
  if(n_k>MAX_K)
    {
      fprintf(stderr,"number of objectives too high. Maximum is currently set at %d\n", MAX_K);
      exit(1);
    }

  if(B<0)
    {
      fprintf(stderr,"B must be a positive integer.\n");
      exit(1);
    }

  if(A>=B)
    {
      fprintf(stderr,"A must be less than B.\n");
      exit(1);
    }

  if(K<0)
    {
      fprintf(stderr,"K must be a positive integer.\n");
      exit(1);
    }

  if(M<0)
    {
      fprintf(stderr,"M must be a non-negative integer\n");
      exit(1);
    }

  if(m<0)
    {
      fprintf(stderr,"m must be a non-negative integer\n");
      exit(1);
    }

  point_in_plane(M, K, m, n_fac);

  //  for(i=0;i<n_fac;i++)
  //  printf("%d %d\n",XYpos[i][0], XYpos[i][1]);
  
  max_dist=0;
  for(i=0;i<n_fac;i++)
    for(j=0;j<n_fac;j++)
      {   
        d1=(int)(sqrt((XYpos[i][0]-XYpos[j][0])*(XYpos[i][0]-XYpos[j][0]) + ((XYpos[i][1]-XYpos[j][1]))*((XYpos[i][1]-XYpos[j][1]))));
        d_matrix[i][j] = d_matrix[j][i] = d1;      
        if(d1>max_dist)
          max_dist=d1;
      }

  max_flow=0;
  for(i=0;i<n_fac;i++)
    for(j=0;j<n_fac;j++)
      {   
        if(i==j)
          f_matrix[0][i][i]=0;
        else
          {
            r1 = RN;
            f1 = (int)(pow(10,((B-A)*r1+A)));
            f_matrix[0][i][j] = f_matrix[0][j][i] = f1;
            if(f1>max_flow)
              max_flow=f1;
            for(k=1;k<n_k;k++)
              {
                if(i==j)
                  f_matrix[k][i][i]=0;
                else
                  {
                    r2=RN;
                    if(r2>overlap) // if we want an un-correlated value
                      {
                        if(f_matrix[0][i][j]==0)                                                
                          {
                            r2=RN*B;
                            fk = (int)(pow(10,r2));
                            f_matrix[k][i][j] = f_matrix[k][j][i] = fk;
                          }
                        else
                          {
                            fk=0;
                            f_matrix[k][i][j] = f_matrix[k][j][i] = fk;
                          }
                      }
                    else
                      {
                        r2 = RN;
                        if(corr[k]>=0)
                          {
                            //val = corr[k]*r1+(1.0-corr[k])*r2;
                            val = correl_val(r1,corr[k]);
                            //fprintf(stderr, "val=%lf\n", val);
                          }
                        else
                          {
                            val = 1.0-correl_val(r1,-corr[k]);
                            //fprintf(stderr, "val=%lf\n", val);
                          }
                        fk = (int)(pow(10,((B-A)*val+A)));
                        f_matrix[k][i][j]=f_matrix[k][j][i]=fk;
                      }
                  }
                if(fk>max_flow)
                  max_flow=fk;
                // fprintf(stderr,"f1 = %d, fk=%d, max_flow=%d\n",f1, fk, max_flow);
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
  printf("facilities = %d objectives = %d max_distances = %d max flows = %d overlap =  %.2f seed = %ld ", n_fac, n_k, max_dist, max_flow, overlap, init_seed);
  for(k=1;k<n_k;k++)
    {
      printf("correlation%d = %.3lf ", k, corr[k]);
    }
  printf("A = %d B = %d M = %d K = %d m = %d\n", A, B, M, K, m);

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


void point_in_plane(const int M, const int K, const int m, const int max_points)
{
  // Makes points in the plane for constructing the distance matrix for real-world-like QAP problems. 
  // The function follows the procedure outlined in Taillard'95 - BibTeX reference at top of program.

  double Theta;
  double theta;
  int R, r, N, i;
  int num_points=0;

  while(num_points<max_points)    
    {      
      Theta = RN*2*PI;
      R = (int)(RN * M); 
      N = (int)(RN * K)+1;
      for(i=0;i<N;i++)
        {
          if(num_points<max_points)
            {
              theta = RN*2*PI;
              r = (int)(RN *m);
              XYpos[num_points][0]=(int)(R * cos(Theta) + r * cos(theta));
              XYpos[num_points][1]=(int)(R * sin(Theta) + r * sin(theta));
              num_points++;
              // printf("%d\n", num_points);
            }
          else
            return;
        }      
    }
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

