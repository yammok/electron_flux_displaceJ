#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "moparam.h"
#include "basis.h"

#define PI M_PI

static double calcElDensity(double *x, int h, int n, 
			               int na, struct stat_t st); 
static double calcAO(int iao, double *x, struct nstat_t nu); 

double calcDen1(int h, int i, int j, struct estat_t el);  
double calcDist(double *d1, double *d2);
void prtVal(double *x, double *j);

void calcDSpDist(struct input_t in, struct stat_t st)
{
   int a;
   double dm; 
   double   x[NXYZ]={0.0},  dx[NXYZ]={0.0};
   double max[NXYZ]={0.0}, min[NXYZ]={0.0};

   for (a=X; a<=Z; a++) {
      min[a]=-in.s.w[a]*0.5+in.s.o[a];
      max[a]= in.s.w[a]*0.5+in.s.o[a];
       dx[a]= in.s.d[a];
   }

   for (x[X]=min[X]; x[X]<=max[X]; x[X]+=dx[X]) {
      for (x[Y]=min[Y]; x[Y]<=max[Y]; x[Y]+=dx[Y]) {
         for (x[Z]=min[Z]; x[Z]<=max[Z]; x[Z]+=dx[Z]) {
            dm=calcElDensity(x,in.hmo,in.nmo,in.nat,st);
	    printf("%9.5lf %9.5lf %9.5lf ", x[X],x[Y],x[Z]); 
	    printf("%17.10lf\n", dm); 
         }
      }
   }
   putchar('\n');
}

static double calcElDensity(double *x, int h, int n, 
			               int na, struct stat_t st)
{
   enum atom_t {A,B}; 
   int iao, jao; 
   double dm=0.0, dmao;
   double ao[2]; 

   for (iao=0; iao<n; iao++) {
      for (jao=0; jao<n; jao++) {
	 ao[A]=calcAO(iao,x,st.nu); 
	 ao[B]=calcAO(jao,x,st.nu); 
         dmao =calcDen1(h,iao,jao,st.el);
	 dm+=dmao*ao[A]*ao[B];
      }
   }
   return dm; 
}

static double calcAO(int iao, double *x, struct nstat_t nu)
{
   int i; 
   double ao=0.0, xrn, c; 
   double rn[NXYZ]={0.0}; 

   for (i=X; i<=Z; i++) { rn[i]=nu.c[i][iao]; }
   xrn=calcDist(x,rn); 
   for (i=0; i<NGAU; i++) { 
      c=2.0*e[i]/PI; c=pow(c,0.75); 
      ao+=d[i]*c*exp(-e[i]*xrn*xrn); 
   }

   return ao; 
}
