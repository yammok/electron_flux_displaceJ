#include <math.h>
#include "moparam.h"

#define PI M_PI
#define TLIM 0.01
#define DT   0.01

double calcDist(double *d1, double *d2)
{
   int i;
   double d=0.0; 

   for (i=X; i<=Z; i++) { d+=(d1[i]-d2[i])*(d1[i]-d2[i]); }

   return sqrt(d); 
}

#if 0
double Fn(int n, double x) 
{
   if (x<=TLIM) { 
      double nd=(double)n; 

      return 1.0/(nd+1.0)-1.0/(nd+3.0)*x; 
   }

   if      (n==0) { return 0.5*sqrt(PI/x)*erf(sqrt(x)); } 
   else if (n==1) { return 0.5*(1.0-exp(-x))/x; }   

   return 0.5*(n*Fn(n-2,x)-exp(-x))/x; 
}
#endif 

double Fn(int n, double x) 
{
   double fn=0.0;
   double t;

   for (t=0.0; t<=1.0; t+=DT) {
      double tn=pow(t,n);
      double t2=pow(t,2);

      fn+=DT*tn*exp(-x*t2);
   }
   return fn;
}
