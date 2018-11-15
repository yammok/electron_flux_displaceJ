#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "moparam.h"
#include "basis.h"
#include "_basis.h"

#define PI   M_PI

static void calcCom(double a, double b, double *r, double *c1, double *c2); 
static inline double grdVGTO(double ea,  double eb,   
		             double *ra, double *rb, double *x, int    ix); 
static inline double grdVGTODt(double ea,  double eb,   
		               double *ra, double *rb, double *va,double *vb,
	     		       double *x,  int    ix); 

double calcDist(double *d1, double *d2); 
double Fn(int n, double a); 

double grdScalor(int i, int j, int ix, double *x, struct nstat_t nu)
{
   int a, b, ixyz; 
   double gs=0.0; 
   double ra[NXYZ]={0.0}, rb[NXYZ]={0.0}; 

   for (ixyz=X; ixyz<=Z; ixyz++) {
      ra[ixyz]=nu.c[ixyz][i]; 
      rb[ixyz]=nu.c[ixyz][j]; 
   }

   for (a=0; a<NGAU; a++) {
      for (b=0; b<NGAU; b++) {
	 gs+=d[a]*d[b]*grdVGTO(e[a],e[b],ra,rb,x,ix); 
      }
   }

   return gs; 
}

static inline double grdVGTO(double ea,  double eb,   
		             double *ra, double *rb, double *x, int    ix)
{
   double gv, xp, xp2, ca, cb, r2; 
   double rp[NXYZ]={0.0}; 
   double ep=ea+eb; 
   double r=calcDist(ra,rb);  

   ca=2.0*ea/PI; ca=pow(ca,0.75); 
   cb=2.0*eb/PI; cb=pow(cb,0.75); 

   calcCom(ea,eb,rp,ra,rb);
   xp =calcDist(x,rp); 
   xp2=xp*xp; 
   r2 =r*r;

   gv=-4.0*PI*ca*cb*(x[ix]-rp[ix])*exp(-ea*eb/ep*r2)*Fn(2,ep*xp2); 

   return gv; 
}

double dgrdScalor(int i, int j, int ix, double *x, struct nstat_t nu)
{
   int a, b, ixyz; 
   double dsdt=0.0; 
   double ra[NXYZ]={0.0}, rb[NXYZ]={0.0}; 
   double va[NXYZ]={0.0}, vb[NXYZ]={0.0}; 

   for (ixyz=X; ixyz<=Z; ixyz++) {
      ra[ixyz]=nu.c[ixyz][i]; 
      rb[ixyz]=nu.c[ixyz][j]; 
      va[ixyz]=nu.v[ixyz][i]; 
      vb[ixyz]=nu.v[ixyz][j]; 
   }

   for (a=0; a<NGAU; a++) {
      for (b=0; b<NGAU; b++) {
         dsdt+=d[a]*d[b]*grdVGTODt(e[a],e[b],ra,rb,va,vb,x,ix); 
      }
   }

   return dsdt; 
}

static inline double grdVGTODt(double ea,  double eb,   
		               double *ra, double *rb, double *va,double *vb,
	     		       double *x,  int    ix) 
{
   int i; 
   double dvdt, xp, xp2, ca, cb, r2; 
   double rv=0.0, vxp=0.0; 
   double rp[NXYZ]={0.0}, vp[NXYZ]={0.0}; 
   double ep=ea+eb; 
   double r=calcDist(ra,rb);

   ca=2.0*ea/PI; ca=pow(ca,0.75); 
   cb=2.0*eb/PI; cb=pow(cb,0.75);

   calcCom(ea,eb,rp,ra,rb);
   calcCom(ea,eb,vp,va,vb);
   xp =calcDist(x,rp); 
   xp2=xp*xp; 
   r2 =r*r;

   for (i=X; i<=Z; i++) {
      rv +=(ra[i]-rb[i])*(va[i]-vb[i]); 
      vxp+= vp[i]*(x[i]-rp[i]); 
   }

   dvdt=-4.0*PI*ca*cb*exp(-ea*eb/ep*r2)
       *(-Fn(2,ep*xp2)*(vp[ix]+2.0*ea*eb/ep*rv*(x[ix]-rp[ix]))
         +2.0*ep*Fn(4,ep*xp2)*vxp*(x[ix]-rp[ix])); 

   return dvdt; 
}

static void calcCom(double a, double b, double *r, double *c1, double *c2)
{
   int i; 

   for (i=X; i<=Z; i++) { r[i]=(a*c1[i]+b*c2[i])/(a+b); }
}
