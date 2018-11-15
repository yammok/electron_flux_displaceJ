#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "moparam.h"
#include "err.h"

#define PI M_PI

static void calcElCurrent(double *j, double *x, 
			int h, int n, int na, struct stat_t st); 
static double calcDotDen1(int h, int n, int iao, int jao, int na, 
			struct estat_t el, struct nstat_t nu); 
static double grdDen(struct estat_t el, int h, int n, 
				 int x, int a, int iao, int jao); 

double calcDen1(int h, int iao, int jao, struct estat_t el); 
double  grdScalor(int i, int j, int ix, double *x, struct nstat_t nu); 
double dgrdScalor(int i, int j, int ix, double *x, struct nstat_t nu); 
void prtVal(double *x, double *j); 

void calcJSpDist(struct input_t in, struct stat_t st)
{
   int a; 
   double   x[NXYZ]={0.0},   d[NXYZ]={0.0}; 
   double max[NXYZ]={0.0}, min[NXYZ]={0.0}; 

   for (a=X; a<=Z; a++) {
      min[a]=-in.s.w[a]/2+in.s.o[a]; 
      max[a]= in.s.w[a]/2+in.s.o[a]; 
	d[a]= in.s.d[a]; 
   }

   for (x[X]=min[X]; x[X]<=max[X]; x[X]+=d[X]) {
      for (x[Y]=min[Y]; x[Y]<=max[Y]; x[Y]+=d[Y]) {
	 for (x[Z]=min[Z]; x[Z]<=max[Z]; x[Z]+=d[Z]) {
	    double j[NXYZ]={0.0}; 

	    calcElCurrent(j,x,in.hmo,in.nmo,in.nat,st); 
	    prtVal(x,j); 
	 }
      }
   }
   putchar('\n'); 
}

static void calcElCurrent(double *j, double *x, 
			int h, int n, int na, struct stat_t st)
{
   int iao, jao, ix; 
   double dmdt, dm; 

   for (iao=0; iao<n; iao++) {
      for (jao=0; jao<n; jao++) {
	 dmdt=calcDotDen1(h,n,iao,jao,na,st.el,st.nu); 
	 for (ix=X; ix<=Z; ix++) {
	    j[ix]+=dmdt*grdScalor(iao,jao,ix,x,st.nu); 
	 }

	 dm=calcDen1(h,iao,jao,st.el); 
	 for (ix=X; ix<=Z; ix++) {
	    j[ix]+=dm*dgrdScalor(iao,jao,ix,x,st.nu); 
	 }
      }
   }
   for (ix=X; ix<=Z; ix++) { j[ix]/=4.0*PI; }
}

double calcDen1(int h, int iao, int jao, struct estat_t el)
{
   int imo;
   double dm=0.0; 

   for (imo=0; imo<h; imo++) { dm+=el.cao[iao][imo]*el.cao[jao][imo]; }
   return 2.0*dm; 
}

static double calcDotDen1(int h, int n, int iao, int jao, int na, 
			struct estat_t el, struct nstat_t nu)
{
   int ia, ix; 
   double dden=0.0; 

   for (ia=0; ia<na; ia++) {
      for (ix=X; ix<=Z; ix++) {
	 dden+=nu.v[ix][ia]*grdDen(el,h,n,ix,ia,iao,jao); 
      }
   }
   return dden; 
}

static double grdDen(struct estat_t el, int h, int n, 
				 int x, int a, int iao, int jao)
{
/* H2 sto-ng only */
   return -2.0*el.dsmo[x][a][0][0]*el.cao[iao][0]*el.cao[jao][0]; 
#if 0
   double prt_s=0.0, prt_u=0.0; 
   int imo, mmo; 

   for (imo=0; imo<h; imo++) {
      for (mmo=0; mmo<n-h; mmo++) {
	 prt_u+=(el.cao[i][mmo]*el.cao[j][imo] 
	       + el.cao[i][imo]*el.cao[j][mmo])*el.u[x][a][mmo][imo]; 
      }
   }

   for (imo=0; imo<h; imo++) {
      for (mmo=0; mmo<h; mmo++) {
	 prt_s-=(el.cao[i][mmo]*el.cao[j][imo] 
	       + el.cao[i][imo]*el.cao[j][mmo])*el.dsmo[x][a][mmo][imo];  
      }
   }
   prt_s*=0.5;

   return 2.0*(prt_u+prt_s);
#endif 
}
