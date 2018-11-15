#include <stdio.h>
#include "moparam.h"

void initMat(double **cao, int m, int n); 

static void initCPHFMat(double ***u, int h, int n, int na); 
static void initGeomMat(struct nstat_t *nu, int na); 

void prtMat(double **mat, int m, int n);

void initMOArray(struct stat_t *s, struct input_t in)
{
   int i; 

   initMat(s->el.cao,in.nmo,in.nmo); 
   for (i=X; i<=Z; i++) { 
      initCPHFMat(s->el.u[i]   ,in.hmo,in.nmo,in.nat); 
      initCPHFMat(s->el.dsmo[i],in.nmo,in.nmo,in.nat); 
   }
   initGeomMat(&s->nu,in.nat); 
} 

static void initGeomMat(struct nstat_t *nu, int na)
{
   int i,j; 

   for (i=0; i<na; i++) { nu->elem[i]=0; }
   for (i=X; i<=Z; i++) {
      for (j=0; j<na; j++) { nu->c[i][j]=0.0; }
      for (j=0; j<na; j++) { nu->v[i][j]=0.0; }
   }
} 

void initMat(double **cao, int m, int n)
{
   int i,j; 

   for (i=0; i<m; i++) {
      for (j=0; j<n; j++) { cao[i][j]=0.0; }
   }
}

static void initCPHFMat(double ***u, int h, int n, int na)
{
   int i,j,k; 

   for (i=0; i<na; i++) {
      for (j=0; j<h; j++) {
	 for (k=0; k<n-h; k++) { u[i][j][k]=0.0; }
      }
   }
}
