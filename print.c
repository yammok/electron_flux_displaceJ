#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include "moparam.h"
#include "err.h"

void prtMat(double **mat, int m, int n)
{
   int i, j; 

   for (i=0; i<m; i++) {
      for (j=0; j<n; j++) { fprintf(stderr, "%17.10lf", mat[i][j]); }
      fprintf(stderr, "\n");
   }
}

void prtVal(double *x, double *j)
{
   int i; 

   for (i=X; i<=Z; i++) { printf("%9.5lf ", x[i]); }
   for (i=X; i<=Z; i++) { printf("%17.10lf ", j[i]); }
   putchar('\n'); 
}

void printPoints(struct space_t s, char *iname)
{
   FILE *fp; 
   char pfile[NBUF]={'\0'}; 
   int a, n[NXYZ]={0}, ntot=0; 
   double   x[NXYZ]={0.0},   d[NXYZ]={0.0};
   double max[NXYZ]={0.0}, min[NXYZ]={0.0};

   for (a=X; a<=Z; a++) {
      min[a]=-s.w[a]/2+s.o[a];
      max[a]= s.w[a]/2+s.o[a];
        d[a]= s.d[a];
   }

   for (x[X]=min[X]; x[X]<=max[X]; x[X]+=d[X]) { n[X]++; }
   for (x[Y]=min[Y]; x[Y]<=max[Y]; x[Y]+=d[Y]) { n[Y]++; }
   for (x[Z]=min[Z]; x[Z]<=max[Z]; x[Z]+=d[Z]) { n[Z]++; }

   for (x[Z]=min[Z]; x[Z]<=max[Z]; x[Z]+=d[Z]) {
      for (x[Y]=min[Y]; x[Y]<=max[Y]; x[Y]+=d[Y]) {
         for (x[X]=min[X]; x[X]<=max[X]; x[X]+=d[X]) {
	    ntot++; 
         }
      }
   }

   strncpy(pfile,iname,NBUF-1); 
   strncat(pfile,".point",strlen(".point"));
   fp=fopen(pfile,"w"); 
   if (fp==NULL) {
      ERRMSG("Failure oppening the point file"); 
      exit(0); 
   }

   fprintf(fp,"%d %d %d\n",n[X],n[Y],n[Z]); 
   fprintf(fp,"%d\n",ntot); 

   fclose(fp); 
}

void printInpData(struct input_t in)
{
   int i; 

   fprintf(stderr, "Type: %s\n", in.type); 
   fprintf(stderr, "nat: %d\n", in.nat);
   fprintf(stderr, "nmo: %d\n", in.nmo); 
   fprintf(stderr, "hmo: %d\n", in.hmo); 

   fprintf(stderr, "%-10s: ", "origin"); 
   for (i=X; i<=Z; i++) { fprintf(stderr, "%14.7lf ", in.s.o[i]); }
   fprintf(stderr, "\n"); 

   fprintf(stderr, "%-10s: ", "mesh"); 
   for (i=X; i<=Z; i++) { fprintf(stderr, "%14.7lf ", in.s.d[i]); }
   fprintf(stderr, "\n"); 

   fprintf(stderr, "%-10s: ", "box size"); 
   for (i=X; i<=Z; i++) { fprintf(stderr, "%14.7lf ", in.s.w[i]); }
   fprintf(stderr, "\n"); 
}

void printGrdMat(double ***mx[NXYZ], int na, int m, int n)
{
   int ia, ix; 

   for (ia=0; ia<na; ia++) {
      fprintf(stderr, "atom %d\n", ia); 
      for (ix=X; ix<=Z; ix++) {
	 fprintf(stderr, "XYZ: %d\n", ix); 
	 prtMat(mx[ix][ia], m, n); 
      }
   }
}

void printStat(struct stat_t s, int na, int nao)
{
   int i, j, k, l;

   fprintf(stderr, "=== Nuclear    status ===\n"); 
   fprintf(stderr, "Position\n"); 
   for (i=0; i<na; i++) { 
      for (j=X; j<=Z; j++) { fprintf(stderr,"%9.5lf ", s.nu.c[j][i]); }
      fprintf(stderr, "\n"); 
   }

   fprintf(stderr, "Velocity\n"); 
   for (i=0; i<na; i++) { 
      for (j=X; j<=Z; j++) { fprintf(stderr, "%9.5lf ", s.nu.v[j][i]); }
      fprintf(stderr, "\n"); 
   }

   fprintf(stderr, "=== electronic status ===\n"); 
   fprintf(stderr, "MO coefficients\n"); 
   for (i=0; i<nao; i++) {
      for (j=0; j<nao; j++) { fprintf(stderr, "%9.5lf ", s.el.cao[j][i]); }
      fprintf(stderr, "\n"); 
   }

   fprintf(stderr, "Overlap derivative matrix (AO)\n"); 
   for (l=0; l<na; l++) {
      for (k=X; k<=Z; k++) {
         fprintf(stderr, "%-3d %3d\n", l,k);
         for (i=0; i<nao; i++) {
            for (j=0; j<nao; j++) { fprintf(stderr, "%9.5lf ", 
				            s.el.dsmo[k][l][i][j]); }
            fprintf(stderr, "\n"); 
         }
      }
   }
}
