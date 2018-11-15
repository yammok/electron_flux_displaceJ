#include <stdio.h>
#include <stdlib.h>
#include "moparam.h"
#include "err.h"

void allocMem(double ***c, int m, int n); 
void freeMem(double ***a, int n); 

static void allocCPHFArray(double ****u, int h, int n, int na); 
static void allocNuMem(struct nstat_t *n, int nu); 

void allocMOArray(struct stat_t *s, struct input_t in)
{
   int i; 

   allocMem(&s->el.cao,in.nmo,in.nmo); 
   for (i=X; i<=Z; i++) {
      allocCPHFArray(&s->el.u[i]   ,in.hmo,in.nmo,in.nat); 
      allocCPHFArray(&s->el.dsmo[i],in.nmo,in.nmo,in.nat); 
   }
   allocNuMem(&s->nu,in.nat); 
}

static void allocNuMem(struct nstat_t *n, int na)
{
   int i; 

   n->elem=malloc(sizeof(double)*na); 
   if (n->elem==NULL) {
      ERRMSG("Failure allocating memory");
      exit(1); 
   }

   for (i=X; i<=Z; i++) { 
      n->c[i]=malloc(sizeof(double) * na); 
      if (n->c[i]==NULL) {
         ERRMSG("Failure allocating memory");
	 goto fail; 
      }
   }

   for (i=X; i<=Z; i++) { 
      n->v[i]=malloc(sizeof(double)*na); 
      if (n->v[i]==NULL) {
         ERRMSG("Failure allocating memory");
	 goto fail; 
      }
   }

   return ; 
fail: 
   for (i=X; i<=Z; i++) {
      free(n->c[i]); 
      free(n->v[i]); 
   }
   free(n->elem); 

   exit(1); 
} 

static void allocCPHFArray(double ****u, int h, int n, int na)
{
   int i; 

   *u=malloc(sizeof(double **)*na); 
   if (*u==NULL) {
      ERRMSG("Failure allocating memory");
      exit(1); 
   }

   for (i=0; i<na; i++) {
      (*u)[i]=malloc(sizeof(double *)*h); 
      if ((*u)[i]==NULL) {
         ERRMSG("Failure allocating memory");
	 goto fail;
      }
   }

   for (i=0; i<na; i++) {
      int j; 

      for (j=0; j<h; j++) {
 	 (*u)[i][j]=malloc(sizeof(double)*(n-h)); 
	 if ((*u)[i][j]==NULL) {
	    ERRMSG("Failure allocating memory");
 	    goto fail;
	 }
      }
   }

   return ; 
fail: 
   for (i=0; i<na; i++) {
      int j; 

      for (j=0; j<h; j++) { free((*u)[i][j]); }	
   }
   for (i=0; i<na; i++) { free((*u)[i]); }
   free(*u); 

   exit(1); 
}

void allocMem(double ***c, int m, int n)
{
   int i; 

   *c=malloc(sizeof(double *) * m); 
   if ((*c)==NULL) {
      ERRMSG("Failure allocating memory");
      goto fail;
   }
   for (i=0; i<m; i++) {
      (*c)[i]=malloc(sizeof(double) * n); 
      if ((*c)[i]==NULL) {
       	 ERRMSG("Failure allocating memory");
         goto fail;
      }
   }

   return ; 
fail:
	freeMem(c,n); 
	exit(1); 
}

void freeMem(double ***a, int n)
{
   int i;

   for (i=0; i<n; i++) { free((*a)[i]); }
   free((*a));
}
