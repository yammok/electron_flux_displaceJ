#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "moparam.h"
#include "err.h"

void allocMem(double ***c, int m, int n);
void freeMem(double ***a, int n);
void initMat(double **cao, int m, int n);
void transToAu(double **v, int n); 

static void loadData(char *s, double *d); 
static void loadMatData(FILE *fp, double **c, int n); 
static void loadDSMO(FILE *fp, double ****dsmo, double **c, int n, int na); 
static void loadNuGeom(int na, FILE *fp, double **c); 
static void loadNuVel(int na, FILE *fp, double **v); 
static void genDSMO(double **dsmo, double **dsao, double **c, int n); 

static double transMat(double **mat, double **c, int n, int lmo, int mmo); 

void loadInpData(struct input_t *in, char *iname, char *type)
{
   FILE *fp; 
   char fname[NBUF]={'\0'}; 
   int i;

   strncpy(fname,iname,NBUF-1); 
   strncat(fname,".cur",strlen(".cur")); 

   fp=fopen(fname, "r");  
   if (fp==NULL) {
      ERRMSG("Failure to open 'cur' file");
      exit(1);
   }

   if      (strcmp(type,"-d")==0) { strncpy(in->type,"density",NBUF-1); }
   else if (strcmp(type,"-f")==0) { strncpy(in->type,"flux",NBUF-1);    }

   for (i=0; ; i++) {
      char s[NBUF]={'\0'}; 

      fgets(s,NBUF-1,fp);
      if (feof(fp)) { break; }
      switch (i) {
      case 0:
         loadData(s,in->s.o); 
         break;
      case 1:
         loadData(s,in->s.d); 
         break;
      case 2:
         loadData(s,in->s.w); 
         break;
      case 3:
         in->hmo=atoi(s);
         break;
      case 4:
         in->nmo=atoi(s);
         break;
      case 5:
         in->nat=atoi(s);
         break;
      default:
         break;
      }
   }
   fclose(fp); 
}

static void loadData(char *s, double *d)
{
   char *sp; 
   int i=0; 

   for (sp=s; ; sp=NULL) {
      char *tp=strtok(sp," "); 

      if (tp==NULL) { break; }
      d[i++]=atof(tp); 
   }
}

void loadElData(struct input_t in, FILE **fp, struct stat_t *stat)
{
   int i, j; 

   loadMatData(fp[MO], stat->el.cao, in.nmo); 
   for (i=X; i<=Z; i++) {
      for (j=0; j<in.nat; j++) {
	 loadMatData(fp[CPHF],stat->el.u[i][j],in.hmo); 
      }
   }
   loadDSMO(fp[DSMO],stat->el.dsmo,stat->el.cao,in.nmo,in.nat); 
} 

static void loadDSMO(FILE *fp, double ****dsmo, double **c, int n, int na)
{
   /* Note: c[iao][imo] */
   double **dsao; 

   allocMem(&dsao,n,n); 

   initMat(dsao,n,n); 
   loadMatData(fp,dsao,n); 
   dsmo[X][0][0][0]=2.0*c[1][0]*c[0][0]*dsao[1][0]; 
   initMat(dsao,n,n); 
   loadMatData(fp,dsao,n); 
   dsmo[Y][0][0][0]=2.0*c[1][0]*c[0][0]*dsao[1][0]; 
   initMat(dsao,n,n); 
   loadMatData(fp,dsao,n); 
   dsmo[Z][0][0][0]=2.0*c[1][0]*c[0][0]*dsao[1][0]; 

   initMat(dsao,n,n); 
   loadMatData(fp,dsao,n); 
   dsmo[X][1][0][0]=2.0*c[1][0]*c[0][0]*dsao[0][1]; 
   initMat(dsao,n,n); 
   loadMatData(fp,dsao,n); 
   dsmo[Y][1][0][0]=2.0*c[1][0]*c[0][0]*dsao[0][1]; 
   initMat(dsao,n,n); 
   loadMatData(fp,dsao,n); 
   dsmo[Z][1][0][0]=2.0*c[1][0]*c[0][0]*dsao[0][1]; 

   freeMem(&dsao,n); 
}

#if 0
static void loadDSMO(FILE *fp, double ****dsmo, double **c, int n, int na)
{
   int i, ix, ia; 
   double **dsao; 

   allocMem(&dsao,n,n); 
   for (ia=0; ia<na; ia++) {
      for (ix=X; ix<=Z; ix++) {
	 initMat(dsao,n,n); 
	 loadMatData(fp,dsao,n); 
	 genDSMO(dsmo[ix][ia],dsao,c,n); 
      }
   }
   freeMem(&dsao,n); 
}
#endif 

static void genDSMO(double **dsmo, double **dsao, double **c, int n)
{
   int imo, jmo; 

   for (imo=0; imo<n; imo++) {
      for (jmo=0; jmo<n; jmo++) {
	 dsmo[imo][jmo]=transMat(dsao,c,n,imo,jmo); 
      }
   }
}

static double transMat(double **mat, double **c, int n, int lmo, int mmo)
{
   double t=0.0;
   int iao, jao; 

   for (iao=0; iao<n; iao++) {
      for (jao=0; jao<n; jao++) { t+=c[iao][lmo]*mat[iao][jao]*c[jao][mmo]; }
   }
   return t; 
}

static void loadMatData(FILE *fp, double **c, int n)
{
   int i; 

   for (i=0; i<n; i++) {
      char s[NBUF]={'\0'}, *sp; 
      int j=0; 

      fgets(s,NBUF-1,fp); 
      for (sp=s; ; sp=NULL) {
 	 char *tp=strtok(sp," "); 

	 if (tp==NULL) { break; }
	 c[i][j++]=atof(tp); 
      }
   }
} 

void loadNuData(struct input_t in, FILE **fp, struct stat_t *stat)
{
   loadNuGeom(in.nat,fp[GEOM],stat->nu.c); 
   loadNuVel(in.nat,fp[VEL],stat->nu.v); 
}

static void loadNuGeom(int na, FILE *fp, double **c)
{
   int i; 

   for (i=0; i<na; i++) {
      char s[NBUF]={'\0'}, *sp; 
      int ixyz=0; 

      fgets(s,NBUF-1,fp); 
      if (feof(fp)) { break; }

      for (sp=s; ; sp=NULL) {
	 char *tp=strtok(sp," "); 

 	 if (tp==NULL) { break; }
	    c[ixyz++][i]=atof(tp); 
      }
   }
}

static void loadNuVel(int na, FILE *fp, double **v)
{
   int i; 

   for (i=0; i<na; i++) {
      char s[NBUF]={'\0'}, *sp; 
      int ixyz=0; 

      fgets(s,NBUF-1,fp); 
      if (feof(fp)) { break; }

      for (sp=s; ; sp=NULL) {
         char *tp=strtok(sp," "); 

	 if (tp==NULL) { break; }
	 v[ixyz++][i]=atof(tp); 
      }
   }
   transToAu(v,na);
}
