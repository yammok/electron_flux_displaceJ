#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "moparam.h"

#define NARGS 3

void loadInpData(struct input_t *in, char *iname, char *type); 
void allocMOArray(struct stat_t *s, struct input_t in); 
void initMOArray(struct stat_t *s , struct input_t in); 
void openMOFiles(FILE **fp, char *iname);
void loadElData(struct input_t in, FILE **fp, struct stat_t *stat); 
void loadNuData(struct input_t i, FILE **fp, struct stat_t *stat); 
void calcJSpDist(struct input_t i, struct stat_t s); 
void calcDSpDist(struct input_t i, struct stat_t s); 
void printPoints(struct space_t s, char *iname);
void printInpData(struct input_t in); 
void printGrdMat(double ***mx[NXYZ], int na, int m, int n); 

static void usage(char *cmd); 

int main(int argc, char *argv[])
{
   char inpname[NBUF]={'\0'}; 
   struct input_t inp; 
   struct stat_t  stat; 
   FILE *fp[FTOT]; 

   if (argc!=NARGS) {
      usage(argv[0]); 
      return 1; 
   }
   else { strncpy(inpname,argv[2],NBUF-1); }

   loadInpData(&inp,inpname,argv[1]); 
   printPoints(inp.s,inpname); 
   allocMOArray(&stat,inp); 
   initMOArray(&stat,inp); 
   openMOFiles(fp,inpname); 

   while (1) {
      loadElData(inp,fp,&stat); 
      loadNuData(inp,fp,&stat); 
      if (feof(fp[MO])) { break; }
      if      (strcmp(inp.type,"flux")==0)    { calcJSpDist(inp,stat); }
      else if (strcmp(inp.type,"density")==0) { calcDSpDist(inp,stat); }
   }

   return 0; 
}

static void usage(char *cmd)
{
   fprintf(stderr,"Usage : %s [option] [cur file]\n",&cmd[2]);
   fprintf(stderr,"Option:\n"); 
   fprintf(stderr,"	-d: one electron density calculation\n"); 
   fprintf(stderr,"	-f: one electron flux density calculation\n"); 
}
