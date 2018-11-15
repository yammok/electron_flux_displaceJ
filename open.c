#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "moparam.h"
#include "err.h"

void openMOFiles(FILE **fp, char *iname)
{
   int i; 

   for (i=MO; i<FTOT; i++) {
      char s[NBUF]={'\0'}; 

      strncpy(s,iname,NBUF-1); 
      switch (i) {
      case MO: 
 	 strncat(s,".mo",strlen(".mo")); 
	 break; 
      case DSMO: 
	 strncat(s,".dsmx",strlen(".dsmx")); 
	 break; 
      case CPHF: 
         strncat(s,".cpmx",strlen(".cpmx")); 
	 break; 
      case GEOM: 
	 strncat(s,".geom",strlen(".geom")); 
	 break; 
      case VEL: 
         strncat(s,".vel",strlen(".vel")); 
	 break; 
      default:
         ERRMSG("Undefined file detected");  
	 exit(1); 
      }

      fp[i]=fopen(s,"r"); 
      if (fp[i]==NULL) { 
	 ERRMSG("Failure opening data files"); 
	 goto fail; 
      }
   }

   return ; 
fail:
   while ((--i)>=0) { fclose(fp[i]); }
   exit(1); 
}
