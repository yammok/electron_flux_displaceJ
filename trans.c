#include "unit.h"
#include "moparam.h"

void transToAu(double **v, int n)
{
   int i, j; 

   for (i=0; i<n; i++) {
      for (j=X; j<=Z; j++) { v[j][i]*=0.1*PLANK/EHAR; }
   }
}
