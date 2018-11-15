enum xyz_t  {X, Y, Z, NXYZ}; 
enum file_t {MO, DSMO, CPHF, GEOM, VEL, FTOT}; 
enum mem_t  {NBUF=4096}; 

struct estat_t {
   double **cao;
   double ***u[NXYZ]; 
   double ***dsmo[NXYZ]; 
}; 

struct nstat_t {
   int *elem; 
   double *c[NXYZ];
   double *v[NXYZ]; 
};

struct stat_t {
   struct estat_t el; 
   struct nstat_t nu; 
};

struct space_t {
   double o[NXYZ]; 
   double d[NXYZ]; 
   double w[NXYZ]; 
}; 

struct input_t {
   char type[NBUF]; 
   struct space_t s; 
   int hmo; 
   int nmo; 
   int nat; 
}; 

void printInpData(struct input_t in);
void prtMat(double **mat, int m, int n); 
void printStat(struct stat_t s, int na, int nao);
