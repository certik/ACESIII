#include "f77_name.h"
#include "f_types.h"


void F77_NAME(proto1, PROTO1)(int *myrank, int *size, int *commonstart, int *end, 
           int *maxa, int *maxb, int *maxc, int *maxd, int *maxi, int *maxj, int *endfunction, int *nshells,
           int *role, double *dbuf, int *ibuf);

           
void F77_NAME(integral, INTEGRAL)
            (int *message0, int *message1, int *message2, int *message3,
             int *message4, int *message5, int *message6, int *message7,
             int *commonstart,   double *dbuf, int *sendnumber, 
             int *intbuffer, int * buffersize, double *dbuf1);
             
             
void F77_NAME(cwork, CWORK) (double *v, double *told, double *tnew, int *ni, int *nj, int *na, int *nb, int *nc, int *nd);
