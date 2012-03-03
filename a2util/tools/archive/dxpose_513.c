
/* THERE IS A LOGIC ERROR SOMEWHERE IN HERE! TEST A 19x4 MATRIX. */

/*
 * This routine transposes a rectangular matrix in situ.
 * The algorithm was blatantly ripped from netlib/TOMS.
 * (cf. http://www.netlib.org/toms/513)
 *
 * ARGUMENTS:
 *   a(,) - the double precision array to transpose
 *   m    - THE 'FASTEST-RUNNING' INDEX
 *          [In Fortran, this is the number of rows in a(,), but
 *           in C, this is the number of columns.]
 *   n    - THE 'SLOWEST-RUNNING' INDEX
 *          [In Fortran, this is the number of columns in a(,), but
 *           in C, this is the number of rows.]
 */

#include <stdio.h>
#include <stdlib.h>

#include "f_types.h"

#define WORK unsigned char /* this MUST be 8 bits */
#define bool short
#define TRUE  1
#define FALSE 0

void
#ifdef C_SUFFIX
  dxpose_513_
#else
  dxpose_513
#endif /* C_SUFFIX */
(double * a, f_int * m, f_int * n)
{

    if ((*m < 2) || (*n < 2)) return;

    if (*m == *n)
    {
        long j;
        for ( j=0; j<(*m-1); j++)
        {
            long i;
            for ( i=(j+1); i<*m; i++)
            {
                double b      = *(a+i+(*m*j));
                *(a+i+(*m*j)) = *(a+j+(*m*i));
                *(a+j+(*m*i)) = b;
            }
        }
    }
    else
    {

        WORK *work; long workdim;

        long ncount;

        long mn, k, i, im;

        {
            size_t size;
            /* workdim = ((*m)+(*n))/2; */
            workdim = (*m)*(*n)/2;
            size = 1 + (workdim/8);
            work = malloc(size);
            if (work == NULL)
            {
                printf("@dmat_xpose: Not enough free system memory.\n");
                exit(-1);
            }
            else
            {
                int z;
                for ( z=0; z<size; z++) *(work+z) &= 0;
                workdim = size*8;
            }
        }

        if ((*m != 2) && (*n != 2))
        {
            /*
             * Calculate the number of fixed points
             * using Euclid's algorithm for GCD(m-1,n-1).
             */
            long ir2  = *m-1;
            long ir1  = *n-1;
            long itmp = 1;
            while (itmp != 0)
            {
                itmp = ir2 % ir1;
                ir2  = ir1;
                ir1  = itmp;
            }
            ncount = ir2+1;
        }
        else
        {
            ncount = 2;
        }

        mn = (*m)*(*n);
        k  = mn-1;
        i  = 1;
        im = *m;
        while (ncount < mn)
        {
            long i1, i1c, kmi;
            double b, c;
            bool still_going = TRUE;

            i1  = i;
            kmi = k-i;
            b   = *(a+i1);
            i1c = kmi;
            c   = *(a+i1c);
            while (still_going)
            {
                long i2, i2c;

                /* leave this line as is since i1/(*n) returns the number
                 * of times *n wholly divides i1
                 */
                i2  = *m*i1 - k*(i1/(*n));
                i2c = k-i2;
                if (i1  <= workdim) *(work+(i1 /8)) |= (1 << (i1  % 8));
                if (i1c <= workdim) *(work+(i1c/8)) |= (1 << (i1c % 8));
                ncount += 2;
                if (i2 == i)
                {
                    *(a+i1)  = b;
                    *(a+i1c) = c;
                    still_going = FALSE;
                }
                else
                {
                    if (i2 == kmi)
                    {
                        double d = b;
                        b        = c;
                        c        = d;
                        *(a+i1)  = b;
                        *(a+i1c) = c;
                        still_going = FALSE;
                    }
                    else
                    {
                        *(a+i1)  = *(a+i2);
                        *(a+i1c) = *(a+i2c);
                        i1  = i2;
                        i1c = i2c;
                    } /* if (i2 == kmi) */
                } /* if (i2 == i) */
            } /* while (still_going) */
            if (ncount < mn)
            {
                still_going = TRUE;
                while (still_going)
                {
                    long i2, i2c, max;

                    max = k-i;
                    i += 1;
                    if (i > max)
                    {
                        printf("\n"
"@dmat_xpose: Congratulations!\n"
"   You have found a circumstance in which the default algorithm for sizing\n"
"a work array has failed. Unfortunately, there is nothing you can do about\n"
"this besides altering the code, recompiling, and reinstalling.\n"
"      i = %D, max = %D\n"
"      workdim = %D\n"
"      mxn = %Dx%D\n",
i,max,workdim,*m,*n
                        );
                        exit(-1);
                    }
                    if (im > k) im -= k;
                    else        im += *m;
                    i2 = im;
                    if (i != i2)
                    {
                        if (i > workdim)
                        {
                            while ((i < i2) && (i2 < max))
                            {
                                i1 = i2;
                                i2 = *m*i1 - k*(i1/(*n));
                            }
                            still_going = (i2 != i);
                        }
                        else
                        {
                            still_going = ((*(work+(i/8)) >> (i % 8)) & 1);
                        }
                    } /* if (i != i2) */
                } /* while (still_going) */
            } /* if (ncount < mn) */
        } /* while (ncount < mn) */

        free(work);
    } /* if (*m == *n) */

    return;
}

