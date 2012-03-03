
/*
 * This routine allocates an array of fiInts Fortran integers and
 * passes back l0 (that's ell-zero), the Fortran index of the first
 * usable address from the anchor address lCore.
 *
 * If malloc fails, the caller should test if (icore(1).eq.0), in which
 * case *fiInts will be unchanged (for sensible error messages). Testing
 * the value of l0 is NOT recommended since the index could be any value.
 *
 * WARNINGS:
 *
 *  - l0 may be negative.
 *
 *  - lCore(l0) does not point to the first address returned by malloc.
 * The first usable address is the first element of the next cache line.
 * To override this (at compile time), simply set the CACHELINE preprocessor
 * directive to some other power of 2 (presumably 8 for double boundaries).
 *
 *  - fiInts is updated to contain the actual number of Fortran integers
 * addressable at lCore(l0). The length of one cache line is added to
 * the byte-length to ensure at least fiInts integers are allocated.
 *
 *  - The value of lCore is updated to contain the actual address of the
 * new array; thus,
 *       call c_free(lCore)
 * does exactly what it looks like it should do provided:
 *  1) c_free dereferences lCore, which is the current implementation
 *  2) the address can fit in 32 bits in a 32-bit binary
 */

#include <unistd.h> /* for NULL   */
#include <stdlib.h> /* for malloc */
#include <stdio.h>  /* for printf */
#ifdef _USE_SBRK
#define GET_HEAP sbrk
#else
#define GET_HEAP malloc
#endif

#include "f77_name.h"
#include "f_types.h"

#ifndef CACHELINE
#   define CACHELINE 256 /* in Bytes */
#endif

void
F77_NAME(aces_malloc,ACES_MALLOC)
(f_int * fiInts, f_int * lCore, f_int * l0)
{
    const int iIntLn = sizeof(f_int);
#ifdef _PTRS_ARE_WORDS
    const int iPtrInt = 1;
#else
    const int iPtrInt = iIntLn;
#endif
    size_t Bytes;
    void * pHeap;

    if (*fiInts <= 0)
    {
     /* define lCore itself (leave *fiInts) */
        *lCore = 0; *l0 = 1;
        return;
    }

    Bytes = (size_t)(*fiInts*iIntLn + CACHELINE);
    if (Bytes <= 0)
    {
        printf("@ACES_MALLOC: size_t overflow, reduce memory request\n");
        exit(1);
    }

    pHeap = GET_HEAP(Bytes);
    if (pHeap)
    {
        const int iTmp = CACHELINE*iPtrInt/iIntLn;
        size_t zPtr = (size_t)pHeap+(iTmp-1);
        zPtr /= iTmp; /* the number of lines beneath the heap's address */
        zPtr *= iTmp; /* the address of the first whole line in the heap */

        *lCore = (f_int)pHeap; /* assign the heap address for freeing */
        *l0 = 1+((f_int *)zPtr-lCore); /* Fortran indexing */
        {
            size_t zMemSize = Bytes/iIntLn;
            zMemSize -= (zPtr-(size_t)pHeap+iPtrInt-1)/iPtrInt;
#ifdef _VERBOSE
            printf("@ACES_MALLOC: allocated >= %li MB of core memory\n",
                   (zMemSize*iIntLn)>>20);
#endif
            *fiInts = zMemSize; /* the number of f_ints at lCore(l0) */
        }

     /* printf("@ACES_MALLOC: %li bytes at %p\n",Bytes,pHeap); */

     /* make sure an f_int can index the first and last integer in the heap */
        {
            size_t first = 1+(zPtr-(size_t)lCore)/iPtrInt;
            size_t last  = first+(size_t)*fiInts-1;
            f_int  fi0   = first;
            f_int  fi1   = last;
            if ((first != (size_t)fi0) || (last != (size_t)fi1))
            {
                printf("\n@ACES_MALLOC: "
                       "The default integer cannot address the heap.\n"
                       "lCore = %p (%li)\n"
                       "&heap = %p (%li)\n\n",
                        (void *)*lCore,(long)*lCore,
                                 pHeap,(long) pHeap
                      );
                free(pHeap);
                exit(1);
            }
        }

        if ((void *)*lCore != pHeap)
        {
            printf("\n@ACES_MALLOC: "
                   "WARNING - the lCore address was overflowed.\n"
                   "              This memory cannot be freed.\n"
                   "lCore = %p (%li)\n"
                   "&heap = %p (%li)\n\n",
                    (void *)*lCore,(long)*lCore,
                             pHeap,(long) pHeap
                  );
        }
    }
    else
    {
     /* define lCore itself (leave *fiInts) */
        *lCore = 0; *l0 = 1;
    }
    return;
}

