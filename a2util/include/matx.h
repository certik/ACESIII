
#ifndef _MATX_H_
#define _MATX_H_ /* MATRIX HANDLES
 *
 * A block-diagonal symmetric matrix can be stored in several different forms.
 * These include:
 *   full:  All elements (including zeroes) are stored.
 *   sqr :  Each square block is stored (both upper and lower triangles), but
 *          none of the zeros.
 *   tri :  Only the upper triangular part of each block is stored.
 *
 * An individual block can be referenced in one of several ways including:
 *   triblk :  One of the triangular blocks.
 *   triseg :  A specific block in a "tri" matrix.
 *   sqrblk :  One of the square blocks.
 *   sqrseg :  A specific block in a "sqr" matrix.
 *   fullseg:  A specific block in a "full" matrix (the zeroes are ignored).
 *
 * To make working with these matrices and blocks easier, the following
 * constants are defined:
 */

#define MAT_FULL     0
#define MAT_SQR      1
#define MAT_TRI      2

#define MAT_FULLSEG  10
#define MAT_SQRBLK   20
#define MAT_SQRSEG   21
#define MAT_TRIBLK   30
#define MAT_TRISEG   31

#endif /* _MATX_H_ */

