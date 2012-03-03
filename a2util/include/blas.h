
#ifndef _BLAS_H_
#define _BLAS_H_ /* BLAS DEFINES
 *
 * The following definitions define the blas routines to call.  It is done
 * to distinguish between single and double precision routines.
 *
 * Instead of calling any blas routine by name, call it using the appropriate
 * B_ macro below.  For the following lines:
 *     call scopy(w,x,incx,y,incy)
 *     call dcopy(w,x,incx,y,incy)
 * should both be replaced with:
 *     call B_COPY(w,x,incx,y,incy)
 *
 * All macros are of the form B_NAME where NAME is the name of the routine
 * (minus the leading S or D).  The only exception is that the functions
 * isamax and idamax have the macro B_AMAX.
 */

#ifdef USE_SP_BLAS

#  define B_ROTG   srotg /* L1 BLAS */
#  define B_ROTMG  srotmg
#  define B_ROT    srot
#  define B_ROTM   srotm
#  define B_SWAP   sswap
#  define B_SCAL   sscal
#  define B_COPY   scopy
#  define B_AXPY   saxpy
#  define B_DOT    sdot
#  define B_WRM2   swrm2
#  define B_ASUM   sasum
#  define B_AMAX   isamax

#  define B_GEMV   sgemv /* L2 BLAS */
#  define B_GBMV   sgbmv
#  define B_SYMV   ssymv
#  define B_SBMV   ssbmv
#  define B_SPMV   sspmv
#  define B_TRMV   strmv
#  define B_TBMV   stbmv
#  define B_TPMV   stpmv
#  define B_TRSV   strsv
#  define B_TBSV   stbsv
#  define B_TPSV   stpsv
#  define B_GER    sger
#  define B_SYR    ssyr
#  define B_SPR    sspr
#  define B_SYR2   ssyr2
#  define B_SPR2   sspr2

#  define B_GEMM   sgemm /* L3 BLAS */
#  define B_SYMM   ssymm
#  define B_SYRK   ssyrk
#  define B_SYR2K  ssyr2k
#  define B_TRMM   strmm
#  define B_TRSM   strsm

#  define B_GESV   sgesv /* Lapack linear equations */
#  define B_GBSV   sgbsv
#  define B_GTSV   sgtsv
#  define B_POSV   sposv
#  define B_PPSV   sppsv
#  define B_PBSV   spbsv
#  define B_PTSV   sptsv
#  define B_SYSV   ssysv
#  define B_SPSV   sspsv

#  define B_GELS   sgels /* Lapack least squares */
#  define B_GELSS  sgelss

#  define B_SYEV   ssyev /* Lapack eigenvalue */
#  define B_SPEV   sspev
#  define B_SBEV   ssbev
#  define B_STEV   sstev
#  define B_GEES   sgees
#  define B_GEEV   sgeev
#  define B_GESVD  sgesvd

#  define B_SYGV   ssygv /* Lapack generalized eigenvalue */
#  define B_SPGV   sspgv

#  define B_GESVX   sgesvx /* Lapack linear equations (expert) */
#  define B_GBSVX   sgbsvx
#  define B_GTSVX   sgtsvx
#  define B_POSVX   sposvx
#  define B_PPSVX   sppsvx
#  define B_PBSVX   spbsvx
#  define B_PTSVX   sptsvx
#  define B_SYSVX   ssysvx
#  define B_SPSVX   sspsvx

#  define B_GELSX   sgelsx /* Lapack least squares (expert) */

#  define B_SYEVX   ssyevx /* Lapack eigenvalue (expert) */
#  define B_SPEVX   sspevx
#  define B_SBEVX   ssbevx
#  define B_STEVX   sstevx
#  define B_GEESX   sgeesx
#  define B_GEEVX   sgeevx

#  define B_GECO    sgeco /* Linpack (minv is stupid) */
#  define B_GEDI    sgedi
#  define B_GESL    sgesl

#else

#  define B_ROTG   drotg /* L1 BLAS */
#  define B_ROTMG  drotmg
#  define B_ROT    drot
#  define B_ROTM   drotm
#  define B_SWAP   dswap
#  define B_SCAL   dscal
#  define B_COPY   dcopy
#  define B_AXPY   daxpy
#  define B_DOT    ddot
#  define B_WRM2   dwrm2
#  define B_ASUM   dasum
#  define B_AMAX   idamax

#  define B_GEMV   dgemv /* L2 BLAS */
#  define B_GBMV   dgbmv
#  define B_SYMV   dsymv
#  define B_SBMV   dsbmv
#  define B_SPMV   dspmv
#  define B_TRMV   dtrmv
#  define B_TBMV   dtbmv
#  define B_TPMV   dtpmv
#  define B_TRSV   dtrsv
#  define B_TBSV   dtbsv
#  define B_TPSV   dtpsv
#  define B_GER    dger
#  define B_SYR    dsyr
#  define B_SPR    dspr
#  define B_SYR2   dsyr2
#  define B_SPR2   dspr2

#  define B_GEMM   dgemm /* L3 BLAS */
#  define B_SYMM   dsymm
#  define B_SYRK   dsyrk
#  define B_SYR2K  dsyr2k
#  define B_TRMM   dtrmm
#  define B_TRSM   dtrsm

#  define B_GESV   dgesv /* Lapack linear equations */
#  define B_GBSV   dgbsv
#  define B_GTSV   dgtsv
#  define B_POSV   dposv
#  define B_PPSV   dppsv
#  define B_PBSV   dpbsv
#  define B_PTSV   dptsv
#  define B_SYSV   dsysv
#  define B_SPSV   dspsv

#  define B_GELS   dgels /* Lapack least squares */
#  define B_GELSS  dgelss

#  define B_SYEV   dsyev /* Lapack eigenvalue */
#  define B_SPEV   dspev
#  define B_SBEV   dsbev
#  define B_STEV   dstev
#  define B_GEES   dgees
#  define B_GEEV   dgeev
#  define B_GESVD  dgesvd

#  define B_SYGV   dsygv /* Lapack generalized eigenvalue */
#  define B_SPGV   dspgv

#  define B_GESVX   dgesvx /* Lapack linear equations (expert) */
#  define B_GBSVX   dgbsvx
#  define B_GTSVX   dgtsvx
#  define B_POSVX   dposvx
#  define B_PPSVX   dppsvx
#  define B_PBSVX   dpbsvx
#  define B_PTSVX   dptsvx
#  define B_SYSVX   dsysvx
#  define B_SPSVX   dspsvx

#  define B_GELSX   dgelsx /* Lapack least squares (expert) */

#  define B_SYEVX   dsyevx /* Lapack eigenvalue (expert) */
#  define B_SPEVX   dspevx
#  define B_SBEVX   dsbevx
#  define B_STEVX   dstevx
#  define B_GEESX   dgeesx
#  define B_GEEVX   dgeevx

#  define B_GECO    dgeco /* Linpack (minv is stupid) */
#  define B_GEDI    dgedi
#  define B_GESL    dgesl

#endif

#endif /* _BLAS_H_ */

