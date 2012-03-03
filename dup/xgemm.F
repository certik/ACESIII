C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine xgemm(transa,transb,
     &                 m,n,k,alpha,a,lda,
     &                             b,ldb,
     &                       beta, c,ldc)

c----------------------------------------------------------------------------
c  This code overrides the normal ACESII serial version of xgemm in the 
c  "tools" library.  Its purpose is to provide a version of xgemm that works 
c  has the proper interface to the calling sequence required in the joda
c  routines, while calling the actual system-level dgemm routine with the
c  properly typed arguments.
c
c  This version of xgemm gets linked in ahead of the tools/xgemm.o.  Then
c  the joda code within the xaces3 executable calls this version of xgemm,
c  which calls the proper system-level dgemm routine.  
c----------------------------------------------------------------------------
 
c ARGUMENT LIST
      character*1 transa, transb
#ifdef ACES64
      integer*8 m, n, k, lda, ldb, ldc
#else
      integer m, n, k, lda, ldb, ldc
#endif
      double precision alpha, beta
c      double precision a(lda,*), b(ldb,*), c(ldc,*)
      double precision a, b, c

#ifdef __crayx1
      integer*8   int_m, int_n, int_k, int_lda, int_ldb, int_ldc
#else
      integer   int_m, int_n, int_k, int_lda, int_ldb, int_ldc
#endif

c ----------------------------------------------------------------------

c This quick return code was taken straight from blas/dgemm. We need it
c for cases when m=lda=0 (e.g., the number of occupied orbitals is zero
c but we didn't check for it). Rather than have dgemm crash, we just
c want to return having done nothing. This also saves us from filtering
c lda, ldb, and ldc through max(1,ld?).

c BEWARE!!! xgemm will, therefore, not crash when the other arguments
c are ill-conditioned.

      if ((m.eq.0).or.
     &    (n.eq.0).or.
     &    (((alpha.eq.(0.0d0)).or.(k.eq.0)).and.(beta.eq.(1.0d0)))
     &   ) return

c ----------------------------------------------------------------------

c   o recast the arguments
      int_m   = m
      int_n   = n
      int_k   = k
      int_lda = lda
      int_ldb = ldb
      int_ldc = ldc

#ifdef USE_SP_BLAS
      call sgemm(transa,transb,int_m,int_n,int_k,
     &           alpha,a,int_lda,
     &                 b,int_ldb,
     &           beta, c,int_ldc)
#else
      call dgemm(transa,transb,int_m,int_n,int_k,
     &           alpha,a,int_lda,
     &                 b,int_ldb,
     &           beta, c,int_ldc)
#endif

      return
      end

