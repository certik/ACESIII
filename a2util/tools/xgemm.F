
cYAU - WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING
c
c    The original xgemm filtered lda, ldb, and ldc through max(1,ld?).
c This subtly changes the behavior of dgemm depending on the other arguments.
c I refused to continue that tradition. There is quite a bit of debugging
c info in this routine which should help in finding the places where the
c leading dimensions are not properly checked.

      subroutine xgemm(transa,transb,
     &                 m,n,k,alpha,a,lda,
     &                             b,ldb,
     &                       beta, c,ldc)

c ARGUMENT LIST
      character*1 transa, transb
      integer m, n, k, lda, ldb, ldc
      double precision alpha, beta
c      double precision a(lda,*), b(ldb,*), c(ldc,*)
      double precision a, b, c

c EXTERNAL FUNCTIONS

c INTERNAL VARIABLES
c Be careful here. If you build your own generic gemm, then make sure
c you use the same integer type. Otherwise, if you link against
c system-BLAS, make sure you cast the arguments correctly.

#ifdef _32BIT_BLAS
      integer*4 int_m, int_n, int_k, int_lda, int_ldb, int_ldc
#else
      integer   int_m, int_n, int_k, int_lda, int_ldb, int_ldc
#endif /* _32BIT_BLAS */

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

#ifdef _DEBUG
      if (.not.(
     &             (
     &                  (transa.eq.'N').or.(transa.eq.'n').or.
     &                  (transa.eq.'T').or.(transa.eq.'t').or.
     &                  (transa.eq.'C').or.(transa.eq.'c')
     &             )
     &             .and.
     &             (
     &                  (transb.eq.'N').or.(transb.eq.'n').or.
     &                  (transb.eq.'T').or.(transb.eq.'t').or.
     &                  (transb.eq.'C').or.(transb.eq.'c')
     &             )
     &         )
     &   ) then
         print *, '@XGEMM: TRANSA and TRANSB must be N, T, or C.'
         print *, '        TRANSA : ', TRANSA, ', TRANSB : ', TRANSB
         call c_exit(1)
      end if
      if (min(lda,ldb,ldc).eq.0) then
         print *, '@XGEMM: A leading dimension cannot be zero.'
         print *, '        M = ',M,', N = ',N,', K = ',K
         print *, '        LDA = ',LDA,', LDB = ',LDB,', LDC = ',LDC
c      o Trap this with a debugger. Hopefully, the optimizer will fail
c        to recognize the denominator is definitely 0.
         lda=1/lda
         ldb=1/ldb
         ldc=1/ldc
         call c_exit(1)
      end if
#endif

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

