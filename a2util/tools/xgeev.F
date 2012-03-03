
cYAU - WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING
c
c    For some reason, which I am not privy to, the original xgeev was not bound
c the same way dgemm is. I moved the original one to libr/mn_geev.f

      subroutine xgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr,
     &                  ldvr, work, lwork, info )
      implicit none

      character*1 jobvl, jobvr
      integer info, lda, ldvl, ldvr, lwork, n
      double precision a(lda,*), vl(ldvl,*), vr(ldvr,*)
      double precision wi(*), work(*), wr(*)

#ifdef _32BIT_BLAS
      integer*4 int_n, int_lda, int_ldvl, int_ldvr, int_lwork, int_info
#else
      integer   int_n, int_lda, int_ldvl, int_ldvr, int_lwork, int_info
#endif /* _32BIT_BLAS */

      if (n.eq.0) return

      int_n     = n
      int_lda   = lda
      int_ldvl  = ldvl
      int_ldvr  = ldvr
      int_lwork = lwork

#ifdef USE_SP_BLAS
      call sgeev(jobvl,jobvr,int_n,a,int_lda,wr,wi,
     &           vl,int_ldvl,vr,int_ldvr,work,int_lwork,int_info)
#else
      call dgeev(jobvl,jobvr,int_n,a,int_lda,wr,wi,
     &           vl,int_ldvl,vr,int_ldvr,work,int_lwork,int_info)
#endif

      info = int_info

      return
      end

