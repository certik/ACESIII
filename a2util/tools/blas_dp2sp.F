#ifdef USE_SP_BLAS
      subroutine daxpy(n,da,dx,incx,dy,incy)
      double precision dx(*),dy(*),da
      integer incx,incy,n
      call saxpy(n,da,dx,incx,dy,incy)
      return
      end

      subroutine dcopy(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
      integer incx,incy,n
      call scopy(n,dx,incx,dy,incy)
      return
      end

      subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,
     &                 work,lwork,info)
      character jobvl,jobvr
      integer info,lda,ldvl,ldvr,lwork,n
      double precision a(lda,*),vl(ldvl,*),vr(ldvr,*),wi(*)
      double precision work(*),wr(*)
      call sgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,
     &           work,lwork,info)
      return
      end

      subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      character transa,transb
      integer lda,ldb,ldc,m,n,k
      double precision alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
      call sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      return
      end

      subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      character trans
      integer incx,incy,lda,m,n
      double precision alpha,beta,a(lda,*),x(*),y(*)
      call sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      return
      end

      subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      integer lda,ldb,n,nrhs,ipiv(*)
      double precision a(lda,*),b(ldb,*)
      call sgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      return
      end

      subroutine dgetrf(m,n,a,lda,ipiv,info)
      integer info,lda,m,n,ipiv(*)
      double precision a(lda,*)
      call sgetrf(m,n,a,lda,ipiv,info)
      return
      end

      subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      integer info,lda,lwork,n,ipiv(*)
      double precision a(lda,*),work(*)
      call sgetri(n,a,lda,ipiv,work,lwork,info)
      return
      end

      subroutine drot(n,dx,incx,dy,incy,c,s)
      double precision dx(*),dy(*),c,s
      integer incx,incy,n
      call srot(n,dx,incx,dy,incy,c,s)
      return
      end

      subroutine dscal(n,da,dx,incx)
      double precision da,dx(*)
      integer incx,n
      call sscal(n,da,dx,incx)
      return
      end

      subroutine dswap(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
      integer incx,incy,n
      call sswap(n,dx,incx,dy,incy)
      return
      end

      subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      character jobz,uplo
      integer info,lda,lwork,n
      double precision a(lda,*),w(*),work(*)
      call ssyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      return
      end

      subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
      character side,uplo
      integer lda,ldb,ldc,m,n
      double precision alpha,beta,a(lda,*),b(ldb,*),c(ldc,*)
      call ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
      return
      end

      subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
      character uplo
      integer lda,n,incx,incy
      double precision alpha,beta,a(lda,*),x(*),y(*)
      call ssymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
      return
      end

      subroutine dsyr(uplo,n,alpha,x,incx,a,lda)
      character uplo
      integer incx,lda,n
      double precision alpha,a(lda,*),x(*)
      call ssyr(uplo,n,alpha,x,incx,a,lda)
      return
      end

      subroutine dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
      character uplo
      integer incx,incy,lda,n
      double precision alpha,a(lda,*),x(*),y(*)
      call ssyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
      return
      end

      double precision function dasum(n,dx,incx)
      double precision dx(*),sasum
      integer incx,n
      dasum=sasum(n,dx,incx)
      return
      end

      double precision function ddot(n,dx,incx,dy,incy)
      double precision dx(*),dy(*),sdot
      integer incx,incy,n
      ddot=sdot(n,dx,incx,dy,incy)
      return
      end

      double precision function dlamch(cmach)
      character cmach
      double precision slamch
      dlamch=slamch(cmach)
      return
      end

      double precision function dnrm2(n,x,incx)
      double precision x(*),snrm2
      integer incx,n
      dnrm2=snrm2(n,x,incx)
      return
      end

      double precision function dsum(n,dx,incx)
      double precision dx(*),ssum
      integer incx,n
      dsum=ssum(n,dx,incx)
      return
      end

      integer function idamax(n,dx,incx)
      double precision dx(*)
      integer incx,n,isamax
      idamax=isamax(n,dx,incx)
      return
      end
#else
#include "return.f"
#endif /* USE_SP_BLAS */
