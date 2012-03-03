
c This transforms a triangular block of a matrix to a square one
      subroutine mat_trans_tri_sqr(tri,sqr,trilen,sqrlen,sqroff)
#include <aces.h>
      implicit none
#include <callstack.com>
      integer trilen,sqrlen,sqroff
      M_REAL tri(trilen),sqr(sqrlen,sqrlen)

      integer i,j,itri,len

      callstack_curr='MAT_TRANS_TRI_SQR'
      len=int(sqrt(real(1+8*trilen))-1)/2
      itri=0
      do i=sqroff,sqroff+len-1
        do j=sqroff,i
          itri=itri+1
          sqr(i,j)=tri(itri)
          sqr(j,i)=tri(itri)
        end do
      end do
      return
      end
