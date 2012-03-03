
c This transforms a square block of a matrix to a triangular one
      subroutine mat_trans_sqr_tri(sqr,tri,sqrlen,trilen,sqroff)
#include <aces.h>
      implicit none
#include <callstack.com>
      integer trilen,sqrlen,sqroff
      M_REAL tri(trilen),sqr(sqrlen,sqrlen)

      integer i,j,itri,len

      callstack_curr='MAT_TRANS_SQR_TRI'
      len=int(sqrt(real(1+8*trilen))-1)/2
      itri=0
      do i=sqroff,sqroff+len-1
        do j=sqroff,i
          itri=itri+1
          tri(itri)=sqr(j,i)
        end do
      end do
      return
      end
