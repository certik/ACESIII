
c This transforms a square block of a matrix
      subroutine mat_trans_sqr_sqr(in,out,inlen,outlen,inoff,outoff)
#include <aces.h>
      implicit none
#include <callstack.com>
      integer inoff,outoff,inlen,outlen
      M_REAL in(inlen,inlen),out(outlen,outlen)

      integer len,i

      callstack_curr='MAT_TRANS_SQR_SQR'
      len=min(inlen,outlen)
      do i=1,len
        call B_COPY(len,in(inoff,inoff+i-1),1,out(outoff,outoff+i-1),1)
      enddo
      return
      end
