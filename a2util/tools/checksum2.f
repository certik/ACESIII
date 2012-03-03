
c This routine prints the sum of elements a(1:n) and their L2 vector norm.
c It also dots them into a sawtooth waveform and prints that product.

      subroutine checksum2(sTag,a,n)
      implicit none

      character*(*) sTag
      integer n
      double precision a(n)

      double precision sm, nm, sw, sc, rs
      integer i

      sm=0.d0
      nm=0.d0
      sw=0.d0
      rs=1.d0/n
      sc=rs
      do i=1,n
         sm=sm+a(i)
         nm=nm+a(i)*a(i)
         sw=sw+a(i)*sc
         sc=sc+rs
      end do
      write(*,'(2a,1x,3e20.14)') '@CHECKSUM2: ',sTag,sm,sqrt(nm),sw

      return
      end

