
      subroutine dmat_prt(dmat,m,n,ld)
      implicit none
      double precision dmat(*)
      integer m, n, ld

      integer i, j, start, cols, doing

      if ((m.eq.0).or.(n.eq.0).or.(ld.lt.m)) return

      start = 1
      cols = n
      do while (cols.ne.0)
         doing = min(3,cols-1)
         do i = 0, m-1
c            print '(4F19.14)', (dmat(start+i+j*ld),j=0,doing)
            print '(4F19.12)', (dmat(start+i+j*ld),j=0,doing)
         end do
         print '(/)'
         doing = doing + 1
         cols = cols - doing
         start = start + ( ld * doing )
c     end do while (cols.ne.0)
      end do

      return
      end

