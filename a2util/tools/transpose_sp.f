
c ---------------------------------------------------------------
c --- "TRANSPOSE" AN UPPER- OR LOWER- PACKED TRIANGULAR ARRAY ---
c ---                ANTHONY YAU, 29 JULY 2000                ---
c ---------------------------------------------------------------

c CALL TRANSPOSE_SP(UPLO,M,DY,IX)
c
c UPLO = { "U", "u", "L", "l" }
c        the condition of DY on entry (upper- or lower-packed)
c
c M    = { integer }
c        the order of the square array
c
c DY   = { double precision array of dimension M*(M+1)/2 }
c        the double precision array to be transposed
c
c IX   = { integer array of dimension M*(M+1)/2 }
c        a work array

      subroutine transpose_sp(uplo,m,dy,ix)
      implicit none

c ARGUMENT LIST
      character*1 uplo
      integer m
      integer ix(m*(m+1)/2)
      double precision dy(m*(m+1)/2)

c INTERNAL VARIABLES
      double precision dtmp
      integer offset, row, rowptr, col, colsum

c DISCARD THE CHAFF
      if (m.lt.0) then
         print *, "@TRANSPOSE_SP: The order of the array to",
     &            " transpose cannot"
         print *, "               be less than zero."
         call c_exit(1)
      end if
      if ((m.eq.0).or.(m.eq.1).or.(m.eq.2)) return
      if (m.eq.3) then
         dtmp =dy(3)
         dy(3)=dy(4)
         dy(4)=dtmp
         return
      end if

c GET TO WORK
c   o If m is 4 then the code is short enough to implement by hand
c     and will be much faster than treating it generally; otherwise,
c     create the integer array to be sorted. This is a list of the
c     destination addresses for the elements in dy.
      if ((uplo.eq."U").or.(uplo.eq."u")) then
         if (m.eq.4) then
            dtmp =dy(3)
            dy(3)=dy(4)
            dy(4)=dy(7)
            dy(7)=dy(8)
            dy(8)=dy(6)
            dy(6)=dy(5)
            dy(5)=dtmp
            return
         end if
         ix(1)=1
         offset=2
         do row = 2,m
            ix(offset)=row
            offset=offset+1
            colsum=0
            do col=1,(row-1)
               colsum=colsum+(m-col)
               ix(offset)=row+colsum
               offset=offset+1
            end do
         end do
      else if ((uplo.eq."L").or.(uplo.eq."l")) then
         if (m.eq.4) then
            dtmp =dy(3)
            dy(3)=dy(5)
            dy(5)=dy(6)
            dy(6)=dy(8)
            dy(8)=dy(7)
            dy(7)=dy(4)
            dy(4)=dtmp
            return
         end if
         offset=1
         do row = 1,(m-1)
            rowptr=ishft((row*(row+1)),-1)
            ix(offset)=rowptr
            offset=offset+1
            colsum=0
            do col=row,(m-1)
               colsum=colsum+col
               ix(offset)=rowptr+colsum
               offset=offset+1
            end do
         end do
         ix(offset)=ishft((m*(m+1)),-1)
      else
         print *, "@TRANSPOSE_SP: An invalid character has been",
     &            " received."
         print *, "               The character ", uplo,
     &            " should be U, u, L, or l."
         call c_exit(1)
      end if

c   o We're done with the counting variables so we can pass them as
c     arguments to idsort.
      colsum=ishft((m*(m+1)),-1)
      offset=2
      call idsort(ix,dy,colsum,offset)

      return
      end

