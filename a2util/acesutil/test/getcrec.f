
      program getcrec
      implicit none

c ARGUMENTS
      integer iLength
      character*(11) szDest

c INTERNAL VARIABLES
      integer i, iNdx, nLeft, iBuf(3), iRecLen

c ----------------------------------------------------------------------

c     This should read 'Hello World '
      iLength = 11
      iBuf(1) = 1819043144
      iBuf(2) = 1867980911
      iBuf(3) =  543452274

c   o iNdx is the moving szDest index
      iNdx = 1

c   o iRecLen is the number of whole integers that contain szDest(1:iLength)
c     (residual handled later)
      iRecLen = iLength/4

c   o the format definition for reading whole integers into szDest
 10   format(a4)

c   o read from iBuf into szDest by whole integers
      if (iRecLen.gt.0) then
         do i = 1, iRecLen
            write(szDest(iNdx:iNdx+4-1),10) iBuf(i)
            iNdx = iNdx + 4
         end do
      end if

c   o pick up the slack
      nLeft = iand(iLength,4-1)
      if (nLeft.ne.0) then
         iRecLen = iRecLen + 1
         do i = 1, nLeft
            write(szDest(iNdx:iNdx),'(a1)') iBuf(iRecLen)
            iBuf(iRecLen) = rshift(iBuf(iRecLen),8)
            iNdx = iNdx + 1
         end do
      end if

      print *, '->',szDest(1:iLength),'<-'

c     end subroutine getcrec
      end

