
      program putcrec
      implicit none

c ARGUMENTS
      integer iLength
      character*(12) szSrc

c INTERNAL VARIABLES
      integer i, iNdx, nLeft, iBuf(3), iRecLen

c ----------------------------------------------------------------------

      iLength = 11
      szSrc   = 'Hello World'

c   o iNdx is the moving szSrc index
      iNdx = 1

c   o iRecLen is the number of whole integers that contain szSrc(1:iLength)
c     (residual handled later)
      iRecLen = iLength/4

c   o the format definition for reading whole integers from szSrc
 10   format(a4)

c   o read from szSrc into iBuf by whole integers
      if (iRecLen.gt.0) then
         do i = 1, iRecLen
            read(szSrc(iNdx:iNdx+4-1),10) iBuf(i)
            iNdx = iNdx + 4
         end do
      end if

c   o pick up the slack (while flushing the integer with spaces)
      nLeft = iand(iLength,4-1)
      if (nLeft.ne.0) then
         iRecLen = iRecLen + 1
         read(szSrc(iNdx:iNdx+nLeft-1),10) iBuf(iRecLen)
      end if

      write(*,*) (iBuf(i),i=1,iRecLen)

      end

