
c This routine prints the chemical system (minus dummy atoms) in the form:
c   A1 X1 Y1 Z1
c   A2 X2 Y2 Z2

      subroutine xyz
      implicit none

      integer maxatm
      parameter(maxatm=100)

      double precision coord(3,maxatm)
      integer nAtoms, i, indx
      character*(5*maxatm) zsymuni

#include "machsp.com"

c ----------------------------------------------------------------------

      call getrec(1,'JOBARC','NATOMS',1,nAtoms)
      if (nAtoms.gt.maxatm) then
         print *, '@XYZ: Assertion failed.'
         print *, '      maxatm = ',maxatm
         print *, '      nAtoms = ',nAtoms
         call aces_exit(1)
      end if

c   o get the Cartesian coordinates (with dummy atoms)
      call getrec(1,'JOBARC','COORD',iintfp*3*nAtoms,coord)

c   o get the atomic symbols
      call getcrec(1,'JOBARC','ZSYM',5*nAtoms,zsymuni)

c   o print the matrix without dummies
      indx = 1
      do i = 1, nAtoms
         if (zsymuni(indx:indx).ne.'X') print '(a2,3(x,f12.8))',
     &      zsymuni(indx:indx+1), coord(1,i), coord(2,i), coord(3,i)
         indx = indx + 5
      end do

      return
c     end subroutine xyz
      end

