
C Shell to drive LAPACK algorithm.

C A    - REAL SYMMETRIC MATRIX TO BE DIAGONALIZED (E.VALUES IN DIAGS AFTERWARDS)
C B    - EIGENVECTORS RETURNED IN COLUMNS
C JUNK - not used
C N    - SIZE OF MATRIX
C SORT - EIGENVECTORS AND EIGENVALUES ARE REORDERED, with eigenvalues:
c      0     - ascending
c      1     - unordered
c      other - descending

#include "aces.h"

      subroutine eig(a,b,junk,n,sort)
      implicit none

      integer junk, n, sort
      double precision a(n,n),b(n,n)

      integer i, j, k, index, ierr, lda
      double precision dtmp, reserve(9)

      double precision tolerance
      parameter (tolerance = 1.0d-5)

c ----------------------------------------------------------------------

      if (n.lt.1) return
      if (n.eq.1) then
         b(1,1) = 1.0d0
         return
      end if

      do j = 1, n
         do i = 1, n
            b(i,j) = a(i,j)
         end do
      end do
#ifdef _ASSERT
      dtmp = 0.d0
      do j = 2, n
         do i = 1, j-1
            dtmp = dtmp + dabs(a(j,i)-a(i,j))
         end do
      end do
      if (dtmp.gt.1.d-10) then
         print *, '@EIG: WARNING - the matrix is not symmetric'
         print *, '      sum of differences is ',dtmp
      end if
#endif
      ierr=0
      if (n.lt.4) then
         call B_SYEV('V','L',n,b,n,a,reserve,9,      ierr)
      else
         call B_SYEV('V','L',n,b,n,a,a(1,2), (n*n-n),ierr)
      end if
      if (ierr.ne.0) then
         write(*,*) '@EIG: Eigenvalue not found. Error code: ',ierr
         call errex
      end if

c   o make first significant element of each eigenvector positive
      do j = 1, n
         index = 0
         do i = 1, n
            if (index.eq.0) then
               if (abs(b(i,j)).gt.tolerance) index = i
            end if
         end do
         if (index.ne.0) then
            if (b(index,j).lt.0.0d0) then
               do i = 1, n
                  b(i,j) = -b(i,j)
               end do
            end if
         end if
      end do

      if ((sort.ne.0).and.(sort.ne.1)) then
         do i = 1, n-1
            do j = i+1, n
               if ( a(i,1) .lt. a(j,1) ) then
                  dtmp   = a(i,1)
                  a(i,1) = a(j,1)
                  a(j,1) = dtmp
                  do k = 1, n
                     dtmp   = b(k,i)
                     b(k,i) = b(k,j)
                     b(k,j) = dtmp
                  end do
               end if
            end do
         end do
      end if

c   o put the eigenvalues along the diagonal of A and zero the rest
      do j = 2, n
         do i = 1, n
            a(i,j) = 0.0d0
         end do
      end do
      do i = 1, n
         a(i,i) = a(i,1)
      end do
      do i = 2, n
         a(i,1) = 0.0d0
      end do

      return
      end

