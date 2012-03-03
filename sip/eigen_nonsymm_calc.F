C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine eigen_nonsymm_calc(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Calculates the eigenvalues and eigenvectors of a 2-d square matrix.
c   The matrix does NOT have to be symmetric. 
c   The matrix is also diagonalized on output.
c
c   The eigenvalues are stored in a common block.  
c   The instruction is executed as follows:
c   
c   execute eig Array1 Array2
c
c   where Array1 is the input matrix that is to be diagonalized.
c         Array2 is the matrix whose columns are the eigenvectors of Array1.
c
c   Both Array1 and Array2 must be "static" arrays.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer*8 iarray, ievec, get_index_from_base
      integer i, j, n
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer ind(mx_array_index)
      integer junk

      common /epsilon/epsilon(1000), epsilonb(1000)
      double precision epsilon, epsilonb

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

c----------------------------------------------------------------------------
c   Locate the data for the input matrix.
c---------------------------------------------------------------------------

      array = op(c_result_array)
      evec_array = op(c_op1_array)
      if (array .eq. 0 .or. evec_array .eq. 0) then
         print *,'Error: EIG routine requires 2 array arguments.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      array_type = array_table(c_array_type, array)
      evec_type  = array_table(c_array_type, evec_array)
      if (array_type .ne. static_array .or. 
     *    evec_type .ne. static_array) then
         print *,'Error: Both arrays in EIG routine must be static.'
         call abort_job()
      endif

#ifdef ALTIX
      dptr = dshptr
#endif

      nindex = array_table(c_nindex, array)
      nindex_evec = array_table(c_nindex, evec_array)
      if (nindex .ne. nindex_evec) then
         print *,'Error: Both arrays in EIG must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex,' indices.'
         print *,'Second array has ',nindex_evec,' indices.'
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray = get_index_from_base(address_table(array), x, 2)
      ievec  = get_index_from_base(address_table(evec_array), x, 2)

      n = index_table(c_index_size, ind(1))  ! pick up length of index

c---------------------------------------------------------------------------
c   Diagonalization.
c---------------------------------------------------------------------------
      call diag_nonsymm(x(iarray), x(ievec), n)

c     call print_eig(x(ievec), n)

      return
      end

      subroutine diag_nonsymm(a,b,n)
      implicit none
      include 'machine_types.h'
      include 'parallel_info.h'

      integer n 
      double precision a(n,n),b(n,n)
      double precision wr(n), wi(n), z(n,n),fv1(n),x 
      integer iv1(n)
      integer i, j, ierr 
      integer*8 arg64(3), iarg64(n)
c
c      do i = 1, n 
c      do j = 1, n 
c      write(6,*)'      a(',i,',',j,')=',a(i,j) 
c      enddo 
c      enddo 
c
      do i = 1, n 
      do j = 1, n 
         b(i,j) = 0.0d0 
      enddo 
      enddo 

c ----------------------------------------------------------------------
c   Call the ACESII "eig" routine.
c-----------------------------------------------------------------------
#ifdef ACES64
      arg64(1) = n
      arg64(2) = 1  
      call rg(arg64(1),arg64(1),a,wr,wi,arg64(2),z,iarg64,fv1,arg64(3))   
      ierr = arg64(3) 
#else
      call rg(n,n,a,wr,wi,1,z,iv1,fv1,ierr)   
#endif
c
c Put the eigenvalues on the diagonal of B. 
c ----------------------------------------- 
c
c      write(6,*)'space dimension = ',n
      do i = 1, n 
         a(i,i) = wr(i) 
c      write(6,*) ' Eig of A, real: ',wr(i), 'img: ',wi(i) 
      enddo 
c
c Put the real part of the eigenvectors in A. 
c ------------------------------------------- 
c
c     write(6,*) ' Eigenvectors of diagonalized array', n, n  
      do i = 1, n 
      do j = 1, n 
         b(i,j) = z(i,j) 
      enddo 
      enddo 
c
c Calculate norm of eigenvectors
c ------------------------------------------- 
      do i = 1, n 
      x=0.0d+0
c
      do j = 1, n 
      x = x + (b(j,i) * z(j,i) )
      enddo 
c
c Normalize the eigenvectors
c ------------------------------------------- 
      do j = 1, n 
      b(j,i) = z(j,i) * (1.0d0/(dsqrt( x ))) 
      enddo 
c
      enddo 
c 
      return
      end

C ***************************************************************
C  STANDARD DIAGONALIZATION PROCEDURE FOR NON-HERMITEAN MATRICES
C ***************************************************************
C
      subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
      integer i, j 
c
#ifdef ACES64
      integer*8 n,nm,is1,is2,ierr,matz
      integer*8 iv1(n)
#else
      integer n,nm,is1,is2,ierr,matz
      integer iv1(n)
#endif
      double precision a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.  complex conjugate
c        pairs of eigenvalues appear consecutively with the
c        eigenvalue having the positive imaginary part first.
c
c        z  contains the real and imaginary parts of the eigenvectors
c        if matz is not zero.  if the j-th eigenvalue is real, the
c        j-th column of  z  contains its eigenvector.  if the j-th
c        eigenvalue is complex with positive imaginary part, the
c        j-th and (j+1)-th columns of  z  contain the real and
c        imaginary parts of its eigenvector.  the conjugate of this
c        vector is the eigenvector for the conjugate eigenvalue.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for hqr
c           and hqr2.  the normal completion code is zero.
c
c        iv1  and  fv1  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c      write(6,*)'rg entered'
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  balanc(nm,n,a,is1,is2,fv1)
      call  elmhes(nm,n,is1,is2,a,iv1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)

      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)

      if (ierr .ne. 0) go to 50
      call  balbak(nm,n,is1,is2,fv1,n,z)

c     write(6,*) ' Z in rg' 
c     do i = 1, n 
c     do j = 1, n 
c        write(6,*) i, j, z(i,j) 
c     enddo 
c     enddo 

   50 return
      end

