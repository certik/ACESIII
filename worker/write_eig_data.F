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
      subroutine write_eig_data(ca, cb, focka, fockb, n, need_predef)
c---------------------------------------------------------------------------
c   Writes eigenvector/eigenvalue data to jobarc file.
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer n
      logical need_predef(*)
      double precision ca(n,n), cb(n,n), focka(n,n), fockb(n,n)
      double precision scr(n)

      integer i, j, ierd

      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      equivalence (ibuf(1), dbuf(1))
      common dbuf
#endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif
      if (intpkg .eq. flocke_package) then

c---------------------------------------------------------------------------
c   Transform the CA and CB data from ERD order/scaling to VMOL order/scaling.
c---------------------------------------------------------------------------

         if (need_predef(1) .or.
     *       need_predef(2)) then
            do j = 1, n
               do i = 1, n
                  scr(i) = ca(i,j)
               enddo

               do i = 1, n
                  ierd = ibuf(ierdind+i-1)
                  ca(ierd,j) = scr(i)/dbuf(iscale_fac+ierd-1)
               enddo
            enddo   
         endif

         if (need_predef(3)) then 
            do j = 1, n
               do i = 1, n
                  scr(i) = cb(i,j)
               enddo

               do i = 1, n
                  ierd = ibuf(ierdind+i-1)
                  cb(ierd,j) = scr(i)/dbuf(iscale_fac+ierd-1)
               enddo
            enddo
         endif   

c---------------------------------------------------------------------------
c   Save the diagonals of the FOCKA array.
c---------------------------------------------------------------------------

         if (need_predef(9)) then
            do i = 1, n
               scr(i) = focka(i,i)
            enddo

c---------------------------------------------------------------------------
c   Unscramble the shells from the binpack'ed order.
c---------------------------------------------------------------------------

            call descramble_coefficients(ca, n, scr, nfps, ixshells, 
     *                                nshells)

            do i = 1, n
               focka(i,i) = scr(i)
            enddo
         endif   

c---------------------------------------------------------------------------
c   Save the diagonals of the FOCKB array.
c---------------------------------------------------------------------------

         if (need_predef(10)) then
            do i = 1, n
               scr(i) = fockb(i,i)
            enddo

c---------------------------------------------------------------------------
c   Unscramble the shells from the binpack'ed order.
c---------------------------------------------------------------------------

            call descramble_coefficients(cb, n, scr, nfps, ixshells,
     *                                nshells)

            do i = 1, n
               fockb(i,i) = scr(i)
            enddo
         endif  ! need_predef
      endif     ! flocke_package

c---------------------------------------------------------------------------
c   Write the eigenvector data to JOBARC.
c---------------------------------------------------------------------------

         if (need_predef(1) .or.
     *       need_predef(2)) then
            call dputrec(1,'JOBARC','SCFEVCA0',n*n, ca)
         endif

         if (need_predef(3)) then
            call dputrec(1,'JOBARC','SCFEVCB0',n*n, cb)
         endif

c---------------------------------------------------------------------------
c   Write the eigenvalue data to JOBARC.
c---------------------------------------------------------------------------

      if (need_predef(9)) then
         do i = 1, n
            scr(i) = focka(i,i)
         enddo

         call dputrec(1, 'JOBARC', 'SCFEVLA0',n, scr)
      endif

c---------------------------------------------------------------------------
c   Write the data to JOBARC.
c---------------------------------------------------------------------------

      if (need_predef(10)) then
         do i = 1, n
            scr(i) = fockb(i,i)
         enddo

         call dputrec(1, 'JOBARC', 'SCFEVLB0',n, scr)
      endif

c---------------------------------------------------------------------------
c   Write out OCCUPYA0, OCCUPYB0.
c---------------------------------------------------------------------------

      call iputrec(1, 'JOBARC', 'OCCUPYA0', 1, nalpha_occupied)
      call iputrec(1, 'JOBARC', 'OCCUPYB0', 1, nbeta_occupied)

c---------------------------------------------------------------------------
c   Write out the SCFOCCA and SCFOCCB values.
c---------------------------------------------------------------------------

      do i = 1, n
         if (i .le. nalpha_occupied) then
            scr(i) = 1.
         else
            scr(i) = 0.
         endif
      enddo

      call dputrec(1, 'JOBARC', 'SCFOCCA', nalpha_occupied, scr)

      do i = 1, n
         if (i .le. nbeta_occupied) then
            scr(i) = 1.
         else
            scr(i) = 0.
         endif
      enddo

      call dputrec(1, 'JOBARC', 'SCFOCCB', nbeta_occupied, scr)


      return
      end

