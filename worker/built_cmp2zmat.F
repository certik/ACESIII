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
C
      Subroutine built_cmp2zmat(Natoms, Iordgp) 
C
      implicit none
      include 'machine_types.h'
C
      Character*4 CMPPGP      
      integer Natoms, Iordgp
      integer Norbits, Ncen
      integer i
      integer itype, L
      Integer NCenW1(Natoms), NCenW2(Natoms), CenMap(Natoms)
      Double precision Coord1(3, Natoms), Coord2(3, IOrdGp)

      Call GetCRec(10, 'JOBARC', 'COMPPTGP', 4, CMPPGP) 
      Call igetrec(10, 'JOBARC', 'COMPNORB', 1, NOrbits)
      if (Norbits .gt. Natoms) then
         print *,'Error: COMPPOPV has ',norbits,' entries in GETREC,',
     *     ' should be less than natoms = ',natoms
         call abort_job()
      endif
 
      Call igetrec(10 , 'JOBARC', 'COMPPOPV', NOrbits, NCenW1) 
      Call igetrec(0 , 'JOBARC', 'COMPMEMB', NCen, NCenW2) 
      if (Ncen .gt. Natoms) then
         print *,'Error: COMPMEMB has ',ncen,' entries in GETREC,',
     *     ' should be less than natoms = ',natoms
         call abort_job()
      endif

      Call igetrec(10 , 'JOBARC', 'COMPMEMB', NCen, NCenW2)
      L = 3 * NCen 
      Call dgetrec( -1, 'JOBARC', 'COORD   ', L, Coord1)
      if (CMPPGP(1:3) .eq. 'D2 ')
     *   CALL igetrec(20,'JOBARC','D2TYPXYZ',1,ITYPE)
C
      Call Remapx(CMPPGP, IORDGP, NCen, NOrbits, Coord1, Coord2,
     &           NCenW2, NCenW1, CenMap, itype)
c
      Call iputrec( 20, 'JOBARC', 'MAP2ZMAT', NCen, CenMap) 
      
      Return
      ENd
