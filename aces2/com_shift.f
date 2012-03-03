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
      subroutine com_shift(q, natoms, atmass, iprnt)
c-----------------------------------------------------------------------------
c   Performs center-of-mass trnanslation of a molecular system.
c-----------------------------------------------------------------------------
      implicit none
      integer natoms, iprnt
      double precision q(3*natoms), atmass(natoms)  

      integer i, j, luout, idegen
      double precision cmx, cmy, cmz, molwt
      double precision cm(3)
C
C TRANSLATE TO CENTER OF MASS
C
      LUOUT = 6
      IDEGEN=0
      IF(IPRNT.GE.3)WRITE(LUOUT,7733)(ATMASS(J),J = 1,NATOMS)
 7733 FORMAT(1X,F15.10)
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      MOLWT=0.D0
      DO 20 I = 1,NATOMS
        CMX = ATMASS(I)*Q(3*I-2)+CMX
        CMY = ATMASS(I)*Q(3*I-1)+CMY
        CMZ = ATMASS(I)*Q(3*I)+CMZ
   20 MOLWT = MOLWT+ATMASS(I)
      IF (MOLWT.LT.1.0D-8) THEN
         WRITE(LUOUT,*) '@com_shift: No real atoms in Z-matrix.'
         CALL ERREX
      END IF
      CM(1) = CMX/MOLWT
      CM(2) = CMY/MOLWT
      CM(3) = CMZ/MOLWT
      DO I = 1,NATOMS
        DO J = 0,2
          Q(3*I-J) = Q(3*I-J)-CM(3-J)
        ENDDO
      ENDDO 
      IF(IPRNT .GE. 4) THEN
           WRITE(LUOUT,*)
     &     'After translation to center of mass coordinates '
           WRITE(LUOUT,80)(Q(I),I = 1,NATOMS)
      ENDIF
   80 FORMAT((4X,3(2X,F16.12)))

      write(6,*) ' @symmetry-i, Coordinates after  COM shift '
      do i=1,natoms
        write(6,'(3F20.12)') q(3*i-2),q(3*i-1),q(3*i)
      enddo

      return 
      end 
