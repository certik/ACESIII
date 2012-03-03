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
      SUBROUTINE ASSIGN_XYZ_LABELS (TOTNO_XYZ, XYZ_LABEL)
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Integer TOTNO_XYZ
C
C Assign a label for xyz coordiantes of each atom
C
#include "mxatms.par"
C
      CHARACTER*4 XYZ_LABEL(TOTNO_XYZ)
      CHARACTER*3 LABEL(3*MXATMS)
C
      do i = 1, min(9,3*MXATMS)
         write(label(i),'(2i1)') 0, i
      end do
      do i = 10, min(99,MXATMS)
         write(label(i),'(i2)') i
      end do
      do i = 100, 3*MXATMS
         write(label(i),'(i3)') i
      end do
C
      DO I=1,TOTNO_XYZ
         XYZ_LABEL(I)='R'//LABEL(I)
      END DO
C            
      RETURN
      END
   
