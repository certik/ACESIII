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
      SUBROUTINE MORSEA(IATOM,JATOM,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(80,80)
      CALL ZERO(VEC,6400)
      VEC(1,1)=0.993509D0
      VEC(17,13)=0.702441D0
      VEC(13,9)=1.224876D0
      VEC(13,1)=0.685881D0
      VEC(13,8)=1.089910D0
      VEC(5,5)=0.927159D0
      VEC(17,5)=0.843132D0
      VEC(17,4)=0.779540D0
      VEC(9,4)=0.962386D0
      VEC(4,1)=0.731181D0
      VEC(8,4)=1.064388D0
      VEC(9,5)=1.143940D0
      VEC(5,1)=0.866977D0    
      VEC(7,5)=1.122298D0
      VEC(8,5)=1.129189D0
      VEC(6,6)=1.078581D0
      VEC(17,6)=0.385271D0
      VEC(6,1)=0.996595D0
      VEC(17,17)=1.077926D0
      VEC(17,9)=1.422952D0
      VEC(7,6)=1.187941D0
      VEC(8,6)=1.238251D0
      VEC(16,6)=0.970701D0
      VEC(17,1)=0.920386D0
      VEC(9,1)=1.196723D0
      VEC(3,3)=0.388674D0
      VEC(3,1)=0.582811D0
      VEC(7,7)=1.296765D0
      VEC(11,11)=0.372353D0
      VEC(11,1)=0.562428D0
      VEC(8,7)=1.316422D0
      VEC(8,8)=1.266555D0
      VEC(8,1)=1.142097D0
      VEC(16,16)=0.870337D0
      VEC(17,14)=0.753723D0
      VEC(14,9)=0.939980D0
      VEC(14,7)=1.008525D0
      VEC(14,8)=1.011156D0
      VEC(16,8)=1.040830D0
      I1=MAX(IATOM,JATOM)
      J1=MIN(IATOM,JATOM)
      A=VEC(I1,J1)
      RETURN
      END
