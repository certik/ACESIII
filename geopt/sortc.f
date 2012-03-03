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

      SUBROUTINE SortC (NAtms, XX, AtNr, Y, NORD, X)
C     X is a scratch array of the same type and dimension as XX and Y
C     AtNr are the atomix numbers of the centers (dummy atom = 0)
C
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR AND ATOMIC NUMBER VECTOR (AtNr)
C
      Implicit double precision (a-h,o-z)
      DIMENSION X(3*NATMS),XX(3*NATMS),Y(3*NATMS),NORD(NATMS)
      Integer AtNr(natms)
C
      ILINE(J)=1+J/3
C
C     SORT ON THE X - IF TWO X'S ARE EQUIVALENT, SORT ON Y AND SO ON.
C     
      DO 80 I=1,3*NATMS
 80      X(I)=XX(I)
C
C     FIRST GIVE DUMMY ATOMS RIDICULOUS SCRATCH COORDINATES - ENSURES
C     THAT THEY WILL WIND UP AT THE BOTTOM OF THE LIST
C
      DO 81 I=1,3*NATMS-2,3
         IF(AtNr(ILINE(I)) .eq. 0)THEN
            DO 82 J=0,2
 82            X(J+I) = -99995.
         ENDIF
 81   CONTINUE
      JK=1
 429  J=1
      DO 96 I=1,3*NATMS-2,3
C
C     CONTINUE WITH SORTING.
C
         IF(X(I)-X(J).GT.1D-6)J=I
         IF(DABS(X(I)-X(J)).LT.1D-6)THEN
            IF(X(I+1)-X(J+1).GT.1D-6)J=I
            IF(DABS(X(I+1)-X(J+1)).LT.1D-6)THEN
               IF(X(I+2)-X(J+2).GT.1D-6)J=I
            ENDIF
         ENDIF
 96   CONTINUE
      DO 93 I=0,2
C
C     Mass-WEIGHT SORTED VECTOR - WILL ZERO ELEMENTS CORRESPONDING
C     TO DUMMY ATOMS SO THEY DONT MUCK UP THE SYMMETRY CHECKING.
C     
         Y(3*JK-2+I)=X(J+I)*AtNr(ILINE(J))
 93      X(J)=-99999.D0
      NORD(JK)=(J+2)/3
      JK=JK+1
      if(jk.eq.NATMS+1)go to 999
      go to 429
 999  Continue
      Return
      end
