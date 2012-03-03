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
      SUBROUTINE EXPND4(VPACK,VFULL,LENGTH,IDIM)
C
C THIS ROUTINE EXPANDS A "PACKED" MATRIX TO FULL FORMAT, WHERE
C   THE PACKED ORDER IS :  V(1,1),V(2,1),...,V(N,1),V(2,2),...
C
C THIS IS USEFUL FOR EXPANDING THE D FUNCTION REPRESENTATIONS TO
C   A FULL 3x3 MATRIX
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DLOC(6),FLOC(10),GLOC(15),DPAK,FPAK,GPAK,DSIZ,FSIZ,GSIZ
      DIMENSION VPACK(1),VFULL(LENGTH*LENGTH)
      DIMENSION DFAC(6),FFAC(10),GFAC(15)
C
      PARAMETER (HALF = 0.5D0)
      PARAMETER (ONE  = 1.0D0)
      PARAMETER (TWO  = 2.0D0)
      PARAMETER (THIRD=1.0D0/3.0D0)
      PARAMETER (SIXTH=1.0D0/6.0D0)
      PARAMETER (FOURTH=0.25D0)
      PARAMETER (TWELFTH=1.0D0/12.0D0)
C   
      DATA DPAK /6/
      DATA FPAK /10/
      DATA GPAK /15/
      DATA DSIZ /9/
      DATA FSIZ /27/
      DATA GSIZ /81/
      DATA DLOC /1,2,3,5,6,9/
      DATA FLOC /1,2,3,5,6,9,14,15,18,27/
      DATA GLOC /1,2,3,5,6,9,14,15,18,27,41,42,45,54,81/
      DATA DFAC /ONE,HALF,HALF,ONE,HALF,ONE/
      DATA FFAC /ONE,THIRD,THIRD,THIRD,SIXTH,THIRD,ONE,THIRD,
     &           THIRD,ONE/
      DATA GFAC /ONE,FOURTH,FOURTH,SIXTH,TWELFTH,
     &           SIXTH,FOURTH,TWELFTH,TWELFTH,FOURTH,
     &           ONE,FOURTH,SIXTH,FOURTH,ONE/
C
      INDX2(I,J)=I+3*(J-1)
      INDX3(I,J,K)=I+3*(J-1)+9*(K-1)
      INDX4(I,J,K,L)=I+3*(J-1)+9*(K-1)+27*(L-1) 
C
C CODE FOR TWO DIMENSIONAL MATRICES (SUCH AS d FUNCTIONS)
C
      IF(IDIM.EQ.2)THEN
       CALL ZERO  (VFULL,DSIZ)
       CALL VECPRD(DFAC,VPACK,VPACK,DPAK)
       CALL SCATTER(DPAK,VFULL,DLOC,VPACK)
       DO 10 I=1,LENGTH
        DO 20 J=I,LENGTH
         X=VFULL(INDX2(J,I))
         VFULL(INDX2(I,J))=X
20      CONTINUE
10     CONTINUE
C
C CODE FOR THREE DIMENSIONAL MATRICES (SUCH AS f FUNCTIONS)
C
      ELSEIF(IDIM.EQ.3)THEN
       CALL ZERO  (VFULL,FSIZ)
       CALL VECPRD (FFAC,VPACK,VPACK,FPAK)
       CALL SCATTER(FPAK,VFULL,FLOC,VPACK)
       DO 110 I=1,LENGTH
        DO 120 J=I,LENGTH
         DO 130 K=J,LENGTH
          X=VFULL(INDX3(K,J,I))
          VFULL(INDX3(J,I,K))=X
          VFULL(INDX3(J,K,I))=X
          VFULL(INDX3(K,I,J))=X
          VFULL(INDX3(I,J,K))=X
          VFULL(INDX3(I,K,J))=X
130      CONTINUE
120     CONTINUE
110    CONTINUE
C
C CODE FOR FOUR DIMENSIONAL MATRICES (g FUNCTIONS)
C
      ELSEIF(IDIM.EQ.4)THEN
       CALL ZERO  (VFULL,GSIZ)
       CALL VECPRD (GFAC,VPACK,VPACK,GPAK)
       CALL SCATTER(GPAK,VFULL,GLOC,VPACK)
       DO 210 I=1,LENGTH
        DO 220 J=I,LENGTH
         DO 230 K=J,LENGTH
CDIR$ NOVECTOR
*VOCL LOOP,SCALAR
          DO 240 L=K,LENGTH
           X=VFULL(INDX4(L,K,J,I))
C
           VFULL(INDX4(I,J,K,L))=X
           VFULL(INDX4(I,J,L,K))=X
           VFULL(INDX4(I,L,J,K))=X
           VFULL(INDX4(I,L,K,J))=X
           VFULL(INDX4(I,K,J,L))=X
           VFULL(INDX4(I,K,L,J))=X
C
           VFULL(INDX4(J,I,K,L))=X
           VFULL(INDX4(J,I,L,K))=X
           VFULL(INDX4(J,L,I,K))=X
           VFULL(INDX4(J,L,K,I))=X
           VFULL(INDX4(J,K,I,L))=X
           VFULL(INDX4(J,K,L,I))=X
C
           VFULL(INDX4(K,I,J,L))=X
           VFULL(INDX4(K,I,L,J))=X
           VFULL(INDX4(K,J,I,L))=X
           VFULL(INDX4(K,J,L,I))=X
           VFULL(INDX4(K,L,I,J))=X
           VFULL(INDX4(K,L,J,I))=X
C
           VFULL(INDX4(L,I,J,K))=X
           VFULL(INDX4(L,I,K,J))=X
           VFULL(INDX4(L,K,I,J))=X
           VFULL(INDX4(L,J,I,K))=X
           VFULL(INDX4(L,J,K,I))=X
C        
240       CONTINUE
230      CONTINUE
220     CONTINUE
210    CONTINUE
C
      ENDIF
C
      RETURN
      END
