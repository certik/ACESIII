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
      SUBROUTINE GETLAMBDA(HESMOD, GRDHES, RFAMAT, NOPT, EVRFAMAT, 
     &                     LMBDAN, LMBDAP, IADITNL, IMODE, LUOUT, 
     &                     TS, RFA, QSTLST_CLIMB)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL TS, RFA, QSTLST_CLIMB
      DOUBLE PRECISION LMBDAN, LMBDAP, LP1, LP2 
C
      DIMENSION HESMOD(NOPT, NOPT), GRDHES(NOPT), EVRFAMAT(NOPT+IADITNL,
     &          NOPT+IADITNL), RFAMAT(NOPT+IADITNL, NOPT+IADITNL)
C
C If this is a simple LST or QST climb, we do not need Lambda(n) parameter.
      IF (.NOT. QSTLST_CLIMB) THEN

      JCT=1
      NDIM = NOPT + IADITNL

      do i = 1, NOPT+IADITNL
      do j = 1, NOPT+IADITNL
         rfamat(i,j)   = 0.d0
         evrfamat(i,j) = 0.d0
      enddo
      enddo

      DO 503 I = 1, NOPT
C
         IF(I .EQ. IMODE) GOTO 503
         RFAMAT(JCT, JCT)  = HESMOD(I,I)
         RFAMAT(JCT, NDIM) = GRDHES(I)
         RFAMAT(NDIM, JCT) = GRDHES(I)
         JCT = JCT + 1
C
 503  CONTINUE
C
C Diagonalize the RFA matrix
C
      CALL EIG(RFAMAT, EVRFAMAT, NDIM, NDIM, 0)
      LMBDAN = RFAMAT(1,1)
C
      WRITE(LUOUT, 2105)LMBDAN
 2105 FORMAT(T3,' Shift Parameter Lambda(n)=',f10.7,'.')
C
C     ENDIF (.NOT. QSTLST_CLIMB)
      ENDIF
C
C Calculate Lambda(P)
C
      IF (TS.OR.QSTLST_CLIMB) THEN
C
         IF (.NOT.QSTLST_CLIMB) THEN
         WRITE(LUOUT, 1651) IMODE
 1651    FORMAT(T3,' Eigenvector following is turned on.',/,t3,
     &          ' Mode ',i3,' will be searched uphill.')
         ELSE
         WRITE(LUOUT, 1652) IMODE
 1652    FORMAT(T3,' QST or LST climbing is turned on.',/,t3,
     &          ' Mode ',i3,' will be searched uphill.')
         ENDIF
C
         LP1 = HESMOD(IMODE,IMODE)*0.5D0
         LP2 = HESMOD(IMODE,IMODE)**2 + 4.D0*GRDHES(IMODE)**2
         LMBDAP = LP1 + 0.5D0*DSQRT(LP2)
C
         WRITE(LUOUT,2106)LMBDAP
 2106    FORMAT(T3,' Shift Parameter Lambda(p)=',f10.7)
C
      ELSE
         LMBDAP = 0.0D00
      ENDIF
C
      RETURN
      END


       
