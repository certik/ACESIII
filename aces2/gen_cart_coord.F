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

      SUBROUTINE GEN_CART_COORD()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      double precision scr(natoms)

#include "mxatms.par"
      LOGICAL LINEAR

#include "coord_dup.com"
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol

      CALL ZERO(Q,3*MXATMS)
      IF (NATOMS.NE.1) Q(1,2)=R(1,2)
C
C*** FIND THE FIRST ATOM (IF ANY) WHICH IS NOT COLLINEAR WITH ATOMS
C 1 AND 2.  ASSUME THAT THE ATOMS ARE COLLINEAR IF THE COSINE OF THE
C ANGLE IS GREATER THAN 0.99999
C
      IF (NATOMS.GT.2) THEN
         LINEAR=.TRUE.
         IND=2
         DO WHILE (LINEAR.AND.IND.LT.NATOMS)
            CCOS=DCOS(R(2,IND+1))
            LINEAR = (DABS(CCOS).GT.0.999999)
            IND=IND+1
         END DO
         IF (LINEAR) THEN
            WRITE(*,'(T3,A)')'@GMETRY: Linear frameworks not permitted.'
            WRITE(*,'(T3,A)')'         Use dummy atoms in the Z-matrix.'
            CALL ERREX
         END IF
c      o put the first IND-1 atoms on the X axis
         Q(1,IND-1)=R(1,IND-1)
         Q(2,IND-1)=0.0D0
         Q(3,IND-1)=0.0D0
         DO K=IND-2,1,-1
            Q(1,K)=Q(1,K+1)-R(1,K+1)
            Q(2,K)=0.0D0
            Q(3,K)=0.0D0
         END DO
c      o hardcode the IND atom
         Q(1,IND)=Q(1,IND-1)-R(1,IND)*CCOS
         Q(2,IND)=R(1,IND)*DSIN(R(2,IND))
         Q(3,IND)=0.0D0
c      o process the remaining atoms
         DO I=IND+1,NATOMS
            COSA=DCOS(R(2,I))
            MB=NCON(2,I)
            MC=NCON(1,I)
            XB=Q(1,MB)-Q(1,MC)
            YB=Q(2,MB)-Q(2,MC)
            ZB=Q(3,MB)-Q(3,MC)
            RBC=1.0D0/DSQRT(XB*XB+YB*YB+ZB*ZB)
            IF (ABS(COSA).GE.0.999999) THEN
C              ATOMS MC, MB, AND (I) ARE COLLINEAR
               RBC=R(1,I)*RBC
               Q(1,I)=Q(1,MC)-XB*RBC
               Q(2,I)=Q(2,MC)-YB*RBC
               Q(3,I)=Q(3,MC)-ZB*RBC
            ELSE
C              THE ATOMS ARE NOT COLLINEAR
               MA=NCON(3,I)
               XA=Q(1,MA)-Q(1,MC)
               YA=Q(2,MA)-Q(2,MC)
               ZA=Q(3,MA)-Q(3,MC)
C
C ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
               XYB=DSQRT(XB*XB+YB*YB)
               K=-1
               IF (XYB.LE.0.1) THEN
                  dTmp =  ZA
                  ZA   = -XA
                  XA   = dTmp
                  dTmp =  ZB
                  ZB   = -XB
                  XB   = dTmp
                  XYB  = DSQRT(XB*XB+YB*YB)
                  K    = 1
               END IF
C
C ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
               COSTH=XB/XYB
               SINTH=YB/XYB
               XPA=XA*COSTH+YA*SINTH
               YPA=YA*COSTH-XA*SINTH
               SINPH=ZB*RBC
               COSPH=DSQRT(DABS(1.0-SINPH*SINPH))
               XQA=XPA*COSPH+ZA*SINPH
               ZQA=ZA*COSPH-XPA*SINPH
C
C ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
               YZA=DSQRT(YPA*YPA+ZQA*ZQA)
               COSKH=YPA/YZA
               SINKH=ZQA/YZA
C
C COORDINATES :-   R=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C NONE ARE NEGATIVE.
C THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
               SINA=DSIN(R(2,I))
               SIND=DSIN(R(3,I))
               COSD=DCOS(R(3,I))
               XD=R(1,I)*COSA
               YD=R(1,I)*SINA*COSD
               ZD=R(1,I)*SINA*SIND
C
C TRANSFORM THE COORDINATES BACK TO THE ORIGIONAL SYSTEM.
C
               YPD=YD*COSKH-ZD*SINKH
               ZPD=ZD*COSKH+YD*SINKH
               XPD=XD*COSPH-ZPD*SINPH
               ZQD=ZPD*COSPH+XD*SINPH
               XQD=XPD*COSTH-YPD*SINTH
               YQD=YPD*COSTH+XPD*SINTH
               IF (K.GE.1) THEN
                  dTmp = -ZQD
                  ZQD  =  XQD
                  XQD  = dTmp
               END IF
               Q(1,I)=XQD+Q(1,MC)
               Q(2,I)=YQD+Q(2,MC)
               Q(3,I)=ZQD+Q(3,MC)
            END IF
C        END DO I=IND+1,NATOMS
         END DO
C     END IF (NATOMS.GT.2)
      END IF

      IF (IPRNT.GE.200) THEN
         WRITE(*,444)
 444     FORMAT(T3,'@GMETRY-I, Cartesian coordinates before scaling:')
         WRITE(*,445) (I,Q(1,I),Q(2,I),Q(3,I),I=1,NATOMS)
 445     FORMAT(T5,I3,1X,F12.6,1X,F12.6,1X,F12.6)
      END IF
c YAU : Why is this here?
      NDIS=0
      IF (NDIS.NE.0) THEN
         dTmp=1.d0/0.529177249d0
c         CALL DSCAL(3*NATOMS,dTmp,Q,1)
         do i = 1, natoms
            q(1,i) = q(1,i) * dtmp
            q(2,i) = q(2,i) * dtmp
            q(3,i) = q(3,i) * dtmp
         enddo

      END IF
      IF (IPRNT.GE.200) THEN
         WRITE(*,446)
 446     FORMAT(T3,'@GMETRY-I, Cartesian coordinates after scaling:')
         WRITE(*,445) (I,Q(1,I),Q(2,I),Q(3,I),I=1,NATOMS)
      END IF
C
C Rotate the frame (supposedly to match GAUSSIAN).
C   X' =  Y
C   Y' = -Z
C   Z' =  X
C
#ifdef _UNICOS
c   o vector shuffle
      call dcopy(natoms,q(3,1),3,scr,1)
      call dcopy(natoms,q(1,1),3,q(3,1),3)
      call dcopy(natoms,q(2,1),3,q(1,1),3)
      do i=1,natoms
         q(2,i) = -scr(i)
      end do
#else
c   o scalar shuffle
      DO I=1,NATOMS
         dTmp   =  Q(3,I)
         Q(3,I) =  Q(1,I)
         Q(1,I) =  Q(2,I)
         Q(2,I) = -dTmp
      END DO
#endif

      IF (IPRNT.GE.20) THEN
         WRITE(*,555)
 555     FORMAT(T3,'@GMETRY-I, Cartesian coordinates from Z-matrix:')
         WRITE(*,445) (I,Q(1,I),Q(2,I),Q(3,I),I=1,NATOMS)
      END IF

      RETURN
      END

