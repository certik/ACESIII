
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR AND ATOMIC MASS VECTOR
C     FOR MASS WEIGHTING.  THE SORT IS FIRST DONE ON THE X COORDINATE,
C     BUT IF IDENTICAL X VALUES ARE ENCOUNTERED, THE Y COORDINATE
C     IS THEN SORTED AS WELL, AND ALSO THE Z IF NECESSARY.

      SUBROUTINE SORT2(XX,Y,NORD,IORDER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"

      DIMENSION XX(3*IORDER),Y(3*IORDER),NORD(IORDER*2)

      DIMENSION X(3*MXATMS),SCR(3*MXATMS)

      TOL=1.d-5

      CALL ZERO(X,3*MXATMS)

c   o sort the X coordinates and mass weight them
      NATOMS=IORDER
      CALL DCOPY(3*IORDER,XX,1,SCR,1)
      IOFF=1
      DO I=1,NATOMS
         NORD(I+NATOMS)=I
         IOFF=IOFF+3
      END DO
      CALL DCOPY(NATOMS,SCR,3,X,1)

      CALL PIKSR3(NATOMS,X,NORD(NATOMS+1))
      DO ITARGET=1,NATOMS
         ISOURCE=NORD(ITARGET+NATOMS)
         IOFFSRC=3*(ISOURCE-1)+1
         IOFFTAR=3*(ITARGET-1)+1
         CALL DCOPY(3,SCR(IOFFSRC),1,Y(IOFFTAR),1)
      END DO

c   o search for clusters of values in the X vector
      IFINDX=1
1     CONTINUE
         NLEFT=NATOMS-IFINDX+1
         ILOC=IFINDNE(NLEFT,X(IFINDX),1,X(IFINDX),TOL)+IFINDX-1
         ICLSIZX=ILOC-IFINDX
         IF (ILOC.NE.IFINDX+1) THEN
c         o located cluster --> sort the Y values
            CALL DCOPY(ICLSIZX,Y(2+3*(IFINDX-1)),3,X,1)
            CALL PIKSR3(ICLSIZX,X,NORD(NATOMS+IFINDX))
            DO ITARGET0=1,ICLSIZX
               ISOURCE=NORD(ITARGET0+NATOMS+IFINDX-1)
               ITARGET=ITARGET0+IFINDX-1
               IOFFSRC=3*(ISOURCE-1)+1
               IOFFTAR=3*(ITARGET-1)+1
               CALL DCOPY(3,SCR(IOFFSRC),1,Y(IOFFTAR),1)
            END DO
c         o search for clusters of values in the Y vector within this X cluster
            IFINDY=1
2           CONTINUE
               NLEFTY=ICLSIZX-IFINDY+1
               ILOC=IFINDNE(NLEFTY,X(IFINDY),1,X(IFINDY),TOL)+IFINDY-1
               ICLSIZY=ILOC-IFINDY
               IF (ICLSIZY.NE.1) THEN
c               o located cluster --> sort the Z values
                  IPOS=IFINDX+IFINDY-1
                  CALL DCOPY(ICLSIZY,Y(3+3*(IPOS-1)),3,X,1)
                  CALL PIKSR3(ICLSIZY,X,NORD(NATOMS+IPOS))
                  DO ITARGET0=1,ICLSIZY
                     ISOURCE=NORD(ITARGET0+NATOMS+IPOS-1)
                     ITARGET=ITARGET0+IPOS-1
                     IOFFSRC=3*(ISOURCE-1)+1
                     IOFFTAR=3*(ITARGET-1)+1
                     CALL DCOPY(3,SCR(IOFFSRC),1,Y(IOFFTAR),1)
                  END DO
               END IF
               IFINDY=IFINDY+ICLSIZY
            IF (IFINDY.LT.ICLSIZX) GOTO 2
         END IF
         IFINDX=IFINDX+ICLSIZX
      IF (IFINDX.LT.NATOMS) GOTO 1

c YAU - why is this here?
      IOFF=1
      DO I=1,NATOMS
         Z=1.0
         CALL DSCAL(3,Z,Y(IOFF),1)
         IOFF=IOFF+3
      END DO

      CALL ICOPY(NATOMS,NORD(NATOMS+1),1,NORD,1)

      RETURN
      END

