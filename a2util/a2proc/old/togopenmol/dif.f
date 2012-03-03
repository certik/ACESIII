      SUBROUTINE DIF
C     
C     This subroutine determines the density
C     difference in the two input cubes. Note
C     that the proram assumes that both the
C coordinates of the atoms and the points
C     are identical.
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      CHARACTER*25 LINE1
      CHARACTER*25 LINE2
C     
      DIMENSION P(6),Q(6),O(6)
C     
C     I apologize for the following series of commands, but I
C couldn't get things to work otherwise.
      I=ISHELL('cp den.cube denone.cube')
      I=ISHELL('cp den2.cube dentwo.cube')
      I=ISHELL('rm den.cube den2.cube')
C     
      OPEN(UNIT=13,FILE='dentwo.cube')
      OPEN(UNIT=12,FILE='denone.cube')
      OPEN(UNIT=14,FILE='dif.cube')
C     
      READ(12,"(A25)") LINE1
      READ(12,"(A25)") LINE2
      WRITE(14,*) LINE1
      WRITE(14,*) "The Delta Density"
C
      READ(12,*) NATOMS,P1,P2,P3
      WRITE(14,"(I5,3 F12.6)") NATOMS,P1,P2,P3
      READ(12,*) IXPTS,P1,P2,P3
      WRITE(14,"(I5,3 F12.6)") IXPTS,P1,P2,P3
      READ(12,*) IYPTS,P1,P2,P3
      WRITE(14,"(I5,3 F12.6)") IYPTS,P1,P2,P3
      READ(12,*) IZPTS,P1,P2,P3
      WRITE(14,"(I5,3 F12.6)") IZPTS,P1,P2,P3
C     
      DO 10 I=1,ABS(NATOMS)
         READ(12,*) INT,P1,P2,P3,P4
         WRITE(14,"(I5,4 F12.6)") INT,P1,P2,P3,P4
 10   CONTINUE
C
      MODVAL=mod(IZPTS,6)
      NGRPS=IZPTS/6
C
      DO 25 K=1,(6+ABS(NATOMS))
         READ(13,*)
 25   CONTINUE
C
      DO 20 I=1,(IXPTS*IYPTS)
C
         DO 30 J=1,NGRPS
            READ(12,*) (P(K),K=1,6)
            READ(13,*) (Q(K),K=1,6)
            DO 40 K=1,6
               O(K)=P(K)-Q(K)
 40         CONTINUE
            WRITE(14,"(6(D13.5))") (O(K),K=1,6)
 30      CONTINUE
C
         IF (MODVAL.NE.0) THEN
            READ(12,*) (P(K),K=1,MODVAL)
            READ(13,*) (Q(K),K=1,MODVAL)
            DO 50 K=1,MODVAL
               O(K)=P(K)-Q(K)
 50         CONTINUE
            WRITE(14,"(6(D13.5))") (O(K),K=1,MODVAL)
         END IF
C
 20   CONTINUE
C
      I=ISHELL('cp denone.cube den.cube')
      I=ISHELL('cp dentwo.cube den2.cube')
      I=ISHELL('rm denone.cube dentwo.cube')
C     
      END SUBROUTINE DIF

