      SUBROUTINE GRID(NPTSX,TVALX,BVALX,NPTSY,TVALY,
     $     BVALY,NPTSZ,TVALZ,BVALZ,PCOORD,
     $     PCRDX,PCRDY,PCRDZ,IBIGLOOP,JBIGLOOP)
c     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION PCOORD((NPTSZ+1)*3)
      DIMENSION PCRDX(NPTSZ+1),PCRDY(NPTSZ+1),
     $     PCRDZ(NPTSZ+1)
C     
C     determine the points on the cube
      DO 10 I=0,NPTSX
         PCRDX(I+1)=(BVALX+((TVALX-BVALX)/REAL(NPTSX))*I)
 10   CONTINUE
      DO 11 I=0,NPTSY
         PCRDY(I+1)=(BVALY+((TVALY-BVALY)/REAL(NPTSY))*I)
 11   CONTINUE
      DO 12 I=0,NPTSZ
         PCRDZ(I+1)=(BVALZ+((TVALZ-BVALZ)/REAL(NPTSZ))*I)
 12   CONTINUE
C     
C     fill in the pcoord vector for this given row of points
      K=1
      DO 20 J=1,(NPTSZ+1)
         PCOORD(K)=PCRDX(IBIGLOOP)
C     WRITE(*,*) "PCOORD(",K,")=",PCOORD(IBIGLOOP)
         K=K+1
         PCOORD(K)=PCRDY(JBIGLOOP)
C     WRITE(*,*) "PCOORD( n ~ K,")=",PCOORD(JBIGLOOP)
         K=K+1
         PCOORD(K)=PCRDZ(J)
C     WRITE(*,*) "PCOORD(",K, n ) = t. ~ PCOORD(K)
         K=K+1
 20   CONTINUE
C     
      return
      END SUBROUTINE GRID

