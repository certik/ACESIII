      SUBROUTINE NORBMAKE(DENRELA,DENRELB,TEMPNOCOEFF,ANOCOEFF,
     $     BNOCOEFF,MOREORD,NBASP,XNOOCC,WORK)
C
c Diagonalize the alpha and beta correlated densities to obtain
C the NOs in the MO basis
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DENRELA(NBASP*NBASP),DENRELB(NBASP*NBASP),
     $     TEMPNOCOEFF(NBASP*NBASP),ANOCOEFF(NBASP*NBASP),
     $     MOREORD(NBASP),XNOOCC(NBASP*NBASP),WORK(3*NBASP-1),
     $     BNOCOEFF(NBASP*NBASP)
C     
      WRITE(*,*) "DENRELA:"
      CALL PRINTMAT(DENRELA,NBASP,NBASP)
      WRITE(*,*) "DENRELB:"
      CALL PRINTMAT(DENRELB,NBASP,NBASP)
      WRITE(*,*)
C
C ALPHA PART
C
      DO 10 I=1,(NBASP*NBASP)
         XNOOCC(I)=DENRELA(I)
 10   CONTINUE
C
C Diagonalize the correlated density matrix
      CALL EIG(XNOOCC,TEMPNOCOEFF,WORK,(NBASP),2)
C     
      WRITE(*,*) "Eigenvalues of alpha density matrix:"
      SUM=0.0
      DO 100 I=1,NBASP
         SUM=SUM+XNOOCC((I-1)*NBASP+I)
         WRITE(*,*) XNOOCC((I-1)*NBASP+I)
 100  CONTINUE
      WRITE(*,*)
      WRITE(*,*) "SUM=",SUM
      WRITE(*,*)
      WRITE(*,*) "Alpha NOs in MO basis in correlated order"
      CALL PRINTMAT(TEMPNOCOEFF,NBASP,NBASP)
C     Reorder the MOs (rows) to get ANOCOEFF
C We are going from Correlated to SCF ordering.
      DO 11 I=1,NBASP
         DO 12 J=0,(NBASP-1)
            ANOCOEFF(MOREORD(I)+J*NBASP)=TEMPNOCOEFF(I+J*NBASP)
 12      CONTINUE
 11   CONTINUE
      WRITE(*,*)
      WRITE(*,*) "Alpha NOs ln MO basis in scf order"
      CALL PRINTMAT(ANOCOEFF,NBASP,NBASP)
C
C     BETA PART
C     
      DO 20 I=1,(NBASP*NBASP)
         XNOOCC(I)=DENRELB(I)
 20   END DO
C     
C     Diagonalize the correlated density matrix
      CALL EIG(XNOOCC,TEMPNOCOEFF,WORK,(NBASP),2)
C
      WRITE(*,*) "Eigenvalues of beta density matrix:"
      SUM=0.0
      DO 200 I=1,NBASP
         SUM=SUM+XNOOCC((I-1)*NBASP+I)
         WRITE(*,*) XNOOCC((I-1)*NBASP+I)
 200  CONTINUE
      WRITE(*,*)
      WRITE(*,*) "SUM=",SUM
      WRITE(*,*)
         WRITE(*,*) "Beta NOs in MO basis in correlated order"
         CALL PRINTMAT(TEMPNOCOEFF,NBASP,NBASP)
C     Reorder the MOs (rows) to get BNOCOEFF
         DO 21 I=1,NBASP
            DO 22 J=0,(NBASP-1)
               BNOCOEFF(MOREORD(I)+J*NBASP)=TEMPNOCOEFF(I+J*NBASP)
 22         CONTINUE
 21      CONTINUE
         WRITE(*,*)
         WRITE(*,*) "Beta NOs in MO basis in scf order"
         CALL PRINTMAT(BNOCOEFF,NBASP,NBASP)
C     
C
         RETURN
         END SUBROUTINE NORBMAKE
         
