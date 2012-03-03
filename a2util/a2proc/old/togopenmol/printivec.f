	SUBROUTINE PRINTIVEC(IVECNAM,LEN)
C       
C       Designed to print out an integer vector
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C       
	DIMENSION IVECNAM(LEN)
C       
	WRITE(*,*) "len=", LEN
	DO 20 I=1,LEN
	   WRITE(*,*) IVECNAM(I)
   20 CONTINUE
C       
	END SUBROUTINE PRINTIVEC
	
