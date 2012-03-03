
C THIS ROUTINE CONSTRUCTS THE FORCE CONSTANT MATRIX FOR SYMMETRY
C BLOCK NDIM FROM NUMERICAL DIFFERENTIATION OF THE ENERGY.

c INPUT
c integer NDIM      :
c double  ENERGY(*) :
c integer INVOP     :
c double  STPSIZ    :
c double  E0        :
c integer NDSCR     :

c OUTPUT
c double FCM(NDIM,NDIM) :
c double DSCR(NDSCR)    :

      SUBROUTINE ENER2FCM(NDIM,ENERGY,FCM,INVOP,STPSIZ,E0,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ENERGY(*),FCM(NDIM,NDIM)
      double precision dscr(ndim,*)

      LOGICAL ONEGRD,PRINTQ

      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA FACT /5.14048D03/

      PRINTQ=(IFLAGS(1).GE.1)

C IF THE REPRESENTATION IS TOTALLY SYMMETRIC, THEN WE HAVE DONE ONLY
C ONE DISPLACEMENT FOR EACH DEGREE OF FREEDOM.

      ONEGRD=(INVOP.GT.0)

      IF (ONEGRD) THEN

c      o evaluate diagonal force constant matrix elements
         DTMP=2.d0/(STPSIZ*STPSIZ)
         E00=-DTMP*E0
         IOFFT=1
         DO IDIM=1,NDIM
            FCM(IDIM,IDIM)=E00+DTMP*ENERGY(IOFFT)
            IOFFT=IOFFT+IDIM+1
         END DO

c      o compute off-diagonal upper triangle force constants
         DTMP=0.5d0*DTMP
         E00=-DTMP*E0
         IOFFT0=2
         DO IDIM1=2,NDIM
            Z=E00-0.5d0*FCM(IDIM1,IDIM1)
            IOFFT=IOFFT0
            DO IDIM2=1,IDIM1-1
               FCM(IDIM2,IDIM1)=Z+DTMP*ENERGY(IOFFT)
     &                           -0.5d0*FCM(IDIM2,IDIM2)
               IOFFT=IOFFT+1
            END DO
            IOFFT0=IOFFT0+IDIM1
         END DO

      ELSE

c      o evaluate diagonal force constant matrix elements
         DTMP=1.d0/(STPSIZ*STPSIZ)
         E00=-2.d0*DTMP*E0
         IOFFT=1
         DO IDIM=1,NDIM
            FCM(IDIM,IDIM)=E00+DTMP*(ENERGY(IOFFT)+ENERGY(IOFFT+1))
            IOFFT=IOFFT+2*(IDIM+1)
         END DO

c      o compute off-diagonal upper triangle force constants
         E00=-DTMP*E0
         DTMP=0.5d0*DTMP
         IOFFT0=3
         DO IDIM1=2,NDIM
            Z=E00-0.5d0*FCM(IDIM1,IDIM1)
            IOFFT=IOFFT0
            DO IDIM2=1,IDIM1-1
               FCM(IDIM2,IDIM1)=Z+DTMP*(ENERGY(IOFFT)+ENERGY(IOFFT+1))
     &                           -0.5d0*FCM(IDIM2,IDIM2)
               IOFFT=IOFFT+2
            END DO
            IOFFT0=IOFFT0+2*IDIM1
         END DO

      END IF

      if (ndscr.lt.ndim*(4+ndim)) then
         print *, '@ENER2FCM: Insufficient memory.'
         print *, '           need ',ndim*(4+ndim),' doubles'
         print *, '           have ',ndscr,' doubles'
         call aces_exit(1)
      end if

c   o reflect upper triangle into lower triangle (make a copy for DSYEV)
      DO IDIM2=1,NDIM-1
         DSCR(IDIM2,1+IDIM2)=FCM(IDIM2,IDIM2)
         DO IDIM1=IDIM2+1,NDIM
            FCM(IDIM1,IDIM2)=FCM(IDIM2,IDIM1)
            DSCR(IDIM1,1+IDIM2)=FCM(IDIM2,IDIM1)
         END DO
      END DO
      DSCR(NDIM,1+NDIM)=FCM(NDIM,NDIM)

      IF (PRINTQ) THEN
         WRITE(6,1010)
1010     FORMAT(T3,'Force constant matrix : ')
         WRITE(6,1011)((I,J,FCM(I,J),J=1,NDIM),I=1,NDIM)
1011     FORMAT(3('[',I3,',',I3,']',1X,F9.6,1X),'[',I3,',',I3,']',1X,
     &          F9.6)
      END IF

c   o diagonalize force constant matrix
      LWORK = ndscr - ndim*(ndim+1)
      CALL DSYEV('N','L',NDIM,DSCR(1,2),NDIM,DSCR(1,1),
     &                   DSCR(1,NDIM+2),LWORK,I)
      if (i.ne.0) then
         print *, '@ENER2FCM: Eigensolver failed.'
         print *, '           dsyev returned error ',i
         call aces_exit(i)
      end if

c   o print the frequencies
      WRITE(6,1000)
1000  FORMAT(T3,'Vibrational frequencies  (cm-1) :')
      DO IMODE=1,NDIM
         X=DSCR(IMODE,1)
         IF (X.LT.0.d0) THEN
            WRITE(6,1001)IMODE,SQRT(-X)*FACT
1001        FORMAT(T3,I5,1X,F10.5,'i')
         ELSE
            WRITE(6,1002)IMODE,SQRT(X)*FACT
1002        FORMAT(T3,I5,1X,F10.5)
         END IF
      END DO

      RETURN
      END

