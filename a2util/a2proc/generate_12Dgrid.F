       SUBROUTINE GENERATE_12DGRID(WORK, ILEFT, NATOMS, NFCT,
     &                             NBFNS, NAOBFNS, LNP1, IPOPF,
     &                             IUATMS, COORD, ALPHA, PCOEF,
     &                             NUNQSHL, NSHL, NANGMOMSHL,
     &                             NCONFUNSHL, NPRIMFUNSHL,
     &                             NOFFSETATMP, NOFFSETATMC,
     &                             NOFFSETPRM, NOFFSETCON, 
     &                             TMP1, TMP2, TMP3, PRDUCINT,
     &                             IATMCHRG, IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*2 GRID_TYPE, ORIENTATION(3)
      CHARACTER*7 DENSITY_TYPE
      LOGICAL D2_GRID, D3_GRID, D2_EXIST, D3_EXIST, SPIN_D, POST_SCF,
     &        TRIPLET_PERT
      DIMENSION IRREP_PERT(8)
      PARAMETER (MAX_GRID_POINTS = 100)
C
#include "flags.h"
#include "machsp.com"
#include "jodaflags.com"                   
#include "mxatms.par"
C#include "istart.com"                      
#include "baslims.par"
C    
      DIMENSION ALPHA(NFCT), IPOPF(IUATMS), COORD(NATOMS*3), 
     &          NSHL(IUATMS), NANGMOMSHL(IUATMS, NUNQSHL),
     &          NCONFUNSHL(IUATMS, NUNQSHL), 
     &          NPRIMFUNSHL(IUATMS,NUNQSHL), PCOEF(NFCT, NBFNS), 
     &          PRDUCINT(NAOBFNS,NAOBFNS), 
     &          WORK(ILEFT), NOFFSETPRM(IUATMS,NUNQSHL),
     &          NOFFSETCON(IUATMS,NUNQSHL), NOFFSETATMP(IUATMS), 
     &          NOFFSETATMC(IUATMS), TMP1(LNP1,LNP1),TMP2(LNP1,LNP1),
     &          TMP3(LNP1,LNP1), ORIGIN_2D(3),
     &          NCUBE(4), ATOM_XYZ(4, MXATMS), IATMCHRG(NATOMS),
     &          QUANTITY(0:4, MAX_GRID_POINTS), STEP_3D(3, 0:3)
 
      DATA D0, D1, D2, D4 /0.0D0, 1.0D0, 2.0D0, 4.0D0/
      DATA LUNIT_2D, LUNIT_3D /10, 20/
      DATA XTANG /0.529177249D0/ 
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "Entering ZMAT READ"
#endif 
          POST_SCF = .FALSE.
      TRIPLET_PERT = .FALSE.
      If (iflags(h_IFLAGS_calc) .GT. 1)    POST_SCF = .TRUE.
      If (iflags(h_IFLAGS_props) .EQ. 9) TRIPLET_PERT = .TRUE.
C

      CALL A2RDZMAT_4GRIDINF(GRID_TYPE, DENSITY_TYPE, SPIN_D, 
     &                       NMBR_OF_PERTS, IPICK_PERT, STEP_SIZE, 
     &                       D2RANGE, D3RANGE, XORIGIN, YORIGIN, 
     &                       ZORIGIN)
C
      
      IF (DENSITY_TYPE .EQ. "1-ORDER" .OR. 
     &    DENSITY_TYPE .EQ. "DEFINED") THEN
         CALL GETREC(0, 'JOBARC', 'NTOTPERT', IRECLEN, NTPERT)
        IF (IRECLEN .GT. 0) THEN
           CALL GETREC(20, 'JOBARC', 'NTOTPERT', 1, NTPERT)
        ELSE 
           NTPERT = 0
           CALL IZERO(IRREP_PERT, 8)
           CALL GETREC (20, 'JOBARC', 'NPERTIRR', 8, IRREP_PERT)
           DO IREPS = 1, 8
              NTPERT = NTPERT + IRREP_PERT(IREPS)
           END DO
        END IF
C    
        IF (NMBR_OF_PERTS .GT. NTPERT .OR. IPICK_PERT .GT. NTPERT) 
     &     THEN
           Write(6,'(t3,6(a))') "Inconsistency in the number of",
     &           " internal and input total perturbations",
     &           " or the selected perturbation is beyond the",
     &           " total number of perturbations: the internal", 
     &           " total number of perturbations and the first", 
     &           " perturbation in the list is choosen."
           NMBR_OF_PERTS = NTPERT
           IPICK_PERT    = 1
        ENDIF
      ENDIF
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Write(6,'(1x,a,a,2I3)') "The total number of perturbations ",
     & "and the perturbation of interest", NMBR_OF_PERTS, 
     & IPICK_PERT
#endif
      CALL A2GET_DEN(WORK, ILEFT, DENSITY_TYPE, SPIN_D, 
     &               NMBR_OF_PERTS, NBFNS, NAOBFNS, ISCF_TDEN,
     &               ICOR_TDEN, ISCF_DDEN, ICOR_DDEN, IBEGIN_P_DENS,
     &               IUHF)
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "Offsets for 0/1-order density"
      Write(*,'(7(1x,i5))') ISCF_TDEN,ICOR_TDEN,ISCF_DDEN,ICOR_DDEN,
     &                       IBEGIN_P_DENS, NMBR_OF_PERTS, IPICK_PERT
      Write(6,*)
#endif
C
      IF (GRID_TYPE .EQ. "2D") THEN
         D2_GRID = .TRUE.
         ORIENTATION(1) = 'xy'
         ORIENTATION(2) = 'xz'
         ORIENTATION(3) = 'yz'
         INQUIRE(FILE="2DPLOT", EXIST=D2_EXIST)
         IF (D2_EXIST) THEN
            OPEN(UNIT=LUNIT_2D, FILE="2DPLOT", FORM="FORMATTED", 
     &           STATUS='OLD')
         ELSE
            OPEN(UNIT=LUNIT_2D, FILE="2DPLOT", FORM="FORMATTED", 
     &           STATUS='NEW')
         ENDIF
          
      ELSE 
         D3_GRID = .TRUE.
         INQUIRE(FILE="3DPLOT", EXIST=D3_EXIST)
         IF (D3_EXIST) THEN
            OPEN(UNIT=LUNIT_3D, FILE="3DPLOT", FORM="FORMATTED", 
     &           STATUS='OLD')
         ELSE
            OPEN(UNIT=LUNIT_3D, FILE="3DPLOT", FORM="FORMATTED", 
     &           STATUS='NEW')
         ENDIF
C
         IOFF = 1
         DO IATOM =1, NATOMS
            IOFF = IOFF + (IATOM - 1)*3
            DO IXYZ = 1, 3
               ATOM_XYZ(IXYZ, IATOM) = COORD(IOFF)
               IOFF = IOFF + 1
            ENDDO
            IOFF = 1
         END DO
C
         NCUBE(1) = 80
         NCUBE(2) = 80
         NCUBE(3) = 80
C
#ifdef _DEBUG_LVLM2
         NCUBE(1) = 2
         NCUBE(2) = 2
         NCUBE(3) = 2
#endif
C
       ENDIF 

C#ifdef _NOSKIP

      IF(D2_GRID) THEN 

        STEP_2D      = STEP_SIZE
        ADD_2D       = D2RANGE
        ORIGIN_2D(1) = XORIGIN/XTANG
        ORIGIN_2D(2) = YORIGIN/XTANG
        ORIGIN_2D(3) = ZORIGIN/XTANG
        NR_STEPS     = 1 + ADD_2D*D2/STEP_2D

        DO IORIENTATION = 1,3

          NI = NR_STEPS
          NJ = NR_STEPS
          NK = NR_STEPS
C
#ifdef _DEBUG_LVLM2
          NI = 5
          NJ = 5
          NK = 5
#endif
          IF (ORIENTATION(IORIENTATION) .EQ. 'xy') NK = 1
          IF (ORIENTATION(IORIENTATION) .EQ. 'xz') NJ = 1
          IF (ORIENTATION(IORIENTATION) .EQ. 'yz') NI = 1

          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A)')
          WRITE(LUNIT_2D,'(A2)') ORIENTATION(IORIENTATION)
          WRITE(LUNIT_2D,'(A)')

          DO I = 1,NI
            DO J = 1,NJ
              DO K = 1,NK

                XPOINT = ORIGIN_2D(1)
     &                   - ADD_2D + (I-1)*STEP_2D
                YPOINT = ORIGIN_2D(2)
     &                   - ADD_2D + (J-1)*STEP_2D
                ZPOINT = ORIGIN_2D(3)
     &                   - ADD_2D + (K-1)*STEP_2D

                IF(ORIENTATION(IORIENTATION) .EQ. 'xy')
     &                                     ZPOINT = ORIGIN_2D(3)
                IF(ORIENTATION(IORIENTATION) .EQ. 'xz')
     &                                     YPOINT = ORIGIN_2D(2)
                IF(ORIENTATION(IORIENTATION) .EQ. 'yz')
     &                                     XPOINT = ORIGIN_2D(1)

                CALL A2EVAL_INT_MAIN(IUATMS, NATOMS, NUNQSHL, LNP1, 
     &                               IPOPF, NSHL, NFCT, NBFNS, 
     &                               NOFFSETATMP, NOFFSETATMC, 
     &                               NOFFSETPRM, NOFFSETCON,
     &                               NANGMOMSHL, NPRIMFUNSHL, 
     &                               NCONFUNSHL, XPOINT, YPOINT, 
     &                               ZPOINT, ALPHA, PCOEF, COORD, 
     &                               PRDUCINT, TMP1, TMP2, TMP3)
#ifdef _DEBUG_LVLM2
      Print*, "Offsets for 0-order density"
      Write(*,'(7(1x,i5))') ISCF_TDEN,ICOR_TDEN,ISCF_DDEN,ICOR_DDEN,
     &                       IBEGIN_P_DENS, NMBR_OF_PERTS, IPICK_PERT
      Write(6,*)
      Print*, "The product integrals"
      CALL OUTPUT(PRDUCINT, 1, NAOBFNS, 1, NAOBFNS, NAOBFNS,
     &            NAOBFNS, 1)

#endif
                CALL A2BUILD_QUANTITY(WORK, ILEFT, PRDUCINT, 
     &                                QUANTITY, SPIN_D, DENSITY_TYPE, 
     &                                NMBR_OF_PERTS, IPICK_PERT,
     &                                NBFNS, NAOBFNS, ISCF_TDEN, 
     &                                ICOR_TDEN, ISCF_DDEN, 
     &                                ICOR_DDEN, IBEGIN_P_DENS, 
     &                                MAX_GRID_POINTS, POST_SCF, 
     &                                GRID_TYPE, 1, IUHF)
C
               XPOINT        = XPOINT*XTANG
               YPOINT        = YPOINT*XTANG
               ZPOINT        = ZPOINT*XTANG
               QUANTITY(0,1) = QUANTITY(0,1)
               QUANTITY(1,1) = QUANTITY(1,1)
              
C
               IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
C
                   IF (POST_SCF) THEN
                       WRITE(LUNIT_2D,'(5F20.10)') XPOINT, YPOINT, 
     &                                             ZPOINT,
     &                                             QUANTITY(0,1),
     &                                             QUANTITY(1,1)
                   ELSE
                       WRITE(LUNIT_2D,'(4F20.10)') XPOINT, YPOINT,
     &                                             ZPOINT,
     &                                             QUANTITY(0,1)
                   ENDIF
C
               ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN
                        WRITE(LUNIT_2D,'(5F20.10)') XPOINT, YPOINT,
     &                                              ZPOINT,
     &                                              QUANTITY(0,1),
     &                                              QUANTITY(1,1)
C
              ELSE IF (DENSITY_TYPE .EQ. "DEFINED") THEN
C
                       WRITE(LUNIT_2D,'(5F20.10)') XPOINT, YPOINT, 
     &                                             ZPOINT, 
     &                                             QUANTITY(0,1)
               ENDIF
  
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF(D3_GRID) THEN 

C calculate cube parameters

C NCUBE(0)  number of atoms
C NCUBE(1)  number of points in the "x"-direction
C NCUBE(2)  number of points in the "y"-direction
C NCUBE(3)  number of points in the "z"-direction

C STEP_3D(*,0)  coordinates of initial point
C STEP_3D(*,1)  step vector for the "x"-direction
C STEP_3D(*,2)  step vector for the "y"-direction
C STEP_3D(*,3)  step vector for the "z"-direction

C ATOM_XYZ(1,I)  x-coordinate of atom I
C ATOM_XYZ(2,I)  y-coordinate of atom I
C ATOM_XYZ(3,I)  z-coordinate of atom I
C ATOM_XYZ(4,I)  charge of atom I

        NCUBE(0) = NATOMS
        ADD_3D   = D3RANGE
        CALL DZERO(STEP_3D, 12)
C
        DO I = 1,3
          STEP_3D(I,0) = ATOM_XYZ(I,1)
          STEP_3D(I,I) = ATOM_XYZ(I,1)
        ENDDO
C
        DO I = 1, NATOMS
           DO M = 1, 3
              STEP_3D(M,0)  = MIN(STEP_3D(M,0),ATOM_XYZ(M,I))
              STEP_3D(M,M)  = MAX(STEP_3D(M,M),ATOM_XYZ(M,I))
           END DO
              ATOM_XYZ(4,I) = DBLE(IATMCHRG(I))
        ENDDO

C add some space, calculate step sizes and convert step vectors to
C angstrom
C
        DO I = 1,3
          STEP_3D(I,0) = STEP_3D(I,0) - ADD_3D
          STEP_3D(I,I) = (STEP_3D(I,I) + (ADD_3D*D2))
          STEP_3D(I,I) = (STEP_3D(I,I) - STEP_3D(I,0))/(NCUBE(I) - 1)
        ENDDO

C Write header

        WRITE(LUNIT_3D,'(A)') 'some'
        WRITE(LUNIT_3D,'(A)') 'text'
        DO I = 0,3
          WRITE(LUNIT_3D,'(I5,4F12.6)') NCUBE(I), (STEP_3D(J,I),J=1,3)
        ENDDO
        DO I = 1,NATOMS
          WRITE(LUNIT_3D,'(I5,4F12.6)')
     &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
        ENDDO

        DO ICUBE = 1,NCUBE(1)
          DO JCUBE = 1,NCUBE(2)
            DO KCUBE = 1,NCUBE(3)

              XPOINT = STEP_3D(1,0) + (ICUBE-1)*STEP_3D(1,1)
     &                              + (JCUBE-1)*STEP_3D(1,2)
     &                          + (KCUBE-1)*STEP_3D(1,3)
              YPOINT = STEP_3D(2,0) + (ICUBE-1)*STEP_3D(2,1)
     &                          + (JCUBE-1)*STEP_3D(2,2)
     &                          + (KCUBE-1)*STEP_3D(2,3)
              ZPOINT = STEP_3D(3,0) + (ICUBE-1)*STEP_3D(3,1)
     &                          + (JCUBE-1)*STEP_3D(3,2)
     &                          + (KCUBE-1)*STEP_3D(3,3)
C
              CALL A2EVAL_INT_MAIN(IUATMS, NATOMS, NUNQSHL, LNP1, IPOPF,
     &                             NSHL, NFCT, NBFNS, NOFFSETATMP,
     &                             NOFFSETATMC, NOFFSETPRM, NOFFSETCON,
     &                             NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL,
     &                             XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF,
     &                             COORD, PRDUCINT, TMP1, TMP2, TMP3)
#ifdef _DEBUG_LVLM2
      Print*, "Offsets for 0-order density"
      Write(*,'(7(1x,i5))') ISCF_TDEN,ICOR_TDEN,ISCF_DDEN,ICOR_DDEN,
     &                       IBEGIN_P_DENS, NMBR_OF_PERTS, IPICK_PERT
      Write(6,*)
      Print*, "The product integrals"
      CALL OUTPUT(PRDUCINT, 1, NAOBFNS, 1, NAOBFNS, NAOBFNS,
     &            NAOBFNS, 1)
      
#endif
               CALl A2BUILD_QUANTITY(WORK, ILEFT, PRDUCINT,
     &                               QUANTITY, SPIN_D, DENSITY_TYPE,
     &                               NMBR_OF_PERTS, IPICK_PERT,
     &                               NBFNS, NAOBFNS, ISCF_TDEN,
     &                               ICOR_TDEN, ISCF_DDEN,
     &                               ICOR_DDEN, IBEGIN_P_DENS,
     &                               MAX_GRID_POINTS, POST_SCF, 
     &                               GRID_TYPE, KCUBE, IUHF)
C
            ENDDO

               XPOINT = XPOINT*XTANG
               YPOINT = YPOINT*XTANG
               ZPOINT = ZPOINT*XTANG
C
               IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
C
                    IF (TRIPLET_PERT) THEN
                        WRITE(LUNIT_3D,'(6E13.5)')  (QUANTITY(0,M),
     &                                               M = 1, NCUBE(3))
                    ELSE
                        WRITE(LUNIT_3D,'(6E13.5)')  (QUANTITY(1,M),
     &                                               M = 1, NCUBE(3))
                    ENDIF
                   
C
               ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN
C
                    IF (TRIPLET_PERT) THEN
                        WRITE(LUNIT_3D,'(6E13.5)') (QUANTITY(1,M),
     &                                              M = 1, NCUBE(3))
                    ELSE
                        WRITE(LUNIT_3D,'(6E13.5)') (QUANTITY(0,M),
     &                                              M = 1, NCUBE(3))
                    END IF

               ELSE IF (DENSITY_TYPE .EQ. "DEFINED") THEN
C
                       WRITE(LUNIT_3D,'(6E13.5)') (QUANTITY(0,M),
     &                                             M = 1, NCUBE(3))
               ENDIF
C
          ENDDO
        ENDDO
C
      ENDIF 
C
C#endif 
#ifdef _DEBUG_LVLM2
       Write(6,*)
       print*, "Entering A2EVAL_PRDCTINT"
       XPOINT=1.0D0
       YPOINT=1.0D0
       ZPOINT=1.0D0
#ifdef _DEBUG_LVL0
       Print*, "Input varaibles:", Iuatms, Natoms, nunqshl, lnp1, nbfns,
     &          nfct 
       Write(6,*)
       write(*, '(a)'), "The exponents"
       Write(*, '(4F10.5)') (alpha(i), i=1, NFCT)
       Write(6,*)
       Write(*, '(a)'), "Cartesian Coords."
       Write(*, '(4F10.5)') (Coord(i), i=1, 3*NATOMS)
       Write(6,*)
       Write(*, '(a)'), "The NSHL array"
       Write(*, '(4I5)') (Nshl(i), i=1, iuatms)
       Write(6,*)
       Write(6,*) "The ipopf array"
       Write(*, '(4I5)') (ipopf(i), i=1, iuatms)
       Write(6,*)
       Write(6, *) "The NOFFSETATMC and NOFFSETATMP arrays"
       Write(*, '(4I5)') (NOFFSETATMC(i), i=1, iuatms)
       Write(*, '(4I5)') (NOFFSETATMP(i), i=1, iuatms)
       write(6,*)
       Write(6,*) "The number of primitive in a shell"
       Write(6, '(4I5)') ((NPRIMFUNSHL(i,j), I=1, iuatms), 
     &                    j=1,NUNQSHL)
       write(6,*)
       Write(*, '(a)') "The Contractions Coef."
       write(*, '(4F10.5)') ((pcoef(i, j), i=1, nfct), 
     &                        j=1,nbfns)
#endif
   
              CALL A2EVAL_INT_MAIN(IUATMS, NATOMS, NUNQSHL, LNP1, IPOPF,
     &                             NSHL, NFCT, NBFNS, NOFFSETATMP,
     &                             NOFFSETATMC, NOFFSETPRM, NOFFSETCON,
     &                             NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL,
     &                             XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF,
     &                             COORD, PRDUCINT, TMP1, TMP2, TMP3)
              CALL A2BUILD_QUANTITY(WORK, ILEFT, PRDUCINT,
     &                              QUANTITY, SPIN_D, DENSITY_TYPE,
     &                              NMBR_OF_PERTS, IPICK_PERT,
     &                              NBFNS, NAOBFNS, ISCF_TDEN,
     &                              ICOR_TDEN, ISCF_DDEN,    
     &                              ICOR_DDEN, IBEGIN_P_DENS,
     &                              MAX_GRID_POINTS, POST_SCF,
     &                              GRID_TYPE, 1, IUHF)      
C     
#endif

      CLOSE (LUNIT_2D, STATUS="KEEP")
      CLOSE (LUNIT_3D, STATUS="KEEP")
C
      Return
      END
		
		
	
