      SUBROUTINE A2BUILD_QUANTITY(WORK, MAXCOR, PRDUTINT, QUANTITY, 
     &                            SPIN_D, DENSITY_TYPE, 
     &                            NMBR_OF_PERTS, IPICK_PERT, NBFNS,
     &                            NAOBFNS, ISCF_TDEN, ICOR_TDEN,  
     &                            ISCF_DDEN, ICOR_DDEN, 
     &                            IBEGIN_P_DENS, MAX_GRID_POINTS, 
     &                            POST_SCF, GRID_TYPE, KCUBE, IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPIN_D, POST_SCF
C
      CHARACTER*7 DENSITY_TYPE
      CHARACTER*2 GRID_TYPE
C
      DIMENSION WORK(MAXCOR), PRDUTINT(NAOBFNS, NAOBFNS),
     &          QUANTITY(0:4, MAX_GRID_POINTS)
      DATA TWO /2.0D0/
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "Offsets for 0/1-order density"
      Write(*,'(7(1x,i5))') ISCF_TDEN,ICOR_TDEN,ISCF_DDEN,ICOR_DDEN,
     &                      IBEGIN_P_DENS, KCUBE, MAX_GRID_POINTS,
     &                      NMBR_OF_PERTS, IPICK_PERT, IBEGIN_P_DENS
      Write(6,*)
#endif
#ifdef _DEBUG_LVLM2

      IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
         IF (.NOT. SPIN_D) THEN
            IOFFTS = ISCF_TDEN
            IF (POST_SCF) IOFFTC = ICOR_TDEN
            Write(6,*)
            Write(6,*) "@-quantity the total density matrices",
     &                 " in AO basis"
             CALL OUTPUT(WORK(IOFFTS), 1, NAOBFNS, 1, NAOBFNS,
     &                  NAOBFNS, NAOBFNS, 1)
             Write(6,*)
             IF (POST_SCF) CALL OUTPUT(WORK(IOFFTC), 1, NAOBFNS,  
     &                                 1, NAOBFNS, NAOBFNS,  
     &                                 NAOBFNS, 1)
         ELSE
            IOFFDS = ISCF_DDEN
            IF (POST_SCF) IOFFDC = ICOR_DDEN
            Write(6,*)
            Write(6,*) "@-quantity the spin density matrices",
     &                 " in AO basis"
             CALL OUTPUT(WORK(IOFFDS), 1, NAOBFNS, 1, NAOBFNS,
     &                  NAOBFNS, NAOBFNS, 1)
             Write(6,*)        
             IF (POST_SCF) CALL OUTPUT(WORK(IOFFDC), 1, NAOBFNS,  
     &                                 1, NAOBFNS, NAOBFNS,  
     &                                 NAOBFNS, 1)                    
         ENDIF
      ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN
            IOFFT = IBEGIN_P_DENS + 2*(IPICK_PERT-1)*
     &              NAOBFNS*NAOBFNS 
            IOFFD = IOFFT + NAOBFNS*NAOBFNS
            ISTOP = IOFFD + NAOBFNS*NAOBFNS
            Write(6,*) "@-quantity the total density matrices",
     &                 " in AO basis", IOFFT, IOFFD
             CALL OUTPUT(WORK(IOFFT), 1, NAOBFNS, 1, NAOBFNS,
     &                  NAOBFNS, NAOBFNS, 1)
            write(6,*)
            Write(6,*) "@-quantity the spin density matrices",
     &                 " in AO basis"
             CALL OUTPUT(WORK(IOFFD), 1, NAOBFNS, 1, NAOBFNS,
     &                  NAOBFNS, NAOBFNS, 1)
      ENDIF
      Write(6,*)
      Write(6,*) "@-quantity the product integral"
      CALL OUTPUT(PRDUTINT, 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
#endif
      IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
         SCF_TOT_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                          WORK(ISCF_TDEN), 1)
         COR_TOT_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                         WORK(ICOR_TDEN), 1)
C
         IF (.NOT. SPIN_D) THEN
C
CSSS            DO IND_CUBE = 1, KCUBE
               QUANTITY(0, KCUBE) = SCF_TOT_DENSITY*TWO
               IF (POST_SCF)  QUANTITY(1, KCUBE) = 
     &                                COR_TOT_DENSITY*TWO
CSSS            END DO
         ENDIF
C
         IF (IUHF .NE. 0 .AND. SPIN_D) THEN
            SCF_DIF_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT,
     &                             1, WORK(ISCF_DDEN), 1)
            COR_DIF_DENSITY = DDOT(NAOBFNS*NAOBFNS, PRDUTINT,
     &                             1, WORK(ICOR_DDEN), 1)
C
CSSS            DO IND_CUBE = 1, KCUBE
               QUANTITY(0, KCUBE) = SCF_DIF_DENSITY*TWO
               IF (POST_SCF)  QUANTITY(1, KCUBE) = 
     &                                 COR_DIFF_DENSITY*TWO
CSSS            END DO
         ENDIF
C
#ifdef _DEBUG_LVLM2
      Print*, "Quantity of interest (scf, correlated total or spin",
     &        " density)"
CSSS      Do IND_CUBE =1, KCUBE
      Write(*,'((1x,2F10.6))') QUANTITY(0, KCUBE), 
     &                         QUANTITY(1,KCUBE)
CSSS      Enddo
#endif
      ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN

          DO IPERT = 1, NMBR_OF_PERTS
             IF (IPERT .EQ. IPICK_PERT) THEN
            
                I_TOT_PDEN_N = IBEGIN_P_DENS + 2*(IPERT-1)*
     &                         NAOBFNS*NAOBFNS
                I_DIF_PDEN_N =  I_TOT_PDEN_N + NAOBFNS*NAOBFNS

                TOT_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                            WORK(I_TOT_PDEN_N), 1)
                DIF_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                            WORK(I_DIF_PDEN_N), 1)
#ifdef _DEBUG_LVLM2
      Write(6,*)
      write(6,'(a,a,I2)'), "Computing the perturbed total/spin", 
     &                    " density for perturbation: ", ipert
      Write(6,*)
      Write(*,'(3(1x,i5),2F10.5)') I_TOT_PDEN_N, I_DIF_PDEN_N,  
     &                             KCUBE, TOT_PDEN_N, DIF_PDEN_N 
      Write(6,*)
#endif
C
CSSS                DO IND_CUBE = 1, KCUBE
                   QUANTITY(0, KCUBE) = TOT_PDEN_N*TWO
                   QUANTITY(1, KCUBE) = DIF_PDEN_N*TWO
CSSS                END DO
C
             END IF
          END DO
c
#ifdef _DEBUG_LVLM2
      Print*, "Quantity of interest (perturbed total and/or spin",
     &        " density)"
      Write(6,*), KCUBE
CSS      Do IND_CUBE =1, KCUBE
      Write(*,'((1x,2F10.6))') QUANTITY(0,KCUBE),
     &                         QUANTITY(1,KCUBE)
CSSS      Enddo 
#endif  
      ELSE IF (DENSITY_TYPE .EQ. "DEFINED") THEN

         OPEN (UNIT=30, FILE='VPOUT', FORM='UNFORMATTED', 
     &         STATUS='OLD')
         IRWND = 0
         IERR  = 0
         IOFF  = 1
C
         DO IPERT = 1, NMBR_OF_PERTS

            CALL SEEKLB ('   DEN  ', IERR, IRWND, 30)
            IF (IERR .NE. 0) CALL ERREX
C
            CALL LOAD_DEFINED(WORK(IOFF), NAOBFNS*NAOBFNS, NAOBFNS)
C
C#ifdef _DEBUG_LVLM2
            CALL HEADER ('Matrix elements of the defined operator', 1
     &                   , 6)
            CALL OUTPUT(WORK(IOFF), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &                  NAOBFNS, 1)
C#endif 
            IRWND = 1
            IOFF  = IOFF + NAOBFNS*NAOBFNS
            
         END DO
C
         DO IPERT = 1, NMBR_OF_PERTS
             IF (IPERT .EQ. IPICK_PERT) THEN

                I_TOT_PDEN_N = IBEGIN_P_DENS + (IPERT-1)*
     &                         NAOBFNS*NAOBFNS

                TOT_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                            WORK(I_TOT_PDEN_N), 1)
#ifdef _DEBUG_LVLM2
      Write(6,*)
      write(6,'(a,a,I2)'), "Computing the "defined"",
     &                    " density for perturbation: ", ipert
      Write(6,*)
      Write(*,'(3(1x,i5),2F10.5)') I_TOT_PDEN_N, KCUBE, TOT_PDEN_N
      Write(6,*)
#endif
C
CSSS                DO IND_CUBE = 1, KCUBE
                   QUANTITY(0, KCUBE) = TOT_PDEN_N*TWO
CSSS                END DO
C
             END IF
          END DO
c
#ifdef _DEBUG_LVL0
      Print*, "Quantity of interest (perturbed total and/or spin",
     &        " density)"
      Write(6,*), KCUBE
CSS      Do IND_CUBE =1, KCUBE
      Write(*,'((1x,2F10.6))') QUANTITY(0,KCUBE),
     &                         QUANTITY(1,KCUBE)
CSSS      Enddo
#endif
      ENDIF
C
      RETURN
      END 
C  
