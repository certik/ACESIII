
      SUBROUTINE A2GET_DEN(WORK, MAXCOR, DENSITY_TYPE, SPIN_D, 
     &                     NMBR_OF_PERTS, NBFNS, NAOBFNS, ISCF_TDEN,
     &                     ICOR_TDEN, ISCF_DDEN, ICOR_DDEN, 
     &                     IBEGIN_P_DENS, IUHF) 
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPHERICAL, SPIN_D
C
      CHARACTER*7 DENSITY_TYPE
      CHARACTER*8 LABELPSCF, LABELDSCF, LABELDENS, LABELSDEN,
     &            RECNAMEA, RECNAMEB
      DIMENSION WORK(MAXCOR)
C
      DATA IONE /1/
C
#include "machsp.com"
#include "syminf.com"
#include "jodaflags.com"
C
      SPHERICAL = (IFLAGS(62) .EQ. 1)
      LABELPSCF ='HFDENSTY'
      LABELDSCF ='HDDENSTY'
      LABELDENS ='TDENSITY'
      LABELSDEN ='DDENSITY'
      NAOBFNS2  = NAOBFNS*NAOBFNS
C
      IBEGIN_P_DENS = IONE
      ISCF_TDEN     = IONE
      ISCF_DDEN     = IONE
      ICOR_TDEN     = IONE
      ICOR_DDEN     = IONE
C
      IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
      
C
C Get the SCF and correlated (relaxed) density matirices in AO basis. For 
C open shell calculations, we also read the density diferences. Note that
C the "relaxed" density read here only contain the correlated contributions.
C
      ISCF_TDEN = IONE
      ICOR_TDEN = ISCF_TDEN + NAOBFNS2
      INEXT     = ICOR_TDEN + NAOBFNS2
C
#ifdef _DEBUG_LVLM2
      Print*, "Offsets for 0-order density"
      Write(*,'(4(1x,i5))') ISCF_TDEN, ICOR_TDEN, INEXT
      Write(6,*)
#endif
      IF (INEXT .GT. MAXCOR) CALL INSMEM("A2GET_DEN", INEXT, MAXCOR)

      CALL GETREC(-20, 'JOBARC', LABELPSCF, NAOBFNS2*IINTFP,
     &             WORK(ISCF_TDEN))
      CALL GETREC(-20, 'JOBARC', LABELDENS, NAOBFNS2*IINTFP, 
     &             WORK(ICOR_TDEN))
C      
      IF (SPIN_D .AND. (IUHF .GT. 0)) THEN
         
          ISCF_DDEN = INEXT
          ICOR_DDEN = ISCF_DDEN + NAOBFNS2
          INEXT     = ICOR_DDEN + NAOBFNS2
#ifdef _DEBUG_LVLM2
      Print*, "Offsets for 0-order density"
      Write(*,'(3(1x,i5))') ISCF_TDEN, ICOR_TDEN, INEXT 
      Write(6,*)
#endif
          IF (INEXT .GT. MAXCOR) CALL INSMEM("A2GET_DEN", INEXT, 
     &                                        MAXCOR)
 
          CALL GETREC(-20, 'JOBARC', LABELDSCF, NAOBFNS2*IINTFP, 
     &                 WORK(ISCF_DDEN))
          CALL GETREC(-20, 'JOBARC', LABELSDEN, NAOBFNS2*IINTFP, 
     &                 WORK(ICOR_DDEN))

      ENDIF
C      
#ifdef _DEBUG_LVLM2
      Write(6,*) "The SCF and Cor. tot. density matrices in AO basis"
      CALL OUTPUT(WORK(ISCF_TDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      CALL OUTPUT(WORK(ICOR_TDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
C
      IF (SPIN_D .AND. IUHF .GT. 0) THEN
      Write(6,*) "The SCF and Cor. spin. density matrices in AO basis"
      CALL OUTPUT(WORK(ISCF_DDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      CALL OUTPUT(WORK(ICOR_DDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      ENDIF
#endif
C
      ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN
C
           MAX_MEM_NEED = 2*NAOBFNS2*NMBR_OF_PERTS + 2*NAOBFNS2
           IF (MAXCOR .LE. MAX_MEM_NEED) CALL INSMEM("A2GET_DEN", 
     &                                                MAX_MEM_NEED, 
     &                                                MAXCOR)
C
           IKEEP_ONE_DENA = IONE 
           IKEEP_ONE_DENB = IKEEP_ONE_DENA + NAOBFNS2
           IKEEP_P_TDEN_N = IKEEP_ONE_DENB + NAOBFNS2
           IBEGIN_P_DENS  = IKEEP_P_TDEN_N  

           DO IPERT = 1, NMBR_OF_PERTS
C
              IKEEP_P_DDEN_N =  IKEEP_P_TDEN_N + NAOBFNS2
              INEXT          =  IKEEP_P_DDEN_N + NAOBFNS2
C
#ifdef _DEBUG_LVLM2
      Print*, "Offsets for 1-order density" 
      Write(*,'(a,1x,i5)')  "The number of pert.", NMBR_OF_PERTS
      Write(*,'(4(1x,i5))') IKEEP_ONE_DENA, IKEEP_ONE_DENB, 
     &                         IKEEP_P_TDEN_N, IKEEP_P_DDEN_N
      Write(6,*)
#endif
              WRITE(RECNAMEA, "(A,I2)") "PTDENA", IPERT
              CALL GETREC(-20, "JOBARC", RECNAMEA, NAOBFNS2*IINTFP, 
     &                    WORK(IKEEP_ONE_DENA))
C
              WRITE(RECNAMEB, "(A,I2)") "PTDENB", IPERT
              CALL GETREC(-20, "JOBARC", RECNAMEB, NAOBFNS2*IINTFP,
     &                     WORK(IKEEP_ONE_DENB))
C
#ifdef _DEBUG_LVLM2
      Write(6,'(1x,a,a,1x,a)') "The record labels ", RECNAMEA, 
     &      RECNAMEB
      Write(6,*)
      Write(6,*) "The  1-st order alpha density matrices in AO basis"
      CALL OUTPUT(WORK(IKEEP_ONE_DENA), 1, NAOBFNS, 1, NAOBFNS,
     &            NAOBFNS, NAOBFNS, 1)
C   
      Write(6,*)
      Write(6,*) "The 1-st order beta density matrices in AO basis"
      CALL OUTPUT(WORK(IKEEP_ONE_DENB), 1, NAOBFNS, 1, NAOBFNS,
     &            NAOBFNS, NAOBFNS, 1)
           
#endif   
              CALL ZERO(WORK(IKEEP_P_TDEN_N), NAOBFNS2)
              CALL ZERO(WORK(IKEEP_P_DDEN_N), NAOBFNS2)
C
C Build perturbed total density for n perturbation 
C
              CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENA), 1, 
     &                   WORK(IKEEP_P_TDEN_N), 1) 
              CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENB), 1, 
     &                   WORK(IKEEP_P_TDEN_N), 1) 
C
C Build perturbed spin density for n perturbation 
C
              CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENA), 1, 
     &                   WORK(IKEEP_P_DDEN_N), 1) 
              CALL SAXPY(NAOBFNS2, -1.0D0, WORK(IKEEP_ONE_DENB), 1, 
     &                   WORK(IKEEP_P_DDEN_N), 1) 
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Write(6,*) "The  1-st order tot. density matrices in AO basis"
      CALL OUTPUT(WORK(IKEEP_P_TDEN_N), 1, NAOBFNS, 1, NAOBFNS, 
     &            NAOBFNS, NAOBFNS, 1)
C    
      Write(6,*)
      Write(6,*) "The 1-st order spin. density matrices in AO basis"
      CALL OUTPUT(WORK(IKEEP_P_DDEN_N), 1, NAOBFNS, 1, NAOBFNS,
     &            NAOBFNS, NAOBFNS, 1)

#endif
              IKEEP_P_TDEN_N =  INEXT
C
           ENDDO
C
C endif for density type
C
      ENDIF
C
C At this point we need to transform the density matrices to Cartesian
C basis in order to directly contract with the Cartesian AO integrals.
C 
C$$$      IF (SPHERICAL) THEN
C$$$         CALL GETREC(20, 'JOBARC','CMP2CART', NBAS*NBASP*IINTFP, 
C$$$     &               SPH2CART)

C$$$      Write(6,*) "The spherical to Cartesian transformation matrix"
C$$$      CALL OUTPUT(SPH2CART, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
C
C$$$         CALL TRANS_SPH2CART(SCFDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$         CALL TRANS_SPH2CART(RELDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$      ENDIF
C
C
      RETURN
      END


