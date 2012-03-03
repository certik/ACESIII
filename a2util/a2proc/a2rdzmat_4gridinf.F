      subroutine A2RDZMAT_4GRIDINF(GRID_TYPE, DENSITY_TYPE, SPIN_D,
     &                             NMBR_OF_PERTS, IPICK_PERT, 
     &                             STEP_SIZE, D2RANGE, D3RANGE, 
     &                             XORIGIN, YORIGIN, ZORIGIN)
C
      implicit double precision (a-h,o-z)
      logical ZMAT_PRESENT, SPIN_D
      character*80 wrk
      CHARACTER*2 GRID_TYPE
      CHARACTER*7 DENSITY_TYPE
      
      INQUIRE(FILE='ZMAT',EXIST=ZMAT_PRESENT)
C
      IF (ZMAT_PRESENT)THEN
         OPEN(UNIT=4,FILE='ZMAT',FORM='FORMATTED',STATUS='OLD')
         rewind(4)
      else
         write(6,'(T3,a)')'@RDZMAT_4GRID.f, ZMAT file not present.'
         call errex
      ENDIF
c
300   READ(4,'(A)', END=900) WRK
      IF (WRK(1:9) .NE.'*GRID_INF') goto 300
C
C The density type can be zeroth order, first-order ...Specify
C the density type as 0-ORDER, 1-ORDER. 
 
      READ(4,'(A)', END=900) WRK
      DENSITY_TYPE = WRK(1:7)
      SPIN_D = .FALSE.
      IF (WRK(9:21) .EQ. "SPIN-DENSITY") SPIN_D = .TRUE.
C
C The number of perturbations. 0 for 1-ORDER and 1-ORDER depends
C on the type of pertubation. For example, 3 for dipole moments.

      READ(4,*,END=900) NMBR_OF_PERTS, IPICK_PERT
C
C The gridtype can be two or three dimensional. Specify the 
C grid type as GRID_TYPE={2D, 3D}

      READ(4,'(A)', END=900) WRK
      GRID_TYPE = WRK(1:2)
C
C The step size for the grid and the range (if range is 
C specifies as x, the grid would be generated -x:x).
C       
      READ(4,*,END=900) STEP_SIZE, D2RANGE, D3RANGE
C
C The orgin of the grid
C       
      READ(4,*,END=900) XORIGIN, YORIGIN, ZORIGIN
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "The user given grid info"
      Write(6, "(t3,a,1x,a,F5.2,1X,i2,1x,i2,1x,5F5.2,1X,l)"), 
     &           GRID_TYPE, DENSITY_TYPE, STEP_SIZE, NMBR_OF_PERTS, 
     &           IPICK_PERT, D2RANGE, D3RANGE, XORIGIN, YORIGIN, 
     &           ZORIGIN, SPIN_D
      Write(6,*)
#endif
C
      GO TO 99
900   WRITE(6,901)
901   FORMAT(T3,'@RDZMAT_4GRID.f, *GRID_INF namelist not found or', 
     &       ' incomplete.')
      GO TO 999 
 999  CONTINUE
      CLOSE(UNIT=4,STATUS='KEEP')
      CALL ERREX
C
  99  CONTINUE
      CLOSE(UNIT=4,STATUS='KEEP')
      
      RETURN
      END
