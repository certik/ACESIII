      subroutine A2get_gridnmlist(GRID_FILE, NMBROF_GRIDPTS)
C
      implicit double precision (a-h,o-z)
      character*80 wrk
      Character*12 GRID_FILE
      
      OPEN(UNIT=4,FILE='ZMAT',FORM='FORMATTED',STATUS='OLD')
      REWIND(4)
C
300   READ(4,'(A)', END=900) WRK
      IF (WRK(1:9) .NE.'*PES_GRID') goto 300
C
C Obtain the number grid points
C
      READ(4,*,END=900) NMBROF_GRIDPTS
C
C The name of the file that the GRID is stored (maximum of 10 
C charaters.

      READ(4,'(A)', END=900) WRK
      GRID_FILE = WRK(1:12)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Print*, "The user given grid info"
      Write(6, "(t3,a,1x,i5)"), 
     &           GRID_FILE, NMBROF_GRIDPTS
      Write(6,*)
#endif
C
      GO TO 99
900   WRITE(6,901)
c
901   FORMAT(T3,'@A2get_gridnmlist, *PES_GRID namelist not found or', 
     &       ' incomplete.')
      GO TO 999 
 999  CONTINUE
C
      CLOSE(UNIT=4,STATUS='KEEP')
      CALL ERREX
C
  99  CONTINUE
      CLOSE(UNIT=4,STATUS='KEEP')
      
      RETURN
      END
