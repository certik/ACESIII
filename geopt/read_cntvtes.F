C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      SUBROUTINE READ_CNTVTES(IBNDTO, IUNIT, NATOMS)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL FLAG
C
#include "mxatms.par"
C
      DIMENSION IBNDTO(NATOMS, NATOMS), NCONPRCNTR(MXATMS),
     &          NCON(MAXCNTVS)

      DATA IONE, IJZERO /1, 0/
   
C Read the connectivites from the file CONNECT if present and
C bypass the automatic genration done in assign_cntvtes.F. 
C This is usefull when dealing with fused rings consists 
C of many atoms.
C
      CALL IZERO(NCON, MAXCNTVS)
      DO IATOMS = 1, NATOMS
         READ(IUNIT, 99) NCONPRCNTR(IATOMS), (NCON(J), J=1, 
     &                   NCONPRCNTR(IATOMS))
         NO_CNTS_PER_CENTER = NCONPRCNTR(IATOMS)
C
         IF (NO_CNTS_PER_CENTER .GE. MAXCNTVS) THEN
            WRITE(6,1000)
            CALL ERREX
         ENDIF
#ifdef _DEBUG_LVL0
         Write(6, "(a,I2)") "The # of connectivites: ",
     &                       NCONPRCNTR(IATOMS)
         Write(6,*)
         Write(6,"(a,10I3)") "The NCON array: ",
     &                     (NCON(J), J=1, NCONPRCNTR(IATOMS))
#endif
      ENDDO
   99 FORMAT (I2, 10(1X,I2))
C
C Create the IBNDTO array. This dicatate everything else that
C proceeds.
C
      DO IATMS = 1, NATOMS
         DO IBNDS = 1, NCONPRCNTR(IATMS)
            IF (NCON(IBNDS) .NE. IJZERO) 
     &          IBNDTO(NCON(IBNDS), IATMS) = IONE
         ENDDO
      ENDDO
C
C Some internal checks. Caution: needs to verify them.
C
      DO IATMS = 1, NATOMS
         DO IBNDS = 1, NCONPRCNTR(IATMS)
            KBNDS = IBNDTO(IBNDS, IATMS)
            DO JBNDS = 1,  NCONPRCNTR(KBNDS)
C
               IF (IBNDTO(JBNDS, KBNDS) .EQ. IATMS) THEN
                   FLAG = .TRUE.
                   WRITE(6, 2000)
                   CALL ERREX
               ENDIF
C
           ENDDO
        ENDDO
      ENDDO
C
CSSS        IF (.NOT. FLAG) THEN
CSSS            NCONPRCNTR(KBNDS) = NCONPRCNTR(KBNDS) + 1
CSSS            IBNDTO(NCONPRCNTR(KBNDS), KBNDS) = IATMS
CSSS        ENDIF
C
 1000 FORMAT("@-Rd_cntvtes - Exceeds the maximum number of
     &        bonds allowed for a center")

 2000 FORMAT("@-Rd_cntvtes - Duplicate connectivity", " I=",I5,
     &      " IB=",I3," J=",I5," JB=",I3,'.')

       RETURN 
       END

