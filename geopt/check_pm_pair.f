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
C
      LOGICAL FUNCTION CHECK_PM_PAIR(IBGN, IMD1, IMD2, IEND, JEND, 
     &                             SAME_PLANE, NATOMS)
      INTEGER NATOMS
      LOGICAL SAME_PLANE(NATOMS), SM_PLANE_BM, SM_PLANE_END
     
CSSS      WRITE(6,*) IBGN, IMD1, IMD2, IEND, JEND

      IF (IEND .EQ. IBGN .OR. IEND .EQ. IMD1 .OR. IEND. EQ. IMD2)
     &    THEN
          CHECK_PM_PAIR = .FALSE. 
          RETURN
      ENDIF
      IF (JEND .EQ. IBGN .OR. JEND .EQ. IMD1 .OR. JEND. EQ. IMD2)
     &    THEN 
          CHECK_PM_PAIR = .FALSE. 
          RETURN
      ENDIF
C
      IF (SAME_PLANE(IBGN) .AND. SAME_PLANE(IMD1) .AND. 
     &    SAME_PLANE(IMD2)) SM_PLANE_BM = .TRUE.
      
      SM_PLANE_END = SAME_PLANE(IEND) .AND.  SAME_PLANE(JEND)
CSSS      WRITE(6,*) SM_PLANE_BM, SM_PLANE_END
     
      IF ((SM_PLANE_BM .AND. SM_PLANE_END) .OR. (.NOT. 
     &     SM_PLANE_BM .AND. .NOT. SM_PLANE_END)) THEN
          CHECK_PM_PAIR = .FALSE.
          RETURN
      ENDIF
      IF ((SM_PLANE_BM .AND. .NOT. SM_PLANE_END) .OR. (.NOT. 
     &    SM_PLANE_BM .AND. SM_PLANE_END)) THEN
          CHECK_PM_PAIR = .TRUE.
          RETURN
      ENDIF
     
      END

