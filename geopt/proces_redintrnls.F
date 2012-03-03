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
      SUBROUTINE PROCESS_REDINTRNLS(Q, REDUNCO, IREDUNCO, IBNDTO, 
     &                              I_UNIQUE, PTGRP, TOTNOFBND, 
     &                              TOTNOFANG, TOTNOFDIH, TOTREDNCO,
     &                              NEQ, IEQUIV, NATOMS, IPLUS_MINUS,
     &                              IPRNT, EPSILON)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
c      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
c      PARAMETER (MXATMS=96, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 900)
#include "mxatms.par"
C
      INTEGER TOTREDNCO, ORDER_PTGRP, TOTNOFBND, TOTNOFANG, 
     &        TOTNOFDIH
      LOGICAL SYMEQVLNT, SAME_BGN, I_UNIQUE(TOTREDNCO), SAME_END,
     &        SAME_MID, EQEV_MID, EQEV_BRDGE, SAME_TWO, ABELIAN,
     &        SAMEDIHS, SAMEBNDS, SAMEANGLS, SAME_PLANE(MXATMS),
     &        PM_PAIR, CHECK_PM_PAIR, PRV_USED
      CHARACTER*4 PTGRP

      DIMENSION IREDUNCO(4, TOTREDNCO), IBNDTO(NATOMS, NATOMS),
     &          Q(3*NATOMS), NEQ(MAXREDUNCO), REDUNCO(TOTREDNCO), 
     &          IEQUIV(3*MXATMS, 3*MXATMS), 
     &          IUNIQUE_OFFSET(MAXREDUNCO), IPLUS_MINUS(MAXREDUNCO),
     &          I_PAIR_FLAG(MAXREDUNCO)
C
      CALL IZERO(IBNDTO, NATOMS*NATOMS) 
      CALL IZERO(IPLUS_MINUS, MAXREDUNCO)
      ABELIAN = .FALSE. 

      IF (PTGRP .EQ. "C2v ") THEN
         ABELIAN     = .TRUE. 
         ORDER_PTGRP = 4
      ELSE IF (PTGRP .EQ. "D2h ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 8 
      ELSE IF (PTGRP .EQ. "D2  ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 4
      ELSE IF (PTGRP .EQ. "C2  ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 2
      ELSE IF (PTGRP .EQ. "C2h ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 4
      ELSE IF (PTGRP .EQ. "C s ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 2 
      ELSE IF (PTGRP .EQ. "C i ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 2 
      ELSE IF (PTGRP .EQ. "C1  ") THEN
         ABELIAN     = .TRUE.
         ORDER_PTGRP = 1  
      ENDIF 
C
C Let's get the permutation of atoms with highest order symmetry 
C operation, and the order of the highest rotation axis.
C
      IF (.NOT. ABELIAN) THEN
         CALL IGETREC(20, "JOBARC", "SYMEQUIV", NATOMS, NEQ)
      ENDIF 
C
CSSS      WRITE(6,*) "The Permu. Array"
CSSS      WRITE(6,*) (NEQ(I), I= 1, NATOMS)
C
      IF (ABELIAN) THEN
         DO IATOM = 1, NATOMS  
            DO JATOM = 1, NATOMS
               IF (IATOM .NE. JATOM) THEN
                  CALL IDNTFY_SYMEQV_ATMS(PTGRP, ORDER_PTGRP, 
     &                                    Q((IATOM - 1)*3 + 1), 
     &                                    Q((JATOM - 1)*3 + 1),
     &                                    SYMEQVLNT)
                  IF (SYMEQVLNT) THEN 
                     IBNDTO(JATOM, IATOM) = 1
                  ENDIF 
               ENDIF
            ENDDO
         ENDDO
C
      ELSE
C
         DO IATOM = 1, NATOMS
C
            KATOM1 = NEQ(IATOM)
C
            DO JATOM = 1, IATOM - 1

               KATOM2 = NEQ(JATOM)

               IBNDTO(IATOM, KATOM1) = 1
               IBNDTO(KATOM1, IATOM) = 1
               IBNDTO(JATOM, KATOM2) = 1
               IBNDTO(KATOM2, JATOM) = 1
C
               IF (IATOM .EQ. KATOM2) THEN
                  IBNDTO(JATOM, KATOM1) = 1
                  IBNDTO(KATOM1, JATOM) = 1
               ENDIF
               IF (JATOM .EQ. KATOM1) THEN
                  IBNDTO(KATOM2, IATOM) = 1
                  IBNDTO(IATOM, KATOM2) = 1
               ENDIF
               IF (IBNDTO(JATOM, KATOM1) .EQ. 1) THEN
                  IBNDTO(KATOM2, IATOM) = 1
                  IBNDTO(IATOM, KATOM2) = 1
               ENDIF
               IF (IBNDTO(JATOM, KATOM2) .EQ. 1) THEN
                   IBNDTO(IATOM, KATOM1) = 1
                   IBNDTO(KATOM1, IATOM) = 1
               ENDIF
               IF (IBNDTO(IATOM, KATOM2) .EQ. 1) THEN
                  IBNDTO(JATOM, KATOM1) = 1
                  IBNDTO(KATOM1, JATOM) = 1
               ENDIF
               IF (IBNDTO(IATOM, KATOM1) .EQ. 1) THEN
                   IBNDTO(JATOM, KATOM2) = 1
                   IBNDTO(KATOM2, JATOM) = 1
               ENDIF
C
            ENDDO 
         ENDDO
C
      ENDIF
C
C Before we identify symmetry equivalent bonds, angles and dihedral
C angles, we need to some way to identify the atoms that can make 
C plus-minus dihedral pairs. At this point the logic is built on
C the condition that dihedral pairs arise when two symmetry identical
C atoms below or above (left or right) of a plane or axis of symmetry
C and three of the atoms that make the dihedral is on the plane or
C axis of symmetry. 
C
      IF (NATOMS .GE. 4 .AND. .NOT. PTGRP .EQ. "C1  " ) THEN
         CALL IGETREC(20, "JOBARC", "SAMEPLNE", NATOMS, NEQ)
cSSS         WRITE(6,*) (NEQ(I), I=1,NATOMS)
         DO IATOM = 1, NATOMS
            IF (NEQ(IATOM) .EQ. IATOM) SAME_PLANE(IATOM) = .TRUE. 
         ENDDO
      ENDIF
C
CSSS      WRITE(6,*) "SAME PLANE ARRAY"
CSSS      WRITE(6,*) (SAME_PLANE(I), I=1, NATOMS)
C
C Identify the symmetry unique RIC bond coordinates.
C
CSSS      WRITE(6,*) "Symmetry Equivalent Atoms"
CSSS      WRITE(6,*) IBNDTO
cSSS      WRITE(6,*) TOTNOFBND, TOTNOFANG, TOTNOFDIH
      CALL IZERO(NEQ, MAXREDUNCO)
      IUNIQUE = 0 
      DO IBONDS = 1, TOTNOFBND
         IBGN = IREDUNCO(1, IBONDS)
         IEND = IREDUNCO(2, IBONDS) 
         DO JBONDS = 1, IBONDS - 1
            SAMEBNDS = (DABS(REDUNCO(IBONDS) - REDUNCO(JBONDS))
     &                  .LE. EPSILON)
            IF (SAMEBNDS) THEN
               IF (IBONDS .NE. JBONDS) THEN
                  SAME_BGN = .FALSE.
                  SAME_END = .FALSE. 
                  SAME_TWO = .FALSE.
                  JBGN = IREDUNCO(1, JBONDS)
                  JEND = IREDUNCO(2, JBONDS)
cSSS                  WRITE(6,*) IBONDS, JBONDS 
cSSS                  IF (IPRNT .GT. 10) WRITE(6,*) IBGN, IEND, JBGN, 
cSSS     &                               JEND, IUNIQUE
                  IF (IBGN .EQ. JBGN) SAME_BGN = .TRUE. 
                  IF ((IBNDTO(JEND, IEND) .NE. 0) .AND. (SAME_BGN))
     &               THEN 
                     IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                        IOFFSET = IUNIQUE_OFFSET(JBONDS)
                     ENDIF
                     NEQ(IOFFSET) = NEQ(IOFFSET) + 1 
                     IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
cSSS            WRITE(6,*) IBONDS, JBONDS, NEQ(JBONDS)
cSSS            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                     GO TO 10 
                  END IF
C
                  IF (IEND .EQ. JEND) SAME_END = .TRUE.
                  IF ((IBNDTO(JBGN, IBGN) .NE. 0) .AND. (SAME_END))
     &               THEN 
                     IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                        IOFFSET = IUNIQUE_OFFSET(JBONDS)
                     ENDIF
                     NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                     IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
cSSS            WRITE(6,*) IBONDS, JBONDS, NEQ(JBONDS)
cSSS            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                     GO TO 10
                  ENDIF
C
                  IF (IBGN .EQ. JEND) SAME_TWO = .TRUE.
                  IF ((IBNDTO(JBGN, IEND) .NE. 0) .AND. (SAME_TWO))
     &            THEN
                     IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                        IOFFSET = IUNIQUE_OFFSET(JBONDS)
                     ENDIF
                     NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                     IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
cSSS            WRITE(6,*) IBONDS, JBONDS, NEQ(IUNIQUE)
cSSS            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                     GO TO 10
                  ENDIF
C
                  IF (IEND .EQ. JBGN) SAME_TWO = .TRUE.
                  IF ((IBNDTO(JEND, IBGN) .NE. 0) .AND. (SAME_TWO))
     &            THEN
                     IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                        IOFFSET = IUNIQUE_OFFSET(JBONDS)
                     ENDIF
                     NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                     IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
cSSS            WRITE(6,*) "I am Here"
cSSS            WRITE(6,*) IBONDS, JBONDS, NEQ(IUNIQUE)
cSSS            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                     GO TO 10
                  ENDIF
C
                  IF ((IBNDTO(IBGN, JBGN) .NE. 0) .AND.
     &               (IBNDTO(IEND, JEND) .NE. 0) .OR.
     &               (IBNDTO(IBGN, JEND) .NE. 0) .AND.  
     &               (IBNDTO(IEND, JBGN) .NE. 0))
     &            THEN
                     IF (IUNIQUE_OFFSET(JBONDS) .NE. 0) THEN
                        IOFFSET = IUNIQUE_OFFSET(JBONDS)
                     ENDIF
                     NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                     IEQUIV(IOFFSET, NEQ(IOFFSET)) = IBONDS
cSSS            WRITE(6,*) IBONDS, JBONDS, NEQ(JBONDS)
cSSS            WRITE(6,*) IEQUIV(JBONDS, NEQ(JBONDS))
C
                     GO TO 10
                  ENDIF
C 
               ENDIF
            ENDIF
C
         END DO
C
         I_UNIQUE(IBONDS) = .TRUE. 
         IUNIQUE = IUNIQUE + 1
         IUNIQUE_OFFSET(IBONDS) = IUNIQUE
C
 10      CONTINUE
C
      END DO
C
C Identify the symmetry unique RIC angle coordinates. Also, this block 
C of code mark the bond angles that are constrained by symmetry. For
C exmple 60 degere bond angles in BF3 are symmetry constrained. 
C 
cSSS      WRITE(6,*) (IUNIQUE_OFFSET(J), J = 1, TOTNOFBND)
      DO IANGL = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
         IBGN = IREDUNCO(1, IANGL)
         IMID = IREDUNCO(2, IANGL)
         IEND = IREDUNCO(3, IANGL)
         DO JANGL = (TOTNOFBND + 1), IANGL - 1 
            SAMEANGLS = (DABS(REDUNCO(IANGL) - REDUNCO(JANGL))
     &                  .LE. EPSILON)
            IF (SAMEANGLS) THEN
cSSS               WRITE(6,*) IANGL, JANGL
               SAME_MID = .FALSE.
               EQEV_MID = .FALSE. 
               JBGN = IREDUNCO(1, JANGL)
               JMID = IREDUNCO(2, JANGL)
               JEND = IREDUNCO(3, JANGL)
               IF (IANGL .NE. JANGL) THEN
                  IF (IMID .EQ. JMID)  SAME_MID = .TRUE.
                     IF (SAME_MID) THEN 
cSSS                     WRITE(6,*) "CHECKING SYM EQV ANGLES"
cSSS                     WRITE(6,*) IBGN, JBGN,IEND, JEND
cSSS                     WRITE(6,*) IBNDTO(IBGN, JBGN), IBNDTO(IEND, JEND),
cSSS     &                          IBNDTO(IEND, JBGN), IBNDTO(IBGN, JEND)
                         IF (IBGN .EQ. JBGN .AND. 
     &                      IBNDTO(IEND, JEND) .NE. 0) THEN 
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*)IANGL,JANGL,NEQ(IOFFSET),IOFFSET
                            GO TO 20
                         ENDIF
                         IF (IEND .EQ. JEND .AND. 
     &                      IBNDTO(IBGN, JBGN) .NE. 0) THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*)IANGL,JANGL,NEQ(IOFFSET),IOFFSET
                            GO TO 20
                         ENDIF
                         IF (IBGN .EQ. JEND .AND. 
     &                      IBNDTO(IEND, JBGN) .NE. 0) THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*)IANGL,JANGL,NEQ(IOFFSET),IOFFSET
                            GO TO 20
                         ENDIF
                         IF (IEND .EQ. JBGN .AND. 
     &                      IBNDTO(IBGN, JEND) .NE. 0) THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*) IANGL,JANGL, NEQ(IOFFSET),IOFFSET
                            GO TO 20
                         ENDIF
                         IF ((IBNDTO(IBGN, JBGN) .NE. 0) .AND. 
     &                      (IBNDTO(IEND, JEND) .NE. 0)  .OR. 
     &                      (IBNDTO(IBGN, JEND) .NE. 0)  .AND.
     &                      (IBNDTO(IEND, JBGN) .NE. 0))
     &                   THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*) IANGL,JANGL,NEQ(IOFFSET),IOFFSET
                            GO TO 20
                         ENDIF
                     ENDIF
C
                     IF (IBNDTO(IMID, JMID) .NE. 0) EQEV_MID = .TRUE. 
                     IF (EQEV_MID) THEN
                        IF (IBGN .EQ. JBGN .AND.
     &                      IBNDTO(IEND, JEND) .NE. 0) THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*)IANGL,JANGL,NEQ(IOFFSET),IOFFSET
                           GO TO 20
                        ENDIF
                        IF (IBGN .EQ. JEND .AND.
     &                      IBNDTO(IEND, JBGN) .NE. 0) THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*)IANGL,JANGL,NEQ(IOFFSET),IOFFSET
                            GO TO 20
                        ENDIF
                        IF (IEND .EQ. JBGN .AND.
     &                      IBNDTO(IBGN, JEND) .NE. 0) THEN
                            IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                               IOFFSET = IUNIQUE_OFFSET(JANGL)
                            ENDIF
                            NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                            IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                            WRITE(6,*) IANGL,JANGL, NEQ(IOFFSET),IOFFSET
                            GO TO 20
                         ENDIF
                         IF ((IBNDTO(IBGN, JBGN) .NE. 0) .AND.
     &                     (IBNDTO(IEND, JEND) .NE. 0)  .OR. 
     &                     (IBNDTO(IBGN, JEND) .NE. 0)  .AND.
     &                     (IBNDTO(IEND, JBGN) .NE. 0)) 
     &                   THEN
                           IF (IUNIQUE_OFFSET(JANGL) .NE. 0) THEN
                              IOFFSET = IUNIQUE_OFFSET(JANGL)
                           ENDIF
                           NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                           IEQUIV(IOFFSET, NEQ(IOFFSET)) = IANGL
cSSS                           WRITE(6,*) "WHY AM I HERE", IANGL, JANGL,
cSSS     &                     NEQ(IOFFSET), IOFFSET
                           GO TO 20
                     ENDIF
C
                  ENDIF
               ENDIF
            ENDIF
C
         ENDDO
C 
         I_UNIQUE(IANGL) = .TRUE.
         IUNIQUE = IUNIQUE + 1
         IUNIQUE_OFFSET(IANGL) = IUNIQUE
cSSS         WRITE(6,*) I_UNIQUE(IANGL), IUNIQUE_OFFSET(IANGL), IANGL
C
 20      CONTINUE
C
      ENDDO 
C
C Identify the symmetry unique RIC dihedral angle coordinates. Also,
C this block of code mark the dihedrals angles that are constrained 
C by symmetry. For exmple 180 dihedral angles and dihedral angle pairs.
C
      DO IDIHS = (TOTNOFANG + TOTNOFBND + 1),  TOTREDNCO
         IBGN = IREDUNCO(1, IDIHS)
         IMD1 = IREDUNCO(2, IDIHS)
         IMD2 = IREDUNCO(3, IDIHS)
         IEND = IREDUNCO(4, IDIHS)
         DO JDIHS = (TOTNOFANG + TOTNOFBND + 1), IDIHS - 1 
            SAMEDIHS = (DABS(REDUNCO(IDIHS) - REDUNCO(JDIHS)) 
     &                  .LE. EPSILON)
CSSS            WRITE(6,*) "THE SAME DIHEDRALS", SAMEDIHS
            IF (SAMEDIHS) THEN
                EQEV_BRDGE = .FALSE.  
                SAME_MID   = .FALSE. 
                JBGN = IREDUNCO(1, JDIHS)
                JMD1 = IREDUNCO(2, JDIHS)
                JMD2 = IREDUNCO(3, JDIHS)
                JEND = IREDUNCO(4, JDIHS)  
CSSS                WRITE(6,*) IPLUS_MINUS(JDIHS), I_PAIR_FLAG(IDIHS),
CSSS     &                     IDIHS, JDIHS
                PRV_USED = (IPLUS_MINUS(JDIHS) .NE. 0 .OR. 
     &                      I_PAIR_FLAG(JDIHS) .NE. 0)
C
                IF (IDIHS .NE. JDIHS .AND. .NOT. PRV_USED) THEN 
C                  
                   PM_PAIR  = .FALSE. 
                   IF ((IMD1 .EQ. JMD1 .AND. IMD2 .EQ. JMD2) .OR. 
     &                (IMD1 .EQ. JMD2 .AND. IMD2. EQ. JMD1)) 
     &                SAME_MID  = .TRUE.
C
                   IF ((IMD1 .EQ. JMD1 .AND. IBNDTO(IMD2, JMD2)
     &                 .NE. 0)) SAME_MID = .TRUE. 
C
                   IF ((IMD1 .EQ. JMD2 .AND. IBNDTO(IMD2, JMD1)    
     &                 .NE. 0)) SAME_MID = .TRUE. 
C
                   IF ((IMD2 .EQ. JMD1 .AND. IBNDTO(IMD1, JMD2)
     &                 .NE. 0)) SAME_MID = .TRUE.
C
                   IF ((IMD2 .EQ. JMD2 .AND. IBNDTO(IMD1, JMD1)
     &                 .NE. 0)) SAME_MID = .TRUE.
C
                   IF ((IBGN .EQ. JBGN .AND. IEND .EQ. JEND .AND. 
     &                  SAME_MID) .OR. (IBGN .EQ. JEND .AND. IEND
     &                  .EQ. JBGN .AND. SAME_MID)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
CSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
CSSS     &                           PM_PAIR, "I GOT HERE"
                      GO TO 30
                   END IF
C   
                   IF ((IBGN .EQ. JBGN .AND. SAME_MID .AND. 
     &                 IBNDTO(IEND, JEND) .NE. 0)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
CSSS                      WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
CSSS     &                           PM_PAIR
                      GO TO 30
                   END IF 
C
                   IF ((IEND .EQ. JEND .AND. SAME_MID .AND. 
     &                 IBNDTO(IBGN, JBGN) .NE. 0)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IEND, IMD1, IMD2,
     &                           IBGN, JBGN, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                      WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                           PM_PAIR
                      GO TO 30
                   END IF 

                   IF ((IBGN .EQ. JEND .AND. SAME_MID .AND.
     &                 IBNDTO(IEND, JBGN) .NE. 0)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JBGN, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                      WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                           PM_PAIR
                      GO TO 30
                   END IF

                   IF ((IEND .EQ. JBGN .AND. SAME_MID .AND.
     &                 IBNDTO(IBGN, JEND) .NE. 0)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IEND, IMD1, IMD2,
     &                           IBGN, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                            PM_PAIR
                       GO TO 30
                   END IF
                   IF ((IBGN .EQ. JBGN .AND. SAME_MID .AND.
     &                  IEND .EQ. JEND)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
CSSS                       WRITE(6,*)IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
CSSS     &                           PM_PAIR
                       GO TO 30
                   END IF
                   IF ((IBGN .EQ. JEND .AND. SAME_MID .AND.
     &                 IEND .EQ. JBGN)) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JBGN, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1  
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                            PM_PAIR
                       GO TO 30
                   END IF
                   IF ((SAME_MID .AND. IBNDTO(IBGN, JBGN)  
     &                  .NE. 0) .AND.
     &                  (IBNDTO(IEND, JEND) .NE. 0) .OR.
     &                  (IBNDTO(IBGN, JEND) .NE. 0) .AND.
     &                  (IBNDTO(IEND, JBGN) .NE. 0))
     &             THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JEND, SAME_PLANE, NATOMS)
                       IF (.NOT. (SAME_PLANE(IBGN) .AND. 
     &                            SAME_PLANE(IEND) .AND.  
     &                            SAME_PLANE(JBGN) .AND. 
     &                            SAME_PLANE(JEND)))
     &                            PM_PAIR = .TRUE.  

cSSS                       IF (.NOT. PM_PAIR) THEN 
cSSS                           PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1,
cSSS     &                               IMD2, IEND, JBGN, 
cSSS     &                               SAME_PLANE, NATOMS)
cSSS                       ENDIF
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
                       
cSSS                       WRITE(6,*)IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                           PM_PAIR,"I AM HERE", IBNDTO
cSSS     &                           (IEND,JEND), IBNDTO(IBGN,JBGN)
                       GO TO 30
                   END IF
C
                   IF ((IBNDTO(IMD1, JMD1) .NE. 0) .AND.
     &                (IBNDTO(IMD2, JMD2) .NE. 0)  .OR.
     &                (IBNDTO(IMD1, JMD2) .NE. 0)  .AND.
     &                (IBNDTO(IMD2, JMD1) .NE. 0)) 
     &                 EQEV_BRDGE = .TRUE. 
cSSS                   WRITE(6,*) "I see a EQUV BRIDGE", EQEV_BRDGE
C
                   IF (EQEV_BRDGE .AND. IEND .EQ. JEND .AND.
     &                 IBNDTO(IBGN, JBGN) .NE. 0) THEN
                       PM_PAIR = CHECK_PM_PAIR(IEND, IMD1, IMD2,
     &                           IBGN, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                            PM_PAIR
                       GO TO 30
                   ENDIF
C
                   IF (EQEV_BRDGE .AND. IBGN .EQ. JBGN .AND.
     &                 IBNDTO(IEND, JEND) .NE. 0) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                 PM_PAIR
                       GO TO 30
                   ENDIF
C   
                   IF (EQEV_BRDGE .AND. IBGN .EQ. JEND .AND.
     &                 IBNDTO(IEND, JBGN) .NE. 0) THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JBGN, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                      ENDIF
                      IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                      ENDIF
                      NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                      IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                      WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                PM_PAIR
                      GO TO 30
                   ENDIF
C
                   IF (EQEV_BRDGE .AND. IEND .EQ. JBGN .AND.
     &                 IBNDTO(IBGN, JEND) .NE. 0) THEN
                       PM_PAIR = CHECK_PM_PAIR(IEND, IMD1, IMD2,
     &                           IBGN, JEND, SAME_PLANE, NATOMS)
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                 PM_PAIR
                       GO TO 30
                   ENDIF
C
                   IF ((EQEV_BRDGE .AND. IBNDTO(IBGN, JBGN)
     &                  .NE. 0) .AND.
     &                  (IBNDTO(IEND, JEND) .NE. 0) .OR. 
     &                  (IBNDTO(IBGN, JEND) .NE. 0) .AND.
     &                  (IBNDTO(IEND, JBGN) .NE. 0)) 
     &             THEN
                       PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1, IMD2,
     &                           IEND, JEND, SAME_PLANE, NATOMS)
                       IF (.NOT. (SAME_PLANE(IBGN) .AND.
     &                            SAME_PLANE(IEND) .AND.
     &                            SAME_PLANE(JBGN) .AND.
     &                            SAME_PLANE(JEND)))
     &                            PM_PAIR = .TRUE.
cSSS                       IF (.NOT. PM_PAIR) THEN
cSSS                           PM_PAIR = CHECK_PM_PAIR(IBGN, IMD1,
cSSS     &                               IMD2, IEND, JBGN, 
cSSS     &                               SAME_PLANE, NATOMS)
cSSS                       ENDIF 
                       IF (IUNIQUE_OFFSET(JDIHS) .NE. 0) THEN
                          IOFFSET = IUNIQUE_OFFSET(JDIHS)
                       ENDIF 
                       IF (PM_PAIR) THEN
                          IPLUS_MINUS(JDIHS) = 1
                          I_PAIR_FLAG(IDIHS) = 1
                          GO TO 40
                       ENDIF
                       NEQ(IOFFSET) = NEQ(IOFFSET) + 1
                       IEQUIV(IOFFSET, NEQ(IOFFSET)) = IDIHS
cSSS                       WRITE(6,*) IDIHS,JDIHS,NEQ(IOFFSET),IOFFSET,
cSSS     &                 PM_PAIR
                       GO TO 30
                   END IF
C
                ENDIF
            ENDIF 
cSSS            I_UNIQUE(IDIHS) = .TRUE.
cSSS            IUNIQUE = IUNIQUE + 1
C
         ENDDO                
C
 40      CONTINUE
         I_UNIQUE(IDIHS) = .TRUE.
         IUNIQUE = IUNIQUE + 1
         IUNIQUE_OFFSET(IDIHS) = IUNIQUE
 30      CONTINUE
      ENDDO 
C
cSSS      WRITE(6,*) IUNIQUE_OFFSET(31), IUNIQUE_OFFSET(32),
cSSS     &           IUNIQUE_OFFSET(33), IUNIQUE_OFFSET(34)
cSSS      WRITE(6,*) I_UNIQUE
cSSS      WRITE(6,*) NEQ
cSSS      DO I = 1, IUNIQUE
cSSS         DO J =1, NEQ(I)
cSSS         WRITE(6,*) I, IEQUIV(I, J)
cSSS         ENDDO
cSSS      ENDDO
C
      RETURN
      END
