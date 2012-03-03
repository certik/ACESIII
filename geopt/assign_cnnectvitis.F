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
      SUBROUTINE ASIGN_CNTVTES(CARTCOORD, BNDLENTHS, IATOMICNMBER, 
     &                         SMOFCOVRADI, NCONPRCNTR, IBNDTO,
     &                         MAXCNTVS, NRATMS, MARK_FRAGMENTS, 
     &                         LENGTH_FRAGMENTS,FRAGSCR,
     &                         BLNGT_IJFRGS)
C
C We need to decide some sort of a connectivity pattern if want to do
C Cartesian optimizations based on RIC. I prefer that the user specify the 
C connectivities, but Gaussian crowd certainly going to cry!
C So Let's use this covalent radii scheme Dr. Fau is talking about. 
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION PNTSEV
      LOGICAL FLAG, NOPClFPS, FRAGMENT, CONN_FILE_EXISTS, REGULAR_BND
      INTEGER AATOM, BATOM, VALENCY, iND, IDUM,
     &        JDUM, CONNECTED_ATOMS
C
C Watson
C
      INTEGER  FRAGSCR (1:NRATMS*NRATMS)
      INTEGER  IPRINT,IERR,MAXFRAG,NFRAG_TMP
C
C Watson
C
CSSS    NEIGHBORS(NRATMS, NRATMS)
C
      DIMENSION CARTCOORD(3*NRATMS), BNDLENTHS(NRATMS, NRATMS), 
     &          IATOMICNMBER(NRATMS), SMOFCOVRADI(NRATMS, NRATMS),
     &          IBNDTO(NRATMS, NRATMS), NCONPRCNTR(NRATMS),
     &          VALENCY(86), MARK_FRAGMENTS(NRATMS, NRATMS), 
     &          LENGTH_FRAGMENTS(NRATMS),BLNGT_IJFRGS(NRATMS)
C      
      DIMENSION VRADII(17)
      DATA VRADII /1.20D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.7D0,
     &             1.55D0, 1.52D0, 1.35D0, 0.0D0, 0.0D0, 0.0D0,
     &             0.0D0, 0.0D0, 1.90D0, 1.85D0, 1.80D0/
      DATA VALENCY /1, 0, 1, 2, 3, 4, 4, 2, 1, 0, 1, 2, 6, 6, 5,
     &              6, 1, 0, 1, 2, 6, 6, 6, 8, 6, 6, 6, 6, 6, 3,
     &              4, 3, 2, 1, 0, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6,
     &              6, 6, 3, 4, 3, 2, 1, 0, 1, 2, 12, 6, 6, 6, 
     &              6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     &              6, 6, 6, 6, 6, 3, 4, 3, 2, 1, 0, 1/
C
      DATA PNTSEV /0.71D0/, ZZERO /0.D0/, HALF/0.50D0/
      DATA IONE /1/, JZERO /0/
      DATA ZIONE /1.D0/
      DATA AMTOBOHR /1.8897265/
C
C
C Calculate internuclear distances; The calculated distances are 
C in Angstroms Call Zero sets up a 0 matrix for bond lenghts
C
      CALL ZERO(BNDLENTHS, NRATMS*NRATMS)
      CALL IZERO(IBNDTO, NRATMS*NRATMS)
C      
      INDEX1 = 1
      INDEX2 = 4
      IUNIT  = 199
C
      INQUIRE(FILE='CONNECT',EXIST=CONN_FILE_EXISTS)
      IF (CONN_FILE_EXISTS) THEN
          OPEN(UNIT=IUNIT,FILE='CONNECT',FORM='FORMATTED',
     &         STATUS='OLD')
#ifdef _DEBUG_LVL0
          Write(6,*) "Entering Read Connectivities"
#endif 
          CALL READ_CNTVTES(IBNDTO, IUNIT, NRATMS)
          CLOSE (IUNIT)
          RETURN
      ENDIF
C
      DO 10 IBNDS = 2, NRATMS
         DO 20 JBNDS = 1, IBNDS - 1
            BNDLENTHS(JBNDS, IBNDS) = DIST(CARTCOORD(INDEX1),
     &                                CARTCOORD(INDEX2))
            BNDLENTHS(IBNDS, JBNDS) = BNDLENTHS(JBNDS, IBNDS)
            INDEX1 = INDEX1 + 3 
 20      CONTINUE
         INDEX2 = INDEX2 + 3 
         INDEX1 = IONE
 10   CONTINUE
C
C Now calculate the scaled sum of the covalent radii for each unique
C pair of atoms
C     
cSSS CALL GETREC(20, 'JOBARC', 'ATOMCHRG', NRATMS, IATOMICNMBER)
C      
      CALL GETCOVLRADI(IATOMICNMBER, SMOFCOVRADI, NRATMS)
C      
C Now assign preliminary connectivites based on simple fact that
C length between two centers is less than or equal to 1.3 times sum of the 
C covalent radii of the centers. If they are connected, assign integer
C 1 for that pair, if not leave it as 0.
C This simple criteria is used by Helgaker and many others to assign 
C connectivities. 
C
      DO IBNDS = 2, NRATMS
         DO JBNDS = 1, (IBNDS - 1)
            IF (IBNDTO(JBNDS,IBNDS).EQ.JZERO.AND.
     &          BNDLENTHS(JBNDS,IBNDS).LE.SMOFCOVRADI(JBNDS,IBNDS)
     &         ) THEN
                IBNDTO(JBNDS,IBNDS) = IONE
#ifdef _DEBUG_LVLM1
            Print*, "Bonded atoms, bond-length, sum of c-radii"
            Print*,jbnds, ibnds, BNDLENTHS(JBNDS, IBNDS)*0.529177249d0, 
     &             SMOFCOVRADI(JBNDS,IBNDS)*0.529177249d0, 
     &             SMOFCOVRADI(IBNDS,JBNDS)*0.529177249d0
#endif
            END IF
         END DO
      END DO
C
C These condition pretty much take care of regular bonds. Now we
C need to consider isolated fragments and connect them.
C Identifying isolated fragments is a graph theory problem and I don't
C claim to be an expert on graph theory. There are efficient algorithms 
C and specialized programs to do this. The following is very simple
C and hopefully will work 90% of the time without heavy penalty. 
C 
C First, count the total number of connections for each atom and
C store the index of the other atoms that are connected to the
C given atom (in NCONPRCNTR, ICONNECTO arrays respectively). 
C Reorder these two arrays in descending order with respect to the
C number of connections.  
C Take a pair of atoms (starting with maximum connectivities) 
C and check whether they are connected. If they are not
C connected then search outwards (at predefined maxlevels
C representing a maximum possible chain length, for our
C cases I am thinking of 5). The searching is done as follows;
C Take the number of connections for the chosen pair and test whether at
C least one of them is connected; if so, then the choosen pair is not
C part of a disconnected fragment. If all of the atoms are not
C connected, then we look at the next level. When we consider all pairs
C and the FRAGMENT is still true then we have to mark that pair as a
C part of an isolated fragment. 
C 
      CALL IZERO(NCONPRCNTR, NRATMS)
      DO IBNDS = 1, NRATMS
         DO JBNDS = 1, NRATMS
            IF (IBNDTO(JBNDS, IBNDS) .NE. JZERO .OR.
     &          IBNDTO(IBNDS, JBNDS) .NE. JZERO) THEN
                NCONPRCNTR(IBNDS) = NCONPRCNTR(IBNDS) + 1
                IF (NCONPRCNTR(IBNDS) .GT. MAXCNTVS) THEN
                    WRITE(6, 9) MAXCNTVS
                    CALL ERREX
                ENDIF
            ENDIF 
         ENDDO 
      ENDDO
  9   FORMAT(T4, 'The maximum no. of connectvities allowed for',
     &       'for a center ', I2, ' exceeded')
C
#ifdef _DEBUG_LVLM1
      Print*, "The number of connectivities per center"
      Print*,  (NCONPRCNTR(I), I=1, NRATMS)
#endif 
C
C Lets identify hydrogen bonds. H-bonds are assigned if H-X distance
C is within the sum of the covalent radii of H and X(=N, O, P, Cl, F,
C P, S) and 0.9 times the Van der Walls radii and Y...H-X angle
C is greater than 90 degrees.
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Print*, "The connectivity assignments before H-bond block"
      Print*, ((ibndto(j, i), j=1, nratms), I=1, nratms)     
      Write(6,*)
#endif
C
      NHBONDS = 0
      CALL DSCAL(NRATMS*NRATMS, 1.0D0/1.3D0, SMOFCOVRADI, 1)
      DO IBNDS = 1, NRATMS
         IF (IATOMICNMBER(IBNDS) .EQ. 1) THEN
            DO JBNDS = 1, (IBNDS - 1)
               NOPClFPS = ((IATOMICNMBER(JBNDS).EQ. 7)  .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 8)  .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 9)  .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 10) .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 15) .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 16) .OR.
     &                    (IATOMICNMBER(JBNDS) .EQ. 17))

C
C The check for regular bond is required in order to avoid treating
C previoulsly assigned h-bond connections to assign new hydrogren
C bonds. 02/22, Ajith Perera.
C
               DIST_IJ =  BNDLENTHS(JBNDS, IBNDS)
               COVR_IJ = SMOFCOVRADI(JBNDS,IBNDS)*1.30D0
               IF (DIST_IJ .LE. COVR_IJ) THEN
                    REGULAR_BND = .TRUE.
               ELSE
                    REGULAR_BND = .FALSE.
               ENDIF
C
               IF (NOPClFPS.AND.IBNDTO(JBNDS,IBNDS).NE.JZERO) THEN
                  DO KBNDS = 1, NRATMS
                     IF (((IATOMICNMBER(KBNDS) .EQ. 7)  .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 8)  .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 9)  .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 10) .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 15) .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 16) .OR.
     &                   (IATOMICNMBER(KBNDS)  .EQ. 17)) .AND.
     &                   ((KBNDS .NE. IBNDS) .AND. (KBNDS .NE.
     &                    JBNDS))) THEN
                        IF (IBNDS .GT. KBNDS) THEN
                           DIST_HY     = BNDLENTHS(KBNDS, IBNDS)
                           SMOFCOVLRAD = SMOFCOVRADI(KBNDS,IBNDS)
                        ELSE
                           DIST_HY = BNDLENTHS(IBNDS, KBNDS)
                           SMOFCOVLRAD = SMOFCOVRADI(IBNDS,KBNDS)
                        END IF
                        CALL CMP_ANGLE(JBNDS, IBNDS, KBNDS,
     &                                 CARTCOORD, ANG_XHY, NRATMS)
                        KATOM = IATOMICNMBER(KBNDS)
                        IATOM = IATOMICNMBER(IBNDS)
                        SMOFVRADII = (VRADII(KATOM)+VRADII(IATOM))
     &                               *0.9D0
                        SMOFVRADII = AMTOBOHR*SMOFVRADII
#ifdef _DEBUG_LVLM1
                  Write(6,"(a)") "Critical parameters for H-bonds"
                  Write(6,"(4(1x,F10.5))") DIST_HY, ANG_XHY , 
     &                          SMOFVRADII, SMOFCOVLRAD
#endif
                        IF ((DIST_HY .GT. SMOFCOVLRAD).AND.
     &                      (DIST_HY .LT. SMOFVRADII).AND.
     &                      (ANG_XHY .GT. 90.0D0)) THEN
                           NHBONDS = NHBONDS + 1
                           IF (IBNDS .GT. KBNDS) THEN
                              IBNDTO(KBNDS, IBNDS) = IONE
                           ELSE
                              IBNDTO(IBNDS, KBNDS) = IONE
                           END IF
                        END IF
                     END IF
                  END DO
               END IF
            END DO
         END IF
      END DO
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a,I4)") "The no. of H bonds", NHBONDS
#endif
C
C#ifdef _NOSKIP
C
C
CSSS      do iatoms=2,nratms
CSSS         do jatoms=1,(iatoms-1)
CSSS            if (ibndto(jatoms,iatoms).eq.1) ibndto(iatoms,jatoms)=1
CSSS
CSSS      call identfy_fragments(nratms,ibndto,connected_atoms,nfrags,
CSSS     &                       length_fragments,mark_fragments)
CSSS          enddo
CSSS      enddo
C
C Watson
C
#ifdef _DEBUG_LVL0
      iprint = 50
#else
      iprint =0
#endif
      
      CALL  FRAG__INIT (nratms,ibndto,iprint,ierr,fragscr,
     +                  nfrags,length_fragments,mark_fragments)
C
C Watson
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,i3)") "The number of isolated fragments per center", 
     &  nfrags
      Write(6,*)
      Write(6,"(a)") "The length_fragment array"
      Write(6,*)
      Write(6,"(4I3)")  (length_fragments(I), I=1, NRATMS)
      Write(6,*)
      Write(6,"(a)") "Fragment atom lables"
      Write(6,*)
      Write(6,"(4I3)") ((mark_fragments(I, J), I=1, nratms),
     &         J=1, NRATMS)
#endif
C
C Assign the primary interfragment connectivity
C
      IF (NFRAGS .GT. 1) THEN
         INDEX = 1
         DO IATMS = 1, NFRAGS
            DO JATMS = NFRAGS, IATMS+1, -1
               ARCH_DIST = 20.0D0
             DIST_IJFRAG =  0.0D0
               DO IFRAGS = 1, LENGTH_FRAGMENTS(IATMS)
                  DO JFRAGS = 1, LENGTH_FRAGMENTS(JATMS)
                     IFRAG_CENT  = MARK_FRAGMENTS(IATMS, IFRAGS)
                     JFRAG_CENT  = MARK_FRAGMENTS(JATMS, JFRAGS)
C
                     DIST_IJFRAG = BNDLENTHS(IFRAG_CENT,
     &                                       JFRAG_CENT)
C               
#ifdef _DEBUG_LVL0
                     Write(6,"(4I4)") IATMS, JATMS, IFRAG_CENT, 
     &                                JFRAG_CENT
                     Write(6,"(2F10.5)") DIST_IJFRAG, ARCH_DIST
#endif

                     IF (DIST_IJFRAG .LE. ARCH_DIST) THEN
                         ARCH_DIST = DIST_IJFRAG
                         I_ARCH    = IFRAG_CENT
                         J_ARCH    = JFRAG_CENT
                     ENDIF
                  ENDDO
               ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a,F10.5)") "The smallest distance between ",
     &                       "fragments: ", ARCH_DIST
      Write(6,"(a,2I3)") "The atom indices of corrs. to smalest: ", 
     &         I_ARCH, J_ARCH
#endif
C
C In order to avoid connecting fragments that are way too far apart
C let's have a minimum reasonble threshold to the distance between the 
C fragment's atoms to a assign a connectivity. The choice for this
C distance must be reevlauted as we learn more about this. At the moment
C the threshold is 2.5 sum of the covalent radii of the fragment atoms.
C This supposedly connect the nearest neighbors. Ajith Perera, 08/2008 
C 
                   BLNGT_IJFRGS(INDEX) = ARCH_DIST
                                 INDEX = INDEX + 1
C
             THRESHOLD_4NEARST = 2.50D0*SMOFCOVRADI(I_ARCH, J_ARCH)
             IF (ARCH_DIST .LE. THRESHOLD_4NEARST) THEN

                IBNDTO(J_ARCH, I_ARCH) = IONE
                IBNDTO(I_ARCH, J_ARCH) = IONE

             ENDIF
C
             ENDDO
         ENDDO
      ENDIF
C
C Assign auxilary interfragment connectivities
C
C#ifdef _NOSKIP
C
      IF (NFRAGS .GT. 1) THEN
         INDEX = 1
         DO IATMS = 1, NFRAGS
            DO JATMS = NFRAGS, IATMS+1, -1
               ARCH_DIST = BLNGT_IJFRGS(INDEX)
               DO IFRAGS = 1, LENGTH_FRAGMENTS(IATMS)
                  DO JFRAGS = 1, LENGTH_FRAGMENTS(JATMS)
                     IFRAG_CENT  = MARK_FRAGMENTS(IATMS, IFRAGS)
                     JFRAG_CENT  = MARK_FRAGMENTS(JATMS, JFRAGS)
C
                     DIST_IJFRAG  = BNDLENTHS(IFRAG_CENT,
     &                                        JFRAG_CENT)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a,F10.5)") "The smallest distance between ",
     &                       "fragments: ", ARCH_DIST
      Write(6,"(a,a,2I3,F10.5)") "The frag. atom indices and the",
     &          " distance: ", IFRAG_CENT, JFRAG_CENT, DIST_IJFRAG
#endif

                     IF (DIST_IJFRAG .LT. 1.30D0*ARCH_DIST
     &                   .OR. DIST_IJFRAG .LT. 2.0D0*
     &                   AMTOBOHR) THEN
                                      I_ARCH    = IFRAG_CENT
                                      J_ARCH    = JFRAG_CENT
                         IBNDTO(J_ARCH, I_ARCH) = IONE
                         IBNDTO(I_ARCH, J_ARCH) = IONE
                     ENDIF
                   ENDDO
                ENDDO
                INDEX = INDEX + 1
             ENDDO
         ENDDO
      ENDIF
C
C#endif  
C#endif
C
C 
C Symmetrize the connectivity array.
C
cSSS      DO 50 IATMS = 1, NRATMS
cSSS         DO 60 IBNDS = 1, NCONPRCNTR(IATMS)
cSSS            KBNDS = IBNDTO(IBNDS, IATMS)
cSSS            FLAG = .FALSE.
cSSS            DO 70 JBNDS = 1, NCONPRCNTR(KBNDS)
cSSS
cSSS                  IF (IBNDTO(JBNDS, KBNDS) .EQ. IATMS) THEN
cSSS                  FLAG = .TRUE.
cSSS                  CALL WRITE(LUOUT, 2000)
cSSS                  CALL ERREX
cSSS               ENDIF
cSSS
cSSS 70         CONTINUE
cSSS 60      CONTINUE
cSSS     
cSSS         IF (.NOT. FLAG) THEN
cSSS            NCONPRCNTR(KBNDS) = NCONPRCNTR(KBNDS) + 1
cSSS            IBNDTO(NCONPRCNTR(KBNDS), KBNDS) = IATMS
cSSS         ENDIF
cSSS
cSSS 50   CONTINUE
cSSS
cSSS 1000 FORMAT("@-assign_connectvities - Exceeds the maximum number of 
cSSS     &        bonds allowed for a center")
cSSS
cSSS 2000 FORMAT("@-assign_connectvitis - Duplicate connectivity", " I=",I5,
cSSS              " IB=",I3," J=",I5," JB=",I3,'.')

      RETURN
      END
 
