         SUBROUTINE  FRAG__WORK
     +
     +                   (NATOMS,
     +                    CONNECT,
     +                    IPRINT,IERR,
     +                    INITFRAG,
     +
     +                              NFRAG,
     +                              FRAGLEN,
     +                              FRAGMAT  )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : FRAG__WORK
C  MODULE      : FRAGMENT DETERMINE
C  MODULE-ID   : FRAG
C  SUBROUTINES : FRAG__LOOP_THROUGH
C
C  DESCRIPTION : This routine grabs an atom and then checks for
C                connections against it to fill in the fragment
C                array.
C
C                  Input:
C
C                    IPRINT       =  level of printing variable
C                    IERR         =  sets the return error
C                    NATOMS       =  total number of atoms
C                    CONNECT      =  connectivity array
C                    INITFRAG     =  array to hold previous
C                                    fragmentation pattern
C
C                  Output:
C
C                    NFRAG        =  total number of fragments
C                    FRAGLEN      =  contains the length of the
C                                    fragments in FRAGMAT
C                    FRAGMAT      =  contains the fragments
C
C
C  AUTHOR      : Thomas Watson Jr.
C------------------------------------------------------------------------
C
C
C             ...Declare variables and include files!
C
C
         IMPLICIT    NONE

         LOGICAL    INCFRAG,HOTEQNATOM
         LOGICAL    ISFRAG,ALREADY

         INTEGER    I,J,K,IHOT,NFRAG
         INTEGER    IFRAG,IATOM,ILOC
         INTEGER    ONE,BOND,HOT,NEW
         INTEGER    IPRINT,IERR,IFRAGLN
         INTEGER    NATOMS,LSCR1
         INTEGER    ZERO,ICOUNT

         INTEGER    FRAGLEN  (1:NATOMS)

         INTEGER    CONNECT  (1:NATOMS,1:NATOMS)
         INTEGER    INITFRAG (1:NATOMS,1:NATOMS)
         INTEGER    FRAGMAT  (1:NATOMS,1:NATOMS)

         PARAMETER  ( ZERO = 0 )
         PARAMETER  ( ONE  = 1 )
C
C
C------------------------------------------------------------------------
C
C
C             ...First we're going to see what is
C                bonded to what.  We still assume
C                there are NATOM fragments, but 
C                we can form this matrix and then
C                search it for overlapping elements.
C
C                This will also assign diagonal 
C                elements of the connectivity array
C                to the atom number, since each atom
C                is initially assumed to be its own 
C                fragment.
C
C                Also, print the array if user wants.
C
C
         DO IFRAG = 1,     NATOMS
         DO IATOM = IFRAG, NATOMS

            BOND = CONNECT (IFRAG,IATOM)
            IF ((BOND .EQ. ONE) .OR. (IFRAG .EQ. IATOM)) THEN
               INITFRAG (IFRAG,IATOM) = IATOM
               FRAGMAT  (IFRAG,IATOM) = IATOM
            ENDIF

         ENDDO
         ENDDO

         IF (IPRINT .GE. 10) THEN
            WRITE (*,*) ' @FRAG__WORK: Connectivity array - '
            DO I = 1, NATOMS
               WRITE (*,'(99I3)') (CONNECT (I,J), J = 1, NATOMS)
            ENDDO
            WRITE (*,*) ''
            WRITE (*,*) ' @FRAG__WORK: First bonding pass - '
            DO I = 1, NATOMS
               WRITE (*,'(99I3)') (INITFRAG (I,J), J = 1, NATOMS)
            ENDDO
            WRITE (*,*) ''
         ENDIF
C
C
C             ...Basically, what we're gonna do here
C                is start with IFRAG equal to 1 and
C                pick the first atom in that fragment.
C
C                Since we know what that atom is bonded
C                do, we know what other atoms are in
C                that fragment.  If we can find what 
C                those other atoms are bonded to, then
C                we can add them in to the fragment.o
C
C                When we've exhausted all other fragments,
C                we're done with IFRAG equal to 1. So
C                we'll increment IFRAG up 1. 
C
C
         NFRAG = NATOMS
         IFRAG = 1

         DO WHILE (NFRAG .NE. 0)

            DO IHOT = IFRAG, NATOMS
C
C
C          ...Pick an atom to search for connections against!
C             If it is equal to the last atom, note that so
C             we can do a second pass in the next routine.
C
C
               HOT = INITFRAG (IFRAG,IHOT)
               HOTEQNATOM = (IHOT .EQ. NATOMS)

               CALL  FRAG__LOOP_THROUGH
     +
     +                    ( NATOMS,IFRAG,
     +                      HOT,IHOT,
     +                      HOTEQNATOM,
     +                      INCFRAG,
     +                      IPRINT,IERR,
     +                      INITFRAG,
     +
     +                            FRAGMAT )
     +
            ENDDO

            IF (IPRINT .GE. 50) THEN
               WRITE (*,*) ' @FRAG_WORK: Debugging INITFRAG - '
               DO I = 1, NATOMS
                  WRITE (*,'(99I3)') (INITFRAG (I,J), J = 1, NATOMS)
               ENDDO
               WRITE (*,*) ''

               WRITE (*,*) ' @FRAG_WORK: Debugging FRAGMAT  - '
               DO I = 1, NATOMS
                  WRITE (*,'(99I3)') (FRAGMAT (I,J), J = 1, NATOMS)
               ENDDO
               WRITE (*,*) ''
            ENDIF
C
C
C          ...Increment the fragment if we have to, and
C             decrement the number of fragments!
C
C
            IF (INCFRAG) THEN
               NFRAG = NFRAG - 1
               IF (IFRAG .NE. NATOMS) IFRAG = IFRAG + 1
            ENDIF
C
C
C          ...Copy the fragmentation matrix into the test
C             matrix to go through another round of checking!
C
C
            DO I = 1, NATOMS
            DO J = 1, NATOMS
               INITFRAG (I,J) = FRAGMAT (I,J)
            ENDDO
            ENDDO

         ENDDO ! WHILE (NFRAG .NE. 0)
C
C
C             ...Now count the total number of fragments!
C                Also, reorder the arrays for ACES format!
C
C
         DO IFRAG = 1, NATOMS
         DO IATOM = 1, NATOMS
            INITFRAG (IFRAG,IATOM) = 0
         ENDDO
         ENDDO

         NFRAG  = 0
         ICOUNT = 0
         DO IFRAG = 1, NATOMS
            IFRAGLN = 0
            ALREADY = .FALSE.
            DO IATOM = 1, NATOMS
               
               ISFRAG = (FRAGMAT (IFRAG,IATOM) .GT. 0)
               IF (ISFRAG) THEN
                  IFRAGLN = IFRAGLN + 1
                  INITFRAG (IFRAG,IFRAGLN) = FRAGMAT (IFRAG,IATOM)
                  FRAGMAT  (IFRAG,IATOM)   = 0

                  IF (.NOT. ALREADY) THEN
                     NFRAG   = NFRAG + 1
                     ALREADY = .TRUE.
                  ENDIF
               ENDIF
            
            ENDDO

            IF (IFRAGLN .GT. 0) THEN
               ICOUNT = ICOUNT + 1
               FRAGLEN (ICOUNT) = IFRAGLN
               DO IATOM = 1, NATOMS
                  FRAGMAT (ICOUNT,IATOM) = INITFRAG (IFRAG,IATOM)
               ENDDO
            ENDIF

         ENDDO
C
C
C          ...Print out the information if user wants it!
C
C
         IF (IPRINT .GE. 10) THEN
            WRITE (*,*) ' @FRAG__WORK: Number of fragments - ',NFRAG

            WRITE (*,*) ' @FRAG__WORK: Fragmentation Array - '
            DO I = 1, NATOMS
               WRITE (*,'(99I3)') (FRAGMAT (I,J), J = 1, NATOMS)
            ENDDO
            WRITE (*,*) ''

            WRITE (*,*) ' @FRAG__WORK: Fragment Lengths - '
            WRITE (*,'(99I3)') (FRAGLEN (I), I = 1, NATOMS)
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END

