         SUBROUTINE  FRAG__LOOP_THROUGH
     +
     +                   (NATOMS,IFRAG,
     +                    HOT,IHOT,
     +                    HOTEQNATOM,
     +                    INCFRAG,
     +                    IPRINT,IERR,
     +                    INITFRAG,
     +                                 FRAGMAT )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : FRAG__LOOP_THROUGH
C  MODULE      : FRAGMENT DETERMINE
C  MODULE-ID   : FRAG
C  SUBROUTINES :
C
C  DESCRIPTION : This flag will check another iteration in the
C                bonding pattern and determine if the current
C                fragment is filled.
C
C                  Input:
C
C                    IPRINT       =  level of printing variable
C                    IERR         =  sets the return error
C                    IFRAG        =  current fragment number
C                    HOT          =  current atom we're checking
C                    IHOT         =  pointer to the current atom
C                    HOTEQNATOM   =  if HOT is the last atom in 
C                                    this fragment iteration
C                    INCFRAG      =  flag that tells if we've
C                                    exhausted the fragment and
C                                    atom space against HOT
C                    NATOMS       =  total number of atoms
C                    INITFRAG     =  array to hold previous
C                                    fragmentation pattern
C
C                  Output:
C
C                    FRAGMAT      =  contains the fragments for
C                                    an iteration
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

         INTEGER    I,J,K,IHOT
         INTEGER    IHOT2,HOT2
         INTEGER    IFRAG,IATOM,ILOC
         INTEGER    ONE,BOND,HOT,NEW
         INTEGER    IPRINT,IERR
         INTEGER    NATOMS,LSCR1
         INTEGER    ZERO

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
C             ...Loop through and move all overlapped 
C                elements into the proper IFRAG.
C
C
         IF (HOT .GT. ZERO) THEN

            DO I = 1, NATOMS
            DO J = 1, NATOMS
               IF (I .NE. IFRAG) THEN
                  NEW = INITFRAG (I,J)
                  IF (NEW .EQ. HOT) THEN
                     DO K = 1, NATOMS
                        IF (INITFRAG (I,K) .GT. 0) THEN
                           FRAGMAT  (IFRAG,K) = INITFRAG (I,K)
                           INITFRAG (I,K) = ZERO
                           FRAGMAT  (I,K) = ZERO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            ENDDO

         ENDIF ! (HOT .GT. ZERO)
C
C
C             ...Now look through for at least one more
C                overlap.  If it exists, we do not increment
C                the fragment.
C
C
         IF (HOTEQNATOM) THEN

            INCFRAG = .TRUE.

            DO IHOT2 = IFRAG, NATOMS
               HOT2 = FRAGMAT (IFRAG,IHOT2)
               IF (HOT2 .GT. ZERO) THEN
                  DO I = 1, NATOMS
                  DO J = 1, NATOMS
                     IF (I .NE. IFRAG) THEN
                        NEW = INITFRAG (I,J)
                        IF (NEW .EQ. HOT2) THEN 
                           INCFRAG = .FALSE.
                           RETURN
                        ENDIF
                     ENDIF
                  ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END

