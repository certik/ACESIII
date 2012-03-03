         SUBROUTINE  FRAG__INIT
     +
     +                   (NATOMS,
     +                    CONNECT,
     +                    IPRINT,IERR,
     +                    INITFRAG,
     +
     +                                 NFRAG,
     +                                 FRAGLEN,
     +                                 FRAGMAT  )
     +
     +             
C------------------------------------------------------------------------
C  OPERATION   : FRAG__INIT
C  MODULE      : FRAGMENT DETERMINE
C  MODULE-ID   : FRAG
C  SUBROUTINES : FRAG__WORK
C                                    
C  DESCRIPTION : Main program that drives the fragmentation pattern.
C
C                  Input:
C
C                    IPRINT       =  level of printing variable
C                    IERR         =  sets the return error
C                    NATOMS       =  total number of atoms
C                    CONNECT      =  connectivity array
C                                    fragmentation pattern
C                    INITFRAG     =  scratch space for routine
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

         INTEGER    I,J,ZERO
         INTEGER    IPRINT,IERR
         INTEGER    NATOMS,NFRAG

         INTEGER    FRAGLEN  (1:NATOMS)

         INTEGER    CONNECT  (1:NATOMS,1:NATOMS)
         INTEGER    FRAGMAT  (1:NATOMS,1:NATOMS)
         INTEGER    INITFRAG (1:NATOMS,1:NATOMS)

         PARAMETER  (ZERO = 0)
C
C
C------------------------------------------------------------------------
C
C
C             ...Return if 1 atom, analyze otherwise!
C
C
         IF (NATOMS .LE. 1) THEN
            NFRAG   = 1
            FRAGLEN ( 1 ) = 1
            FRAGMAT (1,1) = 1
            RETURN
         ENDIF
C
C
C          ...Zero out the necessary arrays!
C
C
         DO I = 1, NATOMS
            FRAGLEN (I) = ZERO
            DO J = 1, NATOMS
               INITFRAG (I,J) = ZERO
               FRAGMAT  (I,J) = ZERO
            ENDDO
         ENDDO
C
C
C          ...Figure out the fragmentation pattern!
C
C
         CALL  FRAG__WORK
     +
     +              (NATOMS,
     +               CONNECT,
     +               IPRINT,IERR,
     +               INITFRAG,
     +
     +                          NFRAG,
     +                          FRAGLEN,
     +                          FRAGMAT  )
     +
C
C
C             ...ready!
C
C
         RETURN
         END
