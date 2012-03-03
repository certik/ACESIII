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
         SUBROUTINE  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +                    ( N,M,
     +                      FIRST,LAST,
     +                      ABSOLUT,INCRESE,
     +                      OFF,
     +                      VECTOR,
     +
     +                              SORT )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__SORT_FLP_VECTOR_ELEMENTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine returns those indices of a vector
C                containing floating point elements, which correspond
C                to increasing (>=) or decreasing (=<) of (absolute)
C                vector values within a selected index range.
C                Schematically this can be shown by the following
C                picture:
C
C
C                    vector indices = 1 2 3 4 5 6 7 8 9 10 11 12
C
C                                           |            |
C                                         first         last
C
C                in which all vector values between index 4 and 10
C                will be sorted. The sorting technique used is
C                quicksort.
C
C
C                  Input:
C
C                    N            =  total # of vector elements
C                    M            =  total # of sorted vector indices
C                                    plus offset.
C                    FIRST,LAST   =  starting and ending indices
C                                    subjected to vector sorting.
C                    ABSOLUT      =  is true, if the absolute values
C                                    of the vector are to be taken
C                                    for ordering.
C                    INCRESE      =  is true, if the ordering is to
C                                    be performed according to
C                                    increasing (absolute) vector
C                                    values.
C                    OFF          =  offset value for placing sorted
C                                    vector indices.
C                    VECTOR (I)   =  the floating point vector elements.
C                                    Range of I is 1 to N.
C
C
C                  Output:
C
C                    SORT (OFF+I) =  contains the I-th increasing or
C                                    decreasing vector index, with I
C                                    ranging from 1 to LAST-FIRST+1.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ABSOLUT
         LOGICAL     INCRESE

         INTEGER     FACTOR
         INTEGER     FIRST,LAST
         INTEGER     I,J,M,N
         INTEGER     JBEG,JEND
         INTEGER     LEFT,RIGHT
         INTEGER     MAXIDX
         INTEGER     OFF,RELOFF
         INTEGER     SORTJ

         INTEGER     SORT  (1:M)

         DOUBLE PRECISION  DIFF
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  VECTOR (1:N)

         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (N.LT.LAST)  THEN
              WRITE (*,*) ' Last order index exceeds vector dimension! '
              WRITE (*,*) ' nlo__sort_flp_vector_elements '
              WRITE (*,*) ' N,LAST = ',N,LAST
              WRITE (1,*) ' Last order index exceeds vector dimension! '
              WRITE (1,*) ' nlo__sort_flp_vector_elements '
              WRITE (1,*) ' N,LAST = ',N,LAST
              STOP
         END IF

         RELOFF = OFF - FIRST + 1
         MAXIDX = RELOFF + LAST

         IF  (MAXIDX.GT.M)  THEN
              WRITE (*,*) ' Sorted index vector dimension too small! '
              WRITE (*,*) ' nlo__sort_flp_vector_elements '
              WRITE (*,*) ' largest index, dimension = ',MAXIDX,M
              WRITE (1,*) ' Sorted index vector dimension too small! '
              WRITE (1,*) ' nlo__sort_flp_vector_elements '
              WRITE (1,*) ' largest index, dimension = ',MAXIDX,M
              STOP
         END IF
C
C
C             ...perform quicksort on vector index values:
C
C                i) based on absolute flp vector values.
C
C
         IF (INCRESE) THEN
             FACTOR = +1
         ELSE
             FACTOR = -1
         END IF

         IF (ABSOLUT) THEN

             DO 100 I = FIRST,LAST
                LEFT = FIRST
                RIGHT = I - 1
  110           IF (LEFT .GT. RIGHT) THEN
                    JBEG = RELOFF + LEFT
                    IF (I .GT. LEFT) THEN
                        JEND = RELOFF + I - 1
                        DO 120 J = JEND,JBEG,-1
                           SORT (J + 1) = SORT (J)
  120                   CONTINUE
                    END IF
                    SORT (JBEG) = I
                ELSE
                    J = (LEFT + RIGHT) / 2
                    SORTJ = SORT (RELOFF + J)
                    DIFF = DABS (VECTOR (I)) - DABS (VECTOR (SORTJ))
                    IF ((FACTOR * DIFF) .LT. ZERO) THEN
                        RIGHT = J - 1
                    ELSE
                        LEFT = J + 1
                    END IF
                    GOTO 110
                END IF
  100        CONTINUE

         ELSE
C
C
C              ii) based on original flp vector values.
C
C
             DO 200 I = FIRST,LAST
                LEFT = FIRST
                RIGHT = I - 1
  210           IF (LEFT .GT. RIGHT) THEN
                    JBEG = RELOFF + LEFT
                    IF (I .GT. LEFT) THEN
                        JEND = RELOFF + I - 1
                        DO 220 J = JEND,JBEG,-1
                           SORT (J + 1) = SORT (J)
  220                   CONTINUE
                    END IF
                    SORT (JBEG) = I
                ELSE
                    J = (LEFT + RIGHT) / 2
                    SORTJ = SORT (RELOFF + J)
                    DIFF = VECTOR (I) - VECTOR (SORTJ)
                    IF ((FACTOR * DIFF) .LT. ZERO) THEN
                        RIGHT = J - 1
                    ELSE
                        LEFT = J + 1
                    END IF
                    GOTO 210
                END IF
  200        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
