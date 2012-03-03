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
         SUBROUTINE  NLO__HANDLE_MATRIX_SECTIONS
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWY,DDCOLY,
     +                      DDTSEC,DDNSEC,
     +                      NSEC,
     +                      SECIDX,SECDIM,SECOFF,
     +                      EXTRACT,
     +                      LTRG,
     +                      X,
     +
     +                              M,
     +                              Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__HANDLE_MATRIX_SECTIONS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : Given a square matrix X divided into sections, the
C                routine extracts or replaces the submatrix Y of
C                sections specified. The sections can be of different
C                lengths but must have the same decomposition along
C                the rows as well the columns of matrix X.
C
C                The following example (this one extracts/replaces
C                a lower triangle from/to X) makes this clear:
C
C
C                If X is:
C
C
C                               A       B       C      D    E
C             offset A ->   -----------------------------------
C                          |*        |     |         |   |     |
C                          |* *      |     |         |   |     |
C                        A |* * *    |     |         |   |     |
C                          |* * * *  |     |         |   |     |
C                          |* * * * *|     |         |   |     |
C             offset B ->   -----------------------------------
C                          |* * * * *|*    |         |   |     |
C                        B |* * * * *|* *  |         |   |     |
C                          |* * * * *|* * *|         |   |     |
C             offset C ->   -----------------------------------
C                          |         |     |         |   |     |
C                          |         |     |         |   |     |
C                        C |         |     |         |   |     |
C                          |         |     |         |   |     |
C                          |         |     |         |   |     |
C             offset D ->   -----------------------------------
C                          |* * * * *|* * *|         |*  |     |
C                        D |* * * * *|* * *|         |* *|     |
C             offset E ->   -----------------------------------
C                          |         |     |         |   |     |
C                        E |         |     |         |   |     |
C                          |         |     |         |   |     |
C                           -----------------------------------
C
C
C
C                then, if the above marked sections A,B and D are
C                extracted from the X matrix, Y is:
C
C
C                                A       B    D
C                            -------------------
C                           |*        |     |   |
C                           |* *      |     |   |
C                         A |* * *    |     |   |
C                           |* * * *  |     |   |
C                           |* * * * *|     |   |
C                            -------------------
C                           |* * * * *|*    |   |
C                         B |* * * * *|* *  |   |
C                           |* * * * *|* * *|   |
C                            -------------------
C                           |* * * * *|* * *|*  |
C                         D |* * * * *|* * *|* *|
C                            -------------------
C
C
C                and if Y is to be replaced into X, then the sections
C                of Y go into the above marked sections A,B and D of X.
C
C
C                Each section is characterized by two numbers:
C
C                     i) offset position inside the X matrix.
C                        This is the first index - 1 of where
C                        the section starts (see above picture
C                        of X). 
C
C                    ii) the # of rows/columns the section is
C                        made of.
C
C
C                  Input:
C
C                    DDROWz,DDCOLz  =  row and column dimensions for
C                                      matrices z=X,Y
C                    DDTSEC         =  dimension for total # of
C                                      sections
C                    DDNSEC         =  dimension for chosen # of
C                                      sections to be extracted or
C                                      replaced
C                    NSEC           =  # of sections to be extracted
C                                      from matrix X
C                    SECIDX         =  section labels indicating
C                                      which sections are going to
C                                      be extracted
C                    SECDIM         =  # of elements in each section
C                    SECOFF         =  offset index in matrix X for
C                                      each section
C                    EXTRACT        =  is true, if an extraction is
C                                      wanted, false if a replacement
C                                      is to be performed
C                    LTRG           =  is true, if only lower triangle
C                                      of X is to be extracted and
C                                      placed in lower triangle of Y
C                    X (if extract) =  matrix X from which the sections
C                                      are going to be extracted in
C                                      case of an extraction
C                    Y (if replace) =  matrix Y which is to be replaced
C
C
C                  Output:
C
C                    M (if extract) =  dimension of M x M extraction
C                                      matrix Y
C                    Y (if extract) =  extracted matrix Y
C                    X (if replace) =  replaced matrix X
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

         LOGICAL     EXTRACT
         LOGICAL     LTRG

         INTEGER     DDROWX,DDCOLX,DDROWY,DDCOLY
         INTEGER     DDNSEC,DDTSEC
         INTEGER     DIM
         INTEGER     I,J,K,L,M,N
         INTEGER     IDX,IXI,IXJ
         INTEGER     NCOL,NROW
         INTEGER     NSEC
         INTEGER     OLDOFF,NEWOFF
         INTEGER     XCOFF,YCOFF,XROFF,YROFF
         INTEGER     XCOL,YCOL

         INTEGER     SECDIM (1:DDTSEC)
         INTEGER     SECIDX (1:DDNSEC)
         INTEGER     SECOFF (1:DDTSEC)

         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)
         DOUBLE PRECISION  Y (1:DDROWY,1:DDCOLY)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine any inconsistencies in the sections
C                provided and check dimensions of X and Y matrices
C                supplied.
C
C
         IF (NSEC.GT.DDNSEC) THEN
             WRITE (*,*) ' Dimensions for chosen # of sects too small! '
             WRITE (*,*) ' nlo__handle_matrix_sections '
             WRITE (*,*) ' DDNSEC,NSEC = ',DDNSEC,NSEC
             WRITE (1,*) ' Dimensions for chosen # of sects too small! '
             WRITE (1,*) ' nlo__handle_matrix_sections '
             WRITE (1,*) ' DDNSEC,NSEC = ',DDNSEC,NSEC
             STOP
         END IF

         M = 0
         N = 0
         OLDOFF = 0
         DO 10 I = 1,NSEC
            IDX = SECIDX (I)
            DIM = SECDIM (IDX)
            M = M + DIM
            NEWOFF = SECOFF (IDX)
            IF (NEWOFF .LT. OLDOFF) THEN
                WRITE (*,*) ' Overlapping sections for extraction! '
                WRITE (*,*) ' Occuring at section # ',I
                WRITE (*,*) ' nlo__handle_matrix_sections '
                WRITE (1,*) ' Overlapping sections for extraction! '
                WRITE (1,*) ' Occuring at section # ',I
                WRITE (1,*) ' nlo__handle_matrix_sections '
                STOP
            END IF
            OLDOFF = NEWOFF + DIM
            N = MAX0 (N,OLDOFF)
   10    CONTINUE

         IF (N.GT.DDROWX .OR. N.GT.DDCOLX) THEN
             WRITE (*,*) ' Dimensions of matrix X too small! '
             WRITE (*,*) ' nlo__handle_matrix_sections '
             WRITE (*,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             WRITE (1,*) ' Dimensions of matrix X too small! '
             WRITE (1,*) ' nlo__handle_matrix_sections '
             WRITE (1,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             STOP
         END IF

         IF (M.GT.DDROWY .OR. M.GT.DDCOLY) THEN
             WRITE (*,*) ' Dimensions of matrix Y too small! '
             WRITE (*,*) ' nlo__handle_matrix_sections '
             WRITE (*,*) ' DDROWY,DDCOLY,M = ',DDROWY,DDCOLY,M
             WRITE (1,*) ' Dimensions of matrix Y too small! '
             WRITE (1,*) ' nlo__handle_matrix_sections '
             WRITE (1,*) ' DDROWY,DDCOLY,M = ',DDROWY,DDCOLY,M
             STOP
         END IF
C
C
C             ...decide on extraction or replacement.
C                Outer loops will run over column sections and
C                elements, inner loops over row sections and
C                elements.
C
C
         IF (EXTRACT) THEN
C
C
C             ...lower triangular extraction.
C
C
             IF (LTRG) THEN

                 YCOFF = 0
                 DO 100 J = 1,NSEC
                    IXJ = SECIDX (J)
                    NCOL = SECDIM (IXJ)
                    XCOFF = SECOFF (IXJ)
                    DO 110 L = 1,NCOL
                       XCOL = XCOFF + L
                       YCOL = YCOFF + L

                       DO 120 K = L,NCOL
                          Y (YCOFF+K,YCOL) = X (XCOFF+K,XCOL)
  120                  CONTINUE

                       YROFF = YCOFF + NCOL
                       DO 130 I = J+1,NSEC
                          IXI = SECIDX (I)
                          NROW = SECDIM (IXI)
                          XROFF = SECOFF (IXI)
                          DO 140 K = 1,NROW
                             Y (YROFF+K,YCOL) = X (XROFF+K,XCOL)
  140                     CONTINUE
                          YROFF = YROFF + NROW
  130                  CONTINUE

  110               CONTINUE
                    YCOFF = YCOFF + NCOL
  100            CONTINUE

             ELSE
C
C
C             ...full matrix extraction.
C
C
                 YCOFF = 0
                 DO 200 J = 1,NSEC
                    IXJ = SECIDX (J)
                    NCOL = SECDIM (IXJ)
                    XCOFF = SECOFF (IXJ)
                    DO 210 L = 1,NCOL
                       XCOL = XCOFF + L
                       YCOL = YCOFF + L

                       YROFF = 0
                       DO 220 I = 1,NSEC
                          IXI = SECIDX (I)
                          NROW = SECDIM (IXI)
                          XROFF = SECOFF (IXI)
                          DO 230 K = 1,NROW
                             Y (YROFF+K,YCOL) = X (XROFF+K,XCOL)
  230                     CONTINUE
                          YROFF = YROFF + NROW
  220                  CONTINUE

  210               CONTINUE
                    YCOFF = YCOFF + NCOL
  200            CONTINUE

             END IF

         ELSE
C
C
C             ...lower triangular replacement.
C
C
             IF (LTRG) THEN

                 YCOFF = 0
                 DO 300 J = 1,NSEC
                    IXJ = SECIDX (J)
                    NCOL = SECDIM (IXJ)
                    XCOFF = SECOFF (IXJ)
                    DO 310 L = 1,NCOL
                       XCOL = XCOFF + L
                       YCOL = YCOFF + L

                       DO 320 K = L,NCOL
                          X (XCOFF+K,XCOL) = Y (YCOFF+K,YCOL)
  320                  CONTINUE

                       YROFF = YCOFF + NCOL
                       DO 330 I = J+1,NSEC
                          IXI = SECIDX (I)
                          NROW = SECDIM (IXI)
                          XROFF = SECOFF (IXI)
                          DO 340 K = 1,NROW
                             X (XROFF+K,XCOL) = Y (YROFF+K,YCOL)
  340                     CONTINUE
                          YROFF = YROFF + NROW
  330                  CONTINUE

  310               CONTINUE
                    YCOFF = YCOFF + NCOL
  300            CONTINUE

             ELSE
C
C
C             ...full matrix replacement.
C
C
                 YCOFF = 0
                 DO 400 J = 1,NSEC
                    IXJ = SECIDX (J)
                    NCOL = SECDIM (IXJ)
                    XCOFF = SECOFF (IXJ)
                    DO 410 L = 1,NCOL
                       XCOL = XCOFF + L
                       YCOL = YCOFF + L

                       YROFF = 0
                       DO 420 I = 1,NSEC
                          IXI = SECIDX (I)
                          NROW = SECDIM (IXI)
                          XROFF = SECOFF (IXI)
                          DO 430 K = 1,NROW
                             X (XROFF+K,XCOL) = Y (YROFF+K,YCOL)
  430                     CONTINUE
                          YROFF = YROFF + NROW
  420                  CONTINUE

  410               CONTINUE
                    YCOFF = YCOFF + NCOL
  400            CONTINUE

             END IF
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
