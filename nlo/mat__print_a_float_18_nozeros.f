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
         SUBROUTINE  MAT__PRINT_A_FLOAT_18_NOZEROS
     +
     +                    ( UNITID,
     +                      TITLE,
     +                      DDROW,DDCOL,
     +                      ROW,COL,
     +                      A )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__PRINT_A_FLOAT_18_NOZEROS
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation prints a two-dimensional matrix A to
C                the unit specified by UNITID in floating point format
C                F26.18 .Values below 5.0D-17 are not printed.
C
C                The following print options are possible:
C
C                    ROW equals zero   =   upper triangle of matrix
C                    COL equals zero   =   lower triangle of matrix
C                    ROW,COL not zero  =   full matrix
C
C
C                The matrix is printed by the following algorithm:
C
C                    a) divide number of columns into blocks of 4.
C                    b) loop over all blocks.
C                            loop over all rows in matrix.
C                                 loop over all 4 columns in block.
C                    c) determine number of columns left.
C                            loop over all rows in matrix.
C                                 loop over all resting columns. 
C                    
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    B
         INTEGER    BLOCKS
         INTEGER    DDROW,DDCOL
         INTEGER    I,J,K,L
         INTEGER    JSTART,JEND
         INTEGER    ROW,COL
         INTEGER    UNITID
         INTEGER    WIDTH

         CHARACTER*(*)      TITLE

         CHARACTER*26       S (1:4)

         DOUBLE PRECISION   LIMIT

         DOUBLE PRECISION   A (1:DDROW,1:DDCOL)

         DATA  LIMIT  /5.0D-17/
         DATA  WIDTH  /4/
C
C
C------------------------------------------------------------------------
C
C
C             ...printout title.
C
C
         WRITE (UNITID,5000)
         WRITE (UNITID, *  ) TITLE
C
C
C             ...immediate return if both lenghts are zero.
C
C
         IF  ( ROW.EQ.0  .AND.  COL.EQ.0 )  RETURN
C
C
C             ...check dimensions.
C
C
         IF  ( ROW .GT. DDROW   .OR.
     +         COL .GT. DDCOL         )  THEN

               WRITE (1,*) ' Dimensions of matrix A too small: '
               WRITE (1,*) ' mat__print_a_float_18_nozeros '
               WRITE (1,*) ' DDROW,DDCOL,ROW,COL = ',
     +                       DDROW,DDCOL,ROW,COL

               WRITE (*,*) ' Dimensions of matrix A too small: '
               WRITE (*,*) ' mat__print_a_float_18_nozeros '
               WRITE (*,*) ' DDROW,DDCOL,ROW,COL = ',
     +                       DDROW,DDCOL,ROW,COL

               STOP

         END IF

         IF  ( ROW.EQ.0 )  THEN
C
C
C             ...upper triangle of square matrix COL*COL is printed.
C
C
               BLOCKS = COL / WIDTH
               JEND   = 0

               DO  100  B = 1,BLOCKS

                   JSTART = JEND + 1
                   JEND   = JEND + WIDTH

                   WRITE (UNITID,6000) (J,J=JSTART,JEND)
                   WRITE (UNITID,7000)

                   DO  200  I = 1,JEND
                       K = 0
                       DO  300  J = JSTART,JEND
                           K = K + 1

                           IF  ( I.GT.J )  THEN
                                S(K) = ' '
                           ELSE IF ( ABS (A(I,J)) .LT. LIMIT )  THEN
                                S(K) = ' '
                           ELSE 
                                WRITE (S(K),8000) A(I,J)
                           END IF

  300                  CONTINUE
                       WRITE (UNITID,9000) I,(S(L),L=1,WIDTH)
  200              CONTINUE

  100          CONTINUE

               JSTART = JEND + 1

               IF  ( JSTART .GT. COL )  RETURN

               WRITE (UNITID,6000) (J,J=JSTART,COL)
               WRITE (UNITID,7000)

               DO  250  I = 1,COL
                   K = 0
                   DO  350  J = JSTART,COL
                       K = K + 1

                       IF  ( I.GT.J )  THEN
                             S(K) = ' '
                       ELSE IF  ( ABS (A(I,J)) .LT. LIMIT )  THEN
                             S(K) = ' '
                       ELSE 
                             WRITE (S(K),8000) A(I,J)
                       END IF

  350              CONTINUE
                   WRITE (UNITID,9000) I,(S(L),L=1,K)
  250          CONTINUE

               RETURN

         END IF

         IF  ( COL.EQ.0 )  THEN
C
C
C             ...lower triangle of square matrix ROW*ROW is printed.
C
C
               BLOCKS = ROW / WIDTH
               JEND   = 0

               DO  400  B = 1,BLOCKS

                   JSTART = JEND + 1
                   JEND   = JEND + WIDTH

                   WRITE (UNITID,6000) (J,J=JSTART,JEND)
                   WRITE (UNITID,7000)

                   DO  500  I = JSTART,ROW
                       K = 0
                       DO  600  J = JSTART,JEND
                           K = K + 1

                           IF  ( I.LT.J )  THEN
                                S(K) = ' '
                           ELSE IF ( ABS (A(I,J)) .LT. LIMIT )  THEN
                                S(K) = ' '
                           ELSE 
                                WRITE (S(K),8000) A(I,J)
                           END IF

  600                  CONTINUE
                       WRITE (UNITID,9000) I,(S(L),L=1,WIDTH)
  500              CONTINUE

  400          CONTINUE

               JSTART = JEND + 1

               WRITE (UNITID,6000) (J,J=JSTART,ROW)
               WRITE (UNITID,7000)

               DO  550  I = JSTART,ROW
                   K = 0
                   DO  650  J = JSTART,ROW
                       K = K + 1

                       IF  ( I.LT.J )  THEN
                             S(K) = ' '
                       ELSE IF  ( ABS (A(I,J)) .LT. LIMIT )  THEN
                             S(K) = ' '
                       ELSE 
                             WRITE (S(K),8000) A(I,J)
                       END IF

  650              CONTINUE
                   WRITE (UNITID,9000) I,(S(L),L=1,K)
  550          CONTINUE

               RETURN

         END IF
C
C
C             ...full matrix ROW*COL is printed.
C
C
         BLOCKS = COL / WIDTH
         JEND   = 0

         DO  700  B = 1,BLOCKS

             JSTART = JEND + 1
             JEND   = JEND + WIDTH

             WRITE (UNITID,6000) (J,J=JSTART,JEND)
             WRITE (UNITID,7000)

             DO  800  I = 1,ROW
                 K = 0
                 DO  900  J = JSTART,JEND
                     K = K + 1

                     IF  ( ABS (A(I,J)) .LT. LIMIT )  THEN
                           S(K) = ' '
                     ELSE 
                           WRITE (S(K),8000) A(I,J)
                     END IF

  900            CONTINUE
                 WRITE (UNITID,9000) I,(S(L),L=1,WIDTH)
  800        CONTINUE

  700    CONTINUE

         JSTART = JEND + 1

         IF  ( JSTART .GT. COL )  RETURN

         WRITE (UNITID,6000) (J,J=JSTART,COL)
         WRITE (UNITID,7000)

         DO  850  I = 1,ROW
             K = 0
             DO  950  J = JSTART,COL
                 K = K + 1

                 IF  ( ABS (A(I,J)) .LT. LIMIT )  THEN
                       S(K) = ' '
                 ELSE 
                       WRITE (S(K),8000) A(I,J)
                 END IF

  950        CONTINUE
             WRITE (UNITID,9000) I,(S(L),L=1,K)
  850    CONTINUE
C
C
C             ...formats of printing.
C
C
 5000    FORMAT  (/)
 6000    FORMAT  (/,1X,4I26)
 7000    FORMAT  ()
 8000    FORMAT  (F26.18)
 9000    FORMAT  (4X,I3,4A26)
C
C
C             ...ready!
C
C
         RETURN
         END
