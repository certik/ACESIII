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
         SUBROUTINE  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      DDWORK,
     +                      ROW,SUM,COL,
     +                      ROWOFF,COLOFF,
     +                      SAVEA,SAVEB,SAVEC,
     +                      A,B,C,
     +                      WORK,
     +
     +                              D )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_ORTHOTRAN_OFFDIAG
C  MODULE      : Matrix
C  MODULE-ID   : MAT
C  DESCRIPTION : This routine generates an offdiagonal block of the
C                orthogonally transformed matrix B:
C
C                              D = A(T) * B * C
C
C                for which we can write for each element:
C
C                          D(ij)  =  sum  A(ki) * B(kl) * C(lj)
C                                     kl
C
C                where the index ranges are:
C
C                           i index    ->   1 to ROW
C                           j index    ->   1 to COL
C                         k=l indices  ->   1 to SUM
C
C                The matrix A(T) is passed as such in argument, that
C                is the array A contains the SUM x ROW A(T) elements.
C                Note, that the working array WORK is only needed in
C                case all three matrices A,B,C need to be saved.
C
C                There are five possible effective paths the routine
C                can take:
C
C
C                   1) If matrix A can be destroyed (SAVEA = .false.):
C                      ----------------------------------------------
C
C                          i) A(T) * B -> store in D (size SUM x ROW)
C                         ii) Copy D to A(T)
C                        iii) A(T) * C -> store in D
C
C
C                   2) If matrix B can be destroyed (SAVEB = .false.):
C                      ----------------------------------------------
C
C                          i) A(T) * B -> store in D (size SUM x ROW)
C                         ii) Copy D to B
C                        iii) B(T) * C -> store in D
C
C
C                   3) If matrix B can be destroyed (SAVEB = .false.):
C                      ----------------------------------------------
C
C                          i) B * C -> store in D (size SUM x COL)
C                         ii) Copy D to B
C                        iii) A(T) * B -> store in D
C
C
C                   4) If matrix C can be destroyed (SAVEC = .false.):
C                      ----------------------------------------------
C
C                          i) B * C -> store in D (size SUM x COL)
C                         ii) Copy D to C
C                        iii) A(T) * C -> store in D
C
C
C                   5) If all matrices A,B,C have to be saved:
C                      --------------------------------------
C
C                          i) B * C -> store in D (size SUM x COL)
C                         ii) A(T) * D -> overwrite in D using WORK
C
C
C                These five main paths need different intermediate sizes
C                of matrix D, hence fitting of a path according to
C                declared D matrix dimensions should always be checked.
C
C                There is also the option of placing the resulting
C                rectangular matrix inside the transmitted D array. The
C                final transformed ROW x COL matrix can be placed
C                starting at a specified row index with offset ROWOFF
C                and at a specified column index with offset COLOFF.
C                Pictorially this looks as follows inside the D array:
C
C
C                                     COLOFF + 1    COLOFF + COL
C
C                                    ---------------------------
C                                   |        .           .      |
C                                   |        .           .      |
C                                   |        .           .      |
C                                   |        .           .      |
C                  ROWOFF + 1   ->  |. . . . |-----------|      |
C                                   |        |/ / / / / /|      |
C                                   |        |/ D block /|      |
C                                   |        |/ / / / / /|      |
C                  ROWOFF + ROW ->  |. . . . |-----------|      |
C                                   |                           |
C                                   |                           |
C                                   |                           |
C                                   |                           |
C                                    ---------------------------
C
C
C                If for some reason one wants the D block to start at
C                the first row, simply pass ROWOFF = 0 in the argument.
C                Same with first column with COLOFF = 0.
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

         LOGICAL     FITBSR,FITBSC,FITDSR,FITDSC
         LOGICAL     OFFSET
         LOGICAL     ROWCOL
         LOGICAL     SAVEA,SAVEB,SAVEC,SAVALL
         LOGICAL     PATH1,PATH2,PATH3,PATH4,PATH5

         INTEGER     DDROWA,DDCOLA
         INTEGER     DDROWB,DDCOLB
         INTEGER     DDROWC,DDCOLC
         INTEGER     DDROWD,DDCOLD
         INTEGER     DDWORK
         INTEGER     I,J,K
         INTEGER     ROW,SUM,COL
         INTEGER     ROWOFF,COLOFF

         DOUBLE PRECISION  X,ZERO

         DOUBLE PRECISION  WORK (1:DDWORK)

         DOUBLE PRECISION  A (1:DDROWA,1:DDCOLA)
         DOUBLE PRECISION  B (1:DDROWB,1:DDCOLB)
         DOUBLE PRECISION  C (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  D (1:DDROWD,1:DDCOLD)

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions for:
C
C                      i) matrix A
C                     ii) matrix B
C                    iii) matrix C
C                     iv) matrix D holding final result
C
C
         IF  (SUM.GT.DDROWA .OR. ROW.GT.DDCOLA)  THEN
              WRITE (*,*) ' Dimensions of matrix A too small! '
              WRITE (*,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (*,*) ' DDROWA,DDCOLA,ROW,SUM = ',
     +                      DDROWA,DDCOLA,ROW,SUM
              WRITE (1,*) ' Dimensions of matrix A too small! '
              WRITE (1,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (1,*) ' DDROWA,DDCOLA,ROW,SUM = ',
     +                      DDROWA,DDCOLA,ROW,SUM
              STOP
         END IF

         IF  (SUM.GT.DDROWB .OR. SUM.GT.DDCOLB)  THEN
              WRITE (*,*) ' Dimensions of matrix B too small! '
              WRITE (*,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (*,*) ' DDROWB,DDCOLB,SUM = ',DDROWB,DDCOLB,SUM
              WRITE (1,*) ' Dimensions of matrix B too small! '
              WRITE (1,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (1,*) ' DDROWB,DDCOLB,SUM = ',DDROWB,DDCOLB,SUM
              STOP
         END IF

         IF  (SUM.GT.DDROWC .OR. COL.GT.DDCOLC)  THEN
              WRITE (*,*) ' Dimensions of matrix C too small! '
              WRITE (*,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (*,*) ' DDROWC,DDCOLC,SUM,COL = ',
     +                      DDROWC,DDCOLC,SUM,COL
              WRITE (1,*) ' Dimensions of matrix C too small! '
              WRITE (1,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (1,*) ' DDROWC,DDCOLC,SUM,COL = ',
     +                      DDROWC,DDCOLC,SUM,COL
              STOP
         END IF

         IF  ((ROW+ROWOFF).GT.DDROWD .OR. (COL+COLOFF).GT.DDCOLD)  THEN
              WRITE (*,*) ' Matrix D cannot hold final result! '
              WRITE (*,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (*,*) ' DDROWD,DDCOLD,ROW+ROWOFF,COL+COLOFF = ',
     +                      DDROWD,DDCOLD,ROW+ROWOFF,COL+COLOFF
              WRITE (1,*) ' Matrix D cannot hold final result! '
              WRITE (1,*) ' mat__c_eq_orthotran_offdiag '
              WRITE (1,*) ' DDROWD,DDCOLD,ROW+ROWOFF,COL+COLOFF = ',
     +                      DDROWD,DDCOLD,ROW+ROWOFF,COL+COLOFF
              STOP
         END IF
C
C
C             ...find the allowed paths. The variables FITDSR and
C                FITDSC indicate, if the intermediate products of
C                sizes SUM x ROW and SUM x COL, respectively, can
C                be fit into the D matrix. The variables FITBSR
C                and FITBSC indicate, if the intermediate products
C                can be copied to the B matrix.
C
C                Note, that when the # of rows is greater than the
C                # of columns, the paths 1) and 2) are not as favorable
C                as paths 3),4) or 5) due to the intermediate storage
C                and copying of D, which has then more elements in the
C                first two paths. Hence we exclude paths 1) and 2)
C                from the list, if either paths 3),4) or 5) are
C                allowed.
C
C
         PATH1  = .FALSE.
         PATH2  = .FALSE.
         PATH3  = .FALSE.
         PATH4  = .FALSE.
         PATH5  = .FALSE.

         SAVALL = SAVEA .AND. SAVEB .AND. SAVEC

         FITBSR =  ROW.LE.DDCOLB
         FITBSC =  COL.LE.DDCOLB
         FITDSR = (SUM.LE.DDROWD .AND. ROW.LE.DDCOLD)
         FITDSC = (SUM.LE.DDROWD .AND. COL.LE.DDCOLD)

         ROWCOL = ROW.GT.COL
         OFFSET = (ROWOFF.NE.0 .OR. COLOFF.NE.0)

         IF  (.NOT.SAVEA .AND. FITDSR             )  PATH1  = .TRUE.
         IF  (.NOT.SAVEB .AND. FITDSR .AND. FITBSR)  PATH2  = .TRUE.
         IF  (.NOT.SAVEB .AND. FITDSC .AND. FITBSC)  PATH3  = .TRUE.
         IF  (.NOT.SAVEC .AND. FITDSC             )  PATH4  = .TRUE.
         IF  (  SAVALL   .AND. FITDSC             )  PATH5  = .TRUE.

         IF (ROWCOL .AND. (PATH3.OR.PATH4.OR.PATH5)) THEN
             PATH1  = .FALSE.
             PATH2  = .FALSE.
         END IF
C
C
C             ...proceed with the path chosen.
C
C                Path # 1: matrix in array A can be destroyed.
C
C                       i) A(T) * B -> store in D (size SUM x ROW)
C                      ii) Copy D to A(T)
C                     iii) A(T) * C -> store in D
C
C
         IF (PATH1) THEN

             CALL    MAT__CT_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWD,DDCOLD,
     +                      ROW,SUM,SUM,
     +                      A,B,
     +
     +                              D )
     +
     +
             CALL    MAT__C_EQ_A_FLOAT
     +
     +                    ( DDROWD,DDCOLD,
     +                      DDROWA,DDCOLA,
     +                      SUM,ROW,
     +                      D,
     +
     +                              A )
     +
     +
             CALL    MAT__C_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      ROW,SUM,COL,
     +                      A,C,
     +
     +                              D )
     +
     +
             IF (OFFSET) THEN
                 DO J = 1,COL
                 DO I = 1,ROW
                    D (ROWOFF+I,COLOFF+J) = D (I,J)
                 END DO
                 END DO
             END IF
C
C
C             ...path # 2: matrix in array B can be destroyed.
C
C                       i) A(T) * B -> store in D (size SUM x ROW)
C                      ii) Copy D to B
C                     iii) B(T) * C -> store in D
C
C
         ELSE IF (PATH2) THEN

             CALL    MAT__CT_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWD,DDCOLD,
     +                      ROW,SUM,SUM,
     +                      A,B,
     +
     +                              D )
     +
     +
             CALL    MAT__C_EQ_A_FLOAT
     +
     +                    ( DDROWD,DDCOLD,
     +                      DDROWB,DDCOLB,
     +                      SUM,ROW,
     +                      D,
     +
     +                              B )
     +
     +
             CALL    MAT__C_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      ROW,SUM,COL,
     +                      B,C,
     +
     +                              D )
     +
     +
             IF (OFFSET) THEN
                 DO J = 1,COL
                 DO I = 1,ROW
                    D (ROWOFF+I,COLOFF+J) = D (I,J)
                 END DO
                 END DO
             END IF
C
C
C             ...path # 3: matrix in array B can be destroyed.
C
C                       i) B * C -> store in D (size SUM x COL)
C                      ii) Copy D to B
C                     iii) A(T) * B -> store in D
C
C
         ELSE IF (PATH3) THEN

             CALL    MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                    ( DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      SUM,SUM,COL,
     +                      B,C,
     +
     +                              D )
     +
     +
             CALL    MAT__C_EQ_A_FLOAT
     +
     +                    ( DDROWD,DDCOLD,
     +                      DDROWB,DDCOLB,
     +                      SUM,COL,
     +                      D,
     +
     +                              B )
     +
     +
             CALL    MAT__C_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWD,DDCOLD,
     +                      ROW,SUM,COL,
     +                      A,B,
     +
     +                              D )
     +
     +
             IF (OFFSET) THEN
                 DO J = 1,COL
                 DO I = 1,ROW
                    D (ROWOFF+I,COLOFF+J) = D (I,J)
                 END DO
                 END DO
             END IF
C
C
C             ...path # 4: matrix in array C can be destroyed.
C
C                       i) B * C -> store in D (size SUM x COL)
C                      ii) Copy D to C
C                     iii) A(T) * C -> store in D
C
C
         ELSE IF (PATH4) THEN

             CALL    MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                    ( DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      SUM,SUM,COL,
     +                      B,C,
     +
     +                              D )
     +
     +
             CALL    MAT__C_EQ_A_FLOAT
     +
     +                    ( DDROWD,DDCOLD,
     +                      DDROWC,DDCOLC,
     +                      SUM,COL,
     +                      D,
     +
     +                              C )
     +
     +
             CALL    MAT__C_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      ROW,SUM,COL,
     +                      A,C,
     +
     +                              D )
     +
     +
             IF (OFFSET) THEN
                 DO J = 1,COL
                 DO I = 1,ROW
                    D (ROWOFF+I,COLOFF+J) = D (I,J)
                 END DO
                 END DO
             END IF
C
C
C             ...path # 4: matrices in arrays A,B,C need to be saved.
C
C                       i) B * C -> store in D (size SUM x COL)
C                      ii) A(T) * D -> overwrite in D using WORK
C
C
         ELSE IF (PATH5) THEN

             IF (ROW.GT.DDWORK) THEN
                 WRITE (*,*) ' Dimensions of working array too small! '
                 WRITE (*,*) ' mat__c_eq_orthotran_offdiag '
                 WRITE (*,*) ' DDWORK,ROW = ',DDWORK,ROW
                 WRITE (1,*) ' Dimensions of working array too small! '
                 WRITE (1,*) ' mat__c_eq_orthotran_offdiag '
                 WRITE (1,*) ' DDWORK,ROW = ',DDWORK,ROW
                 STOP
             END IF

             CALL    MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                    ( DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      SUM,SUM,COL,
     +                      B,C,
     +
     +                              D )
     +
     +
             DO J = 1,COL
                DO I = 1,ROW
                   X = ZERO
                   DO K = 1,SUM
                      X = X + A (K,I) * D (K,J)
                   END DO
                   WORK (I) = X
                END DO
                DO I = 1,ROW
                   D (ROWOFF+I,J) = WORK (I)
                END DO
             END DO

             IF (COLOFF.NE.0) THEN
                 DO J = 1,COL
                 DO I = 1,ROW
                    D (ROWOFF+I,COLOFF+J) = D (ROWOFF+I,J)
                 END DO
                 END DO
             END IF
C
C
C             ...no path possible. Print error message.
C
C
         ELSE
             WRITE (*,*) ' No path possible! Cannot do operation! '
             WRITE (*,*) ' mat__c_eq_orthotran_offdiag '
             WRITE (1,*) ' No path possible! Cannot do operation! '
             WRITE (1,*) ' mat__c_eq_orthotran_offdiag '
             STOP
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
