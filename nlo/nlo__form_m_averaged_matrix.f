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
         SUBROUTINE  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWY,DDCOLY,
     +                      N,M,
     +                      ROWOFF,
     +                      MJUMP,LTRG,ADD,
     +                      X,
     +
     +                              Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_M_AVERAGED_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms the m-averaged matrix of size
C                N/M x N/M from an input matrix of size N x N. By
C                m-averaged it is meant that m-diagonal elements 
C                of X are summed over all different m-values and
C                the average is formed. Pictorially this can be
C                seen as follows. For an example in which m=x,y,z:
C
C
C                x y z x y z 
C                -----------
C             x |a    |d    |
C             y |  b  |  e  |                    (a+b+c)/3  (d+e+f)/3
C             z |    c|    f|   x,y,z averaged   
C                -----------   --------------->
C             x |g    |j    |    2 x 2 matrix
C             y |  h  |  k  |                    (g+h+i)/3  (j+k+l)/3
C             z |    i|    l|
C                -----------
C
C
C                If the x,y,z basis functions are ordered such that
C                each x,y,z are consecutive, then we would have:
C
C
C                x x y y z z 
C                -----------
C             x |a b|   |   |
C             x |c d|   |   |                    (a+e+i)/3  (b+f+j)/3
C                -----------     x,y,z averaged
C             y |   |e f|   |   --------------->
C             y |   |g h|   |     2 x 2 matrix
C                -----------                     (c+g+k)/3  (d+h+l)/3
C             z |   |   |i j|
C             z |   |   |k l|
C                ----------- 
C
C
C
C                Thus the m-average is only defined, if N is divisible
C                by M and can be evaluated in the above two ways,
C                depending on the m-order of the basis set. For each
C                m-averaged matrix element:
C
C
C                 1) Spaced form, taking M subblock diagonal elements
C                    from the X matrix spaced by M:
C
C                                           M
C                       Y (i,j) =  (1/M) * sum  X (M*[i-1]+k,M*[j-1]+k)
C                                           k
C
C                 2) Consecutive form, taking M subblock elements
C                    from the X matrix:
C
C                                           M
C                       Y (i,j) =  (1/M) * sum  X (i+M*[k-1],j+M*[k-1])
C                                           k
C
C                where i,j=1,N/M. In case N is not divisible by M,
C                the routine stops with an error message.
C
C                  Input:
C
C                    DDROWz,DDCOLz  =  declared dimensions of matrices
C                                      z = X and Y
C                                N  =  size of square N x N matrix X
C                                M  =  size of m-space
C                           ROWOFF  =  row offset value for picking
C                                      the N rows of N x N matrix X
C                            MJUMP  =  is .true., if the m values in
C                                      the m-space are ordered such
C                                      that the same m values are
C                                      separated. Hence this needs
C                                      evaluation of the spaced average.
C                                      If false, the consecutive average
C                                      form is generated
C                             LTRG  =  is .true., if only the lower
C                                      triangle of Y is wanted. In this
C                                      case also only the lower triangle
C                                      of X has to be passed. If false,
C                                      the full matrix Y is calculated
C                              ADD  =  is .true., if the averaged matrix
C                                      is to be added to existing values
C                                      in the Y matrix.
C                                X  =  matrix of size N x N from which
C                                      the average is calculated
C
C                  Output:
C
C                                Y  =  averaged matrix of size N/M x N/M
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

         LOGICAL     ADD,LTRG,MJUMP

         INTEGER     DDROWX,DDCOLX,DDROWY,DDCOLY
         INTEGER     I,J,K,L,M,N
         INTEGER     MSTEP,NSIZE
         INTEGER     ROWOFF

         DOUBLE PRECISION  FACTOR
         DOUBLE PRECISION  YIJ
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)
         DOUBLE PRECISION  Y (1:DDROWY,1:DDCOLY)

         PARAMETER   (ZERO = 0.D0 )
         PARAMETER   (ONE  = 1.D0 )
C
C
C------------------------------------------------------------------------
C
C
C             ...check, if N is divisible by M and the dimensions
C                of X and Y.
C
C
         IF  (MOD (N,M).NE.0)  THEN
              WRITE (1,*) ' Cannot form average matrix! '
              WRITE (1,*) ' nlo__form_m_averaged_matrix '
              WRITE (1,*) ' N,M = ',N,M
              WRITE (*,*) ' Cannot form average matrix! '
              WRITE (*,*) ' nlo__form_m_averaged_matrix '
              WRITE (*,*) ' N,M = ',N,M
              STOP
         END IF

         NSIZE = N / M

         IF  ((N+ROWOFF).GT.DDROWX .OR. N.GT.DDCOLX)  THEN
              WRITE (1,*) ' Dimensions of matrix X too small: '
              WRITE (1,*) ' nlo__form_m_averaged_matrix '
              WRITE (1,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
              WRITE (*,*) ' Dimensions of matrix X too small: '
              WRITE (*,*) ' nlo__form_m_averaged_matrix '
              WRITE (*,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
              STOP
         END IF

         IF  (NSIZE.GT.DDROWY .OR. NSIZE.GT.DDCOLY)  THEN
              WRITE (1,*) ' Dimensions of matrix Y too small: '
              WRITE (1,*) ' nlo__form_m_averaged_matrix '
              WRITE (1,*) ' DDROWY,DDCOLY,NSIZE = ',DDROWY,DDCOLY,NSIZE
              WRITE (*,*) ' Dimensions of matrix Y too small: '
              WRITE (*,*) ' nlo__form_m_averaged_matrix '
              WRITE (*,*) ' DDROWY,DDCOLY,NSIZE = ',DDROWY,DDCOLY,NSIZE
              STOP
         END IF
C
C
C             ...everything ok => form the m-average in desired form.
C
C
         FACTOR = ONE / DFLOAT (M)

         IF (MJUMP) THEN
             IF (LTRG) THEN
                 IF (ADD) THEN
C
C
C             ...lower triangle m-average in spaced form, adding
C                result to existing Y matrix.
C
C
                     L = 0
                     DO J = 1,NSIZE
                        K = L + ROWOFF
                        DO I = J,NSIZE
                           YIJ = ZERO
                           DO MSTEP = 1,M
                              YIJ = YIJ + X (K+MSTEP,L+MSTEP)
                           END DO
                           Y (I,J) = Y (I,J) + FACTOR * YIJ
                           K = K + M
                        END DO
                        L = L + M
                     END DO
                 ELSE
C
C
C             ...lower triangle m-average in spaced form, creating
C                new Y matrix.
C
C
                     L = 0
                     DO J = 1,NSIZE
                        K = L + ROWOFF
                        DO I = J,NSIZE
                           YIJ = ZERO
                           DO MSTEP = 1,M
                              YIJ = YIJ + X (K+MSTEP,L+MSTEP)
                           END DO
                           Y (I,J) = FACTOR * YIJ
                           K = K + M
                        END DO
                        L = L + M
                     END DO
                 END IF
             ELSE
                 IF (ADD) THEN
C
C
C             ...full m-average matrix in spaced form, adding result
C                to existing Y matrix.
C
C
                     L = 0
                     DO J = 1,NSIZE
                        K = ROWOFF
                        DO I = 1,NSIZE
                           YIJ = ZERO
                           DO MSTEP = 1,M
                              YIJ = YIJ + X (K+MSTEP,L+MSTEP)
                           END DO
                           Y (I,J) = Y (I,J) + FACTOR * YIJ
                           K = K + M
                        END DO
                        L = L + M
                     END DO
                 ELSE
C
C
C             ...full m-average matrix in spaced form, creating new
C                Y matrix.
C
C
                     L = 0
                     DO J = 1,NSIZE
                        K = ROWOFF
                        DO I = 1,NSIZE
                           YIJ = ZERO
                           DO MSTEP = 1,M
                              YIJ = YIJ + X (K+MSTEP,L+MSTEP)
                           END DO
                           Y (I,J) = FACTOR * YIJ
                           K = K + M
                        END DO
                        L = L + M
                     END DO
                 END IF
             END IF
         ELSE
             IF (LTRG) THEN
                 IF (ADD) THEN
C
C
C             ...lower triangle m-average in consecutive form, adding
C                result to existing Y matrix.
C
C
                     DO J = 1,NSIZE
                     DO I = J,NSIZE
                        K = I + ROWOFF
                        L = J
                        YIJ = ZERO
                        DO MSTEP = 1,M
                           YIJ = YIJ + X (K,L)
                           K = K + NSIZE
                           L = L + NSIZE
                        END DO
                        Y (I,J) = Y (I,J) + FACTOR * YIJ
                     END DO
                     END DO
                 ELSE
C
C
C             ...lower triangle m-average in consecutive form, creating
C                new Y matrix.
C
C
                     DO J = 1,NSIZE
                     DO I = J,NSIZE
                        K = I + ROWOFF
                        L = J
                        YIJ = ZERO
                        DO MSTEP = 1,M
                           YIJ = YIJ + X (K,L)
                           K = K + NSIZE
                           L = L + NSIZE
                        END DO
                        Y (I,J) = FACTOR * YIJ
                     END DO
                     END DO
                 END IF
             ELSE
                 IF (ADD) THEN
C
C
C             ...full m-average matrix in consecutive form, adding
C                result to existing Y matrix.
C
C
                     DO J = 1,NSIZE
                     DO I = 1,NSIZE
                        K = I + ROWOFF
                        L = J
                        YIJ = ZERO
                        DO MSTEP = 1,M
                           YIJ = YIJ + X (K,L)
                           K = K + NSIZE
                           L = L + NSIZE
                        END DO
                        Y (I,J) = Y (I,J) + FACTOR * YIJ
                     END DO
                     END DO
                 ELSE
C
C
C             ...full m-average matrix in consecutive form, creating
C                new Y matrix.
C
C
                     DO J = 1,NSIZE
                     DO I = 1,NSIZE
                        K = I + ROWOFF
                        L = J
                        YIJ = ZERO
                        DO MSTEP = 1,M
                           YIJ = YIJ + X (K,L)
                           K = K + NSIZE
                           L = L + NSIZE
                        END DO
                        Y (I,J) = FACTOR * YIJ
                     END DO
                     END DO
                 END IF
             END IF
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
