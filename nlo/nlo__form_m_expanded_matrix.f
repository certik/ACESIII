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
         SUBROUTINE  NLO__FORM_M_EXPANDED_MATRIX
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWY,DDCOLY,
     +                      M,N,
     +                      ROWOFF,
     +                      MJUMP,LTRG,ZEROY,
     +                      X,
     +
     +                              Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_M_EXPANDED_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms the new m-expanded matrix Y of size
C                N x N from an input matrix X of size M x M. The
C                m-expansion of a matrix is only defined, if N is
C                divisible by M. The following example for m=x,y,z
C                makes the procedure clear:
C
C
C                                                   x y z x y z 
C                                                   -----------
C                                                x |a    |b    |
C                                                y |  a  |  b  |
C                  a  b     x,y,z expanded       z |    a|    b|
C                         ------------------>       ----------- 
C                  c  d      6 x 6 matrix        x |c    |d    |
C                                                y |  c  |  d  |
C                                                z |    c|    d|
C                                                   -----------
C
C                X matrix                            Y matrix
C
C
C                If the x,y,z basis functions are ordered such that
C                each x,y,z are consecutive, then we would have the
C                following m-expansion:
C
C
C                                                   x x y y z z 
C                                                   -----------
C                                                x |a b|   |   |
C                                                x |c d|   |   |
C                                                   -----------
C                  a  b     x,y,z expanded       y |   |a b|   |
C                         ------------------>    y |   |c d|   |
C                  c  d      6 x 6 matrix           -----------
C                                                z |   |   |a b|
C                                                z |   |   |c d|
C                                                   -----------
C
C                X matrix                            Y matrix
C
C
C                Thus the m-expansion can proceed via two possible
C                ways, depending on the m-order of the basis set.
C                For each Y matrix element:
C
C
C                 1) Spaced form:
C
C                       Y (r+[N/M]*[k-1],r+[N/M]*[l-1])  =  X (k,l)
C
C                 2) Consecutive form:
C
C                       Y (k+M*[r-1],l+M*[r-1])  =  X (k,l)
C
C
C                where r=1,[N/M]. In case N is not divisible by M,
C                the routine stops with an error message.
C
C                  Input:
C
C                    DDROWz,DDCOLz  =  declared dimensions of matrices
C                                      z = X and Y
C                              M,N  =  # of rows and columns of matrices
C                                      X and Y
C                           ROWOFF  =  row offset value for placing
C                                      the N rows of N x N matrix Y
C                            MJUMP  =  is true, if the m values in
C                                      the m-space are ordered such
C                                      that the same m values are
C                                      separated. Hence this needs
C                                      evaluation of the spaced
C                                      expansion. If false, the
C                                      consecutive expansion is done
C                             LTRG  =  is true, if only the lower
C                                      triangle of Y is wanted. In this
C                                      case also only the lower triangle
C                                      of X has to be passed. If false,
C                                      the full matrix Y is calculated
C                            ZEROY  =  is true, if the output Y matrix
C                                      needs to be zeroed first. This
C                                      is necessary because not all
C                                      elements of Y are addressed
C                                      during the m-expansion. If false,
C                                      the zeroing of Y is bypassed
C                                X  =  original matrix of size M x M
C
C                  Output:
C
C                                Y  =  m-expanded matrix of sixe N x N
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

         LOGICAL     MJUMP,LTRG,ZEROY

         INTEGER     DDROWX,DDCOLX,DDROWY,DDCOLY
         INTEGER     I,J,K,L,M,N,R
         INTEGER     NM
         INTEGER     ROWOFF

         DOUBLE PRECISION  XKL
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)
         DOUBLE PRECISION  Y (1:DDROWY,1:DDCOLY)

         PARAMETER   (ZERO = 0.D0)
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
              WRITE (*,*) ' Cannot form expanded matrix! '
              WRITE (*,*) ' nlo__form_m_expanded_matrix '
              WRITE (*,*) ' N,M = ',N,M
              WRITE (1,*) ' Cannot form expanded matrix! '
              WRITE (1,*) ' nlo__form_m_expanded_matrix '
              WRITE (1,*) ' N,M = ',N,M
              STOP
         END IF

         IF  (M.GT.DDROWX .OR. M.GT.DDCOLX)  THEN
              WRITE (*,*) ' Dimensions of matrix X too small: '
              WRITE (*,*) ' nlo__form_m_expanded_matrix '
              WRITE (*,*) ' DDROWX,DDCOLX,M = ',DDROWX,DDCOLX,M
              WRITE (1,*) ' Dimensions of matrix X too small: '
              WRITE (1,*) ' nlo__form_m_expanded_matrix '
              WRITE (1,*) ' DDROWX,DDCOLX,M = ',DDROWX,DDCOLX,M
              STOP
         END IF

         IF  ((N+ROWOFF).GT.DDROWY .OR. N.GT.DDCOLY)  THEN
              WRITE (*,*) ' Dimensions of matrix Y too small: '
              WRITE (*,*) ' nlo__form_m_expanded_matrix '
              WRITE (*,*) ' DDROWY,DDCOLY,N = ',DDROWY,DDCOLY,N
              WRITE (1,*) ' Dimensions of matrix Y too small: '
              WRITE (1,*) ' nlo__form_m_expanded_matrix '
              WRITE (1,*) ' DDROWY,DDCOLY,N = ',DDROWY,DDCOLY,N
              STOP
         END IF
C
C
C             ...everything ok => form the m-expansion in desired form,
C                zeroing Y, if wanted.
C
C
         NM = N / M

         IF (LTRG) THEN
C
C
C             ...lower triangle of Y.
C
C
             IF (ZEROY) THEN
                 DO 10 J = 1,N
                 DO 10 I = J,N
                    Y (I+ROWOFF,J) = ZERO
   10            CONTINUE
             END IF

             IF (MJUMP) THEN
C
C
C             ...lower triangle m-expansion in spaced form.
C
C
                 J = 0
                 DO 100 L = 1,M
                    I = ROWOFF
                    DO 110 K = L,M
                       XKL = X (K,L)
                       DO 120 R = 1,NM
                          Y (I+R,J+R) = X (K,L)
  120                  CONTINUE
                       I = I + NM
  110               CONTINUE
                    J = J + NM
  100            CONTINUE
             ELSE
C
C
C             ...lower triangle m-expansion in consecutive form.
C
C
                 I = ROWOFF
                 J = 0
                 DO 200 R = 1,NM
                    DO 210 L = 1,M
                    DO 210 K = L,M
                       Y (I+K,J+L) = X (K,L)
  210               CONTINUE
                    I = I + M
                    J = J + M
  200            CONTINUE
             END IF

         ELSE
C
C
C             ...full matrix Y.
C
C
             IF (ZEROY) THEN
                 DO 20 J = 1,N
                 DO 20 I = 1,N
                    Y (I+ROWOFF,J) = ZERO
   20            CONTINUE
             END IF

             IF (MJUMP) THEN
C
C
C             ...full m-expansion in spaced form.
C
C
                 J = 0
                 DO 300 L = 1,M
                    I = ROWOFF
                    DO 310 K = 1,M
                       XKL = X (K,L)
                       DO 320 R = 1,NM
                          Y (I+R,J+R) = X (K,L)
  320                  CONTINUE
                       I = I + NM
  310               CONTINUE
                    J = J + NM
  300            CONTINUE
             ELSE
C
C
C             ...full m-expansion in consecutive form.
C
C
                 I = ROWOFF
                 J = 0
                 DO 400 R = 1,NM
                    DO 410 L = 1,M
                    DO 410 K = 1,M
                       Y (I+K,J+L) = X (K,L)
  410               CONTINUE
                    I = I + M
                    J = J + M
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
