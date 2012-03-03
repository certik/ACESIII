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
         SUBROUTINE  NLO__FORM_M_EXPANDED_COEFFS
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWT,DDCOLT,
     +                      DDROWC,DDCOLC,
     +                      N,NLM,NL,
     +                      MJUMP,
     +                      X,T,
     +
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_M_EXPANDED_COEFFS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms the new m-expanded coefficient
C                matrix C of size N x NLM from an input m-expanded
C                coefficient matrix of the same size and and m-averaged
C                transformation matrix T of size NL x NL. This routine
C                thus represents a coefficient matrix transformation:
C
C                                C = X * T
C
C                between X and C involving an expanded NLM x NLM
C                transformation matrix T, which is formally obtained
C                by an m-expansion, as shown in the following example
C                for m=x,y,z:
C
C
C                                                   x y z x y z 
C                                                   -----------
C                                                x |a    |b    |
C                              T matrix          y |  a  |  b  |
C                  a  b     x,y,z expanded       z |    a|    b|
C                         ------------------>       ----------- 
C                  c  d      6 x 6 matrix        x |c    |d    |
C                                                y |  c  |  d  |
C                                                z |    c|    d|
C                                                   -----------
C
C
C                If the x,y,z basis functions are ordered such that
C                each x,y,z are consecutive, then we would have the
C                following T matrix expansion:
C
C
C                                                   x x y y z z 
C                                                   -----------
C                                                x |a b|   |   |
C                                                x |c d|   |   |
C                              T matrix             -----------
C                  a  b     x,y,z expanded       y |   |a b|   |
C                         ------------------>    y |   |c d|   |
C                  c  d      6 x 6 matrix           -----------
C                                                z |   |   |a b|
C                                                z |   |   |c d|
C                                                   -----------
C
C
C                Thus the m-expansion is only defined, if NLM is
C                divisible by NL. The coefficient transformation
C                can thus proceed via two possible ways, depending
C                on the m-order of the basis set. For each C matrix
C                coefficient element:
C
C
C                 1) Spaced form, taking M subblock diagonal elements
C                    from the T matrix spaced by M:
C
C                                         NL
C                    C (r,i+M*[l-1])  =  sum  X (r,i+M*[k-1]) * T (k,l)
C                                         k
C
C                 2) Consecutive form, taking M subblock elements
C                    from the T matrix:
C
C                                         NL
C                    C (r,l+NL*[i-1]) =  sum  X (r,k+NL*[i-1]) * T (k,l)
C                                         k
C
C                where i=1,M and r=1,N. In case NLM is not divisible
C                by NL, the routine stops with an error message.
C
C                  Input:
C
C                    DDROWz,DDCOLz  =  declared dimensions of matrices
C                                      z = X,T and C
C                                N  =  # of rows of matrices X and C
C                              NLM  =  # of columns of matrices X and C
C                               NL  =  size of m-averaged transformation
C                                      matrix T
C                            MJUMP  =  is .true., if the m values in
C                                      the m-space are ordered such
C                                      that the same m values are
C                                      separated. Hence this needs
C                                      evaluation of the spaced
C                                      transformation. If false, the
C                                      consecutive transformation is
C                                      performed
C                                X  =  original coefficient matrix of
C                                      size N x NLM
C                                T  =  m-averaged transformation matrix
C                                      of size NL x NL
C
C                  Output:
C
C                                C  =  transformed coefficient matrix
C                                      of size N x NLM
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

         LOGICAL     MJUMP

         INTEGER     DDROWX,DDCOLX,DDROWT,DDCOLT,DDROWC,DDCOLC
         INTEGER     I,J,K,L,M,N,R
         INTEGER     NL,NM,NLM

         DOUBLE PRECISION  TKL

         DOUBLE PRECISION  C (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  T (1:DDROWT,1:DDCOLT)
         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)
C
C
C------------------------------------------------------------------------
C
C
C             ...check, if NLM is divisible by NL and the dimensions
C                of C,T and X.
C
C
         IF  (MOD (NLM,NL).NE.0)  THEN
              WRITE (1,*) ' Cannot form m-expanded coeff matrix! '
              WRITE (1,*) ' nlo__form_m_expanded_coeffs '
              WRITE (1,*) ' NLM,NL = ',NLM,NL
              WRITE (*,*) ' Cannot form m-expanded coeff matrix! '
              WRITE (*,*) ' nlo__form_m_expanded_coeffs '
              WRITE (*,*) ' NLM,NL = ',NLM,NL
              STOP
         END IF

         IF  (N.GT.DDROWC .OR. NLM.GT.DDCOLC)  THEN
              WRITE (1,*) ' Dimensions of matrix C too small: '
              WRITE (1,*) ' nlo__form_m_expanded_coeffs '
              WRITE (1,*) ' DDROWC,DDCOLC,N,NLM = ',DDROWC,DDCOLC,N,NLM
              WRITE (*,*) ' Dimensions of matrix C too small: '
              WRITE (*,*) ' nlo__form_m_expanded_coeffs '
              WRITE (*,*) ' DDROWC,DDCOLC,N,NLM = ',DDROWC,DDCOLC,N,NLM
              STOP
         END IF

         IF  (NL.GT.DDROWT .OR. NL.GT.DDCOLT)  THEN
              WRITE (1,*) ' Dimensions of matrix T too small: '
              WRITE (1,*) ' nlo__form_m_expanded_coeffs '
              WRITE (1,*) ' DDROWT,DDCOLT,NL = ',DDROWT,DDCOLT,NL
              WRITE (*,*) ' Dimensions of matrix T too small: '
              WRITE (*,*) ' nlo__form_m_expanded_coeffs '
              WRITE (*,*) ' DDROWT,DDCOLT,NL = ',DDROWT,DDCOLT,NL
              STOP
         END IF

         IF  (N.GT.DDROWX .OR. NLM.GT.DDCOLX)  THEN
              WRITE (1,*) ' Dimensions of matrix X too small: '
              WRITE (1,*) ' nlo__form_m_expanded_coeffs '
              WRITE (1,*) ' DDROWX,DDCOLX,N,NLM = ',DDROWX,DDCOLX,N,NLM
              WRITE (*,*) ' Dimensions of matrix X too small: '
              WRITE (*,*) ' nlo__form_m_expanded_coeffs '
              WRITE (*,*) ' DDROWX,DDCOLX,N,NLM = ',DDROWX,DDCOLX,N,NLM
              STOP
         END IF
C
C
C             ...everything ok => form the m-expanded coefficients
C                in desired form.
C
C
         NM = NLM / NL

         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( DDROWC,DDCOLC,
     +                N,NLM,
     +
     +                        C )
     +
     +
         IF (MJUMP) THEN
C
C
C             ...m-expanded coefficients in spaced form.
C
C
             I = 0
             DO L = 1,NL
                J = 0
                DO K = 1,NL
                   TKL = T (K,L)
                   DO M = 1,NM
                      J = J + 1
                      DO R = 1,N
                         C (R,I+M) = C (R,I+M) + TKL * X (R,J)
                      END DO
                   END DO
                END DO
                I = I + NM
             END DO
         ELSE
C
C
C             ...m-expanded coefficients in consecutive form.
C
C
             DO L = 1,NL
             DO K = 1,NL
                TKL = T (K,L)
                I = L
                J = K
                DO M = 1,NM
                   DO R = 1,N
                      C (R,I) = C (R,I) + TKL * X (R,J)
                   END DO
                   I = I + NL
                   J = J + NL
                END DO
             END DO
             END DO

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
