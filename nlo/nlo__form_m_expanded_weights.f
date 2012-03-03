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
         SUBROUTINE  NLO__FORM_M_EXPANDED_WEIGHTS
     +
     +                    ( DDVECX,
     +                      DDVECW,
     +                      NLM,NL,
     +                      MJUMP,
     +                      X,
     +
     +                              W )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_M_EXPANDED_WEIGHTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms the m-expanded weight vector W
C                of size NLM from an input m-averaged weight vector X
C                of size NL. The following example for m=x,y,z makes
C                this clear:
C
C
C                          x,y,z expanded          x y z x y z 
C                  a b   ------------------>       a a a b b b
C                         6 element vector
C
C
C                If the x,y,z basis functions are ordered such that
C                each x,y,z are consecutive, then we would have the
C                following W vector expansion:
C
C
C                          x,y,z expanded          x x y y z z 
C                  a b   ------------------>       a b a b a b
C                         6 element vector
C
C
C
C                Thus the m-expansion is only defined, if NLM is
C                divisible by NL. In case NLM is not divisible by NL
C                the routine stops with an error message.
C
C                  Input:
C
C                           DDVECz  =  declared dimensions of vectors
C                                      z = X and W
C                              NLM  =  # of elements of vector W
C                               NL  =  # of elements of vector X
C                            MJUMP  =  is .true., if the m values in
C                                      the m-space are ordered such
C                                      that the same m values are
C                                      separated. Hence this needs
C                                      evaluation of the spaced
C                                      weight expansion. If false, the
C                                      consecutive weight expansion is
C                                      performed
C                                X  =  original m-averaged weight
C                                      vector size NL
C
C                  Output:
C
C                                W  =  expanded weight vector of size
C                                      NLM
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

         INTEGER     DDVECX,DDVECW
         INTEGER     I,K,M
         INTEGER     NL,NM,NLM

         DOUBLE PRECISION  XK

         DOUBLE PRECISION  W (1:DDVECW)
         DOUBLE PRECISION  X (1:DDVECX)
C
C
C------------------------------------------------------------------------
C
C
C             ...check, if NLM is divisible by NL and the dimensions
C                of X and W.
C
C
         IF  (MOD (NLM,NL).NE.0)  THEN
              WRITE (1,*) ' Cannot form m-expanded weight vector! '
              WRITE (1,*) ' nlo__form_m_expanded_weights '
              WRITE (1,*) ' NLM,NL = ',NLM,NL
              WRITE (*,*) ' Cannot form m-expanded weight vector! '
              WRITE (*,*) ' nlo__form_m_expanded_weights '
              WRITE (*,*) ' NLM,NL = ',NLM,NL
              STOP
         END IF

         IF  (NL.GT.DDVECX)  THEN
              WRITE (1,*) ' Dimensions of vector X too small: '
              WRITE (1,*) ' nlo__form_m_expanded_weights '
              WRITE (1,*) ' DDVECX,NL = ',DDVECX,NL
              WRITE (*,*) ' Dimensions of vector X too small: '
              WRITE (*,*) ' nlo__form_m_expanded_weights '
              WRITE (*,*) ' DDVECX,NL = ',DDVECX,NL
              STOP
         END IF

         IF  (NLM.GT.DDVECW)  THEN
              WRITE (1,*) ' Dimensions of vector W too small: '
              WRITE (1,*) ' nlo__form_m_expanded_weights '
              WRITE (1,*) ' DDVECW,NL = ',DDVECW,NL
              WRITE (*,*) ' Dimensions of vector W too small: '
              WRITE (*,*) ' nlo__form_m_expanded_weights '
              WRITE (*,*) ' DDVECW,NL = ',DDVECW,NL
              STOP
         END IF
C
C
C             ...everything ok => form the m-expanded weights
C                in desired form.
C
C
         NM = NLM / NL

         IF (MJUMP) THEN
C
C
C             ...m-expanded weights in spaced form.
C
C
             I = 0
             DO K = 1,NL
                XK = X (K)
                DO M = 1,NM
                   W (I+M) = XK
                END DO
                I = I + NM
             END DO
         ELSE
C
C
C             ...m-expanded weights in consecutive form.
C
C
             DO K = 1,NL
                XK = X (K)
                I = K
                DO M = 1,NM
                   W (I) = XK
                   I = I + NL
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
