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
         SUBROUTINE  CCO__REMAP_1234_TO_1243
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y1243 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_1243
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y1243.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    IJ,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N123,N124
         INTEGER    NKX,NKY,NKLX,NKLY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y1243  (1:N1*N2*N4*N3)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y1243  (1:N1,1:N2,1:N4,1:N3)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N123 = N12 * N3
         N124 = N12 * N4

         NKX = - N12
         NKY = - N124
         DO 10 K  = 1,N3
            NKX = NKX + N12
            NKY = NKY + N124
            NKLX = NKX - N123
            NKLY = NKY - N12
            DO 20 L  = 1,N4
               NKLX = NKLX + N123
               NKLY = NKLY + N12
               DO 30 IJ = 1,N12
                  Y1243 (IJ+NKLY) = X1234 (IJ+NKLX)
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y1243 (I,J,L,K) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_1324
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y1324 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_1324
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y1324.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N13,N123
         INTEGER    NL,NKLX,NKLY,NJKLX,NJKLY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y1324  (1:N1*N3*N2*N4)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y1324  (1:N1,1:N3,1:N2,1:N4)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N13 = N1 * N3
         N123 = N12 * N3

         NL = - N123
         DO 10 L = 1,N4
            NL = NL + N123
            NKLX = NL - N12
            NKLY = NL - N1
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N1
               NJKLX = NKLX - N1
               NJKLY = NKLY - N13
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NJKLY = NJKLY + N13
                  DO 40 I = 1,N1
                     Y1324 (I+NJKLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y1324 (I,K,J,L) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_1423
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y1423 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_1423
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y1423.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N14,N123,N124
         INTEGER    NLX,NLY,NKLX,NKLY,NJKLX,NJKLY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y1423  (1:N1*N4*N2*N3)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y1423  (1:N1,1:N4,1:N2,1:N3)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N14 = N1 * N4
         N123 = N12 * N3
         N124 = N12 * N4

         NLX = - N123
         NLY = - N1
         DO 10 L = 1,N4
            NLX = NLX + N123
            NLY = NLY + N1
            NKLX = NLX - N12
            NKLY = NLY - N124
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N124
               NJKLX = NKLX - N1
               NJKLY = NKLY - N14
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NJKLY = NJKLY + N14
                  DO 40 I = 1,N1
                     Y1423 (I+NJKLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y1423 (I,L,J,K) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_1432
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y1432 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_1432
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y1432.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N14,N123,N134
         INTEGER    NLX,NLY,NKLX,NKLY,NJKLX,NJKLY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y1432  (1:N1*N4*N3*N2)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y1432  (1:N1,1:N4,1:N3,1:N2)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N14 = N1 * N4
         N123 = N12 * N3
         N134 = N14 * N3

         NLX = - N123
         NLY = - N1
         DO 10 L = 1,N4
            NLX = NLX + N123
            NLY = NLY + N1
            NKLX = NLX - N12
            NKLY = NLY - N14
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N14
               NJKLX = NKLX - N1
               NJKLY = NKLY - N134
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NJKLY = NJKLY + N134
                  DO 40 I = 1,N1
                     Y1432 (I+NJKLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y1432 (I,L,K,J) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_2134
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y2134 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_2134
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y2134.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,KL
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N34
         INTEGER    NKL,NJKLX,NIKLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y2134  (1:N2*N1*N3*N4)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y2134  (1:N2,1:N1,1:N3,1:N4)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N34 = N3 * N4

         NKL = - N12
         DO 10 KL = 1,N34
            NKL = NKL + N12
            NJKLX = NKL - N1
            NSAVE = NKL - N2
            DO 20 J = 1,N2
               NJKLX = NJKLX + N1
               NIKLY = NSAVE
               DO 30 I = 1,N1
                  NIKLY = NIKLY + N2
                  Y2134 (J+NIKLY) = X1234 (I+NJKLX)
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y2134 (J,I,K,L) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_2143
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y2143 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_2143
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y2143.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N123,N124
         INTEGER    NLX,NLY,NKLX,NKLY,NJKLX,NIKLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y2143  (1:N2*N1*N4*N3)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y2143  (1:N2,1:N1,1:N4,1:N3)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N123 = N12 * N3
         N124 = N12 * N4

         NLX = - N123
         NLY = - N12
         DO 10 L = 1,N4
            NLX = NLX + N123
            NLY = NLY + N12
            NKLX = NLX - N12
            NKLY = NLY - N124
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N124
               NJKLX = NKLX - N1
               NSAVE = NKLY - N2
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NIKLY = NSAVE
                  DO 40 I = 1,N1
                     NIKLY = NIKLY + N2
                     Y2143 (J+NIKLY) = X1234 (I+NJKLX)
   40            CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y2143 (J,I,L,K) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_2314
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y2314 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_2314
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y2314.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N23,N123
         INTEGER    NL,NKLX,NKLY,NJKLX,NIKLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y2314  (1:N2*N3*N1*N4)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y2314  (1:N2,1:N3,1:N1,1:N4)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N23 = N2 * N3
         N123 = N12 * N3

         NL = - N123
         DO 10 L = 1,N4
            NL = NL + N123
            NKLX = NL - N12
            NKLY = NL - N2
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N2
               NJKLX = NKLX - N1
               NSAVE = NKLY - N23
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NIKLY = NSAVE
                  DO 40 I = 1,N1
                     NIKLY = NIKLY + N23
                     Y2314 (J+NIKLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y2314 (J,K,I,L) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_2341
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y2341 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_2341
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y2341.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N23,N123,N234
         INTEGER    NLX,NLY,NKLX,NKLY,NJKLX,NIKLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y2341  (1:N2*N3*N4*N1)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y2341  (1:N2,1:N3,1:N4,1:N1)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N23 = N2 * N3
         N123 = N12 * N3
         N234 = N23 * N4

         NLX = - N123
         NLY = - N23
         DO 10 L = 1,N4
            NLX = NLX + N123
            NLY = NLY + N23
            NKLX = NLX - N12
            NKLY = NLY - N2
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N2
               NJKLX = NKLX - N1
               NSAVE = NKLY - N234
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NIKLY = NSAVE
                  DO 40 I = 1,N1
                     NIKLY = NIKLY + N234
                     Y2341 (J+NIKLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y2341 (J,K,L,I) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_2413
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y2413 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_2413
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y2413.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N24,N123,N124
         INTEGER    NLX,NLY,NKLX,NKLY,NJKLX,NIKLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y2413  (1:N2*N4*N1*N3)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y2413  (1:N2,1:N4,1:N1,1:N3)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N24 = N2 * N4
         N123 = N12 * N3
         N124 = N12 * N4

         NLX = - N123
         NLY = - N1
         DO 10 L = 1,N4
            NLX = NLX + N123
            NLY = NLY + N1
            NKLX = NLX - N12
            NKLY = NLY - N124
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NKLY = NKLY + N124
               NJKLX = NKLX - N1
               NSAVE = NKLY - N24
               DO 30 J = 1,N2
                  NJKLX = NJKLX + N1
                  NIKLY = NSAVE
                  DO 40 I = 1,N1
                     NIKLY = NIKLY + N24
                     Y2413 (J+NIKLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y2413 (J,L,I,K) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_3124
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y3124 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_3124
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y3124.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N13,N123
         INTEGER    NL,NKLX,NJLY,NJKLX,NIJLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y3124  (1:N3*N1*N2*N4)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y3124  (1:N3,1:N1,1:N2,1:N4)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N13 = N1 * N3
         N123 = N12 * N3

         NL = - N123
         DO 10 L = 1,N4
            NL = NL + N123
            NKLX = NL - N12
            NSAVE = NL - N13
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NJLY = NSAVE
               NJKLX = NKLX - N1
               DO 30 J = 1,N2
                  NJLY = NJLY + N13
                  NJKLX = NJKLX + N1
                  NIJLY = NJLY - N3
                  DO 40 I = 1,N1
                     NIJLY = NIJLY + N3
                     Y3124 (K+NIJLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y3124 (K,I,J,L) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_3214
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y3214 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_3214
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y3214.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N23,N123
         INTEGER    NL,NKLX,NJLY,NJKLX,NIJLY
         INTEGER    NSAVE

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y3214  (1:N3*N2*N1*N4)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y3214  (1:N3,1:N2,1:N1,1:N4)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N23 = N2 * N3
         N123 = N12 * N3

         NL = - N123
         DO 10 L = 1,N4
            NL = NL + N123
            NKLX = NL - N12
            NSAVE = NL - N3
            DO 20 K = 1,N3
               NKLX = NKLX + N12
               NJLY = NSAVE
               NJKLX = NKLX - N1
               DO 30 J = 1,N2
                  NJLY = NJLY + N3
                  NJKLX = NJKLX + N1
                  NIJLY = NJLY - N23
                  DO 40 I = 1,N1
                     NIJLY = NIJLY + N23
                     Y3214 (K+NIJLY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y3214 (K,J,I,L) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
cYAU|         SUBROUTINE  CCO__REMAP_1234_TO_3412
cYAU|     +
cYAU|     +                    ( N1,N2,N3,N4,
cYAU|     +                      X1234,
cYAU|     +
cYAU|     +                             Y3412 )
cYAU|     +
cYAU|C------------------------------------------------------------------------
cYAU|C  OPERATION   : CCO__REMAP_1234_TO_3412
cYAU|C  MODULE      : Coupled Cluster Outcore
cYAU|C  MODULE-ID   : CCO
cYAU|C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
cYAU|C                to a "4-dimensional" array Y3412.
cYAU|C
cYAU|C  AUTHOR      : Norbert Flocke
cYAU|C------------------------------------------------------------------------
cYAU|C
cYAU|C
cYAU|C             ...include files and declare variables.
cYAU|C
cYAU|C
cYAU|         INTEGER    I,J,K,L
cYAU|         INTEGER    N1,N2,N3,N4
cYAU|
cYAU|         DOUBLE PRECISION   X1234  (1:N1*N2,1:N3*N4)
cYAU|         DOUBLE PRECISION   Y3412  (1:N3*N4,1:N1*N2)
cYAU|
cYAU|C         INTEGER    I,J,K,L
cYAU|C         INTEGER    N1,N2,N3,N4
cYAU|C
cYAU|C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
cYAU|C         DOUBLE PRECISION   Y3412  (1:N3,1:N4,1:N1,1:N2)
cYAU|C
cYAU|C
cYAU|C------------------------------------------------------------------------
cYAU|C
cYAU|C
cYAU|C             ...perform the remapping.
cYAU|C
cYAU|C
cYAU|         CALL    MAT__C_EQ_A_TRANSPOSED_FLOAT
cYAU|     +
cYAU|     +                ( N1*N2,N3*N4,
cYAU|     +                  N3*N4,N1*N2,
cYAU|     +                  N1*N2,N3*N4,
cYAU|     +                  X1234,
cYAU|     +
cYAU|     +                          Y3412 )
cYAU|     +
cYAU|C
cYAU|C
cYAU|C             ...the dirty way for checking.
cYAU|C
cYAU|C
cYAU|C         DO 10 I = 1,N1
cYAU|C         DO 10 J = 1,N2
cYAU|C         DO 10 K = 1,N3
cYAU|C         DO 10 L = 1,N4
cYAU|C            Y3412 (K,L,I,J) = X1234 (I,J,K,L)
cYAU|C   10    CONTINUE
cYAU|C
cYAU|C
cYAU|C             ...ready!
cYAU|C
cYAU|C
cYAU|         RETURN
cYAU|         END
         SUBROUTINE  CCO__REMAP_1234_TO_4213
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y4213 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_4213
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y4213.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N24,N123,N124
         INTEGER    NLX,NKY,NKLX,NJKY,NJKLX,NIJKY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y4213  (1:N4*N2*N1*N3)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y4213  (1:N4,1:N2,1:N1,1:N3)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N24 = N2 * N4
         N123 = N12 * N3
         N124 = N12 * N4

         NLX = - N123
         DO 10 L = 1,N4
            NLX = NLX + N123
            NKY = - N124
            NKLX = NLX - N12
            DO 20 K = 1,N3
               NKY = NKY + N124
               NKLX = NKLX + N12
               NJKY = NKY - N4
               NJKLX = NKLX -  N1
               DO 30 J = 1,N2
                  NJKY = NJKY + N4
                  NJKLX = NJKLX + N1
                  NIJKY = NJKY - N24
                  DO 40 I = 1,N1
                     NIJKY = NIJKY + N24
                     Y4213 (L+NIJKY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y4213 (L,J,I,K) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_4312
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y4312 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_4312
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y4312.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N34,N123,N134
         INTEGER    NLX,NKY,NKLX,NJKY,NJKLX,NIJKY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y4312  (1:N4*N3*N1*N2)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y4312  (1:N4,1:N3,1:N1,1:N2)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N34 = N3 * N4
         N123 = N12 * N3
         N134 = N1 * N34

         NLX = - N123
         DO 10 L = 1,N4
            NLX = NLX + N123
            NKY = - N4
            NKLX = NLX - N12
            DO 20 K = 1,N3
               NKY = NKY + N4
               NKLX = NKLX + N12
               NJKY = NKY - N134
               NJKLX = NKLX -  N1
               DO 30 J = 1,N2
                  NJKY = NJKY + N134
                  NJKLX = NJKLX + N1
                  NIJKY = NJKY - N34
                  DO 40 I = 1,N1
                     NIJKY = NIJKY + N34
                     Y4312 (L+NIJKY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y4312 (L,K,I,J) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
         SUBROUTINE  CCO__REMAP_1234_TO_4321
     +
     +                    ( N1,N2,N3,N4,
     +                      X1234,
     +
     +                             Y4321 )
     +
C------------------------------------------------------------------------
C  OPERATION   : CCO__REMAP_1234_TO_4321
C  MODULE      : Coupled Cluster Outcore
C  MODULE-ID   : CCO
C  DESCRIPTION : This routine remaps a "4-dimensional" array X1234
C                to a "4-dimensional" array Y4321.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    I,J,K,L
         INTEGER    N1,N2,N3,N4
         INTEGER    N12,N34,N123,N234
         INTEGER    NLX,NKY,NKLX,NJKY,NJKLX,NIJKY

         DOUBLE PRECISION   X1234  (1:N1*N2*N3*N4)
         DOUBLE PRECISION   Y4321  (1:N4*N3*N2*N1)

C         INTEGER    I,J,K,L
C         INTEGER    N1,N2,N3,N4
C
C         DOUBLE PRECISION   X1234  (1:N1,1:N2,1:N3,1:N4)
C         DOUBLE PRECISION   Y4321  (1:N4,1:N3,1:N2,1:N1)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the remapping.
C
C
         N12 = N1 * N2
         N34 = N3 * N4
         N123 = N12 * N3
         N234 = N2 * N34

         NLX = - N123
         DO 10 L = 1,N4
            NLX = NLX + N123
            NKY = - N4
            NKLX = NLX - N12
            DO 20 K = 1,N3
               NKY = NKY + N4
               NKLX = NKLX + N12
               NJKY = NKY - N34
               NJKLX = NKLX -  N1
               DO 30 J = 1,N2
                  NJKY = NJKY + N34
                  NJKLX = NJKLX + N1
                  NIJKY = NJKY - N234
                  DO 40 I = 1,N1
                     NIJKY = NIJKY + N234
                     Y4321 (L+NIJKY) = X1234 (I+NJKLX)
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
C
C             ...the dirty way for checking.
C
C
C         DO 10 I = 1,N1
C         DO 10 J = 1,N2
C         DO 10 K = 1,N3
C         DO 10 L = 1,N4
C            Y4321 (L,K,J,I) = X1234 (I,J,K,L)
C   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
