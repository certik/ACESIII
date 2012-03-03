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
         SUBROUTINE  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                    ( STAGE,
     +                      BONDSIZE,
     +                      MAXOCC,
     +                      NHCEN,
     +                      NWSTEP,
     +                      WMAX,WMIN,WSTEP,
     +                      WBOND,WSTAR )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NHO_SEARCH_INFO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints information related to the NHO
C                generation process, depending on the generation stage.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         CHARACTER*1    DASH
         CHARACTER*(*)  STAGE

         INTEGER     BONDSIZE
         INTEGER     NHCEN
         INTEGER     NSTEP
         INTEGER     SIZE

         INTEGER     NWSTEP (1:BONDSIZE)

         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  WBMAX,WBMIN,WSMAX,WSMIN
         DOUBLE PRECISION  WMAX,WMIN
         DOUBLE PRECISION  WSTEP
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  WBOND (1:BONDSIZE)
         DOUBLE PRECISION  WSTAR (1:BONDSIZE)

         DATA  DASH  /'-'/
         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...print out info according to stage.
C
C
         IF (STAGE.EQ.'start') THEN

             WRITE (*,9100) ' --- Starting NHO search --- '
             WRITE (1,9100) ' --- Starting NHO search --- '

             WRITE (*,9110) '# cen','wbond range','wstar range'
             WRITE (1,9110) '# cen','wbond range','wstar range'

             DO 10 SIZE = 1,BONDSIZE
                NSTEP = NWSTEP (SIZE)
                WBMAX = MAXOCC
                WBMIN = WMAX - (NSTEP-1)*WSTEP
                WSMAX = WMIN + (NSTEP-1)*WSTEP
                WSMIN = ZERO
                WRITE (*,9120) SIZE,WBMIN,DASH,WBMAX,WSMAX,DASH,WSMIN
                WRITE (1,9120) SIZE,WBMIN,DASH,WBMAX,WSMAX,DASH,WSMIN
   10        CONTINUE

         ELSE IF (STAGE.EQ.'failure') THEN

             WRITE (*,9200) ' Failure at bondsize = ',NHCEN
             WRITE (1,9200) ' Failure at bondsize = ',NHCEN

         ELSE IF (STAGE.EQ.'complete') THEN

             WRITE (*,9300) ' Trying to complete NHO space ... '
             WRITE (1,9300) ' Trying to complete NHO space ... '

         ELSE IF (STAGE.EQ.'failure complete') THEN

             WRITE (*,9400) ' Failure completing NHO space! '
             WRITE (1,9400) ' Failure completing NHO space! '

         ELSE IF (STAGE.EQ.'more') THEN

             WRITE (*,9500) ' Searching for remaining NHOs ... '
             WRITE (1,9500) ' Searching for remaining NHOs ... '

             WRITE (*,9510) '# cen','wbond range','wstar range'
             WRITE (1,9510) '# cen','wbond range','wstar range'

             DO 50 SIZE = 1,BONDSIZE
                WBMAX = WBOND (SIZE)
                WBMIN = WBOND (SIZE) - WSTEP
                WSMAX = WSTAR (SIZE) + WSTEP
                WSMIN = WSTAR (SIZE)
                WRITE (*,9520) SIZE,WBMIN,DASH,WBMAX,WSMAX,DASH,WSMIN
                WRITE (1,9520) SIZE,WBMIN,DASH,WBMAX,WSMAX,DASH,WSMIN
   50        CONTINUE

         ELSE IF (STAGE.EQ.'success') THEN

             WRITE (*,9600) ' NHO space generation successful ! '
             WRITE (1,9600) ' NHO space generation successful ! '
             WRITE (*,9610) ' ----- Final weight limits ----- '
             WRITE (1,9610) ' ----- Final weight limits ----- '

             WRITE (*,9620) '# cen','wbond range','wstar range'
             WRITE (1,9620) '# cen','wbond range','wstar range'

             DO 60 SIZE = 1,BONDSIZE
                WBMAX = MAXOCC
                WBMIN = WBOND (SIZE)
                WSMAX = WSTAR (SIZE)
                WSMIN = ZERO
                WRITE (*,9630) SIZE,WBMIN,DASH,WBMAX,WSMAX,DASH,WSMIN
                WRITE (1,9630) SIZE,WBMIN,DASH,WBMAX,WSMAX,DASH,WSMIN
   60        CONTINUE

         END IF
C
C
C             ...formats used.
C
C
 9100    FORMAT (A29)
 9110    FORMAT (1X,A5,1X, 2(1X,A11,1X))
 9120    FORMAT (2X,I2,3X, 2(1X,F4.2,1X,A1,1X,F4.2,1X))

 9200    FORMAT (A23,I3)
 9300    FORMAT (A34)
 9400    FORMAT (A31)

 9500    FORMAT (A34)
 9510    FORMAT (1X,A5,1X, 2(1X,A11,1X))
 9520    FORMAT (2X,I2,3X, 2(1X,F4.2,1X,A1,1X,F4.2,1X))

 9600    FORMAT (A35)
 9610    FORMAT (A33)
 9620    FORMAT (1X,A5,1X, 2(1X,A11,1X))
 9630    FORMAT (2X,I2,3X, 2(1X,F4.2,1X,A1,1X,F4.2,1X))
C
C
C             ...ready!
C
C
         RETURN
         END
