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

c The purpose of this routine is to set the number of CACHE_RECS and
c the FILE_RECSIZE based on the MEMORY_SIZE.  The goals are to use 5%
c of ICORE for the I/O cache and to have between 32 and 128
c CACHE_RECS.  We will start with a FILE_RECSIZE of 64 KB and vary it
c to match these goals.

c SG 11/15/97

#define _MAX_CACHE_SLOTS 128 /* from acescore/include/cache.com */

      SUBROUTINE FIGIO(TOTMEM,RECNUM,RECSIZ)
      IMPLICIT NONE

      INTEGER TOTMEM, RECNUM, RECSIZ

c This is the default record size, in bytes. It should be some (large)
c multiple of 512.
cYAU
c   If records will be used will with common indices (e.g., a(i)=b(i)+c(i)),
c then these should NOT be large powers of two like 512, 1024, etc. They
c should certainly be aligned on 128 Byte boundaries, but definitely not
c large powers of two. Search online for 'cache optimization pad'.
      INTEGER DEFSIZ
      PARAMETER (DEFSIZ = 64 * 1024)

      INTEGER MEMUSE
      LOGICAL ADJSIZ, ADJNUM, bTmp1, bTmp2

#include "machsp.com"

c Figure out five percent of MEMORY_SIZE.
      MEMUSE = (TOTMEM / 20)
      ADJSIZ = (RECSIZ .LE. 0)
      ADJNUM = (RECNUM .LE. 0)

      IF ((ADJSIZ .OR. ADJNUM) .AND. (MEMUSE .LT. 32)) THEN
         PRINT *, '@FIGIO: MEMORY_SIZE is too small.'
         CALL ERREX
      END IF

c Calculate the number of records with the default record size
      IF (ADJSIZ) RECSIZ = DEFSIZ / IINTLN

      IF (ADJNUM) THEN

         bTmp1 = .true.
         do while (bTmp1)
            RECNUM = MEMUSE / RECSIZ
            IF (RECNUM .LT. 32) THEN
               IF (ADJSIZ) THEN
c               o Reduce record size
                  RECSIZ = RECSIZ / 2
               ELSE
                  RECNUM = 32
                  bTmp1  = .false.
               END IF
            ELSE
               IF (RECNUM .GT. _MAX_CACHE_SLOTS) THEN
                  IF (ADJSIZ) THEN
c                  o Increase record size
                     RECSIZ = RECSIZ * 2
                  ELSE
                     RECNUM = _MAX_CACHE_SLOTS
                     bTmp1  = .false.
                  END IF
               ELSE
                  bTmp1 = .false.
               END IF
            END IF
c        end do while (.bTmp1)
         end do

      ELSE

         IF (ADJSIZ) THEN
            bTmp2 = .true.
            do while (bTmp2)
               IF ((RECNUM * RECSIZ) .GT. MEMUSE) THEN
c               o Reduce record size
                  RECSIZ = RECSIZ / 2
               ELSE
                  IF ((RECNUM * RECSIZ) .LT. MEMUSE / 2) THEN
c                  o Increase record size
                     RECSIZ = RECSIZ * 2
                  ELSE
c                  o quit
                     bTmp2 = .false.
                  END IF
               END IF
c           end do while (bTmp2)
            end do
c        END IF (ADJSIZ)
         END IF

c     END IF (ADJNUM)
      END IF

c   o override values that are too large
      if (recsiz.gt.262144) then
c      o reduce to 1MB (32-bit) or 2MB (64-bit)
         if (.not.adjsiz) then
            print *, '@FIGIO: reducing records to ',3-iintfp,' MB'
         end if
         if (adjnum) then
            recnum = totmem/5242880
         else
c         o This is a bad convention. The user should be allowed only to
c           specify the size of the cache in units of bytes or as a percentage
c           of MEMSIZE.
            print *, '@FIGIO: increasing number of cache slots'
            recnum = recsiz*recnum/262144
         end if
         recsiz = 262144
         recnum = min(recnum,_MAX_CACHE_SLOTS)
      end if

      RETURN
      END

