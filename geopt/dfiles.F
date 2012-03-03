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

C THIS SUBROUTINE DETERMINES THE INFORMATION ABOUT THE
C LOCATION OF VARIOUS FILES (DEFAULT IS OF COURSE
C THE WORKING DIRECTORY) IN AN ACES2 RUN. THIS OPTION
C IS PARTICULARLY USEFUL FOR WORKSTATIONS, BUT HAS
C SOME BENEFITS ON SUPER COMPUTERS WITH TMPDIR DISKS.
C
C IF A FILE IS INTENDED NOT TO BE IN THE WORKING
C DIRECTORY, THE FULL PATH NAME CAN BE GIVEN IN THE
C FOLLOWING WAY IN THE ZMAT FILE PRIOR TO THE TITLE
C
C  %MOINTS=/work/username/test1/MOINTS
C
C WHICH MEANS THAT MOINTS IS NOW IN THE DIRECTORY
C /work/username/test1. THERE IS ALSO THE POSSIBILITY
C TO RENAME FILES BUT THAT SHOULD BE USED WITH
C COME CAUTION, IN ORDER NOT TO MESS UP THINGS.
C
C JG GOETEBORG SEPT/92

c Yau : 01/2003
c This routine was altered to create FILES only if the user specifies
c a file mapping. Oh, and it was rewritten for more robust parsing.

c Yau : 02/2003
c A special file directive SAVEDIR was created to allow the user to
c specify a location (absolute or relative to the ACES working directory)
c in which precious files will be stored for [finite displacement] restarts.

#include "flags.h"

      SUBROUTINE DFILES(IGNORE)
      IMPLICIT NONE

#include "io_units.par"
#include "linelen.par"

      LOGICAL IGNORE

      character*1 achar
      logical leq
      external leq

      integer i, j
      integer iErr, iFlags(100), iFlags2(500)
      integer iStart1, iEnd1, iStart2, iEnd2, LINENO
      logical bVerbose, bExist, bNoTitle, bSaveDir, bBackup, bRestore
      integer iJOBARC
      CHARACTER*(linelen) A, szSaveDir, szJOBARC
      character*1 czNull, czSpace, czTab

c ----------------------------------------------------------------------

      czNull  = achar(0)
      czTab   = achar(9)
      czSpace = achar(32)

c   o assume SAVEDIR=SAVEDIR (bSaveDir is the user-override flag)
      szSaveDir = 'SAVEDIR'//czNull
      bSaveDir  = .false.

c   o remove FILES (it may be contaminated on a restart)
c      INQUIRE(FILE='FILES',EXIST=bExist)
       bexist = .false.
      if (bExist) then
         call f_remove('FILES')
         bVerbose = .false.
      else
         bVerbose = .true.
      end if

c   o open ZMAT
c      INQUIRE(FILE=ZFIL,EXIST=bExist)
      bexist = .true.
      IF (bExist) THEN
         OPEN(UNIT=LUZ,FILE=ZFIL,FORM='FORMATTED',STATUS='OLD',ERR=666)
         REWIND(LUZ)
      ELSE
         WRITE(*,*) '@DFILES: ZMAT file does not exist.'
         CALL ERREX
      END IF

c   o scan the header until the job title is found
c     (i.e., first non-white character is not '%' or '#')
      bNoTitle = .true.
      bExist   = .false.
      LINENO   = 0
      do while (bNoTitle)

         read(LUZ,'(a)',end=9999,err=666) a
         LINENO = LINENO + 1

c      o scan for the first character
         i = 1
         do while ((i.lt.linelen).and.
     &             ((a(i:i).eq.czSpace).or.(a(i:i).eq.czTab))
     &            )
            i = i + 1
         end do

c      o branch on the first character
         if (a(i:i).ne.'%') then

            bNoTitle = (a(i:i).eq.czSpace).or.
     &                 (a(i:i).eq.czTab).or.
     &                 (a(i:i).eq.'#')

c        else if ([file directive])
         else

c         o scan for a '#' comment delim
            j = i + 1
            do while ((j.lt.linelen).and.(a(j:j).ne.'#'))
               j = j + 1
            end do
            if (a(j:j).eq.'#') j = j - 1

c         o skip white space
            i = i + 1
            do while ((i.lt.j).and.
     &                ((a(i:i).eq.czSpace).or.(a(i:i).eq.czTab))
     &               )
               i = i + 1
            end do

c         o find white space or '=' delim
            iStart1 = i
            do while ((i.lt.j).and.
     &                ((a(i:i).ne.czSpace).and.(a(i:i).ne.czTab).and.
     &                 (a(i:i).ne.'='))
     &               )
               i = i + 1
            end do
            iEnd1 = i - 1

c         o skip white space and '=' delim
            do while ((i.lt.j).and.
     &                ((a(i:i).eq.czSpace).or.(a(i:i).eq.czTab).or.
     &                 (a(i:i).eq.'='))
     &               )
               i = i + 1
            end do

c         o find white space up to j
            iStart2 = i
            do while ((i.lt.j).and.
     &                ((a(i:i).ne.czSpace).and.(a(i:i).ne.czTab))
     &               )
               i = i + 1
            end do
            iEnd2 = i - 1

c         o verify integrity of substrings
            if ((iStart1.eq.iEnd1+1).or.
     &          (iStart2.eq.iEnd2+1)    ) then
               print *, '@DFILES: DEFECT IN FILE DIRECTIVE ON LINE ',
     &                  LINENO, ' OF ZMAT FILE'
               call errex
            end if

            if (leq(a(iStart1:iEnd1),'SAVEDIR')) then
               szSaveDir = a(iStart2:iEnd2)//czNull
               bSaveDir  = .true.
            else
c            o create FILES
               if (.not.bExist) then
                  OPEN(UNIT=LuFiles,FILE='FILES',FORM='FORMATTED',
     &                 STATUS='NEW')
                  REWIND(LuFiles)
                  bExist = .true.
                  if (bVerbose) then
                     print '()'
                     print *, '@DFILES: FILE MAPPING FROM ZMAT HEADER:'
                  end if
               end if
c            o write to FILES
               WRITE(LuFiles,'(a,x,a)')A(ISTART1:IEND1),A(ISTART2:IEND2)
               if (bVerbose) then
                  print *, A(ISTART1:IEND1), ' -> ', A(ISTART2:IEND2)
               end if
            end if

c        end if ([not file directive])
         end if

c     end do while (bNoTitle)
      end do

c   o add SAVEDIR to FILES only when the user overrides it
      if (bSaveDir) then
         if (.not.bExist) then
            OPEN(UNIT=LuFiles,FILE='FILES',FORM='FORMATTED',
     &           STATUS='NEW')
            REWIND(LuFiles)
            bExist = .true.
            if (bVerbose) then
               print '()'
               print *, '@DFILES: Archive directory:'
            end if
         end if
         i = 1
         do while (szSaveDir(i:i).ne.czNull)
            i = i + 1
         end do
         i = i - 1
         write(LuFiles,'(a,a)')'SAVEDIR ',szSaveDir(1:i)
         if (bVerbose) then
            print *, 'SAVEDIR == ', szSaveDir(1:i)
         end if
      end if

c   o close FILES
      if (bExist) then
         if (bVerbose) write(*,*)
         WRITE(LuFiles,*)
         CLOSE(LuFiles)
      end if

c   o close ZMAT
      CLOSE(LUZ)

c   o assume the current file set is good
      call gfname('JOBARC',szJOBARC,iJOBARC)
      inquire(file=szJOBARC(1:iJOBARC),exist=ignore)

c   o backup or restore the precious files
c      if (bSaveDir) then
      if (.true.) then
         bRestore = .true.
         bBackup  = .false.
         if (ignore) then
c         o limit coarse backups to clean geometry optimizations and vib freqs
            call aces_ja_init
            call igetrec(1,'JOBARC','IFLAGS',100,iFlags)
            call igetrec(1,'JOBARC','DIRTYFLG',1,i)
            if (i.eq.0) then
c            o JOBARC is clean
               bRestore = .false.
               if (iFlags(h_IFLAGS_restart).ne.0) then
                  call igetrec(1,'JOBARC','IFLAGS2',500,iFlags2)
                  bBackup = (iFlags2(h_IFLAGS2_geom_opt).ne.0).or.
     &                      (iFlags(h_IFLAGS_vib).eq.3)
               end if
            else
c            o JOBARC is dirty
               bRestore = (iFlags(h_IFLAGS_restart).ne.0)
            end if
            call aces_ja_fin
c        end if (JOBARC exists)
         end if
         if (bRestore) then
            i = 1
            do while (szSaveDir(i:i).ne.czNull)
               i = i + 1
            end do
            i = i - 1
            if (ignore) then
               print *, '@DFILES: The current file set is contaminated.'
               print *, '         Joda will attempt to recover from ',
     &                  szSaveDir(1:i)
            end if
            call coarse_restore(szSaveDir,iErr)
            if (iErr.ne.0) call errex()
            if (.not.ignore) then
               inquire(file=szJOBARC(1:iJOBARC),exist=ignore)
               if (ignore) then
                  print *, ' '
                  print *, '@DFILES: The ACES II file set was missing,',
     &                     ' but joda was able to recover'
                  print *, '         a checkpoint from ',szSaveDir(1:i)
                  print *, '         If this calculation should not',
     &                     ' have been restarted, then remove'
                  print *, '         the CURRENT and/or OLD',
     &                     ' directories in the checkpoint directory.'
                  print *, ' '
               end if
            end if
            bBackup = .false.
         end if
         if (bBackup) then
            call coarse_backup(szSaveDir,iErr)
            if (iErr.ne.0) call errex()
         end if
c     end if (bSaveDir)
      end if

      RETURN

 666  print '()'
      print *, '@DFILES: An error was encountered while attempting to ',
     &         'open or read ZMAT.'
      call errex

 9999 print '()'
      print *, '@DFILES: The ZMAT file is empty.'
      call errex

      END

