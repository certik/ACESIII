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
      subroutine joda_main()
      implicit none

      integer pass1, i, istat
      logical fd, geomopt, num_grad

      integer iprt, ierr
      character*80 BasNam

#include "jodaflags.com"
#include "flags.h"
#include "mxatms.par"
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     1     PARNAM(3*MXATMS)

#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      INTEGER  NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT

      integer genby(mxatms), scratch(mxatms)
 
c   o initialize the runtime environment
      call aces_init_rte

c   o initialize the job archive subsystem
c      call aces_ja_init

c      o this is the first joda run
c         call iputrec(1,'JOBARC','FIRSTRUN',1,1)
c         call iputrec(1,'JOBARC','DIRTYFLG',1,0)
c         call iputrec(1,'JOBARC','JODADONE',1,0)
c         call iputrec(1,'JOBARC','FNDFDONE',1,1)

c         iprt = 1 
c         call gtflgs(iprt, ierr, BasNam)   ! decipher flags from ZMAT
c         if (iflags(h_IOPPAR_coordinates) .eq. 1) then 
c            call getxyz()   ! read ZMAT Cartesian geometry into common 
c         else
c            call read_zmat_internal_coords()  ! read ICs from ZMAT.
c         endif
c
c         call mkvmol(q, ' ', natoms, nunique, zsym, iatnum, 
c     *               GenBy, scratch, istat, BasNam)
c         if (istat .ne. 0) then
c            print *,'@JODA_MAIN: MKVMOL returned istat = ',istat
c            call abort_job()
c         endif

c      geomopt=(iflags2(h_IFLAGS2_geom_opt).ne.0)
c      num_grad=(iflags2(h_IFLAGS2_grad_calc).eq.1)
c      fd = (iflags(h_IFLAGS_vib).eq.3.or.(geomopt.and.num_grad))

      call aces2_joda_main()

c   o finalize the job archive subsystem
c      call aces_ja_fin

      return
      end

