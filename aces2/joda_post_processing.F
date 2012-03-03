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
      subroutine joda_post_processing(jodadone)
c----------------------------------------------------------------------------
c   Performs post-gradient/hessian processing which was done in joda in the 
c   previous serial ACES II implementation.  Returns jodadone = 0 if 
c   convergence is not achieved, 1 if the process has converged (or if we
c   have a single-point energy calculation).
c----------------------------------------------------------------------------

      implicit none
      integer jodadone

      integer dummy, nforcecon

c      call aces_ja_init()
c      call iputrec(1,'JOBARC','FIRSTRUN',1,0)
c      call iputrec(1,'JOBARC','DIRTYFLG',1,1)
c      call geopt(jodadone)
c      call aces_ja_fin()

      call aces2_joda_main()
      call aces_ja_init()
      call igetrec(1, 'JOBARC', 'JODADONE',  1, jodadone)
c      if (jodadone .eq. 0) then

c--------------------------------------------------------------------------
c   For VIB_FINDIF=EXACT calculations, joda doesn't set the jodadone flag.
c   Instead, we check to see if the FORCECON record is on JOBARC.
c---------------------------------------------------------------------------

c         call igetrec(0, 'JOBARC', 'FORCECON', nforcecon, dummy)
c         if (nforcecon .gt. 0) jodadone = 1
c         print *,'JODA_POST_PROCESSING: nforcecon, jodadone ',
c     *       nforcecon, jodadone
c      endif

      call aces_ja_fin()
      return
      end
