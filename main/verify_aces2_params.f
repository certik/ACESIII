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
      subroutine verify_aces2_params()
c----------------------------------------------------------------------------
c   Checks for valid combinations of parameters in the *ACES2 section.
c   If an invalid combination (i. e. not supported by ACESIII) is found,
c   an error message is printed and the job is aborted.
c----------------------------------------------------------------------------

      implicit none

      common /flags/ iflags
      common /flags2/ iflags2
      integer iflags(100)
      integer iflags2(500)

      integer vib, geom_opt, calc, opt_method, coords, init_hessian,
     *        hess_update, dropmo, ecp
      integer grad_calc

      call igetrec(1,'JOBARC','IFLAGS',100,iflags)
      call igetrec(1,'JOBARC','IFLAGS2',500,iflags2)

      vib      = iflags(54)
      geom_opt = iflags2(5)
      calc     = iflags(2)
      opt_method   = iflags(47)
      coords       = iflags(68)
      init_hessian = iflags2(8)
      hess_update  = iflags2(7)
      dropmo       = iflags(27)
      ecp          = iflags(71)
      grad_calc    = iflags2(138)

      if (vib .eq. 1 .and. geom_opt .gt. 0) then
         print *,'ERROR: ACESIII does not support VIB=FINDIF ',
     *           'with geometry optimization.'
         print *,'Either remove the geometry optimization ',
     *           'parameters or remove VIB_FINDIF'
         call abort_job()   
      endif

      if (vib .eq. 1 .and. dropmo .gt. 0) then
         print *,'ERROR: ACESIII does not support Hessian ',
     *           'calculations using DROPMO.'
         call abort_job()
      endif

      if (ecp .gt. 0) then
         print *,'ERROR: ACESIII does not support ECP calculations.'
         call abort_job()
      endif

      if (grad_calc .eq. 2 .and. vib .eq. 1) then
         print *,'Error: GRAD_CALC = NUMERICAL and VIB = FINDIF'
         print *,'       Parameter inconsistency'
         call abort_job()
      endif
      return
      end    
