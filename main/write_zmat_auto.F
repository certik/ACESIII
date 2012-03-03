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
      subroutine write_zmat_auto()
c-------------------------------------------------------------------------
c   Uses information from the ZMAT file and ACES II flags to generate
c   automatic ACES III parameters which are written to the file ZMAT.AUTO.
c
c   If a *SIP section already exists in the ZMAT file, no automatic
c   parameters are generated.
c-------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'

#ifdef ACES64
      integer*8 ishell, ierr
      common /flags/ iflags
      integer*8 iflags(100)
      integer*8 iflags2(500)
#else
      integer ishell, ierr
      common /flags/ iflags
      integer iflags(100)
      integer iflags2(500)
#endif

      integer calc, vib, geom_opt, ref, dropmo, excite, instab, props
      integer grad_calc
      integer nerr
      integer ios
      integer str_trimlen
      integer i, j, n, nw, ns

      integer nauto_params
      parameter (nauto_params = 4)

      character*80 line, sip_line
      character*80 keyword
      character*80 auto_keyword(nauto_params)
      logical set_keyword(nauto_params)
      logical sip_written

      auto_keyword(1) = 'MAXMEM'
      auto_keyword(2) = 'COMPANY'
      auto_keyword(3) = 'IOCOMPANY'
      auto_keyword(4) = 'SIAL_PROGRAM'

      nerr = 0

c---------------------------------------------------------------------------
c   Get IFLAGS and IFLAGS2 off of JOBARC.
c---------------------------------------------------------------------------

      call igetrec(1,'JOBARC','IFLAGS',100,iflags)
      call igetrec(1,'JOBARC','IFLAGS2',500,iflags2)

c---------------------------------------------------------------------------
c   Determine type of calculation.
c---------------------------------------------------------------------------

      calc     = iflags(2)
      if (calc .ne. 0 .and. calc .ne. 1 .and. calc .ne. 10 .and.
     *    calc .ne. 22 .and. calc .ne. 6) then
         nerr = nerr + 1
         print *,'ACES III only supports CALC values of SCF, MBPT(2),',
     *           ' LCCSD, CCSD, and CCSD(T).'
      endif

      vib      = iflags(54)
      geom_opt = iflags2(5)
      props    = iflags(18)   ! Watson: get props.
      grad_calc = iflags2(138)

      if ((geom_opt .gt. 0) .and. (calc .eq. 22)) then
        print *,'ACES III does not support CALC=CCSD(T) with geometry ',
     *          'optimization'
        nerr = nerr + 1
      endif

c---------------------------------------------------------------------------
c   Get the REF value.
c---------------------------------------------------------------------------

      ref = iflags(11)
      if (ref .ne. 0 .and. ref .ne. 1) then
         nerr = nerr + 1
         print *,'ACES III only supports REF values of UHF, ',
     *           'RHF for auto parameterization.'
      endif

c---------------------------------------------------------------------------
c   Get the DROPMO flag.
c---------------------------------------------------------------------------

      dropmo = iflags(27)

      if (dropmo .gt. 0) then
         if (geom_opt .gt. 0) then
            if (calc .eq. 1) then
               print *,'ACES III does not support CALC=MBPT(2) with ',
     *                 'DROPMO and geometry optimization'
               nerr = nerr + 1
            endif
         endif
      endif

c---------------------------------------------------------------------------
c   Get the EXCITE flag.
c---------------------------------------------------------------------------

      excite = iflags(87)

      if (excite .gt. 0) then
         if (calc .ne. 10) then
            print *,'ACES III using EXCITE=EOMEE must use ',
     *                 'CALC=CCSD.'
            nerr = nerr + 1
         endif
      endif

c---------------------------------------------------------------------------
c   Get the HFSTABILITY flag.
c---------------------------------------------------------------------------

      instab = iflags(74)

c--------------------------------------------------------------------------
c   Abort job if it did not pass the error checking.
c--------------------------------------------------------------------------

      if (nerr .gt. 0) then
         print *,'Aborting job in write_zmat_auto due to ',nerr,
     *           ' errors'
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Open the files.
c---------------------------------------------------------------------------

      open (unit=22, file='ZMAT.AUTO', ACTION='WRITE',
     *      err = 100, iostat = ios)

  100 continue
      if (ios .ne. 0) then
         print *,'Error: Cannot open ZMAT.AUTO file.'
         call abort_job()
      endif

      open (unit=23, file='ZMAT', status='OLD',
     *      err = 200, iostat = ios)

  200 continue
      if (ios .ne. 0) then
         print *,'Error: Cannot open ZMAT file.'
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Copy everything before a *SIP section to the ZMAT file.
c---------------------------------------------------------------------------

      sip_written = .false.
  300 continue
      line = ' '
      read (23,9000, end = 400) line
      if (line(1:4) .eq. '*SIP') then
         write (22,9000) line   ! copy the '*SIP'
         sip_written = .true.
         go to 400
      endif
      write (22,9000) line   ! copy to ZMAT.AUTO
      go to 300

  400 continue

c--------------------------------------------------------------------------
c   Process the '*SIP' section?
c--------------------------------------------------------------------------

      do i = 1, nauto_params
         set_keyword(i) = .false.
      enddo

      if (sip_written) then
  500    continue
         line = ' '
         read (23,9000,end=600) line
         if (line(1:1) .eq. '*') go to 600   ! new namelist.
         if (line .eq. ' ') go to 600        ! skip blank line

c---------------------------------------------------------------------------
c   Scan the line for '='.
c---------------------------------------------------------------------------

         do i = 1, len(line)
            if (line(i:i) .eq. '=') then
               n = str_trimlen(line(1:i-1))
               keyword = line(1:n)

c----------------------------------------------------------------------------
c   See if 'keyword' is in the list of automatic params.
c----------------------------------------------------------------------------

               do j = 1,nauto_params
                  if (keyword .eq. auto_keyword(j))
     *                set_keyword(j) = .true.
               enddo
               exit
            endif
         enddo

c----------------------------------------------------------------------------
c   Now copy the line to ZMAT.AUTO.
c----------------------------------------------------------------------------

         write(22,9000) line
         go to 500
      endif

  600 continue

c--------------------------------------------------------------------------
c   If the '*SIP' record has not been written, do it now.
c--------------------------------------------------------------------------

      if (.not. sip_written) then
         sip_line = '*SIP'
         write(22,9000) sip_line
      endif

c--------------------------------------------------------------------------
c   Write out lines for any automatic parameters that have not been set in
c   the '*SIP' section.
c-------------------------------------------------------------------------

      if (.not. set_keyword(1)) then
         sip_line = 'MAXMEM=900'
         write (22,9000) sip_line
      endif

      if (.not. set_keyword(2)) then
         nw = (nprocs * 3) / 4
         nw = max(nw, 1)      ! must have at least 1 worker
         sip_line = ' '
         write (sip_line, 9100) nw   ! COMPANY parameter
         write (22,9000) sip_line
      endif

      if (.not. set_keyword(3)) then
         nw = (nprocs * 3) / 4
         nw = max(nw, 1)
         ns = nprocs - nw
         sip_line = ' '
         write (sip_line, 9200) ns   ! IOCOMPANY parameter
         write (22,9000) sip_line
      endif

      if (.not. set_keyword(4)) then
         call write_default_sial_programs(calc,dropmo,ref,geom_opt,
     *                                    vib, excite, instab, props,
     *                                    grad_calc)
      endif

c---------------------------------------------------------------------------
c   Copy all additional lines to ZMAT.AUTO.
c---------------------------------------------------------------------------

  700 continue
      line = ' '
      read (23,9000,end=9900) line
      write (22,9000) line

 9900 continue
      close(22)
      close(23)

 9000 format(a80)
 9100 format('COMPANY   = 1 1 ',i5,' 0')
 9200 format('IOCOMPANY = 2 1 ',i5,' 0')
      return
      end

