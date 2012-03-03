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
      subroutine write_default_sial_programs(calc,dropmo,ref,geom_opt,
     *                                       vib, excite, instab, props,
     *                                       grad_calc)
c-------------------------------------------------------------------------
c   Writes the SIAL_PROGRAM parameters for a default set of SIAL programs
c   determined by the parameters calc, dropmo, ref, geom_opt, and vib.
c   The SIAL_PROGRAM parameters are written to ZMAT.AUTO.
c-------------------------------------------------------------------------
      implicit none
      integer calc,dropmo,ref,geom_opt, vib, excite, instab, props
      integer grad_calc 
      character*80 jobflow
      integer ierr
      integer n, str_trimlen

c---------------------------------------------------------------------------
c   Determine the jobflow required for the combination of parameters.
c---------------------------------------------------------------------------

c      print *,'ref, calc, dropmo, geom_opt, vib ',
c     *   ref, calc, dropmo, geom_opt, vib, excite
      jobflow = 'UNDEFINED'

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'SCF_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'SCF_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'SCF_RHF_GRADIENT'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'SCF_UHF_GRADIENT'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'SCF_RHF_HESSIAN'     ! analytical hessian calc
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 0 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 3) then
         if (grad_calc .eq. 2) then
            jobflow = 'SCF_RHF_ENERGY'    ! numerical gradient calc
         else
            jobflow = 'SCF_RHF_GRADIENT'    ! numerical hessian calc
         endif 
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'SCF_UHF_HESSIAN'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'SCF_UHF_ENERGY'
         else
            jobflow = 'SCF_UHF_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'MP2_RHF_HESSIAN'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'MP2_RHF_ENERGY'
         else
            jobflow = 'MP2_RHF_GRADIENT'
         endif
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 1)) then
         jobflow = 'MP2_UHF_HESSIAN'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'MP2_UHF_ENERGY'
         else
            jobflow = 'MP2_UHF_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'MP2_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .eq. 0 .and. vib .eq. 0) then
         jobflow = 'MP2_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'MP2_RHF_GRADIENT'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 1 .and.
     *    geom_opt .gt. 0) then
         jobflow = 'MP2_UHF_GRADIENT'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_RHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'LCCSD_UHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'LCCSD_RHF_ENERGY'
         else
            jobflow = 'LCCSD_RHF_GRADIENT'
         endif
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'LCCSD_UHF_ENERGY'
         else 
            jobflow = 'LCCSD_UHF_GRADIENT'
         endif  
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         jobflow = 'LCCSD_RHF_DROPMO_GRADIENT'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 6 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         jobflow = 'LCCSD_UHF_DROPMO_GRADIENT'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_DROPMO_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_RHF_ENERGY'
         else
            jobflow = 'CCSD_RHF_GRADIENT'
         endif
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_UHF_ENERGY'
         else
            jobflow = 'CCSD_UHF_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_RHF_DROPMO_ENERGY'
         else
            jobflow = 'CCSD_RHF_DROPMO_GRADIENT'
         endif 
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 10 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .gt. 0 .or. vib .eq. 3)) then
         if (grad_calc .eq. 2) then
            jobflow = 'CCSD_UHF_DROPMO_GRADIENT'
         else
            jobflow = 'CCSD_UHF_DROPMO_GRADIENT'
         endif
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 22 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_TRIPLES_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 22 .and.
     *    dropmo .eq. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_TRIPLES_ENERGY'
      endif

      if (ref .eq. 0 .and.
     *    calc .eq. 22 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_RHF_DROPMO_TRIPLES_ENERGY'
      endif

      if (ref .eq. 1 .and.
     *    calc .eq. 22 .and.
     *    dropmo .ne. 0 .and.
     *    (geom_opt .eq. 0 .and. vib .eq. 0)) then
         jobflow = 'CCSD_UHF_DROPMO_TRIPLES_ENERGY'
      endif

      n = str_trimlen(jobflow)
      if (instab .gt. 0) jobflow((n+1):(n+7)) = '_INSTAB'

      if (excite .eq. 3) then
         if (ref .eq. 0 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0) then
            jobflow = 'EOM_RHF_CCSD_ENERGY'
         endif

         if (ref .eq. 0 .and.
     *       calc .eq. 10 .and.
     *       dropmo .ne. 0) then
            jobflow = 'EOM_RHF_CCSD_DROPMO_ENERGY'
         endif

         if (ref .eq. 1 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0) then
            jobflow = 'EOM_UHF_CCSD_ENERGY'
         endif

         if (ref .eq. 1 .and.
     *       calc .eq. 10 .and.
     *       dropmo .ne. 0) then
            jobflow = 'EOM_UHF_CCSD_DROPMO_ENERGY'
         endif
C
C   Watson
C
         if (ref .eq. 0 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0  .and.
     *       props  .eq. 1)   then
            jobflow = 'EOM_RHF_CCSD_DENS_ENERGY'
         endif

         if (ref .eq. 1 .and.
     *       calc .eq. 10 .and.
     *       dropmo .eq. 0  .and.
     *       props   .eq. 1)  then
            jobflow = 'EOM_UHF_CCSD_DENS_ENERGY'
         endif
C
C   Watson
C
      endif

      n = str_trimlen(jobflow)
      print *,'Using Default jobflow = ',jobflow(1:n)
      call find_default_jobflow(jobflow(1:n), ierr)
      if (ierr .eq. 0) then
         call write_default_jobflow(22, jobflow(1:n))
      else
         print *,'Error: Cannot find jobflow ',jobflow(1:n),
     *           ' in default_jobflows file.'
         call abort_job()
      endif

      return
      end

      subroutine find_default_jobflow(jobflow, ierr)
c---------------------------------------------------------------------------
c   Position default_jobfiles file to desired jobflow.  ierr > 0 indicates
c   the jobflow was not found.
c---------------------------------------------------------------------------
      implicit none
      character*(*) jobflow
      integer ierr

      character*256 fn, envvar
      character*80 line

      integer iunit, ios, lenenv, n, str_trimlen

c--------------------------------------------------------------------------
c   Check for existence of $ACES_EXE_PATH/sial_config
c--------------------------------------------------------------------------

      iunit = 12

      call c_getenv('ACES_EXE_PATH'//char(0), envvar, lenenv, ierr)
      if (ierr .eq. 0) then
         fn = envvar(1:lenenv) // '/sio/default_jobflows'
         open (unit=iunit, file = fn, status = 'OLD',
     *      err=600, iostat = ios)
      endif

  100 continue

c---------------------------------------------------------------------------
c   Read each line, checking for JOBFLOW statements.
c---------------------------------------------------------------------------

      read (iunit, '(A80)', end = 200) line
      if (line(1:7) .eq. 'JOBFLOW') then
         n = str_trimlen(line)
         if (line(9:n) .eq. jobflow) then

c---------------------------------------------------------------------------
c   The jobflow name matches, we have found it.
c---------------------------------------------------------------------------

            ierr = 0
            return
         endif
      endif

      go to 100

  200 continue
      ierr = 1
      return

  600 continue
      ierr = 2
      print *,'Error: Could not open file: ',fn
      return
      end

      subroutine write_default_jobflow(iunit, jobflow)
c---------------------------------------------------------------------------
c   Reads data until an "ENDJOBFLOW" is reached.  Writes the filenames to
c   the ZMAT.AUTO file.
c---------------------------------------------------------------------------

      implicit none
      integer iunit, n, str_trimlen
      character*(*) jobflow
      character*80 line , pline

  100 continue
      read (12, '(A80)', end=200) line
      n = str_trimlen(line)
      if (line(1:10) .eq. 'ENDJOBFLOW') then
         close(12)
         return
      endif

      pline = 'SIAL_PROGRAM=' // line(1:n)
      write(iunit, 9000) pline
      go to 100

  200 continue
      print *,'Error: Missing ENDJOBFLOW statement for jobflow ',
     *   jobflow
      call abort_job()

 9000 format(a80)
      return
      end
