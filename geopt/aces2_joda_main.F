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

#include "flags.h"

      subroutine aces2_joda_main()
      implicit none

#include "fnamelen.par" 
      character*(fnamelen)fname
      integer pass1, i, i_havegeom, ipost_vib, ilength
      logical fd, geomopt, num_grad, vib_frqs, we_havegeom,
     &        optarc_presnt

      integer iMemMin, iMemInc, Popt_num_fd
      parameter (iMemMin=2097152)
      parameter (iMemInc=1048576)

      integer*8 malloced_previously
      integer*8 ixz
      integer*8 get_max_heap_usage
      integer*8 symcor_mem_bytes, maxheap_bytes
      integer*8 gentranmat_mem_bytes
      integer symcor_ints
      integer gentranmat_ints
 
      integer ishell
      external ishell

#include "icore.com"
#include "istart.com"
      logical ignore
      common /restart_com/ ignore
#include "jodaflags.com"
#include "mxatms.par"
      include 'machine_types.h'

      integer memreq, m, ierr
      double precision z(1)
      integer iz(1)
      equivalence (z, iz)
#ifdef ALTIX
      include 'sheap.h'
      pointer (dptr, z)
#else
#ifndef BLUEGENE
      common z
#endif
#endif

#ifdef _TRAPFPE
c   o since JODA does not call crapsi, have it turn on FPEs
      call trapfpe
#endif /* _TRAPFPE */

#ifdef ALTIX
      dptr = dshptr
#endif

       m      = max(mxatms, maxredunco )
       memreq = 15*m*m + 34*m

c   o gather parallel statistics (needed by gfname in ja_init)
c      call aces_com_parallel_aces

c   o parse command line (overwrites rank and number of processes)
c      call parse_cli

c   o evaluate and repair the health of the current file set
      call dfiles(ignore)

c   o initialize the job archive subsystem
      call aces_ja_init

      if (.not.ignore) then
c      o this is the first joda run
         call iputrec(1,'JOBARC','FIRSTRUN',1,1)
         call iputrec(1,'JOBARC','DIRTYFLG',1,0)
         call iputrec(1,'JOBARC','JODADONE',1,0)
         call iputrec(1,'JOBARC','FNDFDONE',1,1)

c------------------------------------------------------------------------------
c   Allocate heap memory for GEOPT.
c------------------------------------------------------------------------------

         call mem_alloc(z, memreq, bytes_per_double, ixz, .true., ierr)
         if (ierr .ne. 0) then
            print *,'ACES2_JODA_MAIN: failure allocating Z array'
            print *,'memreq = ',memreq
            call abort_job()
         endif

         call geopt(z(ixz), memreq)
         call mem_alloc_reset()         ! release the memory.
      else
c      o load flags
         call igetrec(1,'JOBARC','IFLAGS', 100,iflags)
         call igetrec(1,'JOBARC','IFLAGS2',500,iflags2)
      end if
c
      call igetrec(-1,'JOBARC', 'HAVEGEOM',1, i_havegeom)
      call igetrec(-1, 'JOBARC', 'POSTOPFD', 1, popt_num_fd)
c
      call gfname('OPTARC', fname, ilength)
      inquire(file=fname(1:ilength), exist=optarc_presnt)
c
      If (i_havegeom .eq. 1) we_havegeom = .true.
      geomopt=(iflags2(h_IFLAGS2_geom_opt).ne.0) .and.
     &         .not. we_havegeom
      num_grad=(iflags2(h_IFLAGS2_grad_calc).eq.2)
      vib_frqs=(iflags(h_IFLAGS_vib).eq.3)
c
      fd = ((iflags(h_IFLAGS_vib).eq.3 .and. we_havegeom)
     &     .or. num_grad .or. (geomopt.and.num_grad))
c
#ifdef _DEBUG_LVL0
      Print*, "Finite diffs; ignore and geomopt vars:",
     &         ignore, geomopt
      Print*, "The vib calc. related varrs:"
      Print*, "iflags(h_IFLAGS_vib):", iflags(h_IFLAGS_vib)
      Print*, "The finite diffs:", fd
      Print*, "iflags2(h_IFLAGS2_geom_opt):",
     &         iflags2(h_IFLAGS2_geom_opt)
      Print*, "first run of popt num. frq:", popt_num_fd
 
#endif
c
      if (fd) then
c         print *,'Finite differences are not currently supported.'
c         call errex()

c      o allocate memory for symcor
c         icrsiz = iflags(h_IFLAGS_mem)
          
c         do while ((iCore(1).eq.0).and.(icrsiz.gt.iMemMin))
c            call aces_malloc(icrsiz,iCore,i0, symcor_heap)
c            if (symcor_heap.eq.0) icrsiz = icrsiz - iMemInc
c         end do
c         if (symcor_heap.eq.0) then
c            print *, '@JODA: unable to allocate at least ',
c     &               iMemMin,' integers of memory'
c            call aces_exit(1)
c            call errex() 
c         end if
c
          if (ignore) then 
c         o mid-stream -> keep going
            call igetrec(1,'JOBARC','FNDFDONE',1,i)
            if (popt_num_fd .eq. 1) call post_opt_update
            if (i.eq.0) then

c------------------------------------------------------------------------
c   Set up heap for symcor.
c------------------------------------------------------------------------

               maxheap_bytes = get_max_heap_usage()
               call mem_alloc_query(malloced_previously)
               symcor_mem_bytes = maxheap_bytes - malloced_previously
               symcor_ints = symcor_mem_bytes / intsize 
               call mem_alloc(iz, symcor_ints, intsize, ixz, .true., 
     *                     ierr)
               if (ierr .ne. 0) then
                  print *,'SYMCOR MEM ALLOC ERROR'
                  print *,'maxheap_bytes, malloced_previously ', 
     *                     maxheap_bytes, malloced_previously
                  print *,'symcor_ints ',symcor_ints
                  call abort_job() 
               endif 
               call symcor(iz(ixz), symcor_ints)

c-------------------------------------------------------------------------
c   Undo the symcor heap.
c-------------------------------------------------------------------------

               call mem_alloc_reset()
            endif 
         else
c         o first run -> reset finite difference series
            if (geomopt) i=ishell('cp OPTARC OPTARCBK')
            call iputrec(1,'JOBARC','FNDFDONE',1,0)

c------------------------------------------------------------------------
c   Set up heap for symcor.
c------------------------------------------------------------------------

            maxheap_bytes = get_max_heap_usage()
            call mem_alloc_query(malloced_previously)
            symcor_mem_bytes = maxheap_bytes - malloced_previously
            symcor_ints      = symcor_mem_bytes / intsize 
            call mem_alloc(iz, symcor_ints, intsize, ixz, .true., 
     *                     ierr)
            if (ierr .ne. 0) then
               print *,'SYMCOR MEM ALLOC ERROR'
               print *,'maxheap_bytes, malloced_previously ',
     *              maxheap_bytes, malloced_previously
               print *,'symcor_ints ',symcor_ints
               call abort_job() 
            endif 
            call symcor(iz(ixz), symcor_ints)

c-------------------------------------------------------------------------
c   Undo the symcor heap.
c-------------------------------------------------------------------------

            call mem_alloc_reset()
            ignore=.true.
         end if

         call igetrec(1,'JOBARC','PASS1',1,pass1)
         if (pass1.ne.-1) then
c         o vib freqs w/ an grads -OR- geom opts w/ num grads
            call iputrec(1,'JOBARC','FIRSTRUN',1,0)
            call iputrec(1,'JOBARC','DIRTYFLG',1,1)
            if (geomopt.and.pass1.eq.0) then
               i=ishell('cp OPTARCBK OPTARC')

c------------------------------------------------------------------------------
c   Allocate heap memory for GEOPT.
c------------------------------------------------------------------------------

               call mem_alloc(z, memreq, bytes_per_double, ixz, 
     *                        .true., ierr)
               if (ierr .ne. 0) then
                  print *,'ACES2_JODA_MAIN: failure allocating Z array'
                  print *,'memreq = ',memreq
                  call abort_job()
               endif

               call geopt(z(ixz), memreq)
               call mem_alloc_reset()         ! release the memory.

               call igetrec(1,'JOBARC','JODADONE',1,i)
               if (i.ne.1) then
c               o new geom -> reset finite difference series
                  i=ishell('cp OPTARC OPTARCBK')
                  call iputrec(1,'JOBARC','FNDFDONE',1,0)
c------------------------------------------------------------------------
c   Set up heap for symcor.
c------------------------------------------------------------------------

                  maxheap_bytes = get_max_heap_usage()
                  call mem_alloc_query(malloced_previously)
                  symcor_mem_bytes = maxheap_bytes - malloced_previously
                  symcor_ints      = symcor_mem_bytes / intsize 
                  call mem_alloc(iz, symcor_ints, intsize, ixz, .true., 
     *                     ierr)
                  if (ierr .ne. 0) then
                     print *,'SYMCOR MEM ALLOC ERROR'
                     print *,'maxheap_bytes, malloced_previously ',
     *                    maxheap_bytes, malloced_previously
                     print *,'symcor_ints ',symcor_ints
                     call abort_job() 
                  endif 
                  call symcor(iz(ixz), symcor_ints)

c-------------------------------------------------------------------------
c   Undo the symcor heap.
c-------------------------------------------------------------------------

                  call mem_alloc_reset()

c------------------------------------------------------------------------------
c   Allocate heap memory for GEOPT.
c------------------------------------------------------------------------------

                  call mem_alloc(z, memreq, bytes_per_double, ixz,
     *                        .true., ierr)
                  if (ierr .ne. 0) then
                     print *,'ACES2_JODA_MAIN: failure allocating ',
     *                       'Z array'
                     print *,'memreq = ',memreq
                     call abort_job()
                  endif

                  call geopt(z(ixz), memreq)
                  call mem_alloc_reset()         ! release the memory.
               end if
            else if (vib_frqs .or. num_grad) then 
               if (we_havegeom .and. optarc_presnt) 
     &             i=ishell('rm OPTARC')
               call igetrec(1,'JOBARC','FNDFDONE',1,i)
c------------------------------------------------------------------------------
c   Allocate heap memory for GEOPT.
c------------------------------------------------------------------------------

               call mem_alloc(z, memreq, bytes_per_double, ixz,
     *                        .true., ierr)
               if (ierr .ne. 0) then
                  print *,'ACES2_JODA_MAIN: failure allocating Z array'
                  print *,'memreq = ',memreq
                  call abort_job()
               endif

               call geopt(z(ixz), memreq)
               call mem_alloc_reset()         ! release the memory.

               if (pass1.eq.0) call iputrec(1,'JOBARC','JODADONE',1,1)
               call igetrec(1,'JOBARC','FNDFDONE',1,i) 
            end if
         end if

      else

         if (ignore) then
c         o joda has run before
            call iputrec(1,'JOBARC','FIRSTRUN',1,0)
            call iputrec(1,'JOBARC','DIRTYFLG',1,1)
c------------------------------------------------------------------------------
c   Allocate heap memory for GEOPT.
c------------------------------------------------------------------------------

            call mem_alloc(z, memreq, bytes_per_double, ixz,
     *                        .true., ierr)
            if (ierr .ne. 0) then
               print *,'ACES2_JODA_MAIN: failure allocating Z array'
               print *,'memreq = ',memreq
               call abort_job()
            endif

            call geopt(z(ixz), memreq)
            call mem_alloc_reset()         ! release the memory.
         end if

      end if


      maxheap_bytes = get_max_heap_usage()
      call mem_alloc_query(malloced_previously)
      gentranmat_mem_bytes = maxheap_bytes - malloced_previously
      gentranmat_ints = gentranmat_mem_bytes / intsize
      call mem_alloc(iz, gentranmat_ints, intsize, ixz, .true.,
     *               ierr)
      if (ierr .ne. 0) then
         print *,'GENTRANMAT MEM ALLOC ERROR'
         print *,'maxheap_bytes, malloced_previously ',
     *            maxheap_bytes, malloced_previously
         print *,'gentranmat_ints ',gentranmat_ints
         call abort_job()
      endif

      call gentranmat(iz(ixz), gentranmat_ints)

c---------------------------------------------------------------------------
c   Free any heap memory that might have been used.  It must be made
c   available for the ACES III code's use.
c---------------------------------------------------------------------------

      call mem_alloc_reset()

c----------------------------------------------------------------------------
c   For vib_findif=exact, set jodadone flag to 1.
c----------------------------------------------------------------------------

      if (iflags(h_IFLAGS_vib) .eq. 1) then
         call iputrec(1,'JOBARC', 'JODADONE', 1, 1)
      endif

c   o finalize the job archive subsystem
      call aces_ja_fin

      end

