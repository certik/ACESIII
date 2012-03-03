
c This routine reads the command-line arguments and initializes the
c a2proc environment (and any other dependents).

      subroutine a2proc_init(module)
      implicit none

c ARGUMENT LIST
      character*80 module

c INTERNAL VARIABLES
      logical have_command
      character arg*80
      integer*4 i, iargc, num_args

c ----------------------------------------------------------------------

#ifdef _TRACE_A2PROC
      write(*,*) '@A2PROC: Entered file ', __FILE__
#endif

c ----------------------------------------------------------------------

c   o initialize module*80 and have_command
      module = ''
      have_command = .false.

c Make sure the integers are all of the same size. The Sun Pro 6 compilers,
c for example, require the index argument of getarg to be 4 bytes. This may
c be different with other compilers.
      num_args = iargc()
      if (1.gt.num_args) then
         write(*,*) '@INIT_A2PROC: This program requires an argument.\n'
         call exit(1)
      end if

c   o loop over arguments
      do i = 1, num_args
         if (.not.have_command) then
            call getarg(i,arg)

c         o help
            if ((arg(1:4).eq.'help').or.(arg(1:2).eq.'-h')) then
               module(1:4) = 'help'
               have_command = .true.
            end if

c         o MOLDEN
            if (arg(1:6).eq.'molden') then
               module(1:6) = 'molden'
               have_command = .true.
            end if

c         o GENNBO
            if (arg(1:6).eq.'gennbo') then
               module(1:6) = 'gennbo'
               have_command = .true.
            end if

c         o JOB ARCHIVE CONTROLLER
            if (arg(1:7).eq.'rmfiles') then
               module(1:7) = 'rmfiles'
               have_command = .true.
            end if

c         o JOB ARCHIVE CONTROLLER
            if (arg(1:5).eq.'jasum') then
               module(1:5) = 'jasum'
               have_command = .true.
            end if

c         o I/O CONTROLLER
            if (arg(1:5).eq.'iosum') then
               module(1:5) = 'iosum'
               have_command = .true.
            end if

c         o 2-particle tau intermediate (T2+T1T1)
            if (arg(1:4).eq.'tau2') then
               module(1:4) = 'tau2'
               have_command = .true.
            end if

c   o end loop over arguments
c        end if (.not.have_command)
         end if
c     end do i = 1, num_args
      end do

c   o stop if there is no command
      if (.not.have_command) then
         write(*,*) '@INIT_A2PROC: The supplied arguments are not ',
     &              'recognized.\n'
         call exit(1)
      end if

c ----------------------------------------------------------------------

c VERIFY CONSISTENCY

c for most modules, JOBARC must exist
      inquire(file='JOBARC',exist=have_command)
      if (.not.have_command) then
         write(*,*) '@INIT_A2PROC: There is no JOBARC file, which ',
     &              'probably means'
         write(*,*) '              ACES2 has not been run.\n'
         call exit(1)
      end if

#ifdef _DEBUG
#endif

c ----------------------------------------------------------------------

#ifdef _TRACE_A2PROC
      write(*,*) '@A2PROC: Leaving file ', __FILE__
#endif

      return
      end

