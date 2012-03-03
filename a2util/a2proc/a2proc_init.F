
c This routine reads the command-line arguments and initializes the
c a2proc environment (and any other dependents).

      subroutine a2proc_init(module,args,dimargs)
      implicit none

c ARGUMENT LIST
      integer dimargs
      character*80 module
      character*80 args(*)

c INTERNAL VARIABLES
      logical bExist
      integer i, f_iargc, num_args
      integer iJOBARC
      character*80 szJOBARC

c ----------------------------------------------------------------------

#ifdef _TRACE_A2PROC
      write(*,*) '@A2PROC: Entered file ', __FILE__
#endif

c ----------------------------------------------------------------------

      num_args = f_iargc()
      if (1.gt.num_args) then
         write(*,*) '@INIT_A2PROC: This program requires an argument.\n'
         call aces_exit(1)
      end if

c   o get the module name
      module = ' '
      call f_getarg(1,module)

c   o get additional arguments
      do i = 1, min(dimargs,num_args-1)
         call f_getarg(1+i,args(i))
      end do
      if (num_args-1.gt.dimargs) then
         print *, '@A2PROC_INIT: a2proc can only handle the first ',
     &            dimargs,' arguments'
         print *, '              after the module name'
      else
         dimargs = num_args - 1
      end if

c ----------------------------------------------------------------------

c VERIFY CONSISTENCY

c for most modules, JOBARC must exist
      if (module.eq.'help'.or.module.eq.'-h') return
      if (module.eq.'factor') return
      if (module.eq.'PES_scan') return
      call gfname('JOBARC',szJOBARC,iJOBARC)
      inquire(file=szJOBARC(1:iJOBARC),exist=bExist)
      if (.not.bExist) then
         write(*,*) '@INIT_A2PROC: There is no JOBARC file, which ',
     &              'probably means'
         write(*,*) '              ACES2 has not been run.\n'
         call c_exit(1)
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_A2PROC
      write(*,*) '@A2PROC: Leaving file ', __FILE__
#endif

      return
      end

