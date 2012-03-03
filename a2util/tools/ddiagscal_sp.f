
c ---------------------------------------------------------------------
c --- SCALE THE (OFF) DIAGONAL ELEMENTS OF A PACKED SYMMETRIC ARRAY ---
c ---------------------------------------------------------------------

c CALL DDIAGSCAL_SP(UPLO,DORO,M,ALPHA,DX)
c
c UPLO  = { "U", "u", "L", "l" }
c         the condition of the packed array
c
c DORO  = { "D", "d", "O", "o" }
c         controls whether the diagonal or off-diagonal elements are
c         to be scaled
c
c M     = { integer }
c         the order of the symmetric array
c
c ALPHA = { double precision }
c         the scaling factor
c
c DX    = { double precision array with dimension M*(M+1)/2 }
c         the packed array to be scaled

      subroutine ddiagscal_sp(uplo,doro,m,alpha,dx)
      implicit none

c ARGUMENT LIST
      character*1 uplo, doro
      integer m
      double precision alpha, dx(m*(m+1)/2)

c INTERNAL VARIABLES
      logical diag
      integer offset, i, from, to, step
      double precision factor

c DISCARD THE CHAFF AND THE SIMPLE CASES
      if (m.lt.0) then
         print *, "@DDIAGSCAL_SP: The array to be scaled cannot",
     &            " have a dimension"
         print *, "               less than zero."
         call c_exit(1)
      end if
      if (m.eq.0) return
      if (alpha.eq.(1.0D0)) return
      if ((doro.eq."D").or.(doro.eq."d")) then
         diag=.TRUE.
      else if ((doro.eq."O").or.(doro.eq."o")) then
         diag=.FALSE.
      else
         print *, "@DDIAGSCAL_SP: An invalid character has been",
     &            " received for DorO."
         print *, "               The string ",doro," should be",
     &            " D, d, O, or o."
         call c_exit(1)
      end if
      if (m.eq.1) then
         if (diag) dx(1)=dx(1)*alpha
         return
      end if
      if (m.eq.2) then
         if (diag) then
            dx(1)=dx(1)*alpha
            dx(3)=dx(3)*alpha
         else
            dx(2)=dx(2)*alpha
         end if
         return
      end if

c SET INCREMENT PARAMETERS
c   o For the case where m is 3, we can hard-code the results and save
c     some time. I know it is not a lot of time, but it is my call and
c     I made it.
      if ((uplo.eq."U").or.(uplo.eq."u")) then
         if (m.eq.3) then
            if (diag) then
               dx(1)=dx(1)*alpha
               dx(3)=dx(3)*alpha
               dx(6)=dx(6)*alpha
            else
               dx(2)=dx(2)*alpha
               dx(4)=dx(4)*alpha
               dx(5)=dx(5)*alpha
            end if
            return
         end if
         from =  2
         to   =  m+1
         step =  1
      else if ((uplo.eq."L").or.(uplo.eq."l")) then
         if (m.eq.3) then
            if (diag) then
               dx(1)=dx(1)*alpha
               dx(4)=dx(4)*alpha
               dx(6)=dx(6)*alpha
            else
               dx(2)=dx(2)*alpha
               dx(3)=dx(3)*alpha
               dx(5)=dx(5)*alpha
            end if
            return
         end if
         from =  m
         to   =  1
         step = -1
      else
         print *, "@DDIAGSCAL_SP: An invalid character has been",
     &            " received for UpLo."
         print *, "               The string ",uplo," should be",
     &            " U, u, L, or l."
         call c_exit(1)
      end if

c INVERT THE PROBLEM WHEN .NOT.DIAG (WATCH OUT WHEN ALPHA=0.0D0)
c   o Note: I am still looking for a fast way to zero-out a float.
c           If these were integers, we could IAND(number,0) but I
c           cannot find any such intrinsic function for floats.
c     Also, I examined the time cost of filtering low-order matrices
c     out (instead of calling dscal), but the startup cost of even
c     entering this subroutine is much higher than anything you
c     may gain by tweaking the cut-off parameter. KEEP IT SIMPLE!
      if (.not.diag) then
c         if ((alpha.eq.(0.0D0)).or.(m.lt.10)) then
         if (alpha.eq.(0.0D0)) then
            if ((uplo.eq."L").or.(uplo.eq."l")) then
               offset = 1
               do step = m,2,-1
                  do i = 1,(step-1)
                     dx(offset+i)=dx(offset+i)*alpha
                  end do
                  offset=offset+step
               end do
            else
               offset = 1
               do step = 1,(m-1)
                  do i = 1,step
                     dx(offset+i)=dx(offset+i)*alpha
                  end do
                  offset=offset+step+1
               end do
            end if
            return
         end if
         call xscal(ishft((m*(m+1)),-1),alpha,dx,1)
         factor=1/alpha
      else
         factor=alpha
      end if

c GET TO WORK
      offset=1
      do i = from,to,step
         dx(offset)=dx(offset)*factor
         offset=offset+i
      end do

      return
      end

