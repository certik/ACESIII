
c NOTE: This routine was blatantly ripped from netlib/TOMS.
c       I, Anthony Yau, removed all those GOTO statements and
c       replaced them with DO WHILE loops and a minimum of
c       code replication.
c       (cf. http://www.netlib.org/toms/513)

c ---------------------------------------------
c --- TRANSPOSE A RECTANGULAR ARRAY IN SITU ---
c ---------------------------------------------

c CALL DXPOSE(A,M,N,ADIM,WORK,WORKDIM,IERR)
c
c A    = { double precision array of single dimension ADIM }
c        the double precision array to be transposed
c
c M, N = { integer }
c        the dimensions of A(M,N) before transposition
c        stored columnwise
c
c ADIM = { integer }
c        M*N
c
c WORK = { integer array }
c        a work array
c
c WORKDIM = { integer }
c           the dimension of WORK(), recommended to be (M+N)/2
c
c IERR = { integer }
c        a return code
c      = 0  := normal
c      = -1 := ADIM.NE.(M*N)
c      = -2 := WORKDIM.LE.0
c      > 0  := (rare) the final value of i when the search is
c              complete but some loops have not been moved
c              NOTE: work(i) will stay zero for fixed points

      subroutine dxpose(a,m,n,adim,work,workdim,ierr)
      implicit none

c ARGUMENT LIST
      integer m, n, adim
      double precision a(adim)
      integer workdim, work(workdim), ierr

c INTERNAL VARIABLES
      double precision b, c, d
      integer n1, ncount, max, im
      integer j, j1
      integer i, i1, i2, i1c, i2c
      integer k, kmi
      integer itmp, ir1, ir2
      logical still_going

c ----------------------------------------------------------------------

      ierr = 0
      if (adim.ne.m*n)  ierr = -1
      if (workdim.lt.1) ierr = -2
      if ((ierr.ne.0).or.(m.lt.2).or.(n.lt.2)) return

      if (m.eq.n) then
         n1 = n-1
         do i = 1, n1
            j1 = i+1
            do j = j1, n
               i1 = i + (j-1)*n
               i2 = j + (i-1)*m
               b     = a(i1)
               a(i1) = a(i2)
               a(i2) = b
            end do
         end do
         return
      end if

      ncount = 2
      k = adim-1
      do itmp = 1, workdim
         work(itmp) = 0
      end do

      if ((m.gt.2).and.(n.gt.2)) then
c      o calculate the number of fixed points using Euclid's algorithm
c        for GCD(m-1,n-1)
         ir2  = m-1
         ir1  = n-1
         itmp = 1
         do while (itmp.ne.0)
            itmp = mod(ir2,ir1)
            ir2  = ir1
            ir1  = itmp
         end do
         ncount = ir2+1
      end if

c     set initial values for search
      i  = 1
      im = m

c     we return only when ncount.ge.adim (or we error)
      do while (ncount.lt.adim)
         i1  = i
         kmi = k-i
         b   = a(i1+1)
         i1c = kmi
         c   = a(i1c+1)
         still_going = .true.
         do while (still_going)
c         o leave this line as is since k*(i1/n) recasts the product
            i2  = m*i1 - k*(i1/n)
            i2c = k-i2
            if (i1 .le.workdim) work(i1 ) = 1
            if (i1c.le.workdim) work(i1c) = 1
            ncount = ncount + 2
            if (i2.eq.i) then
               a(i1+1)  = b
               a(i1c+1) = c
               still_going = .false.
            else
               if (i2.eq.kmi) then
                  d = b
                  b = c
                  c = d
                  a(i1+1)  = b
                  a(i1c+1) = c
                  still_going = .false.
               else
                  a(i1+1)  = a(i2+1)
                  a(i1c+1) = a(i2c+1)
                  i1  = i2
                  i1c = i2c
c              end if (i2.eq.kmi)
               end if
c           end if (i2.eq.i)
            end if
c        end do while (still_going)
         end do
         if (ncount.lt.adim) then
            still_going = .true.
            do while (still_going)
               max = k-i
               i   = i+1
               if (i.gt.max) then
                  ierr = i
                  return
               end if
               if (im.gt.k) then
                  im = im-k
               else
                  im = im+m
               end if
               i2 = im
               if (i.ne.i2) then
                  if (i.gt.workdim) then
                     if ((i.lt.i2).and.(i2.lt.max)) then
                        do while ((i.lt.i2).and.(i2.lt.max))
                           i1 = i2
                           i2 = m*i1 - k*(i1/n)
                        end do
                     end if
c                     if (i2.eq.i) still_going = .false.
                     still_going = (i2.ne.i)
                  else
c                     if (work(i).eq.0) still_going = .false.
                     still_going = (work(i).ne.0)
                  end if
               end if
c           end do while (still_going)
            end do
c        end if (ncount.lt.adim)
         end if
c     end do while (ncount.lt.adim)
      end do

      return
      end

