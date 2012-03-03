      double precision function dsum(n,dx,incx)
c
c     takes the sum of the values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c     modified 11/4/03, dasum->dsum conversion
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dsum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dx(i)
   10 continue
      dsum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dx(i) + dx(i + 1) + dx(i + 2)
     *  + dx(i + 3) + dx(i + 4) + dx(i + 5)
   50 continue
   60 dsum = dtemp
      return
      end
