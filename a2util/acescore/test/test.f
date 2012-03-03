
      program main
      implicit none

      integer i, j, k1, k2, icomp, ifile, irec

      integer ipack, upackr, upackf
      integer ipack2, upackr2, upackf2

      IPACK(I,J)=IOR(ISHFT(I-49,29),J)
      UPACKR(I) =IAND(I,((2**29)-1))
      UPACKF(I) =IAND(ISHFT(I,-29),7)+49

      ipack2(ifile,irec) = ior(lshift(ifile-49,29),irec)
      upackr2(icomp) = iand(icomp,((2**29)-1))
      upackf2(icomp) = rshift(icomp,29)+49

      do i = 50, 54
         do j = 1, 5
            k1 = ipack(i,j)
            k2 = ipack2(i,j)
            if (k1.ne.k2) print *, 'ipack(',i,',',j,') = ',k1,',',k2
            if (upackr(k1).ne.upackr2(k2))
     &         print *, 'upackr(',k1,') = ',upackr(k1),
     &                '; upackr2(',k2,') = ',upackr2(k2)
            if (upackf(k1).ne.upackf2(k2))
     &         print *, 'upackf(',k1,') = ',upackf(k1),
     &                '; upackf2(',k2,') = ',upackf2(k2),
     &                '; i,j = ',i,j
         end do
      end do

      end

