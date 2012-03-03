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
      subroutine nuclear_nuclear_gradient(hess, nnhess, natoms)
c--------------------------------------------------------------------------
c   Calculates the nuclear-nuclear contribution to the gradient.
c--------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'hess.h'
      include 'dbugcom.h'

      integer i, j, k, natoms
      integer ihess, jhess, khess, iatom, jatom
      double precision x, y, z, r, r3, r5, pkl
      double precision hess(3*natoms,3*natoms)
      double precision NNhess(3*natoms,3*natoms)
      double precision d(natoms,natoms)

      if (natoms .gt. max_centers) then
         print *,'Error: natoms = ',natoms,' is larger than ',
     *      'max_centers ',max_centers
         call abort_job()
      endif
  
      DO I = 1, NATOMS  
C
         DO J = 1, NATOMS  

            if (i .ne. j) then 
C
            X = acenter(I,1) - acenter(J,1)          
            Y = acenter(I,2) - acenter(J,2)          
            Z = acenter(I,3) - acenter(J,3)          
            R = (X**2 + Y**2 + Z**2) 
            D(i,j) = -1.0d0/R**(1.5)  
c
            endif 
C
         ENDDO ! J 
C
      ENDDO ! I 
C
C     ----- FORM CONTRIBUTION TO GRADIENT -----
C
      DO J = 1, NATOMS 
         DO I = 1, 3
            NNgrad(I,J) = 0.0d0 
            DO K = 1, NATOMS 
               IF (J .ne. K) then 
                  PKL = (acenter(J,I)-acenter(K,I))
                  NNgrad(I,J) = NNgrad(I,J)+PKL*D(K,J)*charge(k)
     &                                                *charge(j) 
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 

      if (dbg) then
         do j = 1, natoms
            print *,'NNgrad for atom ',j,' : ',(nngrad(i,j),i=1,3)
         enddo
      endif
C
C Compute the Nuclear-Nuclear contribution to the Hessian.  
C -------------------------------------------------------- 
C
      DO I = 1, NATOMS  
C
         DO J = 1, NATOMS  

            if (i .ne. j) then 
C
            X = acenter(I,1) - acenter(J,1)          
            Y = acenter(I,2) - acenter(J,2)          
            Z = acenter(I,3) - acenter(J,3)          
            R = (X**2 + Y**2 + Z**2) 
            D(i,j) = 1.0d0/R**(2.5)  
c
            endif 
C
         ENDDO ! J 
C
      ENDDO ! I 
C
      DO IATOM = 1, NATOMS 
      DO JATOM = 1, NATOMS 
         DO I = 1, 3
         DO J = 1, 3
c
            ihess = (iatom-1)*3 + i 
            jhess = (jatom-1)*3 + j 
            hess(jhess,ihess) = 0.0d0 
c
         ENDDO ! J 
         ENDDO ! I 
      ENDDO ! JATOM 
      ENDDO ! IATOM 
C
      DO IATOM = 1, NATOMS 
      DO JATOM = 1, NATOMS 
         DO I = 1, 3
         DO J = 1, 3
c
            ihess = (iatom-1)*3 + i 
            khess = (iatom-1)*3 + j
            jhess = (jatom-1)*3 + j 
c 
            if ((i .ne. j) .and. (iatom .ne. jatom)) then 
c
               hess(jhess,ihess) = hess(jhess,ihess) 
     *                           - 3.0d0*charge(iatom)*charge(jatom)   
     *                           * (acenter(iatom,i) - acenter(jatom,i))  
     *                           * (acenter(iatom,j) - acenter(jatom,j))  
     *                           * D(iatom,jatom)  
c
               hess(khess,ihess) = hess(khess,ihess)
     *                           + 3.0d0*charge(iatom)*charge(jatom)
     *                           * (acenter(iatom,i) - acenter(jatom,i))
     *                           * (acenter(iatom,j) - acenter(jatom,j))
     *                           * D(iatom,jatom)
            endif 
c 
         ENDDO ! J 
         ENDDO ! I 
      ENDDO ! JATOM 
      ENDDO ! IATOM 
C
      DO IATOM = 1, NATOMS 
      DO JATOM = 1, NATOMS 
C
         if (iatom .ne. jatom) then
            X = acenter(IATOM,1) - acenter(JATOM,1)          
            Y = acenter(IATOM,2) - acenter(JATOM,2)          
            Z = acenter(IATOM,3) - acenter(JATOM,3)          
            R = (X**2 + Y**2 + Z**2) 
            R5 = 1.0d0/R**(2.5)  
            R3 = 1.0d0/R**(1.5)  
c
            DO I = 1, 3
            DO J = 1, 3
               ihess = (iatom-1)*3 + i 
               jhess = (jatom-1)*3 + j 
               if (i .eq. j) then 
                  hess(jhess,ihess) = hess(jhess,ihess) 
     *                    - 3.0d0*charge(iatom)*charge(jatom)   
     *                    * (acenter(iatom,i) - acenter(jatom,i))  
     *                    * (acenter(iatom,j) - acenter(jatom,j))  
     *                    * R5    
     *                    + charge(iatom)*charge(jatom)*R3   

                  hess(ihess,ihess) = hess(ihess,ihess) 
     *                    + 3.0d0*charge(iatom)*charge(jatom)   
     *                    * (acenter(iatom,i) - acenter(jatom,i))  
     *                    * (acenter(iatom,j) - acenter(jatom,j))  
     *                    * R5    
     *                    - charge(iatom)*charge(jatom)*R3   
               endif 
            ENDDO ! J 
            ENDDO ! I 
         endif ! iatom .ne. jatom
      ENDDO ! JATOM 
      ENDDO ! IATOM 

c--------------------------------------------------------------------------
c   Write out Hessian data
c--------------------------------------------------------------------------
 
      if (dbg) then
        write(6,*) ' '
        write(6,*) ' Nuclear-Nuclear Hessian data '

        do i = 1, NATOMS  
        do j = 1, NATOMS  

c          do ix = 1, 3
c          do jx = 1, 3

c             ihess = (i-1)*3 + ix
c             jhess = (j-1)*3 + jx

              write(6,33) hess((i-1)*3+1,(j-1)*3+1),
     *                   hess((i-1)*3+1,(j-1)*3+2),
     *                   hess((i-1)*3+1,(j-1)*3+3)

              write(6,33) hess((i-1)*3+2,(j-1)*3+1),
     *                   hess((i-1)*3+2,(j-1)*3+2),
     *                   hess((i-1)*3+2,(j-1)*3+3)

              write(6,33) hess((i-1)*3+3,(j-1)*3+1),
     *                   hess((i-1)*3+3,(j-1)*3+2),
     *                   hess((i-1)*3+3,(j-1)*3+3)

c          enddo ! jx = 1, 3
c          enddo ! ix = 1, 3

        enddo ! j
        enddo ! i
      endif

33      format(3F16.8) 
c
c Zero out the Hessian so then the N-N contribution is not included 
C
      DO IATOM = 1, NATOMS 
      DO JATOM = 1, NATOMS 
         DO I = 1, 3
         DO J = 1, 3
 
            ihess = (iatom-1)*3 + i 
            jhess = (jatom-1)*3 + j 
            NNhess(ihess,jhess) = hess(ihess,jhess)
            hess(ihess,jhess) = 0.0d0 
 
         ENDDO ! J 
         ENDDO ! I 
      ENDDO ! JATOM 
      ENDDO ! IATOM 

      return
      end
