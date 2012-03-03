      Subroutine projec_FC(Coords, Hess, AtmMass, Grad, Hess_project, 
     &                     Work, Threshold, Nreals, Move_CMass, 
     &                     Proj_rots, Proj_grads)
C
      Implicit Double Precision (A-H, O-Z)
C
      Double Precision MO_Inertia(3,3)
      Logical Proj_rots, Proj_grads, Move_CMass
      Integer L_Index, R_Index
C
      Dimension Coords(3*Nreals), Grad(3*Nreals), AtmMass(Nreals), 
     &          Hess(3*Nreals, 3*Nreals), Hess_Project(3*Nreals, 
     &          3*Nreals), Work(3*Nreals, 3*Nreals), Asym_ten(3,3,3),
     &          CM(3), Center_Mass(3)
 
      Data Asym_ten/ 0.0D+00,  0.0D+00,  0.0D+00,
     &               0.0D+00,  0.0D+00, -1.0D+00,
     &               0.0D+00,  1.0D+00,  0.0D+00,
     &               0.0D+00,  0.0D+00,  1.0D+00,
     &               0.0D+00,  0.0D+00,  0.0D+00,
     &              -1.0D+00,  0.0D+00,  0.0D+00,
     &               0.0D+00, -1.0D+00,  0.0D+00,
     &               1.0D+00,  0.0D+00,  0.0D+00, 
     &               0.0D+00,  0.0D+00,  0.0D+00  /
C
C Normalize the gradinet vector. 
C
      If (Proj_grads) Then
         Grad_sqr = Ddot(3*Nreals, Grad, 1, Grad, 1)   
         If (Grad_sqr .GT. Threshold)  Call Normal(Grad, 3*Nreals)
      Endif
C

      Write(6,*)
      Write(6,"(a)") "@-Projec_FC At enetry the Coords and grads"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)
      Write(6,*)
      Write(6, "(3F17.13)") (Grad(i),i=1,3*Nreals)

C
C Move to the center of Mass coordinate system.
C
      If (Move_CMass) Then
         Call Dzero(CM, 3)
         TotMass = 0.0D0
         Do Iatoms = 1, Nreals
            Ioff = 3*(Iatoms - 1)
            TotMass = TotMass + AtmMass(Iatoms)
            Do Ixyz = 1, 3
               CM(Ixyz)   = Coords(Ioff + Ixyz)*AtmMass(Iatoms) +
     &                      CM(Ixyz)
            Enddo
         Enddo
C
         Do Ixyz = 1, 3
            If (TotMass .GT. Threshold) Then
                Center_Mass(IxYz) = CM(Ixyz)/TotMass
            Else
                Write(6, "(a)") "@-Project_FC Zero total mass"
                Call Errex
            Endif
         Enddo
C         
         Do Iatoms = 1, Nreals
            Ioff  = 3*(Iatoms - 1)
            Do Ixyz = 1, 3
               Coords(Ioff + Ixyz) = Coords(Ioff + Ixyz) - 
     &                               Center_Mass(Ixyz)
            Enddo
         Enddo
      Endif
C

      Write(6,*)
      Write(6,"(a)") "@-Projec_FC The center of mass Coords"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)

C
C compute the inverse of inertia matrix.
c 
      Call CalMOI(NReals, Coords, AtmMass, MO_Inertia)

      Write(6,*)
      Write(6,"(a)") "@-Projec_FC the moments of inertia matrix"
      Call output(MO_Inertia, 1, 3, 1, 3, 3, 3,1)


      Call Minv(MO_Inertia, 3, 3, Work, Det, 1.0D-8, 0, 1)
C
C Mass weigh the incomming Cartesians 
C
      Ioff = 1
      Do Iatoms = 1, Nreals
          Sqrtmass = Dsqrt(AtmMass(Iatoms))
          Do Ixyz = 1, 3
             Coords(Ioff) = Coords(Ioff)*Sqrtmass
             Ioff = Ioff + 1 
          Enddo
      Enddo    
C

      Write(6,*)
      Write(6,"(a)") "@-Projec_FC Mass weighted center of mass Coords"
      Write(6, "(3F17.13)") (Coords(i),i=1,3*Nreals)
      Write(6,*)
      Write(6,"(a)") "@-Projec_FC Inv. of the moms. of inertia matrix"
      Call output(MO_Inertia, 1, 3, 1, 3, 3, 3,1)

C
C Build the Hessian projector See, Miller, Handy and Adams, JCP, 
C 72, 99, (1980).
C
      Do Iatms = 1, Nreals
         Ioff = 3*(Iatms - 1)
CSSS         Itmp = Max(3*(Iatms - 1), 6*(Iatms -1)-3*Nreals)

         Do Jatms = 1, Iatms
          Joff = 3*(Jatms - 1)
CSSS          Jtmp = Max(3*(Jatms - 1), 6*(Jatms -1)-3*Nreals)





            Do Iz = 1, 3
               L_index = Ioff + Iz
CSSS               Itmp_L = Itmp + Iz
               Jz_end = 3
               If (Iatms .EQ. Jatms)  Jz_end = Iz
C
                 Do Jz = 1, Jz_end
                    R_index = Joff + Jz
CSSS                    Jtmp_R =  Jtmp + JZ
                    Sum = 0.0D0
C
                    If (Proj_rots) Then 

                        Do Ix = 1, 3
 
                           Do Iy = 1, 3
C
                              If (Asym_ten(Ix, Iy, Iz) .NE. 0.0D0) 
     &                            Then
                                   
                                   Do Jx = 1, 3

                                      Do Jy = 1, 3
                                         If (Asym_ten(Jx, Jy, Jz) .NE. 
     &                                       0.0D0) Then
                                             Sum = Sum + 
     &                                             Coords(Ioff + IY)*
     &                                             Coords(Joff + Jy)*
     &                                            MO_inertia(Ix, Jx)*
     &                                          Asym_ten(Ix, Iy, Iz)*
     &                                          Asym_ten(Jx, Jy, Jz)
                                         Endif
                                      Enddo
                                   Enddo 
                              Endif
                           Enddo
                        Enddo
                    Endif
C 
                    Hess_project(L_index, R_index) = Sum
C
                    If (Proj_grads) Hess_project(L_index, R_index) = 
     &                              Hess_project(L_index, R_index) + 
     &                              Grad(R_Index)*Grad(L_index)
C
                    If (Iz .EQ. Jz) Then



                        Hess_project(L_index, R_index) =
     &                  Hess_project(L_index, R_index) + 
     &                  Dsqrt(AtmMass(Iatms)*AtmMass(Jatms))/
     &                  Totmass
                    Endif
C
                 Enddo 
            Enddo 
C
         Enddo
      Enddo
C






C
C Build the projector (I - P)
C
      Do Jdeg = 1, 3*Nreals
          Do Ideg = 1, Jdeg
              Hess_project(Jdeg, Ideg) = -Hess_project(Jdeg, Ideg)
             If (Ideg .Eq. Jdeg)  Hess_project(Jdeg, Ideg) = 1.0D0 +
     &                            Hess_project(Jdeg, Ideg)   
             If (Dabs(Hess_project(Jdeg, Ideg)) .LT. Threshold) 
     &                            Hess_project(Jdeg, Ideg) = 0.0D0
              Hess_project(Ideg, Jdeg) = Hess_project(Jdeg, Ideg)
          Enddo
      Enddo
C 
C

      NX = 3*Nreals
      Write(6,*)
      Write(6,"(a)") "The Hessian projector:(I-P)"
      CALL OUTPUT(Hess_project, 1, NX, 1, Nx, Nx, Nx, 1)
      

C
C
C Project the Hessian (I - P)H(I - P)
C
      Call Xgemm("N", "N", 3*Nreals, 3*Nreals, 3*Nreals, 1.0D0, 
     &            Hess_project, 3*Nreals, Hess, 3*Nreals, 0.0D0, 
     &            Work, 3*Nreals)
C
      Call Xgemm("N", "N", 3*Nreals, 3*Nreals, 3*Nreals, 1.0D0, 
     &            Work, 3*Nreals, Hess_project, 3*Nreals, 0.0D0, 
     &            Hess, 3*Nreals)
C
      Return
      End
     
