      Subroutine Ln_interpol(Coord_K1P, Coord_K1CM, Coord_K1C, 
     &                       Coord_PV, Grad_on_K1P, Grad_on_K1C, 
     &                       Vec_K1P, Vec_K1C, Grad_K1PM, 
     &                       Grad_K1CM, Grad_K1C, Vec_K1C_K1P, 
     &                       Vec_K1Int, Grad_K1Int, AtmMass, 
     &                       Nreals, Ln_Intrp_Tol)

      Implicit Double Precision (A-H, O-Z)
C    
      Double Precision Ln_intrp_Tol, K1P_K1C_Norm, Pang_K1C_K1P
C
      Dimension Coord_K1P(3*Nreals), Coord_K1CM(3*Nreals), 
     &          Grad_on_K1P(3*Nreals), Grad_on_K1C(3*Nreals),
     &          Vec_K1P(3*Nreals), Vec_K1C(3*Nreals),
     &          Vec_K1C_K1P(3*Nreals), Vec_K1Int(3*Nreals),
     &          Grad_K1Int(3*Nreals), Coord_PV(3*Nreals),
     &          Grad_K1CM(3*Nreals), Grad_K1C(3*Nreals),
     &          Grad_K1PM(3*Nreals), Coord_K1C(3*Nreals), 
     &          AtmMass(Nreals)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") "@-Ln_interpol M.W. prev., "
     &                  "current point, N.M.W. current point"
      Write(6,"(a)")  "The prv. point"
      Write(6, "(3F17.13)") (Coord_K1P(i),i=1,3*Nreals)
      Write(6,"(a)")  "The curr. point"
      Write(6, "(3F17.13)") (Coord_K1CM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The curr. point (no. mass weight)"
      Write(6, "(3F17.13)") (Coord_K1C(i),i=1,3*Nreals)
      Write(6,"(a)")  "The M.W. pivot point"
      Write(6, "(3F17.13)") (Coord_PV(i),i=1,3*Nreals)
      Write(6,*) 
      Write(6, "(a,a)") "@-Ln_interpol M.W. grad and vec on, "
     &                  "K1 vector"
      Write(6,"(a)")  "The Grad_on_K1P"
      Write(6, "(3F17.13)") (Grad_on_K1P(i),i=1,3*Nreals)
      Write(6,"(a)")  "Grad_on_K1C"
      Write(6, "(3F17.13)") (Grad_on_K1C(i),i=1,3*Nreals)
      Write(6,*) 
      Write(6, "(a,a)") "@-Ln_interpol M.W. prv. and curr. grad "
      Write(6,"(a)")  "The M.W previous gradient"
      Write(6, "(3F17.13)") (Grad_K1PM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The M.W current gradient"
      Write(6, "(3F17.13)") (Grad_K1CM(i),i=1,3*Nreals)
#endif                        
C
      Do Iatom = 1, Nreals
         Ioff = 1 + 3*(Iatom - 1)
         Call Vec(Coord_PV(Ioff), Coord_K1CM(Ioff), Vec_K1C(Ioff), 0)
         Call Vec(Coord_PV(Ioff), Coord_K1P(Ioff), Vec_K1P(Ioff), 0)
      Enddo
C
      Vec_K1C_Norm  = Ddot(3*Nreals, Vec_K1C, 1, Vec_K1C, 1)
      Vec_K1P_Norm  = Ddot(3*Nreals, Vec_K1P, 1, Vec_K1P, 1)
      K1P_K1C_Norm  = Ddot(3*Nreals, Vec_K1P, 1, Vec_K1C, 1)
C
      Pang_K1C_K1P = (K1P_K1C_Norm/(Dsqrt(Vec_K1C_Norm)*
     &                              Dsqrt(Vec_K1P_Norm)))
      If (Abs(Pang_K1C_K1P) .Gt. 1.0D0) Pang_K1C_K1P = Pang_K1C_K1P/
     &                                             Abs(Pang_K1C_K1P)
      PAng_K1C_K1P = DAcos(Pang_K1C_K1P)
      Theta_prime  = PAng_K1C_K1P

C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The norms of the vectors for theta prime"
      write(6,"(3F17.13)") dsqrt(Vec_K1P_Norm), (K1P_K1C_Norm),
     &                    dsqrt(Vec_K1C_Norm )
      Write(6,*)
      Write(6,"(a,F17.13)") "The Theta_prime", Theta_prime
#endif 
C
      Do Iatom = 1, Nreals
         Ioff = 1 + 3*(Iatom - 1)
         Call Vec(Vec_K1P(Ioff),  Vec_K1C(Ioff),  Vec_K1C_K1P(Ioff),
     &            0)
      Enddo

      G1_prime = Ddot(3*Nreals, Vec_K1C_K1P, 1, Grad_on_K1C, 1)
      G2_Prime = Ddot(3*Nreals, Vec_K1C_K1P, 1, Grad_on_K1P, 1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The norms of the vectors for theta "
      write(6,"(2F17.13)") G1_prime, G2_Prime
#endif
C
      Theta = Theta_prime*G2_prime/(G2_Prime - G1_Prime)
#ifdef  _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a, F17.13)") "The theta", theta
#endif 
C
      If (Theta .LT. 0.0D0 .AND. Theta_prime .LT. Ln_Intrp_Tol) Then
         Write(6, "(a)") "@-Ln_interpol, Interpolation is invlaid"
         Return
      Endif

      If (DAbs(Theta) .GT. Ln_Intrp_Tol .AND. DAbs(Theta_Prime 
     &    - Theta) .GT. Ln_Intrp_Tol) Then
         Write(6, "(a)") "@-Ln_interpol, Interpolation is invlaid"
         Return
      Endif
C
      Ratio_1 = Theta/Theta_prime
      Ratio_2 = Sin(Theta)/Sin(Theta_prime)
      Ratio_3 = Cos(Theta) - Sin(Theta)*Cos(Theta_prime)/
     &                                  Sin(Theta_prime)
C
      Do Ideg = 1, Nreals
         Ioff = 3*(Ideg - 1) 
         Do Ixyz = 1, 3 
            Ioff = Ioff + 1 
            Vec_K1Int(Ioff)  = Ratio_2*Vec_K1C(Ioff) + Ratio_3*
     &                         Vec_K1P(Ioff)
            Grad_K1Int(Ioff) = Ratio_1*Grad_K1CM(Ioff) + (1.0D0 - 
     &                         Ratio_1)*Grad_K1PM(Ioff)
         Enddo
      Enddo 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "@-ln_interpol, geo. and gradient updates"
      Write(6,"(a)")  "The M.W. geo. update"
      Write(6, "(3F17.13)") (Vec_K1int(i),i=1,3*Nreals)
      Write(6,"(a)")  "The M.W. gradient update" 
      Write(6, "(3F17.13)") (Grad_K1int(i),i=1,3*Nreals)
      Write(6,*) 
#endif
C
      Call Daxpy(3*Nreals, 1.0D0, Coord_PV, 1, Vec_K1Int, 1)
      Call Dcopy(3*Nreals, Vec_K1Int, 1, Coord_K1CM, 1)
      Call Dcopy(3*Nreals, Grad_K1Int, 1, Grad_K1CM, 1)
C
      Ioff = 0
      Do Iatom = 1, Nreals
         Do Ixyz = 1, 3
           Ioff = Ioff + 1
            Grad_K1Int(Ioff) = Grad_K1Int(Ioff)*
     &                         Dsqrt(AtmMass(Iatom))
            Vec_K1Int(Ioff)  = Vec_K1Int(Ioff)/
     &                         Dsqrt(AtmMass(Iatom))
         Enddo
      Enddo
C
      Call Dcopy(3*Nreals, Vec_K1Int, 1, Coord_K1C, 1)
      Call Dcopy(3*Nreals, Grad_K1Int, 1, Grad_K1C, 1)
      
      Return
      End
       


