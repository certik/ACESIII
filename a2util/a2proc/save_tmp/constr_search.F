#include <flags.h>
      Subroutine Constr_search(Coords_K0OM, Coords_K1C, Coords_K1PM,
     &                         Vcoords, Coords, Coords_K1CM, Grad_K1C,
     &                         Grad_K1CM, Grad_K0O, A2Hess, Hess_K1C,
     &                         AtmMass, Vec_K1C0, Vec_K1C1, Vec_K0O, 
     &                         Grad_on_K0O, Grad_on_K1C, Work,
     &                         B2ang, Units, Stride, Delta, 
     &                         First_step, New_IRC, Constr_Min,
     &                         Imap, Ncycles, Nreals, Natoms,
     &                         Opt_Tol, Stride_Total, Stride_Tol,
     &                         Ln_Intrp_Tol, Hess_Eval_Tol,
     &                         Lamsearch_Tol, Binsearch_Tol, 
     &                         G_Cutoff, R_Cutoff, E_currnt)
C
      Implicit Double Precision (A-H, O-Z)
C
#include <machsp.com>
#include <jodaflags.com>
C
      Logical Converged, First_step, Powel_upd, Bfgs_upd, Bohr, 
     &        New_IRC
C
      Double Precision Lamsearch_Tol, Ln_Intrp_Tol
C
      Character*5 Constr_Min
      Character*4 Units

      Dimension Coords_K0OM(3*Nreals), Coords_K1C(3*Nreals),
     &          Coords_K1PM(3*Nreals), Grad_K1C(3,Natoms), 
     &          Imap(Natoms), VEc_K1C_stat(6), Grad_K1CM(3,Natoms), 
     &          Grad_K1C_stat(6), AtmMass(Nreals), Vec_K1C0(3*Nreals),
     &          Grad_on_K1C(3, Nreals), Grad_on_K1C_stat(6),
     &          A2Hess(3*Natoms, 3*Natoms), Coords_K1CM(3*Nreals),
     &          Hess_K1C(3*Natoms, 3*Natoms), Work(36*Nreals*Nreals),
     &          Vcoords(3*Natoms), Coords(3*Natoms), 
     &          Vec_K0O(3*Nreals), Grad_on_K0O(3*Nreals), 
     &          Grad_K0O(3*Nreals), Vec_K1C1(3*Nreals)
C
C Compute the gradients at the current point
C
      Convert = B2Ang
      If (Units .EQ. "Bohr") Convert = 1.0D0
C
      Write(6,*) "First_step, New_IRC", First_step, New_IRC, 
     &            Constr_Min
      If (First_step) Then
C
C Except the reference point, all the subsequent points used
C Cartesian coordinates. The reference point coordiante is choosen
C by the user and can be in internal or Cartesians. If internal
C there could be dummy atoms, and any point other than the reference
C point, the dummy atoms have (0.0,0.0,0.0) coordinates and do not
C enter into the calculation. With this mechanism one can eliminate
C the JOBACR IO error resulting from trying to overwrite an
C exsisting record with record of different length.
C
         If (Iflags(54) .NE. 0)  Iflags(54)   = 0
                                 Iflags2(138) = 1
         If (Iflags(3) .EQ. 2)
     &                           Iflags(3) = 1
C
          Call Dzero(Vcoords, 3*Natoms)
          Call Dzero(Coords, 3*Natoms)
          Call Daxpy(3*Nreals, Convert, Coords_K1C, 1, Vcoords, 1)
          Do Iatoms = 1, Natoms
             If (Imap(Iatoms) .NE. 0) Then
                 Ioff = 3*(Imap(Iatoms) - 1) + 1
                 Joff = 3*(Iatoms - 1) + 1
                 Call Dcopy(3, Vcoords(Joff), 1, Coords(Ioff), 1)
             Endif
          Enddo 
          Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP,
     &                Coords)
          Call Putrec(1,'JOBARC','IFLAGS  ',   100, iflags)
          Call Putrec(1,'JOBARC','IFLAGS2 ',  500, iflags2)
C
          Call aces_ja_fin
C
          Call Runit("runaces2b")
C
          Call aces_ja_init

      Else if (.NOT. First_step .AND. New_IRC) Then
C
          Call Dzero(Vcoords, 3*Natoms)
          Call Dzero(Coords, 3*Natoms)
          Call Daxpy(3*Nreals, Convert, Coords_K1C, 1, Vcoords, 1)
          Do Iatoms = 1, Natoms
             If (Imap(Iatoms) .NE. 0) Then
                 Ioff = 3*(Imap(Iatoms) - 1) + 1
                 Joff = 3*(Iatoms - 1) + 1
                 Call Dcopy(3, Vcoords(Joff), 1, Coords(Ioff), 1)
             Endif
          Enddo
          Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP,
     &                Coords)
C
C
          Call a2_reset_jarc
          Call aces_ja_fin
C
          Call Runit("runaces2b")
C
          Call aces_ja_init

      Endif

      Call Getrec(20,'JOBARC','GRADIENT',3*Nreals*IINTFP,
     &               Grad_K1C(1,1))
      Call Getrec(20,'JOBARC','TOTENERG', IINTFP, E_previs)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search Grad_k1C(comp. order)"
      Write(6, "(3F12.8)") ((Grad_K1C(i,j), i=1,3),j=1,Nreals)
      Write(6,*) "@-Constr_search The imap array"
      Write(6, "(5I2)") (Imap(i), i=1, Nreals)
#endif
C
C First convert to ZMAT order and then remove the dummy atoms.
C (Note that the gradient record does not have dummy atom 
C contributions, but in order to map to ZMAT order, we have to
C add dummy contributions (zero) again and then remove).

      If (First_Step .OR. New_IRC) Then
          Call Dzero(grad_K1CM, 3*Natoms)
          Do i=1, Natoms
             k=Imap(i)
             If (k .NE. 0) Then
                Do j=1,3
                   grad_K1CM(j,k)  = Grad_K1C(j,i)
                Enddo
             Endif
          Enddo
C
          Call Dcopy(3*Natoms, grad_K1CM, 1, Grad_K1C, 1)
C
         Do i=1, Natoms
            If (Imap(i) .NE. 0) Then
                Call Dcopy(3, Grad_K1C(1,Imap(i)), 1, 
     &                     grad_K1CM(1,i), 1)
            Endif
         Enddo
         Call Dcopy(3*Nreals, grad_K1CM, 1, Grad_K1C, 1)
C
C During the first step save the gradient of the current point
C in the location reserved for the previous point gradient.
C
         Call Dcopy(3*Nreals, Grad_K1C, 1, Grad_K0O, 1)
      Endif
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search Grad_K1C (ZMAT order)"
      Write(6, "(3F12.8)") ((Grad_K1C(i,j), i=1,3),j=1,Natoms)
#endif
C
C Grad_k1c_stat(1-6) are the largest absolute, the smallest absolute,
C the largest, the largest, the smallest, the norm, and the
C dynamic range
C
      Call Vstat(Grad_K1C, Grad_K1C_stat, 3*Nreals)
      If (Grad_K1C_stat(1) .LT. Opt_Tol .AND. Grad_K1C_stat(5) .LT.
     &    Opt_Tol/3.0D0 .AND. Str_Total .GT. 1.0D0 ) Then
          Converged = .True.
      Else
          Converged = .False.
      Endif
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search Grad_K1C stats"
      Write(6, "(a,2F12.8)") "largest and the norm",Grad_K1C_stat(1),
     &                        Grad_K1C_stat(5)
#endif
C
C      If (First_step .AND. .NOT. New_IRC) Then       
      If (First_step .OR. New_IRC) Then       
         Ioff = 0
         Do Iatom = 1, Nreals
            Do Ixyz = 1, 3
               Grad_K1CM(Ixyz, Iatom) = Grad_K1CM(Ixyz, Iatom)/
     &                                  Dsqrt(AtmMass(Iatom))
               Ioff = Ioff + 1
               Coords_K0OM(Ioff) = Coords_K0OM(Ioff)*
     &                             Dsqrt(AtmMass(Iatom))
               Coords_K1CM(Ioff)= Coords_K1CM(Ioff)*
     &                            Dsqrt(AtmMass(Iatom))
               Coords_K1PM(Ioff) = Coords_K1PM(Ioff)*
     &                             Dsqrt(AtmMass(Iatom))
            Enddo
         Enddo
      Endif

C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") "@-Constr_search Mass weighted grads, prev."
     &                  "current and pivot point"
      Write(6, "(3F12.8)") ((Grad_K1CM(i,j), i=1,3),j=1,Nreals)
      Write(6,"(a)")  "The prv. point"
      Write(6, "(3F12.8)") (Coords_K0OM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The curr. point"
      Write(6, "(3F12.8)") (Coords_K1CM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The curr. point (no. mass weight)"
      Write(6, "(3F12.8)") (Coords_K1C(i),i=1,3*Nreals)
      Write(6,"(a)")  "The pivot. point"
      Write(6, "(3F12.8)") (Coords_K1PM(i),i=1,3*Nreals)
#endif
C
C Define the vectors k+1: (qk+1 - qk)
C 
      Do Iatom = 1, Nreals
         Ioff = 1 + 3*(Iatom - 1) 
         Call Vec(Coords_K1PM(Ioff),  Coords_K1CM(Ioff), 
     &            Vec_K1C0(Ioff), 0)
      Enddo
C
      Call Vstat(Vec_K1C0, Vec_K1C_Stat, 3*Nreals)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search Vec_K1C0"
      Write(6, "(3F12.8)") (VEC_K1C0(i),i=1,3*Nreals)
      Write(6,*) "@-Constr_search Vec_K1C_stat"
      Write(6, "(F12.8)") Vec_K1C_stat(5) 
#endif
C
      Vec_K1C_length = Vec_K1C_Stat(5)*Dsqrt(Dble(3*Nreals))
C
      If (Abs(Vec_K1C_length - 0.50D0*Stride) .GT. Stride_Tol) Then
         Write(6, "(a)") "The Norm criteria is violated"
         Call Errex
      Endif
C
      Call Dcopy(3*Nreals, Vec_K1C0, 1, Vec_K1C1, 1)
      Call Dscal(3*Nreals, 1.0D0/(Vec_K1C_Stat(5)*Dsqrt
     &           (Dble(3*Nreals))), Vec_K1C1, 1)
C
C Gradient along the k+1 vectors.
C
      Grad_Proj_K1C =  Ddot(3*Nreals, Vec_K1C1, 1, Grad_K1CM, 1)
C
      Do Iatoms = 1, Nreals
         Ioff = 3*(Iatoms - 1)
         Do Ixyz = 1, 3 
            Ioff = Ioff + 1
            Grad_on_K1C(Ixyz, Iatoms) = Grad_K1CM(Ixyz, Iatoms) - 
     &                                  Vec_K1C1(Ioff)*Grad_Proj_K1C
         Enddo
      Enddo
C               
      Call Vstat(Grad_on_K1C, Grad_on_K1C_stat, 3*Nreals)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search Vec_K1C0 and Grad_on_K1C and state"
      Write(6,*) "Vec_K1C0"
      Write(6, "(3F12.8)") (Vec_K1C0(i),i=1,3*Nreals) 
      Write(6,*)
      Write(6, "(3F12.8)") ((Grad_on_K1C(i,j), i=1,3),j=1,Nreals)
      Write(6, "(a,2F12.8)") "largest and the norm of Grad_on_K1C",
     &                    Grad_on_K1C_stat(1),
     &                    Grad_on_K1C_stat(5)*Dsqrt(Dble(3*Nreals))
#endif
      If (Ncycles .ge. 2) Then
         Powel_upd = .true.
         Call Dcopy(3*Nreals, Grad_K0O, 1, Work, 1)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") "@-Constr_search entry to Hess. update"
     &                 "current, previous gradinets and geo. update"
      Write(6,*) "The current gradient"
      Write(6, "(3F12.8)") ((Grad_K1C(i,j), i=1,3),j=1,Nreals)
      Write(6,*) "The previous gradient"
      Write(6, "(3F12.8)") (Work(j), j=1,3*Nreals)
      Write(6,*) "The geo. update"
      Write(6,"(3F12.8)") (Work(9*Nreals+i), i=1,3*Nreals)
      Write(6,*) "The restored Hessian"
      Call output(Hess_K1C, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals,
     &            3*Nreals, 1)
#endif
         Iw_off_pg = 1
         Iw_off_dx = Iw_off_pg + 9*Nreals 
         Iw_off_s1 = Iw_off_dx + 12*Nreals  
         Iw_off_s2 = Iw_off_s1 + 12*Nreals + Nreals*Nreals
         Iw_off_ed = Iw_off_s2 + 3*Nreals*Nreals
         If (Iw_off_ed .GT. 36*Nreals*Nreals) Call Insmem(
     &       "@-Constr_search",Iw_off_ed,36*Nreals*Nreals)
C
         If (Powel_upd) Call Powel_update(Grad_K1C, Work(Iw_off_pg),
     &                                    Hess_K1C, Work(Iw_off_s1),
     &                                    Work(Iw_off_dx),   
     &                                    Work(Iw_off_s2), 1.0D0,
     &                                    3*Nreals)
C
#ifdef _DEBUG_LVL0
      if (Powel_upd) then
      Write(6,*) "The powel updated Hessian"
      Call output(Hess_K1C, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals,
     &            3*Nreals, 1)
      endif
#endif
C
         If (Bfgs_upd) Call Bfgs_update(Grad_K1C, Work(Iw_off_pg), 
     &                      Hess_K1C, Work(Iw_off_s1), 
     &                      Work(Iw_off_dx), Work(Iw_off_s2),
     &                      3*Nreals)
C
#ifdef _DEBUG_LVL0
      if (Bfgs_upd) then
      Write(6,*) "The BFGS updated Hessian"
      Call output(Hess_K1C, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals,
     &            3*Nreals, 1)
      endif
#endif   
C
         Call Ln_interpol(Coords_K0OM, Coords_K1CM, Coords_K1C,
     &                    Coords_K1PM, Grad_on_K0O, Grad_on_K1C, 
     &                    Vec_K0O, Vec_K1C, Grad_K0O, Grad_K1CM, 
     &                    Grad_K1C, Work, Work(3*Nreals+1), 
     &                    Work(6*Nreals+1), AtmMass, Nreals, 
     &                    Ln_Intrp_Tol)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a,a,a)") "@-constr_search After linear interpolation",
     &                 " current point and gradients, M.W. first",
     &                 " and then no M.W"
      Write(6,"(a)")  "The M.W. current geometry " 
      Write(6, "(3F12.8)") (Coords_K1CM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The M.W.current gradient"
      Write(6,"(3F12.8)") ((Grad_K1CM(i,j), i=1,3),j=1,Nreals)
      Write(6,"(a)")  "The no M.W. current geometry " 
      Write(6, "(3F12.8)") (Coords_K1C(i),i=1,3*Nreals)
      Write(6,"(a)")  "The no M.W.current gradient"
      Write(6,"(3F12.8)") ((Grad_K1C(i,j), i=1,3),j=1,Nreals)
#endif
      Endif
C
C Save the current coordinates in the work array. The current coordinates
C will get modified in get_lambda. The unmodified current coordinates are
C needed in subsequent steps. 
      
      Call Dcopy(3*Nreals, Coords_K1C, 1, VCoords, 1)
      Call Get_lambda(Grad_K1CM, Vec_K1C0, Work, Work(3*Nreals+1),  
     &                Hess_K1C, A2Hess, AtmMass, Work(6*Nreals+1), 
     &                Work(9*Nreals + 1), Coords_K1C, Stride, Delta, 
     &                Nreals, Hess_Eval_Tol, Lamsearch_Tol, 
     &                Binsearch_Tol)
      Call Dcopy(3*Nreals, Vcoords, 1, Coords_K1C, 1)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") "@-Constr_search The No M. W. updated vector",
     &                 " (Vec_K1C_Updated)"
      Write(6,"(3F12.8)") (Work(6*Nreals+i), i=1,3*Nreals)
      Write(6,*) "@-Constr_search geo. (before updated, no M.W)"
      Write(6,"(3F12.8)") (Coords_K1C(i), i=1,3*Nreals)
#endif
C
C Save the previous point coordinates, gradients and tangent gradients before
C generating the new point.
C
      Call Dcopy(3*Nreals, Coords_K1CM, 1, Coords_K0OM, 1)
      Call Dcopy(3*Nreals, Grad_on_K1C, 1, Grad_on_K0O, 1)
      Call Dcopy(3*Nreals, Grad_K1C, 1, Grad_K0O, 1)
C
C Add the geo. update to the current point and generate a new point.

      Call Daxpy(3*Nreals, 1.0D0, Work(6*Nreals+1), 1, Coords_K1C, 1)
      Call Dcopy(3*Nreals, Work(6*Nreals+1), 1, Work(9*Nreals+1), 1)
C
      Ioff = 0
      Do Iatom = 1, Nreals
         Do Ixyz = 1, 3
               Ioff = Ioff + 1
                Work(6*Nreals+Ioff)= Work(6*Nreals+Ioff)*
     &                               Dsqrt(AtmMass(Iatom))
         Enddo
      Enddo
      Call Daxpy(3*Nreals, 1.0D0, Work(6*Nreals+1), 1, Coords_K1CM, 1)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search The updated no M.W. geo."
      Write(6,"(3F12.8)") (Coords_K1C(i), i=1,3*Nreals)
      Write(6,*) "@-Constr_search The updated M.W. geo."
      Write(6,"(3F12.8)") (Coords_K1CM(i), i=1,3*Nreals)
#endif
C 
      Delta_q = Ddot(3*Nreals, Work(3*Nreals+1), 1, Work(3*Nreals+1),
     &               1)
      Delta_g = Grad_on_K1C_stat(5)*Dsqrt(Dble(3*Nreals))
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Check for convergence"
      Write(6,"(a,2F12.8)") "G_cutoff, R_cutoff",G_cutoff, R_cutoff
      Write(6,"(a,2F12.8)") "Delta_q, Delta_g  ",Delta_q, Delta_g
#endif
      If (Delta_g .LT. G_cutoff .AND. Delta_q .LT. R_cutoff) 
     &   Constr_Min = "Found"       
C
      If (Constr_min .EQ. "Found") Then

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search The IRC is converged "
      Write(*,"(3F12.8)") (Coords_K1C(I), I = 1, 3*Nreals)
#endif
C
         IUnitO = 30
         Open(IUnitO,File='IRC.map',Status='Unknown')
         ierr = 0
         Do while (ierr .eq. 0)
            Read(IUnitO,7000, iostat=ierr) Temp_buff
         Enddo
C
 7000 Format(5X,A5,F15.10,2F18.10,F18.10)
C
         Write(IUnitO,"(3(2X,F12.8))") (Coords_K1C(I), I = 1,3*Nreals)
         Write(IUnitO, "")
         Close(IUnitO)
         Return
      Endif
C
      Call Dzero(Vcoords, 3*Natoms)
      Call Dzero(Coords, 3*Natoms)
      Call Daxpy(3*Nreals, Convert, Coords_K1C, 1, Vcoords, 1)
      Do Iatoms = 1, Natoms
         If (Imap(Iatoms) .NE. 0) Then
             Ioff = 3*(Imap(Iatoms) - 1) + 1
             Joff = 3*(Iatoms - 1) + 1
             Call Dcopy(3, Vcoords(Joff), 1, Coords(Ioff), 1)
         Endif
      Enddo
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-Constr_search updated geo. to JARC"
      Write(6,"(3F12.8)") (Coords(i), i=1,3*Nreals)
#endif    
      Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP,
     &            Coords)
C
      Call a2_reset_jarc
      Call aces_ja_fin
      Call Runit("runaces2b")
      Call aces_ja_init

      Call Getrec(20,'JOBARC','GRADIENT',3*Nreals*IINTFP,
     &            Grad_K1C(1,1))
      Call Getrec(20,'JOBARC','TOTENERG', IINTFP, E_currnt)

      If (E_currnt .GT. E_previs) Then
          Write(6, "(a,a)") "Warning: The energy is rising ", 
     &                      "during the constrained search." 
           E_Previs = E_currnt
      Endif 
C
      Call Dzero(grad_K1CM, 3*Natoms)
      Ioff = 0
      Do i=1, Natoms
         k=Imap(i)
         If (k .NE. 0) Then
            Do j=1,3
               Grad_K1CM(j,k)  = Grad_K1C(j,i)
            Enddo
         Endif
      Enddo
C
      Call Dcopy(3*Natoms, Grad_K1CM, 1, Grad_K1C, 1)
C
      Do i=1, Natoms
         If (Imap(i) .NE. 0) Then
            Call Dcopy(3, Grad_K1C(1,Imap(i)), 1, Grad_K1CM(1,i), 1)
         Endif
      Enddo
C
      Call Dcopy(3*Nreals, grad_K1CM, 1, Grad_K1C, 1)
C
      Do Iatom = 1, Nreals 
         Do Ixyz = 1, 3
            Grad_K1CM(Ixyz, Iatom) = Grad_K1CM(Ixyz, Iatom)/
     &                               Dsqrt(AtmMass(Iatom))
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") "@-Constr_search Grad_K1C at the updated",
     &                 " M. W. gradient. (ZMAT order)."
      Write(6, "(3F12.8)") ((Grad_K1CM(i,j), i=1,3),j=1,Nreals)
#endif
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") "@-Constr_search leaving M. W.,"
     &                  "prev current and pivot point"
      Write(6, "(3F12.8)") (Coords_K0OM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The curr. point"
      Write(6, "(3F12.8)") (Coords_K1CM(i),i=1,3*Nreals)
      Write(6,"(a)")  "The pivot. point"
      Write(6, "(3F12.8)") (Coords_K1PM(i),i=1,3*Nreals)
#endif
C  
      Return
      End 
