      SUBROUTINE GET_REFVIB_DATA(Vcoords, Norm_coords, Vomega, Omega,
     &                           AtmMass, AtmLabel, SymLabel, Coords,   
     &                           Btmp, A2grad, Grad, Hess, Vhess,  
     &                           Imap, Vib_Type, Nreals, Natoms, 
     &                           Nvibs, B2ang, Au2Invcm, Trns_state,
     &                           Mass_Weigh_nm, Mass_Weigh_gr,
     &                           Get_hess, Get_Grad)

      Implicit Double Precision (A-H, O-Z)
C
#include <machsp.com> 
C
      Character*4 SymLabel(3*Natoms)
      Character*5 AtmLabel(Natoms)
      Character*11 Vib_Type(3*Natoms)
      Double Precision Norm_coords
C
      Logical Trns_state, Mass_weigh_nm, Get_Hess, Get_Grad, 
     &        Mass_Weigh_gr
C
      Dimension VCoords(3*Natoms), Norm_coords(9*Nreals*Nreals), 
     &          Vomega(3*Nreals), Omega(3*Nreals), AtmMass(Natoms), 
     &          Coords(3*Nreals), Btmp(Natoms), Imap(Natoms),
     &          Hess(3*Natoms,3*Natoms), Vhess(9*Natoms*Natoms),
     &          Grad(3,Natoms), A2grad(3,Natoms)

C       
      Call Getrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP,
     &            Vcoords)
      Call Getrec(20, "JOBARC", "NORMCORD", 3*Nreals*3*Nreals*
     &             IINTFP, Norm_coords)
      Call Getrec(20, "JOBARC", 'FORCECON', 3*Nreals*IINTFP,
     &            VOmega)
      Call Getcrec(20, "JOBARC", 'VIB_TYPE', 3*Nreals*11*IINTFP, 
     &             Vib_Type)
C
      Ivibs = 0
      Do Imode = 1, 3*Nreals
C
         If (Vomega(Imode) .lt. 0.0D0) Then
            Vomega(Imode) = - DSQRT(DABS(Vomega(Imode)))*Au2Invcm
         Else
            Vomega(Imode) = DSQRT(Vomega(Imode))*Au2Invcm
         Endif
C
         If (Vib_Type(Imode) .EQ. ' VIBRATION ') Then
             Ivibs = Ivibs + 1
             Omega(Ivibs) = Vomega(Imode)
         Endif
C
      Enddo
C
      If (Trns_state .AND. Omega(1) .GE. 0.0D0) Then
          Write(6, "(a)") "The IRC search must start from the TS."
CSSS          Call Errex
      Endif
C
      Call Getcrec(20, "JOBARC", "VIB_SYMS", 3*Nreals*4, SymLabel)
C
      Call Dcopy(Natoms, AtmMass, 1, Btmp, 1)
      Ireal = 0
      Icord = 1
      Jcord = 1
      Do Iatom = 1, Natoms
         Icord = Icord + (Iatom - 1)*3
         If (.Not. (Btmp(Iatom) .lt. 0.50D0)) Then
            Ireal = Ireal + 1
            Jcord = Jcord + (Ireal - 1)*3
            AtmMass(Ireal) = Btmp(Iatom)
            Call Dcopy(3, Vcoords(Icord), 1, Coords(Jcord), 1)
         Endif
         Icord = 1
         Jcord = 1
      Enddo
C
      Ivibs = 0
      Do Imode = 1, 3*Nreals
         If (Vib_Type(Imode) .EQ. ' VIBRATION ') Then
            Ivibs = Ivibs + 1
            SymLabel(Ivibs) = SymLabel(Imode)
         Endif
      Enddo
C
      Ivibs = 0
      Do Imodes = 1, 3*Nreals
C
         If (Vib_Type(Imodes) .EQ. ' VIBRATION ') Then 
            Ivibs = Ivibs + 1
            Joff  = 3*Nreals*(Ivibs  - 1) + 1
            Ioff  = 3*Nreals*(Imodes - 1) + 1
            Call Dcopy(3*Nreals, Norm_coords(Ioff), 1,
     &                 Norm_coords(Joff), 1)
         Endif
      Enddo
C
      If (Mass_weigh_nm) Then 
          Ioff = 1
          Do Imodes = 1, Nvibs
             Do Iatoms = 1, Nreals
                Sqrtmass = Dsqrt(AtmMass(Iatoms))
                Do Ixyz = 1, 3
                   Norm_coords(Ioff) = Norm_coords(Ioff)*Sqrtmass
                   Ioff = Ioff + 1
                Enddo
             Enddo
          Enddo
      Endif
C 
      If (Get_Hess) then
C           
CSSS         Call Getrec(20,'JOBARC','CART_HES',9*Natoms*Natoms*IINTFP,
CSSS     &               VHess)
          Call Getrec(20,'JOBARC','HESSIANM',9*Nreals*Nreals*IINTFP,
     &                VHess)

#ifdef _DEBUG_LVL0
         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian from JARC"
         Call output(Vhess,1,3*nreals,1,3*nreals,3*nreals,3*nreals,1)
#endif
c
c-----Transform to the ZMAT order
C
         Call Zero(Hess, 9*Natoms*Natoms)
         Ioff = 1
         Do Ivmol = 1, Natoms
            Iatmzmat = Imap(Ivmol)
            Do Ixyz = 1, 3
               Icol = Ixyz + (Iatmzmat - 1)*3
               Do Jvmol = 1, Natoms
                  Jatmzmat = Imap(Jvmol)
                  If (Iatmzmat .NE. 0 .AND. Jatmzmat .NE. 0) Then
                      Irow = 1 + (Jatmzmat - 1)*3
                      Call Blkcpy(Vhess(ioff), 3, 1, Hess, 3*Natoms,
     &                            3*Natoms, Irow, Icol)
                      Ioff = Ioff + 3
                  Endif 
               Enddo
            Enddo
         Enddo
#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,*) "@get_ref_vibdata; after reorder"
         Call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
#endif
C
C-----Eleminate the entry for dummy atoms
C
         Ix = 0
         Jx = 1
         Do Jatm = 1, 3*Natoms
            Do Iatm = 1, 3*Natoms
               Ink = 1 + (Iatm - 0.1)/3
               Jnk = 1 + (Jatm - 0.1)/3
               If ((AtmLabel(Ink)(1:5) .NE. 'X    ') .AND.
     &            ( AtmLabel(Jnk)(1:5) .NE. 'X    ')) Then
                   Ix = Ix + 1
                   Hess(Ix, Jx) = Hess(Iatm, Jatm)
                   If (Ix .EQ. 3*Natoms) Ix = 0
                   If (Ix .EQ. 0) Jx = Jx + 1
               Endif
           Enddo
         Enddo
C
      Endif
C     
#ifdef _DEBUG_LVL0
         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian after removing dummies"
         Call output(hess,1,3*nreals,1,3*nreals,3*nreals,3*nreals,1)
#endif
C
      If (Get_Grad) Then

          Call Getrec(20,'JOBARC','GRADIENT',3*Nreals*IINTFP,
     &                A2Grad(1,1))
C
          Call Dzero(Grad, 3*Natoms)

          Do i=1, Natoms
             k=Imap(i)
             If (k .NE. 0) Then
                Do j=1,3
                   Grad(j,k)  = A2grad(j,i)
                Enddo
             Endif
          Enddo                  
C        
          Call Dcopy(3*Natoms, Grad, 1, A2grad, 1)
C
          Do i=1, Natoms
             If (Imap(i) .NE. 0) Then
                 Call Dcopy(3, A2grad(1,Imap(i)), 1, Grad(1,i), 1)
             Endif
          Enddo
C
          If (Mass_Weigh_gr) Then 
             Do Iatom = 1, Nreals
                Do Ixyz = 1, 3
                   Grad(Ixyz, Iatom) = Grad(Ixyz, Iatom)/
     &                                 Dsqrt(AtmMass(Iatom))
                Enddo
             Enddo
          Endif
C
      Endif
C
      Return
      End

