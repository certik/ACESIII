#include <flags.h>
      Subroutine Pes_scan_main(Streaming, Stationary, Drive_IRC)
      Implicit Double Precision (A-H, O-Z) 
C
#include <icore.com>
#include <machsp.com>
#include <mxatms.par>
#include <jodaflags.com>
#
C
      Double Precision Norm_coords, Omega, Lamsearch_Tol, 
     &                 Ln_Intrp_Tol
      Character*80 Fname
      Character*80 tmp
      Character*5 AtmLabel(Mxatms), Constr_Min
      Character*4 SymLabel(3*Mxatms), Units
      Character*12 Grid_file
      Character*11 Vib_Type(3*Mxatms)

C
C     Conversion factors from bohr to angstrom and atomic units to cm-1
C
      Parameter (B2ang=1.889725989, Au2Invcm=5.14048D03)
C
      Logical ZMT_Present, GBS_Present, JARC_Present, JNDX_Present,
     &        Can_Do_Aces, Ref_done, User_input, Direct, 
     &        Grid_file_Present, Drive_IRC, Trans_state, 
     &        Get_Hess, More_pts_left, End_point, First_step,
     &        Forward, Get_grad, Mass_Weigh_nm, Mass_Weigh_gr,
     &        New_IRC

      Logical FlgACESGeom,FlgACESElec,FlgACESGrad,FlgACESHess,
     &        FRQARC_EXIST, Streaming, Stationary, Wrt_extrnls
C 
      Dimension VCoords(Mxatms*3), A2grad(Mxatms*3), Rgrad(Mxatms*3),
     &          A2hess(Mxatms*Mxatms*9), Imap(Mxatms), IOrd(Mxatms),
     &          AtmMass(Mxatms), Rhess(Mxatms*Mxatms*9),
     &          Norm_coords(9*Mxatms*Mxatms), Omega(3*Mxatms),
     &          Vomega(3*Mxatms),Btmp(3*Mxatms), Deltaq(3*Mxatms),
     &          Deltax(3*Mxatms), Coords(Mxatms*3), 
     &          Coords_K0O(Mxatms*3), Coords_K1C(Mxatms*3),  
     &          Coords_K1P(Mxatms*3), Coords_K1CM(Mxatms*3),
     &          Work(36*Mxatms*Mxatms),Grad_K0O(3*Mxatms),
     &          Vec_K1C0(Mxatms*3), Grad_on_K0O(Mxatms*3), 
     &          Grad_on_K1C(Mxatms*3), Vec_K0O(Mxatms*3), 
     &          Steep_Dpath(3*Mxatms), Grad_K1CM(3*Mxatms),
     &          Vec_K1C1(Mxatms*3)
c
c First check whether JOBARC/JAINDX present in the current directory.
c
      Call Gfname('ZMAT',   Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=ZMT_Present)
      Call Gfname('GENBAS', Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=GBS_Present)
  
      Call Gfname('JOBARC', Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=JARC_Present) 
      Call Gfname('JAINDX', Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=JNDX_Present)
      If (ZMT_Present .And. GBS_Present)   Can_Do_Aces = .True.
      If (JARC_Present .And. JNDX_Present) Ref_done    = .True. 
c 
      FlgACESGeom = .TRUE.
      Wrt_Extrnls = .TRUE.
      FlgACESElec = .FALSE.
      FlgACESGrad = .FALSE.
      FlgACESHess = .FALSE.
C

      If (Can_Do_Aces .And. .Not. Ref_done) Then
c 
c There is a reference ZMAT and GENBAS (minimum requiremnt for 
c an ACES calculations) Lets run the reference point. For PES scans
c on a user defined grid, the reference run is a single point 
c energy or gradient calculation. 
c
         Call Runit("runaces2a")
         call aces_init_rte
         call aces_com_parallel_aces
         call aces_ja_init
C
      Else if (Ref_done) Then
C
         Call aces_init_rte
         Call aces_com_parallel_aces
         Call aces_ja_init
C
      Endif 
C
      call getrec(1,'JOBARC','IFLAGS',   100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',  500,iflags2)
      call getrec(1,'JOBARC','NREALATM', 1,  Nreals)
      Call Getrec(1,'JOBARC','ZMATATMS', 1,  Natoms)
      Call Getrec(1,'JOBARC','ATOMMASS', Natoms*IINTFP,AtmMass)
      Call Getrec(1,'JOBARC','UHFRHF  ', 1,iUHF)
      Call Getrec(1,'JOBARC','NIRREP  ', 1,NIrRep)
      Call Getrec(1,'JOBARC','MAP2ZMAT', Natoms,Imap)
      Call Getrec(20, 'JOBARC', 'LINEAR  ', 1, Iamlinear)
      Call Getcrec(1,'JOBARC','ZSYM',5*Natoms,AtmLabel)
C
      If (Iamlinear .EQ. 1) Then
         Nvibs = 3*Nreals - 5
      Else
         Nvibs = 3*Nreals - 6 
      Endif
c
CSSS      Call a2_reset_jarc()
C
      If (streaming) Then
C
C This is for direct dynamics runs. The reference point has alredy
C been run and the geometry of the point that need to be done is 
C pass from the dynamics program. 
c
         Call a2_reset_jarc()
C
#ifdef _DEBUG_LVL0
         Call Getrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP, 
     &               Coords)
         Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
C
         Ioff = 1
         Do Iatm = 1, 2
            Do Ixyz = 1, 3
               Coords(ioff) = (Coords(ioff)) + 0.01
               ioff = ioff + 1
            Enddo
         Enddo
         Write(6,*)
         Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
         Wrt_Extrnl = .True.
#endif
C
C Write the new incomming geometry. The next calculation is performed
C at this geometry.
C
         Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP,
     &               Coords)
         Call aces_ja_fin
C
C Run the point and read the energy and the gradients and Hessians
C if available. 
C
         Call Runit("runaces2b")
         Call aces_ja_init
C
#ifdef _DYNAMIC_MEMORY
         icrsiz = iflags(h_IFLAGS_mem)
         icore(1) = 0
         do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
            call aces_malloc(icrsiz,icore,i0)
            if (icore(1).eq.0) icrsiz = icrsiz - 1000000
         end do
         if (icore(1).eq.0) then
             print *, '@MAIN: unable to allocate at least ',
     &                 1000000,' integers of memory'
            call aces_exit(1)
         end if
         iatmMass = iocc + 16
         icoord   = iatmMass + 3*Natoms*iintfp
         ivgrad   = icoord   + 3*NAtoms*iintfp
         izgrad   = ivgrad   + 3*Nreals*iintfp
         ihess    = izgrad   + 3*Nreals*iintfp
         ifreq    = ihess    + 9*Natoms*Natoms*iintfp
         inex     = ifreq    + 3*Nreals*iintfp
#endif
C
         Call Get_flags(FlgACESElec, FlgACESGrad,  FlgACESHess,
     &                  FlgACESForc)
C
         Call A2read_Jarc(FlgACESGeom, FlgACESElec, FlgACESGrad,
     &                    FlgACESHess, FlgACESForc, Natoms, Nreals, 
     &                    Imap, Coords, A2Ener, AtmMass, AtmLabel, 
     &                    Rgrad, A2grad, Rhess, A2hess, Vomega, 
     &                    Omega, Iuhf, Nirrep, Nvibs, B2ang, 
     &                    Au2Invcm, Wrt_extrnls)
C
         Call a2_reset_jarc()
C
      Else if (Stationary) Then 
C
C This is primarily reserved for running single points calculations 
C on a grid. The grid can be given in a file or can be genrerated.
C The GIRD file (if not automatic) must be formatted file of Cartesian 
C coordinates. 
C
         Call A2get_gridnmlist(Grid_file, Nmbrof_gridpts)

C 1. GIRD can be user defined. If user defined, give the number
C    numbers of atoms (as in the reference ZMAT) whose coordinates
C    are changing (others are kept at the reference values).
C
C 2. Automatically generated on normal modes. The refence point
C    must be a frequency calculations. 
C
         If (Grid_file .NE. "NORMAL_MODES") Then

            Inquire(FILE=Grid_file,EXIST=Grid_file_Present)

            If (Grid_file_Present) Then
                Open(Unit=13, File=Grid_file, Form="FORMATTED",
     &               Status="OLD")
            Else
                Write(6, "(a,1x,a,a)") "The file", Grid_file, 
     &                                     "does not exsist." 
                Call Errex
            Endif 
C
            Call a2_reset_jarc()
            Do Ipts = 1, Nmbrof_gridpts 
C
C Get the point and write the JOBARC record
C
               Read(13, *, End=20)
               Do Iatm = 1, Natoms
                  Ioff = 3*(Iatm - 1)
                  Read(13,*,End=20) (Coords(Ioff + I), I=1, 3)
               Enddo
C       
#ifdef _DEBUG_LVL0
                Write(6,*)
                Write(6, "(a)") "Grid points read from Grid_file"
                Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
#endif
    
               Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP, 
     &                    Coords)
               Call aces_ja_fin
C
C Run the point and read the energy and the gradients 
C
               Call Runit("runaces2b")
               Call aces_ja_init
C
               Call Get_flags(FlgACESElec, FlgACESGrad,  FlgACESHess,
     &                        FlgACESForc)
C
C Extract the energy, gradients and Hessians (if they are in JOBARC).
C
               Call A2read_Jarc(FlgACESGeom, FlgACESElec, 
     &                          FlgACESGrad, FlgACESHess, 
     &                          FlgACESForc, Natoms, 
     &                          Nreals, Imap, Coords, A2Ener, 
     &                          AtmMass, AtmLabel, Rgrad, A2grad, 
     &                          Rhess, A2hess, Vomega, Omega, Iuhf, 
     &                          Nirrep, Nvibs, B2ang, Au2Invcm, 
     &                          Wrt_extrnls)
               Call a2_reset_jarc()
            Enddo
  20        Close(13)
C
         Else 
C
C This block handles scan along normal modes. The requirement is that
C refrence job is a vibrational frequency calculation. The displacement
C size and the length must be provided in PES_SCAN namelist.
C 
            Trans_state   = .False.
            Mass_weigh_nm = .True. 
            Mass_weigh_gr = .False. 
            Get_Hess      = .False.
            Get_Grad      = .False.
            Call Get_refvib_data(Coords, Norm_coords, Vomega, Omega,
     &                           AtmMass, AtmLabel, SymLabel, Coords,  
     &                           Btmp, A2grad, Rgrad, A2Hess, Rhess, 
     &                           Imap, Vib_Type, Nreals, Natoms, 
     &                           Nvibs, B2ang, Au2invcm, Trans_state, 
     &                           Mass_Weigh_nm, Mass_Weigh_gr, 
     &                           Get_hess, Get_Grad)

#ifdef _OBSOLETE
            Call Getrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP, 
     &                  Coords)
            Call Getrec(20, "JOBARC", "NORMCORD", 3*Nreals*3*Nreals*
     &                  IINTFP, Norm_coords)
            Call Getrec(20, "JOBARC", 'FORCECON', 3*Nreals*IINTFP, 
     &                  VOmega)
            Ivib = 0 
            Do Imode = 1, 3*Nreals
C        
               If (Vomega(Imode) .lt. 0.0D0) Then
                   Vomega(Imode) = - DSQRT(DABS(Vomega(Imode)))*
     &                             Au2Invcm
               Else
                   Vomega(Imode) = DSQRT(Vomega(Imode))*Au2Invcm
               Endif
C        
               If (DABS(Vomega(Imode)) .GT. 1.0D0) Then
                  Ivib = Ivib + 1
                  Omega(Ivib) = Vomega(Imode)
               Endif
C     
            Enddo
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
                  Call Dcopy(3, Coords(Icord), 1, Deltax(Jcord), 1)
             Endif
             Icord = 1
             Jcord = 1
            Enddo
C
            Do Iatom = 1, 3*Nreals
               SymLabel(Iatom) = SymLabel(6+Iatom)
            Enddo
            Ioff = (3*Nreals - Nvibs)*3*Nreals
            Do Imodes = 1, NVibs*3*Nreals
               Norm_coords(Imodes) = Norm_coords(Ioff + Imodes)
            Enddo 
C
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
#endif
C
            Call a2_reset_jarc()
C 
#ifdef _DEBUG_LVL0
            Write(6,*) 
            Write(6,*) "The vibrational frequencies and symmetries"
            Write(6, "(4F12.5,1x)") (Omega(I), I = 1, Nvibs)
            Write(6,*)
            Write(6, "(4a4)") (SymLabel(I), I = 1, Nvibs)
            Write(6,*)
            Write(6,*) "The reference geometry"
            Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
            Write(6,*) "Normal modes"
            Do I = 1, Nvibs
               Ioff = (I-1)*3*Nvibs
                  Write(6,"(3F10.5)") (Norm_coords(Ioff + k), k=1, 
     &                                 3*Nreals)
                  Write(6,*)
            Enddo
#endif
C
C The vibrational SCF, averaging etc and many other post electronic
C structure calcualtions depends on having PES on normal modes.
C
C V(q_1,q_2, ...,q_n) = E_0(q_1(0), q_2(0)..) + Sum_i(DE/dq_i)q_i + C                    Sum_ij(DE^2/dq_idq_j)q_iq_j + Sum_ijk(D^3E/dq_idq_jdq_k) C                    q_iq_jq_k + Sum_ijkl(D^4E/dq_idq_jdq_kdq_l)q_iq_jq_kq_l
C
C The curveture dependent dimensionless step size (JCP, 121, 1383. 2004).
C 
           
            Delx = 0.50D0
            Call Stpsize(Omega, Deltaq, Delx, Nvibs)
C
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "Dimensionless step size will be"
         Write(6, "(4F15.5)") (Deltaq(I), I=1, Nvibs)
#endif
C
            Call Dscal(Nvibs, 0.0D0, Deltaq, 1)
            Call Dcopy(3*Nreals, Coords, 1, Deltax, 1)
            Call Run_ref(Deltaq, Deltax, Norm_coords, Nreals, 
     &                   Nvibs) 
C
            Call Get_flags(FlgACESElec, FlgACESGrad,  FlgACESHess,
     &                    FlgACESForc)

            Call A2read_Jarc(FlgACESGeom, FlgACESElec,
     &                       FlgACESGrad, FlgACESHess,
     &                       FlgACESForc, Natoms,
     &                       Nreals, Imap, Coords, A2Ener,
     &                       AtmMass, AtmLabel, Rgrad, A2grad,
     &                       Rhess, A2hess, Vomega, Omega, Iuhf, 
     &                       Nirrep, Nvibs, B2ang, Au2Invcm, 
     &                       Wrt_extrnls)
            Call a2_reset_jarc()

C
CSSS               Call Run_1st(Deltaq, Norm_coords, Nvibs)
       
CSSS               Call Run_2nd(Deltaq, Norm_coords, Nvibs)

         Endif
C
      Else If (Drive_IRC) Then
C
C If the starting point is the transition state, the fist pivot point 
C is half distance from the steepest decent path. The reference calculation
C is a vibrational frequnecy calculation of the transition state.
C
C **** Nature of the starting point must be read from the IRC namelist***
C
           Trans_state   = .True. 
           Mass_weigh_nm = .False. 
           Mass_weigh_gr = .False. 
           Get_Hess      = .True.
           Get_Grad      = .True. 
           Call Get_refvib_data(VCoords, Norm_coords, Vomega, Omega,
     &                          AtmMass, AtmLabel, SymLabel, Coords, 
     &                          Btmp, A2grad, Grad_K1CM, A2hess, Rhess, 
     &                          Imap, Vib_Type, Nreals, Natoms, Nvibs,  
     &                          B2ang, Au2Invcm, Trans_state, 
     &                          Mass_weigh_nm, Mass_weigh_gr, 
     &                          Get_hess, Get_grad)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The Hessian @-PES_scan"
      Call output(Rhess, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals,
     &            3*Nreals, 1)
#endif
#ifdef _DEBUG_LVL0 
            Write(6,*)             
            Write(6,*) "The vibrational frequencies and symmetries"
            Write(6, "(4F12.5,1x)") (Omega(I), I = 1, Nvibs)
            Write(6,*)
            Write(6, "(4a4)") (SymLabel(I), I = 1, Nvibs)
            Write(6,*)
            Write(6,*) "The reference geometry"
            Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
            Write(6,*) "Normal modes"
            Do I = 1, Nvibs
               Ioff = (I-1)*3*Nvibs
                  Write(6,"(3F12.8)") (Norm_coords(Ioff + k), k=1, 
     &                                 3*Nreals)
                  Write(6,*)
            Enddo 
#endif      
C           
           Call Getrec(20,'JOBARC','TOTENERG', IINTFP, E_refrence)
           Call a2_reset_jarc()
C 
C The eigenvector corresponding to the negative eigenvalue (no longer
C mass weighted). Generate the pivot point: that is 1/2stride of 
C the steepest decent eigenvector (when the hessian is availabe). 
C Otherwise choose the largest gradient 
C
           N_Irc_step    = 0
           More_pts_left = .True.
           End_point     = .False.
           First_Step    = .True.
           Forward       = .False.
           Opt_tol       = Iflags2(h_IFLAGS_convergence)
           Max_Cycles    = Iflags(h_IFLAGS_opt_maxcyc)
           Stride        = 0.30D0
           Stride_total  = 0.0D0
           Stride_Tol    = 1.0D-04
           Ln_Intrp_Tol  = 5.0D0
           Hess_Eval_Tol = 1.0D-09
           Delta         = 5.0D-02
           lamsearch_TOl = 1.0D-16
           Binsearch_Tol = 1.0D-08
           Small         = 1.0D-05
           G_Cutoff      = 1.0D-04
           R_Cutoff      = 1.0D-03
           Units         = "Bohr"
           Halfs         = 0.50D0*Stride
C
           Do while (More_pts_left .AND. .NOT. End_point) 
C
              If (First_step) Then
C
                  If (Trans_state) Then 
                     Call Dcopy(3*Nreals, Norm_Coords, 1, 
     &                         Steep_DPath, 1)
C
                     Dpath_Norm = 0.0D0
                     Ioff       = 1
                     Do Iatom = 1, Nreals
                        Sqrtmass = Dsqrt(AtmMass(Iatom))
                        Do Ixyz = 1, 3
                           Dlen       = Steep_DPath(Ioff)*Sqrtmass
                           Dpath_Norm = Dlen*Dlen + Dpath_Norm 
                           Ioff = Ioff + 1 
                        Enddo
                     Enddo
#ifdef _DEBUG_LVL0   
            Write(6,*) 
            Write(6,*) "@-PES_scan The steepest decent path"
            Write(6, "(3F12.8)") (Steep_DPath(i),i=1,3*Nreals)
#endif
                     Dpath_Norm = DSqrt(Dpath_Norm)
C
                     If (Dpath_Norm .LT. Small) Then
                         Write(6,"(a,a,a)") "Norm of the ",
     &                                   " mass-weighted steepest",
     &                                   " decent path is too small."
                          Call Errex
                     Endif
C
                     Call Dscal(3*Nreals, 1.0D0/Dpath_Norm, 
     &                          Steep_DPath, 1)
C
                  Else
                      Grad_norm = Ddot(3*Nreals, Grad_K1C, 1, 
     &                                 Grad_K1C, 1)
                      Grad_norm = Dsqrt(Grad_norm)
                      If (Grad_norm .LT. Small/1.0D03) Then
                         Write(6,"(a,a,a,a)") "Gradient norm is too"
     &                                   " small, can not start from"
     &                                   " stationary point except"
     &                                   " the transition state."
                         Call Errex
                      Endif 
                      Call Dscal(3*Nreals, -Halfs/Grad_norm, 
     &                           Grad_K1C, 1)
                      
                  Endif 
C                  
                  Step = 0.5D0*Stride
C
                  Imax = Idamax(3*Nreals, Steep_Dpath, 1)
C
                  If (Forward .AND. Steep_Dpath(Imax) .LT. 0.0D0) 
     &               Step = -1.0D0*Stride
                  If (.NOT. Forward .AND. Steep_Dpath(Imax) .GT.
     &               0.0D0)
     &               Step = -1.0D0*Stride

C
C _K0O is the pevious point, _K1P is the pivot point and _K1C is the 
C current point.
C
                  Call Dcopy(3*Nreals,        Coords,    1, 
     &                       Coords_K0O, 1)
                  Call Dcopy(3*Nreals,        Coords,    1, 
     &                       Coords_K1P, 1)
                  Call Dcopy(3*Nreals,        Coords,    1, 
     &                       Coords_K1C, 1)
C
                  If (Trans_state) Then
                      Call Daxpy(3*Nreals, Step,   Steep_DPath, 1, 
     &                           Coords_K1P, 1)
                      Call Daxpy(3*Nreals, 2*step, Steep_DPath, 1, 
     &                           Coords_K1C, 1)
                      Call Dcopy(3*Nreals,        Coords_K1C,    1, 
     &                          Coords_K1CM, 1)
                  Else
                      Call Daxpy(3*Nreals, 1.0D0, Grad_K1C, 1, 
     &                           Coords_K1P, 1)
                      Call Daxpy(3*Nreals, 2.0D0, Grad_K1C, 1, 
     &                           Coords_K1C, 1)
                      Call Dcopy(3*Nreals,        Coords_K1C,    1, 
     &                           Coords_K1CM, 1)
                  Endif

#ifdef _DEBUG_LVL0
            Write(6,*) "The current geometry"
            Write(6, "(3F12.8)") (Coords_K1C(i),i=1,3*Natoms)
            Write(6,*) "The pivot point geometry"
            Write(6, "(3F12.8)") (Coords_K1P(i),i=1,3*Natoms)
            Write(6,*) "The previous geometry"
            Write(6, "(3F12.8)") (Coords_K0O(i),i=1,3*Natoms)
#endif  
C
                  Ncycles = 1
                  Constr_Min = "Null "

                  Do While (.NOT. (Constr_Min .EQ. "Found") .AND.
     &                             Ncycles .Lt. Max_Cycles)
C
                           Call Constr_search(Coords_K0O, Coords_K1C, 
     &                                        Coords_K1P, Vcoords,
     &                                        Coords, Coords_K1CM, 
     &                                        A2grad, Grad_K1CM,  
     &                                        Grad_K0O, A2Hess,
     &                                        Rhess, AtmMAss, 
     &                                        Vec_K1C0, Vec_K1C1,
     &                                        Vec_K0O, Grad_on_K0O, 
     &                                        Grad_on_K1C,
     &                                        Work, B2ang, Units,
     &                                        Stride, Delta, 
     &                                        First_step, New_IRC,
     &                                        Constr_Min,
     &                                        Imap, Ncycles, Nreals, 
     &                                        Natoms, Opt_Tol,
     &                                        Stride_Total, 
     &                                        Stride_Tol,  
     &                                        Ln_Intrp_Tol, 
     &                                        Hess_Eval_Tol, 
     &                                        Lamsearch_Tol,
     &                                        Binsearch_Tol,
     &                                        G_cutoff, R_Cutoff,
     &                                        E_currnt) 

#ifdef _DEBUG_LVL0
            Write(6,*) Ncycles, Constr_Min
            Write(6,*) "@-pes_scan, updated geoemtries"
            Write(6,*) "The previous geometry"
            Write(6, "(3F12.8)") (Coords_K0O(i),i=1,3*Nreals)
            Write(6,*) "The current geometry"
            Write(6, "(3F12.8)") (Coords_K1CM(i),i=1,3*Nreals)
            Write(6,*) "The pivot point geometry"
            Write(6, "(3F12.8)") (Coords_K1P(i),i=1,3*Nreals)
            Write(6,*) "The M.W. current gradient "
            Write(6, "(3F12.8)") (Grad_K1CM(i),i=1,3*Nreals)
            Write(6,*) "The N.M.W. current gradient "
            Write(6, "(3F12.8)") (A2grad(i),i=1,3*Nreals)
            Write(6, "(2F12.8)") E_refrence, E_currnt

#endif
                           First_step = .False.
                           New_IRC    = .False.
                           Ncycles    = Ncycles + 1
                  Enddo 
C
                  If (E_currnt .GT. E_refrence) Then
                      Write(6, "(a,a)") "Energy can not rise during a"
     &                              " downhill IRC search."
                      Call Errex
                  Else
                      E_refrence = E_currnt
                  Endif  

                  First_step = .False.
                  New_IRC    = .True.

              Else
C 
C Generate the new pivot point using the gradient. Note that only the
C first step we have an exact Hessian. 
C 
                  Call Dcopy(3*Nreals, Coords_K1C, 1, Coords_K0O, 1)

                  Grad_length = Ddot(3*Nreals, Grad_K1CM, 1, 
     &                               Grad_K1CM, 1)
                  Call Dcopy(3*Nreals, Grad_K1CM, 1, Work, 1)
                  Call Dscal(3*Nreals, 1.0D0/Dsqrt(Grad_length), 
     &                       Work, 1)
#ifdef _DEBUG_LVL0
           Write(6,*) "The normalized M.W. current gradient"
           Write(6, "(3F12.8)") (Work(i),i=1,3*Nreals)
#endif
                  Ioff = 1
                  Do Iatom = 1, Nreals
                     Sqrtmass = Dsqrt(AtmMass(Iatom)) 
                     Do Ixyz = 1, 3
                        Work(Ioff) = Work(Ioff)/Sqrtmass
                        Ioff = Ioff + 1
                     Enddo
                  Enddo
C
                  Call Dscal(3*Nreals, -Halfs, Work, 1)
C
#ifdef _DEBUG_LVL0
           Write(6,*) "The gradient update"
           Write(6, "(3F12.8)") (Work(i),i=1,3*Nreals)
#endif
C
                  Call Dcopy(3*Nreals, Coords_K1C, 1, Coords_K1P, 1)
                  Call Daxpy(3*Nreals, 1.0D0, Work, 1, Coords_K1P, 1)
                  Call Dcopy(3*Nreals, Coords_K1P, 1, Coords_K1C, 1)
                  Call Daxpy(3*Nreals, 1.0D0, Work, 1, Coords_K1C, 1)
                  Call Dcopy(3*Nreals, Coords_K1C, 1, Coords_K1CM, 1)
C
#ifdef _DEBUG_LVL0
            Write(6,*) "The current geometry"
            Write(6, "(3F12.8)") (Coords_K1C(i),i=1,3*Nreals)
            Write(6,*) "The pivot point geometry"
            Write(6, "(3F12.8)") (Coords_K1P(i),i=1,3*Nreals)
            Write(6,*) "previous geometry"
            Write(6, "(3F12.8)") (Coords_K0O(i),i=1,3*Nreals)
#endif
C
                  Ncycles = 1
                  Constr_Min = "Null "

                  Do While (.NOT. (Constr_Min .EQ. "Found") .AND.
     &                             Ncycles .Lt. Max_Cycles)
C
                           Call Constr_search(Coords_K0O, Coords_K1C,
     &                                        Coords_K1P, Vcoords,
     &                                        Coords, Coords_K1CM,
     &                                        A2grad, Grad_K1CM,
     &                                        Grad_K0O, A2Hess,
     &                                        Rhess, AtmMAss,
     &                                        Vec_K1C0, Vec_K1C1, 
     &                                        Vec_K0O, Grad_on_K0O,
     &                                        Grad_on_K1C,
     &                                        Work, B2ang, Units,  
     &                                        Stride, Delta, 
     &                                        First_step, New_IRC,
     &                                        Constr_Min, 
     &                                        Imap, Ncycles, Nreals,
     &                                        Natoms, Opt_Tol,
     &                                        Stride_Total,
     &                                        Stride_Tol,
     &                                        Ln_Intrp_Tol,
     &                                        Hess_Eval_Tol,
     &                                        Lamsearch_Tol,
     &                                        Binsearch_Tol,
     &                                        G_Cutoff, R_Cutoff,
     &                                        E_currnt)
C
                           First_step = .False.
                           New_IRC    = .True.
                           Ncycles    = Ncycles + 1
                  Enddo
C
                  If (E_currnt .GT. E_refrence) Then
                      Write(6, "(a,a)") "Energy can not rise during a"
     &                              " downhill IRC search."
                      Call Errex
                  Else
                      E_refrence = E_currnt
                  Endif
                  STOP
C       
              Endif
C
           Enddo                     
              
      Endif
C
      Return
      End

