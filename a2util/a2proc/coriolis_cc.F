C
      Subroutine Coriolis_cc(Natom,Nreal,Nx,Nxm6,AtmMass,Coord,Bmat,
     &                       Gmat,Hess,vectors,vib_vectors,Btmp,Scr1,
     &                       Scr2,DM_matx,DM_maty,DM_matz,Iztov,iuhf)
C
      Implicit double precision (a-h,o-z)
C
#include "machsp.com"
#include "mxatms.par"
C
      PARAMETER(EPSILON = 1.0D-8)
      Logical FCM_EXIST
      Character*5 AtmLabel(Mxatms)
      Character*4 IrrepLabel(3*Mxatms)
C
      Dimension AtmMass(Natom),Coord(Natom*3),Bmat(Nxm6,Nx),
     &          Gmat(Nx,Nx),Gmat_inv(Nx*Nx),
     &          Hess(Nx,Nx),vectors(Nx*NX),vib_vectors(Nx*Nxm6),
     &          scr1(Nx,Nx),Btmp(Nx*Nxm6),Scr2(Nx,Nx),
     &          DM_matx(Nx,Nx),DM_maty(Nx,Nx),DM_matz(Nx,NX),
     &          Ivtoz(Natom)
C
      Double Precision M_subx(3,3),M_suby(3,3),M_subz(3,3),Mass

C     translation factor from bohr to angstrom
      Parameter (b2ang=1.889725989,wavens=5.14048D03)
C
      Call Getcrec(1,'JOBARC','ZSYM',5*Natom,AtmLabel)
      Call Getcrec(1,'JOBARC','VIB_SYMS',4*Nx,IrrepLabel)
      Call Getrec(20,'JOBARC','ATOMMASS',Natom*IINTFP,Btmp)
      Call Getrec(20,'JOBARC','COORD   ',3*Natom*IINTFP,Coord)
      Call Getrec(20,'JOBARC','BMATRIXC',NX*Nxm6*IINTFP, Bmat) 
      CALL Getrec(20,'JOBARC','CARTHESC',Nx*Nx*IINTFP,Hess)
      Call Getrec(20,'JOBARC','MAP2ZMAT',Natom,Ivtoz)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,"(a,I3)") "The Number of Real Atoms: ", Nreal
      Write(6,"(a,I3)") "The Number of Atoms:", Natom
      Write(6 ,"(1x,a5)") (AtmLabel(i),i=1,Natom)
      Write(6,*) 
      Write(6,*) "The Cartesian Coordinate in comp order"
      Write(6,*)
      ioff = 1
      Do I = 1, Natom
         Write(6, "(3F10.5)") (Coord(ioff+j-1), J=1,3)
         ioff = ioff + 3
      Enddo
      Write(6,*)
      Write(6,*) "The Bmatrix"
      CALL OUTPUT(Bmat, 1, Nxm6, 1, Nx, Nxm6, Nx, 1)
      Write(6,*)
      Write(6,"(a)") "The Hessian matrix"
      CALL OUTPUT(Hess, 1, Nx, 1, Nx, Nx, Nx, 1)
#endif
C Remove dummyies from the atomic mass array
C
      Ireal = 0
      Do Iatom = 1, Natom
         If (.Not. (Btmp(Iatom) .lt. 0.50D0)) Then
            Ireal = Ireal + 1
            AtmMass(Ireal) = Btmp(Iatom)
         Endif
      Enddo 
C
#ifdef _NO_SKIP
      Ioff = 1
      Do Iatom = 1, Natom
         Izmat = Ivtoz(Iatom)
         Joff = (Izmat-1)*3 + 1
         Call Dcopy(3, Coord(Ioff), 1, Btmp(Joff), 1)
         Ioff = Ioff + 3
      Enddo
      Call Dcopy(Nx, Btmp, 1, Coord, 1)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The Cartesian Coordinate in comp order"
      ioff = 1
      Do I = 1, Natom
         Write(6, "(3F10.5)") (Coord(ioff-1+j), j=1,3)
         ioff = ioff + 3
      Enddo
#endif
#endif
C
C Mass weigh the B-Matrix.
C
      Do Iint = 1, Nxm6
         Ixyz  = 0
         Iatom = 0
C
         Do Icart = 1, Nx/3
            Iatom = Iatom + 1
            Do Jatom = 1, 3
               Ixyz = Ixyz + 1
#ifdef _DEBUG_LVLM1
               Write(6,"(2I2,F15.6)") ixyz, Iatom, AtmMass(Iatom)
#endif
               If (AtmMass(Iatom) .eq. 0.0D0) Then 
                  Bmat(Iint, Icart) = 0.0D0
               Else
                  Bmat(Iint, Icart) = (1.0D0/DSQRT(AtmMass(Iatom)))*
     &                              Bmat(Iint, Ixyz)
               Endif
            Enddo
         Enddo
      Enddo
C
C Mass weigh the Hessian and obtain the normal modes.
C
      Do Ideg = 1, Nx
         Do Jdeg = 1, Nx
            Mass = DSQRT(AtmMass(1+(Ideg-1)/3)*AtmMass(1+(Jdeg-1)/3))
            IF (Mass .lt. 1.0D-3) Then
               Hess(Ideg, Jdeg) = 0.0D0
            Else
               Hess(Ideg, Jdeg) = Hess(Ideg, Jdeg)/Mass
            Endif
         Enddo
      Enddo 
C
#ifdef _DEBUG_LVLM1
      Write(6,*) 
      Write(6,"(a)") "The mass weighted Hessian matrix"
      CALL OUTPUT(Hess, 1, Nx, 1, Nx, Nx, Nx, 1)
#endif
      Call Eig(Hess, vectors, Junk, Nx, 0)
      CALL Getrec(20, "JOBARC", "SYM_NMDS", Nx*Nx*IINTFP, Vectors) 
      Call Dcopy(Nxm6*Nx, vectors(6*Nx+1), 1, vib_vectors, 1)
C
      Do I = 1, NX
         IF (Hess(I,I) .Lt. 0.0D0 .and. (Dabs(Hess(I,I)) .gt.
     &       1.0D0)) Then
            Write(6,"(a,a)")" Imaginary eigenvlaues in the Hessian,",
     &                       " modifications to the code" 
            Write(6,"(a)")" need to be made to handle them"
            Call Errex
         Endif
      Enddo

#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,"(a)") "The eig vectors of Hessian"
      CALL OUTPUT(vectors, 1, Nx, 1, Nx, Nx, Nx, 1)
      Write(6,*)
#endif
C
      Call print_coord(Coord, Atmlabel, B2ang, Natom)
      Call Header("The Normal modes", -1, 6)
      CALL OUTPUT(vib_vectors, 1, Nx, 1, Nxm6, Nx, Nxm6, 1)
      Write(6,*)
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "Check for the sign convention using CH3+"
      Write(6,*) "Reverse the natural sign of mode 5" 
      Call Dscal(NX, -1.0D0, vib_vectors(4*Nx +1), 1)
      Call Header("The Normal modes with revese sign", -1, 6)
      CALL OUTPUT(vib_vectors, 1, Nx, 1, Nxm6, Nx, Nxm6, 1)
      Write(6,*)
      Write(6,*) "Check for the orthogonolity of modes"
      Call Xgemm("T", "N", NXM6, NXM6, NX, 1.0D0, vib_vectors, NX,
     &            Vib_vectors, NX, 0.0D0, Vectors, NX)
      CALL OUTPUT(vectors, 1, Nxm6, 1, Nxm6, Nx, Nx, 1)
#endif

      Call Header("The vibrational frequencies in cm^-1", -1, 6)
C
      Nleft = Mod(Nx, 5)
      istart = 1
      Iend   = 5
      Do Ivib = 1, Int(Nx/5)
         Write(6,"(11x,5(a,10x))") (IrrepLabel(I), I=Istart, Iend)
         Write(6,*)
         Write(6,"(5F15.5)") (wavens*DSQRT(DABS(Hess(i,i))),
     &         i=Istart, Iend)
         Write(6,*)
          istart = istart + 5
          Iend   = Iend   + 5
      Enddo
      If (Nleft .Ne. 0) Then
         Istart = Nx - Nleft + 1
         Iend   = Nx
         Write(6,"(11x,5(a,10x))") (IrrepLabel(I),I=Istart, Iend)
         Write(6,*)
         Write(6,"(5F15.5)") (wavens*DSQRT(DABS(Hess(i,i))),
     &         i=istart, iend)
      Endif
C
C Form the G^-1.
C
      Call Dcopy(Nx*NXm6, Bmat, 1, Btmp, 1)
C
      CALL XGEMM('N', 'T', Nxm6, Nxm6, Nx, 1.0D0, Bmat, Nxm6,
     &           Btmp, Nxm6, 0.0D0, Gmat, Nxm6)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,"(a)") "The G-Matrix:BB^t"
      CALL OUTPUT(Gmat, 1, Nxm6, 1, Nxm6, Nxm6, Nxm6, 1)
#endif

C
C The intermediate G matrix created above is linear dependent, hence
C there are eigenvalues that are zero. Invert the non-zero digonal
C elements to built the Lambda(-1) matrix.
C
      CALL EIG(Gmat, scr1, 1, Nxm6, 1)
C
      NULLEVAL = 0
      DO I = 1, Nxm6
         IF (Gmat(I, I) .LE. EPSILON) THEN
             NULLEVAL = NULLEVAL + 1
             Gmat(I, I) = 0.0D0
         ELSE
             Gmat(I, I) = 1.0D0/Gmat(I, I)
         ENDIF
      ENDDO
C
C Built the generalized inverse of G-matrix, what proceeded is generally
C known as singular value decomposition (ineversion of singular matrices)
C G^(-1) = [K L] Lambda^(-1)[K^(t) L^(t)]

      CALL XGEMM('N', 'N', Nxm6, Nxm6,  Nxm6, 1.0D0, Scr1, Nxm6,
     &           Gmat, Nxm6, 0.0D0, Scr2, Nxm6)
C
      CALL XGEMM('N', 'T', Nxm6, Nxm6, Nxm6, 1.0D0, Scr2, Nxm6,
     &            Scr1, Nxm6, 0.0D0, Gmat, Nxm6)
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,"(a)") "The G-1 matrix"
      CALL OUTPUT(Gmat, 1, Nxm6, 1, Nxm6, Nxm6, Nxm6, 1)
#endif
C
C Construct the M_sub and M_Mat for x, y and z directions.
C
      Call zero(M_Subx, 9)
      Call zero(M_suby, 9)
      Call zero(M_Subz, 9)
      
      M_subx(2,3) =  1.0D0
      M_Subx(3,2) = -1.0D0
      M_suby(1,3) = -1.0D0
      M_suby(3,1) =  1.0D0
      M_subz(1,2) =  1.0D0
      M_subz(2,1) = -1.0D0
C
      Call zero(DM_matx, Nx*NX)
      Call zero(DM_maty, Nx*Nx)
      Call zero(DM_matz, Nx*Nx)
   
      Do Iatom = 1, Nreal
         Do Jatom = 1, Nreal
            If (Iatom .EQ. Jatom) Then
                Ioff =  (Iatom-1)*3
                Joff =  (Jatom-1)*3
C
                Do I = 1, 3
                   Ioff = Ioff + 1
                   Do J = 1, 3
                      Joff = Joff + 1
C
                      DM_matx(Ioff, Joff) = M_subx(I, J)
                      DM_maty(Ioff, Joff) = M_suby(I, J)
                      DM_matz(Ioff, Joff) = M_subz(I, J)
C
                 Enddo
                 Joff = (Jatom-1)*3
               Enddo
           Endif
        Enddo
      Enddo
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "The MX Matrix"
      Call output(DM_mATX, 1, Nx, 1, Nx, Nx, Nx, 1)
      Write(6,*)
      Write(6,*) "The MY Matrx"
      Call output(DM_maty, 1, Nx, 1, Nx, Nx, Nx, 1)
      Write(6,*)
      Write(6,*) "The Mz Matrix"
      Call output(DM_matz, 1, Nx, 1, Nx, Nx, Nx, 1)
#endif
C X-Rotation
      call zero(scr2, Nx*Nx)
      CALL XGEMM('N', 'N', Nx, Nxm6, Nx, 1.0D0, DM_matx, Nx,
     &           Vib_Vectors, Nx, 0.0D0, Scr1, Nx)
C
      CALL XGEMM('T', 'N', Nxm6, Nxm6, Nx, 1.0D0, Vib_Vectors,
     &            Nx, Scr1, Nx, 0.0D0, Scr2, Nxm6)
C
      Call Header("The Coriolis Coupling Coefficients, Rotation on X",
     &             -1, 6)
      Call output(Scr2, 1, Nxm6, 1, Nxm6, Nxm6, Nxm6, 1)
      Write(6,*)
C Y-Rotation
      call zero(scr2, Nx*Nx)
      CALL XGEMM('N', 'N', Nx, Nxm6, Nx, 1.0D0, DM_maty, Nx,
     &           Vib_Vectors, Nx, 0.0D0, Scr1, Nx)
C
      CALL XGEMM('T', 'N', Nxm6, Nxm6, Nx, 1.0D0, Vib_Vectors,
     &            Nx, Scr1, Nx, 0.0D0, Scr2, Nxm6)
C
      Call Header("The Coriolis Coupling Coefficients, Rotation on Y", 
     &             -1, 6)
      Call output(Scr2, 1, Nxm6, 1, Nxm6, Nxm6, Nxm6, 1)
C Z-Rotation
      call zero(scr2, Nx*Nx)
      CALL XGEMM('N', 'N', Nx, Nxm6, Nx, 1.0D0, DM_matz, Nx,
     &           Vib_Vectors, Nx, 0.0D0, Scr1, Nx)
C
      CALL XGEMM('T', 'N', Nxm6, Nxm6, Nx, 1.0D0, Vib_Vectors,
     &            Nx, Scr1, Nx, 0.0D0, Scr2, Nxm6)
C
      Call Header("The Coriolis Coupling Coefficients, Rotation on Z", 
     &             -1, 6)
      Call output(Scr2, 1, Nxm6, 1, Nxm6, Nxm6, Nxm6, 1)

      Write(6,*) 
      Write(6,"(T2,a,a)") "Coriolis Coupling constants were",
     &                    " successfully computed" 
      Return
      End
