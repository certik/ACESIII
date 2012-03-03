      SUBROUTINE BUILT_BGMTRX(CARTCOORD, REDUNCO, IREDUNCO,
     &                        TOTREDNCO, TOTNOFBND, TOTNOFANG,
     &                        TOTNOFDIH, NRATMS, 
     &                        BMATRX, GMATRX, EIGVECTORS,
     &                        EPSILON, BTGMIN, DERBMAT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "io_units.par"
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "jodaflags.com"
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
      LOGICAL TS_SEARCH, MN_SEARCH_EXACT_HESS, CMP_DBMAT
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C 
C BTMP : used as temp storage of size 9*MXATMS*MXATMS need
C        to be dynamically allocated and discarded.
      DIMENSION CARTCOORD(3*NRATMS), IREDUNCO(4, MAXREDUNCO),
     &          BMATRX(TOTREDNCO, 3*NRATMS), 
     &          BTMP(3*MXATMS, 3*MXATMS), 
     &          GMATRX(TOTREDNCO, 3*NRATMS),
     &          EIGVECTORS(TOTREDNCO, TOTREDNCO),
     &          REDUNCO(TOTREDNCO)
      DIMENSION DERBMAT(3*NRATMS,3*NRATMS*TOTREDNCO),
     &          BTGMIN(3*NRATMS,TOTREDNCO)
         
      DATA MONE /-1.0/
C  
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
      CALL ZERO(BMATRX, TOTREDNCO*(3*NRATMS))
      CALL ZERO(DERBMAT,9*NRATMS*NRATMS*TOTREDNCO)
C
      CMP_DBMAT = (IFLAGS2(h_IFLAGS2_init_hessian) .GT. 0) .OR.
     &            (IFLAGS2(h_IFLAGS_eval_hess)    .gt. 0)
C
      DO 20 IBNDS = 1, TOTNOFBND
C
         IF (IREDUNCO(2, IBNDS) .NE. 0) THEN  
C
            CALL  BULT_BNDCRD(CARTCOORD, BMATRX, DISTAB,
     &                        IREDUNCO(1, IBNDS), IREDUNCO(2, IBNDS),
     &                        IBNDS, TOTREDNCO, NRATMS)
            IF (CMP_DBMAT) CALL DBND(CARTCOORD, IREDUNCO(1,IBNDS),
     &                               IREDUNCO(2,IBNDS), IBNDS,
     &                               MAXREDUNCO, NRATMS, DERBMAT,
     &                               IREDUNCO, TOTREDNCO)
            REDUNCO(IBNDS) = DISTAB
C
#ifdef _DEBUG_LVL0
            WRITE(6,"(a,F10.5)") "The bond distance =",DISTAB*0.529177249d0
#endif                        
C
         ENDIF
C
 20   CONTINUE
C
      DO IANGS = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)

         IF (IREDUNCO(4, IANGS) .EQ. MONE) THEN
C
C This is for the linear arrangements. 
C
            CALL  BULT_LANGCRD(CARTCOORD, BMATRX, ANGL,
     &                         IREDUNCO(1, IANGS), IREDUNCO(2, IANGS),
     &                         IREDUNCO(3, IANGS), IREDUNCO(4, IANGS), 
     &                         IANGS, NRATMS, TOTREDNCO, EPSILON)
            REDUNCO(IANGS) = ANGL
#ifdef _DEBUG_LVL0
            WRITE(6,"(a,F10.5)") "The bond Angle =", ANGL/DINVPI
#endif
C
C The Derivative of B matrix is zero for the angle 0.0 and 180.00
CSSS            IF (CMP_DBMAT) CALL DEANG(CARTCOORD, IREDUNCO(1,IANGS),
CSSS     &                                IREDUNCO(2,IANGS),
CSSS     &                                IREDUNCO(3,IANGS), BMATRX, NRATMS,
CSSS     &                                TOTREDNCO, IREDUNCO, MAXREDUNCO,
CSSS     &                                IANGS, DERBMAT)
C
         ELSE
C
C This is for the non-linear arrangements. 
C
            CALL  BULT_ANGCRD(CARTCOORD, BMATRX, ANGL,
     &                       IREDUNCO(1, IANGS), IREDUNCO(2, IANGS),
     &                       IREDUNCO(3, IANGS), IANGS, TOTREDNCO,
     &                       NRATMS)
            REDUNCO(IANGS) = ANGL
C
#ifdef _DEBUG_LVL0
            WRITE(6,"(a,F10.5)") "The Bond Angle =", ANGL/DINVPI
#endif
            IF (CMP_DBMAT .AND. ANGL/DINVPI .LE. 175.0D0)
     &                     CALL DEANG(CARTCOORD, IREDUNCO(1,IANGS),     
     &                                IREDUNCO(2,IANGS),
     &                                IREDUNCO(3,IANGS), BMATRX, 
     &                                NRATMS, TOTREDNCO, IREDUNCO, 
     &                                MAXREDUNCO, IANGS, DERBMAT)
         ENDIF
C
      ENDDO
C
 30   CONTINUE
C
#ifdef _DEBUG_LVL0
            WRITE(6,*) "Entering the dihedral angle block"
#endif 
  
      DO 40 IDIHS = (TOTNOFANG + TOTNOFBND + 1),  TOTREDNCO
C
         CALL  BULT_DIHANGCRD(CARTCOORD, BMATRX, DANG,
     &                        IREDUNCO(1, IDIHS), IREDUNCO(2, IDIHS),
     &                        IREDUNCO(3, IDIHS), IREDUNCO(4, IDIHS),
     &                        IDIHS, TOTREDNCO, NRATMS)
         IF (CMP_DBMAT) CALL DDIH(CARTCOORD,
     &                            IREDUNCO(1,IDIHS), IREDUNCO(2,IDIHS),
     &                            IREDUNCO(3,IDIHS), IREDUNCO(4,IDIHS),
     &                            IDIHS, TOTREDNCO, NRATMS, IREDUNCO,
     &                            MAXREDUNCO, DERBMAT)
            REDUNCO(IDIHS) = DANG
C
#ifdef _DEBUG_LVL0
            WRITE(6,"(a,F10.5)") "The Dihedral Angle =", DANG/DINVPI
#endif
C    
 40   CONTINUE
C
C Let's write the DERBMAT to JOBARC
CSSS      PRINT *,'DERBMAT IN BUILT_DERBMAT'
CSSS      CALL OUTPUT(DERBMAT,1,3*NRATMS,1,3*NRATMS*TOTREDNCO,
CSSS     &            3*NRATMS,3*NRATMS*TOTREDNCO,1)
      IF (CMP_DBMAT) THEN
         LENGTH_DERBMAT = 9*NRATMS*NRATMS*TOTREDNCO
         CALL DPUTREC(1,'JOBARC','DERBMAT',LENGTH_DERBMAT,DERBMAT)
      END IF
C
C Let's write the B-matrix to the JOBARC file. 
C
cSSS     CALL DPUTREC(20,'JOBARC', 'BMATRIX ',TOTREDNCO, BMATRX)
C
C Form the G matrix which is required to transfrom the gradient
C from the Cartesian coordinates to redundent internal coordinates.   
C The procedures that are followed here is exactly what is described
C in Pulay et al., 96, 2856, 1992 and Peng et al., 17, 49. 1996.
C Built G = BuB(t) where u is an arbitrary nonsingular matrix, usually
C the unit matrix. The notations are consistent with Pulay et al.
C paper.
C
#ifdef _DEBUG_LVL0
      Write(6,*) "B-Matrix"
      CALL OUTPUT(BMATRX, 1, TOTREDNCO, 1, 3*NRATMS, TOTREDNCO,
     &            3*NRATMS, 1)    
#endif
C
      CALL DCOPY(3*NRATMS*TOTREDNCO, BMATRX, 1, BTMP, 1)
C
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, 3*NRATMS, 1.0D0,
     &           BMATRX, TOTREDNCO, BTMP, TOTREDNCO, 0.0D0, GMATRX,
     &           TOTREDNCO)
C
CSSS      Print*, "The G-Matrix:BB^t", MXATMS, MAXREDUNCO
CSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
CSSS     &            TOTREDNCO, 1)
      LENGMAT=TOTREDNCO*TOTREDNCO
      CALL DPUTREC(20,'JOBARC','G-MATRX ',LENGMAT,GMATRX)
C---DEBUG
cSSS      CALL ZERO(BMATRX, 3*NRATMS*TOTREDNCO)
cSSS      CALL DCOPY(TOTREDNCO*TOTREDNCO, GMATRX, 1, BMATRX, 1)
C---DEBUG
C
C The intermediate G matrix created above is linear dependent, hence
C there are eigenvalues that are zero. Invert the non-zero digonal
C elements to built the Lambda(-1) matrix.
C
      CALL EIG(GMATRX, EIGVECTORS, 1, TOTREDNCO , 1)
C 
CSSS      Print*, "The Eigenvalues of G-Matrix"
CSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
CSSS     &            TOTREDNCO, 1)
CSSS      Print*, "The Eigenvectors of G-Matrix"
CSSS      CALL OUTPUT(EIGVECTORS, 1, TOTREDNCO, 1, TOTREDNCO, 
CSSS     &            TOTREDNCO, TOTREDNCO, 1)
C
      NULLEVAL = 0
      DO I = 1, TOTREDNCO
         IF (GMATRX(I, I) .LE. EPSILON) THEN
             NULLEVAL = NULLEVAL + 1
             GMATRX(I, I) = 0.0D0
         ELSE
             GMATRX(I, I) = 1.0D0/GMATRX(I, I)
         ENDIF
      ENDDO
C  
      NONZERO_TOTREDNCO = TOTREDNCO - NULLEVAL
  
      CALL DPUTREC(20,'JOBARC','REDEVECS',NULLEVAL*TOTREDNCO,
     &            EIGVECTORS)
      CALL IPUTREC(20,'JOBARC','NUMREDCO', 1, NULLEVAL)
C
C Built the generalized inverse of G-matrix, what proceeded is generally
C known as singular value decomposition (ineversion of singular matrices)
C G^(-1) = [K L] Lambda^(-1)[K^(t) L^(t)]
C 
      CALL XGEMM('N', 'N', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
     &           EIGVECTORS, TOTREDNCO, GMATRX, TOTREDNCO, 0.0D0,
     &           BTMP, TOTREDNCO)
CSSS      Print*, "The INVERSE OF B-Matrix:1"
CSSS      CALL OUTPUT(BTMP, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
CSSS     &            TOTREDNCO, 1)
      CALL XGEMM('N', 'T', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
     &           BTMP, TOTREDNCO, EIGVECTORS, TOTREDNCO, 0.0D0,
     &           GMATRX, TOTREDNCO)
      CALL DPUTREC(20,'JOBARC','GI-MATRX',LENGMAT,GMATRX)
C---DEBULG
CSSS      CALL ZERO(BTMP, TOTREDNCO*TOTREDNCO)
CSSS      CALL XGEMM('N', 'N', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
CSSS     &           GMATRX, TOTREDNCO, BMATRX, TOTREDNCO, 0.0D0,
CSSS     &           BTMP, TOTREDNCO)
CSSS      CALL XGEMM('N', 'N', TOTREDNCO, TOTREDNCO, TOTREDNCO, 1.0D0,
CSSS     &           BTMP, TOTREDNCO, BMATRX, TOTREDNCO, 0.0D0,
CSSS     &           GMATRX, TOTREDNCO)

CSSS      Print*, "The CHECK INVERSE OF B-Matrix"
CSSS      CALL OUTPUT(GMATRX, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
CSSS     &            TOTREDNCO, 1)
C---DEBUG
C
C Built the intermediate matrix needed to built G and A matrices. 
C
      CALL XGEMM('N', 'N', TOTREDNCO, 3*NRATMS, TOTREDNCO, 1.0D0,
     &           GMATRX, TOTREDNCO, BMATRX, TOTREDNCO, 0.0D0,
     &           BTMP, TOTREDNCO)

CGAL  Here we build B(t)*G(-)
      CALL XGEMM('T','N',3*NRATMS,TOTREDNCO,TOTREDNCO,1.0D0,
     &            BMATRX,TOTREDNCO,GMATRX,TOTREDNCO,0.0D0,
     &            BTGMIN,3*NRATMS)
CSSS      Write(6,*) "BTGMIN matrix"
CSSS      CALL OUTPUT(BTGMIN,1,3*NRATMS,1,TOTREDNCO
CSSS     &            ,3*NRATMS,TOTREDNCO,1)

      LENBTGMIN=3*NRATMS*TOTREDNCO
      CALL DPUTREC(20,'JOBARC','BTGMIN',LENBTGMIN,BTGMIN)
C  
C Built the the G-MATRIX (G(-1)B)^t; Note the transpose 
C is necessary because the way the transformation is done in 
C CONVQ.f (ACES II). This is the transformation matrix, when
C act on right convert Cartesian gradients to internals, and
C act on left convert internal coordiantes to Cartesian. 
C   
      CALL TRANSP(BTMP, GMATRX, 3*NRATMS, TOTREDNCO)
C
C Make the transpose of the B-Matrix and  pass it out. We 
C need this to transform the Cartesian Hessian matrix to 
C internals.
C
      CALL TRANSP(BMATRX, BTMP, 3*NRATMS, TOTREDNCO)
      CALL DCOPY(3*NRATMS*TOTREDNCO, BTMP, 1, BMATRX, 1)  
C
CSSS
CSSS      CALL OUTPUT(GMATRX, 1, 3*NRATMS, 1, TOTREDNCO, 3*NRATMS,
CSSS     &            TOTREDNCO, 1) 
C
      RETURN
      END
