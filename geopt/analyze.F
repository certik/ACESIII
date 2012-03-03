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
      SUBROUTINE ANALYZE(SCRATCH,A)
C
C ANALYZES INTERNAL COORDINATE STRUCTURE, OBTAINS NUMBER OF
C CONFORMATIONAL DEGREES OF FREEDOM AND CHECKS FOR OBVIOUS
C REDUNDANCIES AND POTENTIAL ERRORS IN THE Z-MATRIX INPUT.
C SCRATCH IS THE SAME VECTOR USED IN SUBROUTINE SYMMETRY.
C THIS REQUIRES THE TRANSFORMATION MATRIX FROM CARTESIAN
C DERIVATIVES TO INTERNAL DERIVATIVES, AND USES A SIMPLE
C TOTALLY SYMMETRIC PROPERTY (CLOSELY RELATED TO THE NUCLEAR
C REPULSION ENERGY) AS A REFERENCE FUNTION.  NORD IS A AN
C INTEGER SCRATCH ARRAY USED ONLY HERE AND IN DEPENDENTS.
C
C This routine is cleaned up to generate messages that
C make sense for redundent internals. Some of the analysis
C done here is not useful for RICs (not valid). 
C Ajith Perera, 04/2006. 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "flags.h"
#include "jodaflags.com"
#include "cbchar.com"
#include "coord.com"
c
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
C
      DIMENSION SCRATCH(6*NATOMS),NORD(6*MXATMS),A(3*NATOMS,NXM6)
      Character*5 ChScr(3*MxAtms)
      LOGICAL RIC_INUSE
      DOUBLE PRECISION ICHGI,ICHGJ
C
      IDIDIT=0
      NDFU=0
      NPFU=0
      ITRASH=0
      IPNDAP=0
      IERR=0
      IERROR=0
      ZREPUL=0.D0
      ICOORE=0
      RIC_INUSE = ((iFlags2(h_IFLAGS2_geom_opt) .EQ. 3  .OR.
     &              iFlags2(h_IFLAGS2_geom_opt) .EQ. 4) .AND.
     &              iFlags(h_IFLAGS_coord)      .EQ. 1)
C
C LOOP THROUGH CARTESIAN COORDINATES AND OBTAIN DERIVATIVES
C OF THE REFERENCE FUNCTION WRT TO CARTESIAN DISPLACEMENTS.  DONE
C ANALYTICALLY.
C
      CALL ZERO(SCRATCH,6*NATOMS)
#ifdef _DEBUG_LVL0
      Write(6,*)
      WRITE(6,*)' R VECTOR IN ANALYZE '
      Write(6,*)
      WRITE(6,'(3f20.12)')(R(I),I=1,NXM6)
14    FORMAT(3(2X,I2))
      Write(6,*)
      WRITE(6,"(6I4)")(IATNUM(J),J=1,NATOMS) 
#endif
      NOPT0=NOPT
      NDERIV=0
      WRITE(6,150)
      WRITE(6,149)
      WRITE(6,150)
      CALL NUCREP(ZREPUL,ICOORE)
C
      IF(ICOORE.NE.0)THEN
       WRITE(6,547)
       WRITE(6,158)
       WRITE(6,150)
      ELSE
       WRITE(6,546)ZREPUL
       call dputrec(1, ' ', 'NUCREP', 1, zrepul )
      ENDIF
C
149   FORMAT(T3,' Analysis of internal coordinates specified by '
     &,'Z-matrix ')
150   FORMAT(80('-'))
158   FORMAT(T3,' Thank you for using the ACES2 program system.')
546   FORMAT(T3,' *The nuclear repulsion energy is ',f10.5,' a.u.')
547   FORMAT(T3,' **PROGRAM ERROR** Problem with coordinate vector.')
C
      CALL VADD(SCRATCH,SCRATCH,Q,3*NATOMS,1.D0)
C
C EVALUATE DERIVATIVES IN CARTESIAN REPRESENTATION.
C DIFFERENTIATED FUNCTION IS:
C
C         SUM    [CHARGE(I)*CHARGE(J)]**0.50/d(I,J)
C         I<J
C (TAKING ROOT REDUCES SIZE OF SCALING FACTOR
C  FOR HEAVY--HEAVY INTERACTIONS)
C
      DO 10 I=1, 3*NATOMS
        IXYZ=MOD(I,3)
        IF(IXYZ.EQ.0)IXYZ=3
        ICUT=3-IXYZ
        IATOMI=(I+2)/3
        IF(IATNUM(IATOMI).EQ.0)GOTO 10
        ICHGI=DFLOAT(IATNUM(IATOMI))**0.50
        JBOT=3-ICUT
        ZARF=0.D0
        DO 50 J=JBOT,3*NATOMS,3
         IF(J.EQ.I)GOTO 50
         IATOMJ=(J+2)/3
         IF(IATNUM(IATOMJ).EQ.0)GOTO 50
         ICHGJ=DFLOAT(IATNUM(IATOMJ))**0.50
         DISIJ=DIST(Q(3*IATOMI-2),Q(3*IATOMJ-2))
         FACTOR=ICHGI*ICHGJ/DISIJ**3
         ZARF=ZARF+FACTOR*(SCRATCH(I)-SCRATCH(J))
 50     CONTINUE
        SCRATCH(3*NATOMS+I)=ZARF
 10   CONTINUE
C 
#ifdef _DEBUG_LVLM1
      Write(6,*)
      WRITE(6,*)' Cartesian derivatives of scaled nulear repulsion.'
      Write(6,*)
      WRITE(6,80)(SCRATCH(3*NATOMS+J),J=1,NXM6)
      Write(6,*)
#endif
C
C NOW DETERMINE TOTAL NUMBER OF DEGREES OF FREEDOM WITHIN THE
C TOTALLY SYMMETRIC SUBSPACE.
C
      ITRNDF=ITNDF()
C
C DERIVATIVE VECTOR NOW IN TOP HALF OF SCRATCH.  TRANSFORM IT
C TO INTERNAL COORDINATES, PLACED IN LOWER PART OF SCRATCH.
C
      CALL ZERO(SCRATCH,3*NATOMS)
      DO 35 J=1,NXM6
      Z1=0.D0
      DO 15 I=1,3*NATOMS
15    Z1=SCRATCH(3*NATOMS+I)*A(I,J)+Z1
35    SCRATCH(J)=Z1
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      WRITE(6,*)' INT COOR DERIVATIVES '
      Write(6,*)
      WRITE(6,80)(SCRATCH(J),J=1,NXM6)
      WRITE(6,*)
#endif
C
C NOW GO THROUGH LIST OF OPTIMIZED PARAMETERS, COUNTING NUMBER OF
C PROBABLE POSITIVE/NEGATIVE DIHEDRAL ANGLE PAIRS.  THEN
C EVENTUALLY (AT THE BOTTOM OF THE 30 DO LOOP)
C DECREMENT NUMBER OF "Z-MATRIX DEGREES OF FREEDOM" (NUMBER
C OF "INDEPENDENTLY" ADJUSTED PARAMETERS) BY THIS NUMBER.
C
      IF ((iFlags2(h_IFLAGS2_geom_opt) .EQ. 1)) THEN
C
C This test is meaningless for Cartesian opts.
C
      DO 423 I=1,NOPT
       ISTP=0
       DO 424 J=I+1,NOPT
        IF(ISTP.NE.0)GOTO 424
        ZIF=DABS(R(ISQUASH(NOPTI(I)))+R(ISQUASH(NOPTI(J))))
        IF(ZIF.LT.1.D-6)THEN
         ZIF2=DABS(SCRATCH(NOPTI(I))+SCRATCH(NOPTI(J)))
         IF(ZIF2.LT.1.D-6.AND.DABS(R(ISQUASH(NOPTI(J))))
     &   .GT.1.D-6)THEN
          IPNDAP=IPNDAP+1
          ISTP=1
         ENDIF
        ENDIF
424    CONTINUE
423   CONTINUE
C
C Endif of iFlags2(h_IFLAGS2_geom_opt) .EQ.1 
C
      ENDIF
C
C SORT INTERNAL COORDINATE DERIVATIVE VECTOR AND THEN ANALYZE IT.
C
      DO 20 I=1,NXM6
      NORD(I)=I
20    SCRATCH(I)=DABS(SCRATCH(I))
      CALL PIKSR2(NXM6,SCRATCH,NORD)
#ifdef _DEBUG_LVLM1
      WRITE(6,*)' SORTED INT COOR DERIVS '
      Write(6,*)
      WRITE(6,80)(SCRATCH(J),J=1,NXM6)
C      WRITE(6,*) (NORD(I), I=1,NXM6)
C      WRITE(6,*) (ISQUASH(I), I=1,NXM6)
      WRITE(6,*)
#endif
      NDF=0
C      IF(SCRATCH(1).GT.1.D-6)NDERIV=1
C
C FIRST CHECK TO SEE IF ALL INTERNAL COORDINATES WITH THE
C SAME NAME ARE APPARENTLY EQUIVALENT.
C
      DO 501 K = 1, NOPT
         DO 502 I = 1, NXM6-1
         IF (VARNAM(ISQUASH(NORD(I))).EQ.PARNAM(K)) THEN
            DO 503 J = I+1, NXM6
            IF (VARNAM(ISQUASH(NORD(J))).EQ.PARNAM(K)) THEN
               DENOM = DABS(SCRATCH(I))
               IF (DENOM.LT.1.D-12) DENOM = 1.D0
               CRIT = DABS(SCRATCH(I)-SCRATCH(J))/DENOM
               IF (CRIT.GE.1.D-4) THEN
                  WRITE(6,155) NORD(I),NORD(J),VARNAM(ISQUASH(NORD(I)))
                  IERROR = 1
               END IF
            END IF
 503        CONTINUE
         END IF
 502     CONTINUE
 501  CONTINUE
C
      IF(IERROR.EQ.1)THEN
        WRITE(6,311)
        WRITE(6,158)
        WRITE(6,150)
        Call ErrEx
      ENDIF
C
155   FORMAT(T3,'*ERROR* Parameters ',i2,' and ',i2,' are not '
     &,'equivalent [both called ',a5,'].')
311   FORMAT(T3,' *Reconstruct Z-matrix and try again.')
C
C NOW LOOP THROUGH PARAMETERS WHICH ARE TO BE OPTIMIZED TO
C MAKE SURE THAT ALL OF THESE CAN HAVE NON-ZERO GRADIENTS.
C
      IF ((iFlags2(h_IFLAGS2_geom_opt) .EQ. 1)) THEN
C
C This test is meaningless for pure Cartesian optmizations.
C
      DO 840 I=1,NOPT
      CHSCR(1)=PARNAM(I)
       IERR=0
       DO 841 J=1,NXM6
       IF(VARNAM(ISQUASH(NORD(J))).EQ.CHSCR(1).AND.
     &    SCRATCH(J).LT.1.D-8)THEN
          WRITE(6,159)VARNAM(ISQUASH(NORD(J))),
     &    NORD(J)
       IERR=IERR+1
       ENDIF
 841   CONTINUE
C
       IF(IERR.GE.1)THEN
       ITRASH=ITRASH+1
       ENDIF
 840   CONTINUE
C
      ENDIF
C
      IBOT=1
      NPNDPC=0
      NOPNPC=0
      DO 30 I=1,NXM6
      IQUIT=0
      DENOM=DABS(SCRATCH(I))
      IF(DENOM.GT.1.D-15)THEN
       ZIFF=DABS(SCRATCH(I+1)-SCRATCH(I))/DENOM
      ELSE
       ZIFF=DABS(SCRATCH(I+1)-SCRATCH(I))
      ENDIF
C
159   FORMAT(T3,'*WARNING* Parameter ',a5,'[',i3,'] cannot have '
     &,'a nonzero derivative,',/,T4,' but is being optimized. '
     &,' Reconstruction of Z-matrix recommended.')
C
      IF ((iFlags2(h_IFLAGS2_geom_opt) .EQ. 1)) THEN
C
C This test is meaningless for full Cartesian opts.
C
      IF(ZIFF.LT.1.D-4)THEN
         IF(DABS(SCRATCH(I)).GT.1.D-6)THEN 
C
C CHECK NOW FOR APPARENT EQUIVALENCES BETWEEN PARAMETERS HAVING
C DIFFERENT NAMES.  IF AN EVENT IS ENCOUNTERED, CHECK FIRST TO
C SEE IF IT IS A POS/NEG DIHEDRAL ANGLE PAIR.  IF IT ISN'T PRINT
C WARNING MESSAGE.  R IS UNSQUASHED IN THIS ROUTINE.
C
        IF(IPRNT.GE.100)WRITE(6,154)
     &     VARNAM(ISQUASH(NORD(I+1))),
     &     VARNAM(ISQUASH(NORD(IBOT)))
C
        IF(VARNAM(ISQUASH(NORD(I+1))).NE.
     &   VARNAM(ISQUASH(NORD(IBOT))))THEN
         CRIT=DABS(R(ISQUASH(NORD(I+1)))+R(ISQUASH(NORD(IBOT))))
         IF(CRIT.LT.1.D-8 .AND.  .NOT. RIC_INUSE)THEN
          WRITE(6,169)VARNAM(ISQUASH(NORD(I+1))),
     &    VARNAM(ISQUASH(NORD(IBOT)))
         ELSE
C
C THIS IS PROBABLY A POS/NEG DIHEDRAL ANGLE PAIR.  THIS IS NOT
C A PROBLEM. EXIT LOOP.
C
         ENDIF
        ENDIF
       ENDIF
C Endif for (iFlags2(h_IFLAGS2_geom_opt) .EQ. 1)
       ENDIF
C
154   FORMAT(2X,A5,2X,A5,2X,'(*')
169   FORMAT(T3,'*WARNING* Parameters ',a5,' and ',a5,' appear '
     &,'to be equivalent.')
C
      ELSE
C
C NEW VARIABLE.  (NPNDAP IS NUMBER
C OF POSITIVE/NEGATIVE DIHEDRAL ANGLE PAIRS WHICH CAN HAVE NONZERO
C DERIVATIVES, NOPNDP IS NUMBER OF THESE THAT ARE OPTIMIZED).
C INCREMENT NDERIV (THE NUMBER OF DISTINCT PARAMETERS WHICH CAN
C HAVE NONVANISHING GRADIENTS).
C
       IBOT=I+1
       IDID=0
       IDID2=0
       IF(SCRATCH(I).GT.1.D-8)NDERIV=NDERIV+1
       IF(IPRNT.GE.100)WRITE(6,157)NDERIV,I
       NORD(NX+NDERIV)=NORD(I)
       IATOM=(NORD(I)+8)/3
       IF(NORD(I).EQ.1)IATOM=2
       IF(IPRNT.GE.100)WRITE(6,*)IATOM,ZIFF,SCRATCH(I)
       IF(ATMASS(IATOM).LT.1.D-3)THEN
C
C DUMMY ATOM.  NOW SEE IF THERE IS AN EQUIVALENT VALUE WHICH CORRESPONDS
C TO ONE OF THE ATOMIC ICs.  RIGHT NOW, THIS INFORMATION IS NOT USED
C EXTENSIVELY, BUT IT MAY PROVE USEFUL WHEN WE FIGURE OUT A WAY TO
C DETERMINE HOW MANY OF THE Z-MATRIX PARAMETERS HAVE TO BE VARIED TO
C RELEASE ALL CONSTRAINTS ON THE OPTIMIZATION (NOT ALWAYS EQUAL TO ITRN
C
        NDFU=NDFU+1
        IPSO=0
        DO 33 J=I-1,1,-1
         IF(IPSO.EQ.1)GOTO 33
         ZIFF2=DABS(SCRATCH(J)-SCRATCH(I))
         IF(ZIFF2.GT.1.D-5)GOTO 33
         JATOM=(NORD(J)+8)/3
         IF(NORD(J).EQ.1)JATOM=2
         IF(ATMASS(JATOM).GT.1.D-3)THEN
          NDF=NDF+1
          IPSO=1
          ZREF=R(NORD(I))
         ENDIF
33      CONTINUE
        GOTO 30
       ENDIF
       IF(SCRATCH(I).LT.1.D-8)GOTO 30
       ZREF=R(NORD(I))
       NDF=NDF+1
       NDFU=NDFU+1
       IF(IPRNT.GE.100)WRITE(6,88)I,NDF
      ENDIF
30    CONTINUE
C
      NOPT0=NOPT-IPNDAP
C
 88   FORMAT(T3,' On cycle ',i2,' NDF incremented to ',i3)
157   FORMAT(T3,' NDERIV incremented to ',i2,' on cycle ',i3)
C
C NOW PRINT OUT WHAT WE HAVE LEARNED.
C
C WRITE(6,170)NDF (WE DON'T KNOW HOW TO DETERMINE THIS YET.)
C170  FORMAT(T3,' *In the Z-matrix, there are ',i3,' parameters which '
C ,'must be optimized.')
C
      IF(ITRNDF.EQ.1)WRITE(6,152)ITRNDF
C
C This is somewhat cheating. In the case of RIC we assume all the
C degrees of freedoms that are optimized belong to the totally 
C symmetric subspace (until we develop a routine similar to ITNDF
C that works with RICs not with Cartesian directions per symmetry
C equivalent atoms. It can be done. All we need is to find out
C how many uniques RICs and then subtract the rotational and vibrational
C degrees of freedoms that transform as totaly symm. Unfotunately, 
C I don't see the value of it since redundent simply implies 
C that we are optimizing more than minimum required. Ajith Perera,04/2006. 
C
      IF(ITRNDF.GT.1 .AND. .NOT. RIC_INUSE) THEN
        WRITE(6,151)ITRNDF
      ELSE
        WRITE(6,151)NOPT
      ENDIF
C
      WRITE(6,171)NOPT
      IF(IPNDAP.GT.0)THEN
         WRITE(6,545)IPNDAP,NOPT0
      ENDIF
      ICON=2
C
      IF (RIC_INUSE) THEN
          ICON =0
          WRITE(6,176) 
      ELSE
C 
         IF(ITRNDF.EQ.NOPT0.AND.ITRASH.EQ.0.AND.NDERIV.EQ.ITRNDF)THEN
            WRITE(6,172)
            ICON=0
         ELSEIF(ITRNDF.EQ.NOPT0.AND.ITRASH.EQ.0.AND.NDERIV.GT.
     &          ITRNDF)THEN
            WRITE(6,142)
            ICON=11
         ELSEIF(ITRNDF.GT.NOPT0)THEN
            ICON=-1
            WRITE(6,173)
         ELSEIF(ITRNDF.LT.NOPT0.AND.NDERIV.EQ.NOPT0)THEN
            ICON=1
            WRITE(6,174)
         ELSEIF(ITRNDF.LT.NOPT0.AND.NDERIV.GT.NOPT0)THEN
            ICON=6
            WRITE(6,178)
         ENDIF
      ENDIF
C
142   FORMAT(T3,' *Although the number of optimized parameters is '
     &,'equal to the true ',/,t4,' number of degrees of freedom, '
     &,' the optimization may be constrained.')
151   FORMAT(T3,' *There are ',i3,' degrees of freedom within the '
     &,'tot. symm. molecular subspace.')
152   FORMAT(T3,' *There is ',i3,' degree of freedom within the '
     &,'tot. symm. molecular subspace.')
171   FORMAT(T3,' *Z-matrix requests optimization of ',i3,
     &' coordinates.')
172   FORMAT(T3,' *The optimization is unconstrained and your Z-matrix '
     &,'is great.')
173   FORMAT(T3,' *The optimization is constrained.')
174   FORMAT(T3,' *The optimization problem may be overdetermined:')
176   FORMAT(T3,' *The RIC optimization is unconstrained and your',  
     &       ' input is great.')
178   FORMAT(T3,' *The optimization problem is poorly defined by '
     &,'your Z-matrix, and may be',/,T4,'overdetermined.')
545   FORMAT(T3,' *There are ',i3,' positive-negative dihedral '
     &,'angle pairs, ',/,t4,' resulting in ',i3,' Z-matrix '
     &,'degrees of freedom.')
C
      IF(ICON.NE.0)THEN
C 
       WRITE(6,180)NDERIV
       WRITE(6,181)(VARNAM(ISQUASH(NORD
     &    (3*NATOMS+J))),NORD(3*NATOMS+J),J=1,NDERIV)
       WRITE(6,183)NOPT
       WRITE(6,181)(PARNAM(K),NOPTI(K),K=1,NOPT)
C
#ifdef _DEBUG_LVLM1
       Print*, "NORD with offset"
       Write(6,*)
       Print*, (NORD(3*NATOMS+I), I=1, NXM6)
       Write(6,*)
#endif
C
180   FORMAT(T3,' *The following ',I3,' parameters can have non-zero',
     &      'derivatives within the ',/,t4, ' totally symmetric ',
     &       'subspace:')
181   FORMAT((T14,6(A5,'[',I3,']',2x)))
183   FORMAT(T3,' *The following ',I3,' parameters are to be '
     & ,'optimized:')
C
C BLOW OFF THIS LOOP IF NPFU.NE.0 -  THIS IS LIKELY TO BE UNUSUAL
C (CORRESPONDING TO OPTIMIZING ONLY ONE PART OF AN EQUIVALENT
C POS/NEG DIHEDRAL PAIR).  WITH MORE THOUGHT, THIS LOOP CAN BE
C ADJUSTED FOR NONZERO NPFUs, BUT LEAVE IT OUT FOR NOW.
C
       IF(NPFU.EQ.0)THEN
        IF(ICON.EQ.11)WRITE(6,198)
        IF(ICON.EQ.-1.AND.NDF.NE.NDFU)WRITE(6,191)ITRNDF-NOPT0
C
        IF(ICON.EQ.11)THEN
         WRITE(6,190)
         WRITE(6,290)
        ENDIF
C
        IF(ICON.EQ.1)THEN
         WRITE(6,193)NOPT0-ITRNDF
        ENDIF
C
        IF(ICON.EQ.6)THEN
         WRITE(6,193)NOPT0-ITRNDF
         WRITE(6,190)
        ENDIF

190   FORMAT(T3,' *Care should be taken that the following coordinates '
     &,'are dependent ',/,t4,' upon the set of optimized parameters.')
193   FORMAT(T3,' *Among the set of internal coordinates listed above, '
     &,'it is possible that ',/,t4,' as many as ',i3,' may be fixed ')
198   FORMAT(T3,' *There is room for improvement in the Z-matrix.')
290   FORMAT(T3,' *If this is true, optimization is unconstrained.')
191   FORMAT(T3,' *Among the set of internal coordinates listed below, '
     &,'at least ',/,t4,i3,' must be allowed to vary.  Those not '
     &,'varied must not be independent',/,t4,' coordinates.')
C
        IF(ICON.EQ.-1.AND.NDF.EQ.NDFU)WRITE(6,192)
        IDEP=0
C
        DO 890 I=1,NDERIV
         IDP=0
         DO 891 J=1,NOPT
891      IF(VARNAM(ISQUASH(NORD(3*NATOMS+I))).NE.PARNAM(J))IDP=IDP+1
         IF(IDP.EQ.NOPT)THEN
          IDEP=IDEP+1
          NORD(IDEP)=NORD(3*NATOMS+I)
          CHSCR(IDEP)=VARNAM(ISQUASH(NORD(3*NATOMS+I)))
         ENDIF
890     CONTINUE
C
        WRITE(6,181)(CHSCR(JX),NORD(JX),JX=1,IDEP)
        WRITE(6,150)
        IDIDIT=1
C
       ENDIF
      ENDIF
      IF(IDIDIT.EQ.0)WRITE(6,150)
C
80    FORMAT(3(2X,F20.14))

 192   FORMAT(T3,' *The following coordinates must be varied in an '
     &,' unconstrained optimization.')
C
      RETURN
      END
C
CSSS129   FORMAT(T3,' *Parameters ',a5,' and ',a5,' are '
CSSS     &     ,'a pos/neg dihedral angle pair.')
CSSS141   FORMAT(T3,' *The optimization will be unconstrained only if '
CSSS     &,'the following parameters',/,t4,' do not represent independent '
CSSS     &,'coordinates.')
CSSS145   FORMAT(T3,' *All of these are optimized, meaning that '
CSSS     &,'that there ',/,t4,' are ',i3,' Z-matrix degrees of '
CSSS     &,'freedom.')
CSSS192   FORMAT(T3,' *The following coordinates must be varied in an '
CSSS     &,' unconstrained optimization.')
CSSS345   FORMAT(T3,' *Of these,',i3,' are optimized, meaning that '
CSSS     &,'that there ',/,t4,' are ',i3,' Z-matrix degrees of '
CSSS     &,'freedom.')
