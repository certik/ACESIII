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
      SUBROUTINE PERTABLE
C
C     GIVES ATOMIC NUMBER AND ATOMIC MASS INFORMATION
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     1     PARNAM(3*MXATMS)
#include "coord.com"
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
 
#include "io_units.par"
 
C
C     Setup the periodic table - we have entered MXTABL elements here
C
      INTEGER MXTABL
      PARAMETER (MXTABL = 103)
CJDW 2/16/95. Out into the unknown. Add elements 87-103 !
      DOUBLE PRECISION ATMSS(MXTABL)
      CHARACTER*2 ATSYM(MXTABL)
C
      DATA ATSYM /'H ','HE',
     &            'LI','BE','B ','C ','N ','O ','F ','NE',
     &            'NA','MG','AL','SI','P ','S ','CL','AR',
     &            'K ','CA',
     &            'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     &                      'GA','GE','AS','SE','BR','KR',
     &            'RB','SR',
     &            'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',
     &                      'IN','SN','SB','TE','I ','XE',
     &            'CS','BA',
     &            'LA','CE','PR','ND','PM','SM','EU','GD',
     &            'TB','DY','HO','ER','TM','YB','LU',
     &                 'HF','TA','W ','RE','OS','IR','PT','AU','HG',
     &                      'TL','PB','BI','PO','AT','RN',
     &            'FR','RA',
     &            'AC','TH','PA','U ','NP','PU','AM','CM',
     &            'BK','CF','ES','FM','MD','NO','LR'               /
      DATA ATMSS / 1.007825D+00 , 4.00260D+00 , 7.01600D+00 ,
     &     9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,
     &     14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,
     &     19.99244D+00 , 22.98980D+00 , 23.98504D+00 ,
     &     26.98153D+00 , 27.97693D+00 , 30.97376D+00 ,
     &     31.97207D+00 , 34.96885D+00 , 39.96238D+00  ,
     &     38.96371D+00 , 39.96259D+00 , 44.95591D+00  ,
     &     47.94795D+00 , 50.94396D+00 , 51.94051D+00  ,
     &     54.93805D+00 , 55.93494D+00 , 58.93320D+00  ,
     &     57.93535D+00 , 62.93960D+00 , 63.92915D+00  ,
     &     68.92558D+00 , 73.92118D+00 , 74.92159D+00  ,
     &     79.91650D+00 , 78.91830D+00 , 83.91150D+00  ,
     &     84.91180D+00 , 87.90560D+00 , 88.90580D+00  ,
     &     89.90470D+00 , 92.90640D+00 , 97.90540D+00  ,
     &     00.00000D+00 , 101.9043D+00 , 102.9055D+00  ,
     &     105.9032D+00 , 106.9050D+00 , 113.9036D+00  ,
     &     114.9041D+00 , 117.9034D+00 , 120.9038D+00  ,
     &     129.9067D+00 , 126.9044D+00 , 131.9042D+00  ,
     &     132.9054D+00 , 137.9052D+00 , 138.9063D+00  ,
     &     139.9054D+00 , 140.9076D+00 , 141.9077D+00  ,
     &     144.9127D+00 , 151.9197D+00 , 152.9212D+00  ,
     &     157.9241D+00 , 158.9253D+00 , 163.9292D+00  ,
     &     164.9303D+00 , 165.9320D+00 , 168.9342D+00  ,
     &     173.9389D+00 , 174.9408D+00 , 179.9465D+00  ,
     &     180.9480D+00 , 183.9509D+00 , 186.9557D+00  ,
     &     191.9615D+00 , 192.9629D+00 , 194.9648D+00  ,
     &     196.9665D+00 , 201.9706D+00 , 204.9744D+00  ,
     &     207.9766D+00 , 208.9804D+00 , 208.9824D+00  ,
     &     209.9875D+00 , 222.0157D+00 , 223.0197D+00  ,
     &     226.0254D+00 , 227.0277D+00 , 232.0381D+00  ,
     &     231.0359D+00 , 238.0508D+00 , 237.0482D+00  ,
     &     244.0642D+00 , 243.0614D+00 , 247.0703D+00  ,
     &     247.0703D+00 , 251.0796D+00 , 252.0829D+00  ,
     &     257.0751D+00 , 258.0986D+00 , 259.1009D+00  ,
     &     260.1053D+00 /

      integer  atomnumb
      external atomnumb
      iBig_A = iachar('A')
      iLit_a = iachar('a')
C
C     For each entry in the Z-matrix, get a mass
C
      IF (IPRNT .GE. 10) WRITE (LUOUT, 9000)
 9000 FORMAT ('PERTABLE: Mass and at. nr. lookup'/
     1     'Line Symbol AtNr   At. Mass')
      NMPROT=0
      IERR=0
      DO 10 I = 1, NATOMS
         IF ('a'.LE.ZSYM(I)(1:1)) THEN
            ZSYM(I)(1:1)=ACHAR(iBig_A+IACHAR(ZSYM(I)(1:1))-iLit_a)
         END IF
         IF ('a'.LE.ZSYM(I)(2:2)) THEN
            ZSYM(I)(2:2)=ACHAR(iBig_A+IACHAR(ZSYM(I)(2:2))-iLit_a)
         END IF
         J = ATOMNUMB(ZSYM(I)(1:2))
         IF (J.NE.0) THEN
            IATNUM(I) = J
            ATMASS(I) = ATMSS(J)
            NMPROT = NMPROT + J
         ELSE
            IF (ZSYM(I).EQ.'X') THEN
               IATNUM(I) = 0
               ATMASS(I) = 0.0
            ELSE IF (ZSYM(I).EQ.'GH') THEN
               IATNUM(I) = 110
               ATMASS(I) = 100.0
            ELSE
               PRINT *, '@PERTABLE: Unknown element ',ZSYM(I)
               IERR = 1
            END IF
         END IF
         IF (IPRNT .GE.10)
     $      WRITE (LuOut,9010) I,ZSYM(I),IATNUM(I),ATMASS(I)
 9010    FORMAT(I3,2X,A5,2X,I5,2X,F10.5)
 10   CONTINUE
      IF (IERR.NE.0) CALL ERREX
      IONE=1
      CALL IPUTREC(20,'JOBARC','NMPROTON',IONE,NMPROT)
      RETURN
      END
