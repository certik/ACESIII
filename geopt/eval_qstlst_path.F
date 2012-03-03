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

      SUBROUTINE EVAL_QSTLST_PATH(QST_TANGENT, LST_TANGENT, WORK,
     &                            LST, QST)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

#include "mxatms.par"
#include "machsp.com"
#include "coord.com"

      LOGICAL XYZIN, NWFINDIF, QST, LST
      DOUBLE PRECISION LST_NORM_SQR, LST_NORM_INV, LST_TANGENT

      COMMON /USINT/NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/XYZIN,NWFINDIF

      DIMENSION WORK(NX*NX), LST_TANGENT(3*NATOMS),
     &          QST_TANGENT(3*NATOMS)
C
C Use Eqn. 1 and 8 in Schlegel to evaluate the QST and LST normalized
C tangent vectors.
C
C Read the guess structures (user provided) for the transition
C state and the products.
C
c      Write(6,*) (Q(I), I = 1, 3*NATOMS)
      IF (XYZIN) THEN
         CALL DGETREC(1,'JOBARC','PRSTRUCT',3*NATOMS,WORK)
      ELSE

         IRECTNT = 1
         IPRDUCT = IRECTNT + 3*NATOMS

         CALL DGETREC(1,'JOBARC','PRSTRUCT',3*NATOMS,
     &               WORK(IPRDUCT))
         CALL CNVRT2RAD(WORK(IPRDUCT),3*NATOMS)
         CALL SQUSH(WORK(IPRDUCT),3*NATOMS)
c         Write(6,*) "Debug-Info: The Product Structure"
c         Write(6,*) (Work(IPRDUCT + I), I= 0,3*NATOMS-1)

         IF (LST) THEN
            CALL DCOPY(3*NATOMS, R, 1, WORK(IRECTNT), 1)
c            Write(6,*) "Debug-Info: Current Structure"
c            Write(6,*) (Work(IRECTNT + I), I= 0, NATOMS)
            CALL XDAXPY(3*NATOMS, -1.0D0, WORK(IRECTNT), 1,
     &                 WORK(IPRDUCT), 1)
c            Write(6,*) "Debug-Info: The Diff Structure"
c            Write(6,*) (Work(IPRDUCT + I), I= 0, NATOMS)

CSSS            CALL DCOPY(NX, WORK(IPRDUCT), 1, LST_TANGENT, 1)
CSSS            CALL DGETREC(20, "JOBARC", "BMATRIX ", 3*NATOMS*NOPT,
CSSS     &                  WORK(1))
CSSS            IOFFSET = 3*NATOMS*NOPT + 1
CSSS            CALL XGEMM("N", "N", NOPT, 1, 3*NATOMS, 1.0D0, WORK,
CSSS     &                  NOPT, LST_TANGENT, 3*NATOMS, 0.0D0, WORK
CSSS     &                  (IOFFSET), NOPT)
CSSS            CALL DCOPY(NOPT, WORK(IOFFSET), 1, LST_TANGENT, 1)

            CALL SYMUNQONLY(WORK(IPRDUCT), LST_TANGENT)
            CALL NORMAL(LST_TANGENT, NOPT)

c            Write(6,*) "Debug-Info: LST tangent"
c            Write(6,*) (LST_TANGENT(I), I= 1, NOPT)

         ELSE IF (QST) THEN

            ICURENT = IPRDUCT + 3*NATOMS

            CALL DGETREC(1,'JOBARC','RXSTRUCT',3*NATOMS,
     &                  WORK(IRECTNT))
            CALL CNVRT2RAD(WORK(IRECTNT), 3*NATOMS)
            CALL SQUSH(WORK(IRECTNT), 3*NATOMS)
            CALL DCOPY(3*NATOMS, R, 1, WORK(ICURENT), 1)

c            Write(6,*) "Debug-Info: Rectant Structure"
c            Write(6,*) (Work(IRECTNT + I), I= 0, 3*NATOMS)
c            Write(6,*) "Debug-Info: Product Structure"
c            Write(6,*) (Work(IPRDUCT + I), I= 0, 3*NATOMS)
c            Write(6,*) "Debug-Info: Current Structure"
c            Write(6,*) (Work(ICURENT + I), I= 0, 3*NATOMS)

            CALL XDAXPY(3*NATOMS, -1.0D0, WORK(ICURENT), 1,
     &                 WORK(IPRDUCT), 1)
            CALL XDAXPY(3*NATOMS, -1.0D0, WORK(ICURENT), 1,
     &                 WORK(IRECTNT), 1)

            QST_NORM_SQR_PX = DDOT(3*NATOMS, WORK(IPRDUCT), 1,
     &                             WORK(IPRDUCT), 1)
            QST_NORM_SQR_RX = DDOT(3*NATOMS, WORK(IRECTNT), 1,
     &                             WORK(IRECTNT), 1)
            IF ((QST_NORM_SQR_PX.LE.1.0D-14).OR.
     &          (QST_NORM_SQR_RX.LE.1.0D-14)    ) THEN
               PRINT*, "AN INPUT ERROR: INSPECT THE INPUT STRUCTURES"
c               CALL ACES_EXIT
                call errex
            END IF

            QST_NORM_PX_INV = 1.0D0/(QST_NORM_SQR_PX)
            QST_NORM_RX_INV = 1.0D0/(QST_NORM_SQR_RX)

            CALL DSCAL(3*NATOMS, QST_NORM_PX_INV, WORK(IPRDUCT), 1)
            CALL DSCAL(3*NATOMS, QST_NORM_RX_INV, WORK(IRECTNT), 1)

            CALL XDAXPY(3*NATOMS, -1.0D0, WORK(IRECTNT), 1,
     &                 WORK(IPRDUCT), 1)

CSSS            CALL DCOPY(3*NATOMS, WORK(IPRDUCT), 1, QST_TANGENT, 1)
CSSS            CALL DGETREC(20, "JOBARC", "BMATRIX ", 3*NOPT*NATOMS,
CSSS     &                  WORK(1))
CSSS            IOFFSET = 3*NATOMS*NOPT + 1
CSSS            CALL XGEMM("N", "N", NOPT, 1, 3*NATOMS, 1.0D0, WORK,
CSSS     &                  NOPT, QST_TANGENT, 3*NATOMS, 0.0D0, WORK
CSSS     &                  (IOFFSET), NOPT)

            CALL SYMUNQONLY(WORK(IPRDUCT), QST_TANGENT)
            CALL NORMAL(QST_TANGENT, NOPT)
c            Write(6,*) "Debug-Info: QST tangent"
c            Write(6,*) (QST_TANGENT(I), I= 1, NOPT)

         END IF
      END IF

CSSS      LST_NORM_SQR = DDOT(3*NATOMS, WORK(IPRDUCT), 1,
CSSS     &                              WORK(IPRDUCT), 1)
CSSS      IF (DSQRT(LST_NORM_SQR).LE.1.0D-14) THEN
CSSS         PRINT*, "AN INPUT ERROR: INSPECT THE INPUT STRUCTURES"
CSSS         CALL ACES_EXIT
CSSS      END IF
CSSS      LST_NORM_INV = 1.0D0/DSQRT(LST_NORM_SQR)
CSSS      CALL DSCAL(3*NATOMS, LST_NORM_INV, WORK(IPRDUCT), 1)

CSSS      ASQURD = (QST_NORM_SQR_RX)*(QST_NORM_SQR_PX)/
CSSS     &         (QST_NORM_SQR_RX + QST_NORM_SQR_PX - 2.0D0*
CSSS     &          DSQRT((QST_NORM_SQR_RX)*(QST_NORM_SQR_PX)))

CSSS      QST_NORM = DDOT(3*NATOMS, WORK(IPRDUCT), 1,
CSSS     &                          WORK(IPRDUCT), 1)
CSSS      CALL DSCAL(3*NATOMS, DSQRT(QST_NORM), WORK(IPRDUCT), 1)

       RETURN
       END

