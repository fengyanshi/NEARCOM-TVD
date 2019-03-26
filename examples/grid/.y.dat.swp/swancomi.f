!
!     SWAN/COMPU   file 6 of 6
!
!     PROGRAM SWANCOMI.FOR
!
!     This file SWANCOMI of the main program SWAN
!     include the subroutines for solving the band matrix.
!
!     Subroutines in this file are:
!
!     CGSTAB   Solve an unsymmetric system of linear equations
!              by the Bi-CGSTAB method. The subroutine contains
!              a number of preconditioners.
!     DAXPY    ?
!     DCOPY    ?
!     DIAG     Makes a diagonal scaling of the matrix in case of
!              a momentum equations, a transport equation, or a
!              pressure equation.
!     DIAGMU   Multiplication of x with the diagonal matrix given in
!              prec. The array prec sould be filled by subroutine
!              diag.f
!     DINVL3   Multiplication of x by L, the preconditioning matrix
!              given in prec.
!     DINVU3   Multiplication of x by U, the preconditioning matrix
!              given in prec.
!     DMLU3    Calculates an upper triangular matrix U and a lower
!              triangular matrix L, which form an incomplete
!              decomposition of A.
!     DRUMA1   determine the right Hand vector b
!     ISSOLV   The subroutine issolv is used to solve an unsymmetric
!              system of equations of the shape A x = f.
!     MKPREC   The subroutine mkprec is used to build a preconditioner.
!     PREVC    Prevc multiplies the vector x with a preconditioner.
!     PRIRES   This is an output subroutine. It prints the norm of
!              the residual
!     SWCOVA2D Compute covariant base vectors in integration points       31.03
!              two-dimensional case                                       31.03
!     SWDISDT2 Distribute diffusion term for tranport equation in R2      31.03
!     SWESSBC  Puts essential boundary conditions in matrix               31.03
!     SWJCTA2D Compute sqrt(g) x contra variant base vectors in           31.03
!              integration point two-dimensional case                     31.03
!     SWTRAD2D Compute contribution of diffusion term in R2 for a         31.03
!              transport equation per integration point                   31.03
!              Compute righthandside                                      31.03
!     SWSOLV   Prepare for ISSOLV                                         31.03
!     VULMAT   Academic test for solver
!     VULMT1   Fills matrix with coefficents of SWAN
!
!
!*******************************************************************
!
      subroutine CGSTAB(n,amat,rhsd,usol,eps1,eps2,itmax,
     &               res,p,rbar,t,s,v,work,icontr,
     &               infmat ,prec, nprec, ndim, nconct,
     &               upperi, loperi, NSTATC, ITSW, ITERSW)                30.72
!
!********************************************************************
!
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!      1.0 : Kees Vuik
!     30.72: IJsbrand Haagsma
!     30.82: IJsbrand Haagsma
!     40.02: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.14: Annette Kieftenburg
!
!  1. Updates
!
!      1.0 , Mar. 94: New subroutine
!     30.72, Feb. 98: Forced at least one iteration of the CGSTAB solver
!     30.82, Sep. 98: Work with double precision to avoid underflows
!     30.82, Jul. 99: Set OMEGA to 1.0 in case of division by 0.
!     40.02, Sep. 00: Modified function DDOT
!     40.14, Mar. 01: Code adjusted to avoid no convergence problem
!                   : (by Kees Vuik)
!     40.13, Nov. 01: argument INFMAT added to call PRIRES
!
!  4. Argument variables
!
!     ITERSW: input  Iteration counter for SWAN
!     ITSW  : input  Time step counter for SWAN
!     NSTATC: input  Indicates stationarity:
!                    =0; stationary computation
!                    =1; nonstationary computation
!
      INTEGER ITERSW, ITSW, NSTATC
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     Solve an unsymmetric system of linear equations by the Bi-CGSTAB
!     method. The subroutine contains a number of preconditioners.
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     preconditioning
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!
!     amat    i   Matrix from the equations to be solved.
!
!     eps1    i   Determine the accuracy of the final approximation. The
!     eps2    i   termination criterion is based upon
!                    ||b-Ax || <  eps1 + eps2 * (b-Ax ) .
!                          k                         0
!
!     icontr i/o  Integer array in which information about the solution
!                 process must be given by the user. On output information
!                 about the solution process is provided. ICONTR must be
!                 filled as follows:
!
!             1 i Parameter controlling the preconditioning. Possible
!                 values are:
!                  0: No preconditioning.
!                  1: Diagonal preconditioner.
!                 -1: Diagonal postconditioner.
!                  2: ILUD preconditioner, with Eisenstat implemen-
!                     tation.
!                 -2: ILUD postconditioner.
!                  3: ILU preconditioner.
!                 -3: ILU postconditioner.
!
!             2 i Control parameter indicating the amount of output
!                 required. Possible values:
!                 <0: No output.
!                  1: Only fatal errors will be printed.
!                  2: Additional information about the iteration is printed.
!                  3: Gives a maximal amount of output concerning the
!                     iteration process.
!
!             3 o At output an error indication is stored. Possible values:
!                 0: No fatal errors.
!                 3: The number of iterations exceeds itmax.
!                 4: Rnorm is less than roundoff error
!
!             4 o Contains at output the actual number of iterations
!                 performed.
!
!             5 i Parameter indicating the file number to which output
!                  should be written.
!
!             6 i Parameter controlling the size of the Krylov-subspaces,
!                 thus determining the restarts of the iteration process.
!                 Only used by GMRES like methods.
!                 Possible values are:
!                 0: Restarts after at most 200 iterations are determined
!                    by the subroutine itself.
!                 1..200: The iterations restart after ICONTR(6) iterations,
!                    with the latest iterate as new starting vector.
!
!
!     infmat  i   Integer array with information of the matrix structure,
!                 to be used in matrix-vector multiplication subroutine.
!
!     itmax   i   The maximum number of iterations to be performed.
!
!     n       i   The number of rows in the matrix A.
!
!     nconct  i   Maximal number of connections in one row of the matrix.
!
!     ndim    i   Integer indicating the amount of unknowns in every grid-
!                 point. In the momentum equations ndim = 2 or 3, whereas
!                 in the pressure and transport equations ndim = 1.
!
!     nprec   i   Number of diagonals which are used in the precond-
!                 itioning.
!
!     p           Work array to store the direction vector.
!
!     prec    i   Array which contains part of the preconditioning matrix
!
!     rbar        Work array to store the quasi-residual vector.
!
!     res     o   Array containing the residual vector.
!
!     rhsd    i   Vector containing the right-hand side vector of the
!                 system of equations.
!
!     s           Work array to store an auxiliary vector.
!
!     t           Work array to store an auxiliary vector.
!
!     usol   i/o  Solution vector of length n. On input the array should
!                 contain a starting vector. At output the array contains
!                 the last iterate, which is an approximation to the
!                 solution of the system.
!
!     v           Work array to store an auxiliary vector.
!
!     work        Work array to store an auxiliary vector. The array
!                 work(.,3) contains the update of the solution usol
!                 during an iteration. If postconditioning is used,
!                 it is first adapted before it is added to usol.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     alpha    Factor alpha in the CGSTAB-process { (rbar,r ) / (rbar,Ap) }.
!                                                          i
!     beta     Factor beta in the CGSTAB-process { ((rbar,r ) / (rbar,r   ))*
!              (alpha/omega) }.
!
!     eps      Required accuracy (including relative and absolute error).
!
!     isqrn    The largest integer smaller than sqrt(n).
!
!     istar    Starting number for iteration
!
!     iter     The current iteration number.
!
!     j        Counting variable.
!                                                       i
!     omega    Factor omega in the CGSTAB-process { (t,s) / (t,t) }.
!
!     rferr    Maximal accuracy possible for this matrix, due to round-off.
!
!     rho      Inner product of M  r and rbar at the previous iteration step
!              (where M stands for the preconditioning matrix).
!
!     rhosta   Rho at the new iteration step.
!
!     rnorm    Norm of the residual vector.
!
!     sigma    Factor sigma in the CGSTAB-process (inner product of rbar
!              and Ap).
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
!     daxpy    Computes a vector + a constant * other vector. This
!              is a BLAS routine.
!
!     dcopy    Copy a vector to another vector. This is a BLAS
!              routine.
!
!     ddot     Computes an inner product. This is a BLAS routine.
!
!     dinvl    Computes the multiplication of the inverse of a lower
!              triangular matrix, stored in arrays prec and amat,
!              and a vector.
!
!     dinvu    Computes the multiplication of the inverse of a upper
!              triangular matrix, stored in arrays prec and amat,
!              and a vector.
!
!
!     dnrm2    Calculates the 2-norm of a given vector. This is
!              a BLAS routine.
!
!     mtvc1    Computes the product of the sparse matrix, stored in AMAT,
!              and a vector.
!
!     prevc    Computes the product of the preconditioner, stored in
!              amat and prec.
!
!     prires   Provides information about the solution process; the amount
!              of output is determined by the value of ICONTR(5).
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     The subroutine Bi-CGSTAB is used to solve an unsymmetric system
!     of linear equations of the shape A.x = f. To obtain a solution
!     either the Bi-CGSTAB method is used, or the preconditioned
!     Bi-CGSTAB method, possibly combined with the Eisenstat implementation.
!
!     Input:
!
!     The elements ICONTR(1:2), ICONTR(5), EPS1, EPS2, ITMAX must have got
!     a value, and the arrays AMAT, RHSD, USOL and in case of preconditioning
!     PREMAT must have been filled by the calling program.
!
!     Output:
!
!     The elements ICONTR(3:4) contain information about the iteration
!     process, the calculated solution vector is stored in USOL, and the
!     corresponding residual vector in RES.
!
!     References:
!
!     H.A. van der Vorst,
!     Bi-CGSTAB: A fast and smoothly converging variant of Bi-CG for the
!     solution of nonsymmetric linear systems,
!     SIAM J. Sci. Stat. Comp., 13, (1992), pp. 631-644.
!
! ==================================================================
!
!     Parameters:
!
      integer          n,itmax,icontr(9), infmat(*),
     &                 nprec, ndim, nconct
      REAL  amat(n,*),rhsd(n),usol(n),eps1,eps2,
     &                 res(n),p(n),rbar(n),t(n),s(n),v(n),
     &                 prec(n,*),work(n,*), upperi(*), loperi(*)
!
!     Local parameters
!
      integer j,isqrn,iter,istar
!
      DOUBLE PRECISION ALPHA, BETA, BNORM, DDOT, DNRM2, EPS, OMEGA        30.82
      DOUBLE PRECISION RFERR, RHO, RHOSTA, RNORM, SIGMA                   30.82
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'cgstab')
!
! **  Set some constants.
!
!      on real underflow ignore
!
      bnorm = dnrm2(n,rhsd,1)
!
      icontr(4) = 0
      rferr = DBLE(n)                                                     30.82
      rferr = sqrt ( rferr )
      isqrn  = nint( rferr )
      rferr  = 0.0
      rho = 1.0
      istar = 1
      beta = 0.0
      alpha = 1.0
      omega = 1.0
!
!     if icontr(1) equals 2 then the linear system is scaled by a diagonal
!     matrix, such that diag(L)=diag(U)=D=I.
!
      if(icontr(1).eq.2) then
         do j = 1, n
            rhsd(j) = rhsd(j)*prec(j,1)
         enddo
      endif
      if(abs(icontr(1)).ge.2) then
           rnorm = dnrm2(n,rhsd,1)*1D-34                                  30.82
           if(rnorm.gt.0.) then
             do j = 1, n
                usol(j) = REAL(DBLE(usol(j)) + rnorm)                     30.82
             enddo
           endif
      endif
!
!     the update solution is set equal to zero
!
      do j = 1, n
         work(j,3) = 0.0
      enddo
!
! **  Compute initial residual vector, quasi-residual vector
! **  and required accuracy.
!
100   continue
      call dcopy(n,work(1,3),1,work(1,2),1)
      if(icontr(1).lt.0) then
         call prevc(n,work(1,2),work(1,3),amat,ndim,
     &              nconct,prec,nprec,infmat,icontr(1))
      endif
!
!     solution usol is adapted
!
      call daxpy(n,1.0, work(1,3),1,usol,1)
!
!     the update solution is set equal to zero
!
      do j = 1, n
         work(j,3) = 0.0
      enddo
!
!     calculation of residual
!
      call druma1( usol, t, amat, n, nconct, infmat,
     &             upperi, loperi)
      do j = 1,n
         work(j,1) = rhsd(j) - t(j)
      enddo
      if(icontr(1).eq.2) then
!         call dinvl(work,res,amat,n/ndim,ndim,
!     1              nconct,prec,nprec,infmat)
      else
         if(icontr(1).gt.0) then
            call prevc(n,work,res,amat,ndim,nconct,
     &                 prec,nprec,infmat,icontr(1))
         else
            call dcopy(n,work,1,res,1)
         endif
      endif
!
!     calculation of rbar
!
      call dcopy(n,res,1,rbar,1)
      rnorm = dnrm2(n,res,1)
      if (istar.eq.1) then
         eps = DBLE(eps1) + rnorm * DBLE(eps2) + bnorm * DBLE(eps2)       30.82
      endif
!
! Make sure that at least one iteration is performed unless rnorm = 0     30.72
!
      if(rnorm .le. eps) then
        IF (NSTATC.GT.0) THEN                                             30.72
!
! Dynamic mode:                                                           30.72
!
! If in the first time-step and first iteration and first cgstab-         30.72
! iteration, then skip the goto 500 statement and perform an iteration    30.72
!
          IF ((ITSW.GT.1).OR.(ITERSW.GT.1).OR.(ICONTR(4).GT.0)) THEN      30.72
            icontr(3)=0
            goto 500
          END IF                                                          30.72
        ELSE                                                              30.72
!
! Stationary mode:                                                        30.72
!
! If in the first iteration and first cgstab-iteration, then skip the     30.72
! goto 500 statement and perform an iteration                             30.72
!
          IF ((ITERSW.GT.1).OR.(ICONTR(4).GT.0)) THEN                     30.72
            icontr(3)=0
            goto 500
          END IF                                                          30.72
        END IF                                                            30.72
        IF (REAL(RNORM).EQ.0.) GOTO 500                                   30.72
      endif
!
      IF ( TESTFL .AND. ITEST .GT. 75 ) THEN
        WRITE(PRTEST,*)
        WRITE(PRTEST,*)
        WRITE(PRTEST,*) ' values before entering iteration loop'
        WRITE(PRTEST,7050) ITER, ISTAR, ITMAX
 7050   FORMAT(' ++1 : ITER ISTAR ITMAX     :',3I5)
        WRITE(PRTEST,7060) RNORM, RFERR                                   30.82
 7060   FORMAT(' ++1 : RNORM RFERR          :',2D12.5)                    30.82
        WRITE(PRTEST,7070) EPS, EPS1, EPS2                                30.82
 7070   FORMAT(' ++1 : EPS, EPS1, EPS2      :',1D12.5, 2E12.5)            30.82
       ENDIF
!
! **  Iteration loop.
!
      do iter = istar,itmax
!
! **     Check for convergence.
!
         if (rnorm .le. eps) then
!           IF ( ITEST .GT. 25 ) WRITE(PRTEST,*) ' rnorm le eps'           30.50
!
! Make sure that at least one iteration is performed                      30.72
!
           IF (NSTATC.GT.0) THEN                                          30.72
!
! Dynamic mode:                                                           30.72
!
! If in the first time-step and first iteration and first cgstab-         30.72
! iteration, then skip the goto 100 statement and perform an iteration    30.72
!
             IF ((ITSW.GT.1).OR.(ITERSW.GT.1).OR.(ICONTR(4).GT.0)) THEN   30.72
               istar = iter
               goto 100
             END IF                                                       30.72
           ELSE                                                           30.72
!
! Stationary mode:                                                        30.72
!
! If in the first iteration and first cgstab-iteration, then skip the     30.72
! goto 100 statement and perform an iteration                             30.72
!
             IF ((ITERSW.GT.1).OR.(ICONTR(4).GT.0)) THEN                  30.72
               istar = iter
               goto 100
             END IF                                                       30.72
           ENDIF                                                          30.72
         end if
!
! **     Check for maximal accuracy possible due to round-off errors.
! **     If rho is smaller, then continuation is useless.
!
         if (rnorm .lt. rferr) then
            icontr(3)=4
            IF (TESTFL .AND. ITEST .GT. 75) WRITE(PRTEST,7075)
     &      ' ++2 : rnorm=', RNORM, ' < rferr=', RFERR                    30.50
7075        FORMAT(2D12.5)                                                30.82
            goto 400                                                      40.14
         end if
!
! **     Print intermediate results and set iteration level.
!
         call prires(' ISCGSTAB ',rnorm,icontr,.false.,INFMAT)            40.13
         icontr(4) = iter
!
! **     Compute rhosta
!
         call dcopy(n,res,1,s,1)
         RHOSTA = DDOT(RBAR, S, N)                                        40.02
!
         if ( abs(rhosta).lt. 1.D-30*abs(rho) .and.                       30.82
     &        iter .gt. 1 ) then
!         *  rhostar too small, restart process
            istar = iter
!
!           CALL MSGERR(3,' ERROR: rhostar too small ')          removed  30.50
            IF ( ITEST .GT. 75 .AND. TESTFL) THEN                         30.82
             WRITE(PRTEST,*) ' rhostar too small'                         30.82
             WRITE(PRTEST,7080) RHOSTA, RHO, BETA, ISTAR                  30.82
7080         FORMAT(' ++3 : RHOSTA RHO BETA ISTAR:',3D12.5,I4)            30.82
            ENDIF                                                         30.82
!
            go to 100
         end if
!
! **     Compute the direction vector p.
!
         if (iter .eq. 1) then
            call dcopy(n,s,1,p,1)
         else
            beta = (rhosta / rho) * (alpha/omega)
            do j = 1,n
               p(j) = REAL(DBLE(res(j)) +
     &                     beta * (DBLE(p(j)) - omega * DBLE(v(j))))      30.82
            enddo
         endif
!
         rho = rhosta
!
        IF ( TESTFL .AND. ITEST .GT. 75 ) THEN                            30.82
          WRITE(PRTEST,7090) RHOSTA, RHO, BETA, ITER                      30.82
7090      FORMAT(' ++4 : RHOSTA RHO BETA ITER:',3D12.5,2I4)               30.82
        ENDIF                                                             30.82
!
! **     Compute (preconditioned) Ap, store in v.
!
         if(icontr(1).gt.0) then
           call druma1( p, work, amat, n, nconct, infmat,
     &             upperi, loperi)
           call prevc(n,work,v,amat,ndim,nconct,
     &                prec,nprec,infmat,icontr(1))
         else
           call prevc(n,p,work,amat,ndim,nconct,
     &                prec,nprec,infmat,icontr(1))
           call druma1( work, v, amat, n, nconct, infmat,
     &             upperi, loperi)
         endif
!
! **     Compute sigma and terminate if less than zero.
!
         SIGMA = DDOT(RBAR, V, N)                                         40.02
         if (sigma .eq. 0.0) then
            icontr(3)=2
           IF ( ITEST .GT. 75 .AND. TESTFL) THEN                          30.82
             WRITE(PRTEST,7190) SIGMA, ICONTR(3)                          30.82
7190         FORMAT(' ++5 : SIGMA ICONTR(3)      :',1D12.5,I4)            30.82
           ENDIF                                                          30.82
            goto 400                                                      40.14
         end if
!
! **     Update the vectors s.
!
         alpha = rho / sigma
         do j = 1,n
            s(j) = REAL(DBLE(res(j)) - alpha * DBLE(v(j)))                30.82
         enddo
!
!
! **     Compute As and store in t.
!
         if(icontr(1).gt.0) then
           call druma1( s, work, amat, n, nconct, infmat,
     &             upperi, loperi)
           call prevc(n,work,t,amat,ndim,nconct,
     &                prec,nprec,infmat,icontr(1))
         else
            call prevc(n,s,work,amat,ndim,nconct,
     &                prec,nprec,infmat,icontr(1))
            call druma1( work, t, amat, n, nconct, infmat,
     &             upperi, loperi)
         endif
!
! **     Compute omega
!
         IF (DNRM2(n,t,1)**2.EQ.0.0) THEN                                 30.82
           OMEGA = 1.0                                                    30.82
         ELSE                                                             30.82
           OMEGA = DDOT(T,S,N)/(dnrm2(n,t,1)**2)                          40.02
         END IF                                                           30.82
!
! **     Update solution and residual, and compute norm.
!
         do j = 1,n
            work(j,3)= REAL(DBLE(work(j,3)) + alpha * DBLE(p(j)) +        30.82
     &                      omega * DBLE(s(j)))                           30.82
            res(j) = REAL(DBLE(s(j)) - omega * DBLE(t(j)))                30.82
!
           IF ( TESTFL .AND. ITEST .GT. 200 ) THEN                        30.82
             WRITE(PRTEST,7191) J, RES(J), S(J), T(J)                     30.82
7191         FORMAT(' ++6 : J, RES(J), S(J), T(J):',I4,1E12.5,2D12.5)     30.82
           ENDIF                                                          30.82
!
         enddo
!
         rnorm = dnrm2(n,res,1)                                           30.82
!
         IF ( TESTFL .AND. ITEST .GT. 75 ) THEN                           30.82
           WRITE(PRTEST,7192) RNORM                                       30.82
7192       FORMAT(' ++7 : RNORM',1D12.5)                                  30.82
         ENDIF                                                            30.82
!
! **     Adapt each sqrt(n) times rferr for accuracy check.
! **     Compute actual residual and subtract the computed RES. The
! **     vector t is now an indication of the possible accuracy.
!
         if (mod(iter,isqrn) .eq. 0) then
            call dcopy(n,work(1,3),1,work(1,2),1)
            if(icontr(1).lt.0) then
               call prevc(n,work(1,2),work(1,3),amat,ndim,
     &                    nconct,prec,nprec,infmat,icontr(1))
            endif
            call daxpy(n,1.0,work(1,3),1,usol,1)
            do j = 1, n
               work(j,3) = 0.0
            enddo
             call druma1( usol, t, amat, n, nconct, infmat,
     &             upperi, loperi)
            do j = 1,n
               work(j,1) = t(j) - rhsd(j)
            enddo
            if(icontr(1).eq.2) then
!               call dinvl(work,t,amat,n/ndim,ndim,
!     1                    nconct,prec,nprec,infmat)
            else
               if(icontr(1).gt.0) then
                  call prevc(n,work,t,amat,ndim,nconct,
     &                       prec,nprec,infmat,icontr(1))
               else
                  call dcopy(n,work,1,t,1)
               endif
            endif
            do j = 1,n
               t(j) = t(j) + res(j)
            enddo
            rferr = DBLE(iter)                                            30.82
            rferr = exp(0.5 * (rferr/DBLE(n)) ** 2) * dnrm2(n,t,1)        30.82
         end if
!
         IF ( TESTFL .AND. ITEST .GT. 75 ) THEN                           30.82
           WRITE(PRTEST,7181) ITER, ISTAR, ITMAX                          30.82
           WRITE(PRTEST,7180) ALPHA, OMEGA, RNORM, RFERR                  30.82
           WRITE(PRTEST,*) ' return to iter '                             30.82
7181       FORMAT(' ++8 : ITER ISTAR ITMAX        :',3I5)                 30.82
7180       FORMAT(' ++8 : ALPHA OMEGA RNORM RFERR :',4D12.5)              30.82
         ENDIF                                                            30.82
!
      enddo
!
! **  The iteration process terminates because the maximum number of
! **  iterations is reached. Set error indicator.
!
      icontr(3)=3
!                                                                         40.14
!     The solution vector is adapted when an error has occured            40.14
!                                                                         40.14
400   continue                                                            40.14
      call dcopy(n,work(1,3),1,work(1,2),1)                               40.14
      if(icontr(1).lt.0) then                                             40.14
         call prevc(n,work(1,2),work(1,3),amat,ndim,                      40.14
     &              nconct,prec,nprec,infmat,icontr(1))                   40.14
      endif                                                               40.14
!                                                                         40.14
!     solution usol is adapted                                            40.14
!                                                                         40.14
      call daxpy(n,1.0, work(1,3),1,usol,1)                               40.14
!
! **  Normal termination of the iteration process.
!
 500  continue
      call prires('ISCGSTAB ',rnorm,icontr,.true.,INFMAT)                 40.13
      call druma1( usol, t, amat, n, nconct, infmat,
     &             upperi, loperi)
      do j = 1,n
         res(j) = rhsd(j) - t(j)
      enddo
      rnorm = dnrm2(n,res,1)
!
      return
      end
!********************************************************************
!                                                                   *
      subroutine daxpy(n, da, dx, incx, dy, incy)
!                                                                   *
!********************************************************************
!
      implicit REAL  (a-h,o-z)
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  Zdenek                                     |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!      1.1 : Zdenek
!     30.82: IJsbrand Haagsma
!
!  1. Updates
!
!      1.1 , Aug. 91: no dimension, cdc, ce
!     30.82, Mar. 99: Replaced artithmic if
!
!
!     overwrite double precision dy with double precision da*dx + dy.
!     for i = 0 to n-1, replace  dy(ly+i*incy) with da*dx(lx+i*incx) +
!       dy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n,
!       and ly is defined in a similar way using incy.
!
!
      REAL  dx(*), dy(*)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'daxpy')
!
      if (n.le.0 .or. da.eq.0.0) return
      if (incx.eq.incy) THEN                                              30.82
        if (incx-1.LT.0) GOTO 10                                          30.82
        if (incx-1.EQ.0) GOTO 30                                          30.82
        if (incx-1.GT.0) GOTO 70                                          30.82
      ENDIF
   10 continue
!
!        code for nonequal or nonpositive increments.
!
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do 20 i=1,n
         dy(iy) = dy(iy) + da*dx(ix)
         ix = ix + incx
         iy = iy + incy
   20 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop so remaining vector length is a multiple of 4.
!
   30 m = mod(n,4)
      if (m.eq.0) go to 50
      do 40 i=1,m
         dy(i) = dy(i) + da*dx(i)
   40 continue
      if (n.lt.4) return
   50 mp1 = m + 1
      do 60 i=mp1,n,4
         dy(i) = dy(i) + da*dx(i)
         dy(i+1) = dy(i+1) + da*dx(i+1)
         dy(i+2) = dy(i+2) + da*dx(i+2)
         dy(i+3) = dy(i+3) + da*dx(i+3)
   60 continue
      return
!
!        code for equal, positive, nonunit increments.
!
   70 continue
      ns = n*incx
      do 80 i=1,ns,incx
         dy(i) = da*dx(i) + dy(i)
   80 continue
      return
      end
!*******************************************************************
!                                                                  *
      subroutine dcopy(n, dx, incx, dy, incy)
!                                                                  *
!*******************************************************************
!
      implicit REAL  (a-h,o-z)
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  Zdenek                                     |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!      1.1 : Zdenek
!     30.82: IJsbrand Haagsma
!
!  1. Updates
!
!      1.1 , Aug. 91: no dimension, cdc, ce
!     30.82, Mar. 99: Replaced artithmic if
!
!     copy double precision dx to double precision dy.
!     for i = 0 to n-1, copy dx(lx+i*incx) to dy(ly+i*incy),
!     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
!     defined in a similar way using incy.
!
!
      REAL  dx(*), dy(*)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'dcopy')
!
      if (n.le.0) return
      if (incx.eq.incy) THEN                                              30.82
        if (incx-1.LT.0) GOTO 10                                          30.82
        if (incx-1.EQ.0) GOTO 30                                          30.82
        if (incx-1.GT.0) GOTO 70                                          30.82
      ENDIF
   10 continue
!
!        code for unequal or nonpositive increments.
!
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do 20 i=1,n
         dy(iy) = dx(ix)
         ix = ix + incx
         iy = iy + incy
   20 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop so remaining vector length is a multiple of 7.
!
   30 m = mod(n,7)
      if (m.eq.0) go to 50
      do 40 i=1,m
         dy(i) = dx(i)
   40 continue
      if (n.lt.7) return
   50 mp1 = m + 1
      do 60 i=mp1,n,7
         dy(i) = dx(i)
         dy(i+1) = dx(i+1)
         dy(i+2) = dx(i+2)
         dy(i+3) = dx(i+3)
         dy(i+4) = dx(i+4)
         dy(i+5) = dx(i+5)
         dy(i+6) = dx(i+6)
   60 continue
      return
!
!        code for equal, positive, nonunit increments.
!
   70 continue
      ns = n*incx
      do 80 i=1,ns,incx
         dy(i) = dx(i)
   80 continue
      return
      end
!*********************************************************************
!                                                                    *
      DOUBLE PRECISION FUNCTION DDOT(DX, DY, N)                          40.02
!                                                                    *
!*********************************************************************
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.02: IJsbrand Haagsma
!
!  1. Updates
!
!     40.02, Sep. 00: New subroutine to replace old one.
!
!  2. Purpose
!
!     Calculates dot product of two vectors of equal length
!
!  3. Method
!
!     Convert vectors to double precision, multiply elements and sum
!
!  4. Modules used
!
!     --
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     DX    : First vector in dot product
!     DY    : Second vector in dot product
!     N     : Vector length
!
      INTEGER, INTENT(IN) :: N
!
      REAL, INTENT(IN)    :: DX(N), DY(N)
!
!  6. Parameter variables (Constants)
!
!     --
!
!  7. Local variables
!
!     IENT  : Number of entries into this subroutine
!
      INTEGER, SAVE       :: IENT = 0
!
!  8. Subroutines used
!
!     STRACE: Tracing routine for debugging
!
!  9. Subroutines calling
!
!     CGSTAB: Solves an unsymmetric system of linear equations by the Bi-CGSTAB method.
!
! 10. Error messages
!
!     --
!
! 11. Remarks
!
!     --
!
! 12. Structure
!
!     Make DX and DY double precision
!     Multiply DX and DY
!     Sum all elements of result
!
! 13. Source text
!
      CALL STRACE(IENT,'DDOT')
!
      DDOT = SUM(DBLE(DX) * DBLE(DY))
!
      RETURN
!
      END FUNCTION DDOT
!*******************************************************************
!                                                                  *
      subroutine diag(amat,n,ndimso,nconct,prec,nprec,infmat)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!
!     programmer       Kees Vuik
!     version 1.3      date 18-11-1992  all the element of prec are filled
!                      also in the virtual points
!     version 1.2      date 04-08-1992  Zdenek:  ndimso instead of ndim
!     version 1.1      date 07-05-1992
!     developed at     Convex, HP700
!
!
! ******************************************************************
!
!                       DESCRIPTION
!
!      Makes a diagonal scaling of the matrix in case of a momentum
!      equations, a transport equation, or a pressure equation.
!
!
! ******************************************************************
!
!                       KEYWORDS
!
!      linear_solver
!      diagonal
!      preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!      amat     i   the coefficient matrix for the momentum equations
!                   or an equation similar to the pressure equation.
!
!      infmat   i   If infmat is 1 we use the momentum equations,
!                   whereas if infmat is larger than or equal to 4
!                   we use equations with a structure similar to the
!                   pressure equation.
!
!      n        i   number of unknowns in the solution vector.
!
!      nconct   i   number of connections in one row of the matrix
!
!      ndimso   i   integer indicating the dimension of the space in
!                   which the problem must be solved (ndimso = 1 or ndim).
!
!      nprec    i   number of diagonals, which are used in the pre-
!                   conditioning. In this subroutine nprec = 1.
!
!      prec     o   the preconditioning matrix.
!
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!      i            loop counter
!
!      j            loop counter
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     To obtain a diagonal scaling of the matrix. This is done for
!     the momentum and pressure equation. The inverse of the main
!     diagonal is stored in the array prec. For transport equations
!     we use the same structure as the pressure equation.
!
! ==================================================================
!
!     input/output parameters
!
      integer          infmat,n,ndimso,nconct,nprec
      REAL  amat(1:n,1:ndimso,1:nconct),
     &                 prec(1:n,1:ndimso,1:nprec)
!
!
!     local parameters
      integer i,j
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'diag')
!
      do j = 1, ndimso
         do i = 1, n
            if(amat(i,j,1).eq.0.0) then
               prec(i,j,1) = 0.0
            else
               prec(i,j,1) = 1.0/amat(i,j,1)
            endif
         enddo
      enddo
      return
      end
!*******************************************************************
!                                                                  *
      subroutine diagmu(n,x,b,prec,nprec)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!
!     programmer       Kees Vuik
!     version 1.0      date 14-05-1992
!     developed at     Convex, HP700
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     Multiplication of x with the diagonal matrix given in prec.
!
!     The array prec sould be filled by subroutine diag.f
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     diagonal
!     preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!     b      o  the resulting vector after multiplication.
!
!     n      i  number of unknowns in the solution vector.
!
!     nprec  i  number of diagonals, which are used in the pre-
!               conditioning. In this subroutine nprec = 1.
!
!     prec   i  the diagonal preconditioning matrix.
!
!     x      i  the original vector.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     i      loop counter.
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     Multiplication of x with the diagonal matrix given in prec.
!     The result is given in b.
!
! ==================================================================
!
!     Parameters
!
      integer          n,nprec
      REAL  prec(1:n,1:nprec),x(n),b(n)
!
!     local parameters
!
      integer i
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'diagmu')
!
      do i = 1, n
         b(i) = x(i)*prec(i,1)
      enddo
      return
      end
!*******************************************************************
!                                                                  *
      subroutine dinvl3(x,b,matrix,n,ndim,nconct,
     &                  prec,nprec,infmat)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!
!     programmer       Kees Vuik
!     version 1.2      date 01-06-1992
!     developed at     Convex, HP700
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     Multiplication of x by L, the preconditioning matrix given in
!     prec.
!
!     In this case we obtain:
!                 -1
!            b = L  x.
!
!     The array prec should be filled by dmlu3.f. This subroutine
!     contains compiler directives to run in vector speed on the
!     Convex.
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     ilu
!     preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!
!     b       o  the result vector, which contains:
!                         -1
!                    b = L  x.
!
!     infmat  i  if infmat(1) is 1 we use the momentum equations,
!                whereas if infmat(1) is larger than or equal to 4
!                we use equations, with a structure similar to the
!                pressure equations.
!
!     matrix  i  the coefficient matrix for the momentum or
!                an equation similar to the pressure equation.
!
!     n       i  number of unknowns in the solution vector.
!
!     nconct  i  number of connections in one row of the matrix
!
!     ndim    i  integer indicating the dimension of the space in
!                which the problem must be solved (ndim = 2 or 3).
!
!     nprec   i  number of diagonals, which are used in the pre-
!                conditioning. In this subroutine nprec = nconct.
!
!     prec    i  the preconditioning matrix.
!
!     x       i  the original vector.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     i       loop counter.
!
!     ii      loop counter.
!
!     j       loop counter.
!
!     jjmax   loop bound.
!
!     jjmin   loop bound.
!
!     nx      number of points in the x-direction.
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     Multiplication of x with L preconditioning matrix given
!     in matrix and prec. In this case we obtain:
!                 -1
!            b = L  x.
!
!     The arrays prec and matrix should be filled by subroutine dmlu3.f
!
!     **************************************************************
!
!     Change1 15-02-1991 C. Vuik.
!
!     The reason of this change is to solve the lower triangular system
!     in vector speed.
!
!     Literature : High performance preconditioning.
!                  H. A. van der Vorst.
!                  SIAM. J. Sci. Stat. Comput. 10 pp.1174-1185 (1989).
!
!     Implementation : Compare subroutine precd in the program ictest
!                      written by H. A. van der Vorst.
!
!     Remark     : In this version the virtual points along the bound-
!                  ary are made zero. Thereafter the loop only runs
!                  over the non-boundary points.
!
!     **************************************************************
!
!     Change2 01-06-1992 C. Vuik
!
!     The virtual points are removed.
!
!     **************************************************************
!
!     Remark
!
!     - at this moment only infmat(1).ge.4 is allowed.
!
!     **************************************************************
!
! ==================================================================
!
!     Parameters
!
      integer          infmat(*),n,ndim,nconct,nprec
      REAL  matrix(1:n,1:ndim,1:nconct),
     &                 prec(1:n,1:ndim,1:nprec),
     &                 x(n,ndim),b(n,ndim)
!
!     local parameters
!
!      integer i,j,ii,nx,jjmin,jjmax
      integer nx
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'dinvl3')
!
      nx = infmat(2)+1
!                             -1
!     The multiplication b = L  x
!
!     First we fill the points on the lower horizontal boundary.
!
      b(1,1) = x(1,1)
      do i = 2, nx
         b(i,1) = x(i,1)-prec(i,1,5)*b(i-1,1)
      enddo
      i = nx+1
      b(i,1) = x(i,1)-prec(i,1,4)*b(2,1)
     &               -prec(i,1,3)*b(1,1)
      i = nx+2
      b(i,1) = x(i,1)-prec(i,1,5)*b(i-1,1)
     &               -prec(i,1,4)*b(3,1)
     &               -prec(i,1,3)*b(2,1)
     &               -prec(i,1,2)*b(1,1)
!
!     change1
!
       do i = nx+3, n
            b(i,1) = x(i,1)-prec(i,1,5)*b(i-1   ,1)
     &                  -prec(i,1,4)*b(i-nx+1,1)
     &                  -prec(i,1,3)*b(i-nx  ,1)
     &                  -prec(i,1,2)*b(i-nx-1,1)
      end do
      end
!*******************************************************************
!                                                                  *
      subroutine dinvu3(x,b,matrix,n,ndimso,nconct,
     &                  prec,nprec,infmat)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!
!     programmer       Kees Vuik
!     version 1.3      date 04-08-1992  Zdenek: ndimso instead of ndim
!     version 1.2      date 01-06-1992
!     developed at     Convex, HP700
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     Multiplication of x by U, the preconditioning matrix given in
!     prec.
!
!     In this case we obtain:
!                 -1
!            b = U  x.
!
!     The array prec should be filled by dmlu3.f. This subroutine
!     contains compiler directives to run in vector speed on the
!     Convex.
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     ilu
!     preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!     b       o  the result vector, which contains:
!                     -1
!                b = U  x.
!
!     infmat  i  if infmat(1) is 1 we use the momentum equations,
!                whereas if infmat(1) is larger than or equal to 4
!                we use equations, with a structure similar to the
!                pressure equation.
!
!     matrix  i  the coefficient matrix for the momentum or
!                an equation similar to the pressure equation.
!
!     n       i  number of unknowns in the solution vector.
!
!     nconct  i  number of connections in one row of the matrix
!
!     ndimso  i  integer indicating the dimension of the space in
!                which the problem must be solved (ndimso = 1 or ndim).
!
!     nprec   i  number of diagonals, which are used in the pre-
!                conditioning. In this subroutine nprec = nconct.
!
!     prec    i  the preconditioning matrix.
!
!     x       i  the original vector.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     i       loop counter.
!
!     ii      loop counter.
!
!     j       loop counter.
!
!     jjmax   loop bound.
!
!     jjmin   loop bound.
!
!     nx      number of points in the x-direction.
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     Multiplication of x with ILU preconditioning matrix given
!     in matrix and prec. In this case we obtain:
!                 -1
!            b = U  x.
!
!     The arrays prec and matrix should be filled by subroutine dmlu3.f
!
!     **************************************************************
!
!     Change1 15-02-1991 C. Vuik.
!
!     The reason of this change is to solve the upper triangular system
!     in vector speed.
!
!     Literature : High performance preconditioning.
!                  H. A. van der Vorst.
!                  SIAM. J. Sci. Stat. Comput. 10 pp. 1174-1185 (1989)
!
!     Implementation : Compare subroutine precd in the program ictest
!                      written by H. A. van der Vorst.
!
!     Remark     : In this version the virtual points along the boundary
!                  are made zero. Thereafter the loop only runs over
!                  the non-boundary points.
!
!     **************************************************************
!
!     Change2 01-06-1992 C. Vuik.
!
!     The virtual points are removed.
!
!     **************************************************************
!
!     Remark
!
!     - at this moment only infmat(1).ge.4 is allowed.
!
! ==================================================================
!
!     Parameters
!
      integer          infmat(*),n,ndimso,nconct,nprec
      REAL  matrix(1:n,1:ndimso,1:nconct),
     &                 prec(1:n,1:ndimso,1:nprec),
     &                 x(n,ndimso),b(n,ndimso)
!
!     local parameters
!
!      integer i,j,ii,nx,jjmin,jjmax
      integer i,nx
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'dinvu3')
!
      nx = infmat(2)+1
!                             -1
!     The multiplication b = U  x
!
!     First we fill the points on the upper horizontal boundary.
!
      b(n,ndimso) = x(n,ndimso)*prec(n,ndimso,1)
      do i = 2, nx-1
!
!       wijziging
!
!       Het statement wat hier weg gecommentarieerd staat is het foute
!       statement: prec(i,ndimso,1) moet prec(n+1-i,ndimso,1) zijn
!
!           b(n+1-i,ndimso) = (x(n+1-i,ndimso)
!       1                   -prec(n+1-i,ndimso,6)*b(n+2-i,ndimso))
!       2                   *prec(i,ndimso,1)
!
!       Hieronder het goede statement
!
        b(n+1-i,ndimso) = (x(n+1-i,ndimso)
     &                    -prec(n+1-i,ndimso,6)*b(n+2-i,ndimso))
     &                    *prec(n+1-i,ndimso,1)
!
      enddo
      i = n-nx+1
            b(i,ndimso) = (x(i,ndimso)-prec(i,ndimso,6)*b(i+1,ndimso)
     &                  -prec(i,ndimso,7)*b(i+nx-1,ndimso))
     &                  *prec(i,ndimso,1)
      i = n-nx
            b(i,ndimso) = (x(i,ndimso)-prec(i,ndimso,6)*b(i+1,ndimso)
     &                  -prec(i,ndimso,7)*b(i+nx-1,ndimso)
     &                  -prec(i,ndimso,8)*b(i+nx  ,ndimso))
     &                  *prec(i,ndimso,1)
         do i = n-nx-1, 1, -1
            b(i,ndimso) = (x(i,ndimso)-prec(i,ndimso,6)*b(i+1,ndimso)
     &                  -prec(i,ndimso,7)*b(i+nx-1,ndimso)
     &                  -prec(i,ndimso,8)*b(i+nx  ,ndimso)
     &                  -prec(i,ndimso,9)*b(i+nx+1,ndimso))
     &                  *prec(i,ndimso,1)
         enddo
      end
!*******************************************************************
!                                                                  *
      subroutine dmlu3(matrix,n,ndim,nconct,prec,nprec,infmat)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!
!     programmer       Kees Vuik
!     version 1.1      date 01-06-1992
!     developed at     Convex, HP700
!     31.03  Annette Kieftenburg
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     Calculates an upper triangular matrix U and a lower triangular
!     matrix L, which form an incomplete decomposition of A.
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     ilu
!     preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!     infmat  i  if infmat(1) is 1 we use the momentum equations,
!                whereas if infmat(1) is larger than or equal to 4
!                we use equations with a structure similar to the
!                pressure equation. infmat(2) is the number of
!                discretization points in the x-direction.
!
!     matrix  i  the coefficient matrix for the momentum equations
!                or an equation similar to the pressure equation.
!
!     n       i  number of unknowns in the solution vector.
!
!     nconct  i  number of connections in one row of the matrix
!
!     ndim    i  integer indicating the dimension of the space in
!                which the problem must be solved (ndim = 2 or 3).
!
!     nprec   i  number of diagonals, which are used in the pre-
!                conditioning. In this subroutine nprec = nconct.
!
!     prec    o  the preconditioning matrix.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     a       this parameter is used to make a combination of ILU
!             and MILU preconditioning.
!
!     i       loop counter.
!
!     n1      n1 is needed to describe the matrix elements.
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     This subroutine calculates an upper triangular matrix U and a lower
!     triangular matrix L, These matrices form an incomplete decomposition
!     for the matrix A using the following rules:
!
!            A = L U - R,
!
!     (a) diag(L) = I,
!
!     (b) The nonzero pattern of L and U are equal to the nonzero pattern
!         of A,
!
!     (c) If a(i,j) # 0 then L*U(i,j) = a(i,j)
!
!     The off-diagonal elements of U are stored in prec(1:n,1:ndim,6:9),
!     whereas inverse(diag(U)) is stored in prec(1:n,1:ndim,1).
!     The off-diagonal elements of L are stored in prec(1:n,1:ndim,2:5),
!
!     **************************************************************
!
!     Remarks
!
!     There are two difficulties to obtain the incomplete decomposition
!
!     - firstly, since the index of the matrix element can be less than
!       1 we discriminate 9 different situations:(n1 = ni+3)
!              i = 1,
!              i = 2,...,n1-1
!              i = n1,
!              i = n1+1,
!              i = n1+2,...,n-n1-2,
!              i = n-n1-1
!              i = n-n1
!              i = n-n1+1,...,n-1
!              i = n
!
!     - secondly, it is known that in the virtual cells the corresponding
!       row and collumn in the matrix consists of zero elements. So it is
!       possible that d(i) = 0. To circumvent this possibility we imple-
!       ment the following check:
!          if(a(i,1).eq.0) then prec(i,1) = 1 else ...
!
! ==================================================================
!
!     Parameters
!
      integer          infmat(*),n,ndim,nconct,nprec
      REAL  matrix(1:n,1:ndim,1:nconct),
     &                 prec(1:n,1:ndim,1:nprec)
!
!     local parameters
!
      integer          i,n1
      REAL  a
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'dmlu3')
!
!     The parameter a contains the value of alpha used in the RILU(alpha)
!     preconditioner. It appears that a good choice for the poisson
!     equation (infmat(1) = 5) is:
!        alpha = 0     for LSETUP=1, and
!        alpha = 0.975 for LSETUP=2.
!
      a = 0.975
      if (infmat(1).EQ.5) then                                            31.03
         LSETUP = infmat(8)                                               31.03
         if ( LSETUP.EQ.1 ) then                                          31.03
            a = 0.0e0
         else
            a = 0.975e0
         end if
      end if
      n1 = infmat(2)+1
!
!     Compute the incomplete decomposition such that the part
!     prec(i,1,*) is filled
!
!     i = 1
!
      i = 1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
         prec(i,1,7) = matrix(i,1,7)
         prec(i,1,8) = matrix(i,1,8)
         prec(i,1,9) = matrix(i,1,9)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*matrix(i+1   ,1,5)
         prec(i+n1-1,1,4) = prec(i,1,1)*matrix(i+n1-1,1,4)
         prec(i+n1  ,1,3) = prec(i,1,1)*matrix(i+n1  ,1,3)
         prec(i+n1+1,1,2) = prec(i,1,1)*matrix(i+n1+1,1,2)
!
!     i = 2, n1-1
!
      do i = 2, n1-1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
         prec(i,1,7) = matrix(i,1,7)-prec(i,1,5)*prec(i-1,1,8)
         prec(i,1,8) = matrix(i,1,8)-prec(i,1,5)*prec(i-1,1,9)
         prec(i,1,9) = matrix(i,1,9)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*matrix(i+1   ,1,5)
         prec(i+n1-1,1,4) = prec(i,1,1)*(matrix(i+n1-1,1,4)
     &                      -prec(i+n1-1,1,3)*prec(i-1,1,6))
         prec(i+n1  ,1,3) = prec(i,1,1)*(matrix(i+n1  ,1,3)
     &                      -prec(i+n1  ,1,2)*prec(i-1,1,6))
         prec(i+n1+1,1,2) = prec(i,1,1)*matrix(i+n1+1,1,2)
      end do
!
!     i = n1
!
      i = n1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
     &                 -prec(i,1,4)*prec(i-n1+1,1,8)
         prec(i,1,7) = matrix(i,1,7)-prec(i,1,5)*prec(i-1,1,8)
         prec(i,1,8) = matrix(i,1,8)-prec(i,1,5)*prec(i-1,1,9)
         prec(i,1,9) = matrix(i,1,9)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*(matrix(i+1   ,1,5)
     &                      -prec(i+1,1,3)*prec(i-n1+1,1,7))
         prec(i+n1-1,1,4) = prec(i,1,1)*(matrix(i+n1-1,1,4)
     &                      -prec(i+n1-1,1,3)*prec(i-1,1,6))
         prec(i+n1  ,1,3) = prec(i,1,1)*(matrix(i+n1  ,1,3)
     &                      -prec(i+n1  ,1,2)*prec(i-1,1,6))
         prec(i+n1+1,1,2) = prec(i,1,1)*matrix(i+n1+1,1,2)
!
!     i = n1+1
!
      i = n1+1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,3)*prec(i-n1  ,1,8)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
     &                 -prec(i,1,3)*prec(i-n1  ,1,9)
     &                 -prec(i,1,4)*prec(i-n1+1,1,8)
         prec(i,1,7) = matrix(i,1,7)-prec(i,1,5)*prec(i-1,1,8)
         prec(i,1,8) = matrix(i,1,8)-prec(i,1,5)*prec(i-1,1,9)
         prec(i,1,9) = matrix(i,1,9)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*(matrix(i+1   ,1,5)
     &                      -prec(i+1,1,2)*prec(i-n1  ,1,8)
     &                      -prec(i+1,1,3)*prec(i-n1+1,1,7))
         prec(i+n1-1,1,4) = prec(i,1,1)*(matrix(i+n1-1,1,4)
     &                      -prec(i+n1-1,1,3)*prec(i-1,1,6))
         prec(i+n1  ,1,3) = prec(i,1,1)*(matrix(i+n1  ,1,3)
     &                      -prec(i+n1  ,1,2)*prec(i-1,1,6))
         prec(i+n1+1,1,2) = prec(i,1,1)*matrix(i+n1+1,1,2)
!
!     i = n1+2, n-n1-1
!
      do i = n1+2, n-n1-1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,2)*(prec(i-n1-1,1,9)
     &                    +a*prec(i-n1-1,1,7))
     &                    -prec(i,1,3)*prec(i-n1  ,1,8)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
     &                 -prec(i,1,3)*prec(i-n1  ,1,9)
     &                 -prec(i,1,4)*prec(i-n1+1,1,8)
         prec(i,1,7) = matrix(i,1,7)-prec(i,1,5)*prec(i-1,1,8)
         prec(i,1,8) = matrix(i,1,8)-prec(i,1,5)*prec(i-1,1,9)
         prec(i,1,9) = matrix(i,1,9)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*(matrix(i+1   ,1,5)
     &                      -prec(i+1,1,2)*prec(i-n1  ,1,8)
     &                      -prec(i+1,1,3)*prec(i-n1+1,1,7))
         prec(i+n1-1,1,4) = prec(i,1,1)*(matrix(i+n1-1,1,4)
     &                      -prec(i+n1-1,1,3)*prec(i-1,1,6))
         prec(i+n1  ,1,3) = prec(i,1,1)*(matrix(i+n1  ,1,3)
     &                      -prec(i+n1  ,1,2)*prec(i-1,1,6))
         prec(i+n1+1,1,2) = prec(i,1,1)*matrix(i+n1+1,1,2)
      end do
!
!     i = n-n1
!
      i = n-n1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,2)*(prec(i-n1-1,1,9)
     &                    +a*prec(i-n1-1,1,7))
     &                    -prec(i,1,3)*prec(i-n1  ,1,8)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
     &                 -prec(i,1,3)*prec(i-n1  ,1,9)
     &                 -prec(i,1,4)*prec(i-n1+1,1,8)
         prec(i,1,7) = matrix(i,1,7)-prec(i,1,5)*prec(i-1,1,8)
         prec(i,1,8) = matrix(i,1,8)-prec(i,1,5)*prec(i-1,1,9)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*(matrix(i+1   ,1,5)
     &                      -prec(i+1,1,2)*prec(i-n1  ,1,8)
     &                      -prec(i+1,1,3)*prec(i-n1+1,1,7))
         prec(i+n1-1,1,4) = prec(i,1,1)*(matrix(i+n1-1,1,4)
     &                      -prec(i+n1-1,1,3)*prec(i-1,1,6))
         prec(i+n1  ,1,3) = prec(i,1,1)*(matrix(i+n1  ,1,3)
     &                      -prec(i+n1  ,1,2)*prec(i-1,1,6))
!
!
!     i = n-n1+1
!
      i = n-n1+1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,2)*(prec(i-n1-1,1,9)
     &                    +a*prec(i-n1-1,1,7))
     &                    -prec(i,1,3)*prec(i-n1  ,1,8)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
         prec(i,1,6) = matrix(i,1,6)
     &                 -prec(i,1,3)*prec(i-n1  ,1,9)
     &                 -prec(i,1,4)*prec(i-n1+1,1,8)
         prec(i,1,7) = matrix(i,1,7)-prec(i,1,5)*prec(i-1,1,8)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*(matrix(i+1   ,1,5)
     &                      -prec(i+1,1,2)*prec(i-n1  ,1,8)
     &                      -prec(i+1,1,3)*prec(i-n1+1,1,7))
         prec(i+n1-1,1,4) = prec(i,1,1)*(matrix(i+n1-1,1,4)
     &                      -prec(i+n1-1,1,3)*prec(i-1,1,6))
!
!
!     i = n-n1+2, n-1
!
      do i = n-n1+2, n-1
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,2)*(prec(i-n1-1,1,9)
     &                    +a*prec(i-n1-1,1,7))
     &                    -prec(i,1,3)*prec(i-n1  ,1,8)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
            prec(i,1,6) = matrix(i,1,6)
     &                 -prec(i,1,3)*prec(i-n1  ,1,9)
     &                 -prec(i,1,4)*prec(i-n1+1,1,8)
!
!     the matrix L
!
         prec(i+1   ,1,5) = prec(i,1,1)*(matrix(i+1   ,1,5)
     &                      -prec(i+1,1,2)*prec(i-n1  ,1,8)
     &                      -prec(i+1,1,3)*prec(i-n1+1,1,7))
      end do
!
!
!     i = n
!
      i = n
!
!     the matrix U
!
         if(matrix(i,1,1).eq.0) then
            prec(i,1,1) = 1.0
         else
            prec(i,1,1) = matrix(i,1,1)
     &                    -prec(i,1,2)*(prec(i-n1-1,1,9)
     &                    +a*prec(i-n1-1,1,7))
     &                    -prec(i,1,3)*prec(i-n1  ,1,8)
     &                    -prec(i,1,4)*(prec(i-n1+1,1,7)
     &                    +a*prec(i-n1+1,1,6)+a*prec(i-n1+1,1,9))
     &                    -prec(i,1,5)*(prec(i-1   ,1,6)
     &                    +a*prec(i-1   ,1,7))
            prec(i,1,1) = 1.0/prec(i,1,1)
         end if
      return
      end
!***********************************************************************
!                                                                      *
      DOUBLE PRECISION FUNCTION DNRM2(N, DX, INCX)                        30.82
!                                                                      *
!***********************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  Zdenek                                     |
!   --|-----------------------------------------------------------|--
!
!
!  0. Authors
!
!      1.0 : C.L. Lawson
!      1.1 : Zdenek
!     30.82: IJsbrand Haagsma
!
!  1. Updates
!
!      1.0 , Jan. 78: New FUNCTION
!      1.1 , Aug. 91: No dimension: CDC, CE
!     30.82, Sep. 98: To avoid errors using the Cray-cf90 compiler
!                     deleted obsolescent fortran90 (assigned goto)
!     30.82, Sep. 98: Works with double precision to avoid underflows
!
!  2. Purpose
!
!     Calculates the Euclidean norm of a vector DX() of length N
!
!  3. Method
!
!        version 1.1    date   22-08-91 (Zdenek: no dimension, cdc, ce)
!
!     euclidean norm of the n-vector stored in dx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1
!
!           c.l.lawson, 1978 jan 08
!
!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  sqrt(u/eps)  over all known machines.
!         cuthi = minimum of  sqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)
!
!     brief outline of algorithm..
!
!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.
!
!     values for cutlo and cuthi..
!     from the environmental parameters listed in the imsl converter
!     document the limiting values are as follows..
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)
!                   thus cutlo = 2**(-51) = 4.44089e-16
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!                   thus cuthi = 2**(63.5) = 1.30438e19
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!                   thus cutlo = 2**(-33.5) = 8.23181d-11
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /
!
!  4. Argument variables
!
!     N     : Length of the vector in DX()
!     INCX  : Stride of the vector stored in DX()
!
      INTEGER N, INCX
!
!     DX    : Array containing the vector
!
      REAL DX(*)
!
!  6. Local variables
!
!     CUTHI :
!     CUTLO :
!     ONE   : the double precision number 1
!     SUM   :
!     SZERO : the single precision number 0
!     XMAX  :
!     ZERO  : the double precision number 0
!
      DOUBLE PRECISION  ONE, SUM, XMAX, ZERO                              30.82
      REAL   CUTHI, CUTLO, HITEST, SZERO                                  30.82
!
      data cutlo, cuthi /8.232e-11,1.304e19/
      data zero, one /0.0,1.0/
      DATA SZERO /0.0/                                                   30.82
!
!     I
!     J
!     NEXT
!     NN
!
      INTEGER I, J, NEXT, NN
!
!     on real underflow ignore
!
      if (n.gt.0) go to 10
      dnrm2 = zero
      go to 140
!
   10 NEXT = 30                                                           30.82
      sum = zero
      nn = n*incx
!                                                 begin main loop
      i = 1
   20 IF (NEXT.EQ.30) GOTO 30                                             30.82
      IF (NEXT.EQ.40) GOTO 40                                             30.82
      IF (NEXT.EQ.70) GOTO 70                                             30.82
      IF (NEXT.EQ.80) GOTO 80                                             30.82
   30 if (abs(dx(i)).gt.cutlo) go to 110
      NEXT = 40                                                           30.82
      xmax = zero
!
!                        phase 1.  sum is zero
!
   40 if (abs(dx(i)).gt.szero) go to 41                                   30.82
      GOTO 130                                                            30.82
   41 continue                                                            30.82
      if (abs(dx(i)).gt.cutlo) go to 110
!
!                                prepare for phase 2.
      NEXT = 70                                                           30.82
      go to 60
!
!                                prepare for phase 4.
!
   50 i = j
      NEXT = 80
      sum = (sum/DBLE(dx(i)))/DBLE(dx(i))                                 30.82
   60 xmax = DBLE(abs(dx(i)))                                             30.82
      go to 90
!
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.
!
   70 if (abs(dx(i)).gt.cutlo) go to 100
!
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.
!
   80 if (abs(dx(i)).le.REAL(xmax)) go to 90                              30.82
      sum = one + sum*(xmax/DBLE(dx(i)))**2                               30.82
      xmax = DBLE(abs(dx(i)))                                             30.82
      go to 130
!
   90 sum = sum + (DBLE(dx(i))/xmax)**2                                   30.82
      go to 130
!
!
!                  prepare for phase 3.
!
  100 sum = (sum*xmax)*xmax
!
!
!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)
!
  110 hitest = cuthi/REAL(n)                                              30.82
!
!                   phase 3.  sum is mid-range.  no scaling.
!
      do 120 j=i,nn,incx
         if (abs(dx(j)).ge.hitest) go to 50
         sum = sum + DBLE(dx(j))**2                                       30.82
  120 continue
      dnrm2 = sqrt(sum)
      go to 140
!
  130 continue
      i = i + incx
      if (i.le.nn) go to 20
!
!              end of main loop.
!
!              compute square root and adjust for scaling.
!
      dnrm2 = xmax*sqrt(sum)
  140 continue
!     on real underflow abort
      return
!
!     End of function DNRM2
!
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine druma1(x,b,matrix,n,nconct,infmat,
     &                  upperi, loperi)
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!  0. Authors
!
!      1.0 : Kees Vuik
!     30.72: IJsbrand Haagsma
!     31.04: Nico Booij
!
!  1. Updates
!
!      1.0 , Jan. 78: New subroutine
!     31.04, Apr. 98: procedure not done for setup, only propagation
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
      integer          infmat(*),n,nconct
      REAL  matrix(1:n,1:nconct),
     &                 x(n),b(n), upperi(*), loperi(*)
      integer i,ni,n1
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'druma1')
!
!
      ni = infmat(2)
      n1 = ni+1
      i = 1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)+
     &          matrix(i,8)*x(i+n1)+
     &          matrix(i,9)*x(i+n1+1)
      do i = 2, n1-1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)+
     &          matrix(i,8)*x(i+n1)+
     &          matrix(i,9)*x(i+n1+1)
      end do
      i = n1
!
!                 temporary test
!
      if (i.lt.2 .or. (i+1).gt.n .or. (i+n1+1).gt.n)
     &             write (prtest,14) i, n1, n
  14  format (' error druma1 ', 6i6)
!
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)+
     &          matrix(i,8)*x(i+n1)+
     &          matrix(i,9)*x(i+n1+1)
      i = n1+1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,3)*x(i-n1)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)+
     &          matrix(i,8)*x(i+n1)+
     &          matrix(i,9)*x(i+n1+1)
      do i = n1+2, n-n1-1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,2)*x(i-n1-1)+
     &          matrix(i,3)*x(i-n1)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)+
     &          matrix(i,8)*x(i+n1)+
     &          matrix(i,9)*x(i+n1+1)
      end do
      i = n-n1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,2)*x(i-n1-1)+
     &          matrix(i,3)*x(i-n1)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)+
     &          matrix(i,8)*x(i+n1)
      i = n-n1+1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,2)*x(i-n1-1)+
     &          matrix(i,3)*x(i-n1)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)+
     &          matrix(i,7)*x(i+n1-1)
      do i = n-n1+2, n-1
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,2)*x(i-n1-1)+
     &          matrix(i,3)*x(i-n1)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)+
     &          matrix(i,6)*x(i+1)
      end do
      i = n
         b(i) = matrix(i,1)*x(i)+
     &          matrix(i,2)*x(i-n1-1)+
     &          matrix(i,3)*x(i-n1)+
     &          matrix(i,4)*x(i-n1+1)+
     &          matrix(i,5)*x(i-1)
!
!     now the periodic part is taken into account (full circle case)
!     The test if ( infmat(1).EQ.4 ) is added when the poisson equation
!     to compute the SETUP was included by Kees Vuik at 31-10-1997.
!     infmat(1) = 4 represents the original equation (wave propagation)
!     infmat(1) = 5 represents the poisson equation (setup)
!
      if ( infmat(1).EQ.4 ) then                                          31.01
         do i = 1, n1
            b(i) = b(i)+upperi(i)*x(n-n1+i)
            b(n-n1+i) = b(n-n1+i)+loperi(i)*x(i)
         end do
      end if                                                              31.01
      return
      end
!************************************************************************
!                                                                       *
      SUBROUTINE ISSOLV(IINSOL  ,RINSOL  ,MATRIX  ,RHSIDE   ,SOLUT    ,
     &                  NUSOL   ,NCONCT  ,INFMAT  ,WORK     ,NWORK    ,
     &                  PRECON  ,NPREC   ,UPPERI  ,LOPERI   ,INOCNV   ,
     &                  ITSW    ,ITERSW  )                                30.72
!                                                                       *
!************************************************************************
!
      USE M_PARALL                                                        40.31

      INCLUDE 'swcomm1.inc'                                               30.80
      INCLUDE 'swcomm3.inc'                                               30.80
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!     | developed at:  Toshiba PC, Convec, HP700                  |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Nov. 97: Moved ERRPTS to other type declarations
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.72, Feb. 98: Modified argument list for update CGSTAB solver
!     34.01, Feb. 99: Changed STOP statements in a MSGERR(4,'message')
!                     calls
!     34.01, Feb. 99: Introducing STPNOW
!     30.80, July 99: KCGRD, ICMAX, ERRPTS, IX, IY, NSTATC removed from argument list
!                     swcomm1 and swcomm3 now included;
!                     error messages modified
!     30.82, Sep. 99: Modified messages in case of non-convergence
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!
!  2. Purpose
!
!     The subroutine issolv is used to solve an unsymmetric system
!     of equations of the shape A x = f.
!
!  3. Method
!
!     Keywords:
!     linear_solver, preconditioning
!
!     At present the following solution methods are available:
!     CGSTAB and the Bi-CGSTAB method.
!
!  4. Argument variables
!
!     ITERSW: input  Iteration counter for SWAN
!     ITSW  : input  Time step counter for SWAN
!     NSTATC: input  Indicates stationarity:
!                    =0; stationary computation
!                    =1; nonstationary computation
!
      INTEGER ITERSW, ITSW
!
!     iinsol   i   reserved
!     infmat   i   reserved
!     matrix   i   reserved
!     nprec    i   reserved
!     nusol    i   reserved
!     nwork    i   reserved
!     precon  i/o  reserved
!     rhside   i   reserved
!     rinsol   i   reserved
!     solut    o   reserved
!     work     o   reserved
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     conv     Logical variable indicating acquired convergence.
!     dt       Timestep.
!     eps1     Required absolute accuracy.
!     eps2     Required reduction of the initial residual.
!     eps3     Required relative error of the solution.
!     i        Local counting variable
!     icontr   Array containing information for the called subroutines.
!     iprint   File number to which output should be written.
!     ipstrt   Used in the call of isextr
!     irrind   Variable indicating possible failures.
!     istop    Variable indicating whether the process must stop.
!     itmax    Maximal number of iteration to be performed.
!     iwarn    Variable indicating whether output must be given.
!     lwork    Extra work space needed for matrix vector multiplication
!     method   Variable determining the method to be used.
!     name     Name of calling subroutine for error messages
!     nconct   number of connections in one row of the matrix.
!     ndim     Dimension of spacec
!     ndimso   number of unknowns per one gridpoint,
!              ndimso = ndim for the momentum equations,
!              ndimso = 1    for all other equations.
!     ni       number of gridpoints in x-direction.
!     nusol2   number of scalar unkowns for the solvers:
!              nusol2 = nusol*ndim   for momentum equation
!                     = nusol        for other equations
!              nusol2 = nusol*ndimso  for any equation
!     rho      Specific mass.
!
!     TSTFLO:
!
      LOGICAL TSTFLO
!
!  7. Common Blocks used
!
!     ---
!
!  8. Subroutines used
!
!     CGSTAB Computes the solution of an unsymmetric system of linear
!            equations. This is the Bi-CGSTAB method.
!     MKPREC Builds the preconditioner and if iinsol(2).eq.2 then
!            mkprec scales the matrix too.
!
      LOGICAL STPNOW                                                      34.01
!
! 10. Error messages
!
!     128 :  Nwork is too small.
!     129 :  This preconditioning is not allowed in this subroutine.
!            (Default values and user values are filled by ISSP06)
!     130 :  No convergence occured in the linear solver.
!
! 12. Structure
!
!     We start with a call to isextr, which initializes the starting
!     vector in work(1:nusol). If we use GMRES and iinsol(2).ne.0
!     then a call of mkprec builds a preconditioner. In the pressure
!     equation this should be done only the first time. Thereafter
!     the system of equations is solved. Finally if the iterative
!     method has converged the solution is copied from work(1:nusol)
!     into solut by a call of isputs.
!
!     Input:
!     The elements IINSOL(1:8,11:14), and RINSOL(1:7) must have got a value,
!     and the arrays GMAT, MATRIX, RHSIDE, INFMAT must have been filled.
!     The parameters ndefgd,nwork, must have an approptiate
!     value.
!
!     Output:
!     The elements IINSOL(9:10) provide information about the solution
!     process. The calculated solution vector is delivered in SOLUT. Additional
!     information, which is meant for research purposes only, might be
!     contained in the array WORK. The parameter nprec and the array
!     precon are filled. If iinsol(2).eq.2 then the matrix is scaled.
!
!     The routines ISEXTR and ISPUTS must be provided by the user.
!
! 13. Source text
!
      integer iinsol(*), nusol ,
     &        infmat(*), nwork, nprec, i   , INOCNV                       30.80
!
!
      REAL    matrix(*), rhside(*),
     &        precon(*), work(nwork), rinsol(*),
     &        solut(*) , upperi(*), loperi(*)
!
!     Local parameters:
!
      integer method, irrind, iprint, itmax, istop, iwarn, nusol2,
     &        icontr(9), lwork,ndim,ndimso,nconct,ni
!     integer ipstrt
      logical conv
!      REAL  eps1, eps2, eps3,rho,dt
      REAL  eps1, eps2
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'issolv')
!
!     *  fill ndim, ndimso, nusol2
!
      if ( infmat(4) .eq. 0 ) then
        ndim = 2
      else
        ndim = 3
      endif
!
      if ( infmat(1) .eq. 1 ) then
         ndimso = ndim
      else
         ndimso = 1
      end if
!
      nusol2 = nusol * ndimso
!
!     **  Initialization of constants.
!
      do i=1, 9
         icontr(i) = 0
      enddo
      iprint    = iinsol(8)
      itmax     = iinsol(6)
      eps1      = rinsol(1)
      eps2      = rinsol(2)
!      eps3      = rinsol(5)
      icontr(1) = iinsol(2)
      icontr(2) = iinsol(3)
      icontr(5) = iprint
      istop     = iinsol(5)
      iwarn     = iinsol(3)
      ni        = infmat(2)
      irrind = 0
      conv = .false.
      method = iinsol(7)
!
!     Preconditioning is allowed at this moment, abs(iconrt(1))
!     should be <= 2 for the momentum equations,
!               <= 3 for the other equations.
!     Moreover, in 3D case abs(icontr(1)) should be
!               <= 2 for the momentum equations,
!               <= 2 for the other equations.
!
      if (infmat(1).eq.1.and.(abs(icontr(1)).gt.2)
     &    .or.abs(icontr(1)).gt.3
     &    .or.ndim.eq.3.and.abs(icontr(1)).gt.2) then
         if (iwarn .ge. 0) then
           IF (INFMAT(1).EQ.4) THEN
             write (PRTEST,4000) ixcgrd(1)-1, iycgrd(1)-1
 4000        format(
     &       ' Error: preconditioner does not exist for grid point ',
     &       2I5)
             CALL MSGERR(4,'preconditioner does not exist')               34.01
           ELSE
             CALL MSGERR(4,'preconditioner for setup does not exist')     30.80
           ENDIF
         endif
         irrind = 1
         icontr(3) = 1
         goto 999
      end if
!
!     If ni = 2 then the ILU-preconditioning applied to the
!     pressure equations can give breakdown, because the
!     pressure equation is singular and for ni = 2 ILU is
!     equal to the LU-decomposition. In such a case we make
!     icontr(1) = -2.
!
      if (infmat(1).eq.4.and.abs(icontr(1)).eq.3
     &   .and.ni.eq.2) icontr(1) = -2
      call mkprec(matrix,nusol2,ndimso,nconct,precon,nprec,
     &            infmat,icontr(1))
      IF (STPNOW()) RETURN                                                34.01
      if (method .eq. 4 .and. iinsol(13).eq.0 ) then
!
!     Solve using Bi-CGSTAB.
!
         lwork = 10 * nusol2
         if (nwork .lt. lwork) then
            if (iwarn .ge. 0) then
              IF (INFMAT(1).EQ.4) THEN
                write (PRTEST,4010) ixcgrd(1)-1, iycgrd(1)-1              30.80
 4010           format(
     &          ' Error: insufficient memory for grid point ', 2I5)       30.80
                CALL MSGERR(4,'memory is too small for CGStab solver')    30.80
              ELSE
                CALL MSGERR(4,'memory is too small for CGStab solver')    30.80
              ENDIF
            endif
            icontr(3) = 2
         else
!
            tstflo = testfl                                               30.82
            testfl = .true.                                               30.82
            call cgstab(nusol2, matrix, rhside, solut, eps1,
     &         eps2, itmax, work(nusol2+1), work(2*nusol2+1),
     &         work(3*nusol2+1), work(4*nusol2+1), work(5*nusol2+1),
     &         work(6*nusol2+1), work(7*nusol2+1), icontr, infmat,
     &         precon, nprec, ndimso, nconct,upperi, loperi,              30.72
     &         NSTATC, ITSW, ITERSW)                                      30.72
            testfl = tstflo                                               30.82
         end if
      end if
!
!     end solve using Bi-CGSTAB.
!
 999  continue
!
!     Check for normal termination.
!
      irrind = irrind + icontr(3)
      if (icontr(3) .eq. 0) conv = .true.
      iinsol(9) = irrind
!
!     Copy solution from work(1:nusol2)
!     into it position in the second half of the space solut(*),
!     in case of convergence.
!
!     Initialisation of flag CSETUP
!
      CSETUP = .TRUE.
!
      if (conv) then
         iinsol(9)  = 0
         iinsol(10) = icontr(4)
      else
         IF (INFMAT(1).EQ.4) THEN                                         30.82
!
!        No convergence in spectral plane
!
           INOCNV = INOCNV + 1                                            30.82
           IF (ERRPTS.GT.0.AND.INODE.EQ.MASTER) THEN                      40.30
              WRITE(ERRPTS,7002) IXCGRD(1)+MXF-1, IYCGRD(1)+MYF-1, 2      40.30 30.80
 7002         FORMAT (I4, 1X, I4, 1X, I2)
           END IF
           IF ((ISTOP .EQ. 0).AND.(ITEST.GE.30)) THEN                     30.82
             WRITE (PRTEST,7005) IXCGRD(1)-1, IYCGRD(1)-1                 30.80
 7005        FORMAT(' No CGStab conv. in gridpoint (', 2I5,
     &       '); -> continue process')                                    40.00
           END IF
         ELSE
!
!        No convergence for setup
!
           CSETUP = .FALSE.                                               30.82
           if ((istop .eq. 0).AND.(ITEST.GE.30)) then                     30.82
             WRITE (PRTEST,7007)                                          30.80
 7007        FORMAT(' No CGStab conv. in comp. of setup')                 30.80
           ENDIF
         ENDIF
      end if
!
      end
!*******************************************************************
!                                                                  *
      subroutine mkprec(matrix,nusol,ndimso,nconct,precon,
     &           nprec,infmat,mkind)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     34.01: IJsbrand Haagsma
!
!  1. Updates
!
!     34.01, Feb. 99: Changed STOP statement in a MSGERR(4,'message')
!                     call
!
!
!     programmer       Kees Vuik
!     version 1.2      date 24-11-1992  Kees V. dmlu works if ndim = 3.
!     version 1.1      date 04-08-1992  Zdenek:  ndimso instead of ndim
!     version 1.0      date 12-05-1992
!     developed at     Convex/HP700
!
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     The subroutine mkprec is used to build a preconditioner.
!
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!     infmat  i   Array which describes the structure of matrix.
!
!     matrix i/o  Double precision array in which the matrix of the
!                 linear system of equations is stored. In the case
!                 of mkind = 2 the matrix is scaled.
!
!     mkind   i   The kind of the preconditioner required.
!
!     nconct  i   The number of non-zero diagonals of MATRIX.
!
!     ndimso   i   The dimension of the space for the solver:
!                 ( ndimso = 1     for noncoupled equations,
!                   ndimso > 1     for coupled equations )
!
!     nprec   o   Maximum number of diagonals in PRECON.
!
!     nusol   i   The length of the solution vector.
!
!     precon  o   Double precision array in which a preconditioning
!                 matrix might be stored, of length NPREC * NUSOL. It is
!                 assumed that PRECON has a similar structure as MATRIX.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     ndim        Dimension of space.
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
!     diag        builds the diagonal preconditioner.
!
!     dmlu        builds the ILUD preconditioner and scales the
!                 matrix for ndim = 2.
!
!     dmlu2       builds the ILUD preconditioner for ndim = 2.
!
!     dmlu3       builds the ILU preconditioner.
!
!     tmlu        builds the ILUD preconditioner for ndim = 3.
!
!     tmlu2       builds the ILUD preconditioner for ndim = 3.
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
! ==================================================================
!
!     Parameters:
!
      integer          nusol,ndimso,nconct,nprec,infmat(*),mkind
      REAL  matrix(*), precon(*)
!
!     Local parameters:
!
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'mkprec')
!
      if(abs(mkind).eq.1) then
!
!        *** build the diagonal scaling preconditioning matrix
!
         nprec = 1
         call diag(matrix,nusol/ndimso,ndimso,nconct,
     &             precon,nprec,infmat(1))
!
!        *** end of the diagonal scaling
!
      endif
      if(abs(mkind).eq.2) then
        CALL MSGERR(4,'This precondioner does not exits')                 34.01
      endif
      if(abs(mkind).eq.3.and.infmat(1).ge.4) then
!
!     *** build an incomplete LU decomposition of A such that
!         the sparseness pattern of L and U is the same as A.
!
         nprec = nconct
         call dmlu3(matrix,nusol/ndimso,ndimso,nconct,
     &              precon,nprec,infmat)
!
!        *** end of the ILU3 decomposition
!
      endif
      end
!*******************************************************************
!                                                                  *
      subroutine prevc(n,x,b,matrix,ndim,nconct,precon,
     &           nprec,infmat,mkind)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!
!     programmer       Kees Vuik
!     version 1.1      date 19-05-1992
!     developed at     Convex,HP700
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     Prevc multiplies the vector x with a preconditioner.
!
! ******************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     preconditioner
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
!
!     b       o   The output vector which is the preconditioner times
!                 the vector x.
!
!     infmat  i   Array which describes the structure of the matrix.
!
!     matrix  i   Double precision array in which the matrix of the
!                 linear system of equations is stored.
!
!     mkind   i   The kind of the preconditioner required.
!
!     n       i   The length of the solution vector.
!
!     nconct  i   The number of non-zero diagonals of MATRIX.
!
!     ndim    i   The dimension of the space (ndim =2 or 3).
!
!     nprec   i   Maximum number of diagonals in PRECON.
!
!     precon  i   Double precision array in which a preconditioning
!                 matrix might be stored, of length NPREC * NUSOL. It is
!                 assumed that PRECON has a similar structure as MATRIX.
!
!     x       i   The input vector.
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
!     i       loop counter
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
!     dcopy       Copy a vector in another vector.
!
!     diagmu      Computes the multiplication of a diagonal preconditioner
!                 stored in array precon and a vector.
!
!     dinvl       Computes the multiplication of the inverse of a lower
!                 triangular matrix stored in arrays matrix and precon
!                 and a vector. Precon should be filled by dmlu.f.
!
!     dinvl2      Computes the multiplication of the inverse of a lower
!                 triangular matrix stored in arrays matrix and precon
!                 and a vector. Precon should be filled by dmlu2.f.
!
!     dinvl3      Computes the multiplication of the inverse of a lower
!                 triangular matrix stored in array precon and a vector.
!                 Precon should be filled by dmlu3.f.
!
!     dinvu       Computes the multiplication of the inverse of an upper
!                 triangular matrix stored in arrays matrix and precon
!                 and a vector. Precon should be filled by dmlu.f.
!
!     dinvu2      Computes the multiplication of the inverse of an upper
!                 triangular matrix stored in arrays matrix and precon
!                 and a vector. Precon should be filled by dmlu2.f.
!
!     dinvu3      Computes the multiplication of the inverse of an upper
!                 triangular matrix stored in array precon and a vector.
!                 Precon should be filled by dmlu3.f.
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
!     The subroutine prevc is used to multiply x with a preconditioner
!     and store the result in b. The choice of the preconditioner is
!     given by mkind.
!
! ==================================================================
!
!     Parameters:
!
      integer          n,ndim,nconct,nprec,infmat(*),mkind
!     INTEGER   I
      REAL  x(*),b(*),matrix(*), precon(*)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'prevc')
!
      if(mkind.eq.0) then
!
!     *** no preconditioner is used
!
         call dcopy(n,x,1,b,1)
      endif
      if(abs(mkind).eq.1) then
!
!     *** use a diagonal scaling preconditioning matrix
!
         call diagmu(n,x,b,precon,nprec)
      endif
      if(abs(mkind).eq.3.and.infmat(1).ge.4) then
!
!     *** use an incomplete LU decomposition of A such that
!         the sparseness pattern of L and U is the same as A.
!
         call dinvl3(x,b,matrix,n/ndim,ndim,
     &               nconct,precon,nprec,infmat)
         call dinvu3(b,b,matrix,n/ndim,ndim,
     &               nconct,precon,nprec,infmat)
      endif
      return
      end
!*******************************************************************
!                                                                  *
      subroutine prires(text,rnorm,icontr,final,INFMAT)                   40.13
!                                                                  *
!*******************************************************************
!
      INCLUDE 'swcomm3.inc'                                               40.13
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!      1.0 : Kees Vuik
!     30.82: IJsbrand Haagsma
!     40.13: Nico Booij
!
!  1. Updates
!
!     30.82, Sep. 98: Work with double precision to avoid underflows
!     30.82, Aug. 99: Explenation of error codes
!     40.13, Mar. 01: error messages reorganized
!            Nov. 01: argument INFMAT added
!
!
!
! ******************************************************************
!
!                       DESCRIPTION
!
!     This is an output subroutine. It prints the norm of the residual
!
! ******************************************************************
!
!                       KEYWORDS
!
!
! ******************************************************************
!
!                       INPUT / OUTPUT   PARAMETERS
!
! deleted by Ris (02-95):  implicit none
!
!
! ******************************************************************
!
!                       COMMON BLOCKS
!
!
! ******************************************************************
!
!                       LOCAL PARAMETERS
!
! ******************************************************************
!
!                       SUBROUTINES CALLED
!
! ******************************************************************
!
!                       I/O
!
! ******************************************************************
!
!                       ERROR MESSAGES
!
! ******************************************************************
!
!                       PSEUDO CODE
!
! ==================================================================
      character *(*) text
      integer icontr(9)
      logical final
      INTEGER, INTENT(IN) :: INFMAT(*)                                    40.13
      DOUBLE PRECISION  rnorm                                             30.82
      integer iter,iprint,ierr
!     integer  iout
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'prires')
!
      iter = icontr(4)
      iprint = icontr(2)
      ierr = icontr(3)
!      iout=icontr(5)

!
      if (iter .eq. 0) then
!         if (iprint .gt. 0) write(iout,100) text
!         if (iprint .eq. 1) write(iout,110) rnorm
!         if (iprint .ge. 2) write(iout,120)
        if (iprint .gt. 0) write(PRTEST,100) text
        if (iprint .eq. 1) write(PRTEST,110) rnorm
        if (iprint .ge. 2) write(PRTEST,120)
 100    format(//,' Output of subroutine ', a,// )
 110    format(' The 2-norm of the initial residual is:',D12.6)           30.82
 120    format(' Iteration number  residual' //)
      end if
      if (iprint .ge. 2) write(PRTEST,130) iter,rnorm
 130  format(i10,D15.2)                                                   30.82
!
      if (final) then
        if (iprint .ge. 0) then
          if (ierr .eq. 3) then                                           40.13
!           convergence error
            IF (INFMAT(1).EQ.4) THEN
!             error occurs during computation of action densities         40.13
              write(PRTEST,140) text,iter,ixcgrd(1)-1, iycgrd(1)-1        40.13
 140          format(1x,a,' No convergence in ',i5,' iterations ',        40.13
     &               'at point:',2i5)                                     40.13
            ELSE
!             error occurs during computation of setup                    40.13
              write(PRTEST,142) text,iter                                 40.13
 142          format(1x,a,' No convergence in ',i5,' iterations ',        40.13
     &               '(setup computation)')                               40.13
            ENDIF                                                         40.13
          endif                                                           40.13
          IF (IERR.EQ.4) WRITE(PRTEST,144)                                30.82
 144        FORMAT(' Rnorm is less than roundoff error')                  30.82
        end if
        if (iprint .ge. 1) write(PRTEST,150) text,iter,rnorm
 150    format(' The number of executed iterations in ',a,
     &         ' is ',i5,/ ' The 2-norm of the residual is ',D12.6)       30.82
      end if
      return
      end
!
!********************************************************************
!                                                                   *
      subroutine vulmat ( n, nconct, a, infmat, upperi, loperi )
!                                                                   *
!********************************************************************
!
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Applied Mathematics and Informatics            |
!     | P.O. Box 356,   2600 AJ  Delft, the Netherlands           |
!     |                                                           |
!     | Programmer  :  C. Vuik                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
      integer n, nconct, infmat(*), i, nx, n1
      REAL  a(n, nconct), upperi(*), loperi(*)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'VULMAT')
!
      nx = infmat(2)
      n1 = nx+1
      do i = 1, n
         a(i,1) = 4.0
      enddo
      do i = 1, n-1
         a(i+1,5) = -1.0
         a(i  ,6) = -1.0
      enddo
      do i = 1, n-n1
         a(i+n1,3) = -1.0
         a(i   ,8) = -1.0
      enddo
      do i = 1, n1
         upperi(i) = -0.5
         loperi(i) = -0.5
      end do
!
      IF ( TESTFL .AND. ITEST .GE. 120) THEN
        WRITE(PRTEST,*) ' subroutine vulmat '
        WRITE(PRTEST,123) NX, N1 , N
 123    FORMAT(' VULMAT : NX  N1  NTOT  :',3I4)
        WRITE(PRTEST,*)
     & '  NPP     (3)      (5)      (1)      (6)      (8) '
        DO IPP = 1 , N
            WRITE(PRTEST,3351) IPP, A(IPP,3),A(IPP,5), A(IPP,1),
     &                         A(IPP,6),A(IPP,8)
3351          FORMAT(I3,5E10.2)
        ENDDO
      ENDIF
      end
!***********************************************************************
!                                                                      *
      SUBROUTINE VULMT1  (NTOT    ,BAND    ,UPPERI  ,LOPERI  ,RHV     ,
     &                    IMATRA  ,IMATLA  ,IMATDA  ,IMATUA  ,IMAT5L  ,
     &                    IMAT6U  ,SECTOR  ,MDC     ,MSC     ,IDDLOW  ,
     &                    IDDTOP  ,ISSTOP  ,IDCMIN  ,IDCMAX  ,ANYBIN  ,
     &                    IDTOT   ,KCGRD   ,ICMAX                     )   30.21
!                                                                      *
!***********************************************************************
!
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
      INTEGER  IS      ,ID      ,MDC     ,MSC     ,NTOT    ,IDDLOW  ,
     &         IDDTOP  ,ISSTOP  ,IDDT    ,IDDL    ,NPP     ,IDSWAN  ,
     &         IDBAND  ,IDTOT   ,ICMAX                                    30.21
!
      INTEGER  SECTOR(MSC)      ,
     &         IDCMIN(MSC)      ,
     &         IDCMAX(MSC)      ,
     &         KCGRD(ICMAX)                                               30.21
!
      REAL     IMATRA(MDC,MSC)              ,
     &         IMATLA(MDC,MSC)              ,
     &         IMATDA(MDC,MSC)              ,
     &         IMATUA(MDC,MSC)              ,
     &         IMAT5L(MDC,MSC)              ,
     &         IMAT6U(MDC,MSC)              ,
     &         BAND(NTOT,9)                 ,
     &         RHV(NTOT)                    ,
     &         UPPERI(*)                    ,
     &         LOPERI(*)
!
      LOGICAL  ANYBIN(MDC,MSC)     ,
     &         PERIOD
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'VULMT1')
!
!     *** The first dimension of array , WORK, PRECON should be > N
!     *** The dimension of EXACT, RHV, SOLUT should be > N
!     *** The dimension of UPPERI, LOPERI should be > MSC
!
!     *** fill coefficients in diagonal matrix , RHV and SOLUT  ***
!
      DO IDDUM = IDDLOW, IDDTOP
!       *** counter for SWAN arrays ***
        IDSWAN = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!       *** counter for arrays of solver ***
        IDBAND = MOD ( IDDUM - IDDLOW , MDC ) + 1
        DO IS = 1, ISSTOP
          NPP = (IDBAND - 1) * ISSTOP + IS
          IF (NPP.LT.1 .OR. NPP.GT.NTOT .OR. IDSWAN.LT.1 .OR.
     &         IDSWAN.GT.MDC .OR. IS.GT.MSC) WRITE (PRTEST,14)
     &         NPP, NTOT, IDBAND, IDSWAN, IDDLOW, IDDTOP, IS, ISSTOP
  14      FORMAT (' error VULMT1 ', 8I6)
          BAND(NPP,1)  = IMATDA(IDSWAN,IS)
          RHV(NPP)     = IMATRA(IDSWAN,IS)
        ENDDO
      ENDDO
!
!     *** fill lower and upper diagonal ( IMAT5L and IMAT6U )       ***
!     *** containing the coefficients of propagation in freq. space ***
!
      DO IDDUM = IDDLOW, IDDTOP
        IDSWAN = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        IDBAND = MOD ( IDDUM - IDDLOW  , MDC ) + 1
        DO IS = 1, ISSTOP
          NPP = (IDBAND - 1) * ISSTOP + IS
          BAND( NPP , 5 )  = IMAT5L(IDSWAN,IS)
          BAND( NPP , 6 )  = IMAT6U(IDSWAN,IS)
        ENDDO
      ENDDO
!
!     *** fill lower and upper diagonal: IMATLA(3) and IMATDA (8) ***
!     *** containing the coefficients propagation in theta space  ***
!
      DO IS = 1, ISSTOP
        PERIOD = .FALSE.
        IF ( SECTOR(IS) .EQ. 2 ) THEN
!
!         *** check is domain in directional space is periodic ***
!
          IF ( ANYBIN(1,IS) .AND. ANYBIN(MDC,IS) .AND.
     &         IDTOT .EQ. MDC ) PERIOD = .TRUE.
!
          DO IDDUM = IDDLOW, IDDTOP
            IDSWAN = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDBAND = MOD ( IDDUM - IDDLOW  , MDC ) + 1
            IDDL   =  MOD ( IDCMIN(IS) - 1 + MDC , MDC ) + 1
            IDDT   =  MOD ( IDCMAX(IS) - 1 + MDC , MDC ) + 1
            NPP = (IDBAND - 1) * ISSTOP + IS
            IF ( IDSWAN .EQ. MDC .AND. PERIOD ) THEN
              BAND(NPP , 3 )  = IMATLA(IDSWAN , IS)
              LOPERI(IS)      = IMATUA(IDSWAN , IS)
            ELSE IF ( IDSWAN .EQ. 1 .AND. PERIOD ) THEN
              BAND(NPP , 8 )  = IMATUA(IDSWAN , IS)
              UPPERI(IS)      = IMATLA(IDSWAN , IS)
            ELSE IF ( IDSWAN .EQ. IDDL ) THEN
              BAND(NPP , 8 )  = IMATUA(IDSWAN , IS)
            ELSE IF ( IDSWAN .EQ. IDDT ) THEN
              BAND(NPP , 3 )  = IMATLA(IDSWAN , IS)
            ELSE
              BAND(NPP , 3 )  = IMATLA(IDSWAN , IS)
              BAND(NPP , 8 )  = IMATUA(IDSWAN , IS)
            ENDIF
          ENDDO
        ELSE IF ( SECTOR(IS) .EQ. 1 .OR. SECTOR(IS) .EQ. 4 ) THEN
!
!         *** minimum counter = 1, maximum counter = MDC ***
!
          IF ( IDDLOW .NE. 1 .OR. IDDTOP .NE. MDC ) THEN
            WRITE(PRTEST,*) 'Error in VULMT1', IS, SECTOR(IS),            20.44
     &                      IDDLOW, IDDTOP                                20.44
          ENDIF
          DO IDDUM = IDDLOW, IDDTOP
            IDSWAN = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDBAND = MOD ( IDDUM - IDDLOW  , MDC ) + 1
            NPP = (IDBAND - 1) * ISSTOP + IS
            IF ( IDDUM .EQ. IDDLOW ) THEN
!             BAND(NPP , 8 )  = IMATUA(IDDUM , IS)             replaced
!             UPPERI(IS)      = IMATLA(IDDUM , IS)
              BAND(NPP , 8 )  = IMATUA(IDSWAN, IS)                        20.44
              UPPERI(IS)      = IMATLA(IDSWAN, IS)                        20.44
            ELSE IF ( IDDUM .EQ. IDDTOP) THEN
!             BAND(NPP , 3 )  = IMATLA(IDDUM , IS)             replaced
!             LOPERI(IS)      = IMATUA(IDDUM , IS)
              BAND(NPP , 3 )  = IMATLA(IDSWAN, IS)                        20.44
              LOPERI(IS)      = IMATUA(IDSWAN, IS)                        20.44
            ELSE
!             BAND(NPP , 3 )  = IMATLA(IDDUM , IS)             replaced
!             BAND(NPP , 8 )  = IMATUA(IDDUM , IS)
              BAND(NPP , 3 )  = IMATLA(IDSWAN, IS)                        20.44
              BAND(NPP , 8 )  = IMATUA(IDSWAN, IS)                        20.44
            ENDIF
          ENDDO
        ENDIF
!
!       *** test output ***
!
        IF ( TESTFL .AND. ITEST .GE. 120) THEN
          WRITE(PRTEST,*)
          WRITE(PRTEST,1056) IS, IDCMIN(IS), IDCMAX(IS), IDDLOW,
     &                       IDDTOP
 1056     FORMAT(' VULMT1: IS IDCMIN IDCMAX IDDLOW IDDTOP:',5I4)
          WRITE(PRTEST,1058) IDDL, IDDT, SECTOR(IS),IDTOT
 1058     FORMAT(' VULMT1: IDDL IDDT SECTOR IDTOT        :',4I4)
          WRITE(PRTEST,1059) ANYBIN(2,IS),ANYBIN(1,IS),
     &                       ANYBIN(MDC,IS),PERIOD
 1059     FORMAT(' VULMT1: ANYBIN  2   1   MDC  PERIOD   :',4L3)
          WRITE(PRTEST,1060) ANYBIN(1,IS),ANYBIN(MDC,IS),
     &                       ANYBIN(MDC-1,IS),PERIOD
 1060     FORMAT(' VULMT1: ANYBIN  1  MDC   MDC-1  PERIOD:',4L3)
        ENDIF
!
      ENDDO
!
!     *** information about SWAN matrices and SOLBAND matrices ***
!
      IF ( TESTFL .AND. ITEST .GE. 120) THEN
        WRITE(PRTEST,*)
        WRITE(PRTEST,*) '  Subroutine VULMT1'
        WRITE(PRTEST,*)
        WRITE(PRTEST,111) KCGRD(1),MDC, MSC                               30.21
 111    FORMAT(' VULMT1 : POINT MDC MSC              :',3I5)
        WRITE(PRTEST,211) IDDLOW, IDDTOP, ISSTOP
 211    FORMAT(' VULMT1 : IDDLOW IDDTOP ISSTOP       :',3I4)
        WRITE(PRTEST,*)
        WRITE(PRTEST,*) ' matrix coefficients in SWAN  '
        WRITE(PRTEST,*)
        WRITE(PRTEST,*)
     & 'IS ID IDB  IMATLA   IMATDA   IMATUA    IMATRA   IMAT5L   IMAT6U'
        DO IDDUM = IDDLOW, IDDTOP
          ID     = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IDBAND = MOD ( IDDUM - IDDLOW  , MDC ) + 1
           DO IS = 1, ISSTOP
            WRITE(PRTEST,2101) IS, ID , IDBAND, IMATLA(ID,IS),
     &                         IMATDA(ID,IS), IMATUA(ID,IS),
     &                         IMATRA(ID,IS), IMAT5L(ID,IS),
     &                         IMAT6U(ID,IS)
2101        FORMAT(3I3,6E10.2)
          ENDDO
          WRITE(PRTEST,*)
        ENDDO
        WRITE(PRTEST,*)
        WRITE(PRTEST,*) ' matrix coefficients for CGSTAB solver'
        WRITE(PRTEST,*)
        WRITE(PRTEST,*)
     & 'IS ID IDB     (3)     (5)       (1)      (6)      (8)      RHV'
        DO IDDUM = IDDLOW, IDDTOP
          IDSWAN = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IDBAND = MOD ( IDDUM - IDDLOW  , MDC ) + 1
          DO IS = 1, ISSTOP
            IPP = ( IDBAND - 1 ) * ISSTOP + IS
            WRITE(PRTEST,3101) IS ,IDSWAN ,IDBAND, BAND(IPP,3),
     &                        BAND(IPP,5), BAND(IPP,1), BAND(IPP,6),
     &                        BAND(IPP,8), RHV(IPP)
3101        FORMAT(3I3,6E10.2)
          ENDDO
          WRITE(PRTEST,*)
        ENDDO
        WRITE(PRTEST,*)'IS ID      LPER     UPER '
        DO IDDUM = IDDLOW, IDDTOP
          IDSWAN = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IDBAND = MOD ( IDDUM - IDDLOW  , MDC ) + 1
          IF ( IDSWAN .EQ. 1  .OR. IDSWAN .EQ. MDC ) THEN
            DO IS = 1, ISSTOP
              IPP = ( IDBAND - 1 ) * ISSTOP + IS
              WRITE(PRTEST,3141) IS ,IDSWAN , LOPERI(IS),
     &                           UPPERI(IS)
3141          FORMAT(2I3,2E10.2)
            ENDDO
          END IF
        ENDDO
      END IF
!
      RETURN
      END
!
!**********************************************************************
!*                                                                    *
      SUBROUTINE SWCOVA2D ( MXC, MYC, XCG, YCG, CVA )
!*                                                                    *
!**********************************************************************
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     31.00  Kees Kassels
!     31.03  Annette Kieftenburg
!
!  1. Updates
!
!  version 1.0    date 09-09-1997
!
!  2. Purpose
!
!    Compute covariant base vectors in integration points
!    two-dimensional case
!
!  3. Method
!
!     The covariant basis vectors are given by:
!      alpha                    beta
!     a       = dx       / d xsi       = cva ( i, j, alpha, beta, k )
!      (beta)     alpha
!
!     i,j refers to the cell index and k (pnttyp) to the position of the point
!     in a cell
!
!     To evaluate the base vectors a central difference formula is applied
!
!     This leads to:
!
!     pnttype                 d/dksi1                  d/dksi2
!
!        1             x(i+1,j)-x(i,j)            (x(i+1,j+1)-x(i+1,j-1))/4+
!                                                 (x(i,j+1)-x(i,j-1))/4
!        2         (x(i+1,j+1)-x(i-1,j+1))/4+      x(i,j+1)-x(i,j)
!                    (x(i+1,j)-x(i-1,j))/4
!
!
!      *-----------*
!      |           |
!      2           |
!      |           |
!      |           |
!      *--- 1------*
!
!
!  4. Argument variables
!
!     CVA       o    Array containing the covariant basis vectors
!                    in 2D we have
!                                                 l
!                    cva(i,j,l,alpha,p) contains a        in cell i,j in
!                                                 (alpha)       point type p
!     MXC       i    Number of points in the x-direction
!     MYC       i    Number of points in the y-direction
!     XCG       i    x-coordinates
!     YCG       i    y-coordinates
!
      IMPLICIT NONE
      INTEGER MXC, MYC
      REAL    XCG(1:MXC,1:MYC),                                           31.03
     &        YCG(1:MXC,1:MYC),
     &        CVA(1:MXC,1:MYC,1:2,1:2,1:2)
!
!  5. Parameter variables
!
!  6. Local variables
!
!     I         General loop variable
!     J         General loop variable
!
      INTEGER I, J
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SETUP2D   Computation of SETUP, the change of waterlevel by waves.
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
! **********************************************************************
!                       KEYWORDS
!    base_vectors
!    interpolation
!
! **********************************************************************
!
!     --- pnttype = 1:    x(i+1,j)-x(i,j)
!                        (x(i+1,j+1)-x(i+1,j-1))/4+(x(i,j+1)-x(i,j-1))/4
!
!         --- whole area
!
      DO J = 1, MYC
         DO I = 1, MXC-1
            CVA(I,J,1,1,1) = XCG(I+1,J)-XCG(I,J)
            CVA(I,J,2,1,1) = YCG(I+1,J)-YCG(I,J)
         END DO
      END DO
!
!     --- inner area
!
      DO J = 2, MYC-1
         DO I = 1, MXC-1
            CVA(I,J,1,2,1) = 0.25E0 * ( XCG(I+1,J+1)-XCG(I+1,J-1)+
     &                                  XCG(I  ,J+1)-XCG(I  ,J-1) )
            CVA(I,J,2,2,1) = 0.25E0 * ( YCG(I+1,J+1)-YCG(I+1,J-1)+
     &                                  YCG(I  ,J+1)-YCG(I  ,J-1) )
         END DO
      END DO
!
!     --- lower boundary
!
      J = 1
      DO I = 1, MXC-1
         CVA(I,J,1,2,1) = 0.25E0 * ( XCG(I+1,J+1)-
     &                              (2E0*XCG(I+1,J)-XCG(I+1,J+1)) +
     &                               XCG(I  ,J+1)-
     &                              (2E0*XCG(I  ,J)-XCG(I  ,J+1)) )
         CVA(I,J,2,2,1) = 0.25E0 * ( YCG(I+1,J+1)-
     &                              (2E0*YCG(I+1,J)-YCG(I+1,J+1)) +
     &                               YCG(I  ,J+1)-
     &                              (2E0*YCG(I  ,J)-YCG(I  ,J+1)) )
      END DO
!
!     --- upper boundary
!
      J = MYC
      DO I = 1, MXC-1
         CVA(I,J,1,2,1) = 0.25E0 * ((2E0*XCG(I+1,J)-XCG(I+1,J-1)) -
     &                               XCG(I+1,J-1)                 +
     &                              (2E0*XCG(I  ,J)-XCG(I  ,J-1)) -
     &                               XCG(I  ,J-1) )
         CVA(I,J,2,2,1) = 0.25E0 * ((2E0*YCG(I+1,J)-YCG(I+1,J-1)) -
     &                               YCG(I+1,J-1)                 +
     &                              (2E0*YCG(I  ,J)-YCG(I  ,J-1)) -
     &                               YCG(I  ,J-1) )
      END DO
!
!     --- pnttype = 2:  (x(i+1,j+1)-x(i-1,j+1))/4+(x(i+1,j)-x(i-1,j))/4
!                        x(i,j+1)-x(i,j)
!
!         --- whole area
!
      DO J = 1, MYC-1
         DO I = 1, MXC
            CVA(I,J,1,2,2) = XCG(I,J+1)-XCG(I,J)
            CVA(I,J,2,2,2) = YCG(I,J+1)-YCG(I,J)
         END DO
      END DO
!
!     --- inner area
!
      DO J = 1, MYC-1
         DO I = 2, MXC-1
            CVA(I,J,1,1,2) = 0.25E0 * ( XCG(I+1,J+1)-XCG(I-1,J+1)+
     &                                  XCG(I+1,J  )-XCG(I-1,J  ) )
            CVA(I,J,2,1,2) = 0.25E0 * ( YCG(I+1,J+1)-YCG(I-1,J+1)+
     &                                  YCG(I+1,J  )-YCG(I-1,J  ) )
         END DO
      END DO
!
!     --- left boundary
!
      I = 1
      DO J = 1, MYC-1
         CVA(I,J,1,1,2) = 0.25E0 * ( XCG(I+1,J+1)-
     &                              (2E0*XCG(I,J+1)-XCG(I+1,J+1)) +
     &                               XCG(I+1,J  )-
     &                              (2E0*XCG(I  ,J)-XCG(I+1,J  )) )
         CVA(I,J,2,1,2) = 0.25E0 * ( YCG(I+1,J+1)-
     &                              (2E0*YCG(I,J+1)-YCG(I+1,J+1)) +
     &                               YCG(I+1,J  )-
     &                              (2E0*YCG(I  ,J)-YCG(I+1,J  )) )
      END DO
!
!     --- right boundary
!
      I = MXC
      DO J = 1, MYC-1
         CVA(I,J,1,1,2) = 0.25E0 * ((2E0*XCG(I,J+1)-XCG(I-1,J+1)) -
     &                               XCG(I-1,J+1)                 +
     &                              (2E0*XCG(I  ,J)-XCG(I-1,J  )) -
     &                               XCG(I-1,J  ) )
         CVA(I,J,2,1,2) = 0.25E0 * ((2E0*YCG(I,J+1)-YCG(I-1,J+1)) -
     &                               YCG(I-1,J+1)                 +
     &                              (2E0*YCG(I  ,J)-YCG(I-1,J  )) -
     &                               YCG(I-1,J  ) )
      END DO
!
      END
!
!**********************************************************************
!*                                                                    *
      SUBROUTINE SWJCTA2D ( MXC, MYC, CVA, JCTA )
!*                                                                    *
!**********************************************************************
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     31.00  Kees Kassels
!     31.03  Annette Kieftenburg
!
!  1. Updates
!
!     version 1.0    date 10-09-1997
!
!  2. Purpose
!
!    Compute sqrt(g) x contra variant base vectors in integration point
!    two-dimensional case
!
!  3. Method
!
!     sqrt(g) x contravariant base vector in 2D is given by:
!
!              (1)     2                    (2)      2
!     sqrt(g) a     = a            sqrt(g) a     = -a
!              1       (2)                  1        (1)
!
!              (1)     1                    (2)      1
!     sqrt(g) a     =-a            sqrt(g) a     =  a
!              2       (2)                  2        (1)
!
!  4. Argument variables
!
!     cva       i    Array containing the covariant basis vectors
!                    in 2D we have
!                                                 l
!                    cva(i,j,l,alpha,p) contains a        in cell i,j in
!                                                 (alpha)       point type p
!     jcta      o    Jacobian times contravariant basis vectors
!                    in point pnttyp=1 base vector 1
!                    in point pnttyp=2 base vector 2
!     MXC       i    Number of points in the x-direction
!     MYC       i    Number of points in the y-direction
!
      IMPLICIT NONE
      INTEGER MXC, MYC
      REAL    CVA(1:MXC,1:MYC,2,2,2)                                      31.03
      REAL   JCTA(1:MXC,1:MYC,2,2)                                        31.03
!
!  5. Parameter variables
!
!  6. Local variables
!
!     I         General loop variable
!     J         General loop variable
      INTEGER I, J
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SETUP2D   Computation of SETUP, the change of waterlevel by waves.
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
! **********************************************************************
!
!                       KEYWORDS
!
!    base_vectors
!    interpolation
! **********************************************************************
!
!     --- alphad = 1, means pointtype = 1, contravariant base vector 1
!
      DO J = 1, MYC
         DO I = 1, MXC-1
            JCTA(I,J,1,1) = CVA(I,J,2,2,1)
            JCTA(I,J,2,1) = -CVA(I,J,1,2,1)
         END DO
      END DO
!
!     --- alphad = 2, means pointtype = 2, contravariant base vector 2
!
      DO J = 1, MYC-1
         DO I = 1, MXC
            JCTA(I,J,1,2) = -CVA(I,J,2,1,2)
            JCTA(I,J,2,2) = CVA(I,J,1,1,2)
         END DO
      END DO
!
      END
!
! **********************************************************************
!                                                                      *
      SUBROUTINE SWTRAD2D ( MXC, MYC, WFRCX, WFRCY, DEPMIN,               31.03
     &                      ALPHAD, DEPTH, CVA, JCTA, CVC, CTC,           31.03
     &                      DTSUM, RHSIDE)
!                                                                      *
! **********************************************************************
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     31.03  Annette Kieftenburg
!
!  1. Updates
!
!     --
!
!  2. Purpose
!
!     Compute contribution of diffusion term in R2 for a transport equation
!     per integration point
!     Compute righthandside
!
!  3. Method (updated...)
!
!      The diffusion terms in the WesBeek discretization are given by
!
!
!                   (a)  (b)     |+e(a)/2
!          sqrt(g) a   .C    T,  |
!                  -    -      b |-e(a)/2
!
!             (b)
!      where C    are the so-called WesBeek vectors defined in the following
!            -
!      way:
!
!      In principle, we want to approximate, for example, the following
!      diffusion term
!
!                (1)    _
!       sqrt(g) a   . k VT  in point (i+1/2,j)
!               -
!            _          1       2                          (1)
!      where VT = (dT/dx , dT/dx  ) and note that sqrt(g) a
!                                                         -
!                          1
!      is continue along xi  = constant.
!
!
!        |(i+1,j)         _        _
!      T |           = S  VT dx =  VT|           . S  dx
!        |(i,j)               -       (i+1/2,j)        -
!
!      (note: S is denoted as a integral!)
!
!      The last integral can be approximated as follows:
!
!
!
!      C    = S  dx  =   a     |
!      -(1)       -      -(1)    (i,j)
!
!      In order to find one more equation, another integration path has
!      to be taken. With this, we can calculate the other covariant WesBeek
!      vector. These are:
!
!
!
!      C    =    a                +  a             +
!      -(2)      -(2) |  (i,j-1 )    -(2) | (i,j )
!
!
!                a                +   a
!                -(2) | (i+1,j-1)     -(2) |(i+1,j)
!
!                                 _
!      Solving two equations for  VT results in
!
!       _                 (a)
!       VT|            = C    T   |
!          (i+1/2,j)     -     ,a  (i+1/2,j)
!                                                  (a)
!      Finally, the contravariant WesBeek vectors C    can be computed
!                                                 -
!      with the following formulae
!
!
!       (1)    1    2       1   T       (2)    1     2     1    T
!      C     = - ( C    , -C   )     C     =  - ( -C   ,  C   )
!      -       C    (2)     (2)      -               (1)    (1)
!
!
!                2    1        2    1
!      with C = C    C     -  C    C
!                (1)   (2)     (1)  (2)
!
!
!  4. Argument variables
!
!     ALPHAD    i    Direction index of integration.
!     CTC       i    Work array containing the contravariant WESBEEK vectors
!     CVA       o    Array containing the covariant basis vectors
!                    in 2D we have
!                                                 l
!                    CVA(i,j,l,alpha,p) contains a        in cell i,j in
!                                                 (alpha)       point type p
!     CVC       i    Work array containing the covariant WESBEEK vectors
!     DEPMIN    i    Minimum depth
!     DEPTH     i    depth direct addressed                               31.03
!     DTSUM     o    Derivative contributions to the matrix
!     JCTA      i    Jacobian times contravariant basis vectors
!                    in point pnttyp=1 base vector 1
!                    in point pnttyp=2 base vector 2
!     MXC       i    Number of points in the x-direction
!     MYC       i    Number of points in the y-direction
!     RHSIDE    o    Righthandside
!     WFRCX     i    force x-component direct addressed                   31.03
!     WFRCY     i    force y-component direct addressed                   31.03
!
      IMPLICIT NONE
      INTEGER ALPHAD, MXC, MYC                                            31.03
      REAL   CTC(1:MXC, 1:MYC,2,2),
     &       CVA(1:MXC, 1:MYC,2,2,2),
     &       CVC(1:MXC, 1:MYC,2,2),
     &    DEPMIN,                                                         31.03
     &     DEPTH(1:MXC, 1:MYC),                                           31.03
     &     DTSUM(1:MXC, 1:MYC,2),
     &      JCTA(1:MXC, 1:MYC,2,2),
     &    RHSIDE(1:MXC, 1:MYC),
     &     WFRCX(1:MXC, 1:MYC),                                           31.03
     &     WFRCY(1:MXC, 1:MYC)                                            31.03

!
!  5. Parameter variables
!
!  6. Local variables
!
!     BETAD     Direction index, usually the direction perpendicular to the
!               integration direction of the boundary integral.
!               For example betad in local point (1,0) is 1 and in (0,1): 2
!     DET       Determinant
!     DX        Point location in cell relative to local
!               coordinates with (0,0) and (2,2) or
!               (0,0,0) and (2,2,2) as opposite cell corners
!     FACT1     Help factor for computing integration path
!     FACT2     Help factor for computing integration path
!     FACT3     Help factor for computing integration path
!     FACT4     Help factor for computing integration path
!     FACTOR    Help factor
!     GAMMAD    Direction index
!     I         General loop variable
!     IP        Local counting variable i-direction
!     IPAREA    Area of integration points where a contribution to
!               the matrix is to be calculated
!     IP1       Shift for boundaries 1 and 2
!     IP2       Shift for boundaries 3 and 4
!     JP        Local counting variable j-direction
!
      INTEGER BETAD, DX(2,2), GAMMAD, I, IP, IPAREA(2,2),
     &        IP1, IP2, JP
      REAL    DET, FACTOR, FACT1, FACT2, FACT3, FACT4                     31.03
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SETUP2D   Computation of SETUP, the change of waterlevel by waves.
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
!                       KEYWORDS
!
!     coefficients
!     discretization
!     transport_matrix
! **********************************************************************
!
!     --- Wesbeek discretization
!
      IF ( ALPHAD.EQ.1 ) THEN
         IPAREA(1,1) = 1
         IPAREA(2,1) = MXC-1
         IPAREA(1,2) = 1
         IPAREA(2,2) = MYC
      ELSE
         IPAREA(1,1) = 1
         IPAREA(2,1) = MXC
         IPAREA(1,2) = 1
         IPAREA(2,2) = MYC-1
      END IF
!
      DX(1,1)      =  0
      DX(1,2)      =  0
      DX(1,ALPHAD) =  1
!
      DO BETAD = 1, 2
         DO GAMMAD = 1, 2
            DX(2,1)      =  0
            DX(2,2)      =  0
            DX(2,GAMMAD) =  1
!
!           --- INNER AREA
!
            IF ( GAMMAD.EQ.ALPHAD ) THEN
!
               DO JP = IPAREA(1,2), IPAREA(2,2)
                  DO IP = IPAREA(1,1), IPAREA(2,1)
                     CVC(IP,JP,BETAD,GAMMAD) =
     &                            CVA(IP,JP,BETAD,GAMMAD,ALPHAD)
                  END DO
               END DO
!
            ELSE
!
            IP1 = 0
            IP2 = 0
            IF ( ALPHAD.EQ.1 ) THEN
               IP1 = 1
            ELSE
               IP2 = 1
            END IF
!
               DO JP = IPAREA(1,2)+IP1, IPAREA(2,2)-IP1
                  DO IP = IPAREA(1,1)+IP2, IPAREA(2,1)-IP2
!
                     FACT1 = 1E0
                     FACT2 = 1E0
                     FACT3 = 1E0
                     FACT4 = 1E0
                     IF(DEPTH(IP+DX(1,1)-DX(2,1),JP+DX(1,2)-DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT1 = 0E0
                     IF(DEPTH(IP+DX(1,1)+DX(2,1),JP+DX(1,2)+DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT2 = 0E0
                     IF(DEPTH(IP        -DX(2,1),JP        -DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT3 = 0E0
                     IF(DEPTH(IP        +DX(2,1),JP        +DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
!
     &               FACT4 = 0E0
!
                     CVC(IP,JP,BETAD,GAMMAD) =
     &            FACT1*CVA(IP+DX(1,1)-DX(2,1),JP+DX(1,2)-DX(2,2),
     &                  BETAD,GAMMAD,3-ALPHAD)
     &
     &           +FACT2*CVA(IP+DX(1,1)        ,JP+DX(1,2)        ,
     &                  BETAD,GAMMAD,3-ALPHAD)
     &
     &           +FACT3*CVA(IP        -DX(2,1),JP        -DX(2,2),
     &                  BETAD,GAMMAD,3-ALPHAD)
     &
     &           +FACT4*CVA(IP,JP,BETAD,GAMMAD,3-ALPHAD)
!
                  END DO
               END DO
!
            END IF
!
!           --- BOUNDARY ONE
!               --- at the boundary we deal with half a cell : fact=5e-1
!
            IF ( ALPHAD.EQ.2 ) THEN
!
               IF ( GAMMAD.NE.ALPHAD ) THEN
!
                  IP = 1
                  DO JP = 1, MYC-1
                     FACT2 = 2E0
                     FACT4 = 2E0
                     IF(
     &               DEPTH(IP+DX(1,1)+DX(2,1),JP+DX(1,2)+DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT2 = 0E0
                     IF(
     &               DEPTH(IP        +DX(2,1),JP        +DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT4 = 0E0
!
                     CVC(IP,JP,BETAD,GAMMAD) =
     &           +FACT2*CVA(IP+DX(1,1)        ,JP+DX(1,2)        ,
     &                  BETAD,GAMMAD,3-ALPHAD)
     &
     &           +FACT4*CVA(IP,JP,BETAD,GAMMAD,3-ALPHAD)
                  END DO
!
               END IF
!
            END IF
!
!           --- BOUNDARY TWO
!
            IF ( ALPHAD.EQ.2 ) THEN
!
               IF ( GAMMAD.NE.ALPHAD ) THEN
!
                  IP = MXC
                  DO JP = 1, MYC-1
                     FACT1 = 2E0
                     FACT3 = 2E0
                     IF(
     &               DEPTH(IP+DX(1,1)-DX(2,1),JP+DX(1,2)-DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT1 = 0E0
                     IF(
     &               DEPTH(IP        -DX(2,1),JP        -DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT3 = 0E0
!
                     CVC(IP,JP,BETAD,GAMMAD) =
     &            FACT1*CVA(IP+DX(1,1)-DX(2,1),JP+DX(1,2)-DX(2,2),
     &                  BETAD,GAMMAD,3-ALPHAD)
!
     &           +FACT3*CVA(IP        -DX(2,1),JP        -DX(2,2),
     &                  BETAD,GAMMAD,3-ALPHAD)
                   END DO
!
               END IF
!
            END IF
!
!           --- BOUNDARY THREE
!
            IF ( ALPHAD.EQ.1 ) THEN
!
               IF ( GAMMAD.NE.ALPHAD ) THEN
!
                  JP = 1
                  DO IP = 1, MXC-1
                     FACT2 = 2E0
                     FACT4 = 2E0
                     IF(
     &               DEPTH(IP+DX(1,1)+DX(2,1),JP+DX(1,2)+DX(2,2))
     &                    .LT.DEPMIN)                                       31.03
     &               FACT2 = 0E0
                     IF(
     &               DEPTH(IP        +DX(2,1),JP        +DX(2,2))
     &                    .LT.DEPMIN)                                       31.03
     &               FACT4 = 0E0
!
                     CVC(IP,JP,BETAD,GAMMAD) =
     &            FACT2*CVA(IP+DX(1,1)        ,JP+DX(1,2)        ,
     &                  BETAD,GAMMAD,3-ALPHAD)
     &
     &           +FACT4*CVA(IP,JP,BETAD,GAMMAD,3-ALPHAD)
                 END DO
!
               END IF
!
            END IF
!
!           --- BOUNDARY FOUR
!
            IF ( ALPHAD.EQ.1 ) THEN
!
               IF ( GAMMAD.NE.ALPHAD ) THEN
!
                  JP = MYC
                  DO IP = 1, MXC-1
                     FACT1 = 2E0
                     FACT3 = 2E0
                     IF(
     &               DEPTH(IP+DX(1,1)-DX(2,1),JP+DX(1,2)-DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT1 = 0E0
                     IF(
     &               DEPTH(IP        -DX(2,1),JP        -DX(2,2))
     &                  .LT.DEPMIN)                                       31.03
     &               FACT3 = 0E0
!
                        CVC(IP,JP,BETAD,GAMMAD) =
     &            FACT1*CVA(IP+DX(1,1)-DX(2,1),JP+DX(1,2)-DX(2,2),
     &                     BETAD,GAMMAD,3-ALPHAD)
     &           +FACT3*CVA(IP        -DX(2,1),JP        -DX(2,2),
     &                     BETAD,GAMMAD,3-ALPHAD)
                   END DO
!
               END IF
!
            END IF
!
         END DO
      END DO
!
!     --- Determine contravariant WESBEEK vectors
!
      DO JP=IPAREA(1,2), IPAREA(2,2)
         DO IP=IPAREA(1,1), IPAREA(2,1)
!
            IF ( DEPTH(IP,JP).LT.DEPMIN .OR.                              31.03
     &           DEPTH(IP+DX(1,1),JP+DX(1,2)).LT.DEPMIN) THEN             31.03
!
!           --- non existent neighbour, no contribution from this
!               integration point, making ctc = 0e0, which gives
!               also dtsum = 0e0
!
               CTC(IP,JP,1,1) =  0E0
               CTC(IP,JP,2,1) =  0E0
               CTC(IP,JP,1,2) =  0E0
               CTC(IP,JP,2,2) =  0E0
!
            ELSE
!
!           --- First, calculate determinant ...
!
               DET = CVC(IP,JP,1,1)*CVC(IP,JP,2,2) -
     &               CVC(IP,JP,2,1)*CVC(IP,JP,1,2)
!
!              --- THEN contravariant WESBEEK vectors
!
               IF ( DET .LE. 0) THEN

                  CTC(IP,JP,1,1) =  0E0
                  CTC(IP,JP,2,1) =  0E0
                  CTC(IP,JP,1,2) =  0E0
                  CTC(IP,JP,2,2) =  0E0
               ELSE
                  CTC(IP,JP,1,1) =  CVC(IP,JP,2,2)/DET
                  CTC(IP,JP,2,1) = -CVC(IP,JP,1,2)/DET
                  CTC(IP,JP,1,2) = -CVC(IP,JP,2,1)/DET
                  CTC(IP,JP,2,2) =  CVC(IP,JP,1,1)/DET
               END IF
!
            END IF
!
         END DO
      END DO
!
!     --- initialize DTSUM
!
      DO I = 1, 2
         DO JP = IPAREA(1,2), IPAREA(2,2)
            DO IP = IPAREA(1,1), IPAREA(2,1)
               DTSUM(IP,JP,I) = 0E0
            END DO
         END DO
      END DO
!
      DO BETAD = 1, 2
!
         DO JP = IPAREA(1,2), IPAREA(2,2)
            DO IP = IPAREA(1,1), IPAREA(2,1)
!
               DO I = 1, 2
                  DTSUM(IP,JP,BETAD) = DTSUM(IP,JP,BETAD) -
     &                                 5E-1*(DEPTH(IP,JP) +
     &                                 DEPTH(IP+DX(1,1),JP+DX(1,2)) ) *   31.03
     &                                 CTC(IP,JP,I,BETAD)*
     &                                 JCTA(IP,JP,I,ALPHAD)
               END DO
!
            END DO
         END DO
!
      END DO
!
!     --- righthandside inner area
!
      IP1 = 0
      IP2 = 0
      IF ( ALPHAD.EQ.1 ) THEN
         IP1 = 1
      ELSE
         IP2 = 1
      END IF
!
      DO JP = IPAREA(1,2)+IP1, IPAREA(2,2)-IP1
         DO IP = IPAREA(1,1)+IP2, IPAREA(2,1)-IP2
!
            IF ( DEPTH(IP,JP).LT.DEPMIN .OR.                              31.03
     &           DEPTH(IP+DX(1,1),JP+DX(1,2)).LT.DEPMIN ) THEN            31.03
!
!           --- non existent neighbour, no contribution from this
!               integration point
!
               FACTOR = 0E0
            ELSE
               FACTOR = 5E-1*(-WFRCX(IP,JP) -                             31.03 -
     &                        WFRCX(IP+DX(1,1),JP+DX(1,2)) ) *            31.03 -
     &                        JCTA(IP,JP,1,ALPHAD) +
     &                  5E-1*(-WFRCY(IP,JP) -                             31.03 -
     &                        WFRCY(IP+DX(1,1),JP+DX(1,2)) ) *            31.03 -
     &                        JCTA(IP,JP,2,ALPHAD)

            END IF
            RHSIDE(IP,JP) = RHSIDE(IP,JP) + FACTOR
            RHSIDE(IP+DX(1,1),JP+DX(1,2)) =
     &         RHSIDE(IP+DX(1,1),JP+DX(1,2)) - FACTOR
         END DO
      END DO
!

!     --- righthandside  boundary : values * 5e-1, because of half cells
!
       IF ( ALPHAD .EQ. 1) THEN
!
        DO JP = 1, MYC, MYC-1
           DO IP = 1, MXC-1

              IF ( DEPTH(IP,JP).LT.DEPMIN .OR.                            31.03
     &             DEPTH(IP+DX(1,1),JP+DX(1,2)).LT.DEPMIN ) THEN          31.03
!
!             --- non existent neighbour, no contribution from this
!                 integration point
!
                 FACTOR = 0E0
              ELSE
                 FACTOR = 25E-2*(-WFRCX(IP,JP) -                          31.03 -
     &                           WFRCX(IP+DX(1,1),JP+DX(1,2)) ) *         31.03 -
     &                           JCTA(IP,JP,1,ALPHAD) +
     &                    25E-2*(-WFRCY(IP,JP) -                          31.03 -
     &                           WFRCY(IP+DX(1,1),JP+DX(1,2)) ) *         31.03 -
     &                           JCTA(IP,JP,2,ALPHAD)
              END IF
!
              RHSIDE(IP,JP) = RHSIDE(IP,JP) + FACTOR
              RHSIDE(IP+DX(1,1),JP+DX(1,2)) =
     &           RHSIDE(IP+DX(1,1),JP+DX(1,2)) - FACTOR
            END DO
        END DO

       ELSE
!
         DO JP = 1, MYC-1
            DO IP = 1, MXC, MXC-1

               IF ( DEPTH(IP,JP).LT.DEPMIN .OR.                           31.03
     &              DEPTH(IP+DX(1,1),JP+DX(1,2)).LT.DEPMIN) THEN          31.03
!
!             --- non existent neighbour, no contribution from this
!                 integration point
!
                  FACTOR = 0E0
               ELSE
                  FACTOR = 25E-2*(-WFRCX(IP,JP) -                         31.03 -
     &                           WFRCX(IP+DX(1,1),JP+DX(1,2)) ) *         31.03 -
     &                           JCTA(IP,JP,1,ALPHAD) +
     &                     25E-2*(-WFRCY(IP,JP) -                         31.03 -
     &                           WFRCY(IP+DX(1,1),JP+DX(1,2)) ) *         31.03 -
     &                           JCTA(IP,JP,2,ALPHAD)
               END IF
               RHSIDE(IP,JP) = RHSIDE(IP,JP) + FACTOR
               RHSIDE(IP+DX(1,1),JP+DX(1,2)) =
     &            RHSIDE(IP+DX(1,1),JP+DX(1,2)) - FACTOR
             END DO
         END DO
      END IF
!
      END
!
!**********************************************************************
!*                                                                    *
      SUBROUTINE SWDISDT2 ( MXC, MYC, DEPTH, DEPMIN, ALPHAD, MATRIX,
     &                      DTSUM )                                       31.03
!*                                                                    *
!**********************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     31.00  Kees Kassels
!     31.03  Annette Kieftenburg
!
!  1. Updates
!
!     version 1.0    date   26-08-1997
!
!  2. Purpose
!
!     Distribute diffusion term for tranport equation in R2
!
!  3. Method
!
!    In DTSUM the contribution to the matrix is stored per
!    integration point. Take for example integration point IP.
!    The value dtsum(i,j,1) has to be given
!    to the matrix points 1 and 2. The value dtsum(i,j,2) has to be given
!    to the matrix points 3, 4, 5 and 6.
!    The value per integration point is used for two cells at the same
!    time, for the left cell with A  + SIG, for the rigth cell
!    with A  - SIG.
!
!                         4           3
!            *-----------*-----------*
!            |           |           |
!            |       _ _ |_ _  _ _ _ |
!            |      |    |    |      |
!            |      |    |2   |      |1
!            *------|----*--- IP-----*
!            |      |    |    |      |
!            |      |_ __|__ _|_ _ _ |
!            |           |           |
!            |           |6          |5
!            *-----------*-----------*
!
!
!    As the unknowns per cell are numbered as follows
!
!             7           8           9
!            *-----------*-----------*
!            |           |           |
!            |       _ _ |_ _  _ _ _ |
!            |      |    |    |      |
!            |5     |    |1   |      |6
!            *------|----*--- IP-----*
!            |      |    |    |      |
!            |      |_ __|__ _|_ _ _ |
!            |           |           |
!            |2          |3          |4
!            *-----------*-----------*
!
!    the contribution for the left cell are going to the points
!    1, 6, 3, 4, 8, 9
!    for the rigth cell
!    5, 1, 2, 3, 6, 8
!
!
!
!  4. Argument variables
!
!     ALPHAD    i    Direction index of integration.
!     DTSUM     i    Derivative contributions to the matrix
!     MATRIX   i/o   Matrix
!     MXC       i    Number of points in the x-direction
!     MYC       i    Number of points in the y-direction
!     DEPMIN    i    Minimum possible depth                               31.03
!
      INTEGER ALPHAD, MXC, MYC
      REAL    DEPMIN,                                                     31.03
     &        DEPTH(1:MXC,1:MYC),                                         31.03
     &        MATRIX(1:MXC,1:MYC,1:9),                                    31.03
     &        DTSUM(1:MXC,1:MYC,1:2)                                      31.03
!
!  5. Parameter variables
!
!  6. Local variables
!
!     BETAD     Direction index, usually the direction perpendicular to the
!               integration direction of the boundary integral.
!               For example betad in local point (1,0) is 1 and in (0,1): 2
!     IP        Local counting variable i-direction
!     ISHIFT    Shift in i-direction
!     JP        Local counting variable j-direction
!     MCT       Integer matrix information
!     OFSTNR    Offset diagonal number of matrix element that
!               corresponds to the current point of the interpolation
!               or extrapolation formula
!     SIG       sign for multiplication (1 or -1)
!
      INTEGER ishift, ip, jp, betad, sig, ofstnr(6), mct(-1:1,-1:1)
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SETUP2D   Computation of SETUP, the change of waterlevel by waves.
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
! **********************************************************************
!                       KEYWORDS
!
!     discretization
!     transport_matrix
! **********************************************************************
!
!                       DATA STATEMENTS
!
      DATA MCT / 2, 3, 4,
     &           5, 1, 6,
     &           7, 8, 9 /
! **********************************************************************
!
      BETAD = 3-ALPHAD
!
      IF ( ALPHAD.EQ.1) THEN
!
!     --- integration point 1
!
         DO JP = 1, MYC
            DO IP = 1,MXC-1
               DO ISHIFT = 0, 1
!
                  OFSTNR(1) = MCT(ISHIFT  , 0)
                  OFSTNR(2) = MCT(ISHIFT-1, 0)
                  OFSTNR(3) = MCT(ISHIFT  , 1)
                  OFSTNR(4) = MCT(ISHIFT-1, 1)
                  OFSTNR(5) = MCT(ISHIFT  ,-1)
                  OFSTNR(6) = MCT(ISHIFT-1,-1)
!
!                 --- boundaries, one side approximation
!
!                     --- lower boundary
!
                  IF ( JP.EQ.1) THEN
                    OFSTNR(5) = OFSTNR(1)
                    OFSTNR(6) = OFSTNR(2)
                  END IF
!
!                 --- upper boundary
!
                  IF ( JP.EQ.MYC) THEN
                    OFSTNR(3) = OFSTNR(1)
                    OFSTNR(4) = OFSTNR(2)
                  END IF
!
!                 --- non existent neighbours
!
                  IF ( JP.LT.MYC ) THEN
                     IF ( DEPTH(IP  ,JP+1).LT.DEPMIN)                     31.03
     &                  OFSTNR(4) = OFSTNR(2)
                     IF ( DEPTH(IP+1,JP+1).LT.DEPMIN)                     31.03
     &                  OFSTNR(3) = OFSTNR(1)
                  END IF
                  IF ( JP.GT.1 ) THEN
                     IF ( DEPTH(IP  ,JP-1).LT.DEPMIN)                     31.03
     &                  OFSTNR(6) = OFSTNR(2)
                     IF ( DEPTH(IP+1,JP-1).LT.DEPMIN)                     31.03
     &                  OFSTNR(5) = OFSTNR(1)
                  END IF
!
                  SIG = -1+2*ISHIFT
                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(1))=
     &                                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(1))
     &                                      + SIG * DTSUM(IP,JP,ALPHAD)
                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(2))=
     &                                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(2))
     &                                      - SIG * DTSUM(IP,JP,ALPHAD)
                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(3))=
     &                                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(3))
     &                                      + SIG * DTSUM(IP,JP,BETAD)
                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(4))=
     &                                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(4))
     &                                      + SIG * DTSUM(IP,JP,BETAD)
                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(5))=
     &                                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(5))
     &                                      - SIG * DTSUM(IP,JP,BETAD)
                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(6))=
     &                                  MATRIX(IP+1-ISHIFT,JP,OFSTNR(6))
     &                                      - SIG * DTSUM(IP,JP,BETAD)
               END DO
            END DO
         END DO
!
      ELSE
!
!     --- integration point 2
!
         DO JP = 1, MYC-1
            DO IP = 1,MXC
               DO ISHIFT = 0, 1
!
                  OFSTNR(1) = MCT( 0,ISHIFT  )
                  OFSTNR(2) = MCT( 0,ISHIFT-1)
                  OFSTNR(3) = MCT( 1,ISHIFT  )
                  OFSTNR(4) = MCT( 1,ISHIFT-1)
                  OFSTNR(5) = MCT(-1,ISHIFT  )
                  OFSTNR(6) = MCT(-1,ISHIFT-1)
!
!                 --- boundaries, one side approximation
!
!                     --- left boundary
!
                  IF ( IP.EQ.1) THEN
                    OFSTNR(5) = OFSTNR(1)
                    OFSTNR(6) = OFSTNR(2)
                  END IF
!
!                 --- right boundary
!
                  IF ( IP.EQ.MXC) THEN
                    OFSTNR(3) = OFSTNR(1)
                    OFSTNR(4) = OFSTNR(2)
                  END IF
!
!                 --- non existent neighbours
!
                  IF ( IP.LT.MXC ) THEN
                     IF ( DEPTH(IP+1,JP  ).LT.DEPMIN)                     31.03
     &                  OFSTNR(4) = OFSTNR(2)
                     IF ( DEPTH(IP+1,JP+1).LT.DEPMIN)                     31.03
     &                  OFSTNR(3) = OFSTNR(1)
                  END IF
                  IF ( IP.GT.1 ) THEN
                     IF ( DEPTH(IP-1,JP  ).LT.DEPMIN)                     31.03
     &                  OFSTNR(6) = OFSTNR(2)
                     IF ( DEPTH(IP-1,JP+1).LT.DEPMIN)                     31.03
     &                  OFSTNR(5) = OFSTNR(1)
                  END IF
!
                  SIG = -1+2*ISHIFT
                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(1))=
     &                                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(1))
     &                                      + SIG * DTSUM(IP,JP,ALPHAD)
                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(2))=
     &                                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(2))
     &                                      - SIG * DTSUM(IP,JP,ALPHAD)
                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(3))=
     &                                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(3))
     &                                      + SIG * DTSUM(IP,JP,BETAD)
                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(4))=
     &                                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(4))
     &                                      + SIG * DTSUM(IP,JP,BETAD)
                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(5))=
     &                                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(5))
     &                                      - SIG * DTSUM(IP,JP,BETAD)
                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(6))=
     &                                  MATRIX(IP,JP+1-ISHIFT,OFSTNR(6))
     &                                      - SIG * DTSUM(IP,JP,BETAD)
               END DO
            END DO
         END DO
!
      END IF
!
      END
!
!**********************************************************************
!                                                                     *
      SUBROUTINE SWESSBC ( MXC, MYC, MATRIX, RHSIDE, SETUP)               31.03
!                                                                     *
!**********************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     31.00  Kees Kassels
!     31.03  Annette Kieftenburg
!
!  1. Updates
!
!     version 1.0    date 11-09-1997
!
!  2. Purpose
!
!     Puts essential boundary conditions in matrix
!
!  3. Method (updated...)
!
!  4. Argument variables
!
!     MXC       i    Number of points in the x-direction
!     MYC       i    Number of points in the y-direction
!     MATRIX    i/o  Matrix
!     RHSIDE    i/o  Righthandside
!     SETUP     i    Unknown to be computed direct addressed
!
      INTEGER MXC, MYC
      REAL    MATRIX(1:MXC,1:MYC,1:9),
     &        RHSIDE(1:MXC,1:MYC),
     &        SETUP (1:MXC,1:MYC)                                         31.03
!
!  5. Parameter variables
!
!  6. Local variables
!
!     I         General loop variable
!     J         General loop variable
!     K         General loop variable
!
      INTEGER I, J, K
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SETUP2D   Computation of SETUP, the change of waterlevel by waves.
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
! **********************************************************************
!
!     --- boundary 1 and 2
!
      DO J = 1,MYC
        DO I = 1,MXC,MXC-1
           DO K = 2, 9
               MATRIX(I,J,K) = 0E0
           END DO
           MATRIX(I,J,1) = 1E0
           RHSIDE(I,J) = SETUP(I,J)
        END DO
      END DO
!
!      --- boundary 3 and 4
!
      DO J = 1,MYC,MYC-1
         DO I = 1,MXC
            DO K = 2, 9
               MATRIX(I,J,K) = 0E0
            END DO
            MATRIX(I,J,1) = 1E0
            RHSIDE(I,J) = SETUP(I,J)
         END DO
      END DO
!
      END
!
! ********************************************************************
!                                                                    *
      SUBROUTINE SWSOLV ( MATRIX, RHSIDE, SETUP, NPOINT,
     &                    WORK, NWORK, ITSW, ITER,
     &                    UPPERI, LOPERI)
!                                                                    *
! ********************************************************************
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
      IMPLICIT NONE
!
      INCLUDE 'swcomm3.inc'
!
!  0. Authors
!
!      1.0 : Kees Vuik
!     30.82: IJsbrand Haagsma
!     31.03: Annette Kieftenburg
!     34.01: Jeroen Adema
!     30.80: Nico Booij
!
!  1. Updates
!
!     version 1.0    date 12-09-97
!     34.01, Feb. 99: Introducing STPNOW
!     30.82, July 99: Corrected argumentlist SWSOLV
!     30.80, July 99: call of ISSOLV modified
!     30.82, Aug. 99: gets the information now from PNUMS file
!
!  2. Purpose
!
!     Prepare for ISSOLV
!
!  3. Method
!
!  4. Argument variables
!
!     ITER      input    Iteration number for SWAN
!     ITSW      input    Time step number
!     LOPERI             only relevant for computation in periodic domain 30.80
!     MATRIX    input    Matrix
!     NPOINT    input    Number of points MXC*MYC
!     NWORK     input    Dimension for work array
!     RHSIDE    input    Righthandside
!     SETUP     i/o      Unknown to be computed direct addressed
!     UPPERI             only relevant for computation in periodic domain 30.80
!     WORK               work array
!
      INTEGER  ITER,ITSW, NPOINT, NWORK                                   31.03
!
      REAL     LOPERI(*), MATRIX(*), RHSIDE(*), SETUP(*), UPPERI(*)       30.82
      REAL     WORK(*)                                                    31.03
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ERRPTS
!     I         General loop variable
!     IINSOL    Integer information for the solver
!     INFMAT    Integer information for the matrix
!     INOCNV
!     NCONCT    Number of connections in a row of the matrix
!     NPREC     Number of diagonals used in the preconditioner
!     RINSOL    Real information for the solver
!
      REAL      RINSOL(7)                                                 31.03
      INTEGER   INFMAT(10), IINSOL(14), I, NPREC, NCONCT
      INTEGER   INOCNV
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!     ISSOLV   The subroutine ISSOLV is used to solve an unsymmetric
!              system of equations of the shape A x = f.
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SETUP2D   Computation of SETUP, the change of waterlevel by waves.
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
! ********************************************************************
!
!     initialization
!
      INOCNV = 0
!
      DO I = 1, 14
         IINSOL(I) = 0
      ENDDO
      DO I = 1, 7
         RINSOL(I) = 0.
      ENDDO
!
      NPREC     = 9
      NCONCT    = 9
!
!     INFMAT(1) = 4 represents the original equation (in spectral domain)
!     INFMAT(1) = 5 represents the Poisson equation (to determine setup)  30.80
!
      INFMAT(1) = 5
!
      INFMAT(2) = MXC-1
      INFMAT(3) = MYC-1
      INFMAT(4) = 0
      INFMAT(5) = NPOINT
      INFMAT(6) = 9
!
!     INFMAT(8) contains LSETUP. This is used in DMLU3 to fill the
!     value of alpha for the RILU(alpha) preconditioner
!
      INFMAT(8) = LSETUP
!
!     input for the solver which should not be changed
!
      IINSOL(1) = 1
      IINSOL(4) = 0
      IINSOL(5) = 0
      IINSOL(7) = 4
      IINSOL(8) = 6
!
!     input for the solver which may be changed
!
!     IINSOL(2) determines the preconditioner
!     Possible values
!       0   no preconditioner
!      -1   diagonal preconditioner
!      -3   ILU preconditioner (in general the best choice)
!
      IINSOL(2) = INT(PNUMS(22))                                          30.82
!
!     IINSOL(3) control parameter for the amount of output
!     Possible values:
!             <0  : No output.
!              0  : Only fatal errors will be printed.
!              1  : Additional information about the iteration
!                   is printed.
!              2  : Gives a maximal amount of output
!                   concerning the iteration process.
!
      IINSOL(3) = INT(PNUMS(24))                                          30.82
!
!     IINSOL(6) Maximal number of iterations to be performed
!               in each of the solution methods.
!
      IINSOL(6) = INT(PNUMS(25))                                          30.82
!
!     RINSOL(2) required accuracy, the iterative method stops
!               if ||r ||  < RINSOL(2) * ||r ||
!                     k  2                  0  2
!
      RINSOL(2) = PNUMS(23)                                                   30.82
!
      CALL ISSOLV(IINSOL, RINSOL, MATRIX, RHSIDE, SETUP,
     &            NPOINT, NCONCT, INFMAT,
     &            WORK(NPREC*NPOINT+1), (NWORK-NPREC)*NPOINT,
     &            WORK,   NPREC,
     &            UPPERI, LOPERI, INOCNV,
     &            ITSW,   ITER)
      IF (STPNOW()) RETURN                                                34.01
!
      END

