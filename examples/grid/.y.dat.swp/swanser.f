!
!     SWAN - SERVICE ROUTINES
!
!  Contents of this file
!
!     READXY
!     REFIXY
!     INFRAM
!     DISTR
!     KSCIP1
!     AC2TST
!     CVCHEK                                                              30.60
!     CVMESH                                                              30.60
!     NEWTON                                                              30.60
!     EVALF                                                               30.60
!     SWOBST                                                              30.60
!     TCROSS2                                                             40.04
!     SWTRCF
!     SSHAPE                                                              40.00
!     SINTRP                                                              40.00
!     HSOBND                                                              32.01
!     CHGBAS                                                              40.00
!     GAMMA                                                               40.00
!     WRSPEC                                                              40.00
!     SWTSTA                                                              40.23
!     SWTSTO                                                              40.23
!     SWPRTI                                                              40.23
!     TXPBLA                                                              40.23
!     INTSTR                                                              40.23
!     NUMSTR                                                              40.23
!     SWCOPI                                                              40.23
!     SWCOPR                                                              40.23
!     SWI2B                                                               40.30
!     SWR2B                                                               40.30
!
!  functions:
!  ----------
!  DEGCNV  (converts from cartesian convention to nautical and            32.01
!           vice versa)                                                   32.01
!  ANGRAD  (converts radians to degrees)                                  32.01
!  ANGDEG  (converts degrees to radians)                                  32.01
!
!  subroutines:
!  ------------
!  HSIBND  (Hs is calculated at those side where a boundary condition     32.01
!           is provided and directly after the 'BOUNDARY' command         32.01
!           and subsequently stored in array HSI)                         32.01
!  HSOBND  (Hs is calculated after a SWAN computation at all sides.       32.01
!           The calculated wave height from SWAN is then compared with    32.01
!           the wave heigth as provided by the user (HSIBND)              32.01
!  SPCVAR  (A variable boundary condition is read along a side of the     32.01
!           computational grid)                                           32.01
!  SPCRD   (Read an input file with 1D- or 2D- spectra)                   32.01
!  INTSPEC (interpolate in geographical and directional space along       32.01
!           a boundary side when two- or more wave spectra are provided)  32.01
!  ROTSPEC (rotate spectra at boundary for interpolation procedure in     32.01
!           directional space)                                            32.01
!  STRSPEC (strech spectrum at boundary for interpolation procedure in    32.01
!           frequency space)                                              32.01
!  SETUPP  (Compute the wave-induced setup for a one-dimensional and a    32.01
!           two-dimensional run. Note that the one-dimensional mode of    32.01
!           swan has been coded in this project (H3268)                   32.01
!  SETUP2D (Computation of SETUP, the change of waterlevel by waves.      31.03
!           A Poisson equation is solved in general coordinates           31.03
!
!
!***********************************************************************
!                                                                      *
      SUBROUTINE READXY (NAMX, NAMY, XX, YY, KONT, XSTA, YSTA)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm1.inc'                                              40.13
      INCLUDE 'swcomm2.inc'                                               30.74
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
!     40.22: John Cazes and Tim Campbell
!     40.13: Nico Booij
!
!  1. UPDATE
!
!       Nov. 1996               offset values are added to standard values
!                               because they will be subtracted later
!     40.13, Nov. 01: a valid value for YY is required if a valid value
!                     for XX has been given; ocpcomm1.inc reactivated
!
!  2. PURPOSE
!
!       Read x and y, initialize offset values XOFFS and YOFFS
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NAMX, NAMY   inp char    names of the two coordinates as given in
!                                the user manual
!       XX, YY       out real    values of x and y taking into account offset
!       KONT         inp char    what to be done if values are missing
!                                see doc. of INDBLE (Ocean Pack doc.)
!       XSTA, YSTA   inp real    standard values of x and y
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
!       INDBLE (Ocean Pack)
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Read x and y in double prec.
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      DOUBLE PRECISION XTMP, YTMP
      CHARACTER  NAMX *(*), NAMY *(*), KONT *(*)
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'READXY')
!
      CALL INDBLE (NAMX, XTMP, KONT, DBLE(XSTA)+DBLE(XOFFS))
      IF (CHGVAL) THEN                                                    40.13
!       a valid value was given for XX                                    40.13
        CALL INDBLE (NAMY, YTMP, 'REQ', DBLE(YSTA)+DBLE(YOFFS))           40.13
      ELSE                                                                40.13
        CALL INDBLE (NAMY, YTMP, KONT, DBLE(YSTA)+DBLE(YOFFS))
      ENDIF                                                               40.13
      IF (.NOT.LXOFFS) THEN
        XOFFS = REAL(XTMP)
        YOFFS = REAL(YTMP)
        LXOFFS = .TRUE.
      ENDIF
      XX = REAL(XTMP-DBLE(XOFFS))
      YY = REAL(YTMP-DBLE(YOFFS))
!
      RETURN
! * end of subroutine READXY  *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE REFIXY (NDS, XX, YY, IERR)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'swcomm2.inc'                                               30.74
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
!     40.22: John Cazes and Tim Campbell
!
!  1. UPDATE
!
!       first version: 10.18 (Sept 1994)
!
!  2. PURPOSE
!
!       initialize offset values XOFFS and YOFFS, and shift XX and YY
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NDS          in  int     file reference number
!       XX, YY       out real    values of x and y taking into account offset
!       IERR         out int     error indicator: IERR=0: no error, =-1: end-
!                                of-file, =-2: read error
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      DOUBLE PRECISION XTMP, YTMP
      REAL             XX, YY
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'REFIXY')
!
      READ (NDS, *, END=10, ERR=20) XTMP, YTMP
      IF (.NOT.LXOFFS) THEN
        XOFFS = REAL(XTMP)
        YOFFS = REAL(YTMP)
        LXOFFS = .TRUE.
      ENDIF
      XX = REAL(XTMP-DBLE(XOFFS))
      YY = REAL(YTMP-DBLE(YOFFS))
!
      IERR = 0
      RETURN
!     end of file
  10  IERR = -1
      RETURN
!     read error
  20  IERR = -2
      RETURN
! * end of subroutine REFIXY  *
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION  INFRAM (XQQ, YQQ)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm1.inc'                                               30.74
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
!  1. UPDATE
!
!     40.22: John Cazes and Tim Campbell
!
!  2. PURPOSE
!
!       Checking whether a point given in frame coordinates is located
!       in the plotting frame (INFRAM = .TRUE.) or not (INFRAM = .FALSE.)
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       XQQ     REAL   input    X-coordinate (output grid) of the point
!       YQQ     REAL   input    Y-coordinate (output grid) of the point
!
!  5. SUBROUTINES CALLING
!
!       SPLSIT, PLNAME
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Give INFRAM initial value true
!       IF XQQ < 0, XQQ > XQLEN, YQQ < 0 OR YQQ > YQLEN, THEN
!           INFRAM = false
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'INFRAM')
!
      INFRAM = .TRUE.
      IF (XQQ .LT.    0.) INFRAM = .FALSE.
      IF (XQQ .GT. XQLEN) INFRAM = .FALSE.
      IF (YQQ .LT.    0.) INFRAM = .FALSE.
      IF (YQQ .GT. YQLEN) INFRAM = .FALSE.
!
      RETURN
! * end of function INFRAM *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE DISTR (CDIR, DIR, COEF, SPCDIR)                          20.43
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!
!  0. Authors
!
!     30.82: IJsbrand Haagsma
!     40.22: John Cazes and Tim Campbell
!
!  1. Updates
!
!      0.1 , Jul. 87: Standard heading added
!      0.2 , Dec. 89: Value for energy outside of the sector changed
!                     from 0. to 1.E-6
!            Oct. 90: Value for energy outside the sector changed to 1.E-10
!                     logical BDIR introduced to take care for case where
!                     none of the values is positive
!     30.82, Oct. 98: Updated description of several variables
!
!  2. PURPOSE
!
!       Computation of the distribution of the wave energy over the
!       sectors, according to the given directional spread.
!
!  3. METHOD
!
!       ---
!
!  4. Argument variables
!
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
!
      REAL    SPCDIR(MDC,6)                                               30.82
!
!       CDIR    REAL   output   array containing the coefficients of
!                               (energy) distribution
!       DIR     REAL   input    main wave direction, in radians
!       COEF    REAL   input    coefficient of the directional distri-
!                               bution (cos**COEF)
!
!  5. SUBROUTINES CALLING
!
!       REINVA (HISWA/SWREAD), STARTB and SWIND (both HISWA/COMPU)
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       For every direction of the grid do
!           If the direction deviates less than PI/2 from the main wave
!            direction, then
!               Compute the coefficient cos**n
!           Else
!               Coefficient is 1.E-10
!       -----------------------------------------------------------------
!       If any of the directions deviated less than PI/2
!       Then Compute the total of the coefficients
!            For every direction of the grid do
!                Divide the fraction of the distribution by the total
!       -----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL BDIR
      REAL CDIR(*)
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT, 'DISTR')
!
      SOMC = 0.
      BDIR = .FALSE.
      DO 10 ID0 = 1, MDC
        TETA = SPCDIR(ID0,1)
        ACOS = COS(TETA-DIR)
        IF (ACOS .GT. 0.) THEN
          BDIR = .TRUE.
          CDIR(ID0) = MAX (ACOS**COEF, 1.E-10)
          SOMC = SOMC + CDIR(ID0)
        ELSE
          CDIR(ID0) = 1.E-10
        ENDIF
  10  CONTINUE
      IF (BDIR) THEN
        CNORM = 1./(SOMC*DDIR)
        DO 20 ID0 = 1, MDC
          CDIR(ID0) = CDIR(ID0) * CNORM
  20    CONTINUE
      ENDIF
!
!     ***** test *****
!      IF(TESTFL .AND. ITEST .GE. 200)
       IF( ITEST .GE. 200)
     &WRITE(PRINTF,6010) CNORM,(CDIR(JJ) , JJ=1,MDC)
 6010   FORMAT (' Test DISTR',F10.3/(10E12.4))
!
      RETURN
! * end of subroutine DISTR *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE KSCIP1 (MMT, SIG, D, K, CG, N, ND)
!                                                                      *
!***********************************************************************
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
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
!
!  0. Authors
!
!     30.81:  Annette Kieftenburg
!     40.22: John Cazes and Tim Campbell
!
!  1. Updates
!
!     Aug. 94, ver. 10.10: arguments N and ND added
!     Dec. 98, ND corrected, argument list adjusted and IMPLICIT NONE added
!
!  2. Purpose
!
!     Interpolation of the wave number, group velocity and N from a
!     table, and calculation of the derivative of N w.r.t. depth (=ND)
!
!  3. Method
!
!     --
!
!  4. Argument variables
!
!     MMT     input    number of points in freq. for which
!
      INTEGER   MMT
!
!     CG      output   group velocity
!     D       input    local depth
!     K       output   wave number
!     N       output   ratio of group and phase velocity
!     ND      output   derivative of N with respect to D
!                      computation must be done
!     SIG     input    rel. frequency for which wave parameters
!                      must be determined
!
      REAL      CG(MMT), D,
     &          K(MMT), N(MMT), ND(MMT), SIG(MMT)
!
!  5. Parameter variables
!
!     NWMAX
!     NWMIN
!
      INTEGER   NWMAX, NWMIN
      PARAMETER (NWMIN = 0, NWMAX = 25)
!
!     DWND
!     RPDW
!
      REAL      DSND, RPDW
      PARAMETER (DSND = 0.1, RPDW=1./DSND)
!
!  6. Local variables
!
!     CGND      dimensionless group velocity
!     CGTB1D    coefficients for calculating group velocity
!     DIFN      dummy variable
!     FAC       factor used for interpolation of tables
!     IENT      number of entries
!     INPW      integer part of WND*RPDW (used for coefficient tables)
!     IS        counter in frequency (sigma-space)
!     KND       dimensionless wave number
!     KTAB1D    coefficients for calculating wave number
!     NTAB1D    coefficients for calculating ratio of group and phase
!               velocity
!     ROOTDG    square root of D/GRAV
!     WFAC      =WND*RPDW
!     WGD       square root of GRAV*D
!     SND       dimensionless frequency
!     WND2      = WND*WND
!
      INTEGER   IENT, INPW, IS
      REAL      CGND, CGTB1D(NWMIN:NWMAX),
     &          DIFN,                                                      30.81
     &          FAC, KND, KTAB1D(NWMIN:NWMAX),
     &          NTAB1D(NWMIN:NWMAX), ROOTDG, SND, WFAC, WGD,
     &          SND2                                                       30.81
!
!  7. Common Blocks used
!
!     --
!
!  8. Subroutines used
!
!     --
!
!  9. Subroutines calling
!
!     SWOEXA, SWOEXF (Swan/Output)
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
!     -----------------------------------------------------------------
!      Compute non-dimensional frequency WND
!      IF WND >= 2.5, then
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND according to
!        deep water theory
!      ELSE IF WND =< 1.e-6
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND
!        according to extremely shallow water
!      ELSE
!        Compute wave number K, group velocity CGO and the ratio of
!        group and phase velocity N by interpolation from 1-dimensio-
!        nal tables. Compute the derivative of N w.r.t. D = ND.
!     -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT /0/
      DATA KTAB1D/0., 0.10016680, 0.20134288, 0.30457264, 0.41096926,
     &    0.52175194, 0.63846463, 0.76255310, 0.89595461, 1.04083824,
     &    1.19967842, 1.37519836, 1.57023048, 1.78743458, 2.02877331,
     &    2.29598522, 2.58901024, 2.90728855, 3.24976063, 3.61523247,
     &    4.00267029, 4.41130257, 4.84030056, 5.29013348, 5.76005936,
     &    6.25002289/
      DATA CGTB1D/1., 0.99501032, 0.98015553, 0.95579094, 0.92251045,
     &    0.88116169, 0.83277500, 0.77876240, 0.72069371, 0.66038823,
     &    0.59983897, 0.54110348, 0.48613495, 0.43655658, 0.39344674,
     &    0.35705268, 0.32704252, 0.30251282, 0.28235823, 0.26552880,
     &    0.25116783, 0.23864383, 0.22752625, 0.21749985, 0.20837778,
     &    0.20001733/
      DATA NTAB1D/1., 0.9966699 , 0.9867362 , 0.9703590 , 0.9478084 ,
     &    0.9194954 , 0.8861622 , 0.8483538 , 0.8071360 , 0.7637302 ,
     &    0.7196137 , 0.6764768 , 0.6361198 , 0.6002431 , 0.5701529 ,
     &    0.5465249 , 0.5291977 , 0.5173482 , 0.5097758 , 0.5052359 ,
     &    0.5026709 , 0.5013000 , 0.5005887 , 0.5002621 , 0.5001116 ,
     &    0.500045/
      IF (LTRACE) CALL STRACE (IENT, 'KSCIP1')
!
      ROOTDG = SQRT(D/GRAV)                                               30.81
      WGD    = ROOTDG*GRAV                                                30.81
      DO 200 IS = 1, MMT
!       WND is dimensionless frequency
        SND = SIG(IS) * ROOTDG
        IF (SND .GE. 2.5) THEN
!     ******* deep water *******
          K(IS)  = SIG(IS) * SIG(IS) / GRAV                                   30.81
          CG(IS) = 0.5 * GRAV / SIG(IS)                                     30.81
          N(IS)  = 0.5
          ND(IS) = 0.
        ELSE IF (SND.LT.1.E-6) THEN
!     *** very shallow water ***                                          30.81
          K(IS)  = SND/D                                                  30.81
          CG(IS) = WGD
          N(IS)  = 1.
          ND(IS) = 0.
        ELSE
          WFAC = SND * RPDW
          INPW = INT(WFAC)
          FAC  = WFAC - FLOAT(INPW)
            KND  = (1.-FAC)*KTAB1D(INPW) + FAC*KTAB1D(INPW+1)
            CGND = (1.-FAC)*CGTB1D(INPW) + FAC*CGTB1D(INPW+1)
            N(IS)  = (1.-FAC)*NTAB1D(INPW) + FAC*NTAB1D(INPW+1)
          K(IS)  = KND/D                                                  30.81
          CG(IS) = CGND * WGD
!
!         Analytical solution:                                            30.81
!           ND(IS) = K(IS)/SINH(2*K(IS)*D)*(1 - KD2/TANH(2*K(IS)*D))      30.81
!                  = (d N / d WND) * (d WND / d D)                        30.81
!                                                                         30.81
          SND2 = SND*SND                                                  30.81
            DIFN = (NTAB1D(INPW+1) - NTAB1D(INPW)) / DSND                 30.81
          ND(IS) = DIFN*0.5*(D*K(IS)*K(IS) + SND2*(1.-SND2)/D)/SND        30.81
        ENDIF
  200 CONTINUE
!
      RETURN
!     end of subroutine KSCIP1 *
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE AC2TST (XYTST, AC2,KGRPNT)
!                                                                      *
!***********************************************************************
!
      USE M_PARALL                                                        40.31
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'                                               30.74
!
!     0. Authors
!
!
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
!
!
      INTEGER   XYTST(*) ,KGRPNT(MXC,MYC)
      REAL      AC2(MDC,MSC,MCGRD)                                        30.21
!.................................................................
      IF ( ITEST .GE. 100 .AND. TESTFL) THEN
        DO II = 1, NPTST
          IX = XYTST(2*II-1)
          IY = XYTST(2*II)
          INDEX = KGRPNT(IX,IY)
          WRITE (PRINTF, 618) IX+MXF-2, IY+MYF-2, KGRPNT(IX,IY)           40.30
 618      FORMAT(/,'Spectrum for test point(index):', 2I5,2X,'(',I5,')')
          DO ID = 1, MDC
            WRITE (PRINTF, 619) (AC2(ID,IS,INDEX), IS=1,MIN(10,MSC))      30.21
 619        FORMAT (10(1X,E12.4))
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
!****************************************************************
!
      SUBROUTINE CVCHEK (KGRPNT, XCGRID, YCGRID)                          30.72
!
!****************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm2.inc'                                               30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'                                               30.74
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.13: Nico Booij
!
!  1. Updates
!
!            May  96: New subroutine
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.13, Mar. 01: messages corrected and extended
!
!  2. Purpose
!
!     Checks whether the given curvilinear grid is correct
!     also set the value of CVLEFT.
!
!  3. Method
!
!     Going around a mesh in the same direction the interior
!     of the mesh must be always in the same side if the
!     coordinates are correct
!
!  4. Argument variables
!
!     KGRPNT: input  Array of indirect addressing
!
      INTEGER KGRPNT(MXC,MYC)                                             30.72
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!
!     5. SUBROUTINES CALLING
!
!        SWRBC
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!     FIRST = True
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             Inmesh = True
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy),
!                    K3 = KGRPNT(ix+1,iy+1)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1),
!                    K3 = KGRPNT(ix,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1),
!                    K3 = KGRPNT(ix,iy)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy),
!                    K3 = KGRPNT(ix+1,iy)
!                 ---------------------------------------------------
!                 If K1>1 and K2>1 and K3>1
!                 Then Det = (xpg(K3)-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (ypg(K3)-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If FIRST
!                      Then Make FIRST = False
!                           If Det>0
!                           Then Make CVleft = False
!                           Else Make CVleft = True
!                      ----------------------------------------------
!                      If ((CVleft and Det<0) or (not CVleft and Det>0))
!                      Then Write error message with IX, IY, ISIDE
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!****************************************************************
!
!
      LOGICAL  FIRST
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'CVCHEK')
!
!     test output
!
      IF (ITEST .GE. 150 .OR. INTES .GE. 30) THEN
        WRITE(PRINTF,186)
 186    FORMAT(/,' ... Subroutine CVCHEK...',
     &  /,2X,'POINT( IX, IY),  INDEX,      COORDX,       COORDY')
        ICON = 0
        DO 5 IIY = 1, MYC
          DO 6 IIX = 1, MXC
            ICON = ICON + 1
            WRITE(PRINTF,7)IIX-1,IIY-1,KGRPNT(IIX,IIY),
     &      XCGRID(IIX,IIY)+XOFFS, YCGRID(IIX,IIY)+YOFFS                  30.72 40.13
 6        CONTINUE
 5      CONTINUE
      ENDIF
 7    FORMAT(4X,I5,1X,I5,3X,I4,5X,F10.2,4X,F10.2)
!
      FIRST = .TRUE.
!
      DO 10 IX = 1,MXC-1
        DO 15 IY = 1,MYC-1
          DO 20 ISIDE = 1,4
            IF (ISIDE .EQ. 1) THEN
              IX1 = IX                                                    40.13
              IY1 = IY                                                    40.13
              IX2 = IX+1                                                  40.13
              IY2 = IY                                                    40.13
              IX3 = IX+1                                                  40.13
              IY3 = IY+1                                                  40.13
            ELSE IF (ISIDE .EQ. 2) THEN
              IX1 = IX+1                                                  40.13
              IY1 = IY                                                    40.13
              IX2 = IX+1                                                  40.13
              IY2 = IY+1                                                  40.13
              IX3 = IX                                                    40.13
              IY3 = IY+1                                                  40.13
            ELSE IF (ISIDE .EQ. 3) THEN
              IX1 = IX+1                                                  40.13
              IY1 = IY+1                                                  40.13
              IX2 = IX                                                    40.13
              IY2 = IY+1                                                  40.13
              IX3 = IX                                                    40.13
              IY3 = IY                                                    40.13
            ELSE IF (ISIDE .EQ. 4) THEN
              IX1 = IX                                                    40.13
              IY1 = IY+1                                                  40.13
              IX2 = IX                                                    40.13
              IY2 = IY                                                    40.13
              IX3 = IX+1                                                  40.13
              IY3 = IY                                                    40.13
            ENDIF
            K1  = KGRPNT(IX1,IY1)                                         40.13
            XC1 = XCGRID(IX1,IY1)                                         40.13 30.72
            YC1 = YCGRID(IX1,IY1)                                         40.13 30.72
            K2  = KGRPNT(IX2,IY2)                                         40.13
            XC2 = XCGRID(IX2,IY2)                                         40.13 30.72
            YC2 = YCGRID(IX2,IY2)                                         40.13 30.72
            K3  = KGRPNT(IX3,IY3)                                         40.13
            XC3 = XCGRID(IX3,IY3)                                         30.72
            YC3 = YCGRID(IX3,IY3)                                         30.72
            DET   = 0.
            IF (K1 .GE. 2 .AND. K2 .GE. 2 .AND. K3 .GE. 2) THEN
              DET = ((XC3 - XC1) * (YC2 - YC1)) -
     &              ((YC3 - YC1) * (XC2 - XC1))
              IF (DET .EQ. 0.) THEN
!               three grid points on one line                             40.13
                CALL MSGERR (2,'3 comp. grid points on one line')         40.13
                WRITE (PRINTF, 112)
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,                  40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,                  40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                   40.13
 112            FORMAT (3(1X, 2I3, 2(1X, F14.4)))                         40.13
              ENDIF
!
              IF (FIRST) THEN
                FIRST = .FALSE.
                IF (DET .GT. 0.) THEN
                  CVLEFT = .FALSE.
                ELSE
                  CVLEFT = .TRUE.
                ENDIF
              ENDIF
              IF (     (      CVLEFT .AND. DET .GT. 0.)
     &            .OR. (.NOT. CVLEFT .AND. DET .LT. 0.)) THEN
!               crossing grid lines in a mesh                             40.13
                CALL MSGERR (2,'Grid angle <0 or >180 degrees')           40.13
                WRITE (PRINTF, 112)
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,                  40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,                  40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                   40.13
              ENDIF
            ENDIF
 20       CONTINUE
 15     CONTINUE
 10   CONTINUE
      RETURN
!     *** end of subroutine CVCHEK ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID ,YCGRID,
     &                   KGRBND)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm2.inc'                                               30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'                                               30.74
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
!     30.72: IJsbrand Haagsma
!     40.00, 40.13: Nico Booij
!     40.02: IJsbrand Haagsma
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, May  98: procedure for points outside grid accelerated
!     40.00, Feb  99: procedure extended for 1D case
!                     XOFFS and YOFFS added in write statements
!     40.02, Mar. 00: Fixed bug that placed dry testpoints outside computational grid
!     40.13, Mar. 01: message "CVMESH 2nd attempt .." suppressed
!
!  2. Purpose
!
!     procedure to find location in curvilinear grid for a point
!     given in problem coordinates
!
!  3. Method
!
!     First attempt: Use Newton-Raphson method to find XC and YC
!     (Note: in the program XC and YC indicate the mesh and position in
!     the mesh) in a few steps; this may be most efficient if a series of
!     points is processed, because the previous point provides a good
!     first estimate.
!
!     The procedure may fail for one reason:
!     a) the number of iterations is larger than a previously set limit,
!        say 10.
!
!     If the Newton-Raphson procedure fails, a second attempt:
!     Scan all meshes of the grid to find whether XP,YP is inside the mesh.
!     A point is assumed to be inside a mesh if it is on the
!     interior side for all 4 sides of the mesh. Here we use common variable
!     CVLEFT (logical): if True interior of the mesh is always on the left
!     going along the mesh in the order: (ix,iy), (ix+1,iy), (ix+1,iy+1),
!     (ix,iy+1), (ix,iy). If it is False the interior is always on the right.
!     Whether a point  is on the left or on the right of a line from  to
!     can be decided by looking at the sign of the determinant.
!       If the point is inside of any mesh of the computational grid then
!     newton-raphson procedure is used again with the pivoting point like
!     first guess.
!
!  4. Argument variables
!
!     XCGRID  input  Coordinates of computational grid in x-direction     30.72
!     YCGRID  input  Coordinates of computational grid in y-direction     30.72
!     XP, YP  input  a point given in problem coordinates
!     XC, YC  outp   same point in computational grid coordinates
!
      REAL     XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                        30.72
      REAL     XP, YP, XC, YC
!
!     KGRPNT   input   array(MXC,MYC)  grid numbers
!                      if KGRPNT <= 1, point is not in comp. grid.
!     KGRBND   input   lists all boundary grid points consecutively
!
      INTEGER  KGRPNT(MXC,MYC), KGRBND(*)                                 40.00
!
!     Local variables
!
!     MXITNR   number of iterations in Newton-Raphson procedure
!     IX, IY   counter of computational grid point
!     K1       address of grid point
!     IXMIN    counter of grid point closest to (XP,YP)
!     IYMIN    counter of grid point closest to (XP,YP)
!     IBND     counter of boundary grid points
!
      INTEGER       :: IX, IY, K1, IXMIN, IYMIN, IBND
      INTEGER, SAVE :: MXITNR = 0
      INTEGER, SAVE :: IENT = 0
!
!     INMESH   if True, point (XP,YP) is inside the computational grid
!     FINDXY   if True, Newton-Raphson propcedure succeeded
!
      LOGICAL  INMESH ,FINDXY
!
!     XCSAVE   XC computed for previous point, used as first guess
!     YCSAVE   YC computed for previous point, used as first guess
!     XPC1     user coordinate of a computational grid point
!     YPC1     user coordinate of a computational grid point
!     XC0      grid coordinate of grid point closest to (XP,YP)
!     YC0      grid coordinate of grid point closest to (XP,YP)
!
      REAL       :: XPC1, YPC1, XC0, YC0
      REAL, SAVE :: XCSAVE=0., YCSAVE=0.                                  40.31
!
!  5. SUBROUTINES CALLING
!
!     SINCMP
!
!  6. SUBROUTINES USED
!
!       NEWTON
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       XCSAVE and YCSAVE are used as first guess of XC and YC
!       values are declared SAVE
!
!  9. STRUCTURE
!
!     --------------------------------------------------------------
!     Make XC=XCSAVE and YC=YCSAVE
!     Determine XC and YC FROM XP AND YP using a Newton-Raphson iteration
!     If (XC and YC were found) then
!       Procedure is ready; Return values of XC and YC
!       Make XCSAVE=XC and YCSAVE=YC
!       return
!     else
!     ---------------------------------------------------------------------
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             Inmesh = True
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy)
!                 ----------------------------------------------------------
!                 If K1>0 and K2>0
!                 Then Det = (xp-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (yp-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If ((CVleft and Det>0) or (not CVleft and Det<0))
!                      Then Make Inmesh = False
!                      Else  Inmesh = true and XC = IX and YC = IY
!                 Else Make Inmesh = False
!             --------------------------------------------------------
!             If Inmesh
!             Then Determine XC and YC using one Newton-Raphson iteration
!                  step
!                  Procedure is ready; Return values of XC and YC
!                  Make XCSAVE=XC and YCSAVE=YC
!     ---------------------------------------------------------------------
!     No mesh is found: Make XC and YC = exception value for XC
!     Return values of XC and YC
!     ---------------------------------------------------------------------
!
!****************************************************************
!
!
      IF (LTRACE) CALL STRACE (IENT,'CVMESH')
!
      IF (ONED) THEN
        CALL NEWT1D  (XP, YP, XCGRID, YCGRID, KGRPNT,                     40.00
     &                MXITNR ,XC ,YC ,FINDXY)
        IF (.NOT.FINDXY) THEN
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                          40.00
          ENDIF
        ENDIF
      ELSE
!       two-dimensional computation
        IF (XCSAVE .LE. MXC .AND. XCSAVE .GE. 0. .OR.
     &      YCSAVE .LE. MYC .AND. YCSAVE .GE. 0.) THEN
          XC     = XCSAVE
          YC     = YCSAVE
!         *** First attempt, to find XC ,YC with Newton-Raphson method**
          MXITNR = 5
          CALL NEWTON  (XP, YP, XCGRID, YCGRID, KGRPNT,                   40.00
     &                  MXITNR ,XC ,YC ,FINDXY, KGRBND)                   40.02
          IF ((ITEST .GE. 150 .OR. INTES .GE. 20) .AND. FINDXY) THEN      40.02
            WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC                   40.03
          ENDIF
 25       FORMAT (' CVMESH: (XP,YP)=','(',E9.2,',',E9.2,
     &            '), (XC,YC)=','(',F9.2,',',F9.2,')')
!
          IF (FINDXY) GO TO 80
        ENDIF
!
        IF (INMESH (XP, YP, XCGRID ,YCGRID, KGRBND)) THEN
!         select grid point closest to (XP,YP)
          DISMIN = 1.E20
          DO 50 IX = 1,MXC
            DO 40 IY = 1,MYC
              K1  = KGRPNT(IX,IY)
              IF (K1.GT.1) THEN
                XPC1 = XCGRID(IX,IY)
                YPC1 = YCGRID(IX,IY)
                DISXY = SQRT ((XP-XPC1)**2 + (YP-YPC1)**2)
                IF (DISXY .LT. DISMIN) THEN
                  IXMIN  = IX
                  IYMIN  = IY
                  DISMIN = DISXY
                ENDIF
              ENDIF
  40        CONTINUE
  50      CONTINUE
!         second attempt using closest grid point as first guess
          MXITNR = 20
          XC0 = REAL(IXMIN)
          YC0 = REAL(IYMIN)
!         ITEST condition changed from 20 to 120                          40.13
          IF (ITEST.GE.120) WRITE (PRTEST, 55) XP+XOFFS ,YP+YOFFS ,       40.13
     &          XC0-1. ,YC0-1.
  55      FORMAT (' CVMESH 2nd attempt, (XP,YP)=','(',E9.2,',',E9.2,
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
          DO KORNER = 1, 4
            IF (KORNER.EQ.1) THEN
              XC = XC0 + 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.2) THEN
              XC = XC0 - 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.3) THEN
              XC = XC0 - 0.2
              YC = YC0 - 0.2
            ELSE
              XC = XC0 + 0.2
              YC = YC0 - 0.2
            ENDIF
            CALL NEWTON  (XP, YP, XCGRID, YCGRID, KGRPNT,                 40.00
     &                    MXITNR ,XC ,YC ,FINDXY, KGRBND)                 40.02
            IF (FINDXY) THEN
              IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC               40.00
              ENDIF
              GOTO 80
            ENDIF
          ENDDO
          WRITE (PRINTF, 75) XP+XOFFS, YP+YOFFS, MXITNR                   40.00
  75      FORMAT (' search of grid coordinates fails for:', 2F10.2,
     &            ' in', I3, ' iterations')
        ELSE
!         scan boundary to see whether the point is close to the boundary
          IX2 = 0
          DO IBND = 1, NGRBND
            IX1 = IX2
            IY1 = IY2
            XP1 = XP2
            YP1 = YP2
            IX2 = KGRBND(2*IBND-1)
            IY2 = KGRBND(2*IBND)
            IF (IX2.GT.0) THEN
              XP2 = XCGRID(IX2,IY2)
              YP2 = YCGRID(IX2,IY2)
              IF (IX1.GT.0) THEN
!               determine relative distance from boundary section
                SLEN2  = (XP2-XP1)**2 + (YP2-YP1)**2
                RELDIS = ((XP-XP1)*(YP2-YP1)-(YP-YP1)*(XP2-XP1)) /
     &                   SLEN2
                IF (ABS(RELDIS).LT.0.2) THEN
!                 determine location on the boundary section
                  RELLOC = ((XP-XP1)*(XP2-XP1)+(YP-YP1)*(YP2-YP1)) /
     &                      SLEN2
                  IF (RELLOC.GE.0. .AND. RELLOC.LE.1.) THEN
                    XC = FLOAT(IX1) + RELLOC * FLOAT(IX2-IX1) - 1.
                    YC = FLOAT(IY1) + RELLOC * FLOAT(IY2-IY1) - 1.
                    IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                      WRITE(PRINTF, 65) XP+XOFFS, YP+YOFFS, XC, YC        40.00
  65                  FORMAT (' CVMESH: (XP,YP)=','(',E9.2,',',E9.2,
     &                        ') is on the boundary, (XC,YC)=(',
     &                        F9.2,',',F9.2,')')
                    ENDIF
                    GOTO 80
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                          40.00
  85        FORMAT (' CVMESH: (XP,YP)=','(',E9.2,',',E9.2,
     &              ') is outside grid')
          ENDIF
        ENDIF
      ENDIF                                                               40.00
  80  XCSAVE = XC
      YCSAVE = YC
      RETURN
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION INMESH (XP, YP, XCGRID ,YCGRID, KGRBND)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'swcomm2.inc'                                               40.03
      INCLUDE 'swcomm3.inc'
      INCLUDE 'ocpcomm4.inc'                                              40.03
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
!     Nico Booij
!
!  1. Updates
!
!       New function for curvilinear version (ver. 40.00). May '98
!       40.03, Dec 99: test output added; commons swcomm2 and ocpcomm4 added
!
!  2. Purpose
!
!       procedure to find whether a given location is
!       in the (curvilinear) computational grid
!
!  3. Method  suggested by Gerbrant van Vledder
!
!       draw a line from the point (XP,YP) in vertical direction
!       determine the number of crossings with the boundary of the
!       grid; if this number is even the point is outside
!
!  4. Argument variables
!
!
!     KGRBND   int  input   array containing boundary grid points
!
      INTEGER  KGRBND(*)
!
!     XP, YP    real, input   a point given in problem coordinates
!     XCGRID    real, input   array(IX,IY) x-coordinate of a grid point
!     YCGRID    real, input   array(IX,IY) y-coordinate of a grid point
!
      REAL     XCGRID(MXC,MYC) ,YCGRID(MXC,MYC),
     &         XP, YP
!
!  5. Parameter variables
!
!  6. Local variables
!
!     NUMCRS   number of crossings with boundary outline
!
      INTEGER  NUMCRS, IX1, IX2, IY2
      REAL     XP1, XP2, YP1, YP2, YPS
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!       CVMESH
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     --------------------------------------------------------------
!     numcros = 0
!     For all sections of the boundary do
!         determine coordinates of end points (XP1,YP1) and (XP2,YP2)
!         If (XP1<XP and XP2>XP) or (XP1>XP and XP2<XP)
!         then If not (YP1<YP and YP2<YP)
!                   if YPS>YP
!                   then numcros = numcros + 1
!     ---------------------------------------------------------------
!     If numcros is even
!     Then Inmesh = False
!     Else Inmesh = True
!     ---------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE (IENT,'INMESH')
!
      IF (XP.LT.XCGMIN .OR. XP.GT.XCGMAX .OR.
     &    YP.LT.YCGMIN .OR. YP.GT.YCGMAX) THEN
        IF (ITEST.GE.70) WRITE (PRTEST, 22) XP+XOFFS, YP+YOFFS,           40.03
     &    XCGMIN+XOFFS, XCGMAX+XOFFS, YCGMIN+YOFFS, YCGMAX+YOFFS
  22    FORMAT (1X, 2E12.4, ' is outside region ', 4E12.4)
        INMESH = .FALSE.
        GOTO 90
      ENDIF
!
      IF (NGRBND.LE.0) THEN
        CALL MSGERR (3, 'grid outline not yet determined')
        RETURN
      ENDIF
!
      NUMCRS = 0
      IX2 = 0
!     loop over the boundary of the computational grid
      DO IBND = 1, NGRBND
        IX1 = IX2
        XP1 = XP2
        YP1 = YP2
        IX2 = KGRBND(2*IBND-1)
        IY2 = KGRBND(2*IBND)
        IF (IX2.GT.0) THEN
          XP2 = XCGRID(IX2,IY2)
          YP2 = YCGRID(IX2,IY2)
          IF (ITEST.GE.180) WRITE (PRTEST, 28) XP2+XOFFS,                 40.03
     &    YP2+YOFFS
  28      FORMAT (' boundary point ', 2E12.4)
          IF (IX1.GT.0) THEN
            IF (((XP1.GT.XP).AND.(XP2.LT.XP)).OR.
     &          ((XP1.LT.XP).AND.(XP2.GT.XP))) THEN
              IF (YP1.GT.YP .OR. YP2.GT.YP) THEN
!               determine y-coordinate of crossing point
                YPS = YP1 + (XP-XP1) * (YP2-YP1) / (XP2-XP1)
                IF (YPS.GT.YP) THEN
                  NUMCRS = NUMCRS + 1
                  IF (ITEST.GE.70) WRITE (PRTEST, 32) NUMCRS,             40.03
     &            XP+XOFFS, YP+YOFFS, YPS+YOFFS                           40.03
  32              FORMAT (' crossing ', I1, ' point ', 3E12.4)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!     point is inside the grid is number of crossings is odd
      IF (MOD(NUMCRS,2) .EQ. 1) THEN
        INMESH = .TRUE.
      ELSE
        INMESH = .FALSE.
      ENDIF
  90  RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWTON (XP, YP, XCGRID, YCGRID, KGRPNT,                  40.00
     &                   MXITNR, XC, YC, FIND, KGRBND)                    40.02
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'                                               30.74
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     30.80, Oct. 98: computation of update of XC,YC modified to avoid
!                     division by 0
!
!  2. Purpose
!
!     Solve eqs. and find a point  (XC,YC) in a curvilinear grid (compt.
!     grid) for a given point (XP ,YP) in a cartesian grid (problem coord).
!
!  3. Method
!
!     In this subroutine the next equations are solved :
!
!                  @XP             @XP
!     XP(xc,yc) +  --- * @XC   +   --- * @YC  - XP(x,y) = 0
!                  @XC             @YC
!
!                  @YP             @YP
!     YP(xc,yc) +  --- * @XC   +    --- * @YC  - YP(x,y) = 0
!                  @XC             @YC
!
!     In the subroutine, next notation is used for the previous eqs.
!     XVC       + DXDXC * DXC   + DXDYC * DYC - XP  = 0.
!     YVC       + DYDXC * DXC   + DYDYC * DYC - YP  = 0.
!
!
!  4. Argument variables
!
! i   KGRBND: Grid adresses of the boundary points                        40.02
! i   KGRPNT: Grid adresses                                               40.00
! i   MXITNR: Maximum number of iterations                                30.82
!
      INTEGER KGRBND(*), KGRPNT(MXC,MYC), MXITNR                          40.02
!
!   o XC    : X-coordinate in computational coordinates                   30.82
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   XP    : X-coordinate in problem coordinates                         30.82
!   o YC    : Y-coordinate in computational coordinates                   30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
! i   YP    : Y-coordinate in problem coordinates                         30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                     30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                     30.82
!
!   o FIND  : Whether XC and YC are found                                 30.82
!
      LOGICAL FIND                                                        30.82
!
!  6. SUBROUTINES USED
!
!     STRACE
!
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0
      REAL, SAVE    :: TOLDC = 0.001
      IF (LTRACE) CALL STRACE (IENT,'NEWTON')
!
      DXC    = 1000.
      DYC    = 1000.
      TOLDC  = 0.001
      FIND   = .FALSE.
!
      IF (ITEST .GE. 200) THEN
        WRITE(PRINTF,*) ' Coordinates in subroutine NEWTON '
        DO J = 1, MYC
          DO I = 1, MXC
            WRITE(PRINTF,30) I ,J ,XCGRID(I,J) ,YCGRID(I,J)               30.72
          ENDDO
        ENDDO
      ENDIF
 30   FORMAT(2(2X,I5),2(2X,E12.4))
!
      DO 14 K = 1 ,MXITNR
!       *** If the guess point (XC,YC) is outside of compt. ***
!       *** grid, put that point in the closest boundary    ***
        IF (XC .LT. 1. ) XC = 1.
        IF (YC .LT. 1. ) YC = 1.
        IF (XC .GT. MXC) XC = FLOAT(MXC)
        IF (YC .GT. MYC) YC = FLOAT(MYC)
        I1   = INT(XC)                                                    40.00
        J1   = INT(YC)
        IF (I1 .EQ. MXC) I1 = I1 - 1
        IF (J1 .EQ. MYC) J1 = J1 - 1
        I2  = I1 + 1
        J2  = J1 + 1
        FJ1 = FLOAT(J1)
        FI1 = FLOAT(I1)
        FJ2 = FLOAT(J2)
        FI2 = FLOAT(I2)
        IF (KGRPNT(I1,J1).LE.1 .OR. KGRPNT(I2,J1).LE.1 .OR.               40.00
     &      KGRPNT(I1,J2).LE.1 .OR. KGRPNT(I2,J2).LE.1) THEN              40.00
          FIND = .FALSE.
          RETURN                                                          40.00
        ENDIF
!
        XVC   = (YC-FJ1)*((XC-FI1)*XCGRID(I2,J2)  +
     &                    (FI2-XC)*XCGRID(I1,J2)) +
     &          (FJ2-YC)*((XC-FI1)*XCGRID(I2,J1)  +
     &                    (FI2-XC)*XCGRID(I1,J1))
        YVC   = (YC-FJ1)*((XC-FI1)*YCGRID(I2,J2)  +
     &                    (FI2-XC)*YCGRID(I1,J2)) +
     &          (FJ2-YC)*((XC-FI1)*YCGRID(I2,J1)  +
     &                    (FI2-XC)*YCGRID(I1,J1))
        DXDXC = (YC -FJ1)*(XCGRID(I2,J2) - XCGRID(I1,J2)) +
     &          (FJ2-YC )*(XCGRID(I2,J1) - XCGRID(I1,J1))
        DXDYC = (XC -FI1)*(XCGRID(I2,J2) - XCGRID(I2,J1)) +
     &          (FI2-XC )*(XCGRID(I1,J2) - XCGRID(I1,J1))
        DYDXC = (YC -FJ1)*(YCGRID(I2,J2) - YCGRID(I1,J2)) +
     &          (FJ2-YC )*(YCGRID(I2,J1) - YCGRID(I1,J1))
        DYDYC = (XC -FI1)*(YCGRID(I2,J2) - YCGRID(I2,J1)) +
     &          (FI2-XC )*(YCGRID(I1,J2) - YCGRID(I1,J1))
!
        IF (ITEST .GE. 150)
     &    WRITE(PRINTF,35) K, XC-1., YC-1., XP, YP, XVC, YVC              40.00
 35     FORMAT(' NEWTON  iter=', I2, ' (XC,YC)=', 2(1X,F10.2),/,          40.00
     &         ' (XP,YP)=', 2(1X,F10.2),
     &         '  X,Y(XC,YC) = ', 2(1X,F10.2))
        IF (ITEST .GE. 180) WRITE(PRINTF,36)
     &     XCGRID(I1,J1), XCGRID(I1,J2), XCGRID(I2,J1), XCGRID(I2,J2),
     &     YCGRID(I1,J1), YCGRID(I1,J2), YCGRID(I2,J1), YCGRID(I2,J2),
     &                     DXDXC, DXDYC, DYDXC, DYDYC                     40.00
 36     FORMAT(' NEWTON grid coord:', 8(1x, F10.0), /
     &         '        deriv=', 4(1X,F10.2))                             40.00
!
!       *** If the accuracy is reached stop the iteration,  ***
        IF (ABS(DXC) .LE. TOLDC .AND. ABS(DYC) .LE. TOLDC) THEN
!
          FIND = .TRUE.
          XC = XC -1.
          YC = YC -1.
          RETURN
        ENDIF
!
!       *** the derivated terms of the eqs. are evaluated and  ***
!       *** the eqs. are solved                                ***
        DDEN = DXDXC*DYDYC - DYDXC*DXDYC                                  30.80
        DXP  = XP - XVC                                                   30.80
        DYP  = YP - YVC                                                   30.80
        DXC  =  (DYDYC*DXP - DXDYC*DYP) / DDEN                            30.80
        DYC  = (-DYDXC*DXP + DXDXC*DYP) / DDEN                            30.80
!
        XC = XC + DXC
        YC = YC + DYC
!
        IF (ITEST .GE. 120 .OR. INTES .GE. 50 .OR. IOUTES .GE. 50)
     &    WRITE(PRINTF,42) DXC, DYC, XC-1., YC-1.                         40.00
 42     FORMAT(' (DXC,DYC)=', 2(1X,F10.2), ' (XC,YC)=', 2(1X,F10.2))      40.00
!
 14   CONTINUE
      RETURN
!     *** end of subroutine NEWTON ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWT1D (XP, YP, XCGRID, YCGRID, KGRPNT,                  40.00
     &                   MXITNR, XC, YC, FIND)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm2.inc'                                               30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!  0. Authors
!
!     40.00, 40.13: Nico Booij
!
!  1. Updates
!
!     40.00, Feb. 99: New (adaptation from subr NEWTON for 1D case)
!     40.13, Feb. 01: DX and DY renamed to DELX and DELY (DX and DY are
!                     common var.); error in expression for RS corrected
!                     PRINTF replaced by PRTEST in test output
!
!  2. Purpose
!
!     Solve eqs. and find a point  (XC,YC) in a curvilinear 1D grid (compt.
!     grid) for a given point (XP ,YP) in a cartesian grid (problem coord).
!
!  3. Method
!
!     In this subroutine the step on the computational grid is selected
!     for which
!
!           (X-X1).(X2-X1)
!     0 <= --------------- <= 1
!          (X2-X1).(X2-X1)
!
!     where X, X1 and X2 are vectors; X corresponds to (Xp,Yp)
!     X1 and X2 are two neighbouring grid points
!
!  4. Argument variables
!
! i   KGRPNT: Grid adresses                                               40.00
! i   MXITNR: Maximum number of iterations                                30.82
!
      INTEGER KGRPNT(MXC,MYC), MXITNR                                     40.00
!
!   o XC    : X-coordinate in computational coordinates                   30.82
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   XP    : X-coordinate in problem coordinates                         30.82
!   o YC    : Y-coordinate in computational coordinates                   30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
! i   YP    : Y-coordinate in problem coordinates                         30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                     30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                     30.82
!
!   o FIND  : Whether XC and YC are found                                 30.82
!
      LOGICAL FIND                                                        30.82
!
!     Local variables:

      REAL :: DELX, DELY   ! grid line                                    40.13
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'NEWT1D')
!
      IF (ITEST .GE. 120) THEN                                            40.13
        WRITE(PRTEST,*) ' Coordinates in subroutine NEWT1D '              40.13
        DO I = 1, MXC
          WRITE(PRTEST,30) I, XCGRID(I,1)+XOFFS ,YCGRID(I,1)+YOFFS        40.13
        ENDDO
      ENDIF
 30   FORMAT(2X,I5,2(2X,E12.4))
!
      FIND = .FALSE.
      DO 40 IX = 2 ,MXC
        IF (KGRPNT(IX-1,1).GT.1) THEN
          X1 = XCGRID(IX-1,1)
          Y1 = YCGRID(IX-1,1)
        ELSE
          GOTO 40
        ENDIF
        IF (KGRPNT(IX,1).GT.1) THEN
          X2 = XCGRID(IX,1)
          Y2 = YCGRID(IX,1)
        ELSE
          GOTO 40
        ENDIF
!       both ends of the step are valid grid points
!       now verify whether projection of (Xp,Yp) is within the step
        DELX = X2 - X1                                                    40.13
        DELY = Y2 - Y1                                                    40.13
        RS = ((XP - X1) * DELX + (YP - Y1) * DELY) /                      40.13
     &              (DELX * DELX + DELY * DELY)                           40.13
        IF (RS.GE.0. .AND. RS.LE.1.) THEN
          FIND = .TRUE.
          XC = REAL(IX-2) + RS                                            40.00
          YC = 0.
          GOTO 50
        ENDIF
  40  CONTINUE
  50  RETURN
!     *** end of subroutine NEWT1D ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE EVALF (XC ,YC ,XVC ,YVC ,XCGRID ,YCGRID)                 30.72
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!
!  2. Purpose
!
!     Evaluate the coordinates (in problem coordinates) of point (XC,YC)
!     given in compuational coordinates.
!
!  3. Method
!
!     Bilinear interpolation
!
!  4. Argument variables
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!       XC, YC      real, outp    point in computational grid coordinates
!       XVC, YCV    real, OUTP    same point  but in problem coordinates
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'EVALF')
!
      I  = INT(XC)
      J  = INT(YC)
!
!     *** If the guess point (XC,YC) is in the boundary   ***
!     *** where I = MXC or/and J = MYC the interpolation  ***
!     *** is done in the mesh with pivoting point         ***
!     *** (MXC-1, J) or/and (I,MYC-1)                     ***
!
      IF (I .EQ. MXC) I = I - 1
      IF (J .EQ. MYC) J = J - 1
      T = XC - FLOAT(I)
      U = YC - FLOAT(J)
!     *** For x-coord. ***
      P1 = XCGRID(I,J)                                                    30.72
      P2 = XCGRID(I+1,J)                                                  30.72
      P3 = XCGRID(I+1,J+1)                                                30.72
      P4 = XCGRID(I,J+1)                                                  30.72
      XVC = (1.-T)*(1.-U)*P1+T*(1.-U)*P2+T*U*P3+(1.-T)*U*P4
!     *** For y-coord. ***
      P1 = YCGRID(I,J)                                                    30.72
      P2 = YCGRID(I+1,J)                                                  30.72
      P3 = YCGRID(I+1,J+1)                                                30.72
      P4 = YCGRID(I,J+1)                                                  30.72
      YVC = (1.-T)*(1.-U)*P1+T*(1.-U)*P2+T*U*P3+(1.-T)*U*P4
      RETURN
!     *** end of subroutine EVALF ***
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOBST (XCGRID, YCGRID, KGRPNT, CROSS)                   40.31 30.70
!                                                                      *
!***********************************************************************

      USE M_OBSTA                                                         40.31

      IMPLICIT NONE                                                       40.04
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!     30.70
!     30.72  IJsbrand Haagsma
!     30.74  IJsbrand Haagsma
!     40.04  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.31  Marcel Zijlema
!
!  1. Updates
!
!     30.70, Feb. 98: check if neighbouring point is a true grid point
!                     loop over grid points moved from calling routine into this
!                     argument list changed
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.04, Nov. 99: IMPLICIT NONE added, header updated
!                   : Removed include files that are not used
!     40.28, Feb. 02: Adjustments for extended REFLECTION option
!     40.31, Oct. 03: changes w.r.t. obstacles
!
!  2. Purpose
!
!     Obtains all the data required to find obstacles and
!     use subroutine TCROSS2 to find them
!
!  3. Method
!
!  4. Argument variables
!
!     CROSS   output Array which contains 0's if there is no
!                    obstacle crossing
!                    if an obstacle is crossing between the
!                    central point and its neighbour CROSS is equal
!                    to the number of the obstacle
!     KGRPNT  input  Indirect addressing for computational grid points
!     XCGRID  input  Coordinates of computational grid in x-direction     30.72
!     YCGRID  input  Coordinates of computational grid in y-direction     30.72
!
      INTEGER KGRPNT(MXC,MYC), CROSS(2,MCGRD)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)                            30.72
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ICC     index
!     ICGRD   index
!     IENT    number of entries of this subroutine
!     ILINK   indicates which link is analyzed: 1 -> neighbour in x
!                                               2 -> neighbour in y
!     IX      counter of gridpoints in x-direction
!     IY      counter of gridpoints in y-direction
!     JJ      counter for number of obstacles
!     JP      counter for number of corner points of obstacles
!     NUMCOR  number of corner points of obstacle
!     X1, Y1  user coordinates of one end of grid link
!     X2, Y2  user coordinates of other end of grid link
!     X3, Y3  user coordinates of one end of obstacle side
!     X4, Y4  user coordinates of other end of obstacle side
!
      INTEGER    ICC, ICGRD, IENT, ILINK, IX, IY, JJ, JP
      INTEGER    NUMCOR
      REAL       X1, X2, X3, X4, Y1, Y2, Y3, Y4
      LOGICAL    XONOBST                                                  40.04
      TYPE(OBSTDAT), POINTER :: COBST                                     40.31
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!     Function TCROSS2                                                    40.04
!     STRACE
!
      LOGICAL   TCROSS2                                                   40.04
!
!  9. Subroutines calling
!
!     SWPREP
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!       ----------------------------------------------------------------
!       Read number of obstacles from array OBSTA
!       For every obstacle do
!           Read number of corners of the obstacle
!           For every corner of the obstacle do
!               For every grid point do
!                   CALL FUNCTION TCROSS2 to search if there is crossing  40.04
!                   point                                                 40.04
!                   between the line of two points of the stencil and the
!                   line of the corners of the obstacle.
!                   If there is crossing point then
!                   then CROSS(link,kcgrd) = number of the crossing obstacle
!                   else CROSS(link,kcgrd) = 0
!       ----------------------------------------------------------------
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWOBST')
!
      IF (NUMOBS .GT. 0) THEN
!       NUMOBS is the number of obstacles ***
        COBST => FOBSTAC                                                  40.31
        DO 120 JJ = 1, NUMOBS
!         number of corner points of the obstacle
          NUMCOR = COBST%NCRPTS
          IF (ITEST.GE. 120) THEN
            WRITE(PRINTF,50) JJ, NUMCOR
 50         FORMAT( ' Obstacle number : ', I4,'  has ',I4,' corners')
          ENDIF
!         *** X1 X2 X3 ETC. are the coordinates of point according ***
!         *** with the scheme in the subroutine TCROSS2 header     ***   40.04
          X3 = COBST%XCRP(1)                                              40.31
          Y3 = COBST%YCRP(1)                                              40.31
          IF (ITEST.GE. 120)  WRITE(PRINTF,30) 1,X3,Y3
          DO 110 JP = 2, NUMCOR
             X4 = COBST%XCRP(JP)                                          40.31
             Y4 = COBST%YCRP(JP)                                          40.31
           IF (ITEST.GE. 120) WRITE(PRINTF,30) JP,X4,Y4
  30        FORMAT(' Corner number:', I4,'    XP: ',E10.4,' YP: ',E11.4)
            DO 100 IX = 1, MXC
              DO 90 IY = 1, MYC
                ICC = KGRPNT(IX,IY)
                IF (ICC .GT. 1) THEN
                  X1 = XCGRID(IX,IY)                                      30.72
                  Y1 = YCGRID(IX,IY)                                      30.72
!
!                 *** "ILINK" indicates which link is analyzed. Initial  ***
!                 *** neighbour in x , second link with neighbouring in y***
                  DO 80 ILINK = 1, 2
                    IF (ILINK.EQ.1 .AND. IX.GT.1) THEN
                      X2    = XCGRID(IX-1,IY)                             30.72
                      Y2    = YCGRID(IX-1,IY)                             30.72
                      ICGRD = KGRPNT(IX-1,IY)                             30.70
                    ELSE IF (ILINK.EQ.2 .AND. IY.GT.1) THEN
                      X2    = XCGRID(IX,IY-1)                             30.72
                      Y2    = YCGRID(IX,IY-1)                             30.72
                      ICGRD = KGRPNT(IX,IY-1)                             30.70
                    ELSE
                      ICGRD = 0
                    ENDIF
                    IF (ICGRD.GT.1) THEN                                  30.70
!
!                     *** All links are analyzed in each point otherwise the   ***
!                     *** boundaries can be excluded                           ***
!
                      IF (TCROSS2(X1, X2, X3, X4, Y1, Y2, Y3, Y4,         40.04
     &                           XONOBST)) THEN                           40.04
                        CROSS(ILINK,ICC) = JJ
                      ENDIF
                    ENDIF
  80              CONTINUE
                ENDIF
  90          CONTINUE
 100        CONTINUE
            X3 = X4
            Y3 = Y4
 110      CONTINUE
          IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT
          COBST => COBST%NEXTOBST
 120    CONTINUE
      ENDIF
!
      RETURN
! * end of subroutine SWOBST *
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION TCROSS2(X1, X2, X3, X4, Y1, Y2, Y3, Y4, X1ONOBST)
!                                                                      *
!***********************************************************************
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
      IMPLICIT NONE                                                       40.04
!
      INCLUDE 'ocpcomm4.inc'
!
!  0. Authors
!
!     40.00  Gerbrant van Vledder
!     40.04  Annette Kieftenburg
!
!  1. Updates
!
!       30.70, Feb 98: argument list simplified
!                      subroutine changed into logical function
!       40.00, Aug 98: division by zero prevented
!       40.04, Aug 99: method corrected, IMPLICIT NONE added, XCONOBST added,
!                      introduced TINY and EPSILON (instead of comparing to 0)
!                      replaced 0 < LMBD,MIU by  0 <= LMBD,MIU
!                      XCONOBST added to argument list
!
!  2. Purpose
!
!       Find if there is an obstacle crossing the stencil in used
!
!  3. Method
!
!     For the next situation (A, B and C are the points in the stencil,
!     D and E  are corners of the obstacle
!
!
!      obstacle --> D(X3,Y3)
!                    *
!                     *
!                      *
!        (X2,Y2)        * (XC,YC)
!            B-----------@--------------------------A (X1,Y1)
!                        ^*                         /
!                   _____| *                       /
!                  |        *                     /
!                  |         *                   /
!         crossing point      *                 /
!                              *               /
!                               E             /
!                              (X4,Y4)       /
!                                           C
!
!
!       The crossing point (@) should be found solving the next eqs.
!       for LMBD and MIU.
!
!       | XC |    | X1 |           | X2 - X1 |
!       |    | =  |    | +  LMBD * |         |
!       | YC |    | Y1 |           | Y2 - Y1 |
!
!
!       | XC |    | X3 |           | X4 - X3 |                            40.04
!       |    | =  |    | +  MIU  * |         |
!       | YC |    | Y3 |           | Y4 - Y3 |                            40.04
!
!
!     If solution exist and (0 <= LMBD <= 1 and 0 <= MIU <= 1)            40.04
!     there is an obstacle crossing the stencil
!
!  4. Argument variables
!
!     X1, Y1  inp    user coordinates of one end of grid link
!     X2, Y2  inp    user coordinates of other end of grid link
!     X3, Y3  inp    user coordinates of one end of obstacle side
!     X4, Y4  inp    user coordinates of other end of obstacle side
!     X1ONOBST outp   boolean which tells whether (X1,Y1) is on obstacle
!
      REAL       EPS, X1, X2, X3, X4, Y1, Y2, Y3, Y4
      LOGICAL    X1ONOBST
!
!  5. Parameter variables
!
!  6. Local variables
!
!     A,B,C,D    dummy variables
!     DIV1       denominator of value of LMBD (or MIU)
!     E,F        dummy variables
!     IENT       number of entries of function TCROSS2
!     LMBD       coefficient in vector equation for stencil points (or obstacle)
!     MIU        coefficient in vector equation for obstacle (or stencil points)
!
      INTEGER    IENT
      REAL       A, B, C, D, DIV1, E, F, LMBD, MIU
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWOBST
!     SWTRCF
!     OBSTMOVE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     Calculate MIU and LMBD
!     If 0 <= MIU, LMBD <= 1                                              40.04
!     Then TCROSS2 is .True.
!     Else TCROSS2 is .False.
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'TCROSS2')
!
      EPS = EPSILON(X1)*SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))             40.04
      IF (EPS ==0.) EPS = TINY(X1)                                        40.04
      A    = X2 - X1
!     A not equal to zero
      IF (ABS(A) .GT. TINY(X1)) THEN                                      40.04
        B    = X4 - X3
        C    = X3 - X1
        D    = Y2 - Y1
        E    = Y4 - Y3
        F    = Y3 - Y1
      ELSE
!       exchange MIU and LMBD                                             40.04
        A    = X4 - X3
        B    = X2 - X1
        C    = X1 - X3
        D    = Y4 - Y3
        E    = Y2 - Y1
        F    = Y1 - Y3
      ENDIF
      DIV1 = ((A*E) - (D*B))                                              40.00
!
!     DIV1 = 0 means that obstacle is parallel to line through            40.04
!     stencil points, or (X3,Y3) = (X4,Y4);                               40.04
!     A = 0 means trivial set of equations X4= X3 and X2 =X1              40.04
!
      IF ((ABS(DIV1).LE.TINY(X1)) .OR.                                    40.04
     &    (ABS(A).LE.TINY(X1))) THEN                                      40.04
        MIU = -1.                                                         40.00
        LMBD = -1.                                                        40.04
      ELSE                                                                40.00
        MIU  = ((D*C) - (A*F)) / DIV1
        LMBD = (C + (B*MIU)) / A
      END IF                                                              40.00
!
      IF (MIU  .GE. 0. .AND. MIU  .LE. 1. .AND.                           40.04
     &    LMBD .GE. 0. .AND. LMBD .LE. 1.) THEN                           40.04
!
!       Only (X1,Y1) is checked, because of otherwise possible double     40.04
!       counting                                                          40.04
        IF ((LMBD.LE.EPS .AND. ABS(X2-X1).GT.EPS).OR.                     40.04
     &      (MIU .LE.EPS .AND. ABS(X2-X1).LE.EPS))THEN                    40.04
          X1ONOBST = .TRUE.                                               40.04
        ELSE                                                              40.04
          X1ONOBST = .FALSE.                                              40.04
        ENDIF                                                             40.04
!
!       *** test output ***
        IF (ITEST .GE. 120) THEN
          WRITE(PRINTF,70)X1,Y1,X2,Y2,X3,Y3,X4,Y4
  70      FORMAT(' Obstacle crossing  :',/,
     &    ' Coordinates of comp grid points and corners of obstacle:',/,
     &    ' P1(',E10.4,',',E10.4,')',' P2(',E10.4,',',E10.4,')',/,
     &    ' P3(',E10.4,',',E10.4,')',' P4(',E10.4,',',E10.4,')')
        ENDIF
!
        TCROSS2 = .TRUE.
      ELSE
        TCROSS2 = .FALSE.
      ENDIF
!
!     End of subroutine TCROSS2
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      RECURSIVE SUBROUTINE OBSTMOVE (XCGRID, YCGRID, KGRPNT)              40.31
!                                                                      *
!***********************************************************************

      USE M_OBSTA                                                         40.31

      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm3.inc'
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
!     40.09  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.31  Marcel Zijlema
!
!  1. Updates
!
!     40.09, July 00: new subroutine
!     40.28, Feb. 02: Adjustments for extended REFLECTION option
!     40.31, Oct. 03: changes w.r.t. obstacles
!
!  2. Purpose
!
!     Move OBSTACLE points (X3,Y3) and (X4,Y4) a bit if computational gridcell
!     (X1,Y1) is on the OBSTACLE line piece.
!
!  3. Method
!
!     Add EPS*(dY,-dX) to OBSTACLE line piece coordinates so that movement of
!     these OBSTACLE points is perpendicular to the direction
!
!  4. Argument variables
!
!     KGRPNT  input  Indirect addressing for computational grid points
!     XCGRID  input  Coordinates of computational grid in x-direction
!     YCGRID  input  Coordinates of computational grid in y-direction
!
      INTEGER KGRPNT(MXC,MYC)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)
!
!  5. Parameter variables
!
!  6. Local variables
!
!     DISTA    distance between (X1,Y1) and (X2A,Y2A)
!     DISTB    distance between (X1,Y1) and (X2 ,Y2 )
!     DXA, DYA difference X1 - X2A respectively Y1 - Y2A
!     DXB, DYB difference X2 - X1 respectively Y2 - Y1
!     DXO, DYO difference X4 - X3 respectively Y4 - Y3
!     DXYO     distance between (X3,Y3) and (X4,Y4)
!     EPS      multiplication factor
!     ICC     index
!     ICGRD   index
!     IENT    number of entries of this subroutine
!     ILINK   indicates which link is analyzed: 1 -> neighbour in x
!                                               2 -> neighbour in y
!     IX      counter of gridpoints in x-direction
!     IY      counter of gridpoints in y-direction
!     JJ      counter for number of obstacles
!     JP      counter for number of corner points of obstacles
!     MOVED    boolean which tells whether OBSTACLE has been moved
!     NUMCOR  number of corner points of obstacle
!     X1, Y1   computational grid coordinates of one end of grid link
!     X2, Y2   computational grid coordinates of other end of grid link
!              i.e. neighbouring point of X1,Y1 associated with linknumber
!     X2A,Y2A  other neighbouring point of X1,Y1 associated with other
!              linknumber, or if invalid: = X2,Y2
!     X3, Y3   user coordinates of one end of obstacle side
!     X4, Y4   user coordinates of other end of obstacle side
!     XCONOBST boolean variable which tells whether XC in on OBSTACLE
!              line piece
!     XYEPS    displacement factor relative to local computational grid
!
      INTEGER    ICC, ICGRD, IENT, ILINK, IX, IY, JJ, JP
      INTEGER    NUMCOR
      REAL       DISTA, DISTB, DXA, DXB, DYA, DYB,
     &           DXO, DXYO, DYO, EPS, X1, X2, X2A, X3, X4,
     &           XYEPS, Y1, Y2, Y2A, Y3, Y4
      LOGICAL    MOVED, XCONOBST
      TYPE(OBSTDAT), POINTER :: COBST                                     40.31
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!     STRACE
!     Function TCROSS2 indicates whether there is an obstacle crossing the used
!                      stencil or not
!     MSGERR
!
      LOGICAL    TCROSS2
!
!  9. Subroutines calling
!
!     SWPREP
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!       ----------------------------------------------------------------
!       Read number of obstacles from array OBSTA
!       For every obstacle do
!           Read number of corners of the obstacle
!           For every corner of the obstacle do
!               For every grid point do
!                   CALL FUNCTION TCROSS2 to search if there is crossing point
!                   between the line of two points of the stencil and the
!                   line of the corners of the obstacle.
!                   If there is crossing point then
!                     If there computational gridpoint is on the obstacle
!                     then move corner points of obstacle perpendicular
!                     to obstacle with factor
!                     (XYEPS*DYO/DXYO,-XYEPS*DXO/DXYO)
!                     Moved = .True.
!                   If Moved then call OBSTMOVE again to check whether
!                   there are still computational grid points on OBSTACLE
!       ----------------------------------------------------------------
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'OBSTMOVE')
      MOVED = .FALSE.
      EPS =1.E-2
!
      IF (NUMOBS .GT. 0) THEN
        COBST => FOBSTAC                                                  40.31
        DO JJ = 1, NUMOBS
!         number of corner points of the obstacle
          NUMCOR = COBST%NCRPTS
          IF (ITEST.GE. 120) THEN
            WRITE(PRINTF,50) JJ, NUMCOR
 50         FORMAT( ' Obstacle number : ', I4,'  has ',I4,' corners')
          ENDIF
!         *** X1 X2 X3 ETC. are the coordinates of point according ***
!         *** with the scheme in the subroutine TCROSS2 header      ***
          X3 = COBST%XCRP(1)                                              40.31
          Y3 = COBST%YCRP(1)                                              40.31
          DO JP = 2, NUMCOR
            X4 = COBST%XCRP(JP)                                           40.31
            Y4 = COBST%YCRP(JP)                                           40.31
            IF (ITEST.GE. 120) WRITE(PRINTF,30) JP,X4,Y4
  30        FORMAT(' Corner number:', I4,'    XP: ',E10.4,' YP: ',E11.4)
            DO IX = 1, MXC
              DO IY = 1, MYC
                ICC = KGRPNT(IX,IY)
                IF (ICC .GT. 1) THEN
                  X1 = XCGRID(IX,IY)
                  Y1 = YCGRID(IX,IY)
!
!                 *** "ILINK" indicates which link is analyzed. Initial  ***
!                 *** neighbour in x , second link with neighbouring in y***
                  DO ILINK = 1, 2
                    IF (ILINK.EQ.1 .AND. IX.GT.1) THEN
                      X2    = XCGRID(IX-1,IY)
                      Y2    = YCGRID(IX-1,IY)
                      IF (IY.GT.1) THEN
                        X2A    = XCGRID(IX,IY-1)
                        Y2A    = YCGRID(IX,IY-1)
                      ELSE
                        X2A    = X2
                        Y2A    = Y2
                      ENDIF
                      ICGRD = KGRPNT(IX-1,IY)
                    ELSE IF (ILINK.EQ.2 .AND. IY.GT.1) THEN
                      X2    = XCGRID(IX,IY-1)
                      Y2    = YCGRID(IX,IY-1)
                      IF (IX.GT.1) THEN
                       X2A    = XCGRID(IX-1,IY)
                       Y2A    = YCGRID(IX-1,IY)
                      ELSE
                       X2A    = X2
                       Y2A    = Y2
                      ENDIF
                      ICGRD = KGRPNT(IX,IY-1)
                    ELSE
                     ICGRD = 0
                    ENDIF
                    IF (ICGRD.GT.1) THEN
!
!                     *** All links are analyzed in each point otherwise the   ***
!                     *** boundaries can be excluded                           ***
!
                      IF (TCROSS2(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XCONOBST))THEN
                        IF (XCONOBST) THEN
                        DXA=(X1-X2A)
                        DYA=(Y1-Y2A)
                        DXB=(X2-X1)
                        DYB=(Y2-Y1)
                        DISTA=SQRT(DXA*DXA+DYA*DYA)
                        DISTB=SQRT(DXB*DXB+DYB*DYB)
                        XYEPS = EPS*MIN(DISTA,DISTB)
                        DXO = X4-X3
                        DYO = Y4-Y3
                        DXYO = SQRT(DXO*DXO+DYO*DYO)
!                       -DXO/DXYO and DYO/DXYO are used (instead of -DXO
!                       and DYO) because otherwise displacement is dependent
!                       on length of OBSTACLE line piece
                        COBST%XCRP(JP-1) = X3 + XYEPS*DYO/DXYO            40.31
                        COBST%YCRP(JP-1) = Y3 - XYEPS*DXO/DXYO            40.31
                        COBST%XCRP(JP  ) = X4 + XYEPS*DYO/DXYO            40.31
                        COBST%YCRP(JP  ) = Y4 - XYEPS*DXO/DXYO            40.31
                        CALL MSGERR (1, 'Obstacle points moved')
                        WRITE(PRINTF, 17) X3, Y3, X4, Y4,
     &                        X3+XYEPS*DYO/DXYO, Y3-XYEPS*DXO/DXYO,
     &                        X4+XYEPS*DYO/DXYO, Y4-XYEPS*DXO/DXYO,X1,Y1
  17                    FORMAT ('OBSTACLE POINTS (', F11.2, ',',  F11.2,
     &                         '), and (', F11.2,',',  F11.2,'),',
     &                         'moved to: (',  F11.2,',',
     &                         F11.2,'), and (', F11.2,',', F11.2,
     &                         '), because OBSTACLE line piece ',
     &                         'was on computational grid point (',
     &                         F11.2,',', F11.2,').')
                        X3 = X3 + XYEPS * DYO/DXYO
                        Y3 = Y3 - XYEPS * DXO/DXYO
                        X4 = X4 + XYEPS * DYO/DXYO
                        Y4 = Y4 - XYEPS * DXO/DXYO
                        MOVED = .TRUE.
                        ENDIF
                      ENDIF
                    ENDIF
                  END DO
                ENDIF
              END DO
            END DO
            X3 = X4
            Y3 = Y4
          END DO
          IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT
          COBST => COBST%NEXTOBST
        END DO
      ENDIF
      IF (MOVED) CALL OBSTMOVE(XCGRID, YCGRID, KGRPNT)                    40.31
      RETURN
! * end of subroutine OBSTMOVE *
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWTRCF (CROSS, WLEV2, CHS,                               40.31 40.00
     &                   LINK, OBREDF                                     40.03
     &                  ,AC2,    IMATRA, KGRPNT, XCGRID,                  40.09
     &                   YCGRID, CAX,    CAY,    RDX,    RDY,    ANYBIN,  40.09
     &                   SPCSIG)                                          40.28
!                                                                      *
!***********************************************************************

      USE M_OBSTA                                                         40.31
      USE M_PARALL                                                        40.31

      IMPLICIT NONE                                                       40.09
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'                                               30.74
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
!     30.70
!     40.03  Nico Booij
!     40.08  Erick Rogers
!     40.09  Annette Kieftenburg
!     40.14  Annette Kieftenburg
!     40.18  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.30  Marcel Zijlema
!     40.31  Marcel Zijlema
!
!  1. Updates
!
!     30.70, Feb. 98: water level (WLEV2) replaced depth
!                     incident wave height introduced using argument
!                     CHS (sign. wave height in whole comput. grid)
!     40.03, Jul. 00: LINK1 and LINK2 in argumentlist replaced by LINK
!     40.09, Nov. 99: IMPLICIT NONE added, Method corrected
!                     Reflection option for obstacle added
!     40.14, Dec. 00: Reflection call corrected: reduced to neighbouring
!                     linepiece of obstacle (bug fix 40.11D)
!            Jan. 01: Constant waterlevel taken into account as well (bug fix 40.11E)
!     40.18, Apr. 01: Scattered reflection against obstacles added        40.18
!     40.28, Dec. 01: Frequency dependent reflection added
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.31, Oct. 03: changes w.r.t. obstacles
!
!  2. Purpose
!
!      take the value of transmission coefficient given
!      by the user in case  obstacle TRANSMISSION
!      or
!      compute the transmision coeficient in case obstacle DAM
!      based on Goda (1967) [from Seelig (1979)].
!      if reflections are turned on, calculate sourceterm in              40.09
!      subroutine REFLECT                                                 40.09
!
!  3. Method
!
!     Calculate transmission coefficient based on Goda (1967)             40.09
!     from Seelig (1979)                                                  40.09
!     Kt = 0.5*(1-sin {pi/(2*alpha)*(WATHIG/Hi +beta)})
!     where
!     Kt         transmission coefficient
!
!     alpha,beta coefficients dependent on structure of obstacle
!                and waves
!     WATHIG     = F = h-d is the freeboard of the dam, where h is the    40.09
!                crest level of the dam above the reference level and d   40.09
!                is the mean water level (relative to reference level)    40.09
!     Hi         incident (significant) wave height                       40.09
!                                                                         40.09
!     If reflection are switched on and obstacle is not exactly on line   40.09
!     of two neighbouring gridpoints, calculate reflections.              40.09
!
!  4. Argument variables
!
!     AC2      input     Action density array                             40.09
!     ANYBIN   input     Set a particular bin TRUE or FALSE depending on  40.09
!                        SECTOR                                           40.09
!     CAX      input     Propagation velocity                             40.09
!     CAY      input     Propagation velocity                             40.09
!     CHS      input     Hs in all computational grid points
!     CROSS    input..   Array which contains 0's if there is no
!                        obstacle crossing
!                        if an obstacle is crossing between the
!                        central point and its neighbour CROSS is equal
!                        to the number of the obstacle
!     IMATRA   inp/outp  Coefficients of right hand side of matrix        40.09
!                        equation                                         40.09
!     KCGRD    input     Grid address of points of computational stencil
!     LINK     input     indicates whether link in stencil                40.03
!                        crosses an obstacle                              40.03
!     OBREDF   output    Array of action density reduction coefficients
!                        (reduction at the obstacle)
!     RDX,RDY  input     Array containing spatial derivative coefficients 40.09
!     WLEV2    input     Water level in grid points
!
      INTEGER  CROSS(2,MCGRD)                                             40.00
      INTEGER  KGRPNT(MXC,MYC)
      INTEGER  LINK(2)                                                    40.03
      REAL     CHS(MCGRD), OBREDF(MDC,MSC,2), WLEV2(MCGRD)                30.70
      REAL     :: AC2(MDC,MSC,MCGRD)                                      40.09 40.22
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL     :: CAX(MDC,MSC,MICMAX), CAY(MDC,MSC,MICMAX)                40.09 40.22
      REAL     :: IMATRA(MDC,MSC), RDX(10), RDY(10)                       40.09 40.22 40.08
      REAL     :: SPCSIG(MSC)                                             40.28
      LOGICAL  :: ANYBIN(MDC,MSC)                                         40.09
!
!  5. Parameter variables
!
!  6. Local variables
!
!     AC2REF   reflected action density spectrum
!     ALOW     Lower limit for FVH
!     BUPL     Upper limit for FVH
!     FD1      Coeff. for freq. dep. reflection: vertical displacement    40.28
!     FD2      Coeff. for freq. dep. reflection: shape parameter          40.28
!     FD3      Coeff. for freq. dep. reflection: directional coefficient  40.28
!     FD4      Coeff. for freq. dep. reflection: bending point of freq.   40.28
!     FVH      WATHIG/Hsin (= F/Hi in formulation of Goda/Seelig
!              (1967/1979))
!     HGT      elevation of top of obstacle above reference level
!     HSIN     significant wave height in whole computational grid
!     ID       counter in directional space
!     IENT     number of entries of this subroutine
!     IS       counter in frequency space
!     ITRAS    indicates kind of obstacle: 0 -> transm
!                                          1 -> dam
!     JP       counter for number of corner points of obstacles
!     KDIF     fraction of energy reflected in a diffuse manner           40.28
!     LOOP     indicates which link is analyzed: 1 -> neighbour in x
!                                                2 -> neighbour in y
!     LREFDIFF inp  Indicates whether reflected energy should be          40.18
!                   reflected (1) or not (0)                              40.18
!     LRFRD     Indicates whether frequency dependent reflection is       40.28
!               active (#0.) or not (=0.)                                 40.28
!     NMPO     link number
!     NUMCOR   number of corner points of obstacle
!     OBET     user defined coefficient (beta) in formulation of
!              Goda/Seelig (1967/1979)
!     OBHKT    transmission coefficient according to Goda/Seelig (1967/1979)
!     OGAM     user defined coefficient (alpha) in formulation of
!              Goda/Seelig (1967/1979)
!     POWD     inp  User defined power of diffusive redist. function      40.28
!     POWS     inp  User defined power of scattered redistr.function      40.18
!     REFLCOEF reflection coefficient in terms of action density
!     SQRTREF  dummy variable
!     SQRTTRC  dummy variable
!     TRCF     transmission coefficient in terms of action density
!              (user defined or calculated (in terms of waveheigth))
!     X1, Y1   user coordinates of one end of grid link
!     X2, Y2   user coordinates of other end of grid link
!     X3, Y3   user coordinates of one end of obstacle side
!     X4, Y4   user coordinates of other end of obstacle side
!     XCGRID   Coordinates of computational grid in x-direction
!     XONOBST  Indicates whether computational point (X1,Y1) is on        40.14
!              obstacle                                                   40.14
!     YCGRID   Coordinates of computational grid in y-direction
!     WATHIG   freeboard of the dam (= HGT-waterlevel)
!
      INTEGER    ID, IENT, ITRAS, IS, JP, LOOP,
     &           NUMCOR, NMPO
      INTEGER    LREFDIFF, LRFRD                                          40.31
      REAL       ALOW, BUPL, FVH, HGT, HSIN, OBET, OBHKT,
     &           POWS, OGAM, REFLCOEF,                                    40.18 40.09
     &           TRCF, X1, X2, X3, X4, Y1, Y2, Y3, Y4, WATHIG
      REAL       KDIF, FD1, FD2, FD3, FD4, POWD                           40.31 40.28
      REAL       AC2REF(MDC,MSC), SQRTREF, SQRTTRC                        40.09
      LOGICAL    EXC, XGTL, XONOBST                                       40.14
      REAL       ICGRD, XCGRID(MXC,MYC), YCGRID(MXC,MYC)
      INTEGER    ICC, JJ
      TYPE(OBSTDAT), POINTER :: COBST                                     40.31
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!     OBSTLINE
!     REFLECT
!     Function TCROSS2                                                    40.14
!
      LOGICAL TCROSS2                                                     40.14
!
!  9. Subroutines calling
!
!     ACTION
!     SWOMPU                                                              30.70
!
! 10. Error messages
!
! 11. Remarks
!
!     Here the formulation of the transmission coefficients concerns the  40.09
!     ratio of action densities!                                          40.09
!
! 12. Structure
!
!     calculate transmission coefficients                                 40.09
!     if activated: calculate reflection source terms (if computational   40.09
!                   point not exactly on linepiece of obstacle)           40.09
!
! 13. Source text
! ======================================================================
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWTRCF')
!
      DO 30 LOOP = 1 ,2
        TRCF = 1.
        NMPO = LINK(LOOP)
        HSIN = CHS(KCGRD(LOOP+1))                                         40.03
        IF (NMPO .EQ. 0) GO TO 40
        COBST => FOBSTAC                                                  40.31
        DO JJ = 1, NMPO-1                                                 40.31
           IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT                      40.31
           COBST => COBST%NEXTOBST                                        40.31
        END DO                                                            40.31
        ITRAS  = COBST%TRTYPE                                             40.31
        IF (ITRAS .EQ. 0) THEN
!         User defined transmission coefficient concerns ratio of         40.09
!         waveheights, so:                                                40.09
          SQRTTRC = COBST%TRCOEF(1)                                       40.31
          TRCF = SQRTTRC * SQRTTRC                                        40.09
        ELSE IF (ITRAS .EQ. 1) THEN
          HGT    =  COBST%TRCOEF(1)                                       40.31
          OGAM   =  COBST%TRCOEF(2)                                       40.31
          OBET   =  COBST%TRCOEF(3)                                       40.31
!         level of dam above the water:
          WATHIG =  HGT - WLEV2(KCGRD(1)) - WLEV                          40.14 30.70
          IF (HSIN .LT. 0.1E-4) HSIN = 0.1E-4
!
!         *** Here the transmission coeff. is that of Goda and Seelig ***
          FVH  = WATHIG/HSIN
          ALOW = -OBET-OGAM                                               40.09
          BUPL = OGAM-OBET
!
          IF (FVH.LT.ALOW ) FVH = ALOW
          IF (FVH.GT.BUPL ) FVH = BUPL
          OBHKT = 0.5*(1.0-SIN(PI*(FVH+OBET)/(2.0*OGAM)))
          IF (TESTFL) WRITE (PRTEST, 20) IXCGRD(1)+MXF-2,                 40.30
     &           IYCGRD(1)+MYF-2,                                         40.30
     &           LOOP, HGT, WATHIG, HSIN, OBHKT                           40.01
  20      FORMAT (' test SWTRCF ', 2X, 3I5, ' dam level=', F6.2,
     &            ' depth=', F6.2, ' Hs=', F6.2, ' transm=', F6.3)        40.01
          IF (TESTFL .AND. ITEST.GE.140) WRITE (PRTEST, 22)
     &           OGAM, OBET, ALOW, BUPL, FVH                              40.01
  22      FORMAT (8X, 6E12.4)                                             40.01
!
!         Formulation of Goda/Seelig concerns ratio of waveheights.       40.09
!         Here we use action density so:                                  40.09
          TRCF = OBHKT*OBHKT
        ENDIF
!
!     *** REFLECTION ****
!     *** X1 X2 X3 ETC. are the coordinates of point according ***
!     *** with the scheme in the function TCROSS2 header       ***        40.04
        IF ( COBST%RFTYP1.GT.0. ) THEN                                    40.31
!         Reflections are activated                                       40.09
          SQRTREF  = COBST%RFCOEF(1)                                      40.31
          REFLCOEF = SQRTREF * SQRTREF                                    40.09
          LREFDIFF = COBST%RFTYP2                                         40.31
          POWS     = COBST%RFCOEF(2)                                      40.31
          POWD     = COBST%RFCOEF(3)                                      40.31
          KDIF     = COBST%RFCOEF(4)                                      40.31
          FD1      = COBST%RFCOEF(5)                                      40.31
          FD2      = COBST%RFCOEF(6)                                      40.31
          FD3      = COBST%RFCOEF(7)                                      40.31
          FD4      = COBST%RFCOEF(8)                                      40.31
          LRFRD    = COBST%RFTYP3                                         40.31
          X3 = COBST%XCRP(1)                                              40.31
          Y3 = COBST%YCRP(1)                                              40.31
          NUMCOR = COBST%NCRPTS                                           40.31
          DO JP = 2, NUMCOR                                               40.09
            X4 = COBST%XCRP(JP)                                           40.31
            Y4 = COBST%YCRP(JP)                                           40.31
            ICC = KCGRD(1)                                                40.09
            IF (ICC .GT. 1) THEN                                          40.09
              X1 = XCGRID(IXCGRD(1),IYCGRD(1))                            40.09
              Y1 = YCGRID(IXCGRD(1),IYCGRD(1))                            40.09
!                                                                         40.09
              ICGRD = 0                                                   40.09
              IF (KGRPNT(IXCGRD(LOOP+1),IYCGRD(LOOP+1)).GT.1) THEN        40.09
                X2    = XCGRID(IXCGRD(LOOP+1),IYCGRD(LOOP+1))             40.09
                Y2    = YCGRID(IXCGRD(LOOP+1),IYCGRD(LOOP+1))             40.09
                ICGRD = KCGRD(LOOP+1)                                     40.09
                IF (ICGRD.GT.1) THEN                                      40.09
                  IF (TCROSS2(X1, X2, X3, X4, Y1, Y2, Y3, Y4,             40.14
     &                        XONOBST)) THEN                              40.14
                    CALL OBSTLINE(X1,Y1,X2,Y2,X3,Y3,X4,Y4,XGTL,EXC)       40.09
                    CALL REFLECT(AC2, AC2REF, IMATRA, X1, Y1, X2, Y2,     40.09
     &                         X3, Y3, X4, Y4, XGTL, EXC, CAX,            40.09
     &                         CAY, RDX, RDY, LOOP,                       40.09
     &                         REFLCOEF, LREFDIFF, POWS, ANYBIN,          40.18 40.09
     &                         KDIF, POWD,                                40.28
     &                         LRFRD, SPCSIG, FD1, FD2, FD3, FD4)         40.28
                  ENDIF                                                   40.14
                ENDIF                                                     40.09
              ENDIF                                                       40.09
            ENDIF                                                         40.09
            X3 = X4                                                       40.09
            Y3 = Y4                                                       40.09
          END DO                                                          40.09
        END IF
!
        IF (ITEST .GE. 120)  WRITE (PRTEST,10)
     &  IXCGRD(1)-1, IYCGRD(1)-1, NMPO, TRCF
  10    FORMAT(' SWTRCF: Point=', 2I5, ' NMPO  = ', I5, ' transm ',       40.03
     &  F8.3)                                                             40.03
  40    DO IS = 1, MSC                                                   040697
          DO ID = 1, MDC
            OBREDF(ID,IS,LOOP) = TRCF
          ENDDO
        ENDDO
  30  CONTINUE
      RETURN
!     * end of SUBROUTINE SWTRCF
      END
!
!***********************************************************************  40.09
!                                                                      *  40.09
      SUBROUTINE OBSTLINE (X1, Y1, X2, Y2, X3, Y3, X4, Y4, XGTL, EXC)     40.09
!                                                                      *  40.09
!***********************************************************************  40.09
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
      IMPLICIT NONE                                                       40.09
!                                                                         40.09
      INCLUDE 'ocpcomm4.inc'                                              40.09
!                                                                         40.09
!  0. Authors                                                             40.09
!                                                                         40.09
!     40.09  Annette Kieftenburg                                          40.09
!                                                                         40.09
!  1. Updates                                                             40.09
!                                                                         40.09
!     40.09, Dec 99: Function created                                     40.09
!                                                                         40.09
!  2. Purpose                                                             40.09
!                                                                         40.09
!     Find out whether vector (X1,Y1) lies 'above' the line piece         40.09
!     through (X3,Y3) and (X4,Y4)                                         40.09
!                                                                         40.09
!  3. Method                                                              40.09
!                                                                         40.09
!     Calculate coefficients a and b in equation Y = A*X +B               40.09
!     Check whether Y1 is greater than A*X1 + B                           40.09
!     Then XGTL is true (else false)                                      40.09
!                                                                         40.09
!  4. Argument variables                                                  40.09
!                                                                         40.09
!     X1, Y1  inp    user coordinates of one end of grid link             40.09
!     X3, Y3  inp    user coordinates of one end of obstacle side         40.09
!     X4, Y4  inp    user coordinates of other end of obstacle side       40.09
!     EXC     outp   indicates whether X4 = X3, which results in 'excep-  40.09
!                    tional' situation (line parallel to y-axis)          40.09
!     XGTL    outp   indicates whether (X1,Y1) is situated 'above'        40.09
!                    linepiece (X3,Y3) (X4,Y4)                            40.09
!                                                                         40.09
      REAL       X1, X2, X3, X4, Y1, Y2, Y3, Y4                           40.09
      LOGICAL    XGTL, EXC                                                40.09
!                                                                         40.09
!  5. Parameter variables                                                 40.09
!                                                                         40.09
!  6. Local variables                                                     40.09
!                                                                         40.09
!     A,B        dummy variables                                          40.09
!     IENT       number of entries of subroutine OBSTLINE                 40.09
!     EPS        small real                                               40.09
!     RES        residual                                                 40.09
!                                                                         40.09
      INTEGER    IENT                                                     40.09
      REAL A, B, RES, EPS                                                 40.09
!                                                                         40.09
!  7. Common Blocks used                                                  40.09
!                                                                         40.09
!  8. Subroutines used                                                    40.09
!                                                                         40.09
!  9. Subroutines calling                                                 40.09
!                                                                         40.09
!     SWTRCF                                                              40.09
!                                                                         40.09
! 10. Error messages                                                      40.09
!                                                                         40.09
! 11. Remarks                                                             40.09
!                                                                         40.09
! 12. Structure                                                           40.09
!                                                                         40.09
!     If .NOT. (denominator of A = denominator of B = 0)                  40.09
!              (i.e. parallel to y-axis)                                  40.09
!     Then                                                                40.09
!       Calculate coefficients A and B in y = Ax + B                      40.09
!       EXC is .False.                                                    40.09
!       RESidual  = Y1 - (A*X1+B)                                         40.09
!       If RESidual <= Eps                                                40.09
!       If Residual > 0                                                   40.09
!       Then XGTL is .True.                                               40.09
!       Else XGTL is .False.                                              40.09
!     Else                                                                40.09
!       EXC is .True.                                                     40.09
!                                                                         40.09
! 13. Source text                                                         40.09
! ======================================================================  40.09
      SAVE IENT                                                           40.09
      DATA IENT/0/                                                        40.09
      IF (LTRACE) CALL STRACE (IENT,'OBSTLINE')                           40.09
!                                                                         40.09
      EPS = TINY(X1)                                                      40.18 40.09
      XGTL= .FALSE.                                                       40.09
!                                                                         40.09
      IF ( .NOT.( ABS(X4-X3).LE.EPS ) )THEN                               40.09
        A = (Y4-Y3)/(X4-X3)                                               40.09
        B = (X4*Y3-Y4*X3)/(X4-X3)                                         40.09
        EXC = .FALSE.                                                     40.09
!                                                                         40.09
        RES = Y1 - (A*X1+B)                                               40.09
        IF ( RES .GE. 0. ) THEN                                           40.09
          XGTL = .TRUE.                                                   40.09
        ELSE                                                              40.09
          XGTL = .FALSE.                                                  40.09
        ENDIF                                                             40.09
!                                                                         40.09
      ELSE                                                                40.09
        A = -9999.                                                        40.09
        B = -9999.                                                        40.09
        EXC = .TRUE.                                                      40.09
      END IF                                                              40.09
!                                                                         40.09
      RETURN                                                              40.09
!     End of subroutine OBSTLINE                                          40.09
      END                                                                 40.09
!
!************************************************************************ 40.09
!                                                                       * 40.09
      SUBROUTINE REFLECT (AC2, AC2REF, IMATRA, X1, Y1, X2, Y2, X3, Y3,    40.09
     &                    X4, Y4, XGTL, EXC, CAX, CAY, RDX, RDY,          40.09
     &                    LOOP, REF0, LREFDIFF, POWS, ANYBIN,             40.18 40.09
     &                    KDIF, POWD,                                     40.28
     &                    LRFRD,SPCSIG,FD1,FD2,FD3,FD4)                   40.28
!                                                                       * 40.09
!************************************************************************ 40.09
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
!                                                                         40.09
      IMPLICIT NONE                                                       40.09
!                                                                         40.09
      INCLUDE 'swcomm3.inc'                                               40.09
      INCLUDE 'ocpcomm4.inc'                                              40.18
!  0. Authors                                                             40.09
!                                                                         40.09
!     40.09  Annette Kieftenburg                                          40.09
!     40.18  Annette Kieftenburg                                          40.18
!     40.28  Annette Kieftenburg                                          40.28
!     40.38  Annette Kieftenburg                                          40.38
!                                                                         40.09
!  1. Updates                                                             40.09
!                                                                         40.09
!     40.09, Nov. 99: Subroutine created                                  40.09
!     40.18, Apr. 01: Scattered reflection against obstacles added        40.18
!     40.28, Dec. 01: Frequency dependent reflection added                40.28
!     40.38, Feb. 02: Diffuse reflection against obstacles added          40.38
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent   40.08
!                     with other subroutines                              40.08
!                                                                         40.09
!  2. Purpose                                                             40.09
!                                                                         40.09
!     Computation of REFLECTIONS near obstacles                           40.09
!                                                                         40.09
!  3. Method                                                              40.09
!                                                                         40.09
!     Determine the angle of the obstacle,                                40.09
!     Determine the angles between which reflections should be taken      40.09
!     into account                                                        40.09
!     If option is on: determine expression for frequency dependency      40.28
!     Determine redistribution function                                   40.18
!     Determine reflected action density (corrected for angle obstacle    40.09
!     and if option is on: redistribute energy)                           40.18
!     If coefficient is Kd not equal to zero: determine diffuse           40.38
!     reflection component                                                40.38
!     Add reflected spectrum to right hand side of matrix equation        40.09
!                                                                         40.09
!  4. Modules used                                                        40.18
!                                                                         40.18
!     --                                                                  40.18
!                                                                         40.18
!  5. Argument variables                                                  40.09
!                                                                         40.09
!     AC2      inp  Nonstationary case) action density as function        40.09
!                   of D,S,X,Y at time T+DT                               40.09
!     AC2REF   outp (Nonstationary case) reflected action density         40.09
!                   as function of D,S,X,Y at time T+DT                   40.09
!     CAX      inp  Propagation velocity                                  40.09
!     CAY      inp  Propagation velocity                                  40.09
!     EXC      inp  Indicates whether X4 = X3, which results in exception 40.09
!                   situation (line parallel to y-axis)                   40.09
!     FD1      inp  Coeff. freq. dep. reflection: vertical displacement   40.28
!     FD2      inp  Coeff. freq. dep. reflection: shape parameter         40.28
!     FD3      inp  Coeff. freq. dep. reflection: directional coefficient 40.28
!     FD4      inp  Coeff. freq. dep. reflection: bending point of freq.  40.28
!     IMATRA   i/o  Right hand side of matrix equation                    40.09
!     KDIF     inp  Fraction of energy reflected in a diffuse manner      40.38
!     LOOP     inp  Indicates which link is analyzed: 1 -> neighbour in x 40.09
!                                                     2 -> neighbour in y 40.09
!     LREFDIFF inp  Indicates whether reflected energy should be          40.18
!                   reflected (1) or not (0)                              40.18
!     LRFRD    inp  Indicates whether frequency dependent reflection is   40.28
!                   active (#0.) or not (=0.)                             40.28
!     POWD     inp  User defined power of diffuse redistr. function       40.38
!     POWS     inp  User defined power of scattered redistr. function     40.18
!     REF0     inp  User defined reflection coefficient                   40.18
!                   w.r.t. waveheight (0<=REF0<=1)                        40.18
!     RDX,RDY  inp  Array containing spatial derivative coefficients      40.09
!     SPCSIG   inp  Relative frequency (= 2*PI*Freq.)                     40.28
!     XGTL     inp  Indicates whether (X1,Y1) is situated 'above'         40.09
!                   linepiece (X3,Y3) (X4,Y4)                             40.09
!     X1, Y1   inp  Coordinates of computational grid point under         40.09
!                   consideration                                         40.09
!     X2, Y2   inp  Coordinates of computational grid point neighbour     40.09
!     X3, Y3   inp  User coordinates of one end of obstacle side          40.09
!     X4, Y4   inp  User coordinates of other end of obstacle side        40.09
!                                                                         40.09
      REAL       :: AC2(MDC,MSC,MCGRD), AC2REF(MDC,MSC)                   40.18
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL       :: CAX(MDC,MSC,MICMAX), CAY(MDC,MSC,MICMAX)              40.18 40.22
      REAL       :: IMATRA(MDC,MSC)                                       40.18
      REAL       :: RDX(10), RDY(10)                                      40.18 40.08
      REAL       :: FD1, FD2, FD3, FD4, SPCSIG(MSC)                       40.28
      REAL       :: REF0                                                  40.18
      REAL       :: X1, X2, X3, X4, Y1, Y2, Y3, Y4                        40.18
      LOGICAL    :: EXC, XGTL, ANYBIN(MDC,MSC)                            40.18
      INTEGER    :: LOOP                                                  40.18
      REAL       :: POWS                                                  40.18
      REAL       :: KDIF, POWD                                            40.38
      INTEGER    :: LREFDIFF, LRFRD
!                                                                         40.18
      INTENT (IN)     AC2, CAX, CAY,                                      40.18
     &                FD1, FD2, FD3, FD4, LRFRD,                          40.28
     &                RDX, RDY, REF0, SPCSIG, X1, X2,                     40.18
     &                X3, X4, Y1, Y2, Y3, Y4, EXC, XGTL, ANYBIN, LOOP,    40.18
     &                KDIF, POWD,                                         40.38
     &                POWS, LREFDIFF                                      40.18
      INTENT (IN OUT) IMATRA                                              40.18
      INTENT (OUT)    AC2REF                                              40.18
!                                                                         40.09
!  6. Parameter variables                                                 40.09
!                                                                         40.09
!     EPS2    constant used for redistribution function                   40.18
!                                                                         40.18
      REAL    :: EPS2                                                     40.18
!                                                                         40.18
!  7. Local variables                                                     40.09
!                                                                         40.09
!     AC2NEW  new action density spectrum (after optional redistribution) 40.18
!     ACSRED  redistributed action density spectrum                       40.18
!     ACDREF  diffusively reflected action density spectrum               40.28
!     ACTOTTO total action density directed towards obstacle              40.28
!     BETA    direction of obstacle                                       40.09
!     CORID   correction integer for angle of reflected wave              40.09
!     C1      counter for mirrored spectrum                               40.09
!     C2,C3   counters for mirrored spectrum (for obstacle with possitive 40.09
!             resp. negative angle)                                       40.09
!     dBETA   'residual' of (integer) correction CORID,                   40.09
!             used as weight factor                                       40.09
!     DELTA   direction of line piece (X1,Y1) (X4,Y4)                     40.09
!     DUMD    dummy variable for nearest integer for angle delta          40.09
!     DUMG    dummy variable for nearest integer for angle delta          40.09
!     DUMN1   dummy variable                                              40.18
!     DUMN2   dummy variable                                              40.18
!     EPS     gridsize dependent dummy (small) variable                   40.18
!     FD      frequency dependence function                               40.28
!     GAMMA_  direction of line piece (X1,Y1) (X3,Y3)                     40.09
!     ID,IS   counters for directions and frequencies                     40.09
!     IDD,IDG nearest integer for angle DELTA resp. GAMMA                 40.09
!     IDJ     counter for directions                                      40.18
!     IDMAX   maximum of IDG,IDD                                          40.09
!     IDMIN   minimum of IDG,IDD                                          40.09
!     IDR     counter for directions                                      40.18
!     IDRJ    counter for directions                                      40.18
!     IENT    number of entries of this subroutine                        40.09
!     KSCAT   fraction of energy reflected in a scattered manner          40.28
!     MAXID   absolute maximum for counter (ID) for width                 40.18
!     MMRW    maximum or minimum counter for P                            40.18
!                         (modulus of reflected width with minus sign)    40.18
!     NB      angle of normal of obstacle with angle BETA                 40.28
!     NORM    norm used to normalize redistribution function P            40.18
!     P       normalized redistribution function                          40.18
!     PMRW    maximum or minimum counter for P                            40.18
!                         (modulus of reflected width with minus sign)    40.18
!     R       reflection coefficient matrix                               40.09
!             distributed REFLECTED energy directed TOWARDS the obstacle  40.18
!     SUMP    total sum of redistribution function P                      40.18
!     WIDTH   counter used to determine width of redistribution function  40.18
!     WRES    resulting width of redistribution function P                40.18
!                                                                         40.09
      REAL    :: AC2NEW(MDC,MSC)                                          40.18
      REAL    :: ACSRED(MDC,MDC,MSC)                                      40.18
      REAL    :: ACDREF(MDC,MSC)                                          40.28
      REAL    :: ACTOTTO(MSC)                                             40.28
      REAL    :: BETA, dBETA, DELTA                                       40.09
      REAL    :: DUMN1, DUMN2                                             40.18
      REAL    :: DUMA, DUMB, FD(MSC)                                      40.28
      REAL    :: EPS, GAMMA_                                              40.09
      REAL    :: NORM, P(MDC,MDC)                                         40.18
      REAL    :: KSCAT, NB, PD(MDC), SUMPD                                40.28
      REAL    :: R(MDC,MSC)                                               40.09
      INTEGER :: C1, C2, C3, CORID, DUMD, DUMG, ID, IDD, IDG, IDMAX       40.09
      INTEGER :: IDMIN                                                    40.09
      INTEGER :: IDJ, IDR, IDRJ, MAXID,  WIDTH, WRES                      40.18
      INTEGER :: IENT, IS                                                 40.09
      LOGICAL :: STPNOW                                                   40.09
!                                                                         40.18
!  8. Subroutines used                                                    40.09
!                                                                         40.18
!     GAMMLN  natural logaritm of the (standard) GAMMA function           40.18
!                                                                         40.18
      REAL    :: GAMMLN                                                   40.18
!                                                                         40.09
!  9. Subroutines calling                                                 40.09
!                                                                         40.09
!     SWTRCF                                                              40.09
!                                                                         40.09
! 10. Error messages                                                      40.09
!                                                                         40.09
!     if obstacle linepiece is of length < EPS                            40.09
!                                                                         40.09
! 11. Remarks                                                             40.09
!                                                                         40.09
!    -In case the obstacle cuts exactly through computational grid point, 40.09
!     the obstacle should be moved a bit with subroutine OBSTMOVE.        40.09
!    -The order of magnitude of the obstacle linepiece is assumed to be   40.09
!     'long enough' compared to grid resolution (> 0.5*sqrt(dx^2+dy^2))   40.09
!     (if this restriction is violated, the reflections due to an obsta-  40.09
!     cle of one straight line can be very different from a similar line  40.09
!     consisting of several pieces (because only the directions of the    40.09
!     spectrum that are directed towards the obstacle linepiece are       40.09
!     reflected).                                                         40.09
!    -There should be only one intersection per computational gridcell.   40.09
!     Therefore it is better to avoid sharp edges in obstacles.           40.18
!                                                                         40.09
! 12. Structure                                                           40.09
!                                                                         40.09
!     Determine angle of obstacle, Beta                                   40.09
!       If Beta < - PI/2 then Beta + PI                                   40.09
!       If Beta >   PI/2 then Beta - PI                                   40.09
!     Calculate correction on counter for direction CORID and weight      40.09
!     factor dBeta                                                        40.09
!     Determine maximum and minimum angle for which reflection should     40.09
!     be taken into account                                               40.09
!                                                                         40.18
!     In option is on: determine redistribution function                  40.18
!                      determine redistibuted action density spectrum     40.18
!                      (first cut off directions, then redistribute)      40.18
!     Else cut off directions not pointing towards obstacle               40.18
!                                                                         40.18
!     Determine reflected action density                                  40.09
!     by reversing counter and correcting the angle with 2*Beta           40.09
!     i.e. correct counter by NINT(2*Beta) and take weighted average      40.09
!                                                                         40.09
!     Add reflected spectrum to right hand side of matrix equation        40.09
!                                                                         40.09
! 13. Source text                                                         40.09
!                                                                         40.18
      SAVE     IENT                                                       40.09
      DATA     IENT /0/                                                   40.09
      CALL STRACE (IENT, 'REFLECT')                                       40.09
!                                                                         40.09
!     initialization                                                      40.09
      DO IS = 1, MSC                                                      40.09
        FD(IS) = 0.                                                       40.28
        DO ID = 1, MDC                                                    40.09
          R(ID,IS) = 0.                                                   40.09
          AC2REF(ID,IS) = 0.                                              40.09
          AC2NEW(ID,IS) = 0.                                              40.18
          ACDREF(ID,IS) = 0.                                              40.28
          DO IDJ = 1, MDC                                                 40.18
            ACSRED(ID,IDJ,IS) = 0.                                        40.18
          END DO                                                          40.18
        END DO                                                            40.09
      END DO                                                              40.09
!                                                                         40.09
      EPS = EPSILON(X1)*SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))             40.09
      IF (EPS ==0.) EPS = TINY(X1)                                        40.09
      EPS2  = 0.01                                                        40.18
!                                                                         40.09
!     Determine angle of obstacle BETA, and related                       40.09
!                                                                         40.09
      IF (.NOT. ((ABS(Y4-Y3).LE.EPS) .AND. (ABS(X4-X3).LE.EPS)) ) THEN    40.09
        BETA = ATAN2((Y4-Y3),(X4-X3))                                     40.09
      ELSE                                                                40.09
        CALL MSGERR (3, 'Obstacle contains line piece of length 0!')      40.09
      END IF                                                              40.09
      IF (BETA.LE.-PI/2.) THEN                                            40.09
        BETA = BETA + PI                                                  40.09
      ELSE IF (BETA.GT.PI/2.) THEN                                        40.09
        BETA = BETA - PI                                                  40.09
      END IF                                                              40.09
!                                                                         40.09
      CORID = NINT(2.*BETA/DDIR)                                          40.09
      dBETA = (2.*BETA)/DDIR-CORID                                        40.09
!                                                                         40.28
      IF ((EXC .AND. (X1 .GT. X3)).OR.(.NOT.XGTL))  THEN                  40.28
        NB   = BETA + 3*PI/2                                              40.28
      ELSE                                                                40.28
        NB = BETA + PI/2                                                  40.28
      ENDIF                                                               40.28
!                                                                         40.09
!     Wave components that are directed to obstacle are lying between     40.09
!     GAMMA_ and DELTA.                                                   40.09
!     In case MDC is odd and the absolute angle (DELTA or GAMMA_) is      40.09
!     larger then PI/2 a trick is used to find the correct  nearest       40.09
!     integer.                                                            40.09
!                                                                         40.09
        GAMMA_ = ATAN2((Y3-Y1),(X3-X1))                                   40.09
!        IF ( (MOD(MDC,2).EQ.0) .OR.                                       40.09  40.18
!     &     ((MOD(MDC,2).EQ.1).AND.(ABS(GAMMA_).LT.PI/2.)) ) THEN          40.09  40.18
!          DUMG = NINT(GAMMA_/DDIR)                                        40.09  40.18
!        ELSE                                                              40.09  40.18
!          DUMG = SIGN(1.,GAMMA_)*(MDC+1)/2 -                              40.09  40.18
!     &           NINT((SIGN(1.,GAMMA_)*PI-GAMMA_)/DDIR)                   40.09  40.18
!        ENDIF                                                             40.09  40.18
        DUMG = NINT((GAMMA_+0.5*DDIR)/DDIR)                               40.18
        IF (DUMG.GT.0) THEN                                               40.09
          IDG = DUMG                                                      40.09
        ELSE                                                              40.09
          IF (DUMG.NE.0) THEN                                             40.09
            IDG = MDC+DUMG                                                40.18
          ELSE                                                            40.09
            IF (GAMMA_.LE.0.) THEN                                        40.09
              IDG = MDC                                                   40.09
            ELSE                                                          40.09
              IDG= 1                                                      40.09
            ENDIF                                                         40.09
          ENDIF                                                           40.09
        ENDIF                                                             40.09
!                                                                         40.09
        DELTA = ATAN2((Y4-Y1),(X4-X1))                                    40.09
!        IF ( (MOD(MDC,2).EQ.0) .OR.                                       40.09  40.18
!     &     ((MOD(MDC,2).EQ.1).AND.(ABS(DELTA).LT.PI/2.)) ) THEN           40.09  40.18
!          DUMD = NINT(DELTA/DDIR)                                         40.09  40.18
!        ELSE                                                              40.09  40.18
!          DUMD = SIGN(1.,DELTA)*(MDC+1)/2 -                               40.09  40.18
!     &           NINT((SIGN(1.,DELTA)*PI-DELTA)/DDIR)                     40.09  40.18
!        ENDIF                                                             40.09  40.18
        DUMD = NINT((DELTA+0.5*DDIR)/DDIR)                                40.18
        IF (DUMD.GT.0) THEN                                               40.09
          IDD = DUMD                                                      40.09
        ELSE                                                              40.09
          IF (DUMD.NE.0) THEN                                             40.09
            IDD = MDC+DUMD                                                40.18
          ELSE                                                            40.09
            IF (DELTA.LE.0.) THEN                                         40.09
              IDD = MDC                                                   40.09
            ELSE                                                          40.09
              IDD= 1                                                      40.09
            ENDIF                                                         40.09
          ENDIF                                                           40.09
        ENDIF                                                             40.09
!                                                                         40.09
      IF (IDG .GE. IDD) THEN                                              40.09
        IDMIN =IDD                                                        40.09
        IDMAX =IDG                                                        40.09
      ELSE                                                                40.09
        IDMIN =IDG                                                        40.09
        IDMAX =IDD                                                        40.09
      ENDIF                                                               40.09
!                                                                         40.09
!     For those directions that are pointed towards the obstacle the      40.09
!     reflection coefficients have a certain value. For the other         40.09
!     directions it is set to zero.                                       40.09
!                                                                         40.09
      IF (IDMAX .NE. IDMIN) THEN                                          40.09
        IF ( ( (EXC .AND. X3.GT.X1)  .OR.                                 40.18 40.09
     &         (.NOT.XGTL .AND. BETA.LE.0.).OR.                           40.09
     &         (XGTL .AND. BETA.GT.0.) )                                  40.09
     &       .AND.                                                        40.09
     &       ( (GAMMA_*DELTA.LT.0.) .OR.                                  40.09
     &       (GAMMA_*DELTA.EQ.0. .AND.(GAMMA_.GT.0. .OR. DELTA.GT.0.)) )  40.09
     &     ) THEN                                                         40.09
!                                                                         40.09
          DO IS = 1, MSC                                                  40.09
          DUMA = FD3*(SPCSIG(IS)/PI2-FD4)                                 40.28
          DUMB = REF0*FD2                                                 40.28
          IF ((LRFRD.EQ.1) .AND.(DUMA /= 0. .OR. DUMB/=0.)) THEN          40.28
            FD(IS)=FD1+FD2/PI*ATAN2((PI*DUMA),DUMB)                       40.28
              IF(FD(IS) > 1.) FD(IS) = 1.                                 40.28
              IF(FD(IS) < 0.) FD(IS) = 0.                                 40.28
            ELSE                                                          40.28
              FD(IS) = 1.                                                 40.28
            ENDIF                                                         40.28
            DO ID = 1, MDC                                                40.09
              IF ((ID.GE.IDMAX).OR.(ID.LE.IDMIN)) THEN                    40.09
                R(ID,IS) = REF0* FD(IS)                                   40.28 40.09
              ELSE                                                        40.09
                R(ID,IS) = 0.                                             40.09
              END IF                                                      40.09
            END DO                                                        40.09
          END DO                                                          40.09
        ELSE                                                              40.09
          DO IS = 1, MSC                                                  40.09
          DUMA = FD3 * (SPCSIG(IS)/PI2-FD4)                               40.28
          DUMB = REF0 * FD2                                               40.28
          IF ((LRFRD.EQ.1) .AND.(DUMA/= 0. .OR. DUMB/=0.)) THEN           40.28
            FD(IS)=FD1+FD2/PI*ATAN2((PI*DUMA),DUMB)                       40.28
              IF(FD(IS) > 1.) FD(IS) =1.                                  40.28
              IF(FD(IS) < 0.) FD(IS) =0.                                  40.28
            ELSE                                                          40.28
              FD(IS) = 1.                                                 40.28
            ENDIF                                                         40.28
            DO ID = 1, MDC                                                40.09
              IF ((ID.GE.IDMIN).AND.(ID.LE.IDMAX)) THEN                   40.09
                R(ID,IS) = REF0 * FD(IS)                                  40.28 40.09
              ELSE                                                        40.09
                R(ID,IS) = 0.                                             40.09
              END IF                                                      40.09
            END DO                                                        40.09
          END DO                                                          40.09
        END IF                                                            40.09
      ELSE                                                                40.09
        DO IS = 1, MSC                                                    40.28
          DUMA = FD3*(SPCSIG(IS)/PI2-FD4)                                 40.28
          DUMB = REF0*FD2                                                 40.28
          IF ((LRFRD.EQ.1) .AND.(DUMA/= 0. .OR. DUMB/=0.) )THEN           40.28
            FD(IS)=FD1+FD2/PI*ATAN2((PI*DUMA),DUMB)                       40.28
            IF(FD(IS) > 1.) FD(IS) = 1.                                   40.28
            IF(FD(IS) < 0.) FD(IS) = 0.                                   40.28
          ELSE                                                            40.28
            FD(IS) = 1.                                                   40.28
          ENDIF                                                           40.28
          R(IDMIN,IS) = REF0 * FD(IS)                                     40.28 40.09
        END DO                                                            40.09
      END IF                                                              40.09
!                                                                         40.18
!     Determine the width of the redistribution function                  40.18
      IF ((LREFDIFF.EQ.1) .AND. (KDIF < 1.)) THEN                         40.18 40.28
        MAXID  = INT(MDC/4)                                               40.18
        WIDTH = 0                                                         40.18
        DO ID = 1, MAXID                                                  40.18
          IF (((COS(ID*DDIR))**POWS).GE.EPS2 ) THEN                       40.18
            WIDTH = WIDTH +1                                              40.18
          END IF                                                          40.18
        END DO                                                            40.18
        WRES = WIDTH                                                      40.18
!                                                                         40.18
!       Initialize redistribution function                                40.18
!                                                                         40.18
        IF (WRES > 0) THEN                                                40.18
          DUMN1 = (1+POWS)/2.                                             40.18
          DUMN2 = 1+POWS/2.                                               40.18
          NORM = SQRT(PI)*EXP(GAMMLN(DUMN1)-GAMMLN(DUMN2))                40.18
          DO IDRJ = 1, MDC                                                40.18
            DO IDR = 1, MDC                                               40.18
              P(IDR,IDRJ) = 0.                                            40.18
            END DO                                                        40.18
          END DO                                                          40.18
!                                                                         40.18
!                                                                         40.18
!         Use property of cosine function to determine maximum and        40.18
!         minimum angle of redistribution                                 40.18
          DO IDR = 1, MDC                                                 40.18
            DO IDRJ = 1,MDC                                               40.18
              IF (COS(((IDR-IDRJ)*DDIR)) > 0.) THEN                       40.18
                P(IDR,IDRJ) = DDIR*(COS(((IDR-IDRJ)*DDIR))**POWS)/NORM    40.18
              ELSE                                                        40.18
                P(IDR,IDRJ) =0.                                           40.18
              ENDIF                                                       40.18
            END DO                                                        40.18
          END DO                                                          40.18
        ENDIF                                                             40.18
!                                                                         40.18
!       To avoid that the redistributed REFLECTed action density is       40.18
!       directed TOWARDS the obstacle, those directions are exluded.      40.18
!       Use the normal direction NB.                                      40.18
        DO IDR = 1, MDC                                                   40.18
          DO IDJ = 1, MDC                                                 40.18
            IF (COS(((IDR-IDJ)*DDIR)-NB) < 0.) THEN                       40.18
              P(IDR,IDJ) =0.                                              40.18
            ENDIF                                                         40.18
           END DO                                                         40.18
        END DO                                                            40.18
      ENDIF                                                               40.18
!                                                                         40.28
!     Determine the diffusive redistribution function                     40.38
      IF ((LREFDIFF.EQ.1) .AND. (KDIF > 0.)) THEN                         40.38
        SUMPD = 0.                                                        40.38
        DO IDRJ = 1,MDC                                                   40.38
          PD(IDRJ) = (COS((IDRJ-0.5)*DDIR-NB))**POWD                      40.38
          IF (PD(IDRJ) < EPS2)  THEN                                      40.38
            PD(IDRJ) = 0.                                                 40.38
          ELSE                                                            40.38
            SUMPD = SUMPD + PD(IDRJ)                                      40.38
          ENDIF                                                           40.38
        END DO                                                            40.38
        DO IDRJ = 1,MDC                                                   40.38
          PD(IDRJ) = PD(IDRJ)/SUMPD                                       40.38
        END DO                                                            40.38
      ENDIF                                                               40.38
!                                                                         40.18
      IF ((LREFDIFF.EQ.1) .AND. (WRES>0.) .AND. (KDIF<1.)) THEN           40.18 40.38
!     If reflectionis not specular                                        40.18
!     Determine redistributed action density spectrum per direction,      40.18
!     cutting of directions not directed towards the obstacle             40.18
          DO IS = 1, MSC                                                  40.18
            DO IDR = 1, MDC                                               40.18
              DO IDJ = 1, MDC                                             40.18
              ACSRED(IDR,IDJ,IS) = AC2(IDR,IS,KCGRD(1))*R(IDR,IS)         40.18
     &                                          *P(IDR,IDJ)               40.18
              END DO                                                      40.18
            END DO                                                        40.18
          END DO                                                          40.18
!         Sum for all components IDR to determine new action density      40.18
!         in component IDJ                                                40.18
          DO IS = 1, MSC                                                  40.18
            DO IDJ = 1, MDC                                               40.18
            AC2NEW(IDJ,IS) = 0.                                           40.18
              DO IDR = 1, MDC                                             40.18
              AC2NEW(IDJ,IS) = AC2NEW(IDJ,IS)+ACSRED(IDR,IDJ,IS)          40.18
              END DO                                                      40.18
            END DO                                                        40.18
          END DO                                                          40.18
      ELSE                                                                40.18
!       Else determine 'new action density' to be used for reflection by  40.18
!       only cutting of directions that are not directed towards the      40.18
!       obstacle                                                          40.18
        DO IS = 1, MSC                                                    40.18
          ACTOTTO(IS) = 0.                                                40.28
          DO ID = 1, MDC                                                  40.18
            AC2NEW(ID,IS) = AC2(ID,IS,KCGRD(1))*R(ID,IS)                  40.18
          END DO                                                          40.18
        END DO                                                            40.18
      ENDIF                                                               40.18
!                                                                         40.38
!     Determine total energy directed towards obstacle                    40.38
      IF (KDIF > 0.) THEN                                                 40.38
        DO IS = 1, MSC                                                    40.38
          ACTOTTO(IS) = 0.                                                40.38
          DO ID = 1, MDC                                                  40.38
            ACTOTTO(IS) = ACTOTTO(IS) + AC2(ID,IS,KCGRD(1)) * R(ID,IS)    40.38
          END DO                                                          40.38
        END DO                                                            40.38
      ENDIF                                                               40.38
!                                                                         40.09
      IF ((LREFDIFF.EQ.1) .AND. (KDIF > 0.)) THEN                         40.28 40.38
        DO IS = 1, MSC                                                    40.28
          DO ID = 1, MDC                                                  40.28
            ACDREF(ID,IS) = ACTOTTO(IS) * PD(ID)                          40.28
          END DO                                                          40.28
        END DO                                                            40.28
      ENDIF                                                               40.28
      KSCAT = 1. - KDIF                                                   40.28
!                                                                         40.09
!     Determine reflected action denstity spectrum by reversing and       40.09
!     correcting for obstacle angle in this counter (use C1), take        40.09
!     weighted average using dBETA and neighbour in theta direction (C2   40.09
!     for positive BETA or C3 for negative BETA)                          40.09
!     Add reflected spectrum to right hand side of matrix equation        40.09
!     (for every (ID,IS)).                                                40.18
!                                                                         40.09
      DO IS = 1, MSC                                                      40.09
        DO ID = 1, MDC                                                    40.09
          IF (ANYBIN(ID,IS)) THEN                                         40.09
            C1 = MOD((MDC-ID+1+CORID),MDC)                                40.09
            IF (C1.LE.0) C1 = C1+MDC                                      40.09
!                                                                         40.09
            IF (BETA .LT. 0.) THEN                                        40.09
              C2 = MOD((MDC-ID+1+CORID+1),MDC)                            40.09
              IF (C2.LE.0) C2 = C2+MDC                                    40.09
              IF (dBETA .GE. 0.) THEN                                     40.09
                AC2REF(ID,IS) =                                           40.09
     &            (RDX(LOOP)*CAX(ID,IS,1) + RDY(LOOP)*CAY(ID,IS,1))       40.09
     &            * (  ( (1.-dBETA)*AC2NEW(C1,IS)                         40.18
     &            + dBETA*AC2NEW(C2,IS) ) * KSCAT                         40.18  40.38
     &            + ACDREF(ID,IS) * KDIF  )                               40.38
              ELSE                                                        40.09
                AC2REF(ID,IS) =                                           40.09
     &            (RDX(LOOP)*CAX(ID,IS,1) + RDY(LOOP)*CAY(ID,IS,1))       40.09
     &            * (  ( (1.+dBETA)*AC2NEW(C1,IS)                         40.18
     &            - dBETA*AC2NEW(C2,IS) ) * KSCAT                         40.18  40.38
     &            + ACDREF(ID,IS) * KDIF  )                               40.38
              ENDIF                                                       40.09
              IMATRA(ID,IS) = IMATRA(ID,IS) + AC2REF(ID,IS)               40.09
            ELSE                                                          40.09
              C3 = MOD((MDC-ID+1+CORID-1),MDC)                            40.09
              IF (C3.LE.0) C3 = C3+MDC                                    40.09
              IF (dBETA .GE. 0.) THEN                                     40.09
                AC2REF(ID,IS) =                                           40.09
     &            (RDX(LOOP)*CAX(ID,IS,1) + RDY(LOOP)*CAY(ID,IS,1))       40.09
     &            * (  ( (1.-dBETA)*AC2NEW(C1,IS)                         40.18
     &            + dBETA*AC2NEW(C3,IS) ) * KSCAT                         40.18  40.38
     &            + ACDREF(ID,IS) * KDIF  )                               40.38
              ELSE                                                        40.09
                AC2REF(ID,IS) =                                           40.09
     &            (RDX(LOOP)*CAX(ID,IS,1) + RDY(LOOP)*CAY(ID,IS,1))       40.09
     &            * (  ( (1.+dBETA)*AC2NEW(C1,IS)                         40.18
     &            - dBETA*AC2NEW(C3,IS) ) * KSCAT                         40.18  40.38
     &            + ACDREF(ID,IS) * KDIF  )                               40.38
              ENDIF                                                       40.09
              IMATRA(ID,IS) = IMATRA(ID,IS) + AC2REF(ID,IS)               40.09
            END IF                                                        40.09
          END IF                                                          40.09
        END DO                                                            40.09
      END DO                                                              40.09
!                                                                         40.18
      IF (STPNOW()) RETURN                                                40.09
      RETURN                                                              40.09
!     End of subroutine REFLECT                                           40.09
      END                                                                 40.09
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SSHAPE (ACLOC, SPCSIG, SPCDIR, FSHAPL, DSHAPL)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!            Roeland Ris
!            Roberto Padilla
!     30.73: Nico Booij
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.02: IJsbrand Haagsma
!
!  1. Updates
!
!            Dec. 92: new for SWAN
!            Dec. 96: option MEAN freq. introduced see LOGPM
!     30.73, Nov. 97: revised in view of new boundary treatment
!     30.82, Sep. 98: Added error message in case of non-convergence
!     30.80, Oct. 98: correction suggested by Mauro Sclavo, and renames
!                     computation of tail added to improve accuracy
!     30.82, Oct. 98: Updated description of several variables
!     40.02, Oct. 00: Modified test write statement to avoid division by MS=0
!
!  2. Purpose
!
!     Calculating of energy density at boundary point (x,y,sigma,theta)
!
!  3. Method (updated...)
!
!     see: M. Yamaguchi: Approximate expressions for integral properties
!          of the JONSWAP spectrum; Proc. JSCE, No. 345/II-1, pp. 149-152,
!          1984.
!
!     computation of mean period: see Swan system documentation
!
!  4. Argument variables
!
!   o ACLOC : Energy density at a point in space
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL    ACLOC(MDC,MSC)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.82
!
! i   DSHAPL: Directional distribution
! i   FSHAPL: Shape of spectrum:
!             =1; Pierson-Moskowitz spectrum
!             =2; Jonswap spectrum
!             =3; bin
!             =4; Gauss curve
!             (if >0: period is interpreted as peak per.
!              if <0: period is interpreted as mean per.)
!
      INTEGER FSHAPL, DSHAPL                                              40.00
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!     LSHAPE   absolute value of FSHAPL
!
      INTEGER  ID, IS, LSHAPE
!
!     PKPER    peak period                                                30.80
!     APSHAP   aux. var. used in computation of spectrum
!     AUX1     aux. variable
!     AUX2          ,,
!     AUX3          ,,
!     COEFF    coefficient for behaviour around the peak (Jonswap)
!     CPSHAP   aux. var. used in computation of spectrum
!     CTOT     total energy
!     CTOTT    total energy (used for comparison)
!     DIFPER   aux. var. used to select bin closest to given frequency
!     MPER
!     MS       power in directional distribution
!     RA       action density
!     SALPHA
!     SF       frequency (Hz)
!     SF4      SF**4
!     SF5      SF**5
!     FPK      frequency corresponding to peak period (1/PKPER)           30.80
!     FPK4     FPK**4
!     SYF      peakedness parameter
!
      REAL     APSHAP, AUX1, AUX2, AUX3
      REAL     COEFF ,SYF   ,MPER  ,CTOT  ,CTOTT,PKPER  ,DIFPER
      REAL     MS
      REAL     RA    ,SALPHA,SF   ,SF4   ,SF5   ,FPK   ,FPK4
!
!     LOGPM    indicates whether peak or mean frequency is used
!     DVERIF   logical used in verification of incident direction
!
      LOGICAL  LOGPM, DVERIF                                              40.00
!
!  7. Common Blocks used
!
!     PSHAPE   coefficients of spectral distribution (see remarks)
!     SPPARM   array containing integral wave parameters (see remarks)
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
! 10. Error messages
!
! 11. Remarks
!
!     PSHAPE(1): SY0, peak enhancement factor (gamma) in Jonswap spectrum
!     PSHAPE(2): spectral width in case of Gauss spectrum in rad/s
!
!     SPPARM    real     input    incident wave parameters (Hs, Period,
!                                 direction, Ms (dir. spread))
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by the user (either peak or mean)      30.80
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!
!     ---------------------------------------------------------------------
!
!     In the case of a JONSWAP spectrum the initial conditions are given by
!                   _               _       _       _       _
!                  |       _   _ -4  |     |       | S - S   |
!             2    |      |  S  |    |     |       |      p  |
!          a g     |      |  _  |    |  exp|-1/2 * |________ |* 2/pi COS(T-T  )
! E(S,D )= ___  exp|-5/4 *|  S  |    | G   |       | e * S   |              wi
!      wa    5     |      |   p |    |     |_      |_     p _|
!           S      |      |_   _|    |
!                  |_               _|
!
!   where
!         S   : rel. frequency
!
!         D   : Dir. of wave component
!          wa
!
!         a   : equili. range const. (Phillips' constant)
!         g   : gravity aceleration
!
!         S   : Peak frequency
!          p
!
!         G   : Peak enhancement factor
!         e   : Peak width
!
!         T   : local wind direction
!          wi
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       case shape
!       =1:   calculate value of Pierson-Moskowitz spectrum
!       =2:   calculate value of Jonswap spectrum
!       else: Give error message because of wrong shape
!       ----------------------------------------------------------------
!       if LOGPM is True
!       then calculate average period
!            if it differs from given average period
!            then recalculate peak period
!                 restart procedure to compute spectral shape
!       ----------------------------------------------------------------
!       for all spectral bins do
!            multiply all action densities by directional distribution
!       ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE(IENT,'SSHAPE')
!
      IF (ITEST.GE.80) WRITE (PRTEST, 8) FSHAPL, DSHAPL,
     &      (SPPARM(JJ), JJ = 1,4)
   8  FORMAT (' entry SSHAPE ', 2I3, 4E12.4)
      IF (FSHAPL.LT.0) THEN
        LSHAPE = - FSHAPL
        LOGPM  = .FALSE.
      ELSE
        LSHAPE = FSHAPL
        LOGPM  = .TRUE.
      ENDIF
!
      IF (SPPARM(1).LE.0.)                                                40.31
     &   CALL MSGERR(1,' sign. wave height at boundary is not positive')  40.31
!
      PKPER = SPPARM(2)
      ITPER = 0
      IF (LSHAPE.EQ.3) THEN
!       select bin closest to given period
        DIFPER = 1.E10
        DO IS = 1, MSC
          IF (ABS(PKPER - PI2/SPCSIG(IS)) .LT. DIFPER) THEN
            ISP = IS
            DIFPER = ABS(PKPER - PI2/SPCSIG(IS))
          ENDIF
        ENDDO
      ENDIF
!
!     compute spectral shape using peak period PKPER                      30.80
!
 100  FPK  = (1./PKPER)                                                   30.80
      FPK4 = FPK**4
      IF (LSHAPE.EQ.1) THEN
        SALPHA = ((SPPARM(1) ** 2) * (FPK4)) * 5. / 16.
      ELSE IF (LSHAPE.EQ.2) THEN
!       *** SALPHA = alpha*(grav**2)/(2.*pi)**4)
        SALPHA = (SPPARM(1)**2 * FPK4) /
     &             ((0.06533*(PSHAPE(1)**0.8015)+0.13467)*16.)
      ENDIF
!
      CTOTT = 0.
      DO 300 IS = 1, MSC                                                  30.80
!
        IF (LSHAPE.EQ.1) THEN
!         *** LSHAPE = 1 : Pierson and Moskowitz ***
          SF = SPCSIG(IS) / PI2
          SF4 = SF**4
          SF5 = SF**5
          RA = (SALPHA/SF5)*EXP(-(5.*FPK4)/(4.*SF4))/(PI2*SPCSIG(IS))
          ACLOC(MDC,IS) = RA
        ELSE IF (LSHAPE.EQ.2) THEN
!         *** LSHAPE = 2 : JONSWAP ***
          SF = SPCSIG(IS)/(PI2)
          SF4 = SF**4
          SF5 = SF**5
          CPSHAP = 1.25 * FPK4 / SF4
          IF (CPSHAP.GT.10.) THEN                                         30.50
            RA = 0.
          ELSE
            RA = (SALPHA/SF5) * EXP(-CPSHAP)
          ENDIF
          IF (SF .LT. FPK) THEN
            COEFF = 0.07
          ELSE
            COEFF = 0.09
          ENDIF
          APSHAP =  0.5 * ((SF-FPK) / (COEFF*FPK)) **2
          IF (APSHAP.GT.10.) THEN                                         30.50
            SYF = 1.
          ELSE
            PPSHAP = EXP(-APSHAP)
            SYF = PSHAPE(1)**PPSHAP
          ENDIF
          RA = SYF*RA/(SPCSIG(IS)*PI2)
          ACLOC(MDC,IS) = RA
          IF (ITEST.GE.120) WRITE (PRTEST, 112)
     &                 SF, SALPHA, CPSHAP, APSHAP, SYF, RA
 112      FORMAT (' SSHAPE freq. ', 8E12.4)
        ELSE IF (LSHAPE.EQ.3) THEN
!
!         *** all energy concentrated in one BIN ***
!
          IF (IS.EQ.ISP) THEN
            ACLOC(MDC,IS) = ( SPPARM(1)**2 ) /
     &                     ( 16. * SPCSIG(IS)**2 * FRINTF )
          ELSE
            ACLOC(MDC,IS) = 0.
          ENDIF
        ELSE IF (LSHAPE.EQ.4) THEN
!
!         *** energy Gaussian distributed (wave-current tests) ***
!
          AUX1 = SPPARM(1)**2 / ( 16.* SQRT (PI2) * PSHAPE(2))
          AUX2 = ( SPCSIG(IS) - ( PI2 / PKPER ) )**2
          AUX3 = 2. * PSHAPE(2)**2
          RA = AUX1 * EXP ( -1. * AUX2 / AUX3 ) / SPCSIG(IS)
          ACLOC(MDC,IS) = RA
        ELSE
          IF (IS.EQ.1) THEN
            CALL MSGERR (2,'Wrong type for frequency shape')
            WRITE (PRINTF, *) ' -> ', FSHAPL, LSHAPE
          ENDIF
        ENDIF
        IF (ITEST.GE.10)
     &        CTOTT = CTOTT + FRINTF * ACLOC(MDC,IS) * SPCSIG(IS)**2
 300  CONTINUE
      IF (ITEST.GE.10) THEN
        IF (SPPARM(1).GT.0.01) THEN
          HSTMP = 4. * SQRT(CTOTT)
          IF (ABS(HSTMP-SPPARM(1)) .GT. 0.1*SPPARM(1))
     &    WRITE (PRINTF, 303) SPPARM(1), HSTMP
 303      FORMAT (' SSHAPE, deviation in Hs, should be ', F8.3,
     &            ', calculated ', F8.3)
        ENDIF
      ENDIF
!
!     if mean frequency was given recalculate PKPER and start anew
!
      IF (.NOT.LOGPM .AND. ITPER.LT.10) THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.
        AM1 = 0.
        DO IS = 1, MSC
          AS2 = ACLOC(MDC,IS) * (SPCSIG(IS))**2
          AS3 = AS2 * SPCSIG(IS)
          AM0 = AM0 + AS2
          AM1 = AM1 + AS3
        ENDDO
!       contribution of tail to total energy density
        PPTAIL = PWTAIL(1) - 1.                                           30.80
        APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))              30.80
        AM0 = AM0 * FRINTF + APTAIL * AS2                                 30.80
        PPTAIL = PWTAIL(1) - 2.                                           30.80
        EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))              30.80
        AM1 = AM1 * FRINTF + EPTAIL * AS3                                 30.80
!       Mean period:
        IF ( AM1.NE.0. ) THEN                                             40.31
           MPER = PI2 * AM0 / AM1
        ELSE                                                              40.31
           CALL MSGERR(3, ' first moment is zero in calculating the')     40.31
           CALL MSGERR(3, ' spectrum at boundary using param. bc.')       40.31
        END IF                                                            40.31
        IF (ITEST.GE.80) WRITE (PRTEST, 72) ITPER, SPPARM(2), MPER,
     &          PKPER
  72    FORMAT (' SSHAPE iter=', I2, '  period values:', 3F7.2)
        IF (ABS(MPER-SPPARM(2)) .GT. 0.01*SPPARM(2)) THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPARM(2) / MPER) * PKPER                              30.80
          GOTO 100
        ENDIF
      ENDIF
!
      IF (ITPER.GE.10) THEN
        CALL MSGERR(3, 'No convergence calculating the spectrum')         30.82
        CALL MSGERR(3, 'at the boundary using parametric bound. cond.')   30.82
      ENDIF
!
!       now introduce distribution over directions
!
      ADIR = PI * DEGCNV(SPPARM(3)) / 180.                                40.00
      IF (DSHAPL.EQ.1) THEN
        DSPR = PI * SPPARM(4) / 180.
        MS = MAX (DSPR**(-2) - 2., 1.)
      ELSE
        MS = SPPARM(4)
      ENDIF
      IF (MS.LT.12.) THEN
        CTOT = (2.**MS) * (GAMMA(0.5*MS+1.))**2 / (PI * GAMMA(MS+1.))
      ELSE
        CTOT =  SQRT (0.5*MS/PI) / (1. - 0.25/MS)
      ENDIF
      IF (ITEST.GE.100) THEN
        ESOM = 0.
        DO IS = 1, MSC
          ESOM = ESOM + FRINTF * SPCSIG(IS)**2 * ACLOC(MDC,IS)
        ENDDO
        WRITE (PRTEST, *) ' SSHAPE dir ', 4.*SQRT(ABS(ESOM)),
     &        SPPARM(1), CTOT, MS, GAMMA(0.5*MS+1.), GAMMA(MS+1.),
     &        CTOT                                                        40.02
      ENDIF
      DVERIF = .FALSE.
      CTOTT = 0.
      DO ID = 1, MDC
        ACOS = COS(SPCDIR(ID,1) - ADIR)
        IF (ACOS .GT. 0.) THEN
          CDIR = CTOT * MAX (ACOS**MS, 1.E-10)
          IF (.NOT.FULCIR) THEN
            IF (ACOS .GE. COS(DDIR)) DVERIF = .TRUE.
          ENDIF
        ELSE
          CDIR = 0.
        ENDIF
        IF (ITEST.GE.10) CTOTT = CTOTT + CDIR * DDIR
        IF (ITEST.GE.100) WRITE (PRTEST, 360) ID,SPCDIR(ID,1),CDIR
 360    FORMAT (' ID Spcdir Cdir: ',I3,3(1X,E10.4))
        DO IS = 1, MSC
          ACLOC(ID,IS) = CDIR * ACLOC(MDC,IS)
        ENDDO
      ENDDO
      IF (ITEST.GE.10) THEN
        IF (ABS(CTOTT-1.) .GT. 0.1) WRITE (PRINTF, 363) CTOTT
 363    FORMAT (' SSHAPE, integral of Cdir is not 1, but:', F6.3)
      ENDIF
      IF (.NOT.FULCIR .AND. .NOT.DVERIF)
     &   CALL MSGERR (1, 'incident direction is outside sector')
!
      RETURN
!
! End of subroutine SSHAPE
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE SINTRP (W1, W2, FL1, FL2, FL, SPCDIR, SPCSIG)
!                                                                  *
!*******************************************************************
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
!
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  Weimin Luo, Roeland Ris, Nico Booij            |
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
!     30.73: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.00: Nico Booij
!
!  1. Updates
!
!     30.01, Jan. 96: New subroutine for SWAN Ver. 30.01
!     30.73, Nov. 97: revised
!     40.00, Apr. 98: procedure to maintain peakedness introduced
!     30.82, Oct. 98: Update description of several variables
!     30.82, Oct. 98: Made arguments in ATAN2 double precision to prevent
!                     underflows
!
!  2. Purpose
!
!     interpolation of spectra
!
!  3. Method (updated...)
!
!     linear interpolation with peakedness maintained
!     interpolated average direction and frequency are determined
!     average direction and frequency of interpolated spectrum are determ.
!     shifts in frequency and direction are determined from spectrum 1 and
!     2 to the interpolated spectrum
!     bilinear interpolation in spectral space is used to calculate
!     contributions from spectrum 1 and 2.
!     in full circle cases interpolation crosses the boundary 0-360 degr.
!
!  4. Argument variables
!
!   o FL    : Interpolated spectrum.
! i   FL1   : Input spectrum 1.
! i   FL2   : Input spectrum 2.
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
! i   W1    : Weighting coefficient for spectrum 1.
! i   W2    : Weighting coefficient for spectrum 2.
!
      REAL    FL1(MDC,MSC), FL2(MDC,MSC), FL(MDC, MSC)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.82
      REAL    W1, W2
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!
      INTEGER  ID, IS
!
!     DOADD    indicates whether or not values have to be added
!
      LOGICAL  DOADD
!
!     ATOT1    integral over spectrum 1
!     ATOT2    integral over spectrum 2
!     AXTOT1   integral over x-component of spectrum 1
!     AXTOT2   integral over x-component of spectrum 2
!     AYTOT1   integral over y-component of spectrum 1
!     AYTOT2   integral over y-component of spectrum 2
!     ASTOT1   integral over Sigma * spectrum 1
!     ASTOT2   integral over Sigma * spectrum 2
!     ASIG1    average Sigma of spectrum 1
!     ASIG2    average Sigma of spectrum 2
!     DELD1    difference in direction between spectrum 1 and
!              the interpolated spectrum in number of directional steps
!     DELD2    same for spectrum 2
!     DELSG1   shift in frequency between spectrum 1 and interpolated
!              spectrum in number of frequency steps
!     DELSG2   same for spectrum 2
!
      REAL     ATOT1,  ATOT2,  AXTOT1, AXTOT2, AYTOT1, AYTOT2,
     &         ASTOT1, ASTOT2
      REAL     ASIG1,  ASIG2
      REAL     DELD1,  DELD2,  DELSG1, DELSG2
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!      SNEXTI, RBFILE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!      -----------------------------------------------------------------
!      If W1 close to 1
!      Then copy FL from FL1
!      Else If W2 close to 1
!           Then copy FL from FL2
!           Else determine total energy in FL1 and FL2
!                If energy of FL1 = 0
!                Then make FL = W2 * FL2
!                Else If energy of FL2 = 0
!                     Then make FL = W1 * FL1
!                     Else determine average direction of FL1 and FL2
!                          make ADIR = W1 * ADIR1 + W2 * ADIR2
!                          determine average frequency of FL1 and FL2
!                          make ASIG = W1 * ASIG1 + W2 * ASIG2
!                          determine directional shift from FL1
!                          determine directional shift from FL2
!                          determine frequency shift from FL1
!                          determine frequency shift from FL2
!                          For all spectral components do
!                              compose FL from components of FL1 and FL2
!      -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINTRP')
!
!     interpolation of spectra
!     ------------------------
!
      IF (W1.GT.0.99) THEN
        DO 101 ID=1,MDC
          DO 102 IS=1,MSC
            FL(ID,IS) = FL1(ID,IS)
 102      CONTINUE
 101    CONTINUE
      ELSE IF (W1.LT.0.01) THEN
        DO 201 ID=1,MDC
          DO 202 IS=1,MSC
            FL(ID,IS) = FL2(ID,IS)
 202      CONTINUE
 201    CONTINUE
      ELSE
        ATOT1  = 0.
        ATOT2  = 0.
        AXTOT1 = 0.
        AXTOT2 = 0.
        AYTOT1 = 0.
        AYTOT2 = 0.
        ASTOT1 = 0.
        ASTOT2 = 0.
        DO 301 ID=1,MDC
          DO 302 IS=1,MSC
            ATOT1  = ATOT1  + FL1(ID,IS)
            AXTOT1 = AXTOT1 + FL1(ID,IS) * SPCDIR(ID,2)
            AYTOT1 = AYTOT1 + FL1(ID,IS) * SPCDIR(ID,3)
            ASTOT1 = ASTOT1 + FL1(ID,IS) * SPCSIG(IS)
            ATOT2  = ATOT2  + FL2(ID,IS)
            AXTOT2 = AXTOT2 + FL2(ID,IS) * SPCDIR(ID,2)
            AYTOT2 = AYTOT2 + FL2(ID,IS) * SPCDIR(ID,3)
            ASTOT2 = ASTOT2 + FL2(ID,IS) * SPCSIG(IS)
 302      CONTINUE
 301    CONTINUE
        IF (ATOT1.LT.1.E-9) THEN
          DO 401 ID=1,MDC
            DO 402 IS=1,MSC
              FL(ID,IS) = W2*FL2(ID,IS)
 402        CONTINUE
 401      CONTINUE
        ELSE IF (ATOT2.LT.1.E-9) THEN
          DO 501 ID=1,MDC
            DO 502 IS=1,MSC
              FL(ID,IS) = W1*FL1(ID,IS)
 502        CONTINUE
 501      CONTINUE
        ELSE
!         determine interpolation factors in Theta space
          AXTOT  = W1 * AXTOT1 + W2 * AXTOT2
          AYTOT  = W1 * AYTOT1 + W2 * AYTOT2
          IF (ITEST.GE.80) THEN
            WRITE (PRTEST, 509)  ATOT1, ATOT2,                            40.02
     &            AXTOT, AXTOT1, AXTOT2, AYTOT, AYTOT1, AYTOT2
 509        FORMAT (' SINTRP factors ', 8E11.4, /, 15X, 4F7.3)
          ENDIF
!         DELD1 is the difference in direction between spectrum 1 and
!         the interpolated spectrum in number of directional steps
          DELD1  = REAL(ATAN2(DBLE(AXTOT*AYTOT1 - AYTOT*AXTOT1),          30.82
     &                        DBLE(AXTOT*AXTOT1 + AYTOT*AYTOT1))) / DDIR  30.82
!         DELD2 is the difference between spectrum 2 and
!         the interpolated spectrum
          DELD2  = REAL(ATAN2(DBLE(AXTOT*AYTOT2 - AYTOT*AXTOT2),          30.82
     &                        DBLE(AXTOT*AXTOT2 + AYTOT*AYTOT2))) / DDIR  30.82
          IDD1A  = NINT(DELD1)
          RDD1B  = DELD1 - REAL(IDD1A)
          IF (RDD1B .LT. 0.) THEN
            IDD1A = IDD1A - 1
            RDD1B = RDD1B + 1.
          ENDIF
          IDD1B  = IDD1A + 1
          RDD1B  = W1 * RDD1B
          RDD1A  = W1 - RDD1B
          IDD2A  = NINT(DELD2)
          RDD2B  = DELD2 - REAL(IDD2A)
          IF (RDD2B .LT. 0.) THEN
            IDD2A = IDD2A - 1
            RDD2B = RDD2B + 1.
          ENDIF
          IDD2B  = IDD2A + 1
          RDD2B  = W2 * RDD2B
          RDD2A  = W2 - RDD2B
!
!         determine interpolation factors in Sigma space
          ASIG1  = ASTOT1 / ATOT1
          ASIG2  = ASTOT2 / ATOT2
          ATOT   = W1 * ATOT1  + W2 * ATOT2
          ASTOT  = W1 * ASTOT1 + W2 * ASTOT2
          ASIG   = ASTOT / ATOT
!
!         DELSG1 is shift in frequency between spectrum 1 and interpolated
!         spectrum in number of frequency steps
          DELSG1 = ALOG (ASIG1 / ASIG) / FRINTF
          IDS1A  = NINT(DELSG1)
          RDS1B  = DELSG1 - REAL(IDS1A)
          IF (RDS1B .LT. 0.) THEN
            IDS1A = IDS1A - 1
            RDS1B = RDS1B + 1.
          ENDIF
          IDS1B  = IDS1A + 1
          RDS1A  = 1. - RDS1B
!
!         DELSG2 is shift in frequency between spectrum 2 and interpolated
!         spectrum in number of frequency steps
          DELSG2 = ALOG (ASIG2 / ASIG) / FRINTF
          IDS2A  = NINT(DELSG2)
          RDS2B  = DELSG2 - REAL(IDS2A)
          IF (RDS2B .LT. 0.) THEN
            IDS2A = IDS2A - 1
            RDS2B = RDS2B + 1.
          ENDIF
          IDS2B  = IDS2A + 1
          RDS2A  = 1. - RDS2B
!         test output
          IF (ITEST.GE.80) THEN
            WRITE (PRTEST, 510) ATOT, ATOT1, ATOT2,
     &            AXTOT, AXTOT1, AXTOT2, AYTOT, AYTOT1, AYTOT2,
     &            DELD1, DELD2, DELSG1, DELSG2
 510        FORMAT (' SINTRP factors ', 9E11.4, /, 15X, 4F7.3)
            WRITE (PRTEST, 512) IDS1A, RDS1A, IDS1B, RDS1B,
     &            IDS2A, RDS2A, IDS2B, RDS2B,
     &            IDD1A, RDD1A, IDD1B, RDD1B,
     &            IDD2A, RDD2A, IDD2B, RDD2B
 512        FORMAT (' SINTRP ', 8(I2, F7.3))
          ENDIF
!
          DO 601 ID=1,MDC
            DO 602 IS=1,MSC
              FL(ID,IS) = 0.
 602        CONTINUE
 601      CONTINUE
          DO 611 ID=1,MDC
            DOADD = .TRUE.
            ID1A = ID + IDD1A
            IF (FULCIR) THEN
              IF (ID1A.LT.1)   ID1A = ID1A + MDC
              IF (ID1A.GT.MDC) ID1A = ID1A - MDC
            ELSE
              IF (ID1A.LT.1)   DOADD = .FALSE.
              IF (ID1A.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 612 IS = MAX(1,1-IDS1A), MIN(MSC,MSC-IDS1A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1A * RDS1A * FL1(ID1A,IS+IDS1A)
 612          CONTINUE
              DO 613 IS = MAX(1,1-IDS1B), MIN(MSC,MSC-IDS1B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1A * RDS1B * FL1(ID1A,IS+IDS1B)
 613          CONTINUE
            ENDIF
 611      CONTINUE
          DO 621 ID=1,MDC
            DOADD = .TRUE.
            ID1B = ID + IDD1B
            IF (FULCIR) THEN
              IF (ID1B.LT.1)   ID1B = ID1B + MDC
              IF (ID1B.GT.MDC) ID1B = ID1B - MDC
            ELSE
              IF (ID1B.LT.1)   DOADD = .FALSE.
              IF (ID1B.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 622 IS = MAX(1,1-IDS1A), MIN(MSC,MSC-IDS1A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1B * RDS1A * FL1(ID1B,IS+IDS1A)
 622          CONTINUE
              DO 623 IS = MAX(1,1-IDS1B), MIN(MSC,MSC-IDS1B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1B * RDS1B * FL1(ID1B,IS+IDS1B)
 623          CONTINUE
            ENDIF
 621      CONTINUE
          DO 631 ID=1,MDC
            DOADD = .TRUE.
            ID2A = ID + IDD2A
            IF (FULCIR) THEN
              IF (ID2A.LT.1)   ID2A = ID2A + MDC
              IF (ID2A.GT.MDC) ID2A = ID2A - MDC
            ELSE
              IF (ID2A.LT.1)   DOADD = .FALSE.
              IF (ID2A.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 632 IS = MAX(1,1-IDS2A), MIN(MSC,MSC-IDS2A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2A * RDS2A * FL2(ID2A,IS+IDS2A)
 632          CONTINUE
              DO 633 IS = MAX(1,1-IDS2B), MIN(MSC,MSC-IDS2B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2A * RDS2B * FL2(ID2A,IS+IDS2B)
 633          CONTINUE
            ENDIF
 631      CONTINUE
          DO 641 ID=1,MDC
            DOADD = .TRUE.
            ID2B = ID + IDD2B
            IF (FULCIR) THEN
              IF (ID2B.LT.1)   ID2B = ID2B + MDC
              IF (ID2B.GT.MDC) ID2B = ID2B - MDC
            ELSE
              IF (ID2B.LT.1)   DOADD = .FALSE.
              IF (ID2B.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 642 IS = MAX(1,1-IDS2A), MIN(MSC,MSC-IDS2A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2B * RDS2A * FL2(ID2B,IS+IDS2A)
 642          CONTINUE
              DO 643 IS = MAX(1,1-IDS2B), MIN(MSC,MSC-IDS2B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2B * RDS2B * FL2(ID2B,IS+IDS2B)
 643          CONTINUE
            ENDIF
 641      CONTINUE
        ENDIF
      ENDIF
!
!     Test output
      IF (ITEST.GE.80) THEN
        A1 = 0.
        A2 = 0.
        AA = 0.
        DO 801 ID=1,MDC
          DO 802 IS=1,MSC
            A1 = MAX(A1,FL1(ID,IS))
            A2 = MAX(A2,FL2(ID,IS))
            AA = MAX(AA,FL(ID,IS))
 802      CONTINUE
 801    CONTINUE
        WRITE (PRTEST, *) ' SINTRP, maxima ', A1, A2, AA
      ENDIF
!
      RETURN
!  end of subroutine of SINTRP
      END
!***********************************************************************
!                                                                      *
      REAL FUNCTION DEGCNV (DEGREE)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'swcomm3.inc'                                               30.74
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees from nautical to carthesian or vice versa.
!
!  3. METHOD
!
!       DEGCNV = 180 + dnorth - degree
!
!  4. PARAMETERLIST
!
!       DEGCNV      direction in carthesian or nautical degrees.
!       DEGREE      direction in nautical or carthesian degrees.
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!           Nautical convention           Cartesian convention
!
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
!  9. STRUCTURE
!
!     ---------------------------------
!     IF (NAUTICAL DEGREES) THEN
!       CONVERT DEGREES
!     IF (DEGREES > 360 OR < 0) THEN
!       CORRECT DEGREES WITHIN 0 - 360
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'DEGCNV')
!
      IF ( BNAUT ) THEN
        DEGCNV = 180. + DNORTH - DEGREE
      ELSE
        DEGCNV = DEGREE
      ENDIF
!
      IF (DEGCNV .GE. 360.) THEN
        DEGCNV = MOD (DEGCNV, 360.)
      ELSE IF (DEGCNV .LT. 0.) THEN
        DEGCNV = MOD (DEGCNV, 360.) + 360.
      ELSE
!       DEGCNV between 0 and 360; do nothing
      ENDIF
!
!
!     *** end of subroutine DEGCNV ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      REAL FUNCTION ANGRAD (DEGREE)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'swcomm3.inc'                                               30.74
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees to radians
!
!  3. METHOD
!
!       ANGRAD = DEGREE * PI / 180
!
!  4. PARAMETERLIST
!
!       ANGRAD      radians
!       DEGREE      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[radian] = ANGLE[degrees} * PI / 180
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGRAD')
!
      ANGRAD = DEGREE * PI / 180.
!
!
!     *** end of subroutine ANGRAD ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      REAL FUNCTION ANGDEG (RADIAN)
!                                                                      *
!***********************************************************************
!
      INCLUDE 'swcomm3.inc'                                               30.74
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform radians to degrees
!
!  3. METHOD
!
!       ANGDEG = RADIAN * 180 / PI
!
!  4. PARAMETERLIST
!
!       RADIAN      radians
!       ANGDEG      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[degrees] = ANGLE[radians} * 180 / PI
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGDEG')
!
      ANGDEG = RADIAN * 180. / PI
!
!
!     *** end of subroutine ANGDEG ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE HSOBND (AC2   ,SPCSIG,HSIBC ,KGRPNT)                     40.00
!                                                                      *
!***********************************************************************
!
      USE M_PARALL                                                        40.31
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!     30.72: IJsbrand Haagsma
!     32.01: Roeland Ris
!     30.70: Nico Booij
!     40.00: Nico Booij
!
!  1. Updates
!
!     32.01, Sep. 97: new for SWAN
!     30.72, Jan. 98: Changed number of elements for HSI to MCGRD
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: structure scheme corrected
!     40.00, Mar. 98: integration method changed (as in SNEXTI)
!                     structure corrected
!
!  2. Purpose
!
!     Compare computed significant wave height with the value of
!     the significant wave height as predescribed by the user. If
!     the values differ more than e.g. 10 % give an error message
!     and the gridpoints where the error has been located
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: input  Relative frequencies in computational domain in      30.72
!                    sigma-space                                          30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!       REALS:
!       ------
!       AC2        action density
!       HSI        significant wave height at boundary (using SWAN
!                  resolution (has thus not to be equal to the WAVEC
!                  significant wave height )
!       ETOT       total energy in a gridpoint
!       DS         increment in frequency space
!       DDIR       increment in directional space
!       HSC        computed wave height after SWAN computation
!       EFTAIL     contribution of tail to spectrum
!
!       INTEGERS:
!       ---------
!       KGRPNT     values of grid indices
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       TRACE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ------------------------------------------------------------------
!     for all computational grid points do                                30.70
!         if HSI is non-zero
!         then compute Hs from action density array
!              if relative difference is large than HSRERR
!              then write error message
!    -------------------------------------------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      REAL      AC2(MDC,MSC,MCGRD) ,HSIBC(MCGRD)                          30.72
!
      REAL      ETOT  ,HSC                                                40.00
!
      INTEGER   ID    ,IS     ,IX     ,IY    ,INDX
!
      LOGICAL   HSRR
!
      INTEGER   KGRPNT(MXC,MYC)
!
      SAVE IENT, HSRR
      DATA IENT/0/, HSRR/.TRUE./
      CALL STRACE (IENT, 'HSOBND')
!
!     *** initializing ***
!
      HSRR = .TRUE.                                                       40.03
!
      DO IY = MYC, 1, -1
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
          IF ( HSIBC(INDX) .GT. 1.E-25 ) THEN
!           *** compute Hs for boundary point (without tail) ***
            ETOT  = 0.
            DO ID = 1, MDC
              DO IS = 1, MSC                                              40.00
                ETOT = ETOT + SPCSIG(IS)**2 * AC2(ID,IS,INDX)             40.00
              ENDDO
            ENDDO
            IF (ETOT .GT. 1.E-8) THEN
              HSC = 4. * SQRT(ETOT*FRINTF*DDIR)
            ELSE
              HSC = 0.
            ENDIF
            HSREL = (HSIBC(INDX) - HSC) / HSIBC(INDX)                     40.00
            IF (HSREL .GT. HSRERR) THEN
              IF ( HSRR ) THEN
                WRITE (PRINTF,*) ' ** WARNING : ',
     &             'Differences in wave height at the boundary'
                WRITE (PRINTF,802) HSRERR
 802            FORMAT (' Relative difference between input and ',
     &          'computation >= ', F6.2)
                WRITE (PRINTF,*) '                        Hs[m]',
     &                           '      Hs[m]      Hs[-]'
                WRITE (PRINTF,*) '    ix    iy  index   (input)',
     &                           ' (computed) (relative)'
                WRITE (PRINTF,*) ' ----------------------------',
     &                           '----------------------'
                HSRR = .FALSE.
              ENDIF
              WRITE (PRINTF,'(3(1x,I5),3(1x,F10.2))')
     &                 IX+MXF-1, IY+MYF-1, INDX, HSIBC(INDX), HSC, HSREL
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      WRITE(PRINTF,*)
!
      IF ( ITEST .GE. 150 ) THEN
        WRITE(PRINTF,*) 'Values of wave height at boundary (HSOBND)'
        WRITE(PRINTF,*) '------------------------------------------'
        DO IY = MYC, 1, -1
          WRITE (PRINTF,'(13F8.3)') ( HSIBC(KGRPNT(IX,IY)), IX=1 , MXC)
        ENDDO
      ENDIF
!
!     *** end of subroutine HSOBND ***
!
      RETURN
      END
!
!********************************************************************
!                                                                   *
      SUBROUTINE SETUPP (KGRPNT, MSTPDA, SETPDA, AC2, DEP2, DEPSAV,
     &                   SETUP2, WFORCX, WFORCY, XCGRID, YCGRID,
     &                   SPCSIG, SPCDIR, ITSW, ITER, UPPERI, LOPERI)      30.82
!                                                                   *
!********************************************************************
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
!  add pass block - fyshi

      USE PASS
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
!
!  0. Authors
!
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     31.03: Annette Kieftenburg
!     31.04: Nico Booij
!     32.01: Roeland Ris
!     32.03: IJsbrand Haagsma
!     34.01: Jeroen Adema
!
!  1. Updates
!
!     32.01, Sept 97: New Subroutine
!     32.03, Feb. 98: Comma added in FORMAT to prevent compilation error
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: transformation of radiation stress in 1D case
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Dec. 98: Argument list KSCIP1 adjusted
!     34.01, Feb. 99: Introducing STPNOW
!     30.82, July 99: Corrected argumentlist SETUPP and SETUP2D
!     30.82, July 99: Corrected argumentlist KSCIP1
!
!  2. Purpose
!
!     Computes the forces/(RHO*GRAV) responsible for the SETUP
!     and adds the SETUP to the depth
!
!  3. Method
!
!     The wave-induced setup is calculated for the one-dimensional
!     mode of SWAN using the following equation:
!
!        d Sxx                d eta
!        ----- +  ( d + eta ) ----- = 0
!         d x                  d x
!
!     This equation is numerically integrated in this subroutine
!     using an explicit numerical scheme in geograhical space.
!
!  4. Argument variables
!
!     AC2       input  (Nonstationary case) action density                31.03
!                                  as function of D,S,X,Y at time T+DT    31.03
!     DEPSAV    input    depth following from bottom and water level
!     DEP2      i/o      total depth, including SETUP
!                        on entry: includes previous estimate of SETUP
!                        on exit:  includes new estimate of SETUP
!     ITER      input    iteration counter
!     KGRPNT    input    indirect addresses for grid points
!     LOPERI                                                              30.82
!     MSTPDA    input    number of (aux.) data per grid point
!                            value is set at 10 in swancom1.ftn
!     SETPDA    i/o      (aux.) data for computation of Setup
!                        1: Depth, 2: (previous estimate of) Setup,
!                        3: x-comp of force, 4: y-comp of force,
!                        5: rad. stress comp. RSxx, 6: RSxy,
!                        7: RSyy
!                        SETPDA(*,*,5..MSTPDA) is used as work array
!     SETUP2    output   SETUP in grid points, using indirect addresses
!     SPCDIR    input    (*,1); spectral directions (radians)             30.82
!                        (*,2); cosine of spectral directions             30.82
!                        (*,3); sine of spectral directions               30.82
!                        (*,4); cosine^2 of spectral directions           30.82
!                        (*,5); cosine*sine of spectral directions        30.82
!                        (*,6); sine^2 of spectral directions             30.82
!     SPCSIG    input    Relative frequencies in computational domain     30.72
!                        in sigma-space                                   30.72
!     UPPERI                                                              30.82
!     WFORCX    output   Force x-component
!     WFORCY    output   Force y-component
!     XCGRID    input    Coordinates of computational grid in x-direction 30.72
!     YCGRID    input    Coordinates of computational grid in y-direction 30.72
!
      INTEGER ITER, ITSW, MSTPDA, KGRPNT(MXC,MYC)
!
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    DEPSAV(MCGRD)
      REAL    LOPERI(*)                                                   30.82
      REAL    SETPDA(MXC,MYC,MSTPDA)
      REAL    SETUP2(MCGRD)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    UPPERI(*)                                                   30.82
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)                            30.72
      REAL    WFORCX(MCGRD)  , WFORCY(MCGRD)                              31.04
!
!  5. Parameter variables
!
!  6. Local variables
!
!     CG          group velocity
!     CK          CGO*KWAVE
!     DDET        determinant
!     DEPMAX      maximum depth
!     DEPLOC      local depth
!     DIX         di/dx
!     DIY         di/dy
!     DJX         dj/dx
!     DJY         dj/dy
!     DP1         depth in point i in 1-D case
!     DP2         depth in point i+1 in 1-D case
!     DS2         square of mesh length in x- or y-direction
!     DXI         dx/di
!     DXJ         dx/dj
!     DYI         dy/di
!     DYJ         dy/dj
!     ELOC        local energy
!     ETA1        setup in  point i in 1-D case
!     ETA2        setup in  point i+1 in 1-D case
!     ID          counter in directional space
!     IDPMXX      identifier for x-coordinate of location with maximum depth
!     IDPMXY      identifier for y-coordinate of location with maximum depth
!     IENT        number of entries
!     INDX        address of grid point
!     INDXL       address of neighbouring grid point
!     IS          counter in frequency space
!     IX          counter in x-direction
!     IXLO        counter in x-direction for neighbouring grid point
!     IXUP        counter in x-direction for neighbouring grid point
!     IY          counter in y-direction
!     IYLO        counter in y-direction for neighbouring grid point
!     IYUP        counter in y-direction for neighbouring grid point
!     K           wavenumber
!     LINK        counter for neighbouring grid points
!     N           CGroup/CPhase
!     ND          derivative of N with respect to depth
!     NEIGHB      boolean variable indicating whether neighbouring point is wet
!     RRDI        1/number of steps in i-direction
!     RRDJ        1/number of steps in j-direction
!     RSXX        xx-component of the radiation stress
!     RSXXI       derivative of RSXX in i-direction
!     RSXXJ       derivative of RSXX in j-direction
!     RSXY        xy-component of the radiation stress
!     RSXYI       derivative of RSXY in i-direction
!     RSXYJ       derivative of RSXY in j-direction
!     RSYY        yy-component of the radiation stress
!     RSYYI       derivative of RSYY in i-direction
!     RSYYJ       derivative of RSYY in j-direction
!     S_UPCOR     total setupcorrection factor (user defined and S_UPDP)
!     S_UPDP      setup at location with maximum depth, before correction
!     SIG         dummy variable for frequency
!     SXX1        radiation stress in  point i in 1-D case
!     SXX2        radiation stress in  point i+1 in 1-D case
!
      INTEGER  IDPMXX, IDPMXY,                                            31.03
     &         ID, IENT, INDX, INDXL, IS, IX,
     &         IXLO, IXUP, IY, IYLO, IYUP, LINK
!
      REAL     CK, DDET, DEPLOC, DEPMAX, DIX, DIY, DJX, DJY,              31.04
     &         DP1, DP2, DS2, DXI, DXJ, DYI, DYJ, ELOC, ETA1, ETA2,
     &         RRDI, RRDJ, RSXX,  RSXXI, RSXXJ, RSXY,
     &         RSXYI, RSXYJ, RSYY, RSYYI, RSYYJ,
     &         S_UPCOR,                                                   30.82
     &         S_UPDP,                                                    31.03
     &         SXX1  ,SXX2
      REAL     CG(1), K(1), N(1), ND(1), SIG(1)                           30.82
!
      LOGICAL  NEIGHB
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!     STRACE
!     KSCIP1      Calculates KWAVE, CGO
!     SETUP2D     Computation of SETUP, the change of waterlevel by waves.
!                 A Poisson equation is solved in general coordinates
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SWCOMPU
!
! 10. Error messages
!
!     If setup is unequal to zero in dry point
!
! 11. Remarks
!
! 12. Structure
!
!     ---------------------------------------------------------
!     For all grid points do
!         If depth > DEPMIN
!         Then Integrate over spectrum to find RSxx, RSxy, RSyy
!     ---------------------------------------------------------
!     If one-dimensional mode of SWAN
!     Then Calculate Setup in all grid points
!     Else Call SETUP2 to compute setup in all grid points
!     ---------------------------------------------------------
!     S_updp is setup in (first) deepest point
!     Add user defined correction to setup                                30.82
!     For all grid points do
!         If dep2 > DEPMIN
!            copy Setup - S_upcor to SETUP2
!         If dep2 > DEPMIN
!            compute new value for DEP2
!     ---------------------------------------------------------
!     For all grid points do
!         If depth < DEPMIN
!         Then If water level + setup in neighbouring point above
!                   bottom level in this point
!              Then make depth equal to neighbouring water level
!                   + SETUP - bottom level in point itself
!     ---------------------------------------------------------
!
! 13. Source text
!
!***********************************************************************

      SAVE     IENT
      DATA     IENT /0/
      CALL STRACE (IENT, 'SETUPP')
!
      DEPMAX = 0.                                                         31.03
      IDPMXX = 0                                                          31.03
      IDPMXY = 0                                                          31.03
!                                                                         31.03
!     Initializing SETPDA, WFORCX and WFORCY arrays                       31.03
!                                                                         31.03
      DO IY = 1, MYC                                                      31.03
        DO IX = 1, MXC                                                    31.03
          INDX = KGRPNT(IX,IY)
          SETPDA(IX,IY,1) = -9.                                           31.03
          SETPDA(IX,IY,2) = 0.                                            31.03
          SETPDA(IX,IY,3) = 0.                                            31.03
          SETPDA(IX,IY,4) = 0.                                            31.03
          SETPDA(IX,IY,5) = 0.                                            31.03
          SETPDA(IX,IY,6) = 0.                                            31.03
          SETPDA(IX,IY,7) = 0.                                            31.03
          WFORCX(INDX) = 0.                                               31.03
          WFORCY(INDX) = 0.                                               31.03
         ENDDO                                                            31.03
      ENDDO                                                               31.03
!
      DO IY = 1, MYC
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
          IF (INDX.GT.1) THEN
            IF (DEP2(INDX).GT.DEPMIN) THEN
!
!             In dry points, even after inundation, setup = 0,            31.03
!             so there is no need to set SETPDA to the last estimate.     31.03
!
              SETPDA(IX,IY,1) = DEP2(INDX)                                31.03
              SETPDA(IX,IY,2) = SETUP2(INDX)                              31.03
!                                                                         31.03
!             Seek deepest point.                                         31.03
!                                                                         31.03
              IF (DEPSAV(INDX).GT.DEPMAX) THEN                            31.03
                DEPMAX = DEPSAV(INDX)                                     31.03
                IDPMXX = IX                                               31.03
                IDPMXY = IY                                               31.03
              ENDIF                                                       31.03
!
!             compute radiation stress components RSXX, RSXY and RSYY
!
              RSXX = 0.
              RSXY = 0.
              RSYY = 0.
              DEPLOC = SETPDA(IX,IY,1)                                    31.03
              DO IS = 1, MSC
                SIG(1) = SPCSIG(IS)                                       30.82
                CALL KSCIP1 (1,SIG,DEPLOC,K,CG,N,ND)                      30.82
                CK = CG(1) * K(1)                                         30.82
                DO ID = 1, MDC
                  ELOC = SIG(1) * AC2(ID,IS,INDX)                         30.82
!                                  -                                      31.03
!                                  |{cos(Theta)}^2         for i = 4      31.03
!                 SPCDIR(ID,i) is <| sin(Theta)cos(Theta)  for i = 5      31.03
!                                  |{sin(Theta)}^2         for i = 6      31.03
!                                  -                                      31.03
!                                                                         31.03
                  RSXX = RSXX + (CK*SPCDIR(ID,4)+CK - SIG(1)/2.) * ELOC   30.82
                  RSXY = RSXY + CK*SPCDIR(ID,5) * ELOC                    31.03
                  RSYY = RSYY + (CK*SPCDIR(ID,6)+CK - SIG(1)/2.) * ELOC   30.82
                ENDDO
              ENDDO
!
!             store radiation stress components in array SETPDA
!
!             DDIR   is width of directional band
!             FRINTF is frequency integration factor df/f
!
              IF (ONED) THEN                                              30.70
!               transform to computational direction
                SETPDA(IX,IY,5) = DDIR * FRINTF *                         31.04
     &                           ((COSPC*RSXX + SINPC*RSXY) * COSPC +     30.70
     &                            (COSPC*RSXY + SINPC*RSYY) * SINPC)      30.70
              ELSE                                                        30.70
                SETPDA(IX,IY,5) = RSXX * DDIR * FRINTF
                SETPDA(IX,IY,6) = RSXY * DDIR * FRINTF
                SETPDA(IX,IY,7) = RSYY * DDIR * FRINTF
              ENDIF                                                       30.70
            ENDIF                                                         31.03
          ENDIF
        ENDDO
      ENDDO

!
      IF ( ONED ) THEN
!
!       *** compute on the basis of the radiation stresses the setup ***
!       *** output is new setup = SETPDA(1,1,2)
!
        DO IY = 1, MYC
!         *** boundary condition ***
          SETPDA(1,IY,2) = 0.
          ETA2 = 0.
          DO IX = 1, MXC-1
            DP2  = SETPDA(IX+1,IY,1)
            IF ( DP2 .GT. 0. ) THEN
              DP1  = SETPDA(IX  ,IY,1)
              ETA1 = SETPDA(IX  ,IY,2)
              SXX1 = SETPDA(IX  ,IY,5)                                    31.04
              SXX2 = SETPDA(IX+1,IY,5)                                    31.04
              ETA2 = ETA1 + ( SXX1 - SXX2 ) / ( 0.5 * ( DP2 + DP1 ) )
              SETPDA(IX+1,IY,2) = ETA2
            ENDIF
          ENDDO
        ENDDO
!
      ELSE
!
!       compute forces by taking derivative of radiation stress
!
        DO IY = 1, MYC
          DO IX = 1, MXC
            INDX = KGRPNT(IX,IY)                                          31.03
            DEPLOC = SETPDA(IX,IY,1)
            IF (DEPLOC.LE.DEPMIN
     &          .OR. INDX.EQ.1) THEN
              GOTO 700
            ENDIF
            IF (IX.EQ.1) THEN
              IXLO = 1
              IXUP = 2
            ELSE IF (IX.EQ.MXC) THEN
              IXLO = MXC-1
              IXUP = MXC
            ELSE
              IXLO = IX-1
              IXUP = IX+1
            ENDIF
            IF (SETPDA(IXLO,IY,1).LE.DEPMIN) IXLO = IX                    31.03
            IF (SETPDA(IXUP,IY,1).LE.DEPMIN) IXUP = IX                    31.03
            IF (IXLO.EQ.IXUP) THEN
              RRDI = 1e-20
            ELSE
              RRDI = 1. / REAL(IXUP-IXLO)
            ENDIF
            IF (IY.EQ.1) THEN                                             31.03
              IYLO = 1
              IYUP = 2
            ELSE IF (IY.EQ.MYC) THEN
              IYLO = MYC-1
              IYUP = MYC
            ELSE
              IYLO = IY-1
              IYUP = IY+1
            ENDIF
            IF (SETPDA(IX,IYLO,1).LE.DEPMIN) IYLO = IY                    31.03
            IF (SETPDA(IX,IYUP,1).LE.DEPMIN) IYUP = IY                    31.03
            IF (IYLO.EQ.IYUP) THEN
              RRDJ = 1e-20
            ELSE
              RRDJ = 1. / REAL(IYUP-IYLO)
            ENDIF
!
!           determine (x,y) derivatives w.r.t. i and j
!
            DXI = RRDI * (XCGRID(IXUP,IY)-XCGRID(IXLO,IY))
            DYI = RRDI * (YCGRID(IXUP,IY)-YCGRID(IXLO,IY))
            DXJ = RRDJ * (XCGRID(IX,IYUP)-XCGRID(IX,IYLO))
            DYJ = RRDJ * (YCGRID(IX,IYUP)-YCGRID(IX,IYLO))
!
            RSXXI = RRDI * (SETPDA(IXUP,IY,5)-SETPDA(IXLO,IY,5))
            RSXXJ = RRDJ * (SETPDA(IX,IYUP,5)-SETPDA(IX,IYLO,5))          31.03
            RSXYI = RRDI * (SETPDA(IXUP,IY,6)-SETPDA(IXLO,IY,6))
            RSXYJ = RRDJ * (SETPDA(IX,IYUP,6)-SETPDA(IX,IYLO,6))          31.03
            RSYYI = RRDI * (SETPDA(IXUP,IY,7)-SETPDA(IXLO,IY,7))
            RSYYJ = RRDJ * (SETPDA(IX,IYUP,7)-SETPDA(IX,IYLO,7))          31.03
!
            IF (IXLO.EQ.IXUP.AND.IYLO.EQ.IYUP) THEN                       31.03
!             point surrounded by dry points                              31.03
              DIX  = 0.                                                   31.04
              DIY  = 0.                                                   31.04
              DJX  = 0.                                                   31.04
              DJY  = 0.                                                   31.04
            ELSE IF (IXLO.EQ.IXUP) THEN                                   31.03
!             no forces in i-direction                                    31.03
              DS2  = DXJ**2 + DYJ**2                                      31.04
              DIX  = 0.                                                   31.04
              DIY  = 0.                                                   31.04
              DJX  = DXJ/DS2                                              31.04
              DJY  = DYJ/DS2                                              31.04
            ELSE IF (IYLO.EQ.IYUP) THEN                                   31.03
!             no forces in j-direction                                    31.03
              DS2  = DXI**2 + DYI**2                                      31.04
              DIX  = DXI/DS2                                              31.04
              DIY  = DYI/DS2                                              31.04
              DJX  = 0.                                                   31.04
              DJY  = 0.                                                   31.04
            ELSE                                                          31.03
!             coefficients for transformation from
!             (i,j)-gradients to (x,y)-gradients
              DDET = DXI*DYJ - DXJ*DYI
              DIX  =  DYJ / DDET
              DIY  = -DXJ / DDET
              DJX  = -DYI / DDET
              DJY  =  DXI / DDET
            ENDIF                                                         31.04
!
!           force: spatial gradients of rad. stress
            SETPDA(IX,IY,3) =
     &             -(RSXXI*DIX + RSXXJ*DJX + RSXYI*DIY + RSXYJ*DJY)       31.03
            SETPDA(IX,IY,4) =
     &             -(RSXYI*DIX + RSXYJ*DJX + RSYYI*DIY + RSYYJ*DJY)       31.03
!
!           store force values in array for output                        31.04
            WFORCX(INDX) = SETPDA(IX,IY,3)
            WFORCY(INDX) = SETPDA(IX,IY,4)
! pass wave forcing to master (turns out they are in real x y directions)
! -fyshi
            Pass_Wave_Fx(IX,IY)=SETPDA(IX,IY,3)
            Pass_Wave_Fy(IX,IY)=SETPDA(IX,IY,4)
!
 700      CONTINUE
          ENDDO                                                               31.03
        ENDDO
!
        write(*,*)'Passing wave force to master'
!
!       call Setup2d to compute SETUP
!
!       no nesting:  LSETUP = 1
!       nesting:     LSETUP = 2
!
         CALL SETUP2D (XCGRID, YCGRID, SETPDA(1,1,3),SETPDA(1,1,4),
     &                 SETPDA(1,1,1),SETPDA(1,1,2), UPPERI, LOPERI,       30.82
     &                 MSTPDA-4, SETPDA(1,1,5), ITSW, ITER)               31.04
         IF (STPNOW()) RETURN                                             34.01
      ENDIF
!
      IF (LSETUP.EQ.1) THEN                                               31.04
!       Set setup to 0 for deepest point. (This is allowed because the    31.03
!       solution of a Poisson equation + constant is again a solution of  31.03
!       the same Poisson equation)                                        31.03
        S_UPDP = SETPDA(IDPMXX, IDPMXY,2)                                 31.03
        S_UPCOR = S_UPDP - PSETUP(2)                                      30.82
        DO IY = 1, MYC                                                    31.03
          DO IX = 1, MXC                                                  31.03
             INDX = KGRPNT(IX,IY)                                         31.03
             IF (INDX.GT.1) THEN                                          31.03
               IF (DEP2(INDX).GT.DEPMIN) THEN                             31.03
                 SETUP2(INDX) = SETPDA(IX,IY,2) - S_UPCOR                 30.82
               ELSE
                 SETUP2(INDX) = SETPDA(IX,IY,2)                           31.03
                 IF (ABS(SETUP2(INDX)).GE.1e-7) THEN                      31.03
                   WRITE (PRINTF,*) 'Setup =', SETUP2(INDX),              31.03
     &              'in dry point (', IX,',', IY,') !'                    31.03
                 ENDIF                                                    31.03
               END IF                                                     31.03
             END IF                                                       31.03
          ENDDO                                                           31.03
        ENDDO                                                             31.03
      ENDIF                                                               31.04
!
!     include computed SETUP in depth
!
      DO IY = 1, MYC
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
          IF (INDX.GT.1) THEN
              DEP2(INDX) = DEPSAV(INDX) + SETUP2(INDX)
          ENDIF
        ENDDO
      ENDDO
!
!     check whether dry points should be inundated
!
      DO IY = 1, MYC
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
!         Note:    KGRPNT(.,.) = 1 means a permanently dry point!         31.03
          IF (INDX.GT.1) THEN
            IF (DEP2(INDX).LE.DEPMIN) THEN
              DO LINK = 1, 4
                NEIGHB = .TRUE.
                IF (LINK.EQ.1) THEN
                  IF (IX.EQ.1) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX-1,IY)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ELSE IF (LINK.EQ.2) THEN
                  IF (IY.EQ.1) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX,IY-1)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ELSE IF (LINK.EQ.3) THEN
                  IF (IX.EQ.MXC) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX+1,IY)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ELSE IF (LINK.EQ.4) THEN
                  IF (IY.EQ.MYC) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX,IY+1)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ENDIF
                IF (NEIGHB) THEN                                          31.03
                  IF (DEPSAV(INDX) + SETUP2(INDXL) .GT. DEPMIN) THEN      31.03
                    SETUP2(INDX) = SETUP2(INDXL)                          31.03
                    DEP2(INDX) = DEPSAV(INDX) + SETUP2(INDXL)             31.03
                  ENDIF                                                   31.03
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
      RETURN
!     end of subroutine SETUPP
      END
!
!************************************************************************
!                                                                       *
      SUBROUTINE SETUP2D ( XCGRID, YCGRID, WFRCX, WFRCY, DEPTH,           31.03
     &                     SETUP, UPPERI, LOPERI,                         30.82
     &                     NWKARR, WKARR, ITSW, ITER)                     31.03
!                                                                       *
!************************************************************************
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
      IMPLICIT NONE
!
      INCLUDE 'swcomm3.inc'                                               30.74
!
!  0. Authors
!
!     31.03  Annette Kieftenburg
!     34.01: Jeroen Adema
!
!  1. Updates
!
!     34.01, Feb. 99: Introducing STPNOW
!     30.82, July 99: Corrected argumentlist SETUP2D and SWSOLV
!
!  2. Purpose
!
!     Computation of SETUP, the change of waterlevel by waves.
!     A Poisson equation is solved in general coordinates
!
!  3. Method
!
!  4. Argument variables
!
!     DEPTH     input   Depth
!     ITER      input   Iteration number
!     LOPERI                                                              30.82
!     NWKARR    input   Dimension for work array
!     SETUP     output  Unknown set-up; to be computed indirect addressed
!     UPPERI                                                              30.82
!     WFRCX     input   force x-component
!     WFRCY     input   force y-component
!     WKARR             work array
!     XCGRID    input   x-coordinates
!     YCGRID    input   y-coordinates
!
      INTEGER  ITER, ITSW, NWKARR                                         31.03
      REAL DEPTH(1:MXC,1:MYC),                                            31.03
     &     LOPERI(*),                                                     30.82
     &     SETUP(1:MXC,1:MYC),                                            31.03
     &     UPPERI(*),                                                     30.82
     &     WFRCX(1:MXC,1:MYC),                                            31.03
     &     WFRCY(1:MXC,1:MYC),                                            31.03
     &     WKARR(1:MXC*MYC*NWKARR),
     &     XCGRID(1:MXC,1:MYC),                                           31.03
     &     YCGRID(1:MXC,1:MYC)                                            31.03
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ALPHAD    Direction index of integration.
!     I         General loop variable
!     IENT      Number of entries
!     IPCTC     Starting address of array CTC    in work array WKARR
!     IPCVA     Starting address of array CVA    in work array WKARR
!     IPCVC     Starting address of array CVC    in work array WKARR
!     IPDTSUM   Starting address of array DTSUM  in work array WKARR
!     IPJCTA    Starting address of array JCTA   in work array WKARR
!     IPMATRIX  Starting address of array MATRIX in work array WKARR
!     IPRHSIDE  Starting address of array RHSIDE in work array WKARR
!     NPOINT    Number of points MXC*MYC
!
      INTEGER ALPHAD, I, IENT, IPCTC, IPCVA, IPCVC,  IPDTSUM, IPJCTA,
     &        IPMATRIX, IPRHSIDE, NPOINT
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!     SWCOVA2D  Computes covariant base vectors in integration points
!               two-dimensional case
!     SWJCTA2D  Computes the Jacobian of the transformation, sqrt(g),
!               * contra variant base vectors in
!               integration point two-dimensional case
!     SWTRAD2D  Computes contribution of diffusion term in R2 for
!               a transport equation per integration point
!               Compute righthandside
!     SWDISDT2  Distributes diffusion term for tranport equation in R2
!     SWESSBC   Puts essential boundary conditions in matrix
!     SWSOLV    Solves system of equations
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SETUPP    Computes the forces/(RHO*GRAV) responsible for the SETUP  31.03
!               and adds the SETUP to the depth                           31.03
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
! ======================================================================
      SAVE     IENT
      DATA     IENT /0/
      CALL STRACE (IENT, 'SETUP2D')
!
      NPOINT   = MXC*MYC
      IPMATRIX = 1
      IPRHSIDE = IPMATRIX + NPOINT*9
      IPCVA    = IPRHSIDE + NPOINT
      IPJCTA   = IPCVA    + NPOINT*8
      IPCVC    = IPJCTA   + NPOINT*4
      IPCTC    = IPCVC    + NPOINT*4
      IPDTSUM  = IPCTC    + NPOINT*4
!
!     --- covariant base vectors
!
      CALL SWCOVA2D( MXC, MYC, XCGRID, YCGRID, WKARR(IPCVA) )
!
!     --- jacobian times contravariant base vectors
!
      CALL SWJCTA2D( MXC, MYC, WKARR(IPCVA), WKARR(IPJCTA)  )
!
!     --- initialize  matrix and righthandside

      DO I = 1, 9*NPOINT
          WKARR(IPMATRIX-1+I) = 0E0
      END DO
      DO I = 1, NPOINT
          WKARR(IPRHSIDE-1+I) = 0E0
      END DO
!
!     --- initialize  matrix and righthandside

      DO I = 1, 9*NPOINT
          WKARR(IPMATRIX-1+I) = 0E0
      END DO
      DO I = 1, NPOINT
          WKARR(IPRHSIDE-1+I) = 0E0
      END DO
!
!     --- build matrix and righthandside
!
      DO ALPHAD = 1,2
!
         CALL SWTRAD2D( MXC, MYC, WFRCX, WFRCY,                           31.04
     &                  DEPMIN, ALPHAD, DEPTH,                            31.04
     &                  WKARR(IPCVA), WKARR(IPJCTA),WKARR(IPCVC),
     &                  WKARR(IPCTC), WKARR(IPDTSUM),
     &                  WKARR(IPRHSIDE) )
!
         CALL SWDISDT2( MXC, MYC, DEPTH, DEPMIN, ALPHAD,                  31.04
     &                  WKARR(IPMATRIX), WKARR(IPDTSUM) )
!
      END DO
!
!     --- essential boundary conditions
!
      IF ( LSETUP .EQ. 2) THEN
!
         CALL SWESSBC(MXC, MYC, WKARR(IPMATRIX), WKARR(IPRHSIDE),
     &                SETUP )                                             31.04
!
      END IF
!
!     --- solve system of equations  {WKARR(IPCVA) is used as work array}
!
      CALL SWSOLV ( WKARR(IPMATRIX), WKARR(IPRHSIDE), SETUP,              31.04
     &              NPOINT, WKARR(IPCVA),                                 31.04
     &              NWKARR-11, ITSW, ITER,
     &              UPPERI, LOPERI)                                       30.82
      IF (STPNOW()) RETURN                                                34.01
      RETURN                                                              31.04
      END
!
!*****************************************************************
!                                                                *
      SUBROUTINE CHGBAS (X1, X2, PERIOD, Y1, Y2, N1, N2,
     &                   ITEST, PRTEST)
!                                                                *
!*****************************************************************
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  G. van Vledder, N. Booij                       |
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
!  0. Update history
!
!       ver 20.48: also accomodates periodic variables such as directions
!
!  1. Purpose
!
!       change x-basis of a discretized y-function
!
!  2. Method
!
!     A piecewise constant representation of the functions is assumed
!
!     first boundaries of a cell in X1 are determined
!     then it is determined whether there are overlaps with cells
!     in X2. if so Y1*common length is added to Y2
!     Finally Y2 values are divided by cell lengths
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
!     X1       i    ra   x-coordinates of input grid
!     X2       i    ra   x-coordinates of output grid
!     PERIOD   i    r    period, i.e. x-axis is periodic if period>0
!                        e.g. spectral directions
!     Y1       i    ra   function values of input grid
!     Y2       o    ra   function values of output grid
!     N1       i    i    number of x-values of input grid
!     N2       i    i    number of x-values of output grid
!
!  4. Subroutines used
!
!     ---
!
!  5. Error messages
!
!  6. Remarks
!
!       Cell boundaries in X1 are: X1A and X1B
!       X2 is assumed to be monotonically increasing; this is checked
!       X1 is assumed to be monotonous but not necessarily increasing
!
!  7. Structure
!
!       ------------------------------------------------------------------
!       Make all values of Y2 = 0
!       For each cell in X1 do
!           determine boundaries of cell in X1
!           --------------------------------------------------------------
!           For each cell in X2 do
!               determine overlap with cell in X1; limits: RLOW and RUPP
!               add to Y2: Y1 * length of overlapping interval
!       ------------------------------------------------------------------
!       For each cell in X2 do
!           divide Y2 value by cell length
!       ------------------------------------------------------------------
!
!  8. Source text
!
      INTEGER  I1, I2, N1, N2, ITEST, PRTEST
      REAL     X1(N1), Y1(N1), X2(N2), Y2(N2), PERIOD
      REAL     X1A, X1B, X2A, X2B, RLOW, RUPP
      LOGICAL  TWICE
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'CHGBAS')
!
!     initialize output data
!
      DO I2 = 1, N2
        Y2(I2) = 0.
      ENDDO
      DO I2 = 2, N2
        IF (X2(I2).LE.X2(I2-1))
     &    CALL MSGERR (2, 'subr. CHGBAS: values of X2 not increasing')
      ENDDO
!     boundaries of the range in X2
      X2LO  = 1.5 * X2(1)  - 0.5 * X2(2)
      X2HI  = 1.5 * X2(N2) - 0.5 * X2(N2-1)
      TWICE = .FALSE.
!
!     loop over cells in X1
!
      DO 300 I1 = 1, N1
        IF (ABS(Y1(I1)) .LT. 1.E-20) GOTO 300
!
!       determine cell boundaries in X1
!
        IF (I1.EQ.1) THEN
          X1A = 1.5 * X1(1) - 0.5 * X1(2)
        ELSE
          X1A = 0.5 * (X1(I1) + X1(I1-1))
        ENDIF

        IF (I1.EQ.N1) THEN
          X1B = 1.5 * X1(N1) - 0.5 * X1(N1-1)
        ELSE
          X1B = 0.5 * (X1(I1) + X1(I1+1))
        ENDIF
!
!       swap X1A and X1B if X1A > X1B
!
        IF (X1A.GT.X1B) THEN
          RR  = X1A
          X1A = X1B
          X1B = RR
        ENDIF

        IF (PERIOD.LE.0.) THEN
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
        ELSE
!         X is periodic; move interval in X1 if necessary
          TWICE = .FALSE.
          IADD = 0
  60      IF (X1A.GT.X2HI) THEN
            X1A = X1A - PERIOD
            X1B = X1B - PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)
     &         CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 60
          ENDIF
  70      IF (X1B.LT.X2LO) THEN
            X1A = X1A + PERIOD
            X1B = X1B + PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)
     &           CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 70
          ENDIF
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1A.LT.X2LO .AND. X1A+PERIOD.LT.X2HI) TWICE = .TRUE.
        ENDIF
!
!       loop over cells in X2
!
 100    DO 200 I2 = 1, N2

          IF (I2.EQ.1) THEN
            X2A = X2LO
          ELSE
            X2A = 0.5 * (X2(I2) + X2(I2-1))
          ENDIF

          IF (I2.EQ.N2) THEN
            X2B = X2HI
          ELSE
            X2B = 0.5 * (X2(I2) + X2(I2+1))
          ENDIF
!
!         (RLOW,RUPP) is overlapping interval of (X1A,X1B) and (X2A,X2B)
!
          IF (X1A.LT.X2B) THEN
            RLOW = MAX (X1A, X2A)
          ELSE
            GOTO 200
          ENDIF

          IF (X1B.GT.X2A) THEN
            RUPP = MIN (X1B, X2B)
          ELSE
            GOTO 200
          ENDIF

          IF (RUPP.LT.RLOW) THEN
            CALL MSGERR (3, 'interpolation error')
            WRITE (PRTEST, 140) I1, X1A, X1B, I2, X2A, X2B
 140        FORMAT (' I, XA, XB ', 2(I3, 2(1X,E12.4)))
          ELSE
            Y2(I2) = Y2(I2) + Y1(I1) * (RUPP-RLOW)
          ENDIF
 200    CONTINUE
!
!       Cell in X1 covers both ends of sector boundary
        IF (TWICE) THEN
          X1A = X1A + PERIOD
          X1B = X1B + PERIOD
          TWICE = .FALSE.
          GOTO 100
        ENDIF
 300  CONTINUE
!
      DO I2 = 1, N2
        IF (I2.EQ.1) THEN
          CELLEN = X2(2) - X2(1)
        ELSE IF (I2.EQ.N2) THEN
          CELLEN = X2(N2) - X2(N2-1)
        ELSE
          CELLEN = 0.5 * (X2(I2+1) - X2(I2-1))
        ENDIF
!       divide Y2 by cell length
        Y2(I2) = Y2(I2) / CELLEN
      ENDDO
      IF (ITEST.GE.160) THEN
        WRITE (PRTEST, 84) N1, N2
  84    FORMAT (' test CHGBAS ', 2I5)
        WRITE (PRTEST, 85) (X1(II), II = 1, N1)
        WRITE (PRTEST, 85) (Y1(II), II = 1, N1)
        WRITE (PRTEST, 85) (X2(II), II = 1, N2)
        WRITE (PRTEST, 85) (Y2(II), II = 1, N2)
  85    FORMAT (10 (1X,E10.3))
      ENDIF
!
      RETURN
      END
!
!********************************************************************
!                                                                   *
      REAL FUNCTION GAMMA(XX)
!                                                                   *
!********************************************************************
!
!   Updates
!     ver 30.70, Oct 1997 by N.Booij: new subroutine
!
!   Purpose
!     Compute the transcendental function Gamma
!
!   Subroutines used
!     GAMMLN  (Numerical Recipes)
!
      REAL XX, YY, ABIG                                                   40.00
      SAVE IENT, ABIG
      DATA IENT /0/, ABIG /30./
      CALL STRACE (IENT, 'GAMMA')
      YY = GAMMLN(XX)
      IF (YY.GT.ABIG) YY = ABIG
      IF (YY.LT.-ABIG) YY = -ABIG
      GAMMA = EXP(YY)
      RETURN
      END
!********************************************************************
!                                                                   *
      FUNCTION GAMMLN(XX)
!                                                                   *
!********************************************************************
!
!   Method:
!     function is copied from: Press et al., "Numerical Recipes"
!
      DOUBLE PRECISION  COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE WRSPEC (NREF, ACLOC)
!                                                                      *
!***********************************************************************

      USE OUTP_DATA                                                       40.13

      IMPLICIT NONE                                                       40.13

      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm3.inc'                                               30.74
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
!     40.00, 40.13: Nico Booij
!
!  1. UPDATE
!
!     new subroutine, update 40.00
!     40.03, Mar. 00: precision increased; 2 decimals more in output table
!     40.13, July 01: variable format using module OUTP_DATA
!
!  2. PURPOSE
!
!     Writing of action density spectrum in Swan standard format
!
!  3. METHOD
!
!
!  4. Argument variables
!
!       NREF    int    input    unit ref. number of output file
!       ACLOC   real   local    2-D spectrum or source term at one
!                               output location

      INTEGER, INTENT(IN) :: NREF
      REAL, INTENT(IN)    :: ACLOC(1:MDC,1:MSC)

!  5. Parameter variables
!
!  6. Local variables
!
!       ID      counter of spectral directions
!       IS      counter of spectral frequencies
!
      INTEGER :: ID, IS
!
!       EFAC    multiplication factor written to file
!
      REAL    :: EFAC
!
!  7. Common Blocks used
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       determine maximum value of ACLOC
!       if maximum = 0
!       then write 'ZERO' to file
!       else write 'FACTOR'
!            determine multiplication factor, write this to file
!            write values of ACLOC/factor to file
!       ----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0                                           40.13
      IF (LTRACE) CALL STRACE (IENT, 'WRSPEC')
!
!     first determine maximum energy density
      EFAC = 0.
      DO ID = 1, MDC
        DO IS = 1, MSC
          IF (ACLOC(ID,IS).GE.0.) THEN
            EFAC = MAX (EFAC, ACLOC(ID,IS))
          ELSE
            EFAC = MAX (EFAC, 10.*ABS(ACLOC(ID,IS)))
          ENDIF
        ENDDO
      ENDDO
      IF (EFAC .LE. 1.E-10) THEN
        WRITE (NREF, 12) 'ZERO'                                           40.00
  12    FORMAT (A4)
      ELSE
        EFAC = 1.01 * EFAC * 10.**(-DEC_SPEC)                             40.13
!       factor PI/180 introduced to account for change from rad to degr
!       factor 2*PI to account for transition from rad/s to Hz
        WRITE (NREF, 95) EFAC * 2. * PI**2 / 180.
  95    FORMAT ('FACTOR', /, E18.8)                                       40.13
        DO IS = 1, MSC
!         write spectral energy densities to file
          WRITE (NREF, FIX_SPEC) (NINT(ACLOC(ID,IS)/EFAC), ID=1,MDC)      40.13
        ENDDO
      ENDIF
      RETURN
!     end of subroutine WRSPEC
      END
!****************************************************************
!
      SUBROUTINE SWTSTA (ITIMER)
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'timecomm.inc'
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Aug. 02: New subroutine
!
!  2. Purpose
!
!     Start timing
!
!  3. Method
!
!     Get cpu and wall-clock times and store
!
!  4. Argument variables
!
!     ITIMER      number of timer to be used
!
      INTEGER :: ITIMER
!
!  6. Local variables
!
!     C     :     clock count of processor
!     I     :     index in LISTTM, loop variable
!     IFOUND:     index in LISTTM, location of ITIMER
!     IFREE :     index in LISTTM, first free position
!     M     :     maximum clock count
!     R     :     number of clock counts per second
!     TIMER :     current real cpu-time
!     TIMER1:     current cpu-time used
!     TIMER2:     current wall-clock time used
!
      INTEGER          :: I, IFOUND, IFREE
      INTEGER          :: C, R, M
      REAL             :: TIMER
      DOUBLE PRECISION :: TIMER1, TIMER2
!
!  8. Subroutines used
!
!     CPU_TIME         Returns real value from cpu-time clock
!     SYSTEM_CLOCK     Returns integer values from a real-time clock
!
!  9. Subroutines calling
!
!     SWMAIN, SWCOMP, SWOMPU
!
! 12. Structure
!
!     Get and store the cpu and wall-clock times
!
! 13. Source text
!

!
!     --- check whether a valid timer number is given
!
      IF (ITIMER.LE.0 .OR. ITIMER.GT.NSECTM) THEN
         WRITE(PRINTF,*) 'SWTSTA: ITIMER out of range: ',
     &                   ITIMER, 1, NSECTM
         STOP
      END IF
!
!     --- check whether timing for ITIMER was started already,
!         also determine first free location in LISTTM
!
      IFOUND=0
      IFREE =0
      I     =0
 100  IF (I.LT.LASTTM .AND. (IFOUND.EQ.0 .OR. IFREE.EQ.0)) THEN
         I=I+1
         IF (LISTTM(I).EQ.ITIMER) THEN
            IFOUND=I
         END IF
         IF (IFREE.EQ.0 .AND. LISTTM(I).EQ.-1) THEN
            IFREE =I
         END IF
         GOTO 100
      END IF

      IF (IFOUND.EQ.0 .AND. IFREE.EQ.0 .AND. LASTTM.LT.MXTIMR) THEN
         LASTTM=LASTTM+1
         IFREE =LASTTM
      END IF
!
!     --- produce warning if found in the list
!
      IF (IFOUND.GT.0) THEN
         WRITE(PRINTF,*)
     &      'SWTSTA: warning: previous timing for section ',
     &      ITIMER,' not closed properly/will be ignored.'
      END IF
!
!     --- produce error if not found and no free position available
!
      IF (IFOUND.EQ.0 .AND. IFREE.EQ.0) THEN
         WRITE(PRINTF,*)
     &      'SWTSTA: maximum number of simultaneous timers',
     &      ' exceeded:',MXTIMR
         STOP
      END IF
!
!     --- register ITIMER in appropriate location of LISTTM
!
      IF (IFOUND.EQ.0) THEN
         IFOUND=IFREE
      END IF
      LISTTM(IFOUND)=ITIMER
!
!     --- get current cpu/wall-clock time and store in TIMERS
!
      TIMER1=0D0
      CALL CPU_TIME (TIMER)
      TIMER1=DBLE(TIMER)
      CALL SYSTEM_CLOCK (C,R,M)
      TIMER2=DBLE(C)/DBLE(R)

      TIMERS(IFOUND,1)=TIMER1
      TIMERS(IFOUND,2)=TIMER2

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWTSTO (ITIMER)
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'timecomm.inc'
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Aug. 02: New subroutine
!
!  2. Purpose
!
!     Stop timing
!
!  3. Method
!
!     Get cpu and wall-clock times and store
!
!  4. Argument variables
!
!     ITIMER      number of timer to be used
!
      INTEGER :: ITIMER
!
!  6. Local variables
!
!     C     :     clock count of processor
!     I     :     index in LISTTM, loop variable
!     IFOUND:     index in LISTTM, location of ITIMER
!     M     :     maximum clock count
!     R     :     number of clock counts per second
!     TIMER :     current real cpu-time
!     TIMER1:     current cpu-time used
!     TIMER2:     current wall-clock time used
!
      INTEGER          :: I, IFOUND
      INTEGER          :: C, R, M
      REAL             :: TIMER
      DOUBLE PRECISION :: TIMER1, TIMER2
!
!  8. Subroutines used
!
!     CPU_TIME         Returns real value from cpu-time clock
!     SYSTEM_CLOCK     Returns integer values from a real-time clock
!
!  9. Subroutines calling
!
!     SWMAIN, SWCOMP, SWOMPU
!
! 12. Structure
!
!     Get and store the cpu and wall-clock times
!
! 13. Source text
!

!
!     --- check whether a valid timer number is given
!
      IF (ITIMER.LE.0 .OR. ITIMER.GT.NSECTM) THEN
         WRITE(PRINTF,*) 'SWTSTO: ITIMER out of range: ',
     &                   ITIMER, 1, NSECTM
         STOP
      END IF
!
!     --- check whether timing for ITIMER was started already,
!         also determine first free location in LISTTM
!
      IFOUND=0
      I     =0
 100  IF (I.LT.LASTTM .AND. IFOUND.EQ.0) THEN
         I=I+1
         IF (LISTTM(I).EQ.ITIMER) THEN
            IFOUND=I
         END IF
         GOTO 100
      END IF
!
!     --- produce error if not found
!
      IF (IFOUND.EQ.0) THEN
         WRITE(PRINTF,*)
     &      'SWTSTO: section ',ITIMER,' not found',
     &      ' in list of active timings'
         STOP
      END IF
!
!     --- get current cpu/wall-clock time
!
      TIMER1=0D0
      CALL CPU_TIME (TIMER)
      TIMER1=DBLE(TIMER)
      CALL SYSTEM_CLOCK (C,R,M)
      TIMER2=DBLE(C)/DBLE(R)
!
!     --- calculate elapsed time since start of timing,
!         store in appropriate location in DCUMTM,
!         increment number of timings for current section
!
      DCUMTM(ITIMER,1)=DCUMTM(ITIMER,1)+(TIMER1-TIMERS(IFOUND,1))
      DCUMTM(ITIMER,2)=DCUMTM(ITIMER,2)+(TIMER2-TIMERS(IFOUND,2))
      NCUMTM(ITIMER)  =NCUMTM(ITIMER)+1
!
!     --- free appropriate location of LISTTM,
!         adjust last occupied position of LISTTM
!
      IF (IFOUND.GT.0) THEN
         LISTTM(IFOUND)=-1
      END IF
 200  IF (LASTTM.GT.1 .AND. LISTTM(LASTTM).EQ.-1) THEN
         LASTTM=LASTTM-1
         GOTO 200
      END IF
      IF (LISTTM(LASTTM).EQ.-1) LASTTM=0

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWPRTI
!
!****************************************************************
!
      USE M_PARALL                                                        40.31

      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'timecomm.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Aug. 02: New subroutine
!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!
!  2. Purpose
!
!     Print timings info
!
!  6. Local variables
!
!     IDEBUG:     level of timing output requested:
!                 0 - no output for detailed timings
!                 1 - aggregate output for detailed timings
!                 2 - complete output for all detailed timings
!     IENT  :     number of entries
!     J     :     loop counter
!     K     :     loop counter
!     TABLE :     array for computing aggregate cpu- and wallclock-times
!
      INTEGER          :: IENT, J, K, IDEBUG
      DOUBLE PRECISION :: TABLE(30,2)
      PARAMETER (IDEBUG=0)
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWMAIN (in SWANMAIN)
!
! 12. Structure
!
!     Compile table with overview of cpu/wall clock time used in
!     important parts of SWAN and write to PRINT file
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWPRTI')
!
!     --- compile table with overview of cpu/wall clock time used in
!         important parts of SWAN and write to PRINT file
!
      IF ( ITEST.GE.1 .OR. IDEBUG.GE.1 ) THEN
!
!        --- initialise table to zero
!
         DO K = 1, 30
            DO J = 1, 2
               TABLE(K,J) = 0D0
            END DO
         END DO
!
!        --- compute times for basic blocks
!
         DO J = 1, 2
!
!           --- total run-time
!
            TABLE(1,J) = DCUMTM(1,J)
!
!           --- initialisation, reading, preparation:
!
            DO K = 2, 7
               TABLE(2,J) = TABLE(2,J) + DCUMTM(K,J)
            END DO
!
!           --- domain decomposition:
!
            TABLE(2,J) = TABLE(2,J) + DCUMTM(211,J)
            TABLE(2,J) = TABLE(2,J) + DCUMTM(212,J)
            TABLE(2,J) = TABLE(2,J) + DCUMTM(201,J)
!
!           --- total calculation including communication:
!
            TABLE(3,J) = TABLE(3,J) + DCUMTM(8,J)
!
!           --- output:
!
            TABLE(5,J) = TABLE(5,J) + DCUMTM(9,J)
!
!           --- exchanging data:
!
            TABLE(7,J) = TABLE(7,J) + DCUMTM(213,J)
!
!           --- solving system:
!
            TABLE(9,J) = TABLE(9,J) + DCUMTM(119,J)
            TABLE(9,J) = TABLE(9,J) + DCUMTM(120,J)
!
!           --- global reductions:
!
            TABLE(10,J) = TABLE(10,J) + DCUMTM(202,J)
!
!           --- collecting data:
!
            TABLE(11,J) = TABLE(11,J) + DCUMTM(214,J)
!
!           --- setup:
!
            TABLE(12,J) = TABLE(12,J) + DCUMTM(106,J)
!
!           --- propagation velocities:
!
            TABLE(14,J) = TABLE(14,J) + DCUMTM(111,J)
            TABLE(14,J) = TABLE(14,J) + DCUMTM(113,J)
            TABLE(14,J) = TABLE(14,J) + DCUMTM(114,J)
!
!           --- x-y advection:
!
            TABLE(15,J) = TABLE(15,J) + DCUMTM(140,J)
!
!           --- sigma advection:
!
            TABLE(16,J) = TABLE(16,J) + DCUMTM(141,J)
!
!           --- theta advection:
!
            TABLE(17,J) = TABLE(17,J) + DCUMTM(142,J)
!
!           --- wind:
!
            TABLE(18,J) = TABLE(18,J) + DCUMTM(132,J)
!
!           --- whitecapping:
!
            TABLE(19,J) = TABLE(19,J) + DCUMTM(133,J)
!
!           --- bottom friction:
!
            TABLE(20,J) = TABLE(20,J) + DCUMTM(130,J)
!
!           --- wave breaking:
!
            TABLE(21,J) = TABLE(21,J) + DCUMTM(131,J)
!
!           --- quadruplets:
!
            TABLE(22,J) = TABLE(22,J) + DCUMTM(135,J)
!
!           --- triads:
!
            TABLE(23,J) = TABLE(23,J) + DCUMTM(134,J)
!
!           --- limiter:
!
            TABLE(24,J) = TABLE(24,J) + DCUMTM(122,J)
!
!           --- rescaling:
!
            TABLE(25,J) = TABLE(25,J) + DCUMTM(121,J)
!
!           --- reflections:
!
            TABLE(26,J) = TABLE(26,J) + DCUMTM(136,J)

         END DO
!
!        --- add up times for some basic blocks
!
         DO J = 1, 2
!
!           --- total calculation:
!
            TABLE(3,J) = TABLE(3,J) - TABLE( 7,J)
            TABLE(3,J) = TABLE(3,J) - TABLE(10,J)
            IF ( TABLE(3,J).LT.0D0 ) TABLE(3,J) = 0D0
!
!           --- total communication:
!                * exchanging data
!                * global reductions
!                * collecting data
!
            TABLE(4,J) = TABLE(4,J) + TABLE( 7,J)
            TABLE(4,J) = TABLE(4,J) + TABLE(10,J)
            TABLE(4,J) = TABLE(4,J) + TABLE(11,J)
!
!           --- total propagation:
!                * velocities and derivatives
!
            TABLE(6,J) = TABLE(6,J) + TABLE(14,J)
            TABLE(6,J) = TABLE(6,J) + TABLE(15,J)
            TABLE(6,J) = TABLE(6,J) + TABLE(16,J)
            TABLE(6,J) = TABLE(6,J) + TABLE(17,J)
!
!           --- sources:
!                * wind, whitecapping, friction, breaking,
!                * quadruplets, triads, limiter, rescaling,
!                * reflections
!
            DO K = 18, 26
               TABLE(8,J) = TABLE(8,J) + TABLE(K,J)
            END DO
!
!           --- other computing:
!
            TABLE(13,J) = TABLE(13,J) + TABLE( 3,J)
            TABLE(13,J) = TABLE(13,J) - TABLE( 6,J)
            TABLE(13,J) = TABLE(13,J) - TABLE( 8,J)
            TABLE(13,J) = TABLE(13,J) - TABLE( 9,J)
            TABLE(13,J) = TABLE(13,J) - TABLE(12,J)
            IF ( TABLE(13,J).LT.0D0 ) TABLE(13,J) = 0D0

         END DO
!
!        --- print CPU-times used in important parts of SWAN
!
         WRITE(PRINTF,'(/)')
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,111) INODE
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,112) INODE
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,115) INODE,'total time:'       ,(TABLE(1,J),J=1,2)
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,115) INODE,'total pre-processing:',
     &                                                (TABLE(2,j),j=1,2)
         WRITE(PRINTF,115) INODE,'total calculation:',(TABLE(3,j),j=1,2)
         WRITE(PRINTF,115) INODE,'total communication:',
     &                                                (TABLE(4,j),j=1,2)
         WRITE(PRINTF,115) INODE,'total post-processing:',
     &                                                (TABLE(5,j),j=1,2)
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,113) INODE
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,115) INODE,'calc. propagation:',(TABLE(6,j),j=1,2)
         WRITE(PRINTF,115) INODE,'exchanging data:'  ,(TABLE(7,j),j=1,2)
         WRITE(PRINTF,115) INODE,'calc. sources:'    ,(TABLE(8,j),j=1,2)
         WRITE(PRINTF,115) INODE,'solving system:'   ,(TABLE(9,j),j=1,2)
         WRITE(PRINTF,115) INODE,'reductions:'      ,(TABLE(10,j),j=1,2)
         WRITE(PRINTF,115) INODE,'collecting data:' ,(TABLE(11,j),j=1,2)
         WRITE(PRINTF,115) INODE,'calc. setup:'     ,(TABLE(12,j),j=1,2)
         WRITE(PRINTF,115) INODE,'other computing:' ,(TABLE(13,j),j=1,2)
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,114) INODE
         WRITE(PRINTF,110) INODE
         WRITE(PRINTF,115) INODE,'prop. velocities:',(TABLE(14,j),j=1,2)
         WRITE(PRINTF,115) INODE,'x-y advection:'   ,(TABLE(15,j),j=1,2)
         WRITE(PRINTF,115) INODE,'sigma advection:' ,(TABLE(16,j),j=1,2)
         WRITE(PRINTF,115) INODE,'theta advection:' ,(TABLE(17,j),j=1,2)
         WRITE(PRINTF,115) INODE,'wind:'            ,(TABLE(18,j),j=1,2)
         WRITE(PRINTF,115) INODE,'whitecapping:'    ,(TABLE(19,j),j=1,2)
         WRITE(PRINTF,115) INODE,'bottom friction:' ,(TABLE(20,j),j=1,2)
         WRITE(PRINTF,115) INODE,'wave breaking:'   ,(TABLE(21,j),j=1,2)
         WRITE(PRINTF,115) INODE,'quadruplets:'     ,(TABLE(22,j),j=1,2)
         WRITE(PRINTF,115) INODE,'triads:'          ,(TABLE(23,j),j=1,2)
         WRITE(PRINTF,115) INODE,'limiter:'         ,(TABLE(24,j),j=1,2)
         WRITE(PRINTF,115) INODE,'rescaling:'       ,(TABLE(25,j),j=1,2)
         WRITE(PRINTF,115) INODE,'reflections:'     ,(TABLE(26,j),j=1,2)

      END IF

      IF ( IDEBUG.GE.2 ) THEN
         WRITE(PRINTF,120) INODE
         DO J = 1, NSECTM
            IF (NCUMTM(J).GT.0)
     &         WRITE(PRINTF,121) INODE,J,DCUMTM(J,1),DCUMTM(J,2),
     &                           NCUMTM(J)
         END DO
      END IF

  110 FORMAT(i3,1x,'#')
  111 FORMAT(i3,' # Details on timings of the simulation:')
  112 FORMAT(i3,1x,'#',26x,'cpu-time',1x,'wall-clock')
  113 FORMAT(i3,' # Splitting up calc. + comm. times:')
  114 FORMAT(i3,' # Overview source contributions:')
  115 FORMAT(i3,1x,'#',1x,a22,2f11.2)

  120 FORMAT(/,i3,' #    item     cpu-time    real time     count')
  121 FORMAT(i3,1x,'#',4x,i4,2f13.4,i10)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE TXPBLA(TEXT,IF,IL)
!
!****************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     determines the position of the first and the last non-blank
!     (or non-tabulator) character in the text-string
!
!  4. Argument variables
!
!     IF          position of the first non-blank character in TEXT
!     IL          position of the last non-blank character in TEXT
!     TEXT        text string
!
      INTEGER IF, IL
      CHARACTER*(*) TEXT
!
!  6. Local variables
!
!     FOUND :     TEXT is found or not
!     ITABVL:     integer value of tabulator character
!     LENTXT:     length of TEXT
!
      INTEGER LENTXT, ITABVL
      LOGICAL FOUND
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
!DOS      ITABVL = ICHAR('	')
      ITABVL = ICHAR('	')
      LENTXT = LEN (TEXT)
      IF = 1
      FOUND = .FALSE.
  100 IF (IF .LE. LENTXT .AND. .NOT. FOUND) THEN
         IF (.NOT. (TEXT(IF:IF) .EQ. ' ' .OR.
     &              ICHAR(TEXT(IF:IF)) .EQ. ITABVL)) THEN
            FOUND = .TRUE.
         ELSE
            IF = IF + 1
         ENDIF
         GOTO 100
      ENDIF
      IL = LENTXT + 1
      FOUND = .FALSE.
  200 IF (IL .GT. 1 .AND. .NOT. FOUND) THEN
         IL = IL - 1
         IF (.NOT. (TEXT(IL:IL) .EQ. ' ' .OR.
     &              ICHAR(TEXT(IL:IL)) .EQ. ITABVL)) THEN
            FOUND = .TRUE.
         ENDIF
         GOTO 200
      ENDIF

      RETURN
      END
!****************************************************************
!
      CHARACTER*20 FUNCTION INTSTR ( IVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Convert integer to string
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!
      INTEGER IVAL
!
!  6. Local variables
!
!     CVAL  :     character represented an integer of mantisse
!     I     :     counter
!     IPOS  :     position in mantisse
!     IQUO  :     whole quotient
!
      INTEGER I, IPOS, IQUO
      CHARACTER*1, ALLOCATABLE :: CVAL(:)
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IPOS = 1
 100  CONTINUE
      IF (IVAL/10**IPOS.GE.1.) THEN
         IPOS = IPOS + 1
         GO TO 100
      END IF
      ALLOCATE(CVAL(IPOS))

      DO I=IPOS,1,-1
         IQUO=IVAL/10**(I-1)
         CVAL(IPOS-I+1)=CHAR(INT(IQUO)+48)
         IVAL=IVAL-IQUO*10**(I-1)
      END DO

      WRITE (INTSTR,*) (CVAL(I), I=1,IPOS)

      RETURN
      END
!****************************************************************
!
      CHARACTER*20 FUNCTION NUMSTR ( IVAL, RVAL, FORM )
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Convert integer or real to string with given format
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!     FORM        given format
!     RVAL        real to be converted
!
      INTEGER   IVAL
      REAL      RVAL
      CHARACTER FORM*20
!
!  6. Local variables
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IF ( IVAL.NE.INAN ) THEN
         WRITE (NUMSTR,FORM) IVAL
      ELSE IF ( RVAL.NE.RNAN ) THEN
         WRITE (NUMSTR,FORM) RVAL
      ELSE
         NUMSTR = ''
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOPI ( IARR1, IARR2, LENGTH )
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Copies integer array IARR1 to IARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IARR1       source array
!     IARR2       target array
!     LENGTH      array length
!
      INTEGER LENGTH
      INTEGER IARR1(LENGTH), IARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
      INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOPI')

!     --- check array length

      IF ( LENGTH.LE.0 ) THEN
         CALL MSGERR( 3, 'Array length should be positive' )
      END IF

!     --- copy elements of array IARR1 to IARR2

      DO 10 I = 1, LENGTH
         IARR2(I) = IARR1(I)
  10  CONTINUE

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOPR ( ARR1, ARR2, LENGTH )
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Copies real array ARR1 to ARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     ARR1        source array
!     ARR2        target array
!     LENGTH      array length
!
      INTEGER LENGTH
      REAL    ARR1(LENGTH), ARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
      INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOPR')

!     --- check array length

      IF ( LENGTH.LE.0 ) THEN
         CALL MSGERR( 3, 'Array length should be positive' )
      END IF

!     --- copy elements of array ARR1 to ARR2

      DO 10 I = 1, LENGTH
         ARR2(I) = ARR1(I)
  10  CONTINUE

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWI2B ( IVAL, BVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May 03: New subroutine
!
!  2. Purpose
!
!     Calculates 32-bit representation of an integer number
!
!  3. Method
!
!     The representation of an integer number is divided into 4 parts
!     of 8 bits each, resulting in 32-bit word in memory. Generally,
!     storage words are represented with bits counted from the right,
!     making bit 0 the lower-order bit and bit 31 the high-order bit,
!     which is also the sign bit.
!
!     The integer number is always an exact representation of an
!     integer of value positive, negative, or zero. Each bit, except
!     the leftmost bit, corresponds to the actual exponent as power
!     of two.
!
!     For representing negative numbers, the method called
!     "excess 2**(m - 1)" is used, which represents an m-bit number by
!     storing it as the sum of itself and 2**(m - 1). For a 32-bit
!     machine, m = 32. This results in a positive number, so the
!     leftmost bit need to be reversed. This method is identical to the
!     two's complement method.
!
!     An example:
!
!        the 32-bit representation of 5693 is
!
!        decimal    :     0        0       22       61
!        hexidecimal:     0        0       16       3D
!        binary     : 00000000 00000000 00010110 00111101
!
!        since,
!
!        5693 = 2^12 + 2^10 + 2^9 + 2^5 + 2^4 + 2^3 + 2^2 + 2^0
!
!  4. Argument variables
!
!     BVAL        a byte value as a part of the representation of
!                 integer number
!     IVAL        integer number
!
      INTEGER BVAL(4), IVAL
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!     IQUOT :     auxiliary integer with quotient
!     M     :     maximal exponent number possible (for 32-bit machine, m=32)
!
      INTEGER I, IENT, IQUOT, M
!
! 12. Structure
!
!     initialise 4 parts of the representation
!     if integer < 0, increased it by 2**(m-1)
!     compute the actual part of the representation
!     determine the sign bit
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/, M/32/
      IF (LTRACE) CALL STRACE (IENT,'SWI2B')

!     --- initialise 4 parts of the representation

      DO 10 I = 1, 4
         BVAL(I) = 0
 10   CONTINUE

      IQUOT = IVAL

!     --- if integer < 0, increased it by 2*(m-1)

      IF ( IVAL.LT.0 ) IQUOT = IQUOT + 2**(M-1)

!     --- compute the actual part of the representation

      DO 20 I = 4, 1, -1
         BVAL(I) = MOD(IQUOT,256)
         IQUOT = INT(IQUOT/256)
 20   CONTINUE

!     --- determine the sign bit

      IF ( IVAL.LT.0 ) BVAL(1) = BVAL(1) + 128

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWR2B ( RVAL, BVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May 03: New subroutine
!
!  2. Purpose
!
!     Calculates 32-bit representation of a floating-point number
!
!  3. Method
!
!     The representation of a floating-point number is divided into 4
!     parts of 8 bits each, resulting in 32-bit word in memory.
!     Generally, storage words are represented with bits counted from
!     the right, making bit 0 the lower-order bit and bit 31 the
!     high-order bit, which is also the sign bit.
!
!     The floating-point number is a processor approximation. Its format
!     has an 8-bit biased exponent and a 23-bit fraction or mantissa. The
!     leftmost bit is the sign bit which is zero for plus and 1 for
!     minus. The biased exponent equals the bias and the actual exponent
!     (power of two) of the number. For a 32-bit machine, bias=127.
!
!     Furthermore, the floating-point number is usually stored in the
!     normalized form, i.e. it has a binary point to the left of the
!     mantissa and an implied leading 1 to the left of the binary point.
!     Thus, if X is a floating-point number, then it is calculated as
!     follows:
!
!         X = (-1)**sign bit 1.fraction * 2**(biased exponent-bias)
!
!     There are several exceptions. Let a fraction, biased exponent
!     and sign bit be denoted as F, E and S, respectively. The following
!     formats adhere to IEEE standard:
!
!     S = 0, E = 00000000 and F  = 00 ... 0 : X = 0
!     S = 0, E = 00000000 and F <> 00 ... 0 : X = +0.fraction * 2**(1-bias)
!     S = 1, E = 00000000 and F <> 00 ... 0 : X = -0.fraction * 2**(1-bias)
!     S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
!     S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
!     S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
!
!     A NaN (Not a Number) is a value reserved for signalling an
!     attempted invalid operation, like 0/0. Its representation
!     equals the representation of +Inf plus 1, i.e. 2**31 - 2**23 + 1
!
!     An example:
!
!        the 32-bit representation of 23.1 is
!
!        decimal    :    65      184      204      205
!        hexidecimal:    41       B8       CC       CD
!        binary     : 01000001 10111000 11001100 11001101
!
!        since,
!
!        23.1 = 2^4 + 2^2 + 2^1 + 2^0 + 2^-4 + 2^-5 + 2^-8 + 2^-9 +
!               2^-12 + 2^-13 + 2^-16 + 2^-17 + 2^-19
!
!        so that the biased exponent = 4 + 127 = 131 = 10000011 = E
!        and the sign bit = 0 = S. The remaining of the 32-bit word is
!        the fraction, which is
!
!    3 2 1 0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19
!
! F= 0 1 1 1  0  0  0  1  1  0  0  1  1   0   0   1   1   0   0   1   1   0   1
!
!  4. Argument variables
!
!     BVAL        a byte value as a part of the representation of
!                 floating-point number
!     RVAL        floating-point number
!
      INTEGER BVAL(4)
      REAL    RVAL
!
!  6. Local variables
!
!     ACTEXP:     actual exponent in the representation
!     BEXPO :     biased exponent
!     BIAS  :     bias (for 32-bit machine, bias=127)
!     EXPO  :     calculated exponent of floating-point number
!     FRAC  :     fraction of floating-point number
!     I     :     loop counter
!     IENT  :     number of entries
!     IPART :     i-the part of the representation
!     IQUOT :     auxiliary integer with quotient
!     LEADNR:     leading number of floating-point number
!     LFRAC :     length of fraction in representation
!           :     (for 32-bit machine, lfrac=23)
!     RFRAC :     auxiliary real with fraction
!
      INTEGER ACTEXP, EXPO, I, IENT, IPART, LFRAC, BIAS, BEXPO,
     &        LEADNR, IQUOT
      REAL FRAC, RFRAC
!
! 12. Structure
!
!     initialise 4 parts of the representation and biased exponent
!     determine leading number and fraction
!     do while leading number >= 1
!        calculate positive exponent as power of two
!     or while fraction > 0
!        calculate negative exponent as power of two
!     end do
!     compute the actual part of the representation
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/, LFRAC/23/, BIAS/127/
      IF (LTRACE) CALL STRACE (IENT,'SWR2B')

!     --- initialise 4 parts of the representation and biased exponent

      DO 10 I = 1, 4
         BVAL(I) = 0
 10   CONTINUE
      BEXPO  = -1

!     --- determine leading number and fraction

      IF ( ABS(RVAL).LT.1.E-7 ) THEN
         LEADNR = 0
         FRAC   = 0.
      ELSE
         LEADNR = INT(ABS(RVAL))
         FRAC   = ABS(RVAL) - REAL(LEADNR)
      END IF

 20   IF ( LEADNR.GE.1 ) THEN

!        --- calculate positive exponent as power of two

         EXPO  = 0
         IQUOT = LEADNR
 30      IF ( IQUOT.GE.2 ) THEN

              IQUOT = INT(IQUOT/2)
              EXPO  = EXPO + 1

         GOTO 30
         END IF

      ELSE IF ( FRAC.GT.0. ) THEN

!        --- calculate negative exponent as power of two

         EXPO = 0
         RFRAC = FRAC
 40      IF ( RFRAC.LT.1. ) THEN

              RFRAC = RFRAC * 2.
              EXPO  = EXPO - 1

         GOTO 40
         END IF

      ELSE

         GOTO 50

      END IF

!     --- compute the actual part of the representation

      IF ( BEXPO.EQ.-1 ) THEN

!        --- determine biased exponent

         BEXPO = EXPO + BIAS

!        --- the first seven bits of biased exponent belong
!            to first part of the representation

         BVAL(1) = INT(BEXPO/2)

!        --- determine the sign bit

         IF ( RVAL.LT.0. ) BVAL(1) = BVAL(1) + 128

!        --- the eighth bit of biased component is the leftmost
!            bit of second part of the representation

         BVAL(2) = MOD(BEXPO,2)*2**7
         IPART = 2

      ELSE

!        --- compute the actual exponent of bit 1 in i-th part of
!            the representation

         ACTEXP = (IPART-2)*8 + 7 - BEXPO + BIAS + EXPO
         IF ( ACTEXP.LT.0 ) THEN
            ACTEXP = ACTEXP + 8
            IPART = IPART + 1
            IF ( IPART.GT.4 ) GOTO 50
         END IF
         BVAL(IPART) = BVAL(IPART) + 2**ACTEXP

      END IF

      IF ( EXPO.LT.(BEXPO-BIAS-LFRAC) ) GOTO 50
      LEADNR = LEADNR - 2.**EXPO
      IF ( EXPO.LT.0 ) FRAC = FRAC - 2.**EXPO

      GOTO 20
 50   CONTINUE

      RETURN
      END
