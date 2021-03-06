!
!     SWAN/COMPU   file 4 of 6
!
!
!     PROGRAM SWANCOM4.FOR
!
!
!     This file SWANCOM4 of the main program SWAN
!     include the next subroutines
!
!     *** nonlinear 4 wave-wave interactions ***
!
!     BND4WW (determine boundary for arrays for 4 wave-wave interactions
!             to allocate some memory in WAREA)
!     FAC4WW (compute the constants for the nonlinear wave
!             interactions)
!     RANGE4 (compute the counters for the different types of
!             computations for the nonlinear wave interactions)
!     SWSNL1 (nonlinear four wave interactions; semi-implicit and computed
!             for all bins that fall within a sweep with DIA technique.
!             Interaction are calculated per sweep)
!     SWSNL2 (nonlinear four wave interactions; fully explicit and computed
!             for all bins that fall within a sweep with DIA technique.
!             Interaction are calculated per sweep)
!     SWSNL3 (calculate nonlinear four wave interactions fully explicitly
!             for the full circle per iteration by means of DIA approach
!             and store results in auxiliary array MEMNL4)
!     SWSNL4 (calculate nonlinear four wave interactions fully explicitly
!             for the full circle per iteration by means of MDIA approach
!             and store results in auxiliary array MEMNL4)
!     FILNL3 (fill main diagonal and right-hand side of the system with
!             results of array MEMNL4)
!
!     *** nonlinear 3 wave-wave interactions ***
!
!     STRIAD (nonlinear three wave interactions using a semi-
!             implicit scheme)
!     STRIAN (triad-wave interactions calculated with the Lumped Triad
!             Approximation of Eldeberky, 1996)
!     RIAM_SLW (calculate nonlinear four wave interactions by means
!               of the exact FD-RIAM technique)
!
!----------------------------------------------------------------------
!
!******************************************************************
!
      SUBROUTINE BND4WW (MSCMAX,MDCMAX,SPCSIG              )              34.00
!
!******************************************************************
!
      INCLUDE 'swcomm3.inc'                                               34.00
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!
!  2. Purpose
!
!     compute the array size for the nonlinear 4 wave
!     interactions in order to allocate some memory in
!     the WAREA
!
!  3. Method
!
!     For the method see comment at subroutine SWSNL1
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     5. SUBROUTINES CALLING
!
!        ----
!
!     6. SUBROUTINES USED
!
!        none
!
!     7. ERROR MESSAGES
!
!        none
!
!     8. REMARKS
!
!        ---
!
!     9. STRUCTURE
!
!        ---
!
!     10. SOURCE TEXT
!
!***********************************************************************
!
      INTEGER      IDP   ,IDP1  ,IDM   ,IDM1  ,MSC2  ,                    34.00
     &             MSC1  ,ISP   ,ISP1  ,ISM   ,ISM1  ,ISLOW ,ISHGH ,
     &             MSCMAX,MDCMAX,                                         34.00
     &             IDLOW ,IDHGH
!
      REAL         LAMBDA,LAMM2 ,LAMP2 ,DELTH3,AUX1  ,DELTH4,CIDP  ,
     &             CIDM  ,XIS   ,XISLN ,RATE                              34.00
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'BND4WW')                             34.00
!
!     *** compute the auxiliary array boundaries for the 4 wave ***
!     *** interactions                                          ***
!
      LAMBDA = PQUAD(1)                                                   34.00
      LAMM2  = (1.-LAMBDA)**2
      LAMP2  = (1.+LAMBDA)**2
      DELTH3 = ACOS( (LAMM2**2+4.-LAMP2**2) / (4.*LAMM2) )
      AUX1   = SIN(DELTH3)
      DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
      CIDP  = ABS(DELTH4/DDIR)                                            34.00
      IDP   = INT(CIDP)
      IDP1  = IDP + 1
!
      CIDM  = ABS(DELTH3/DDIR)                                            34.00
      IDM   = INT(CIDM)
      IDM1  = IDM + 1
!
      MSC2   = INT ( FLOAT(MSC) / 2.0 )
      MSC1   = MSC2 - 1
!
      XIS    = SPCSIG(MSC2) / SPCSIG(MSC1)
      XISLN  = LOG( XIS )
      ISP    = INT( LOG(1.+LAMBDA) / XISLN )
      ISP1   = ISP + 1
      ISM    = INT( LOG(1.-LAMBDA) / XISLN )
      ISM1   = ISM - 1
!
!     *** Range of array size and calculations ***
!
      ISLOW =  1  + ISM1
      ISHGH = MSC + ISP1 - ISM1
!C      IDLOW = 1 - MAX(IDM1,IDP1)
!C      IDHGH = MDC + MAX(IDM1,IDP1)
      IDLOW = 1   - MDC - MAX (IDM1, IDP1)
      IDHGH = MDC + MDC + MAX (IDM1, IDP1)
      MSC4MI = ISLOW
      MSC4MA = ISHGH
      MDC4MI = IDLOW
      MDC4MA = IDHGH
      MSCMAX = MSC4MA - MSC4MI + 1
      MDCMAX = MDC4MA - MDC4MI + 1
!
!     *** Test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 50 ) THEN
        WRITE(PRINTF,*) ' subroutine BND4WW:'
        RATE = 180. / PI
        WRITE(PRINTF,*) ' DELTH3 DELTH4 DDIR     :',DELTH3*RATE,          34.00
     &  DELTH4*RATE,DDIR*RATE                                             34.00
        WRITE(PRINTF,*) ' ISM ISM1 ISP ISP1      :',ISM,ISM1,ISP,ISP1
        WRITE(PRINTF,*) ' IDM IDM1 IDP IDP1      :',IDM,IDM1,IDP,IDP1
        WRITE(PRINTF,*) ' MDC MSC  XIS XISLN     :',MDC,MSC,XIS,XISLN
        WRITE(PRINTF,*) ' ISLOW,ISHGH,IDLOW,IDHGH:',ISLOW,ISHGH,IDLOW
     &                 ,IDHGH
        WRITE(PRINTF,*) ' S4MI S4MA D4MI D4MA    :',MSC4MI,MSC4MA,
     &                    MDC4MI,MDC4MA
        WRITE(PRINTF,*) ' MSCMAX MDCMAX          :',MSCMAX, MDCMAX
      END IF
!
!     End of the subroutine BND4WW
      RETURN
      END
!
!******************************************************************
!
      SUBROUTINE FAC4WW (ITER ,XIS   ,SNLC1 ,                             34.00
     &                  DAL1  ,DAL2  ,DAL3         ,SPCSIG,               34.00
     &                  WWINT ,WWAWG ,WWSWG                )              40.17 34.00
!
!******************************************************************

      USE M_SNL4                                                          40.17
!
      INCLUDE 'swcomm3.inc'                                               34.00
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.17: IJsbrand Haagsma
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implementation of Multiple DIA
!
!  2. Purpose :
!
!     Calculate interpolation constants for Snl.
!
!  3. Method :
!
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     INTEGERS:
!     ---------
!     ITER              Iteration number
!     MSC2,MSC1         Auxiliary variables
!     MSC,MDC           Maximum counters in spectral space
!     IDP,IDP1          Positive range for ID
!     IDM,IDM1          Negative range for ID
!     ISP,ISP1          idem for IS
!     ISM,ISM1          idem for IS
!     ISCLW,ISCHG       Minimum and maximum counter for discrete
!                       computations in frequency space
!     ISLOW,ISHGH       Minimum and maximum range in frequency space
!     IDLOW,IDHGH       idem in directional space
!     IS                Frequency counter
!     MSC4MI,MSC4MA     Array dimensions in frequency space
!     MDC4MI,MDC4MA     Array dimensions in direction space
!
!     REALS:
!     ------
!     LAMBDA            Coefficient set 0.25
!     GRAV              Gravitational acceleration
!     SNLC1             Coefficient for the subroutines SWSNLn
!     LAMM2,LAMP2
!     DELTH3,DELTH4     Angles between the interacting wavenumbers
!     DAL1,DAL2,DAL3    Coefficients for the non linear interactions
!     CIDP,CIDM
!     WIDP,WIDP1,WIDM,WIDM1  Weight factors
!     WISP,WISP1,WISM,WISM1  idem
!     AWGn              Interpolation weight factors
!     SWGn              Quadratic interpolation weight factors
!     XIS,XISLN         Difference between succeeding frequencies
!     PI                3.14
!     FREQ              Auxiliary frequency to fill scaling array
!     DDIR,RADE         band width in directional space and factor        34.00
!
!     ARRAYS
!     ------
!     AF11    1D   Scaling frequency
!     WWINT   1D   counters for 4WAVE interactions
!     WWAWG   1D   values for the interpolation
!     WWSWG   1D   vaules for the interpolation
!
!     WWINT ( 1 = IDP    WWAWG ( = AGW1    WWSWG ( = SWG1
!             2 = IDP1           = AWG2            = SWG2
!             3 = IDM            = AWG3            = SWG3
!             4 = IDM1           = AWG4            = SWG4
!             5 = ISP            = AWG5            = SWG5
!             6 = ISP1           = AWG6            = SWG6
!             7 = ISM            = AWG7            = SWG7
!             8 = ISM1           = AWG8 )          = SWG8  )
!             9 = ISLOW
!             10= ISHGH
!             11= ISCLW
!             12= ISCHG
!             14= IDLOW
!             15= IDHGH
!             16= MSC4MI
!             17= MSC4MA
!             18= MDC4MI
!             19= MDC4MA
!             20= MSCMAX
!             21= MDCMAX )
!
!  4. Subroutines used :
!
!     ---
!
!  5. Called by :
!
!     ---
!
!  6. Error messages :
!
!     ---
!
!  9. Source code :
!
!     -----------------------------------------------------------------
!     Calculate :
!       1. counters for frequency and direction for NL-interaction
!       2. weight factors
!       3. the minimum and maximum counter in IS and ID space
!       4. the interpolation weights
!       5. the quadratic interpolation rates
!       6. fill the array for the frequency**11
!     ----------------------------------------------------------
!
!****************************************************************
!
      INTEGER     ITER  ,MSC2  ,MSC1  ,IS    ,IDP   ,IDP1  ,              34.00
     &            IDM   ,IDM1  ,ISP   ,ISP1  ,ISM   ,ISM1  ,              34.00
     &            ISLOW ,ISHGH ,ISCLW ,ISCHG ,IDLOW ,IDHGH ,
     &            MSCMAX,MDCMAX                                           34.00
!
      REAL        SNLC1 ,LAMM2 ,LAMP2 ,DELTH3,                            40.17 34.00
     &            AUX1  ,DELTH4,DAL1  ,DAL2  ,DAL3  ,CIDP  ,WIDP  ,
     &            WIDP1 ,CIDM  ,WIDM  ,WIDM1 ,XIS   ,XISLN ,WISP  ,
     &            WISP1 ,WISM  ,WISM1 ,AWG1  ,AWG2  ,AWG3  ,AWG4  ,
     &            AWG5  ,AWG6  ,AWG7  ,AWG8  ,SWG1  ,SWG2  ,SWG3  ,
     &            SWG4  ,SWG5  ,SWG6  ,SWG7  ,SWG8  ,FREQ  ,              34.00
     &            RADE                                                    34.00
!
      REAL       WWAWG(*)               ,                                 40.17
     &           WWSWG(*)
!
      INTEGER    WWINT(*)
!
      SAVE IENT
      DATA IENT/0/

      IF (LTRACE) CALL STRACE (IENT,'FAC4WW')

      IF (ALLOCATED(AF11)) DEALLOCATE(AF11)                               40.17

!     *** Compute frequency indices                               ***
!     *** XIS is the relative increment of the relative frequency ***
!
      MSC2   = INT ( FLOAT(MSC) / 2.0 )
      MSC1   = MSC2 - 1
      XIS    = SPCSIG(MSC2) / SPCSIG(MSC1)                                30.72
!
!     *** set values for the nonlinear four-wave interactions ***
!
      SNLC1  = 1. / GRAV**4                                               40.17 34.00
!
      LAMM2  = (1.-PQUAD(1))**2                                           40.17
      LAMP2  = (1.+PQUAD(1))**2                                           40.17
      DELTH3 = ACOS( (LAMM2**2+4.-LAMP2**2) / (4.*LAMM2) )
      AUX1   = SIN(DELTH3)
      DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
      DAL1   = 1. / (1.+PQUAD(1))**4                                      40.17
      DAL2   = 1. / (1.-PQUAD(1))**4                                      40.17
      DAL3   = 2. * DAL1 * DAL2
!
!     *** Compute directional indices in sigma and theta space ***
!
      CIDP   = ABS(DELTH4/DDIR)                                           40.00
      IDP   = INT(CIDP)
      IDP1  = IDP + 1
      WIDP   = CIDP - REAL(IDP)
      WIDP1  = 1.- WIDP
!
      CIDM   = ABS(DELTH3/DDIR)                                           40.00
      IDM   = INT(CIDM)
      IDM1  = IDM + 1
      WIDM   = CIDM - REAL(IDM)
      WIDM1  = 1.- WIDM
      XISLN  = LOG( XIS )
!
      ISP    = INT( LOG(1.+PQUAD(1)) / XISLN )                            40.17
      ISP1   = ISP + 1
      WISP   = (1.+PQUAD(1) - XIS**ISP) / (XIS**ISP1 - XIS**ISP)          40.17
      WISP1  = 1. - WISP
!
      ISM    = INT( LOG(1.-PQUAD(1)) / XISLN )                            40.17
      ISM1   = ISM - 1
      WISM   = (XIS**ISM -(1.-PQUAD(1))) / (XIS**ISM - XIS**ISM1)         40.17
      WISM1  = 1. - WISM
!
!     *** Range of calculations ***
!
      ISLOW =  1  + ISM1
      ISHGH = MSC + ISP1 - ISM1
      ISCLW =  1
      ISCHG = MSC - ISM1
!C      IDLOW = 1 - MAX(IDM1,IDP1)
!C      IDHGH = MDC + MAX(IDM1,IDP1)
      IDLOW = 1 - MDC - MAX(IDM1,IDP1)
      IDHGH = MDC + MDC + MAX(IDM1,IDP1)
!
      MSC4MI = ISLOW
      MSC4MA = ISHGH
      MDC4MI = IDLOW
      MDC4MA = IDHGH
      MSCMAX = MSC4MA - MSC4MI + 1
      MDCMAX = MDC4MA - MDC4MI + 1
!
!     *** Interpolation weights ***
!
      AWG1   = WIDP  * WISP
      AWG2   = WIDP1 * WISP
      AWG3   = WIDP  * WISP1
      AWG4   = WIDP1 * WISP1
!
      AWG5   = WIDM  * WISM
      AWG6   = WIDM1 * WISM
      AWG7   = WIDM  * WISM1
      AWG8   = WIDM1 * WISM1
!
!     *** quadratic interpolation ***
!
      SWG1   = AWG1**2
      SWG2   = AWG2**2
      SWG3   = AWG3**2
      SWG4   = AWG4**2
!
      SWG5   = AWG5**2
      SWG6   = AWG6**2
      SWG7   = AWG7**2
      SWG8   = AWG8**2
!
!     *** fill the arrays *
!
      WWINT(1) = IDP
      WWINT(2) = IDP1
      WWINT(3) = IDM
      WWINT(4) = IDM1
      WWINT(5) = ISP
      WWINT(6) = ISP1
      WWINT(7) = ISM
      WWINT(8) = ISM1
      WWINT(9) = ISLOW
      WWINT(10)= ISHGH
      WWINT(11)= ISCLW
      WWINT(12)= ISCHG
      WWINT(13)= IDLOW
      WWINT(14)= IDHGH
      WWINT(15)= MSC4MI
      WWINT(16)= MSC4MA
      WWINT(17)= MDC4MI
      WWINT(18)= MDC4MA
      WWINT(19)= MSCMAX
      WWINT(20)= MDCMAX
!
      WWAWG(1) = AWG1
      WWAWG(2) = AWG2
      WWAWG(3) = AWG3
      WWAWG(4) = AWG4
      WWAWG(5) = AWG5
      WWAWG(6) = AWG6
      WWAWG(7) = AWG7
      WWAWG(8) = AWG8
!
      WWSWG(1) = SWG1
      WWSWG(2) = SWG2
      WWSWG(3) = SWG3
      WWSWG(4) = SWG4
      WWSWG(5) = SWG5
      WWSWG(6) = SWG6
      WWSWG(7) = SWG7
      WWSWG(8) = SWG8

      ALLOCATE (AF11(MSC4MI:MSC4MA))                                      40.17

!     *** Fill scaling array (f**11)                     ***
!     *** compute the radian frequency**11 for IS=1, MSC ***
!
      DO 100 IS=1, MSC
        AF11(IS) = ( SPCSIG(IS) / ( 2. * PI ) )**11                       30.72
 100  CONTINUE
!
!     *** compute the radian frequency for the IS = MSC+1, ISHGH ***
!
      FREQ   = SPCSIG(MSC) / ( 2. * PI )                                  30.72
      DO 110 IS = MSC+1, ISHGH
        FREQ   = FREQ * XIS
        AF11(IS) = FREQ**11
 110  CONTINUE
!
!     *** compute the radian frequency for IS = 0, ISLOW ***
!
      FREQ   = SPCSIG(1) / ( 2. * PI )                                    30.72
      DO 120 IS = 0, ISLOW, -1
        FREQ   = FREQ / XIS
        AF11(IS) = FREQ**11
 120  CONTINUE
!
!     *** test output ***
!
      IF (ISLOW .LT. MSC4MI .OR. ISHGH .GT. MSC4MA .OR.
     &    IDLOW .LT. MDC4MI .OR. IDHGH .GT. MDC4MA) THEN
        WRITE (PRINTF,900) IXCGRD(1), IYCGRD(1),
     &                     ISLOW, ISHGH, IDLOW, IDHGH,
     &                     MSC4MI,MSC4MA, MDC4MI, MDC4MA
 900    FORMAT ( ' ** Error : array bounds and maxima in subr FAC4WW, ',
     &           ' point ', 2I5,
     &         /,'            ISL,ISH : ',2I4, '   IDL,IDH : ',2I4,
     &         /,'            SMI,SMA : ',2I4, '   DMI,DMA : ',2I4)
      ENDIF
!
      IF (ITEST .GE. 40) THEN
        RADE = 360.0 / ( 2. * PI )
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' FAC4WW subroutine '
        WRITE(PRINTF,9000) DELTH4*RADE, DELTH3*RADE, DDIR*RADE, XIS
 9000   FORMAT (' THET3 THET4 DDIR XIS  :',4E12.4)
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE(PRINTF,9012) WIDP, WIDP1, WIDM, WIDM1
 9012   FORMAT (' WIDP WIDP1 WIDM WIDM1 :',4E12.4)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9014) WISP, WISP1, WISM, WISM1
 9014   FORMAT (' WISP WISP1 WISM WISM1 :',4E12.4)
        WRITE(PRINTF,9016) ITER, ISCLW, ISCHG
 9016   FORMAT (' ITER ICLW ICHG        :',3I5)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW,IDHGH
 9015   FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of FAC4WW
      END
!
!******************************************************************
!
      SUBROUTINE RANGE4 (WWINT ,IDDLOW,IDDTOP)                            40.00
!
!******************************************************************
!
      INCLUDE 'swcomm3.inc'                                               40.00
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
!  0. Authors
!
!     40.00: Nico Booij
!     40.10: IJsbrand Haagsma
!
!  1. Updates
!
!     40.10, Mar 00: Made modification for exact quadruplets
!
!  2. Purpose :
!
!     calculate the minimum and maximum counters in frequency and
!     directional space which fall with the calculation for the
!     nonlinear wave-wave interactions.
!
!  3. Method :  review for the counters :
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |    ^    |
!             ISLOW   1                     MSC  |  ISHGH
!                     ^                          |
!                     |                          |
!                    ISCLW                     ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!
!       The directional counters depend on the numerical method that
!       is used.
!
!  4. Parameters :
!
!     INTEGER
!     -------
!     IQUAD         Counter for 4 wave interactions
!     ISLOW,ISHGH   Minimum and maximum counter in frequency space
!     ISCLW,ISCHG   idem for discrete computations
!     IDLOW,IDHGH   Minimum and maximum counters in directional space
!     MSC,MDC       Range of the original arrays
!     ISM1,ISP1,
!     IDM1,IDP1     see subroutine FAC4WW
!     IDDLOW        minimum counter of the bin that is propagated
!                   within a sweep
!     IDDTOP        minimum counter of the bin that is propagated
!                   within a sweep
!
!     array:
!     ------
!     WWINT         counters for the nonlinear interactions
!
!     WWINT ( 1  = IDP      2  = IDP1     3  = IDM     4  = IDM1
!             5  = ISP      6  = ISP1     7  = ISM     8  = ISM1
!             9  = ISLOW    10 = ISHGH    11 = ISCLW   12 = ISCHG
!             13 = IDLOW    14 = IDHGH    15 = MSC4MI  16 = MSC4MA
!             17 = MDC4MI   18 = MDC4MA
!             19 = MSCMAX   20 = MDCMAX )
!
!  5. Subroutines used :
!
!     ---
!
!  6. Called by :
!
!     SOURCE
!
!  7. Error messages :
!
!     ---
!
!  9. Source code :
!
!     -----------------------------------------------------------------
!     Calculate :
!       In absence of a current there are always four sectors
!         equal 90 degrees within a sweep that are propagated
!         Extend the boundaries to calculate the source term
!       In presence of a current and if IDTOT .eq. MDC then calculate
!         boundaries for calculation of interaction using the
!         unfolded area.
!     ----------------------------------------------------------
!
!****************************************************************
!
      INTEGER     IDDLOW,IDDTOP                                           40.00
!
      INTEGER     WWINT(*)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'RANGE4')
!
!     *** Range in frequency domain ***
!
      WWINT(9)  =  1  + WWINT(8)
      WWINT(10) = MSC + WWINT(6) - WWINT(8)
      WWINT(11) =  1
      WWINT(12) = MSC - WWINT(8)
!
!     *** Range in directional domain ***
!
      IF ( IQUAD .LT. 3 .AND. IQUAD .GT. 0 ) THEN                         40.10
!       *** counters based on bins which fall within a sweep ***
        WWINT(13) = IDDLOW - MAX( WWINT(4), WWINT(2) )
        WWINT(14) = IDDTOP + MAX( WWINT(4), WWINT(2) )
      ELSE
!       *** counters initially based on full circle ***
        WWINT(13) = 1   - MAX( WWINT(4), WWINT(2) )
        WWINT(14) = MDC + MAX( WWINT(4), WWINT(2) )
      END IF
!
!     *** error message ***
!
      IF (WWINT(9)  .LT. WWINT(15) .OR. WWINT(10) .GT. WWINT(16) .OR.
     &    WWINT(13) .LT. WWINT(17) .OR. WWINT(14) .GT. WWINT(18) ) THEN
        WRITE (PRINTF,900) IXCGRD(1), IYCGRD(1),
     &                     WWINT(9) ,WWINT(10) ,WWINT(13) ,WWINT(14),
     &                     WWINT(15),WWINT(16) ,WWINT(17) ,WWINT(18)
 900    FORMAT ( ' ** Error : array bounds and maxima in subr RANGE4, ',
     &           ' point ', 2I5,
     &         /,'            ISL,ISH : ',2I4, '   IDL,IDH : ',2I4,
     &         /,'            SMI,SMA : ',2I4, '   DMI,DMA : ',2I4)
        IF (ITEST.GE.50) WRITE (PRTEST, 901) MSC, MDC, IDDLOW, IDDTOP
 901    FORMAT (' MSC, MDC, IDDLOW, IDDTOP: ', 4I5)
      ENDIF
!
!     test output
!
      IF (TESTFL .AND. ITEST .GE. 60) THEN
        WRITE(PRTEST,911) WWINT(4), WWINT(2), WWINT(8), WWINT(6)
 911    FORMAT (' RANGE4: IDM1 IDP1 ISM1 ISP1    :',4I5)
        WRITE(PRTEST,916) WWINT(11), WWINT(12), IQUAD
 916    FORMAT (' RANGE4: ISCLW ISCHG IQUAD      :',3I5)
        WRITE (PRTEST,917) WWINT(9), WWINT(10), WWINT(13), WWINT(14)
 917    FORMAT (' RANGE4: ISLOW ISHGH IDLOW IDHGH:',4I5)
        WRITE (PRTEST,919) WWINT(15), WWINT(16), WWINT(17), WWINT(18)
 919    FORMAT (' RANGE4: MS4MI MS4MA MD4MI MD4MA:',4I5)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of RANGE4
      END
!
!********************************************************************
!
      SUBROUTINE SWSNL1 (WWINT   ,WWAWG   ,WWSWG   ,                      34.00
     &                   IDCMIN  ,IDCMAX  ,UE      ,SA1     ,             40.17
     &                   SA2     ,DA1C    ,DA1P    ,DA1M    ,DA2C    ,
     &                   DA2P    ,DA2M    ,SPCSIG  ,SNLC1   ,KMESPC  ,    30.72
     &                   FACHFR  ,ISSTOP  ,DAL1    ,DAL2    ,DAL3    ,
     &                   SFNL    ,DSNL    ,DEP2    ,AC2     ,IMATDA  ,
     &                   IMATRA  ,PLNL4S  ,PLNL4D  ,                      34.00
     &                   IDDLOW  ,IDDTOP  )                               34.00
!
!********************************************************************

      USE M_SNL4                                                          40.17
!
      INCLUDE 'swcomm3.inc'                                               34.00
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.13: Nico Booij
!     40.17: IJsbrand Haagsma
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implentation of Multiple DIA
!     40.23, Aug. 02: some corrections
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988),
!     including the diagonal term for the implicit integration.
!
!     The interactions are calculated for all bin's that fall
!     within a sweep. No additional auxiliary array is required (see
!     SWSNL3)
!
!  3. Method
!
!     Discrete interaction approximation.
!
!     Since the domain in directional domain is by definition not
!     periodic, the spectral space can not beforehand
!     folded to the side angles. This can only be done if the
!     full circle has to be calculated
!
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |    ^    |
!             ISLOW   1                     MSC  |    ISHGH
!                     ^                          |
!                     |                          |
!                    ISCLW                     ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  IDDLOW - MAX(IDM1,IDP1)
!     IDHGH =  IDDTOP + MAX(IDM1,IDP1)
!
!     For the meaning of the counters on the right hand side of the
!     above equations see section 4.
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     Data in PARAMETER statements :
!     ----------------------------------------------------------------
!       DAL1    Real  LAMBDA dependend weight factors (see FAC4WW)
!       DAL2    Real
!       DAL3    Real
!       ITHP, ITHP1, ITHM, ITHM1, IFRP, IFRP1, IFRM, IFRM1
!               Int.  Counters of interpolation point relative to
!                     central bin, see figure below (set in FAC4WW).
!       NFRLOW, NFRHGH, NFRCHG, NTHLOW, NTHHGH
!               Int.  Range of calculations, see section 2.
!       AF11    R.A.  Scaling array (Freq**11).
!       AWGn    Real  Interpolation weights, see numbers in fig.
!       SWGn    Real  Id. squared.
!       UE      R.A.  "Unfolded" spectrum.
!       SA1     R.A.  Interaction constribution of first and second
!       SA2     R.A.    quadr. respectively (unfolded space).
!       DA1C, DA1P, DA1M, DA2C, DA2P, DA2M
!               R.A.  Idem for diagonal matrix.
!       PERCIR        full circle or sector
!     ----------------------------------------------------------------
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*************************************************************
!
      INTEGER   IS     ,ID     ,I      ,J      ,                          34.00
     &          ISHGH  ,IDLOW  ,ISP    ,ISP1   ,IDP    ,IDP1   ,
     &          ISM    ,ISM1   ,IDHGH  ,IDM    ,IDM1   ,ISCLW  ,
     &          ISCHG  ,IDDLOW ,IDDTOP                                    34.00
!
      REAL      X      ,X2     ,CONS   ,FACTOR ,SNLCS1 ,SNLCS2 ,SNLCS3,
     &          E00    ,EP1    ,EM1    ,EP2    ,EM2    ,SA1A   ,SA1B  ,
     &          SA2A   ,SA2B   ,KMESPC ,FACHFR ,AWG1   ,AWG2   ,AWG3  ,
     &          AWG4   ,AWG5   ,AWG6   ,AWG7   ,AWG8   ,DAL1   ,DAL2  ,
     &          DAL3   ,SNLC1  ,SWG1   ,SWG2   ,SWG3   ,SWG4   ,SWG5  ,
     &          SWG6   ,SWG7   ,SWG8           ,JACOBI ,SIGPI             34.00
!
      REAL      AC2(MDC,MSC,MCGRD)                    ,
     &          DEP2(MCGRD)                           ,
     &          UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,
     &          SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          IMATDA(MDC,MSC)                       ,
     &          IMATRA(MDC,MSC)                       ,
     &          PLNL4S(MDC,MSC,NPTST)                 ,                   40.00
     &          PLNL4D(MDC,MSC,NPTST)                 ,
     &          WWAWG(*)                              ,
     &          WWSWG(*)
!
      INTEGER   IDCMIN(MSC)        ,
     &          IDCMAX(MSC)        ,
     &          WWINT(*)
!
      LOGICAL   PERCIR
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL1')
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
      SWG1 = WWSWG(1)
      SWG2 = WWSWG(2)
      SWG3 = WWSWG(3)
      SWG4 = WWSWG(4)
      SWG5 = WWSWG(5)
      SWG6 = WWSWG(6)
      SWG7 = WWSWG(7)
      SWG8 = WWSWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
          DA1C(IS,ID) = 0.
          DA1P(IS,ID) = 0.
          DA1M(IS,ID) = 0.
          DA2C(IS,ID) = 0.
          DA2P(IS,ID) = 0.
          DA2M(IS,ID) = 0.
          DSNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***         40.17
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** check whether the spectral domain is periodic in ***
!     *** directional space and if so, modify boundaries   ***
!
      PERCIR = .FALSE.
      IF ( IDDLOW .EQ. 1 .AND. IDDTOP .EQ. MDC ) THEN
!       *** periodic in theta -> spectrum can be folded    ***
!       *** (can only be present in presence of a current) ***
        IDCLOW = 1
        IDCHGH = MDC
        IIID   = 0
        PERCIR = .TRUE.
      ELSE
!       *** different sectors per sweep -> extend range with IIID ***
        IIID   = MAX ( IDM1 , IDP1 )
        IDCLOW = IDLOW
        IDCHGH = IDHGH
      ENDIF
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
      DO IDDUM = IDLOW - IIID, IDHGH + IIID
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS = 1, MSC
          UE(IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI        30.72
        ENDDO
      ENDDO
!
!     *** set values in area 2 for IS > MSC+1  ***
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW - IIID , IDHGH + IIID
          UE (IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
      DO IS = ISCLW, ISCHG
        DO ID = IDCLOW, IDCHGH
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
!
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         *** Contribution to interactions                          ***
!         *** CONS is the shallow water factor for the NL interact. ***
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               40.17
          SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         40.17
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               40.17
          SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         40.17
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR
 9005       FORMAT (' FACTOR               : ',E12.4)
          END IF
!
          DA1C(IS,ID) = CONS * AF11(IS) * ( SA1A + SA1B )
          DA1P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM1 ) * PQUAD(2)       40.23
          DA1M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP1 ) * PQUAD(2)       40.23
!
          DA2C(IS,ID) = CONS * AF11(IS) * ( SA2A + SA2B )
          DA2P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM2 ) * PQUAD(2)       40.23
          DA2M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP2 ) * PQUAD(2)       40.23
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
      IF ( PERCIR ) THEN
        DO ID = 1, IDHGH - MDC
          ID0   = 1 - ID
          DO IS = ISCLW, ISCHG
            SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
            SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
            DA1C(IS,MDC+ID) = DA1C(IS,  ID   )
            DA1P(IS,MDC+ID) = DA1P(IS,  ID   )
            DA1M(IS,MDC+ID) = DA1M(IS,  ID   )
            DA2C(IS,MDC+ID) = DA2C(IS,  ID   )
            DA2P(IS,MDC+ID) = DA2P(IS,  ID   )
            DA2M(IS,MDC+ID) = DA2M(IS,  ID   )
!
            SA1 (IS,  ID0 ) = SA1 (IS, MDC+ID0)
            SA2 (IS,  ID0 ) = SA2 (IS, MDC+ID0)
            DA1C(IS,  ID0 ) = DA1C(IS, MDC+ID0)
            DA1P(IS,  ID0 ) = DA1P(IS, MDC+ID0)
            DA1M(IS,  ID0 ) = DA1M(IS, MDC+ID0)
            DA2C(IS,  ID0 ) = DA2C(IS, MDC+ID0)
            DA2P(IS,  ID0 ) = DA2P(IS, MDC+ID0)
            DA2M(IS,  ID0 ) = DA2M(IS, MDC+ID0)
          ENDDO
        ENDDO
      ENDIF
!
!     *** Put source term together (To save space I=IS and J=ID ***
!     *** is used)                                              ***
!
      PI3   = (2. * PI)**3
      DO I = 1, ISSTOP
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = IDCMIN(I), IDCMAX(I)
          ID = MOD ( J - 1 + MDC , MDC ) + 1
          SFNL(I,ID) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
          DSNL(I,ID) =   - 2. * ( DA1C(I,J) + DA2C(I,J) )
     &        + SWG1 * ( DA1P(I-ISP1,J-IDP1) + DA2P(I-ISP1,J+IDP1) )
     &        + SWG2 * ( DA1P(I-ISP1,J-IDP ) + DA2P(I-ISP1,J+IDP ) )
     &        + SWG3 * ( DA1P(I-ISP ,J-IDP1) + DA2P(I-ISP ,J+IDP1) )
     &        + SWG4 * ( DA1P(I-ISP ,J-IDP ) + DA2P(I-ISP ,J+IDP ) )
     &        + SWG5 * ( DA1M(I-ISM1,J+IDM1) + DA2M(I-ISM1,J-IDM1) )
     &        + SWG6 * ( DA1M(I-ISM1,J+IDM ) + DA2M(I-ISM1,J-IDM ) )
     &        + SWG7 * ( DA1M(I-ISM ,J+IDM1) + DA2M(I-ISM ,J-IDM1) )
     &        + SWG8 * ( DA1M(I-ISM ,J+IDM ) + DA2M(I-ISM ,J-IDM ) )
!
!         *** store results in IMATDA and IMATRA ***
!
          IF(TESTFL) THEN
            PLNL4S(ID,I,IPTST) = SFNL(I,ID) / SIGPI                       40.00
            PLNL4D(ID,I,IPTST) = -1. * DSNL(I,ID) / PI3                   40.00
          END IF
!
          IMATRA(ID,I) = IMATRA(ID,I) + SFNL(I,ID) / SIGPI
          IMATDA(ID,I) = IMATDA(ID,I) - DSNL(I,ID) / PI3
!
          IF(ITEST.GE.90 .AND. TESTFL) THEN
            WRITE(PRINTF,9006) I,J,SFNL(I,ID),DSNL(I,ID),
     &       SPCSIG(I)                                                    30.72
 9006       FORMAT (' IS ID SFNL DSNL SPCSIG:',2I4,3E12.4)                30.72
          END IF
!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 50 .AND. TESTFL) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' SWSNL1 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW,IDHGH
 9015   FORMAT (' ISLOW ISHGH IDLOW IDHG:',4I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, IDDLOW, IDDTOP
 9016   FORMAT (' ICLW ICHG IDDLOW IDDTO:',2I5)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC, FACHFR, PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,9023) JACOBI
 9023   FORMAT (' JACOBI                :',E12.4)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of the subroutine SWSNL1
      END
!
!*******************************************************************
!
      SUBROUTINE SWSNL2 (IDDLOW  ,IDDTOP  ,WWINT   ,                      34.00
     &                   WWAWG   ,UE      ,SA1     ,ISSTOP  ,             40.17
     &                   SA2     ,SPCSIG  ,SNLC1   ,DAL1    ,DAL2    ,    30.72
     &                   DAL3    ,SFNL    ,DEP2    ,AC2     ,KMESPC  ,
     &                                              IMATDA  ,IMATRA  ,    40.23 34.00
     &                   FACHFR  ,PLNL4S           ,IDCMIN  ,IDCMAX  )    34.00
!
!*******************************************************************

      USE M_SNL4                                                          40.17
!
      INCLUDE 'swcomm3.inc'                                               34.00
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.13: Nico Booij
!     40.17: IJsbrand Haagsma
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implemented Multiple DIA
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!
!  3. Method
!
!     Discrete interaction approximation.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW = IDDLOW - MAX(IDM1,IDP1)
!     IDHGH = IDDTOP + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS     ,ID     ,I      ,J                      ,ISHGH  ,  34.00
     &          ISSTOP ,ISP    ,ISP1   ,IDP    ,IDP1   ,ISM    ,ISM1   ,
     &          IDM    ,IDM1   ,ISCLW  ,ISCHG  ,                          34.00
     &                  IDLOW  ,IDHGH  ,IDDLOW ,IDDTOP ,IDCLOW ,IDCHGH    34.00
!
      REAL      X      ,X2     ,CONS   ,FACTOR ,SNLCS1 ,SNLCS2 ,SNLCS3 ,
     &          E00    ,EP1    ,EM1    ,EP2    ,EM2    ,SA1A   ,SA1B   ,
     &          SA2A   ,SA2B   ,KMESPC ,FACHFR ,AWG1   ,AWG2   ,AWG3   ,
     &          AWG4   ,AWG5   ,AWG6   ,AWG7   ,AWG8   ,DAL1   ,DAL2   ,
     &          DAL3           ,JACOBI ,SIGPI                             34.00
!
      REAL      AC2(MDC,MSC,MCGRD)                    ,                   30.21
     &          DEP2(MCGRD)                           ,                   30.21
     &          UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,
     &          SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA)   ,
     &          IMATRA(MDC,MSC)                       ,
     &          IMATDA(MDC,MSC)                       ,                   40.23
     &          PLNL4S(MDC,MSC,NPTST)                 ,                   40.00
     &          WWAWG(*)
!
      INTEGER   WWINT(*)         ,
     &          IDCMIN(MSC)      ,
     &          IDCMAX(MSC)
!
      LOGICAL   PERCIR
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL2')
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***         40.17
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** check whether the spectral domain is periodic in ***
!     *** direction space and if so modify boundaries      ***
!
      PERCIR = .FALSE.
      IF ( IDDLOW .EQ. 1 .AND. IDDTOP .EQ. MDC ) THEN
!       *** periodic in theta -> spectrum can be folded  ***
!       *** (can only occur in presence of a current)    ***
        IDCLOW = 1
        IDCHGH = MDC
        IIID   = 0
        PERCIR = .TRUE.
      ELSE
!       *** different sectors per sweep -> extend range with IIID ***
        IIID   = MAX ( IDM1 , IDP1 )
        IDCLOW = IDLOW
        IDCHGH = IDHGH
      ENDIF
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
      DO IDDUM = IDLOW - IIID , IDHGH + IIID
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS = 1, MSC
          UE(IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI        30.72
        ENDDO
      ENDDO
!
!     *** set values in the areas 2 for IS > MSC+1 ***
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW - IIID , IDHGH + IIID
          UE (IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
      DO IS = ISCLW, ISCHG
        DO ID = IDCLOW , IDCHGH
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
!
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         *** Contribution to interactions                          ***
!         *** CONS is the shallow water factor for the NL interact. ***
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               40.17
          SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         40.17
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               40.17
          SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         40.17
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR ,ISLOW
 9005       FORMAT (' FACTOR ISLOW         : ',E12.4,I4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
      IF ( PERCIR ) THEN
        DO ID = 1, IDHGH - MDC
          ID0   = 1 - ID
          DO IS = ISCLW, ISCHG
            SA1 (IS,MDC+ID) = SA1 (IS ,  ID    )
            SA2 (IS,MDC+ID) = SA2 (IS ,  ID    )
            SA1 (IS,  ID0 ) = SA1 (IS , MDC+ID0)
            SA2 (IS,  ID0 ) = SA2 (IS , MDC+ID0)
          ENDDO
        ENDDO
      ENDIF
!
!     ***  Put source term together (To save space I=IS and J=ID ***
!     ***  is used)                                              ***
!
      DO I = 1, ISSTOP
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = IDCMIN(I), IDCMAX(I)
          ID = MOD ( J - 1 + MDC , MDC ) + 1
          SFNL(I,ID) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!         *** store results in rhv ***
!
          IF(TESTFL) PLNL4S(ID,I,IPTST) =  SFNL(I,ID) / SIGPI             40.00
          IMATRA(ID,I) = IMATRA(ID,I) + SFNL(I,ID) / SIGPI                40.00
!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 40 .AND. TESTFL) THEN
        WRITE(PRINTF,*) ' SWSNL2 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISHGH, IDDLOW, IDDTOP
 9015   FORMAT (' ISHG IDDLOW IDDTOP    :',3I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, IDLOW, IDHGH
 9016   FORMAT (' ICLW ICHG IDLOW IDHGH :',4I5)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC, FACHFR,PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,9023) JACOBI,ISLOW
 9023   FORMAT (' JACOBI  ISLOW         :',E12.4,I4)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of SWSNL2
      END
!
!************************************************************
!

      SUBROUTINE SWSNL3 (                  WWINT   ,WWAWG   ,             40.17
     &                   UE      ,SA1     ,SA2     ,SPCSIG  ,SNLC1   ,    40.17
     &                   DAL1    ,DAL2    ,DAL3    ,SFNL    ,DEP2    ,    40.17
     &                   AC2     ,KMESPC  ,MEMNL4  ,FACHFR           )    40.17
!
!*******************************************************************

      USE M_SNL4                                                          40.17
!
      INCLUDE 'swcomm3.inc'                                               40.17
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.17: IJsbrand Haagsma
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implemented Multiple DIA
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle (option if a current is present). Note: using
!     this subroutine requires an additional array with size
!     (MXC*MYC*MDC*MSC). This requires more internal memory but can
!     speed up the computations sigificantly if a current is present.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated an a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     ICMAX : number of points in computational stencil
!     KCGRD : grid address of points of computational stencil
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
      INTEGER WWINT(*)
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!     WWAWG : weight coefficients for the quadruplet interactions
!
      REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1                     40.17
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    MEMNL4(MDC,MSC,MCGRD)
      REAL    SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    WWAWG(*)
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,
     &          ISHGH   ,IDLOW   ,IDHGH   ,ISP     ,ISP1    ,
     &          IDP     ,IDP1    ,ISM     ,ISM1    ,IDM     ,IDM1    ,
     &          ISCLW   ,ISCHG
!
      REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS2  ,
     &          SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     ,
     &          SA1A    ,SA1B    ,SA2A    ,SA2B    ,
     &          AWG1    ,AWG2    ,AWG3    ,AWG4    ,AWG5    ,AWG6    ,
     &          AWG7    ,AWG8    ,
     &          JACOBI  ,SIGPI
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL3')
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***         40.17
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )               30.21
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** extend the area with action density at periodic boundaries ***
!
      DO IDDUM = IDLOW, IDHGH
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS=1, MSC
          UE (IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI       30.72
        ENDDO
      ENDDO
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW, IDHGH
          UE(IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate (unfolded) interactions ***
!     *** Energy at interacting bins        ***
!
      DO IS = ISCLW, ISCHG
        DO ID = 1, MDC
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         Contribution to interactions
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               40.17
          SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         40.17
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               40.17
          SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         40.17
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR,JACOBI
 9005       FORMAT (' FACTOR JACOBI        : ',2E12.4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles -> domain in theta is ***
!     *** periodic                                               ***
!
      DO ID = 1, IDHGH - MDC
        ID0   = 1 - ID
        DO IS = ISCLW, ISCHG
          SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
          SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
          SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
          SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
        ENDDO
      ENDDO
!
!     *** Put source term together (To save space I=IS and ***
!     *** J=MDC is used)  ----                             ***
!
      DO I = 1, MSC
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = 1, MDC
          SFNL(I,J) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!         *** store value in auxiliary array and use values in ***
!         *** next four sweeps (see subroutine FILSNL3)        ***
!
          MEMNL4(J,I,KCGRD(1)) = SFNL(I,J) / SIGPI                        30.21
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 50 .AND. TESTFL) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' SWSNL3 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW, IDHGH
 9015   FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, JACOBI
 9016   FORMAT (' ICLW ICHG JACOBI      :',2I5,E12.4)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC,FACHFR,PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,*)
!
!       *** value source term in every bin ***
!
        IF(ITEST.GE. 150 ) THEN
          DO I=1, MSC
            DO J=1, MDC
              WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),      30.21
     &                           SPCSIG(I)                                30.72
 2006         FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)           30.72
            ENDDO
          ENDDO
        END IF
      END IF
!
      RETURN
!
      END SUBROUTINE SWSNL3
!
!*******************************************************************
!
      SUBROUTINE SWSNL4 (WWINT   ,WWAWG   ,
     &                   SPCSIG  ,SNLC1   ,
     &                   DAL1    ,DAL2    ,DAL3    ,DEP2    ,
     &                   AC2     ,KMESPC  ,MEMNL4  ,FACHFR  ,
     &                   IDIA    ,ITER    )
!
!*******************************************************************

      USE M_SNL4

      INCLUDE 'swcomm3.inc'
      INCLUDE 'swcomm4.inc'
      INCLUDE 'ocpcomm4.inc'
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
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
!     40.17: IJsbrand Haagsma
!
!  1. Updates
!
!     40.17, Dec. 01: New Subroutine based on SWSNL3
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle (option if a current is present). Note: using
!     this subroutine requires an additional array with size
!     (MXC*MYC*MDC*MSC). This requires more internal memory but can
!     speed up the computations sigificantly if a current is present.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated an a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     ICMAX : number of points in computational stencil
!     KCGRD : grid address of points of computational stencil
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
      INTEGER WWINT(*)
      INTEGER IDIA
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!     WWAWG : weight coefficients for the quadruplet interactions
!
      REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    MEMNL4(MDC,MSC,MCGRD)
      REAL    SPCSIG(MSC)
      REAL    WWAWG(*)
!
!  6. Local variables
!
      REAL, ALLOCATABLE :: SA1(:,:), SA2(:,:), SFNL(:,:), UE(:,:)
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,
     &          ISHGH   ,IDLOW   ,IDHGH   ,ISP     ,ISP1    ,
     &          IDP     ,IDP1    ,ISM     ,ISM1    ,IDM     ,IDM1    ,
     &          ISCLW   ,ISCHG
!
      REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS2  ,
     &          SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     ,
     &          SA1A    ,SA1B    ,SA2A    ,SA2B    ,
     &          AWG1    ,AWG2    ,AWG3    ,AWG4    ,AWG5    ,AWG6    ,
     &          AWG7    ,AWG8    ,
     &          JACOBI  ,SIGPI
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL4')

      ALLOCATE(SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
      ALLOCATE(SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
      ALLOCATE(SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
      ALLOCATE(UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** extend the area with action density at periodic boundaries ***
!
      DO IDDUM = IDLOW, IDHGH
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS=1, MSC
          UE (IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI
        ENDDO
      ENDDO
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW, IDHGH
          UE(IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate (unfolded) interactions ***
!     *** Energy at interacting bins        ***
!
      DO IS = ISCLW, ISCHG
        DO ID = 1, MDC
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         Contribution to interactions
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * CNL4_1(IDIA)
          SA1B   = SA1A - EP1*EM1*DAL3 * CNL4_2(IDIA)
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * CNL4_1(IDIA)
          SA2B   = SA2A - EP2*EM2*DAL3 * CNL4_2(IDIA)
!

          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR,JACOBI
 9005       FORMAT (' FACTOR JACOBI        : ',2E12.4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles -> domain in theta is ***
!     *** periodic                                               ***
!
      DO ID = 1, IDHGH - MDC
        ID0   = 1 - ID
        DO IS = ISCLW, ISCHG
          SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
          SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
          SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
          SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
        ENDDO
      ENDDO
!
!     *** Put source term together (To save space I=IS and ***
!     *** J=MDC is used)                                   ***
!
      FAC = 1.                                                            40.17

      DO I = 1, MSC
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = 1, MDC
          SFNL(I,J) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!         *** store value in auxiliary array and use values in ***
!         *** next four sweeps (see subroutine FILSNL3)        ***
!
          IF (IDIA.EQ.1) THEN
            MEMNL4(J,I,KCGRD(1)) = FAC * SFNL(I,J) / SIGPI
          ELSE
            MEMNL4(J,I,KCGRD(1)) = MEMNL4(J,I,KCGRD(1)) +
     &                             FAC * SFNL(I,J) / SIGPI
          END IF
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 50 .AND. TESTFL) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' SWSNL4 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW, IDHGH
 9015   FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, JACOBI
 9016   FORMAT (' ICLW ICHG JACOBI      :',2I5,E12.4)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC,FACHFR,PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,*)
!
!       *** value source term in every bin ***
!
        IF(ITEST.GE. 150 ) THEN
          DO I=1, MSC
            DO J=1, MDC
              WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),
     &                           SPCSIG(I)
 2006         FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)
            ENDDO
          ENDDO
        END IF
      END IF
!
      DEALLOCATE (SA1, SA2, SFNL, UE)

      RETURN
!
      END SUBROUTINE SWSNL4
!
!*******************************************************************
!
      SUBROUTINE FILNL3 (MDC     ,MSC     ,IDCMIN  ,IDCMAX  ,IMATRA  ,
     &                                              IMATDA  ,AC2     ,    40.23
     &                   MEMNL4  ,PLNL4S  ,ISSTOP  ,KCGRD   ,MCGRD   ,
     &                   ICMAX                                       )    30.21
!
!*******************************************************************
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
!  0. Authors
!
!     ??.??: ?  ?
!
!  1. Updates
!
!
!  2. Purpose
!
!     Fill the IMATRA/IMATDA arrays with the nonlinear wave-wave interaction
!     source term for a gridpoint ix,iy per sweep direction
!
!  3. Method
!
!
!  4. Argument variables
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!     Do for every frequency and spectral direction within a sweep
!         fill IMATRA/IMATDA
!     -------------------------------------------
!     End of FILNL3
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,MSC     ,MDC     ,ISSTOP   ,MCGRD   ,
     &          ICMAX                                                     30.21
!
      INTEGER   KCGRD(ICMAX)                                              30.21
!
      REAL      IMATRA(MDC,MSC)           ,
     &          IMATDA(MDC,MSC)           ,                               40.23
     &          AC2(MDC,MSC,MCGRD)        ,                               40.23
     &          PLNL4S(MDC,MSC,NPTST)     ,                               40.00
     &          MEMNL4(MDC,MSC,MCGRD)                                     30.21
!
      INTEGER   IDCMIN(MSC)         ,
     &          IDCMAX(MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'FILNL3')
!
      DO 990 IS=1, ISSTOP
        DO 980 IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IF(TESTFL) PLNL4S(ID,IS,IPTST) = MEMNL4(ID,IS,KCGRD(1))         40.00
          IMATRA(ID,IS) = IMATRA(ID,IS) + MEMNL4(ID,IS,KCGRD(1))          30.21
  980   CONTINUE
  990 CONTINUE
!
      IF ( TESTFL .AND. ITEST.GE.50 ) THEN
        WRITE(PRINTF,9000) IDCMIN(1),IDCMAX(1),MSC,ISSTOP
 9000   FORMAT(' FILNL3: ID_MIN ID_MAX MSC ISTOP :',4I6)
        IF ( ITEST .GE. 100 ) THEN
          DO IS=1, ISSTOP
            DO IDDUM = IDCMIN(IS), IDCMAX(IS)
              ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
              WRITE(PRINTF,6001) IS,ID,MEMNL4(ID,IS,KCGRD(1))             30.21
 6001         FORMAT(' FILNL3: IS ID MEMNL()          :',2I6,E12.4)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
!
      RETURN
!     End of FILNL3
      END
!
!********************************************************************
!
      SUBROUTINE STRIAD (AC2     ,DEP2    ,CGO     ,IMATRA  ,KWAVE   ,
     &                   HS      ,IDDLOW  ,IDDTOP  ,
     &                   SPCSIG  ,SMEBRK  ,IMATDA  ,PLTRI   ,URSELL  )    40.03
!
!********************************************************************
!
      INCLUDE 'swcomm3.inc'                                               30.80
      INCLUDE 'swcomm4.inc'                                               30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: Y. Eldeberky, R.C. Ris                       |
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     30.50         : EXP is renamed TRIEXP to avoid confusion with
!                     standard function EXP; EXP is made FALSE
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.80, Aug. 99: include file swcomm3.inc added
!     40.03, Apr. 00: computation of Ursell moved to subr SDISPA
!     40.23, Sep. 02: small correction (0.1 replaced by PTRIAD(5))
!
!  2. Purpose
!
!     Calculates triad wave interaction based on a Boussinesq-type
!     wave equation
!
!  3. Method
!
!     Original code from Y. Eldeberky (17-03-1995). Subroutine has been
!     recoded and optimized by R.C. Ris (12-09-1995). The procedure has
!     been modified such that an explicit or an implicit scheme for the
!     triad interaction source term can be used.
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     INTEGER
!     -------
!     i             Frequency counter at central bin
!     j             Frequency counter at (higher and lower) harmonic bin
!     MSC,MDC       Range of the original arrays
!     IDDLOW/IDDTOP minimum, maximum counter in directional space
!     IRES          counter that represents the resonance condition in
!                   discrete freq. grid. Is a constant for logaritmic
!                   distribution
!     ISMAX         maximum counter for which the triad formulation
!                   is valid and should be applied
!
!     REALS:
!     ------
!     DPI           two times PI
!     Wi, Wj        radian frequency at central and interacting freq.
!                   bin, respectively
!     WNi,WNj       wavenumber at central and interacting freq.
!                   bin, respectively
!     CGi,CGj       Group velocity at central and interacting freq.
!                   bin, respectively
!     SMEBRK        mean frequency according to first oder moment
!                   (see subroutine SDISPA)
!     XISTRI        Rate between two succeding freq. counters
!
!     ARRAY:
!     ------
!     DEP2          depth at ix,iy
!     KWAVE         wave number
!     CGO           group velocity
!     AC2           action density as function of d,s,x,y at time t
!     E             energy density as functiion of d,f,x,y
!     IMATRA        right hand vector
!     IMATDA        diagonal
!
!  4. Subroutines used :
!
!     ---
!
!  5. Called by :
!
!     SOURCE
!
!  6. Error messages :
!
!     ---
!
!  7. Source code :
!
!     -----------------------------------------------------------------
!     Calculate :
!
!     -----------------------------------------------------------
!     End of the subroutine STRIAD
!     ----------------------------------------------------------
!
!*************************************************************
!
      INTEGER   IS      ,ID      ,                           IDDLOW  ,
     &          IDDTOP  ,         IDDUM   ,IRES    ,I       ,J       ,
     &          ISMAX   ,I1      ,I2                                      30.21
!
      REAL               DEP     ,DEP_2   ,DEP_3   ,         DPI     ,
     &          FT      ,RHV_i   ,RHV_j   ,DIA_i   ,DIA_j   ,B       ,
     &          ALPHA   ,BETA    ,XISTRI  ,SMEBRK  ,         BIPH    ,
     &          SINBPH  ,JACi    ,JACj    ,Wi      ,Wj      ,WNi     ,
     &          WNj     ,CGi     ,CGj     ,HS      ,TMN
!
      REAL  ::  AC2(MDC,MSC,MCGRD)
      REAL  ::  DEP2(MCGRD)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  ::  KWAVE(MSC,MICMAX)                                         40.22
      REAL  ::  IMATRA(MDC,MSC)
      REAL  ::  IMATDA(MDC,MSC)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  ::  CGO(MSC,MICMAX)                                           40.22
      REAL  ::  PLTRI(MDC,MSC,NPTST)
      REAL  ::  URSELL(MCGRD)                                             40.03
      REAL  ::  E(1:MSC)                                                  40.22
!
      LOGICAL   TRIEXP
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'STRIAD')
!
!     *** initialization of variables ***
!
      DPI   = 2. * PI
      DEP   = DEP2(KCGRD(1))                                              30.21
      DEP_2 = DEP**2
      DEP_3 = DEP**3
      B     = 1. / 15.
      IF (ITRIAD.EQ.1) THEN
!       implicit computation
        TRIEXP = .FALSE.
      ELSE
!       explicit computation
        TRIEXP = .TRUE.
      ENDIF
!
!     *** determine resonance condition and determine maximum  ***
!     *** discrete counter for 3. "peak" frequency in Hz       ***
!
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XISTRI = SPCSIG(I2) / SPCSIG(I1)                                    30.72
      IRES   = NINT ( LOG( 2.) / LOG ( XISTRI ) )
!
      ISMAX = 1                                                           20.88
      DO IS = 1, MSC
        IF ( SPCSIG(IS) .LT. ( PTRIAD(2) * SMEBRK) ) THEN                 30.72
          ISMAX = IS
        ENDIF
      ENDDO
      ISMAX = MAX ( ISMAX , IRES + 1 )
!
      TMN = DPI / SMEBRK
      URSELL(KCGRD(1)) = MIN ( URSELL(KCGRD(1)) , 10. )
      IF (URSELL(KCGRD(1)) .LT. PTRIAD(5)) THEN
        BIPH = 0.0
      ELSE
        BIPH = - DPI / 8. * ( LOG10(URSELL(KCGRD(1))) + 1.)
      ENDIF
      SINBPH = ABS( SIN(BIPH) )
!
      IF ( URSELL(KCGRD(1)) .GE. PTRIAD(5) ) THEN
        DO IDDUM = IDDLOW, IDDTOP
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
          IF ( ITEST .GE. 10 .AND. TESTFL) THEN
            WRITE(PRINTF,1001) ID
 1001       FORMAT (' directional counter for interacting triad: ',2I4)
          ENDIF
!
!         *** initialize array with E(f) for a direction considered ***
!
          DO IS = 1, MSC
            E(IS)  = AC2(ID,IS,KCGRD(1)) * DPI * SPCSIG(IS)               30.72
          ENDDO
!
!         *** (j) corresponds to central grid point           ***
!         *** (i) corresponds to grid point lower harmonic    ***
!         ***                                                 ***
!         *** for log. distr.:                                ***
!         ***                                                 ***
!         ***       j =  i + ires                             ***
!         ***      <------------>                             ***
!         ***                                                 ***
!         ***     i             j                             ***
!         ***  ---+-------------+------------+--------        ***
!         ***    fp/2           fp                            ***
!         ***                                                 ***
!         ***                                                 ***
!         ***  start at i=1   ------> ISMAX-IRES              ***
!         ***                                                 ***
!         ***  at central point (j) :                         ***
!         ***  ----------------------                         ***
!         ***  Wj    radian frequency                         ***
!         ***  WNj   wave number                              ***
!         ***  CGj   group velocity                           ***
!         ***  JACj  Jacobian function A(s) = E(f) / 2.pi.s   ***
!         ***                                                 ***
!
          DO I = 1, ISMAX-IRES
!           *** j = central , i = harmonic ***
            J   = I + IRES
            Wi  = SPCSIG(I)                                               30.72
            Wj  = SPCSIG(J)                                               30.72
            WNi = KWAVE(i,1)
            WNj = KWAVE(j,1)
            CGi = CGO(i,1)
            CGj = CGO(j,1)
            JACi = DPI * Wi
            JACj = DPI * Wj
!
            ALPHA = 4. * WNi**2 *
     &              ( 0.5 + ( Wi**2 / ( WNi**2 * GRAV * DEP ) ) )
!
            BETA  = -2. * WNj * ( GRAV * DEP +
     &                            2. * B * GRAV * DEP_3 * WNj**2 -
     &                          ( B + 1./3. ) * Wj**2 * DEP_2   )
!
!           *** constant FT, PTRIAD(1) controls the intensity ***
!
            FT = PTRIAD(1) * CGj * SINBPH * ( GRAV * ALPHA / BETA )**2
!
!           *** explicit or implicit calculation of source term ***
!
            RHV_i = 0.
            RHV_j = 0.
            DIA_i = 0.
            DIA_j = 0.
!
            IF ( TRIEXP ) THEN
!
!             *** explicit calculation ***
!
              RHV_j = FT * (     (Wj/WNj) * E(i) * E(i)  -
     &                      2. * (Wi/WNi) * E(j) * E(i)    )
              RHV_i = RHV_j
!
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
              ENDIF
!
!             *** multiply source term Se(f) with Jacobian ***
!
              IMATRA(ID,i) = IMATRA(ID,i) - RHV_i / JACi
              IMATRA(ID,j) = IMATRA(ID,j) + 0.5 * RHV_j / JACj
!
            ELSE
!
!             *** semi-implicit calculation ***
!
!             *** use explicit scheme for point J since there is ***
!             *** grow at the higher frequencies                 ***
!
              RHV_j = FT * (     (Wj/WNj) * E(i) * E(i)  -
     &                      2. * (Wi/WNi) * E(j) * E(i)    )
!
!             *** source term at point i : E(i) is unknown           ***
!             *** at this point there is dissipation -> use implicit ***
!             *** scheme                                             ***
!
              DIA_i = FT * (     (Wj/WNj) * E(i)         -
     &                      2. * (Wi/WNi) * E(j)           )
!
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
                DIA_j = 0.
                DIA_i = 0.
              ENDIF
!
!             *** source term at point j : E(j) is unknown ***   deleted
!
!              RHV_j = FT * (       (Wj/WNj) * E(i) * E(i)  )     deleted
!              DIA_j = FT * ( -2. * (Wi/WNi) *      * E(i)  )     deleted
!              IMATRA(ID,j) = IMATRA(ID,j) + 0.5 * RHV_j / JACj   deleted
!              IMATDA(ID,j) = IMATDA(ID,j) - 0.5 * DIA_j          deleted
!
              IMATRA(ID,j) = IMATRA(ID,j) + 0.5 * RHV_j / JACj
              IMATDA(ID,i) = IMATDA(ID,i) + DIA_i
!
            ENDIF
!
!           *** store results in array for plot of triad term ***
!
            IF( TESTFL ) THEN
              PLTRI(ID,I,IPTST) = PLTRI(ID,I,IPTST)                       40.00
     &           - RHV_i / JACi - DIA_i * AC2(ID, i ,KCGRD(1))
!
              PLTRI(ID,J,IPTST) = PLTRI(ID,J,IPTST)                       40.00
     &        + 0.5 * RHV_j / JACj + 0.5 * DIA_j * AC2(ID,j,KCGRD(1))     30.21
            END IF
          ENDDO
!
!         *** test output for a particular direction ***
!
          IF ( ITEST .GE. 10 .AND. TESTFL) THEN
            WRITE(PRINTF,1002) ALPHA, BETA, FT, i, j
 1002       FORMAT (' STRIAD: ALPH BETA FT i j :',3E12.4,2X,2I3)
            WRITE(PRINTF,1003) Wi, WNi, CGi, JACi
 1003       FORMAT (' STRIAD: Wi WNi CGi JACi  :',4E12.4)
            WRITE(PRINTF,1004) Wj, WNj, CGj, JACj
 1004       FORMAT (' STRIAD: Wj WNj CGj JACj  :',4E12.4)
            WRITE(PRINTF,1005) RHV_i, RHV_j, DIA_i, DIA_j
 1005       FORMAT (' STRIAD: RHV_i,j DIA_i,j  :',4E12.4)
          ENDIF
!
        ENDDO
      ENDIF
!
!     *** test output ***
!
      IF ( ITEST .GE. 10 .AND. TESTFL) THEN
         WRITE(PRINTF,2000) KCGRD(1), IRES, ISMAX
 2000    FORMAT (' STRIAD: POINT IRES ISMAX :',3I5)
         WRITE(PRINTF,2001) GRAV, DEP, DEP_2, DEP_3
 2001    FORMAT (' STRIAD: G DEP DEP2 DEP3  :',4E12.4)
         WRITE(PRINTF,2002) DPI, PTRIAD(1), B, URSELL(KCGRD(1))
 2002    FORMAT (' STRIAD: DPI PAR B URSELL :',4E12.4)
         WRITE(PRINTF,2003) SMEBRK, TMN, HS, BIPH
 2003    FORMAT (' STRIAD: SMEBRK TMN HS BIPH:',4E12.4)
      ENDIF
!
!     *** end of subroutine STRIAD ***
!
      RETURN
      END
!
!********************************************************************
!
      SUBROUTINE STRIAN (AC2   ,DEP2  ,CGO   ,IMATRA,KWAVE ,
     &                   HS    ,IDDLOW,IDDTOP,
     &                   SPCSIG,SMEBRK,IMATDA,PLTRI ,URSELL        )      40.03
!
!********************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'swcomm3.inc'                                               40.00
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     40.00, 40.13: Nico Booij
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!            Jan. 97: New subroutine (Roeland Ris)
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.81, Aug. 98: Modified calculation of the Ursell number
!     30.82, Sep. 98: Declaration of argument variable reorganised
!     40.03, Apr. 00: computation of Ursell moved to subr SDISPA
!     40.13, July 01: var. coeff. PTRIAD(3) [urslim] introduced
!                     var. coeff. PTRIAD(4) [urcrit] introduced
!                     array E(:) is made allocatable
!     40.23, Sep. 02: PTRIAD(3) should be PTRIAD(5)
!
!  2. Purpose
!
!     In this subroutine the triad-wave interactions are calculated
!     with the Lumped Triad Approximation of Eldeberky (1996). His
!     expression is based on a parametrization of the bi-phase (in
!     terms of the Ursell number), is directionally uncoupled and
!     takes into account for the self-self interactions only.
!
!     For a full description of the equations reference is made
!     to Eldeberky (1996). Here only the main equations are given.
!
!  3. Method
!
!     The Ursell number is given by (computed in routine SINTGRL):
!
!                  g         Hs Tm**2         g          Hs
!     UR =  --------------- --------- = ---------- -------------
!           8 sqrt(2) pi**2   d**2      2 sqrt (2) sig_m**2 d**2
!
!     in which:
!
!     d      depth (m)
!     g      gravitational acceleration (m/s^2)
!     Hs     significant wave height (m)
!     sig_m  mean wave frequency (in rad/s) according to m1/m0
!            (personal communication, Y. Eldeberky, 1997)
!     Tm     mean wave period (s)
!
!     The biphase BIPH is given by (see eq. 3.19, Eldeberky, 1996 ):
!
!                                  0.2
!     BIPH = - pi/2 + pi/2 tanh ( ----- )
!                                   Ur
!
!     The source term is (see eq. 7.25):
!
!             +      -
!     S(p) = S(p) + S(p)
!
!     in which
!      +
!     S(p) = a Cp Cg,p (R(p/2,p/2))**2 sin (|BIPH|) ( E(p/2)**2 -
!
!                          2 E(p) E(p/2)  )
!      -          +
!     S(p) = - 2 S(2p)
!
!     in which:
!     a    tunable coefficient
!     p    frequency
!
!     The value for the interaction coefficient R(p/2,p/2) is
!     not given here (see eq 7.26 in Eldeberky, 1996)
!
!     Note that a slightly adapted formulation of the LTA is used in
!     in the SWAN model:
!
!     - only positive contributions to higher harmonics are considered
!       here (no energy is transferred to lower harmonics).
!
!     - the mean frequency in the expression of the Ursell number
!       is calculated according to the first order moment over the
!       zeroth order moment (personnal communucation, Y.Eldeberky, 1997)
!
!     - the interactions are calculated up to 2.5 times the mean
!       frequency only.
!
!     - to avoid expensive interpolations between discrete bins the
!       interactions are taken at discrete gridpoints only. Since the
!       grid is logarimically distributed in frequency space the
!       number of bins between the central bin and lower harmonic bin
!       is constant (equal to IRES)
!
!     Moreover:
!
!     - the interactions are calculated in terms of energy density
!       instead of action density. So the action density spectrum
!       is firstly converted to the energy density grid, then the
!       interactions are calculated and then the spectrum is converted
!       to the action density spectrum again.
!
!     - To ensure numerical stability the positive contributions
!       of S(p) are calculated with an explicit scheme and the
!       negative contributions with an implicit scheme. Optionally
!       a fully explicit scheme is available.
!
!  4. Argument variables
!
!     ICMAX                                                               30.82
!     IDDLOW: Minimum counter in directional space (see subr. COUNT)??    30.82
!     IDDTOP: Maximum counter in dierctional space (see subr. COUNT)??    30.82
!     ITRIAD                                                              30.82
!     KCGRD                                                               30.82
!     MCGRD                                                               30.82
!     MDC   : Size of array in t-direction??                              30.82
!     MSC   : Size of array in s-direction??                              30.82
!     MTRIAD: Size of array containing triad-coefficients                 30.82
!
      INTEGER        IDDLOW, IDDTOP                                       30.82
!
!     AC2   : Action density as function of d,s,x,y at time t             30.82
!     CGO   : Group velocity                                              30.82
!     DEP2  : Depth at gridpoint ix,iy (obtained from SWANCOM1??)         30.82
!     GRAV                                                                30.82
!     HS                                                                  30.82
!     IMATRA: Right hand vector of matrix??                               30.82
!     IMATDA: Diagonal of matrix??                                        30.82
!     KWAVE : Wave number                                                 30.82
!     PI    : Constant Pi                                                 30.82
!     PLTRI : Values of the triad source terms in TEST-points             30.82
!     PTRIAD: Tunable coefficients for nonlinear sourceterms              30.82
!     SMEBRK: Mean frequency (see subroutine SDISPA)                      30.82
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL :: HS      ! significant wave height                           30.82
      REAL :: SMEBRK  ! average (angular) frequency                       30.82
      REAL :: AC2(MDC,MSC,MCGRD)                                          30.82
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL :: CGO(MSC,MICMAX)                                             30.82 40.22
      REAL :: DEP2(MCGRD)                                                 30.82
      REAL :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)                            30.82
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL :: KWAVE(MSC,MICMAX)                                           30.82 40.22
      REAL :: PLTRI(MDC,MSC,NPTST),                 SPCSIG(MSC)           30.82
      REAL :: URSELL(MCGRD)                                               40.03
!
!  6. Local variables
!
!     I     : Frequency counter at bin of lower harmonic (f/2)            30.82
!             and subscript for variable??                                30.82
!     I1    : Temporary variable to determine the resonance conditions??  30.82
!     I2    : Temporary variable to determine the resonance conditions??  30.82
!     ID    : Grid counter in spectral space (direction)                  30.82
!     IENT  : Number of entries into this subroutine                      30.82
!     II    : Counter                                                     30.82
!     IRES  : Number of bins between central and lower harm. freq.        30.82
!             (is constant for a logaritmic frequency distribution)       30.82
!     IS    : Grid counter in spectral space (frequency)                  30.82
!     ISMAX : Maximum of the counter in spectral space (frequency) for    30.82
!             which the triad interactions are calculated (cut-of)        30.82
!     IX    : Grid counter in geographical space (x-direction)            30.82
!     IY    : Grid counter in geographical space (y-direction)            30.82
!     J     : Frequency counter at central gridpoint (f)                  30.82
!             and subscript for variable??                                30.82
!
      INTEGER I, I1, I2, ID, IENT, II, IRES, IS, ISMAX, J                 30.82
!
!     AUX1  : Temporary variables                                         30.82
!     AUX2  : Temporary variables                                         30.82
!     BIPH  : Parameterized bi-phase of the spectrum                      30.82
!     CGI   : Group velocity at interacting bin                           30.82
!     CGJ   : Group velocity at central frequency bin                     30.82
!     CI    : Phase velocity at interacting bin                           30.82
!     CJ    : Phase velocity at central frequency bin                     30.82
!     DEP   : =DEP2(IX,IY); Depth at gridpoint (IX,IY)                    30.82
!     DEP_2 : =DEP**2;                                                    30.82
!     DEP_3 : =DEP**3;                                                    30.82
!     DPI   : =2*PI;                                                      30.82
!     E     : Energy density as function of location, freq. and direction 30.82
!     FT    : Coefficient??                                               30.82
!     JACI       Jacobian function at bin i and j respectively   ??       30.82
!     JACJ       equal to :  A(s) = E(s) / 2.pi.s                ??       30.82
!     URSELL: Ursell number                                               30.82
!     WI    : Radian frequency at interacting bin                         30.82
!     WJ    : Radian frequency at central frequency bin                   30.82
!     WNI   : Wavenumber at interacting bin                               30.82
!     WNJ   : Wavenumber at central freqency bin                          30.82
!     RINT  : Interaction coefficient                                     30.82
!     XISTRI: Rate between two succeding freq. counters??                 30.82
!
      REAL    AUX1, AUX2, BIPH, CGI, CGJ, CI, CJ, DEP, DEP_2, DEP_3, DPI  30.82
      REAL :: E(1:MSC)                                                    40.22
      REAL    FT, JACI, JACJ, WI, WJ, WNI, WNJ, RINT, XISTRI              30.82
!
!     TRIEXP: Explicit integration scheme is used for calculating triads  30.82
!             (semi-implicit when false)                                  30.82
!
      LOGICAL TRIEXP                                                      30.82
!
!  7. SUBROUTINES USED
!
!     ---
!
!  8. SUBROUTINES CALLING
!
!     SOURCE
!
!  9. ERROR MESSAGES
!
!     ---
!
! 10. REMARKS
!
!     ---
!
! 11. STRUCTURE
!
!     -----------------------------------------------------------------
!     - Initialize variables
!     - Determine resonance condition and the maximum discrete freq.
!       for which the interactions are calculated.
!     -----------------------------------------------------------------
!     If Ursell number larger than 0.1 compute interactions
!       Calculate Biphase
!       Do for each direction
!         Convert action density to energy density
!         Do for all frequencies
!           Calculate interaction coefficient and interaction factor
!           Compute interactions and store results in matrix
!         Enddo
!         -------------------------------------------------------------
!         If testfl store results in test array
!       Enddo
!       ---------------------------------------------------------------
!
! 13. Source text
!
      REAL      RHV_i ,RHV_j ,DIA_i ,DIA_j , SINBPH
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'STRIAN')
!
!     *** initialization of variables ***
!
      DPI   = 2. * PI
      DEP   = DEP2(KCGRD(1))
      DEP_2 = DEP**2
      DEP_3 = DEP**3
      IF (ITRIAD.EQ.3) THEN
!       implicit computation
        TRIEXP = .FALSE.
      ELSE
!       explicit computation
        TRIEXP = .TRUE.
      ENDIF
!
      E(:) = 0.                                                           40.22
!
!     *** determine "resonance condition" and determine maximum    ***
!     *** discrete counter (ISMAX) at 2.5 times the mean frequency ***
!     *** value of 2.5 can be varied using PTRIAD(2)               ***
!
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XISTRI = SPCSIG(I2) / SPCSIG(I1)                                    30.72
      IRES   = NINT ( LOG( 2.) / LOG ( XISTRI ) )
!
      ISMAX = 1
      DO IS = 1, MSC
       IF ( SPCSIG(IS) .LT. ( PTRIAD(2) * SMEBRK) ) THEN                  30.72
          ISMAX = IS
        ENDIF
      ENDDO
      ISMAX = MAX ( ISMAX , IRES + 1 )
!
!     Note: Ursell number is computed in routine SINTGRL                  40.13
!
      IF ( URSELL(KCGRD(1)) .GE. PTRIAD(5) ) THEN                         40.13
!
!       *** calculate Biphase ***
!
        BIPH   = (0.5*PI)*(TANH(PTRIAD(4)/URSELL(KCGRD(1)))-1)            40.13
        SINBPH = ABS( SIN(BIPH) )
!
        DO II = IDDLOW, IDDTOP
          ID = MOD ( II - 1 + MDC , MDC ) + 1
!
!         *** initialize array with E(f) for the direction considered ***
!
          DO IS = 1, MSC
            E(IS)  = AC2(ID,IS,KCGRD(1)) * DPI * SPCSIG(IS)               30.72
          ENDDO
!
!         *** (i) corresponds to grid point lower harmonic    ***
!         *** (j) corresponds to central grid point           ***
!         ***                                                 ***
!         *** for log. distribution:                          ***
!         ***                                                 ***
!         ***       j =  i + ires                             ***
!         ***      <------------>                             ***
!         ***                                                 ***
!         ***     i             j                             ***
!         ***  ---+-------------+------------+--------        ***
!         ***    fp/2           fp                            ***
!         ***                                                 ***
!         ***  start at i=1   ------> ISMAX-IRES              ***
!         ***                                                 ***
!
          DO I = 1, ISMAX-IRES
!           *** j = central bin, i = bin at lower harmonic ***
            J   = I + IRES
            Wi  = SPCSIG(I)                                               30.72
            Wj  = SPCSIG(J)                                               30.72
            WNi = KWAVE(i,1)
            WNj = KWAVE(j,1)
            Ci  = Wi / WNi
            Cj  = Wj / WNj
            CGi = CGO(i,1)
            CGj = CGO(j,1)
            JACi = DPI * Wi
            JACj = DPI * Wj
!
            AUX1 = WNi**2 * ( GRAV * DEP + 2. * Ci**2 )
            AUX2 = WNj * DEP * ( GRAV * DEP +
     &                          (2./15.) * GRAV * DEP_3 * WNj**2 -        30.82
     &                          (2./ 5.) * Wj**2 * DEP_2 )                30.82
            RINT = AUX1 / AUX2
!
!           *** constant FT, PTRIAD(1) controls the intensity ***
!
            FT = PTRIAD(1) * Cj * CGj * RINT**2 * SINBPH
!
!           *** explicit or implicit calculation of source term ***
!
            RHV_i = 0.
            RHV_j = 0.
            DIA_i = 0.
            DIA_j = 0.
!
            IF ( TRIEXP ) THEN
!
!             *** explicit calculation ***
!
              RHV_j = FT * (  E(i) * E(i)  - 2. * E(j) * E(i) )
              RHV_i = RHV_j
!
!             *** consider only transfer from lower to higher harm ***
!
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
              ENDIF
!
!             *** multiply source term Se(f) with Jacobian ***
!
              IMATRA(ID,i) = IMATRA(ID,i) - 2. * RHV_i / JACi
              IMATRA(ID,j) = IMATRA(ID,j) + RHV_j / JACj
!
            ELSE
!
!             *** semi-implicit calculation ***
!
!             *** use explicit scheme for point J since there is ***
!             *** grow at the higher frequencies                 ***
!
              RHV_j = FT * ( E(i) * E(i) - 2. * E(j) * E(i) )
!
!             *** source term at point i : E(i) is unknown           ***
!             *** at this point there is dissipation -> use implicit ***
!             *** scheme                                             ***
!
              DIA_i = FT * ( E(i) - 2. * E(j) )
!
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
                DIA_j = 0.
                DIA_i = 0.
              ENDIF
!
              IMATRA(ID,j) = IMATRA(ID,j) + RHV_j / JACj
              IMATDA(ID,i) = IMATDA(ID,i) + 2. * DIA_i
!
            ENDIF
!
!           *** store results in array for plot of triad term ***
!
            IF( TESTFL ) THEN
              PLTRI(ID,I,IPTST) = PLTRI(ID,I,IPTST)                       40.00
     &           - 2. * RHV_i / JACi - 2. * DIA_i * AC2(ID, i ,KCGRD(1))
!
              PLTRI(ID,J,IPTST) = PLTRI(ID,J,IPTST)                       40.00
     &          + RHV_j / JACj  + DIA_j * AC2(ID, j ,KCGRD(1))
            END IF
          ENDDO
!
!         *** test output for a particular direction ***
!
          IF ( ITEST .GE. 5 .AND. TESTFL) THEN
            WRITE(PRINTF,1002) AUX1, AUX2, FT, i, j
 1002       FORMAT (' STRIAN: ALPH BETA FT i j  :',3E12.4,2X,2I3)
            WRITE(PRINTF,1003) Wi, WNi, Ci, CGi, JACi
 1003       FORMAT (' STRIAN: i, W WN C CG JAC:',5E12.4)
            WRITE(PRINTF,1004) Wj, WNj, Cj, CGj, JACj
 1004       FORMAT (' STRIAN: j, W WN C CG JAC:',5E12.4)
            WRITE(PRINTF,1005) RINT, RHV_i, RHV_j, DIA_i, DIA_j
 1005       FORMAT (' STRIAN: R RHVi,j DIAi,j :',5E12.4)
          ENDIF
!
        ENDDO
      ENDIF
!
!     *** test output ***
!
      IF ( ITEST .GE. 5 .AND. TESTFL) THEN
         WRITE(PRINTF,2000) KCGRD(1), IRES, ISMAX
 2000    FORMAT (' STRIAN: KCGRD IRES ISMAX  :',4I4)
         WRITE(PRINTF,2001) GRAV, DEP, DEP_2, DEP_3
 2001    FORMAT (' STRIAN: G DEP DEP2 DEP3   :',4E12.4)
         WRITE(PRINTF,2002) DPI, PTRIAD(1), PTRIAD(2), URSELL(KCGRD(1))
 2002    FORMAT (' STRIAN: DPI P(1) P(2) URSELL  :',4E12.4)
         WRITE(PRINTF,2003) SMEBRK, HS, BIPH, ABS(SIN(BIPH))
 2003    FORMAT (' STRIAN: SMEBRK HS B SIN(B):',4E12.4)
      ENDIF

      RETURN
      END subroutine STRIAN                                               40.13
!
!********************************************************************
!
!     <<  Numerical Computations of the Nonlinear Energy Transfer
!           of Gravity Wave Spectra in Finite Water Depth  >>
!
!                  << developed by Noriaki Hashimoto >>
!
!        References:  N. Hashimoto et al. (1998)
!                     Komatsu and Masuda (2000) (in Japanese)
!
      SUBROUTINE RIAM_SLW(LMAX,N,N2,G,H,DQ,DQ2,DT,DT2,W,P,ACT,SNL,mint)   40.17

!     LMAX
!     N     : Number of directional bins                                  40.17
!     N2
!     G     : Gravitational acceleration                                  40.17
!     H     : Depth                                                       40.17
!     DQ
!     DQ2
!     DT    : size of the directional bins (Delta Theta)                  40.17
!     DT2
!     W     : discretised frequency array                                 40.17
!     P     : Density of water                                            40.17
!     ACT   : Action density                                              40.17
!     SNL   : Quadruplet source term                                      40.17
!     mint

!     IW4   : counter for the 4th frequency                               40.17
!     MMM
!     W4    : frequency of the fourth quadruplet component                40.17
!     AK4   : wavenumber of the fourth quadruplet component               40.17
!     DNA   : Coefficient in eq. 17 of Hashimoto (98)                     40.17
!             = 2 w4 k4 / Cg(k4)                                          40.17
!     IW3L
!     II
!     JJ
!     DI
!     DJ
!     CGK4  : group velocity for the fourth quadruplet component          40.17

      PARAMETER(MAX=1000000)                                              40.17

      REAL      :: W(LMAX), ACT(LMAX,N), SNL(LMAX,N)                      40.17
      DIMENSION SSS(MAX),II(3,MAX),JJ(3,MAX),DI(3,MAX),DJ(3,MAX)


!     Initialisation of the quadruplet source term                        40.17

      SNL(:,:) = 0.                                                       40.17
!
!     =================
      DO 130 IW4=1,LMAX
!     =================
!
      MMM=0
      W4=W(IW4)

!     WAVE converts nondimensional Sqr(w)d/g to nondimensional kd

      AK4=WAVE(W4**2*H/G)/H

!     CGCMP computes group velocity

      CALL CGCMP(G,AK4,H,CGK4)

!     Calculates the coefficient in equation (17) of Hashimoto (1998)

      DNA=2.*W4*AK4/CGK4

      IW3L=MAX0(1,NINT(IW4-ALOG(3.)/DQ))

!     -----------------------------------------------------
      CALL PRESET(LMAX,n,IW3L,IW4,N2,G,H,W4,AK4,W,DQ,DQ2,DT,DT2,          40.17
     &                                 SSS,II,JJ,DI,DJ,MMM,P)             40.17
!     -----------------------------------------------------
!
!
!     ==============
      DO 130 IT4=1,N
!     ==============
!
!     ==============
      DO 120 M=1,MMM  !    (W1 - W3 - W4 - W2)
!     ==============
!
      K1=II(1,M)+IW4
      K2=II(2,M)+IW4
      K3=II(3,M)+IW4
!
      if(k1.lt.1) go to 120                                               40.17
      IF(k2.gt.lmax) GOTO 120                                             40.17
!
      GG=SSS(M)*DNA
!
      M1P= JJ(1,M)+IT4
      M1N=-JJ(1,M)+IT4
      M2P= JJ(2,M)+IT4
      M2N=-JJ(2,M)+IT4
      M3P= JJ(3,M)+IT4
      M3N=-JJ(3,M)+IT4
!
      M1P=MOD(M1P,N)
      M1N=MOD(M1N,N)
      M2P=MOD(M2P,N)
      M2N=MOD(M2N,N)
      M3P=MOD(M3P,N)
      M3N=MOD(M3N,N)
!
      IF(M1P.LT.1) M1P=M1P+N
      IF(M1N.LT.1) M1N=M1N+N
      IF(M2P.LT.1) M2P=M2P+N
      IF(M2N.LT.1) M2N=M2N+N
      IF(M3P.LT.1) M3P=M3P+N
      IF(M3N.LT.1) M3N=M3N+N
!
      A4=ACT(IW4,IT4)
!
      IF(MINT.EQ.1) THEN
!
        K01=0
        K02=0
        K03=0
        M01=0
        M02=0
        M03=0
        IF(DI(1,M).LE.1.) K01=1
        IF(DI(2,M).LE.1.) K02=1
        IF(DI(3,M).LE.1.) K03=1
        IF(DJ(1,M).LE.1.) M01=-1
        IF(DJ(2,M).LE.1.) M02=-1
        IF(DJ(3,M).LE.1.) M03=-1
!
        K11=K1+K01
        K21=K2+K02
        K31=K3+K03
        IF(K11.GT.LMAX) K11=LMAX
        IF(K21.GT.LMAX) K21=LMAX
        IF(K31.GT.LMAX) K31=LMAX
!
        K111=K1+K01-1
        K211=K2+K02-1
        K311=K3+K03-1
        IF(K111.LT.1) K111=1
        IF(K211.LT.1) K211=1
        IF(K311.LT.1) K311=1
!
        M1P1=M1P+M01
        M2P1=M2P+M02
        M3P1=M3P+M03
        IF(M1P1.LT.1) M1P1=N
        IF(M2P1.LT.1) M2P1=N
        IF(M3P1.LT.1) M3P1=N
!
        M1N1=M1N-M01
        M2N1=M2N-M02
        M3N1=M3N-M03
        IF(M1N1.GT.N) M1N1=1
        IF(M2N1.GT.N) M2N1=1
        IF(M3N1.GT.N) M3N1=1
!
        M1P11=M1P+M01+1
        M2P11=M2P+M02+1
        M3P11=M3P+M03+1
        IF(M1P11.GT.N) M1P11=1
        IF(M2P11.GT.N) M2P11=1
        IF(M3P11.GT.N) M3P11=1
!
        M1N11=M1N-M01-1
        M2N11=M2N-M02-1
        M3N11=M3N-M03-1
        IF(M1N11.LT.1) M1N11=N
        IF(M2N11.LT.1) M2N11=N
        IF(M3N11.LT.1) M3N11=N
!
        A1P=(DI(1,M)*(DJ(1,M)*ACT(K11, M1P1)+ACT(K11, M1P11))
     &              + DJ(1,M)*ACT(K111,M1P1)+ACT(K111,M1P11))
     &              /((1.+di(1,m))*(1.+dj(1,m)))
!
        A1N=(DI(1,M)*(DJ(1,M)*ACT(K11, M1n1)+ACT(K11, M1n11))
     &              + DJ(1,M)*ACT(K111,M1n1)+ACT(K111,M1n11))
     &              /((1.+di(1,m))*(1.+dj(1,m)))
!
        A2P=(DI(2,M)*(DJ(2,M)*ACT(K21, M2P1)+ACT(K21, M2P11))
     &              + DJ(2,M)*ACT(K211,M2P1)+ACT(K211,M2P11))
     &              /((1.+di(2,m))*(1.+dj(2,m)))
!
        A2N=(DI(2,M)*(DJ(2,M)*ACT(K21, M2n1)+ACT(K21, M2n11))
     &              + DJ(2,M)*ACT(K211,M2n1)+ACT(K211,M2n11))
     &              /((1.+di(2,m))*(1.+dj(2,m)))
!
        A3P=(DI(3,M)*(DJ(3,M)*ACT(K31, M3P1)+ACT(K31, M3P11))
     &              + DJ(3,M)*ACT(K311,M3P1)+ACT(K311,M3P11))
     &              /((1.+di(3,m))*(1.+dj(3,m)))
!
        A3N=(DI(3,M)*(DJ(3,M)*ACT(K31, M3n1)+ACT(K31, M3n11))
     &              + DJ(3,M)*ACT(K311,M3n1)+ACT(K311,M3n11))
     &              /((1.+di(3,m))*(1.+dj(3,m)))
!
      ELSE
!
        IF(K1.LT.1) K1=1
        IF(K2.LT.1) K2=1
        IF(K3.LT.1) K3=1
!
        A1P=ACT(K1,M1P)
        A1N=ACT(K1,M1N)
        A2P=ACT(K2,M2P)
        A2N=ACT(K2,M2N)
        A3P=ACT(K3,M3P)
        A3N=ACT(K3,M3N)
!
      ENDIF
!
      W1P2P=A1P+A2P
      W1N2N=A1N+A2N
      S1P2P=A1P*A2P
      S1N2N=A1N*A2N
      W3P4=A3P+A4
      W3N4=A3N+A4
      S3P4=A3P*A4
      S3N4=A3N*A4
!
      XP=S1P2P*W3P4-S3P4*W1P2P
      XN=S1N2N*W3N4-S3N4*W1N2N
!
      SNL(K1,M1P)=SNL(K1,M1P)-XP*GG
      SNL(K1,M1N)=SNL(K1,M1N)-XN*GG
      SNL(K2,M2P)=SNL(K2,M2P)-XP*GG
      SNL(K2,M2N)=SNL(K2,M2N)-XN*GG
      SNL(K3,M3P)=SNL(K3,M3P)+XP*GG
      SNL(K3,M3N)=SNL(K3,M3N)+XN*GG
      SNL(IW4,IT4)=SNL(IW4,IT4)+(XP+XN)*GG
!
! ============
  120 CONTINUE
! ============
!
! ============
  130 CONTINUE
! ============
!
      RETURN
      END

      SUBROUTINE PRESET(LMAX,N,IW3L,IW4,N2,G,H,W4,AK4,W,DQ,DQ2,DT,DT2,    40.17
     &                                         SSS,II,JJ,DI,DJ,M,P)       40.17

!     T1A   : Angle between theta_1 and theta_a
!     T2A   : Angle between theta_2 and theta_a
!     T34   : Angle between theta_3 and theta_4
!     WA    : wa = w1 + w2 = w3 + w4 (resonance conditions)
!     AKA   : absolute value of ka = k1 + k2 = k3 + k4 (resonance cond.)
!     TA    : Angle theta_a representing vector ka

      PARAMETER (MAX=1000000)                                             40.17
      REAL      :: W(LMAX)                                                40.17
      DIMENSION SSS(MAX),II(3,MAX),JJ(3,MAX),DI(3,MAX),DJ(3,MAX)          40.17
!
!     ================
      DO 120 IT34=1,N2
!     ================
!
      T34=DT*(IT34-1)
      DT3=DT
      IF(IT34.EQ.1.OR.IT34.EQ.N2) DT3=DT2
!
!     ===================
      DO 110 IW3=IW3L,IW4
!     ===================
!
      W3=W(IW3)
      AK3=WAVE(W3**2*H/G)/H
!
      WA=W3+W4
      AKA=SQRT(AK3*AK3+AK4*AK4+2.*AK3*AK4*COS(T34))
      TA=ATAN2(AK3*SIN(T34),AK3*COS(T34)+AK4)
      R=SQRT(G*AKA*TANH(AKA*H/2.))/WA-1./SQRT(2.)
!
      TL=0.
      ACS=AKA/(2.*WAVE(WA**2*H/(4.*G))/H)
      ACS=SIGN(1.,ACS)*MIN(1.,ABS(ACS))
      IF(R.LT.0.) TL=ACOS(ACS)
      IT1S=NINT(TL/DT)+1
!
!     ===================
      DO 100 IT1A=IT1S,N2
!     ===================
!
      T1A=DT*(IT1A-1)
!
      DT1=DT
      IF(IT1A.EQ.1.OR.IT1A.EQ.N2) DT1=DT2
!
!     ------------------------------------------------------------
      CALL FINDW1W2(LMAX,IW3,W,G,H,W1,W2,W3,WA,AK1,AK2,AKA,T1A,T2A,IND)   40.17
!     ------------------------------------------------------------
      IF(IND.LT.0) GO TO 100
!
      IF(W1.LE.W3.AND.W3.LE.W4.AND.W4.LE.W2) THEN
!       ----------------------------------------------------
        CALL KERNEL(n,G,H,W1,W2,W3,W4,AK1,AK2,AK3,AK4,AKA,
     &       T1A,T2A,T34,TA,DQ,DT,DT1,DT3,SSS,II,JJ,DI,DJ,M,P)            40.17
!       ----------------------------------------------------
      ENDIF
!
! ============
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
! ============
!
      RETURN
      END

      SUBROUTINE FINDW1W2(LMAX,                                           40.17
     &                    IW3,W,G,H,W1,W2,W3,WA,AK1,AK2,AKA,T1A,T2A,IND)  40.17
      REAL  :: W(LMAX)                                                    40.17
!
      IND=0
      EPS=0.0005
!
      X1=0.005
      X2=W3
!
!     ---------------------------------------
  110 CALL FDW1(G,H,AKA,T1A,WA,X1,X2,X,EPS,M)
!     ---------------------------------------
!
      IF(M.EQ.0) GO TO 999
!
      W1=X
!
      AK1=WAVE(W1**2*H/G)/H
      AK2=SQRT(AKA*AKA+AK1*AK1-2.*AKA*AK1*COS(T1A))
      W2=SQRT(G*AK2*TANH(AK2*H))
      IF(W1.GT.W2) GO TO 999
      T2A=ATAN2(-AK1*SIN(T1A),AKA-AK1*COS(T1A))
      RETURN
!
  999 IND=-999
      RETURN
      END

      SUBROUTINE FDW1(G,H,AKA,T1A,WA,X1,X2,X,EPS,M)
!
      M=0
      F1=FUNCW1(G,H,X1,AKA,T1A,WA)
      F2=FUNCW1(G,H,X2,AKA,T1A,WA)
   20 IF(F1*F2) 18,19,21
   18 M=M+1
      X=X2-((X2-X1)/(F2-F1)*F2)
      IF((ABS(X-X2)/ABS(X))-EPS) 22,23,23
   23 IF((ABS(X-X1)/ABS(X))-EPS) 22,24,24
   22 RETURN
   24 F=FUNCW1(G,H,X,AKA,T1A,WA)
      FM=F*F1
      IF(FM) 25,22,26
   25 X2=X
      F2=F
      GO TO 20
   26 X1=X
      F1=F
      GO TO 20
   19 M=M+1
      IF(F1) 17,16,17
   16 X=X1
      GO TO 22
   17 X=X2
      GO TO 22
   21 M=0
      RETURN
      END

      FUNCTION FUNCW1(G,H,X,AKA,T1A,WA)
!
      AK1=WAVE(X**2*H/G)/H
      AK2=SQRT(AKA*AKA+AK1*AK1-2.*AKA*AK1*COS(T1A))
      FUNCW1=WA-SQRT(G*AK1*TANH(AK1*H))-SQRT(G*AK2*TANH(AK2*H))
!
      RETURN
      END

      SUBROUTINE KERNEL(n,G,H,W1,W2,W3,W4,AK1,AK2,AK3,AK4,AKA,
     &             T1A,T2A,T34,TA,DQ,DT,DT1,DT3,SSS,II,JJ,DI,DJ,M,P)      40.17
      PARAMETER (MAX=1000000)
      DIMENSION SSS(MAX),II(3,MAX),JJ(3,MAX),DI(3,MAX),DJ(3,MAX)
!
      PI=3.141592654
      PI2=2.*PI
!
!     ------------------------
      CALL CGCMP(G,AK1,H,CGK1)
      CALL CGCMP(G,AK2,H,CGK2)
      CALL CGCMP(G,AK3,H,CGK3)
!     ------------------------

!     Calculate denominator S (eq. 14, Hashimoto (1998))

      SS0=ABS(1.+CGK2/CGK1*(AK1-AKA*COS(T1A))/AK2)
      IF(SS0.LE.1.E-7) RETURN

!                        k1 k3 w3     G
!     Kernel function  -------------  - dth3 dth1 dOm    (eq. 17)
!                      Cg(k1) Cg(k3)  S

!     Where S defined in eq. 14 and G defined in eq. 5 (Hashimoto, 1998)

      CF=AK1*AK3*W3/(CGK1*CGK3)/SS0
     &   *9.*PI*G*G/(4.*P*P*W1*W2*W3*W4)*DT1*DT3*DQ

!     --------------------------------------
      CALL NONKDD(G,H,W3,W4,-W1,AK3,AK4,AK1,
     &                T34,0.,T1A+TA+PI,DDD1)
!     --------------------------------------
      GGG1=CF*DDD1*DDD1
!
!     --------------------------------------
      CALL NONKDD(G,H,W3,W4,-W1,AK3,AK4,AK1,
     &               T34,0.,-T1A+TA+PI,DDD2)
!     --------------------------------------
      GGG2=CF*DDD2*DDD2
!
      I1=1-intw(w1/w4,dq)
      I2=1-intw(w2/w4,dq)
      I3=1-intw(w3/w4,dq)
!
      DI1=weigw(w1/w4,dq)
      DI2=weigw(w2/w4,dq)
      DI3=0.
!
!     ----------------------------
!
      TT1=T1A+TA
      TT2=T2A+TA
      TT3=T34
!
      M=M+1
!
      SSS(M)=GGG1
      II(1,M)=I1
      II(2,M)=I2
      II(3,M)=I3
      DI(1,M)=DI1
      DI(2,M)=DI2
      DI(3,M)=DI3
      JJ(1,M)=intt(tt1,dt,pi2,n)-1
      JJ(2,M)=intt(tt2,dt,pi2,n)-1
      JJ(3,M)=intt(tt3,dt,pi2,n)-1
      DJ(1,M)=weigt(tt1,dt,pi2,n)
      DJ(2,M)=weigt(tt2,dt,pi2,n)
      DJ(3,M)=0.
!
!     ----------------------------
!
      TT1=-T1A+TA
      TT2=-T2A+TA
      TT3=T34
!
      M=M+1
!
      SSS(M)=GGG2
      II(1,M)=I1
      II(2,M)=I2
      II(3,M)=I3
      DI(1,M)=DI1
      DI(2,M)=DI2
      DI(3,M)=DI3
      JJ(1,M)=intt(tt1,dt,pi2,n)-1
      JJ(2,M)=intt(tt2,dt,pi2,n)-1
      JJ(3,M)=intt(tt3,dt,pi2,n)-1
      DJ(1,M)=weigt(tt1,dt,pi2,n)
      DJ(2,M)=weigt(tt2,dt,pi2,n)
      DJ(3,M)=0.
!
      RETURN
      END

      FUNCTION INTW(W,DQ)
      AA=ALOG(W)/DQ
      INTW=1+NINT(-AA)
      RETURN
      END

      FUNCTION WEIGW(W,DQ)
      AA=ALOG(W)/DQ
      L=1+INT(-AA)
      A=DQ*(1-L)-ALOG(W)
      IF(A.NE.0.) THEN
        B=ALOG(W)+DQ*L
        WEIGW=B/A
      ELSE
        WEIGW=1000.
      ENDIF
      RETURN
      END

      FUNCTION INTT(T,DT,PI2,MM)
      T=AMOD(T,PI2)
      IF(T.LT.0.) T=T+PI2
      INTT=NINT(T/DT)+1
      IF(INTT.GT.MM) INTT=INTT-MM
      IF(INTT.LT.1)  INTT=INTT+MM
      RETURN
      END

      FUNCTION WEIGT(T,DT,PI2,MM)
      T=AMOD(T,PI2)
      IF(T.LT.0.) T=T+PI2
      N=INT(T/DT)+1
      IF(N.GT.MM) N=N-MM
      IF(N.LT.1)  N=N+MM
      C=T-DT*(N-1)
      IF(C.NE.0.) THEN
        D=DT*N-T
        WEIGT=D/C
      ELSE
        WEIGT=1000.
      ENDIF
      RETURN
      END

      SUBROUTINE NONKDD(G,H,W1,W2,W3,AK1,AK2,AK3,T1,T2,T3,DDD)

!     See Herterich and Hasselmann (1980; p. 223, eq. B1)

      CALL NONKD(G,H,W1,W2,W3,AK1,AK2,AK3,T1,T2,T3,DD1)
      CALL NONKD(G,H,W2,W3,W1,AK2,AK3,AK1,T2,T3,T1,DD2)
      CALL NONKD(G,H,W3,W1,W2,AK3,AK1,AK2,T3,T1,T2,DD3)
!
      DDD=(DD1+DD2+DD3)/3.
!
      RETURN
      END

      SUBROUTINE NONKD(G,H,W1,W2,W3,AK1,AK2,AK3,T1,T2,T3,DDD)

!     GG    : Square of acceleration of gravity (=g**2)
!     G2    : Twice the square of acceleration of gravity (=2g**2)
!     WP23  : Sum of w2 and w3 (=w2+w3)
!     W123  : Sum of w1, w2 and w3 (=w1+w2+w3)
!     AKX1  : Size of the x-component of wavenumber k1
!     AKY1  : Size of the y-component of wavenumber k1
!     AKX2  : Size of the x-component of wavenumber k2
!     AKY2  : Size of the y-component of wavenumber k2
!     AKX3  : Size of the x-component of wavenumber k3
!     AKY3  : Size of the y-component of wavenumber k3
!     AK23  : Dot product of wavenumbers k2 and k3
!     W23   : frequency corresponding to wavenumber k2+k3


      DDD=0.
      GG=G*G
      G2=2.*G*G
!
      AKX1=AK1*COS(T1)
      AKY1=AK1*SIN(T1)
!
      AKX2=AK2*COS(T2)
      AKY2=AK2*SIN(T2)
!
      AKX3=AK3*COS(T3)
      AKY3=AK3*SIN(T3)
!
      AK23=SQRT((AKX2+AKX3)**2+(AKY2+AKY3)**2)
      W23=SQRT(G*AK23*TANH(AK23*H))
!
      WP23=W2+W3
      CF0=W23**2-WP23**2
      IF(CF0.EQ.0.) RETURN
!
      W123=W1+W2+W3
      AKXY23=AKX2*AKX3+AKY2*AKY3
      AKXY123=AKX1*(AKX2+AKX3)+AKY1*(AKY2+AKY3)
      CF1=0.
      IF(AK1*H.LE.10.) CF1=AK1*AK1/COSH(AK1*H)**2
!
      CF2=0.
      CF3=0.
      IF(AK3*H.LE.10.) CF2=W2*AK3*AK3/COSH(AK3*H)**2
      IF(AK2*H.LE.10.) CF3=W3*AK2*AK2/COSH(AK2*H)**2

!     See Herterich and Hasselmann (1980), first eq. after eq. B2
!     Equation is divided by I, to avoid complex numbers

      DD=WP23*(AK2*AK3*TANH(AK2*H)*TANH(AK3*H)-AKXY23)-(CF2+CF3)/2.

!     See Herterich and Hasselmann (1980), second eq. after eq. B2

      EE=(AKXY23-W2*W3*(W2**2+W3**2+W2*W3)/GG)/(2.*G)


      CF4=0.
      IF(AK23*H.LE.10.) CF4=W1*AK23**2/COSH(AK23*H)**2

!     Equation B1 of Herterich and Hasselmann (1980), where DD is
!     multiplied by I, to compensate for the earlier division

      DDD=-DD/CF0*(2.*W123*(W1**2*W23**2/GG-AKXY123)-CF4-WP23*CF1)
     &    +DD*W1*(W1*W1+W23*W23)/GG
     &    +EE*(W1**3*WP23/G-G*AKXY123-G*CF1)
     &    +W1*AKXY23*(W123*(W2**2+W3**2)+W2*W3*WP23)/G2
     &    -W1*W2*W2*AK3*AK3*(W123+W3)/G2
     &    -W1*W3*W3*AK2*AK2*(W123+W2)/G2
!
      RETURN
      END

      SUBROUTINE CGCMP(G,AK,H,CG)

!     Calculates group velocity Cg based on depth and wavenumber
!     Includes a deep water limit for kd > 10.

!     G     : Gravitational acceleration
!     AK    : Wave number
!     H     : Depth
!     CG    : Group velocity

!     AKH   : Depth x Wave number (kd)
!     RGK   : Square root of gk
!     SECH2 : The square of the secant Hyperbolic of kd
!             = 1 - TANH(kd)**2

!     Calculation of group velocity Cg

      AKH=AK*H
      IF(AKH.LE.10.) THEN                                                 40.33

!     Shallow water:

!     Cg = (1/2) (1 + (2 kd) / Sinh (2 kd)) (w / k)
!        = (1/2) (1 + (kd SECH2) / Tanh (kd)) (w / k)
!        = (1/2) (1 + (kd SECH2) / Tanh (kd)) (Sqrt (gk Tanh(kd)) / k)
!        = (1/2) (Sqrt (gk) (Sqrt(Tanh (kd)) + kd SECH2 / Sqrt (Tanh (kd)) / k
!        = (1/2) (Sqrt (k) g (Tanh (kd) + kd SECH2) / (k Sqrt (g Tanh (kd)))
!        = (g Tanh(kd) + gkd SECH2) / (2 Sqrt(gk Tanh(kd))


        SECH2=1.-TANH(AKH)**2
        CG=(G*AKH*SECH2+G*TANH(AKH))/(2.*SQRT(G*AK*TANH(AKH)))
      ELSE                                                                40.17

!     Deep water:

!     Cg = w / (2 k)
!        = g / (2 Sqrt (gk))

        RGK=SQRT(G*AK)
        CG=G/(2.*RGK)
      ENDIF                                                               40.17
!
      RETURN

      END SUBROUTINE CGCMP

      SUBROUTINE DCGCMP(G,AK,H,DCG)
!
      AKH=AK*H
      IF(AKH.GT.10.) GO TO 100
      SECH2=1-TANH(AKH)**2
      SECH3=SECH2*SQRT(SECH2)
      DCG=G*H*SECH3*(COSH(AKH)-AKH*SINH(AKH))/SQRT(G*AK*TANH(AKH))
     &   -(G*AKH*SECH2+G*TANH(AKH))**2/(4.*SQRT(G*AK*TANH(AKH))**3)
      RETURN
!
  100 RGK=SQRT(G*AK)
      DCG=-G*G/(4.*RGK**3)
!
      RETURN
      END

      REAL FUNCTION WAVE(D)

!     Transforms nondimensional Sqr(w)d/g into nondimensional kd using the
!     dispersion relation and an iterative method

      IF(D-10.) 2,2,1
    1 XX=D
      GO TO 6
    2 IF(D-1.) 3,4,4
    3 X=SQRT(D)
      GO TO 5
    4 X=D
    5 COTHX=1.0/TANH(X)
      XX=X-(X-D*COTHX)/(1.+D*(COTHX**2-1.))
      E=1.-XX/X
      X=XX
      IF(ABS(E)-0.0005) 6,5,5
    6 WAVE=XX
      RETURN
      END
