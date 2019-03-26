!                 SWAN COMMON BLOCKS, file 1 of 4
!
!     Contents of this file
!
!     COMMON / SWNAME / *** names and other character data
!     COMMON / SWUITV / *** information for output
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
!     30.72 IJsbrand G. Haagsma
!     32.01 Roeland Ris & Cor van der Schelde
!     40.00 Nico Booij
!
!  1. Updates
!
!     30.72, Sept 97: Created from the old swancomm.ftn file
!     32.01, Jan. 98: NMOVAR increased by 1 for wave-setup
!     40.00, Apr. 98: variable OVKEYW added for all quantities
!                     NMOVAR increased
!                     quantities TSEC and TIME added
!                     array OUTPAR added
!                     char. string CHTIME added
!     40.03, May  00: output variables URSELL, ASTD added
!
!  2. Purpose
!
!     Common variables used by the subroutines in SWAN.
!     This file needs to be included at the proper places with the
!     INCLUDE 'file.inc' compiler directive
!
!  3. Method
!
!     Defines common blocks, with the needed parameter variables.
!     All variables should be declared explicitely
!
!  4. Argument variables
!
!     --
!
!  5. Parameter variables
!
!        [???] default value
!
! NMOVAR [ 48] maximum number of output variables                         40.03
! MOUTPA [  5] number of output parameters
!
      INTEGER             NMOVAR, MOUTPA                                  40.00
      PARAMETER           (NMOVAR = 48, MOUTPA=5)                         40.03
!
!  6. Local variables
!
!     --
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
!     --
!
! 10. Error messages
!
!     --
!
!-------------------------------------------------------------------------------
!
!     COMMON / SWNAME /  *** names and other character data ***
!
!     Common variables
!
!        [???????] default value, [CALCULAT] is calculated
!
! CHTIME [ '    '] character string representation of date-time of computation 40.00
! FNEST [CALCULAT] BCF (=BLANK), name of nest file
! OVKEYW(NMOVAR)   keyword identifying output quantity in a SWAN command  40.00
!  ( 1) [    'XP']
!  ( 2) [    'YP']
!  ( 3) [  'DIST']
!  ( 4) [   'DEP']
!  ( 5) [   'VEL']
!  ( 6) [  'UBOT']
!  ( 7) [  'DISS']
!  ( 8) [    'QB']
!  ( 9) [   'LEA']
!  (10) [    'HS']
!  (11) [  'TM01']
!  (12) [   'RTP']
!  (13) [   'DIR']
!  (14) [   'PDI']
!  (15) [   'TDI']
!  (16) [  'DSPR']
!  (17) [  'WLEN']
!  (18) [  'STEE']
!  (19) [   'TRA']
!  (20) [   'FOR']
!  (21) ['      ']
!  (22) ['      ']
!  (23) ['      ']
!  (24) [    'XC']
!  (25) [    'YC']
!  (26) [  'WIND']
!  (27) [   'FRC']
!  (28) [ 'RTM01']
!  (29) ['      ']
!  (30) [   'DHS']
!  (31) ['DRTM01']
!  (32) [  'TM02']
!  (33) [  'FSPR']
!  (34) [  'URMS']
!  (35) [  'UFRI']
!  (36) [  'ZLEN']
!  (37) [  'TAUW']
!  (38) [ 'CDRAG']
!  (39) [ 'SETUP']                                                        40.00
!  (40) [  'TIME']                                                        40.00
!  (41) [  'TSEC']                                                        40.00
! OVSNAM(NMOVAR)   short name of output quantity
!  ( 1) [    'Xp']
!  ( 2) [    'Yp']
!  ( 3) [  'Dist']
!  ( 4) [ 'Depth']
!  ( 5) [   'Vel']
!  ( 6) [  'Ubot']
!  ( 7) ['Dissip']
!  ( 8) [    'Qb']
!  ( 9) [  'Leak']
!  (10) [    'Hs']
!  (11) [  'Tm01']
!  (12) [ 'Tpeak']
!  (13) [   'Dir']
!  (14) [ 'PkDir']
!  (15) [  'TDir']
!  (16) [  'Dspr']
!  (17) [  'Wlen']
!  (18) ['Steepn']
!  (19) ['Transp']
!  (20) ['WForce']
!  (21) ['AcDens']
!  (22) ['EnDens']
!  (23) [   'Aux']
!  (24) [    'Xc']
!  (25) [    'Yc']
!  (26) [ 'Windv']
!  (27) ['FrCoef']
!  (28) [ 'RTm01']
!  (29) ['EnDens']
!  (30) [   'dHs']
!  (31) [   'dTm']
!  (32) [  'Tm02']
!  (33) [  'FSpr']
!  (34) [  'Urms']
!  (35) [ 'Ufric']
!  (36) [  'Zlen']
!  (37) [  'TauW']
!  (38) [ 'Cdrag']
!  (39) [ 'Setup']
!  (40) [  'Time']                                                        40.00
!  (41) [  'Tsec']                                                        40.00
! OVLNAM(NMOVAR)   long name of output quantity
!  ( 1) [                        'X user coordinate']
!  ( 2) [                        'Y user coordinate']
!  ( 3) [              'distance along output curve']
!  ( 4) [                                    'Depth']
!  ( 5) [                         'Current velocity']
!  ( 6) [           'Orbital velocity at the bottom']
!  ( 7) [                       'Energy dissipation']
!  ( 8) [                  'Fraction breaking waves']
!  ( 9) [     'Energy leak over spectral boundaries']
!  (10) [                  'Significant wave height']
!  (11) [             'Average absolute wave period']
!  (12) [                              'Peak period']
!  (13) [                   'Average wave direction']
!  (14) [    'direction of the peak of the spectrum']
!  (15) [        'direction of the energy transport']
!  (16) [                    'directional spreading']
!  (17) [                      'Average wave length']
!  (18) [                           'Wave steepness']
!  (19) [                    'Wave energy transport']
!  (20) [       'Wave driven force per unit surface']
!  (21) [                  'spectral action density']
!  (22) [                  'spectral energy density']
!  (23) [                       'auxiliary variable']
!  (24) [          'X computational grid coordinate']
!  (25) [          'Y computational grid coordinate']
!  (26) [    'Wind velocity at 10 m above sea level']
!  (27) [              'Bottom friction coefficient']
!  (28) [             'Average relative wave period']
!  (29) [ 'energy density integrated over direction']
!  (30) [      'difference in Hs between iterations']
!  (31) [      'difference in Tm between iterations']
!  (32) [                     'Zero-crossing period']
!  (33) [         'Frequency spectral width (Kappa)']
!  (34) [    'RMS of orbital velocity at the bottom']
!  (35) [                        'Friction velocity']
!  (40) [                                'Date-time']                     40.00
!  (41) [      'Time in seconds from reference time']                     40.00
!  (36) ['Zero velocity thickness of boundary layer']
!  (37) [                                     '    ']
!  (38) [                         'Drag coefficient']
!  (39) [                       'Setup due to waves']
! OVUNIT(NMOVAR)   unit of of output quantity
!  ( 1) [CALCULAT] =UL
!  ( 2) [CALCULAT] =UL
!  ( 3) [CALCULAT] =UL
!  ( 4) [CALCULAT] =UH
!  ( 5) [CALCULAT] =UV
!  ( 6) [CALCULAT] =UV
!  ( 7) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  ( 8) [     ' ']
!  ( 9) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (10) [CALCULAT] =UH
!  (11) [CALCULAT] =UT
!  (12) [CALCULAT] =UT
!  (13) [CALCULAT] =UDI
!  (14) [CALCULAT] =UDI
!  (15) [CALCULAT] =UDI
!  (16) [CALCULAT] =UDI
!  (17) [CALCULAT] =UL
!  (18) [     ' ']
!  (19) [   'm2s'] changed to 'W/m' (=UP), if INRHOG=1
!  (20) [CALCULAT] =UF
!  (21) [   'm2s'] changed to 'Js/m2', if INRHOG=1
!  (22) [    'm2'] changed to 'J/m2', if INRHOG=1
!  (23) [     ' ']
!  (24) [     ' ']
!  (25) [     ' ']
!  (26) [CALCULAT] =UV
!  (27) [     ' ']
!  (28) [CALCULAT] =UT
!  (29) [    'm2'] changed to 'J/m2', if INRHOG=1
!  (30) [CALCULAT] =UH
!  (31) [CALCULAT] =UT
!  (32) [CALCULAT] =UT
!  (33) [     ' ']
!  (34) [CALCULAT] =UV
!  (35) [CALCULAT] =UV
!  (36) [CALCULAT] =UL
!  (37) [     ' ']
!  (38) [     ' ']
!  (39) [     'm']
!  (40) [     ' ']                                                        40.00
!  (41) [     's']                                                        40.00
! SNAME [        ] name of output point set
! UAP   [  'W/m2'] unit of dissipation
! UD    [        ] not used
! UDI   [  'degr'] unit of direction
! UDL   [  'm2/s'] unit of dissipation
! UET   [  'm3/s'] unit of energy transport, and wave force
! UF    [  'N/m2'] unit of pressure or shear stress (force per area)
! UH    [     'm'] unit of vertical length
! UL    [     'm'] unit of horizontal length
! UP    [   'W/m'] unit of energy flux density
! UST   [ 'm2/s2'] not used
! UT    [   'sec'] unit of time (change to s ???)
! UV    [   'm/s'] unit of velocity
!
      CHARACTER*20        CHTIME                                          40.00
      CHARACTER*36        FBCL,        FBCR,        FNEST
      CHARACTER*8         OVKEYW(NMOVAR)                                  40.00
      CHARACTER*40        OVLNAM(NMOVAR)
      CHARACTER*6         OVSNAM(NMOVAR)
      CHARACTER*16        OVUNIT(NMOVAR)
      CHARACTER*8         SNAME
      CHARACTER*6         UAP,         UD,          UDI,         UDL
      CHARACTER*6         UET,         UF,          UH,          UL
      CHARACTER*6         UP,          UST,         UT,          UV
!
      COMMON / SWNAME /   FNEST,       FBCL,        FBCR,        SNAME
      COMMON / SWNAME /   OVKEYW                                          40.00
      COMMON / SWNAME /   OVSNAM,      OVLNAM,      OVUNIT,      UH
      COMMON / SWNAME /   UV,          UT,          UL,          UET
      COMMON / SWNAME /   UDI,         UST,         UF,          UP
      COMMON / SWNAME /   UAP,         UDL,         UD,          CHTIME   40.00
!
!-------------------------------------------------------------------------------
!
!     COMMON / SWUITV /  *** information for output ***
!
!     Common variables
!
!        [???] default value
!
! ERRPTS [   ]   unit ref. number of file containing coord. of "problem points"
!                =16, unit reference number of the file
! INRHOG [  0]   indicates the choice for output based on "variance" or "true energy"
!                =0, output based on variance
!                =1, output based on true energy
! IUBOTR [  0]   set to 1, when IVTYPE=6 or 18
! OVSVTY(NMOVAR) type of the output variable
!                =1, scalar
!                =2, angle
!                =3, vector
!                =4, tensor
!                =5, fully spectral quantity
!                =6, directional spectral quantity
!  ( 1)  [  1]   Xp
!  ( 2)  [  1]   Yp
!  ( 3)  [  1]   Dist
!  ( 4)  [  1]   Depth
!  ( 5)  [  3]   Vel
!  ( 6)  [  1]   Ubot
!  ( 7)  [  1]   Dissip
!  ( 8)  [  1]   Qb
!  ( 9)  [  1]   Leak
!  (10)  [  1]   Hs
!  (11)  [  1]   Tm01
!  (12)  [  1]   Tpeak
!  (13)  [  2]   Dir
!  (14)  [  2]   PkDir
!  (15)  [  2]   TDir
!  (16)  [  1]   Dspr
!  (17)  [  1]   Wlen
!  (18)  [  1]   Steepn
!  (19)  [  3]   Transp
!  (20)  [  3]   WForce
!  (21)  [  5]   AcDens
!  (22)  [  5]   EnDens
!  (23)  [  1]   Aux
!  (24)  [  1]   Xc
!  (25)  [  1]   Yc
!  (26)  [  3]   Windv
!  (27)  [  1]   FrCoef
!  (28)  [  1]   RTm01
!  (29)  [  5]   EnDens
!  (30)  [  1]   dHs
!  (31)  [  1]   dTm
!  (32)  [  1]   Tm02
!  (33)  [  1]   FSpr
!  (34)  [  1]   Urms
!  (35)  [  1]   Ufric
!  (36)  [  1]   Zlen
!  (37)  [  1]   TauW
!  (38)  [  1]   Cdrag
!  (39)  [  1]   Setup
! SPCPOW [  1]   power in expression for computation of average frequency
!                =power; set by command SET ... [power] ... (not documented)
!
      INTEGER             ERRPTS
      INTEGER             INRHOG,      IUBOTR
      INTEGER             OVSVTY(NMOVAR),           SPCPOW
!
!
!        [???????] default value, [CALCULAT] is calculated
!
! AKPOWR [       ] power in expression for computation of average wave number
!                  =kpower; set by command SET ... [kpower] ... (not documented)
! ALCQ   [       ] angle between x-axes of computational grid and output frame
! ALPQ   [       ] angle between x-axes of user coord. system and output frame
! COSCQ  [       ] cos of ALCQ
! COSPQ  [       ] cos of ALPQ
! DXK    [       ] mesh size of output frame
! DYK    [       ] mesh size of output frame
! OVEXCV(NMOVAR)   exception value for output quantity
!  ( 1)  [ -1.E10] Xp
!  ( 2)  [ -1.E10] Yp
!  ( 3)  [  -999.] Dist
!  ( 4)  [  -1.E5] Depth
!  ( 5)  [     0.] Vel
!  ( 6)  [   -10.] Ubot
!  ( 7)  [   -10.] Dissip
!  ( 8)  [    -1.] Qb
!  ( 9)  [   -10.] Leak
!  (10)  [   -10.] Hs
!  (11)  [   -10.] Tm01
!  (12)  [   -10.] Tpeak
!  (13)  [  -999.] Dir
!  (14)  [  -999.] PkDir
!  (15)  [  -999.] TDir
!  (16)  [   -10.] Dspr
!  (17)  [   -10.] Wlen
!  (18)  [    -1.] Steepn
!  (19)  [     0.] Transp
!  (20)  [     0.] WForce
!  (21)  [  -100.] AcDens
!  (22)  [  -100.] EnDens
!  (23)  [ -1.E10] Aux
!  (24)  [ -1.E10] Xc
!  (25)  [ -1.E10] Yc
!  (26)  [     0.] Windv
!  (27)  [    -1.] FrCoef
!  (28)  [   -10.] RTm01
!  (29)  [  -100.] EnDens
!  (30)  [   -10.] dHs
!  (31)  [   -10.] dTm
!  (32)  [   -10.] Tm02
!  (33)  [    -1.] FSpr
!  (34)  [   -10.] Urms
!  (35)  [   -10.] Ufric
!  (36)  [    -1.] Zlen
!  (37)  [   -10.] TauW
!  (38)  [    -1.] Cdrag
!  (39)  [    -9.] Setup
! OVLEXP(NMOVAR)   lower expected limit of output quantity
!  ( 1)  [ -1.E10] Xp
!  ( 2)  [ -1.E10] Yp
!  ( 3)  [     0.] Dist
!  ( 4)  [  -100.] Depth
!  ( 5)  [    -2.] Vel
!  ( 6)  [     0.] Ubot
!  ( 7)  [     0.] Dissip
!  ( 8)  [     0.] Qb
!  ( 9)  [     0.] Leak
!  (10)  [     0.] Hs
!  (11)  [     0.] Tm01
!  (12)  [     0.] Tpeak
!  (13)  [     0.] Dir
!  (14)  [     0.] PkDir
!  (15)  [     0.] TDir
!  (16)  [     0.] Dspr
!  (17)  [     0.] Wlen
!  (18)  [     0.] Steepn
!  (19)  [   -10.] Transp
!  (20)  [   -10.] WForce
!  (21)  [     0.] AcDens
!  (22)  [     0.] EnDens
!  (23)  [ -1.E10] Aux
!  (24)  [  -100.] Xc
!  (25)  [  -100.] Yc
!  (26)  [   -50.] Windv
!  (27)  [     0.] FrCoef
!  (28)  [     0.] RTm01
!  (29)  [     0.] EnDens
!  (30)  [     0.] dHs
!  (31)  [     0.] dTm
!  (32)  [     0.] Tm02
!  (33)  [     0.] FSpr
!  (34)  [     0.] Urms
!  (35)  [     0.] Ufric
!  (36)  [     0.] Zlen
!  (37)  [     0.] TauW
!  (38)  [     0.] Cdrag
!  (39)  [    -1.] Setup
! OVLLIM(NMOVAR)   lower limit of validity of output quantity
!  ( 1)  [ -1.E10] Xp
!  ( 2)  [ -1.E10] Yp
!  ( 3)  [     0.] Dist
!  ( 4)  [  -1.E4] Depth
!  ( 5)  [  -100.] Vel
!  ( 6)  [     0.] Ubot
!  ( 7)  [     0.] Dissip
!  ( 8)  [     0.] Qb
!  ( 9)  [     0.] Leak
!  (10)  [     0.] Hs
!  (11)  [     0.] Tm01
!  (12)  [     0.] Tpeak
!  (13)  [     0.] Dir
!  (14)  [     0.] PkDir
!  (15)  [     0.] TDir
!  (16)  [     0.] Dspr
!  (17)  [     0.] Wlen
!  (18)  [     0.] Steepn
!  (19)  [  -100.] Transp
!  (20)  [  -1.E5] WForce
!  (21)  [     0.] AcDens
!  (22)  [     0.] EnDens
!  (23)  [ -1.E10] Aux
!  (24)  [ -1000.] Xc
!  (25)  [ -1000.] Yc
!  (26)  [  -100.] Windv
!  (27)  [     0.] FrCoef
!  (28)  [     0.] RTm01
!  (29)  [     0.] EnDens
!  (30)  [     0.] dHs
!  (31)  [     0.] dTm
!  (32)  [     0.] Tm02
!  (33)  [     0.] FSpr
!  (34)  [     0.] Urms
!  (35)  [     0.] Ufric
!  (36)  [     0.] Zlen
!  (37)  [     0.] TauW
!  (38)  [     0.] Cdrag
!  (39)  [    -1.] Setup
! OVHEXP(NMOVAR)   upper expected limit of output quantity
!  ( 1)  [  1.E10] Xp
!  ( 2)  [  1.E10] Yp
!  ( 3)  [  1.E10] Dist
!  ( 4)  [   100.] Depth
!  ( 5)  [     2.] Vel
!  ( 6)  [     1.] Ubot
!  ( 7)  [   100.] Dissip
!  ( 8)  [   100.] Qb
!  ( 9)  [   100.] Leak
!  (10)  [    10.] Hs
!  (11)  [   100.] Tm01
!  (12)  [   100.] Tpeak
!  (13)  [   360.] Dir
!  (14)  [   360.] PkDir
!  (15)  [   360.] TDir
!  (16)  [    60.] Dspr
!  (17)  [   200.] Wlen
!  (18)  [    0.1] Steepn
!  (19)  [    10.] Transp
!  (20)  [    10.] WForce
!  (21)  [   100.] AcDens
!  (22)  [   100.] EnDens
!  (23)  [  1.E10] Aux
!  (24)  [   100.] Xc
!  (25)  [   100.] Yc
!  (26)  [    50.] Windv
!  (27)  [     1.] FrCoef
!  (28)  [   100.] RTm01
!  (29)  [   100.] EnDens
!  (30)  [     1.] dHs
!  (31)  [     2.] dTm
!  (32)  [   100.] Tm02
!  (33)  [     1.] FSpr
!  (34)  [     1.] Urms
!  (35)  [     1.] Ufric
!  (36)  [     1.] Zlen
!  (37)  [     1.] TauW
!  (38)  [     1.] Cdrag
!  (39)  [     1.] Setup
! OVULIM(NMOVAR)   upper limit of validity
!  ( 1)  [  1.E10] Xp
!  ( 2)  [  1.E10] Yp
!  ( 3)  [  1.E10] Dist
!  ( 4)  [   1.E4] Depth
!  ( 5)  [   100.] Vel
!  ( 6)  [    10.] Ubot
!  ( 7)  [  1000.] Dissip
!  ( 8)  [     1.] Qb
!  ( 9)  [  1000.] Leak
!  (10)  [   100.] Hs
!  (11)  [  1000.] Tm01
!  (12)  [  1000.] Tpeak
!  (13)  [   360.] Dir
!  (14)  [   360.] PkDir
!  (15)  [   360.] TDir
!  (16)  [   360.] Dspr
!  (17)  [  1000.] Wlen
!  (18)  [     1.] Steepn
!  (19)  [   100.] Transp
!  (20)  [   1.E5] WForce
!  (21)  [  1000.] AcDens
!  (22)  [  1000.] EnDens
!  (23)  [  1.E10] Aux
!  (24)  [  1000.] Xc
!  (25)  [  1000.] Yc
!  (26)  [   100.] Windv
!  (27)  [     1.] FrCoef
!  (28)  [  1000.] RTm01
!  (29)  [  1000.] EnDens
!  (30)  [   100.] dHs
!  (31)  [   100.] dTm
!  (32)  [  1000.] Tm02
!  (33)  [     1.] FSpr
!  (34)  [    10.] Urms
!  (35)  [    10.] Ufric
!  (36)  [     1.] Zlen
!  (37)  [    10.] TauW
!  (38)  [     1.] Cdrag
!  (39)  [     1.] Setup
! SINCQ  [       ] sin of ALCQ
! SINPQ  [       ] sin of ALPQ
! XQLEN  [       ] length of x-side of output frame
! XQP    [       ] x-coordinate (user coord.) of origin of output frame
! YQLEN  [       ] length of y-side of output frame
! YQP    [       ] y-coordinate (user coord.) of origin of output frame
!
      REAL                AKPOWR,      ALCQ,        ALPQ,        COSCQ
      REAL                COSPQ,       DXK,         DYK
      REAL                OVEXCV(NMOVAR),           OVLEXP(NMOVAR)
      REAL                OVLLIM(NMOVAR),           OVHEXP(NMOVAR)
      REAL                OVULIM(NMOVAR),           SINCQ,       SINPQ
      REAL                XPQ,         XQLEN,       XQP,         YPQ
      REAL                YQLEN,       YQP,         OUTPAR(MOUTPA)        40.00
!
      COMMON / SWUITV /   ALCQ,        COSCQ,       SINCQ,       IUBOTR
      COMMON / SWUITV /   INRHOG,      ERRPTS,      DXK,         DYK
      COMMON / SWUITV /   ALPQ,        COSPQ,       SINPQ,       XQP
      COMMON / SWUITV /   YQP,         XQLEN,       YQLEN,       OVSVTY
      COMMON / SWUITV /   OVLLIM,      OVULIM,      OVLEXP,      OVHEXP
      COMMON / SWUITV /   OVEXCV,      SPCPOW,      AKPOWR
      COMMON / SWUITV /   XPQ,         YPQ,         OUTPAR                40.00
!
