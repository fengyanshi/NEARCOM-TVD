!                 SWAN COMMON BLOCKS, file 4 of 4
!
!     Contents of this file
!
!     COMMON / SWTEST / *** information for test output
!     COMMON / TP_SWTEST /  *** thread private components of SWTEST       40.22
!     COMMON / SWPROP / *** 3rd order propagation and shperical coordinates
!     COMMON / TP_SWPROP /  *** thread private components of SWPROP       40.22
!     COMMON / TP_WCAP /  *** thread private whitecapping data            40.22
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
!     30.72: IJsbrand G. Haagsma
!     40.00: Nico Booij
!     40.22: John Cazes and Tim Campbell                                  40.22
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Sept 97: Created from the old swancomm.ftn file
!     30.72, Mar. 98: Introduced UNDFLW in / SWTEST /
!     40.00, Apr. 98: IFS1D and IFS2D added (output of source terms in test
!                     points
!     40.00, Nov. 98: IFPAR added (output of wave parameters in test points)
!     40.22, Sep. 01: Moved IPTST and TESTFL into thread private common   40.22
!                     block TP_SWTEST for OpenMP                          40.22
!     40.22, Sep. 01: Moved COSLAT and PROPSL into thread private common  40.22
!                     block TP_SWPROP for OpenMP                          40.22
!     40.23, Nov. 02: parameter PROPFL added meant for flux-limiting      40.23
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
!     --
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
!     COMMON / SWTEST /  *** information for test output ***
!
!     Common variables
!
!        [???] default value
!
! ICOTES [  0] minimum value for ITEST,
!              =30; for command COTES, not documented in the manual
!              =cotes; for command COTES [cotes], not documented
! INTES  [  0] testing parameter,
!              =30; for command INTE, not documented in the manual
!              =intes; for command INTE [intes], not documented
! IOUTES [  0] minimum value for ITEST,
!              =30; for command OUTE, not documented in the manual
!              =itest; for command OUTE [itest], not documented
! IPTST  [   ] sequence number of a test point
! IFPAR  [  0] unit ref. number for output of parameters in test points   40.00
! IFS1D  [  0] unit ref. number for output of 1D spectra of source terms  40.00
!              if used, value is made non-zero by subr FOR
! IFS2D  [  0] unit ref. number for output of 2D spectra of source terms  40.00
!              if used, value is made non-zero by subr FOR
! LXDMP  [ -1] grid counter for a test point in the x-direction,
!              =ix; set by command TEST ... POI <[ix] [iy]>
! LYDMP  [  0] grid counter for a test point in the y-direction,
!              =iy; set by command TEST ... POI <[ix] [iy]>
! MAXMES [200] not used
!              =maxmes; set by command SET [maxmes]
! NEGMES [  0] not used
! NPTST  [  0] number of test points; set by command TEST                 40.00
! NPTSTA [  1] number of test points, equal to MAX(1,NPTST)               40.00
!
      INTEGER             ICOTES,      INTES,       IOUTES,      IPTST
      INTEGER             IFPAR,       IFS1D,       IFS2D                 40.00
      INTEGER             LXDMP,       LYDMP,       MAXMES,      NEGMES
      INTEGER             NPTST,       NPTSTA
!
! TESTFL [FAL] test output must/must not be made, mainly for testpoints
!
      LOGICAL             TESTFL
!
! UNDFLW [      1.E-15] Small number to prevent underflows                30.72
!
      REAL                UNDFLW                                          30.72
!
      COMMON / SWTEST /   LXDMP,       LYDMP,       NEGMES,      MAXMES
      COMMON / SWTEST /   NPTST,       NPTSTA                             40.22
      COMMON / SWTEST /   INTES,       ICOTES,      IOUTES,      UNDFLW   30.72
      COMMON / SWTEST /   IFPAR,       IFS1D,       IFS2D                 40.00

      COMMON / TP_SWTEST /  IPTST,       TESTFL                           40.22
!$OMP THREADPRIVATE(/TP_SWTEST/)                                          40.22

!
!  0. Authors
!
!     33.08: Erick Rogers, Jim Kaihatu
!     33.09: Nico Booij
!
!  1. Updates
!
!     33.08, June 98: new common area
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
!     --
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
!     COMMON / SWPROP /  *** higher order propagation and spherical coord ***  33.08
!
!     Common variables
!
!        [???] default value
!
!     PROJ_METHOD  [0]  projection method; 0=(quasi-)Cartesian,
!                       1=uniform Mercator (only spherical coordinates)
!     PROPFL  [cal]  indicates whether flux-limiting in spectral space is used
!     PROPSC  [cal]  indicates which numerical scheme is to be used for spatial propagation.
!                    1=first order (BSBT), 2=SORDUP, 3=3rd order (S&L)
!     PROPSS  [  2]  indicates which numerical scheme is to be used in stationary computations
!                    1=first order (BSBT), 2=SORDUP
!     PROPSN  [  3]  indicates which numerical scheme is to be used in nonstationary computations
!                    1=first order (BSBT), 3=3rd order (S&L)
!     PROPSL  [cal]  indicates which numerical scheme is used locally.
!     KSPHER  [  0]  indicates whether spherical coordinates are used, and which projection method
!                    0=Cartesian coordinates, >0=spherical coordinates
!     KREPTX  [  0]  if >0, the domain repeats itself in x-direction (primarily intended for
!                    propagation around the globe)
!
      INTEGER             PROPSC,    PROPSL,    KSPHER,      KREPTX
      INTEGER             PROPSS,    PROPSN,    PROJ_METHOD
      INTEGER             PROPFL                                          40.23
!
!     WAVAGE  [ 0.]  indicates "wave age" parameter (used in countering Garden Sprinkler Effect
!                    in subroutine SANDL)
!     REARTH  [6E6]  radius of the earth
!     LENDEG  [1E5]  length of a degree NS
!     COSLAT  [ 1.]  cos of latitude; =1 for Cartesian coordinates
!
      REAL                WAVAGE,    REARTH,    LENDEG
      REAL                COSLAT(10)                                         33.09
!
!
      COMMON / SWPROP /   PROPSC,    PROPSS,    PROPSN                    40.22
      COMMON / SWPROP /   PROPFL                                          40.23
      COMMON / SWPROP /   WAVAGE                                             33.08
      COMMON / SWPROP /   KSPHER,    REARTH,    LENDEG,      KREPTX          33.08
      COMMON / SWPROP /   PROJ_METHOD                                     40.22

      COMMON / TP_SWPROP /   COSLAT,PROPSL                                40.22
!$OMP THREADPRIVATE(/TP_SWPROP/)                                          40.22

!
!
!  0. Authors
!
!     40.22: John Cazes and Tim Campbell                                  40.22
!
!  1. Updates
!
!     Sep. 01: New include file for variables moved from module           40.22
!              swmod to thread private common block TP_WCAP for           40.22
!              OpenMP.  OpenMP support of modules is unclear.             40.22
!
!  2. Purpose
!
!     Common variables used by the subroutines SINTGRL and SWCAP          40.22
!     in SWAN that must be thread private for OpenMP.                     40.22
!
!  3. Method
!
!     Thread private common block                                         40.22
!
!  4. Modules used
!
!     ---
!
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     ---
!
!  7. Local variables
!
!     ACTOT  : Total action density per gridpoint
!     EDRKTOT: Zeroth moment of energy / SQRT(wavenumber)
!     EKTOT  : Zeroth moment of energy * wavenumber
!     ETOT1  : First moment of the energy density
!     ETOT2  : Second moment of the energy density
!     ETOT4  : Fourth moment of the energy density
!     KM_WAM : Mean average wavenumber according to the WAM-formulation
!     KM01   : Mean average wavenumber according to first order moment
!     SIGM_10: Mean frequency according to zeroth order moment
!     SIGM01 : Mean frequency according to first order moment
!
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     SINTGRL: Calculating integral paramters
!     SWCAP   : Calculating whitecapping source term
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
!     ---
!
! 13. Source text
!

      REAL :: ACTOT , EDRKTOT , EKTOT , ETOT1 , ETOT2 , ETOT4             40.22
      REAL :: KM_WAM, KM01 , SIGM_10 , SIGM01                             40.22

      COMMON /TP_WCAP/ ACTOT , EDRKTOT , EKTOT , ETOT1 , ETOT2 , ETOT4,   40.22
     &              KM_WAM, KM01 , SIGM_10 , SIGM01                       40.22
!$OMP THREADPRIVATE(/TP_WCAP/)                                            40.22
