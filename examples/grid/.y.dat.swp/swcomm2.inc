!                 SWAN COMMON BLOCKS, file 2 of 4
!
!     Contents of this file
!
!     COMMON / SWGRID / *** location and dimensions of input grids
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
!     40.03, 40.13: Nico Booij
!
!  1. Updates
!
!     30.72, Sept 97: Created from the old swancomm.ftn file
!     40.00, Jan. 98: NBFILS, NBSPEC, NBGRPT added
!                     common areas TIMFIL, INPGRS, INPFRM added
!                     NSTATM replaces NONSTA (moved to swcomm3.inc)
!            Feb. 99: IDYNCU, IDYNFR, IDYNWI, IDYNWL removed;
!                     DYNDEP added
!     40.03, May  00: air-sea temperature difference added
!     40.13, Jan. 01: COSPG is initialised at 1 (corresponding to ALPG)
!                     two comments corrected
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
! NUMGRD [ 10] maximum number of grids                                    40.03
!
      INTEGER             NUMGRD
      PARAMETER           (NUMGRD = 10)                                   40.03
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
!     COMMON / SWGRID /  *** location and dimensions of input grids ***
!
!     Common variables
!
! The following input grids are known in SWAN:
! grid(1): depth
! grid(2): current velocity x-component
! grid(3): current velocity y-component
! grid(4): friction coefficient
! grid(5): wind velocity x-component
! grid(6): wind velocity y-component
! grid(7): water level
! grid(8): vegetation density (not used yet)
! grid(9): air-sea temperature difference                                 40.03
!
!                [??] default value
!
! ICOND          [ 0] initial conditions
!                     =0 when mode stationary, or no initial conditions needed
!                     =1 when mode non-stationary and initial conditions should
!                        be calculated
! IGTYPE(NUMGRD) [  ] =0 when grid has constant values
!                     =1 when grid is regular
!                     =2 when grid is curvilinear
! LEDS(NUMGRD)   [ 0] =0 when values have not been read
!                     =1 if values were read
! MXG(NUMGRD)    [ 0] number of meshes in x-direction
! MYG(NUMGRD)    [ 0] number of meshes in y-direction
! NESRUN         [ 0] indicator (=1) for a nested run
! NWAMN          [ 0] indicator (=1) for a WAM-nested run
! OPTG           [ 1] type of the computational grid
!                     =1 when regular
!                     =2 when irregular, but rectangular (not used)
!                     =3 when curvilinear
! NBFILS         [ 0] number of boundary condition files                  40.00
! NBSPEC         [ 0] number of boundary spectra                          40.00
! NBGRPT         [ 0] number of comp. grid points for which boundary      40.00
!                     condition holds
!
      INTEGER             ICOND
      INTEGER             IGTYPE(NUMGRD),           LEDS(NUMGRD)
      INTEGER             MXG(NUMGRD), MYG(NUMGRD), NESRUN,      NWAMN
      INTEGER             OPTG,        NBFILS,      NBSPEC,      NBGRPT   40.00
!
!               [???] default value
!
! CVLEFT        [   ] the (curv.lin.) comput. grid is left/right-oriented
! DYNDEP        [FAL] is True if depth varies with time                   40.00
! LXOFFS        [FAL] offset values were/were not initialized already
! VARAST        [FAL] air-sea temp. diff. is/is not variable over space   40.03
! VARFR         [FAL] friction coefficient is/is not variable over space
! VARWI         [FAL] wind velocity is/is not variable over space
! VARWLV        [FAL] water level is/is not variable over space
!
      LOGICAL             CVLEFT,      LXOFFS,      VARFR                 40.00
      LOGICAL             VARWI,       VARWLV,      DYNDEP                40.00
      LOGICAL             VARAST                                          40.03
!
!               [??????] default value
!
!               [CALCUL] calculated value
!
! ALPG(NUMGRD)  [    0.] direction of the x-axis w.r.t. user coordinates
! COSPG(NUMGRD) [    1.] cos of ALPG                                      40.13
! COSVC         [CALCUL] =COS(-ALPG(2)),
!                        cos of angle of current input grid w.r.t. user coordinates 40.13
! COSWC         [CALCUL] =COS(-ALPG(5)),
!                        cos of angle of wind input grid w.r.t. user coordinates 40.13
! DXG(NUMGRD)   [    0.] mesh size of input grid in x-direction
! DYG(NUMGRD)   [    0.] mesh size of input grid in y-direction
! EXCFLD(NUMGRD)[-1.E20] exception values for input grids
! RDTIM         [      ] =0 when in stationary mode
!                        =1/DT when in non-stationary mode
! SINPG(NUMGRD) [    0.] sin of ALPG
! SINVC         [CALCUL] =SIN(-ALPG(2)),
!                        sin of angle of current input grid w.r.t. comp. grid
! SINWC         [CALCUL] =SIN(-ALPG(5)),
!                        sin of angle of wind input grid w.r.t. comp. grid
! STAGX(NUMGRD) [    0.] staggering of curv.lin. inp. grid w.r.t. comp. grid in X
! STAGY(NUMGRD) [    0.] staggering of curv.lin. inp. grid w.r.t. comp. grid in Y
! XPG(NUMGRD)   [    0.] x of origin
! XOFFS         [    0.] offset value in x
!                        (from user coord. system to internal coord. system)
! YPG(NUMGRD)   [    0.] y of origin
! YOFFS         [    0.] offset value in y
!                        (from user coord. system to internal coord. system)
!
      REAL                ALPG(NUMGRD),             COSPG(NUMGRD)
      REAL                COSVC,       COSWC,       DXG(NUMGRD)
      REAL                DYG(NUMGRD), EXCFLD(NUMGRD)
      REAL                RDTIM,       SINPG(NUMGRD)
      REAL                SINVC,       SINWC,       STAGX(NUMGRD)
      REAL                STAGY(NUMGRD),            XOFFS
      REAL                XPG(NUMGRD), YOFFS,       YPG(NUMGRD)
!
!
      COMMON / SWGRID /   XPG,         YPG,         ALPG,        COSPG
      COMMON / SWGRID /   SINPG,       DXG,         DYG,         MXG
      COMMON / SWGRID /   MYG,         LEDS,        IGTYPE,      VARFR
      COMMON / SWGRID /   VARWI,       COSVC,       SINVC,       COSWC
      COMMON / SWGRID /   SINWC,       XOFFS,       YOFFS,       LXOFFS
      COMMON / SWGRID /   VARWLV,      DYNDEP,      NESRUN,      NWAMN    40.00
      COMMON / SWGRID /   OPTG,        STAGX,       STAGY,       CVLEFT
      COMMON / SWGRID /   RDTIM,       ICOND,       EXCFLD,      NBFILS
      COMMON / SWGRID /   NBSPEC,      NBGRPT,      VARAST                40.03
!
!     commons added version 40.00
!
!     IFLDYN     if =0: data is stationary, if =1: nonstationary
!     IFLIDL     lay-out in input file
!     IFLIFM     format identifier
!     IFLNHF     number of heading lines per file
!     IFLNHD     number of heading lines per input field
!     IFLNDS     unit ref number of data file
!     IFLNDF     unit ref number of namelist file
!
      INTEGER    IFLDYN
      INTEGER    IFLIDL, IFLIFM, IFLNHF, IFLNHD, IFLNDS, IFLNDF
!
!     IFLBEG     begin time of data on file
!     IFLINT     time interval of data on file
!     IFLEND     end time of data on file
!     IFLTIM     time of last reading
!     IFLFAC     multiplication factor
!
      REAL       IFLBEG, IFLINT, IFLEND, IFLTIM
      REAL       IFLFAC
!
!
      COMMON /INPGRS/ IFLIDL(NUMGRD), IFLIFM(NUMGRD), IFLNHF(NUMGRD)
      COMMON /INPGRS/ IFLNHD(NUMGRD), IFLFAC(NUMGRD), IFLNDS(NUMGRD)
      COMMON /INPGRS/ IFLNDF(NUMGRD), IFLDYN(NUMGRD), IFLTIM(NUMGRD)
      COMMON /INPGRS/ IFLBEG(NUMGRD), IFLINT(NUMGRD), IFLEND(NUMGRD)
!
!     IFLFRM     format string
!
      CHARACTER IFLFRM *40
      COMMON /INPFRM/ IFLFRM(NUMGRD)
!
