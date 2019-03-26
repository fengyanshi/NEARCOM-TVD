!                 SWAN COMMON BLOCKS, file 3 of 4
!
!     Contents of this file
!
!     COMMON / COMPDA / *** pointers for data arrays on computational grid
!     COMMON / SWCOMG / *** location and dimensions of computational grid
!     COMMON / TP_SWCOMG /  *** thread private components of SWCOMG       40.22
!     COMMON / SWFYSP / *** physical parameters
!     COMMON / SWNUMS / *** information related to the numerical scheme
!     COMMON / TP_SWNUMS /  *** thread private components of SWNUMS       40.22
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
!     30.75: Nico Booij
!     30.82: IJsbrand Haagsma
!     32.01: Roeland Ris & Cor van der Schelde
!     32.06: Roeland Ris
!     34.00: Nico Booij
!     40.00: Nico Booij
!     40.03: Nico Booij
!     40.14: Annette Kieftenburg
!     40.22: John Cazes and Tim Campbell
!     40.31: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Sept 97: Created from the old swancomm.ftn file
!     32.01, Jan. 98: MXAUX increased by 3 for the h3268 project
!     32.01, Jan. 98: Introduced JSETUP & JDPSAV in / COMPDA /
!     32.01, Jan. 98: Introduced LSETUP in / SWNUMS /
!     30.75, Mar. 98: SPDIR2 added to / SWCOMG /
!     31.04, Apr. 98: JWFRCX and JWFRCY added to / COMPDA /
!     34.00, May  98: array PQUAD added for Alkyon
!                     MWCAP increased to 15
!     40.00, Jan. 98: JBFILS, JBSPEC, JBTIME, JBGRID, JBSLOC, JBSDIR,
!                     JBSFRQ, JBSAUX added to common COMPDA
!                     spectral shape parameters added: FSHAPE, DSHAPE, PSHAPE
!                     array SPPARM added: holds integral parameters of
!                     incident spectrum
!            Mar. 98: JHSIBC added (stores Hs from boundary condition)
!            May  98: IXCGRD, IYCGRD, KCGRD added
!            June 98: NGRBND, XCGMIN etc. added to /SWCOMG/
!                     JGRBND added to /COMPDA/
!                     BNDCHK added to /SWNUMS/
!                     NSTATC replaces ISTAT, NSTATM replaces NONSTA
!            Febr 99: BRESCL added (rescaling on/off)
!     32.06, June 99: IGEN added
!     30.80, Aug. 99: URSELL added
!     30.82, Aug. 99: MNUMS increased by 5 (accuracy for SETUP)
!     30.82, Aug. 99: Array PSETUP added
!     30.82, Sep. 99: Flag CSETUP added
!     40.07, Nov. 99: ACUPDA added
!     33.10, Apr. 00: MNUMS increased by 5 (num. diffusion SORDUP scheme)
!     40.03, Apr. 00: Ursell number array added to Compda
!     40.03, May  00: air-sea temperature difference added
!     40.22, Sep. 01: Split MSWMAT into two arrays MSWMATR(real)          40.22
!                     and MLSWMAT(logical                                 40.22
!     40.22, Sep. 01: Moved IXCGRD, IYCGRD, and KCGRD into thread         40.22
!                     private common block TP_SWCOMG for OpenMP           40.22
!     40.22, Sep. 01: Moved ICMAX and CSETUP into thread                  40.22
!                     private common block TP_SWNUMS for OpenMP           40.22
!     40.31, Dec. 03: removing POOL mechanism
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
!        [???] value
!
! MBOT   [ 10] dimension of array PBOT
! MICMAX [ 10] max. number of points in comput. stencil                   40.00
! MNUMS  [ 30] dimension of array PNUMS                                   33.10
! MQUAD  [ 10] dimension of array PQUAD                                   40.07
! MSHAPE [  2] dimension of array MSETUP                                  30.82
! MSHAPE [  5] dimension of array PSHAPE                                  40.00
! MSPPAR [  5] dimension of array SPPARM                                  40.00
! MSURF  [ 10] dimension of array PSURF
! MTRIAD [ 10] dimension of array PTRIAD
! MWCAP  [ 15] dimension of array PWCAP                                   34.00
! MWIND  [ 40] dimension of array PWIND
!
!
      INTEGER             MBOT,        MICMAX,      MNUMS,       MQUAD
      INTEGER             MSETUP,      MSHAPE,      MSPPAR,      MSURF    30.82
      INTEGER             MTRIAD,      MWCAP,       MWIND                 40.31 30.82
!
      PARAMETER           (MBOT   = 10)
      PARAMETER           (MICMAX = 10)                                   40.00
      PARAMETER           (MNUMS  = 30)                                   33.10
      PARAMETER           (MQUAD  = 10)                                   40.07
      PARAMETER           (MSETUP =  2)                                   30.82
      PARAMETER           (MSHAPE =  5)                                   40.00
      PARAMETER           (MSPPAR =  5)                                   40.00
      PARAMETER           (MSURF  = 10)
      PARAMETER           (MTRIAD = 10)
      PARAMETER           (MWCAP  = 15)                                   34.00
      PARAMETER           (MWIND  = 40)
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
!  COMMON / COMPDA / *** pointers for data arrays on computational grid ***
!
!     Common variables
!
! J?????  sequence number of variable ????? within array *;
!         i.e. variable ????? in COMPDA is found in COMPDA(IX,IY,J?????)
!         or variable ????? in ARRAY is found in ARRAY(J?????)
!
! M?????  number of subarrays within an array
!
!        [???] default value, [CAL]:calculated
!
! JABIN  [  7] within array SWMATR
! JABLK  [  8] within array SWMATR
! JAOLD  [ 12] within array SWMATR
! JASTD2 [   ] new air-sea temp. diff. within array COMPDA                40.03
! JASTD3 [   ] last read air-sea temp. diff. within array COMPDA          40.03
! JBTIME [   ] ??                                                         40.00
! JCDRAG [  1] drag coefficient within array COMPDA,
!              set by command WCAP JANS ...
! JDHS   [  6] wave height correction within array COMPDA
!              (difference in Hs between last two iterations)
! JDIS0  [  9] within array SWMATR
! JDIS1  [ 10] within array SWMATR
! JDISS  [  2] dissipation within array COMPDA
! JDPSAV [  1] saved depth (for setup) within array COMPDA                32.01
! JDP1   [  7] old depth within array COMPDA
! JDP2   [  8] new depth within array COMPDA
! JDP3   [ 15] last read depth within array COMPDA
! JDTM   [ 20] wave period correction within array COMPDA                 40.00
!              (difference in average wave period between last two iterations)
! JFRC2  [  1] friction coefficcient within array COMPDA                  40.00
!              set by command READ FR ...
! JFRC3  [  1] friction coefficcient within array COMPDA                  40.00
!              set by command READ FR ...
! JHS    [  1] significant wave height Hs within array COMPDA
! JHSIBC [ 25] significant wave height from boundary condition in array COMPDA 40.00
! JLEAK  [ 21] "leak" within array COMPDA                                 40.00
!              (refractive energy tranport over sector boundaries)
! JLEK1  [ 11] within array SWMATR
! JMAT5  [  5] within array SWMATR
! JMAT6  [  6] within array SWMATR
! JMATD  [  1] within array SWMATR
! JMATL  [  3] within array SWMATR
! JMATR  [  2] within array SWMATR
! JMATU  [  4] within array SWMATR
! JP4D   [  7] within array SWTSDA, quadruplet interactions
! JP4S   [  6] within array SWTSDA, quadruplet interactions
! JPBTFR [  4] within array SWTSDA, bottom friction
! JPTRI  [  8] within array SWTSDA, triad interactions
! JPWBRK [  5] within array SWTSDA, surf breaking
! JPWCAP [  3] within array SWTSDA, white capping
! JPWNDA [  1] within array SWTSDA, wind source term part A
! JPWNDB [  2] within array SWTSDA, wind source term part B
! JQB    [  4] fraction of breaking waves within array COMPDA
! JSETUP [  1] setup values within array COMPDA                           32.01
! JSTP   [  5] steepness within array COMPDA
! JTAUW  [  1] TauW within array COMPDA,
!              set by command WCAP JANS ...
! JUBOT  [  3] bottom orbital velocity within array COMPDA
! JURSEL [  1] Ursell number as used in Triad computation                 40.03
! JUSTAR [  1] friction velocity within array COMPDA,
!              set by command WCAP JANS ...
! JVX1   [  9] x of old current velocity within array COMPDA
! JVX2   [ 11] x of new current velocity within array COMPDA
! JVX3   [ 13] x of last read current velocity within array COMPDA
! JVY1   [ 10] y of old current velocity within array COMPDA
! JVY2   [ 12] y of new current velocity within array COMPDA
! JVY3   [ 14] y of last read current velocity within array COMPDA
! JWFRCX [  0] within array COMPDA: x-comp. wave induced force            31.04
! JWFRCY [  0] within array COMPDA: y-comp. wave induced force            31.04
! JWLV1  [ 22] old water level within array COMPDA                        40.00
! JWLV2  [ 24] new water level within array COMPDA                        40.00
! JWLV3  [ 23] last read water level within array COMPDA                  40.00
! JWX2   [ 16] x of new wind velocity within array COMPDA                 40.00
! JWX3   [ 18] x of last read wind velocity within array COMPDA           40.00
! JWY2   [ 17] y of new wind velocity within array COMPDA                 40.00
! JWY3   [ 19] y of last read wind velocity within array COMPDA           40.00
! JZEL   [  1] roughness within array COMPDA,
!              set by command WCAP JANS ...
!
! MCMVAR [ 25] within array COMPDA,                                       40.00
!              =MCMVAR+1, for command READ BOT ... (add JDP2) (obsolete)
!              =MCMVAR+2, for command READ FR ... (add JFRC2, JFRC3)      40.00
!              =MCMVAR+1, for command READ WL ... (add JWLV2) (obsolete)
!              =MCMVAR+2, for command WIND ... (add JWX2, JWY2) (obsolete)
!              =MCMVAR+2, for command READ WI ... (add JWX2, JWY2) (obsolete)
!              =MCMVAR+2, for command READ CUR ... (add JVX2, JVY2) (obsolete)
!              =MCMVAR+2, for command READ COOR ... (add JCOOX, JCOOY)
!              =MCMVAR+4, for command WCAP JANS ... (add JCDRAG, JTAUW,
!               JUSTAR, JZEL)
!              =MCMVAR+1, for command TRIAD .......(add Ursell)           40.03
! MSWMAT [ 12] within array SWMATR
! MTSVAR [  8] within array TESTDA
!
      INTEGER             JABIN,       JABLK,       JAOLD                 40.31
      INTEGER             JASTD2,      JASTD3                             40.14 40.03
      INTEGER             JBTIME,      JCDRAG,      JDHS,        JDIS0    40.31
      INTEGER             JDIS1,       JDISS,       JDPSAV,      JDP1     32.01
      INTEGER             JDP2                                            32.01
      INTEGER             JDP3,        JDTM,        JFRC2,       JFRC3    40.31
      INTEGER             JHS,         JHSIBC,      JLEAK                 40.31 40.00
      INTEGER             JLEK1,       JMAT5,       JMAT6,       JMATD
      INTEGER             JMATL,       JMATR,       JMATU
      INTEGER             JP4D,        JP4S,        JPBTFR,      JPTRI
      INTEGER             JPWBRK,      JPWCAP,      JPWNDA,      JPWNDB
      INTEGER             JQB,         JSETUP,      JSTP,        JTAUW    32.01
      INTEGER             JUBOT,       JUSTAR,      JVX1,        JVX2
      INTEGER             JVX3,        JVY1,        JVY2,        JVY3
      INTEGER             JWFRCX,      JWFRCY                             40.31 31.04
      INTEGER             JWLV1,       JWLV2,       JWLV3
      INTEGER             JWX2,        JWX3                               40.00
      INTEGER             JWY2,        JWY3,        JZEL
      INTEGER             MCMVAR,                   MTSVAR,      JURSEL   40.03 40.22
      INTEGER             MSWMATR,     MLSWMAT                            40.22
!
      COMMON / COMPDA /   MCMVAR,      JHS,         JDISS
      COMMON / COMPDA /   JUBOT,       JQB,         JSTP,        JDHS
      COMMON / COMPDA /   JDP1,        JVX1,        JVY1,        JDP2
      COMMON / COMPDA /   JVX2,        JVY2,        JFRC2,       JFRC3    40.00
      COMMON / COMPDA /   JWX2,        JWY2,        JWLV1,       JWLV2    40.31
      COMMON / COMPDA /   MTSVAR,      JPWNDA,      JPWNDB,      JPWCAP
      COMMON / COMPDA /   JPBTFR,      JPWBRK,      JP4S,        JP4D
      COMMON / COMPDA /   JPTRI,       JDTM                               40.31 40.22
      COMMON / COMPDA /   JMATD,       JMATR,       JMATL,       JMATU
      COMMON / COMPDA /   JMAT5,       JMAT6,       JABIN,       JABLK
      COMMON / COMPDA /   JDIS0,       JDIS1,       JLEK1,       JAOLD
      COMMON / COMPDA /   JLEAK,       JWLV3,       JVX3,        JVY3
      COMMON / COMPDA /   JWX3,        JWY3,        JDP3
      COMMON / COMPDA /   JSETUP,      JDPSAV                             40.31 40.00
      COMMON / COMPDA /   JWFRCX,      JWFRCY,      JUSTAR,      JZEL     40.00
      COMMON / COMPDA /   JTAUW,       JCDRAG,      JBTIME,      JHSIBC   40.31
      COMMON / COMPDA /   JURSEL,      JASTD2,      JASTD3                40.14 40.03
      COMMON / COMPDA /   MSWMATR,     MLSWMAT                            40.22
!
!-------------------------------------------------------------------------------
!
!     COMMON / SWCOMG /  *** location and dimensions of computational grid ***
!
!     Common variables
!
! J?????  sequence number of variable ????? within array *;
!         or variable ????? in ARRAY is found in ARRAY(J?????)
!
! M?????  number of ?????
!
!        [???] default value, [CAL]:calculated
!
! ICOMP  [  1] unused
! IXCGRD [   ] IX of points of computational stencil                      40.00
! IYCGRD [   ] IY of points of computational stencil                      40.00
!
! KCGRD  [   ] grid address of points of computational stencil            40.00
! MCGRD  [  1] number of wet grid points of the computational grid
! MDC    [  0] grid points in theta-direction of computational grid
! MDC4MA [CAL] =IDHGH; some counter for quadruplet interactions. Stored in WWINT(18)
! MDC4MI [CAL] =IDLOW; some counter for quadruplet interactions. Stored in WWINT(17)
! MMCGR  [CAL] =MXC*MYC; grid points in computational grid
! MSC    [  0] grid points in sigma-direction of computational grid
! MSC4MA [CAL] =ISHGH; some counter for quadruplet interactions. Stored in WWINT(16)
! MSC4MI [CAL] =ISLOW; some counter for quadruplet interactions. Stored in WWINT(15)
! MTC    [  1] computational time steps
! MXC    [  0] grid points in x-direction of computational grid
! MYC    [  0] grid points in y-direction of computational grid
!
! NGRBND [cal] number of grid points on computational grid boundary       40.00
! NX     [CAL] =MXC-1; only used locally. Equal to MXS
! NY     [CAL] =MYC-1; only used locally. Equal to MYS
!
      INTEGER             IXCGRD(MICMAX), IYCGRD(MICMAX), KCGRD(MICMAX)   40.00
      INTEGER             ICOMP,       MCGRD
      INTEGER             MDC,         MDC4MA,      MDC4MI,      MMCGR
      INTEGER             MSC,         MSC4MA,      MSC4MI,      MTC
      INTEGER             MXC,         MYC,         NX,          NY
      INTEGER             NGRBND                                          40.00
!
!        [???] default value, [CAL]:calculated
!
! FULCIR [TRU] spectral directions cover the full/part of circle
!
      LOGICAL             FULCIR
!
! *CP and *PC  to transform between computational- and user-coordinates
!
!        [???] default value, [CAL]:calculated
!
! ALCP   [CAL] =-ALPC;
!              direction of user coordinates w.r.t. computational coordinates
! ALPC   [ 0.] direction of x-axis of computational grid w.r.t. user coordinates
! COSPC  [CAL] =COS(ALPC)
! DX     [CAL] =XCLEN/MXS; mesh size in x-direction of computational grid
! DY     [CAL] =YCLEN/MYS; mesh size in y-direction of computational grid
! DDIR   [CAL] =(SPDIR2-SPDIR1)/MDC;
!              mesh size in theta-direction of computational grid
! DXRP   [   ] unused
! DYRP   [   ] unused
! FRINTF [CAL] =ALOG(SHIG/SLOW)/(MSC-1); frequency integration factor (df/f)
!              (integral over frequency of G(f) = SUM_j (sigma_j*Gj*FRINTF) )
! FRINTH [CAL] =SQRT(SFAC); frequency mesh boundary factor
!              (mesh in frequency space runs from sigma/FRINTH to sigma*FRINTH)
! SHIG   [CAL] =2*PI*FRHIG; highest spectral value of sigma
! SINPC  [CAL] =SIN(ALPC)
! SLOW   [CAL] =2*PI*FRLOW; lowest spectral value of sigma
! SPDIR1 [ 0.] represents first spectral direction (if FULCIR=.FALSE.)
! SPDIR2 [   ] represents second spectral direction (if FULCIR=.FALSE.)   30.75
! XCGMAX [cal] maximum x-coordinate of computational grid points          40.00
! XCGMIN [cal] minimum x-coordinate of computational grid points          40.00
! XCLEN  [   ] length of computational grid in x-direction
! XCP    [CAL] =-XPC*COSPC-YPC*SINPC;
!              origin of user coordinates w.r.t. computational coordinates
! XPC    [ 0.] x coordinate of origin of computational grid
! YCGMAX [cal] maximum y-coordinate of computational grid points          40.00
! YCGMIN [cal] minimum y-coordinate of computational grid points          40.00
! YCLEN  [   ] length of computational grid in y-direction
! YCP    [CAL] =XPC*SINPC-YPC*COSPC;
!              origin of user coordinates w.r.t. computational coordinates
! YPC    [ 0.] y coordinate of origin of computational grid
!
      REAL                ALCP,        ALPC,        COSPC,       DDIR
      REAL                DX,          DXRP,        DY,          DYRP
      REAL                FRINTF,      FRINTH,      SHIG,        SINPC
      REAL                SLOW,        SPDIR1,      SPDIR2,      XCLEN    30.75
      REAL                XCP,         XPC,         YCLEN,       YCP
      REAL                YPC
      REAL                XCGMIN,      XCGMAX,      YCGMIN,      YCGMAX   40.00
!
      COMMON / SWCOMG /   ICOMP,       XPC,         YPC,         ALPC
      COMMON / SWCOMG /   COSPC,       SINPC,       XCLEN,       YCLEN
      COMMON / SWCOMG /   MTC,         MXC,         MYC,         MDC
      COMMON / SWCOMG /   MSC,         SLOW,        SHIG,        DX
      COMMON / SWCOMG /   DY,          DDIR,        NX,          NY
      COMMON / SWCOMG /   XCP,         YCP,         ALCP,        DXRP
      COMMON / SWCOMG /   DYRP,        MSC4MI,      MSC4MA,      MDC4MI
      COMMON / SWCOMG /   MDC4MA,      FRINTF,      FRINTH,      MMCGR
      COMMON / SWCOMG /   FULCIR,      SPDIR1
      COMMON / SWCOMG /   MCGRD,       SPDIR2                             30.75
      COMMON / SWCOMG /   XCGMIN,      XCGMAX,      YCGMIN,      YCGMAX   40.00
      COMMON / SWCOMG /   NGRBND                                          40.00

      COMMON / TP_SWCOMG /  IXCGRD,  IYCGRD,  KCGRD                       40.22
!$OMP THREADPRIVATE(/TP_SWCOMG/)                                          40.22
!
!-------------------------------------------------------------------------------
!
!     COMMON / SWFYSP /  *** physical parameters ***
!
!     Common variables
!
!        [????????????] default value
!
! CASTD  [          0.] (constant) air-sea temperature difference         40.03
! DEGRAD [  CALCULATED] =PI/180; constant to transform degrees to radians
! DNORTH [         90.] direction of the North w.r.t. x-axis of user coordinates
!                       =nor; set by command SET ... [nor] ...
! GRAV   [        9.81] acceleration due to gravity
!                       =grav; set by command SET ... [grav] ...
! PI     [3.1415926536] circular constant
! PI2    [  CALCULATED] =2*PI;
! PWTAIL(10)            coefficients to calculate tail of the spectrum
!    (1) [          4.] tail power of energy density spectrum as function of freq.
!                       =5.; for command GEN3 JANS ...
!                       =5.; for command WCAP JANS ...,
!                       not documented in manual
!                       =5.; for command GROWTH G3 JANS ...,
!                       not documented in manual
!                       =pwtail; set by command SET ... [pwtail]
!    (2) [         2.5] energy spectrum w.r.t. wave number, not used
!    (3) [  CALCULATED] =PWTAIL(1)+1; action density spectrum w.r.t. frequency
!    (4) [          3.] action density spectrum w.r.t. wave number, not used
!    (5) [  CALCULATED] =1./(PWTAIL(1)*(1.+PWTAIL(1)*(FRINTH-1.)));
!                       tail factor for the action integral
!    (6) [  CALCULATED] =1./((PWTAIL(1)-1.)*(1.+(PWTAIL(1)-1.)*(FRINTH-1.)));
!                       tail factor for the energy integral
!    (7) [  CALCULATED] =1./((PWTAIL(1)-2.)*(1.+(PWTAIL(1)-2.)*(FRINTH-1.)));
!                       tail factor for the first moment of energy
!    (8) [  CALCULATED] =1./((PWTAIL(1)-3.)*(1.+(PWTAIL(1)-3.)*(FRINTH-1.)));
!                       tail factor for the second moment of energy
! RHO    [       1025.] density of the water
!                       =rho; set by command SET ... [rho] ...
! WLEV   [          0.] water level
!                       =level; set by command SET [level] ...
!
      REAL                CASTD                                             40.03
      REAL                DEGRAD,      DNORTH,      GRAV,        PI
      REAL                PI2,         PWTAIL(10),  RHO,         WLEV
!
      COMMON / SWFYSP /   GRAV,        WLEV,        PI,          PI2
      COMMON / SWFYSP /   RHO,         DEGRAD,      DNORTH,      PWTAIL
      COMMON / SWFYSP /   CASTD                                             40.03
!
!-------------------------------------------------------------------------------
!
!     COMMON / SWNUMS /  *** information related to the numerical scheme ***
!
!     Common variables
!
!        [???] default value
!
! DSHAPE [  2] indicates option for computation of directional distribution
!              in the spectrum (boundary spectra etc.)
!              =1: directional spread in degrees is given
!              =2: power of COS is given
!
! FSHAPE [  2] indicates option for computation of frequency distribution
!              in the spectrum (boundary spectra etc.)
!              =2: Jonswap(default set by subr SWINIT),
!              =1: Pierson-Moskowitz, =3: bin, =4: Gaussian (set by command
!              BOUNshape ..)
! IBOT   [  0] indicator bottom friction:
!              =0; no bottom friction dissipation
!              =1; set by command FRIC JON  ..., Jonswap bottom friction model
!              =2; for command FRIC COLL ..., Collins bottom friction model
!              =3; for command FRIC MAD  ..., Madsen bottom friction model
! ICMAX  [  3] number of points in computational stencil
! ICOR   [   ] not used
! ICUR   [  0] indicates presence of currents:
!              =0; no currents
!              =1; for command READ CUR ..., currents are present
! IDBR   [  1] not used
! IDIF   [  0] not used
! IGEN   [  3] indicates the generation mode                              32.06
!              =1; for command GEN1 ...,
!              =2; for command GEN2 ...,
!              =3; for command GEN3 ...,
! IINC   [  0] not used
! IPRE   [   ] not used
! IQUAD  [  2] indicates quadruplet interaction term:
!              =0; for command OFF QUAD
!              =0; for command GEN1 ...,
!              =0; for command GEN2 ...,
!              =0; for command GROWTH G1 ... (not documented in manual),
!              =0; for command GROWTH G2 ... (not documented in manual),
!              quadruplets are inactive
!              =1; quadruplets are calculated semi implicit per sweep direction
!              =2; for command GEN3 ...,
!              =2; for command QUAD, not documented in the manual
!              =2; set when IWIND=3 or 4 and ICUR=0 in SUBR ERRCHK,
!              quadruplets are calculated fully explicit per sweep direction
!              =3; set when IWIND=3 or 4 and ICUR=1 in SUBR ERRCHK,
!              quadruplets are calculated fully explicit per iteration
!              =iquad; set by command GEN3 ... QUAD [iquad] ...,
! IREFR  [  1] indicates refraction effect:
!              =0; for command OFF REF, refraction is inactive
!              =1; refraction is active
! ISURF  [  1] indicates surf breaking (shallow water) term:
!              =0; for command OFF BRE, surf breaking is inactive
!              =1; for command BRE CON ..., surf breaking with constant parameter
!              =2; for command BRE VAR ..., surf breaking according to Nelson
! ITERMX [   ] maximum number of iterations:
!              is set equal to MXITST in case of stationary computations     40.03
!              is set equal to MXITNS in case of nonstationary computations  40.03
! ITFRE  [  1] indicator for transport of action in frequency space
!              =0; for command OFF FSH, frequency shifting inactive
!              =1; frequency shifting active
! ITRIAD [  0] indicates triad interaction term:
!              =0; triads are inactive
!              =1; for command TRI DTA IMP ..., not documented in manual
!              =2; for command TRI DTA EXP ..., not documented in manual
!              =3; for command TRI [trfac] [cutfr], as in manual
!              =3; for command TRI LTA IMP ..., not documented in manual
!              =4; for command TRI LTA EXP ..., not documented in manual
! ITRSY  [  1] not used
! IWCAP  [  1] indicates whitecapping:
!              =0; for command GEN1 ...,
!              =0; for command GEN2 ...,
!              =0; for command OFF WCAP, no whitecapping
!              =1; for command GEN3 KOM ...,
!              =1; for command WCAP KOM ..., not documented in manual,
!              standard WAM formulation (Komen et al.; 1984)
!              =2; for command GEN3 JANS ...,
!              =2; for command WCAP JANS ..., not documented in manual,
!              according to Janssen (1989, 1991)
!              =3; for command WCAP LHIG ..., not documented in manual,
!              according to Longuet-Higgins (1967), Yuan et al. (1986)
!              =4; for command WCAP BJ ..., not documented in manual,
!              according to Battjes & Janssen (1978)
!              =5; for command WCAP KBJ ..., not documented in manual,
!              combined formulation of Komen (1) and Battjes & Janssen (4)
! IWIND  [  0] indicates presence of wind, and type of source term used:
!              =0; no wind
!              =1; for command GEN1 ..., if wind is made active,
!              =1; for command GROWTH G1 ..., not documented in manual,
!              1st generation source term
!              =2; for command GEN2 ..., if wind is made active,
!              =2; for command GROWTH G2 ..., not documented in manual,
!              2nd generation source term (as in Dolphin)
!              =3; for command WIND ..., if IWIND still was 0, else unchanged,
!              =3; for command GEN3 KOM ..., if wind is made active,
!              =3; for command GROWTH G3 KOM ..., not documented in manual,
!              3rd generation source term (Snyder)
!              =4; for command GEN3 JANS ..., if wind is made active,
!              =4; for command GROWTH G3 JANS ..., not documented in manual,
!              source term by P. Janssen (1989, 1991)
!              =5; for command GEN3 YAN ..., if wind is made active,
!              =5; for command GROWTH G3 YAN ..., not documented in manual,
!              not documented in manual
! LSETUP [  0] =0; setup is not calculated                                32.01
!              =1; setup is calculated                                    32.01
!              =2; setup is calculated with the boundary conditions from  32.01
!              a nest file                                                32.01
! MXITST [ 15] max. number of iterations in stationary computations       40.03
! MXITNS [  1] max. number of iterations in nonstationary computations    40.03
! NCOR   [  1] not used
! NSTATC [  1] indicates stationarity of computation:
!              =0; stationary computation
!              =1; nonstationary computation
! NSTATM  [-1] 0: stationary mode, 1: nonstationary mode, -1: unknown
! NUMOBS [   ] number of obstacles
! NCOMPT [   ] number of COMPUTE commands
!
      INTEGER             DSHAPE,      FSHAPE                             40.00
      INTEGER             IBOT,        ICMAX,       ICOR,        ICUR
      INTEGER             IDBR,        IDIF,        IGEN,        IINC     32.06
      INTEGER             IPRE,        IQUAD,       IREFR,       ISURF    32.06
      INTEGER             ITERMX,      ITFRE,       ITRIAD,      ITRSY
      INTEGER             IWCAP,       IWIND,       LSETUP                32.01
      INTEGER             MXITST,      MXITNS,                   NCOR     40.03
      INTEGER             NSTATC,      NSTATM,      NUMOBS,      NCOMPT   40.00
!
! ACUPDA [true] indicates whether or not action densities are to be updated
!               during computation                                        40.07
!
      LOGICAL             ACUPDA                                          40.07
!
!        [???????] default value, [CALCULAT] is calculated
!
! DEPMIN [   0.05] threshold depth (to prevent zero divisions)
!                  =depmin; set by command SET ... [depmin] ...
! HSRERR [    0.1] The error margin allowed between pre-scribed and calculated Hs
!                  at the up-wave boundary. If exceeded a warning is produced
! PBOT(MBOT)       coefficients for the bottom friction models
!  ( 1)  [     0.] =cfc; set by command FRIC COL [cfw] [cfc], (Collins equation),
!                  not documented in the manual
!  ( 2)  [  0.015] =cfw; set by command FRIC COL [cfw], (Collins equation),
!                  also used as CFW in the source code
!  ( 3)  [  0.067] =cfjon; set by command FRIC JON [cfjon], (Jonswap equation)
!  ( 4)  [  -0.08] =mf; value cannot be changed, (Madsen equation)
!  ( 5)  [   0.05] bottom roughness length scale, (Madsen equation),
!                  =kn; set by command FRIC JON [kn],
!                  also used as AKN in the source code
! PNUMS(MNUMS)     numerical coefficients
!                  accuracy criterion:
!  ( 1)  [   0.03] relative error in Hs and Tm01
!                  =drel; set by command NUM ACCUR [drel] ...
!  ( 2)  [   0.03] absolute error in Hs (m)
!                  =dhabs; set by command NUM ACCUR ... [dhabs] ...
!  ( 3)  [    0.3] absolute error in Tm01 (s)
!                  =dtabs; set by command NUM ACCUR ... [dtabs] ...
!  ( 4)  [  97.00] percentage of wet grid points were absolute and relative
!                  accuracy has been reached
!                  =npnts; set by command NUM ACCUR ... [npnts] ...
!  ( 5)  [     0.] not used
!                  diffusion schemes:
!  ( 6)  [    0.5] numerical diffusion over theta,
!                  =cdd; set by command NUM DIR [cdd]
!  ( 7)  [    0.5] numerical diffusion over sigma,
!                  =css; set by command NUM SIGIM [css] ...
!  ( 8)  [     1.] numerical scheme in frequency space
!                  =1.; for command NUM SIGIM, implicit scheme
!                  =2.; for command NUM SIGEX, explicit scheme,
!                  with CFL criterion
!                  =3.; for command NUM FIL, explicit scheme,
!                  without CFL criterion, not documented in the manual
!  ( 9)  [   0.01] diffusion coefficient for explicit scheme,
!                  =diffc; set by command NUM FIL [diffc],
!                  not documented in the manual
!                  linear CGSTAB solver schemes:
!  (10)  [    -3.] value for the pre-conditioner (=IINSOL(2)),
!                  = 0.; no pre-conditioner
!                  =-1.; diagonal pre-conditioner
!                  =-3.; ILU pre-conditioner
!                  =prec; set by command NUM SIGIM ... [prec] ...
!  (11)  [  1.E-6] termination criterion for iterative solver (=RINSOL(1)),
!                  =eps1; set by command NUM SIGIM ... [eps1] ...
!  (12)  [  1.E-4] termination criterion for iterative solver (=RINSOL(2)),
!                  =eps2; set by command NUM SIGIM ... [eps2] ...
!                  (||Ax-b|| < eps1 + eps2 * ||Ax-b|| + eps2 * ||b||)
!  (13)  [     0.] output for CGSTAB solver (=IISOL(3))
!                  =outp; set by command NUM SIGIM ... [outp] ...
!                  <0. no output
!                  =0. only fatal errors are printed
!                  =1. additional information about the iteration is printed
!                  =2. maximal output regarding the iteration process
!  (14)  [    20.] maximum number of iterations in each of the solution methods
!                  =niter; set by command NUM SIGIM ... [niter], (=IINSOL(6))
!  (15)  [   0.01] global error in Hs
!  (16)  [   0.01] global error in Tm01
!  (17)  [     0.] coefficient for limitation of Ctheta (not used currently)
!  (18)  [    0.8] limitation on Froude number (current velocity is reduced if
!                  larger than CGMAX=PNUMS(18)*SQRT(GRAV*DEPW)),
!                  =froudmax; set by command SET ... [froudmax] ...,
!                  not documented in the manual
!  (19)  [CALCULA] =0.5*SQRT(2.); CFL criterion for explicit scheme in
!                  frequency space, manual mentions default = 0.7,
!                  =cfl; set by command NUM SIGEX [cfl]
!  (20)  [  1.E20] maximum growth in spectral bin,
!                  =0.1; for command WIND ..., if IWIND=3 or 4
!                  =0.1; for command GEN3
!                  =0.1; for command TRI
!                  =0.1; for command QUAD,
!                  not documented in the manual
!                  =1.E20; for command OFF QUAD
!                  =limiter; set by command GEN3 ... QUAD ... [limiter]
!                  =limiter; set by command NUM ACCUR ... [itermax] [limiter],
!                  not documented in the manual
!                  =limiter; set by command QUAD [iquad] [limiter],
!                  not documented in the manual
!  (21)            not used
!  (22)            =prec for SETUP; value for the pre-conditioner (=IINSOL(2)),
!                  = 0.; no pre-conditioner
!                  =-1.; diagonal pre-conditioner
!                  =-3.; ILU pre-conditioner
!  (23)            =eps2 for SETUP; termination criterion for iterative solver
!                  (=RINSOL(2)), set by command NUM SETUP ... [eps2] ...
!                  (||Ax-b|| < eps2 * ||Ax-b|| + eps2 * ||b||)
!  (24)            =outp for SETUP; output for CGSTAB solver (=IISOL(3))
!                  <0. no output
!                  =0. only fatal errors are printed
!                  =1. additional information about the iteration is printed
!                  =2. maximal output regarding the iteration process
!  (25)            =niter for SETUP
!                  maximum number of iterations in each of the solution methods
! PQUAD(MQUAD)     coefficients for quadruplet interaction                34.00
!  ( 1)  [   0.25] lambda in eq. B29 of user manual                       34.00
!  ( 2)  [   3.E7] coefficient of interactions                            34.00
!  ( 3)  [    5.5] coefficient for shallow water interactions             34.00
!  ( 4)  [  0.833] coefficient for shallow water interactions             34.00
!  ( 5)  [  -1.25] coefficient for shallow water interactions             34.00
! PSETUP(MSETUP)                                                          30.82
!  ( 1)  [    0.0] not used
!  ( 2)  [    0.0] user defined level for correction of the setup         30.82
! PSHAPE(MSHAPE)   coefficients for calculation of spectrum from integral
!                  parameters                                             40.00
!  ( 1)  [    3.3] peak enhancement factor of Jonswap spectrum
!                  =gamma; set by command BOUN SHAPE .. JON [gamma]
!  ( 2)  [    0.1] width of Gaussian spectrum (in Hz)
!                  =sigfr; set by command BOUN SHAPE .. GAU [sigfr]
!                  after reading the value is converted to radians/second
!                  by a 2 pi multiplication
!
! PSURF(MSURF)     surf breaking coefficients
!  ( 1)  [    1.0] coef. for determining rate of dissipation, (Battjes Jansen),
!                  =1.5; for command BRE VAR
!                  =alpha, set by command BRE CON [alpha] ...
!                  =alpha, set by command BRE VAR [alpha] ...
!  ( 2)  [   0.73] breaker parameter
!                  =gamma, set by command BRE CON ... [gamma]
!  ( 4)  [       ] the min. value of the breaker parameter of Nelson
!                  =0.55; set by command BRE VAR
!                  =gammin; set by command BRE VAR ... [gammin] ...
!  ( 5)  [       ] the max. value of the breaker parameter of Nelson
!                  =0.81; set by command BRE VAR
!                  =gammax; set by command BRE VAR ... [gammax] ...
!  ( 6)  [       ] breaker parameter for negative bottom slopes
!                  =0.73; set by command BRE VAR
!                  =gamneg; set by command BRE VAR ... [gamneg] ...
!  ( 7)  [       ] proportionality coefficient in expression of Nelson
!                  =0.88; set by command BRE VAR
!                  =coeff1; set by command BRE VAR ... [coeff1] ...
!  ( 8)  [       ] coefficient in the exp in the expression of Nelson
!                  =0.012; set by command BRE VAR
!                  =coeff2; set by command BRE VAR ... [coeff2]
! PTRIAD(MTRIAD)
!  ( 1)  [    0.5] controls the proportionality coefficient,
!                  =0.25; for command TRIAD LTA
!                  =0.5; for command TRIAD DTA
!                  =trfac; set by command TRIAD ... [trfac] ...
!  ( 2)  [       ] controls the maximum frequency considered in the comp.,
!                  =2.5; for command TRIAD ...
!                  =cutfr; set by command TRIAD ... [cutfr]
!  ( 3)  [    0.1] controls above which Ursell number the quadruplets (and
!                  thus the limiter) are switched off
!                  =ursell; set by command LIM [ursell]
!                  =urcrit; set by command TRIAD ... [urcrit]
!  ( 4)  [    0.2] controls above which Ursell number the triads are
!                  calculated
!                  =urslim; set by command TRIAD ... [urslim]
! PWCAP(MWCAP)     whitecapping coefficients
!  ( 1)  [2.36E-5] coefficient for Komen et al. (1984),
!                  ALFAWC (Emperical coefficient)
!                  =cds2; set by command GEN3 KOM [cds2] ...,
!                  =cds2; set by command WCAP KOM [cds2] ...,
!                  not documented in the user manual
!  ( 2)  [3.02E-3] coefficient for Komen et al. (1984),
!                  ALFAPM (Alpha of Pierson Moskowitz frequency)
!                  =stpm; set by command GEN3 KOM ... [stpm],
!                  =stpm; set by command WCAP KOM ... [stpm],
!                  not documented in the user manual
!  ( 3)  [    4.5] coeff. for Janssen (1989,1991), acc. to Komen et al. (1994),
!                  CFJANS (cds coefficient)
!                  =cds1; set by command GEN3 JANS [cds1] ...,
!                  =cds1; set by command WCAP JANS [cds1] ...,
!                  not documented in the user manual
!  ( 4)  [    0.5] (=DELTA)
!                  =delta; set by command GEN3 JANS ... [delta],
!                  =delta; set by command WCAP JANS ... [delta],
!                  not documented in the user manual
!  ( 5)  [     1.] coefficient of Longuet Higgins,
!                  =cflhig; set by command WCAP LHIG [cflhig],
!                  not documented in the user manual
!  ( 6)  [   0.88] GAMBTJ (Steepness limited wave breaking)
!                  =bjstp; set by command WCAP BJ [bjstp] ...
!                  =bjstp; set by command WCAP KBJ [bjstp] ...
!  ( 7)  [     1.] Alpha in Battjes/Janssen,
!                  =bjalf; set by command WCAP BJ ... [bjalf]
!                  =bjalf; set by command WCAP KBJ ... [bjalf] ...
!  ( 8)  [   0.75] numerical diffusion over sigma
!                  =kconv; set by command WCAP KBJ ... [kconv]
! PWIND(MWIND)     wind growth term coefficients
!  ( 1)  [   188.] controls linear wave growth,
!                  =cf10; set by command GEN1 [cf10] ...
!                  =cf10; set by command GEN2 [cf10] ...
!                  =cf10; set by command GROWTH G1 [cf10] ...,
!                  not documented in the user manual
!                  =cf10; set by command GROWTH G2 [cf10] ...,
!                  not documented in the user manual
!  ( 2)  [   0.59] controls the exponential wave growth,
!                  =cf20; set by command GEN1 ... [cf20] ...
!                  =cf20; set by command GEN2 ... [cf20] ...
!                  =cf20; set by command GROWTH G1 ... [cf20] ...,
!                  not documented in the user manual
!                  =cf20; set by command GROWTH G2 ... [cf20] ...,
!                  not documented in the user manual
!  ( 3)  [   0.12] controls the exponential wave growth,
!                  =cf30; set by command GEN1 ... [cf30] ...
!                  =cf30; set by command GEN2 ... [cf30] ...
!                  =cf30; set by command GROWTH G1 ... [cf30] ...,
!                  not documented in the user manual
!                  =cf30; set by command GROWTH G2 ... [cf30] ...,
!                  not documented in the user manual
!  ( 4)  [   250.] controls the dissipation rate,
!                  =cf40; set by command GEN1 ... [cf40] ...
!                  =cf40; set by command GEN2 ... [cf40] ...
!                  =cf40; set by command GROWTH G1 ... [cf40] ...,
!                  not documented in the user manual
!                  =cf40; set by command GROWTH G2 ... [cf40] ...,
!                  not documented in the user manual
!  ( 5)  [ 0.0023] controls the spectral energy of the limit spectrum
!                  =cf50; set by command GEN2 ... [cf50] ...
!                  =cf50; set by command GROWTH G2 ... [cf50] ...,
!                  not documented in the user manual
!  ( 6)  [ -0.223] controls the spectral energy of the limit spectrum
!                  =cf60; set by command GEN2 ... [cf60] ...
!                  =cf60; set by command GROWTH G2 ... [cf60] ...,
!                  not documented in the user manual
!  ( 7)  [     0.] cf70, not used
!  ( 8)  [  -0.56] cf80, not used
!  ( 9)  [CALCULA] density air / density water (=RHOAW),
!                  =PWIND(16)/RHO
!                  =rhoaw; set by command GROWTH G1 ... [rhoaw] ...,
!                  not documented in the user manual
!                  =rhoaw; set by command GROWTH G2 ... [rhoaw] ...,
!                  not documented in the user manual
!     PWIND(13) = PMLM     0.13     (  )
!  (10)  [ 0.0036] limit energy Pierson Moskowitz spectrum
!                  =edmlpm; set by command GEN1 ... [edmlpm] ...
!                  =edmlpm; set by command GEN2 ... [edmlpm] ...
!                  =edmlpm; set by command GROWTH G1 ... [edmlpm] ...,
!                  not documented in the user manual
!                  =edmlpm; set by command GROWTH G2 ... [edmlpm] ...,
!                  not documented in the user manual
!  (11)  [0.00123] drag coefficient
!                  =cdrag; set by command GEN1 ... [cdrag] ...
!                  =cdrag; set by command GEN2 ... [cdrag] ...
!                  =cdrag; set by command GROWTH G1 ... [cdrag] ...,
!                  not documented in the user manual
!                  =cdrag; set by command GROWTH G2 ... [cdrag] ...,
!                  not documented in the user manual
!  (12)  [    1.0] minimum wind velocity, relative to current at 10 m above msl
!                  =umin; set by command GEN1 ... [umin] ...
!                  =umin; set by command GEN2 ... [umin] ...
!                  =umin; set by command GROWTH G1 ... [umin] ...,
!                  not documented in the user manual
!                  =umin; set by command GROWTH G2 ... [umin] ...,
!                  not documented in the user manual
!  (13)  [   0.13] coefficient that determines the Pierson Moskowitz spectrum
!                  =cfpm; set by command GEN1 ... [cfpm]
!                  =cfpm; set by command GEN2 ... [cfpm]
!                  =cfpm; set by command GROWTH G1 ... [cfpm],
!                  not documented in the user manual
!                  =cfpm; set by command GROWTH G2 ... [cfpm],
!                  not documented in the user manual
!  (14)  [   0.01] (=ALPHA) Alpha, according to Janssen (1991) wave growth model
!  (15)  [   0.41] (=XKAPPA)carnock: Kappa
!  (16)  [   1.28] density of the air (=RHOA)
!  (17)  [CALCULA] =RHO, density of the water (=RHOW)
!  (31)  [     0.] proportionality coefficient in the wave growth term of
!                  Caveleri and Malanotte,
!                  =0.0015; for command GEN3 ... AGROW
!                  =a; set by command GEN3 ... AGROW [a]
! BRESCL  [True]   rescaling on/off                                       40.00
! SIGMAG [    0.1] width of the Gaussian frequency spectrum in Hz
!                  =0.01; for command BOU STAT ... GAU
!                  =sigfr; set by command BOU STAT ... GAU [sigfr]
!                  after reading the value is converted to radians/second
!                  by a 2 pi multiplication
!
! SPPARM           integral parameters used for computation of incident spectrum
!  ( 1)  [    ---] significant wave height
!  ( 2)  [    ---] wave period (peak or mean)
!  ( 3)  [    ---] average wave direction
!  ( 4)  [    ---] directional distribution coefficient
!
! SY0    [    3.3] peak enhancement parameter of the JONSWAP spectrum,
!                  =gamma; set by command BOU STAT ... JON [gamma]
! U10    [     0.] wind velocity
!                  =vel; set by command WIND [vel] ...
! WDIC   [CALCULA] =PI2*((WDIP/PI2-NINT(WDIP/PI2)),
! WDIP   [     0.] wind direction with respect to problem coordinates
!                  =dir; set by command WIND ... [dir]
!
      REAL                DEPMIN,      PBOT(MBOT),  PNUMS(MNUMS)
      REAL                PSETUP(MSETUP),           PSHAPE(MSHAPE)        30.82
      REAL                PSURF(MSURF),             PTRIAD(MTRIAD)
      REAL                PWCAP(MWCAP),             PWIND(MWIND)
      REAL                SIGMAG
      REAL                SPPARM(MSPPAR),SY0,         U10,         U10C   40.00
      REAL                WDIC,          WDIP,        HSRERR
      REAL                PQUAD(MQUAD)                                    30.80
      REAL                RCOMPT(10,5)
!
! BNAUT  [.FALSE.] Indicates whether nautical or cartesian directions are used
! BNDCHK [.TRUE.]  indicates whether computed Hs on boundary must be compared
!                  with value entered as boundary condition               40.00
! ONED   [.FALSE.] Indicates whether the calculation should be performed in 1D-mode
! CSETUP [.TRUE. ] Incicates whether solver for setup has converged       30.82
!
      LOGICAL             BNDCHK,      BNAUT,       ONED,        BRESCL   40.00
      LOGICAL             CSETUP                                          30.82
!
      COMMON / SWNUMS /   NCOR,        IWCAP,       IPRE,        ICOR
      COMMON / SWNUMS /   IBOT,        ICUR,        IDBR
      COMMON / SWNUMS /   IDIF,        IINC,        ITRIAD,      IREFR
      COMMON / SWNUMS /   ISURF,       ITRSY,       IWIND,       IQUAD
      COMMON / SWNUMS /                ITERMX,      NSTATC,      NSTATM   40.00 40.22
      COMMON / SWNUMS /   U10,         U10C
      COMMON / SWNUMS /   WDIP,        WDIC,        DEPMIN,      PWCAP
      COMMON / SWNUMS /   PBOT,        PTRIAD,      PNUMS,       PSURF
      COMMON / SWNUMS /   PWIND,       SY0,         SIGMAG,      ITFRE
      COMMON / SWNUMS /   NUMOBS,      LSETUP,      BNDCHK,      HSRERR   40.00
      COMMON / SWNUMS /   FSHAPE,      DSHAPE,      PSHAPE,      SPPARM   40.00
      COMMON / SWNUMS /   BNAUT,       ONED,        PQUAD,       BRESCL   40.00
      COMMON / SWNUMS /   IGEN,                     PSETUP                40.03 40.22
      COMMON / SWNUMS /   ACUPDA,      MXITST,      MXITNS                40.07
      COMMON / SWNUMS /   NCOMPT,      RCOMPT
!
      COMMON / TP_SWNUMS /   ICMAX, CSETUP                                 40.22
!$OMP THREADPRIVATE(/TP_SWNUMS/)                                           40.22
