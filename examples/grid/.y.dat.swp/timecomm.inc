!                     TIME RELATED COMMON BLOCKS
!
!     Contents of this file
!
!     COMMON / TIMCOM /  *** Time related variables for the computation
!     COMMON / TIMRED /  *** Time related variables for nested runs
!     COMMON / TIMCPU /  *** contains the cpu and wall-clock times
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
!     40.23 Marcel Zijlema
!
!  1. Updates
!
!     30.72, Nov. 97: Created to make INCLUDE file from the various
!                     appearances of these COMMONS in SWAN
!     40.00, July 98: common area TIMFIL moved to swcomm2
!                     reason: it is related to input grids
!                     REFTIM mover to Ocean Pack common OCPCOMM1.INC
!     40.23, Aug. 02: Introduction new common block TIMCPU containing
!                     cpu and wall-clock times spent in different
!                     sections of SWAN
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
!     COMMON / TIMCOM /  *** Time related variables for the computation
!
!     Common variables
!
!           [????] default value
!
! TINIC     [    ] Start time and date of the computation (in seconds since
!                  the reference day (REFDAY))
!                  =tbegc; set by command COMP [tbegc] ...
! DT        [    ] Time step of the computation (in seconds)
!                  =deltc; set by command COMP ... [deltc] ...
!                  =deltc*60; set by command COMP ... [deltc] MI ...
!                  =deltc*60*60; set by command COMP ... [deltc] HR ...
!                  =deltc*60*60*24; set by command COMP ... [deltc] DA ...
! TFINC     [    ] End time and date of the computation (in seconds since
!                  the reference day (REFDAY))
! TIMCO     [    ] Time and date of the computation during the simulation (in
!                  seconds since the reference day (REFDAY))
!
      REAL                TINIC,       DT,          TFINC,       TIMCO
!
      COMMON / TIMCOM /   TINIC,       DT,          TFINC,       TIMCO
!
!-------------------------------------------------------------------------------
!
!     COMMON / TIMFIL /  *** Time related variables for the grids
!
!     Common variables
!
!           [????] default value
!
! INTECU    [    ] Time step between non-stationary input conditions
!                  for currents (in seconds)
!                  =deltinp; set by command INP CUR ... [deltinp] S/MI/HR/DA ...
! INTEFR    [    ] Time step between non-stationary input conditions
!                  for bottom friction (in seconds)
!                  =deltinp; set by command INP FR ... [deltinp] S/MI/HR/DA ...
! INTEWI    [    ] Time step between non-stationary input conditions
!                  for wind (in seconds)
!                  =deltinp; set by command INP WI ... [deltinp] S/MI/HR/DA ...
! INTEWL    [    ] Time step between non-stationary input conditions
!                  for water levels (in seconds)
!                  =deltinp; set by command INP WLEV ... [deltinp] S/MI/HR/DA ...
!
!
!           [????] default value
!
! TBEGCU    [    ] Start time for the non-stationary input conditions
!                  for currents (in seconds since the reference day (REFDAY))
!                  =tbeginp; set by command INP CUR ... [tbeginp] ...
! TBEGFR    [    ] Start time for the non-stationary input conditions
!                  for bottom friction (in seconds since the reference day (REFDAY))
!                  =tbeginp; set by command INP FR ... [tbeginp] ...
! TBEGWI    [    ] Start time for the non-stationary input conditions
!                  for wind (in seconds since the reference day (REFDAY))
!                  =tbeginp; set by command INP WI ... [tbeginp] ...
! TBEGWL    [    ] Start time for the non-stationary input conditions
!                  for water levels (in seconds since the reference day (REFDAY))
!                  =tbeginp; set by command INP WLEV ... [tbeginp] ...
! TENDCU    [    ] End time for the non-stationary input conditions
!                  for currents (in seconds since the reference day (REFDAY))
!                  =tendinp; set by command INP CUR ... [tendinp] ...
! TENDFR    [    ] End time for the non-stationary input conditions
!                  for bottom friction (in seconds since the reference day (REFDAY))
!                  =tendinp; set by command INP FR ... [tendinp] ...
! TENDWI    [    ] End time for the non-stationary input conditions
!                  for wind (in seconds since the reference day (REFDAY))
!                  =tendinp; set by command INP WI ... [tendinp] ...
! TENDWL    [    ] End time for the non-stationary input conditions
!                  for water levels (in seconds since the reference day (REFDAY))
!                  =tendinp; set by command INP WLEV ... [tendinp] ...
! TIMCU     [    ] Last time that non-stationary input conditions has been read
!                  for currents (in seconds since the reference day (REFDAY))
! TIMFR     [    ] Last time that non-stationary input conditions has been read
!                  for bottom friction (in seconds since the reference day (REFDAY))
! TIMWI     [    ] Last time that non-stationary input conditions has been read
!                  for wind (in seconds since the reference day (REFDAY))
! TIMWL     [    ] Last time that non-stationary input conditions has been read
!                  for water levels (in seconds since the reference day (REFDAY))
!
!
!
!-------------------------------------------------------------------------------
!
!     COMMON / TIMRED /  *** Time related variables for nested runs
!
!     Common variables
!
!           [????] default value
!
! IFACMX    [   1] ????
!                  =ifacmx; set by command BOU (STAT) NE ... [ifacmx] ...
!                  (NOT documented)
! IFACMY    [   1] ????
!                  =ifacmy; set by command BOU (STAT) NE ... [ifacmy]
!                  (NOT documented)
!
      INTEGER             IFACMX,      IFACMY
!
!           [????] default value
!
! BEGBOU    [    ] Start time for the non-stationary boundary conditions
!                  (in seconds since the reference day (REFDAY)) in the case
!                  of nested runs. Read from the nest file
! TIMERB    [    ] Last time that non-stationary boundary conditions has been read
!                  (in seconds since the reference day (REFDAY))) in the case
!                  of nested runs. Read from the nest file
! TINTBO    [    ] Time step between non-stationary boundary conditions (in seconds)
!                  in the case of nested runs. Read from the nest file
!
      REAL                BEGBOU,      TIMERB,      TINTBO
!
      COMMON / TIMRED /   BEGBOU,      TINTBO,      TIMERB,      IFACMX
      COMMON / TIMRED /   IFACMY
!
!-------------------------------------------------------------------------------
!
!     COMMON / TIMCPU /  *** contains the cpu and wall-clock times
!
!     Common variables
!
! DCUMTM    [ 600] cumulative time; columns 1,2: cpu-time, wall-clock time
! LASTTM    [   1] last occupied position in LISTTM, 0 if all positions in
!                  LISTTM are free
! LISTTM    [  10] list of section numbers for all running/active timers
!                  in array TIMERS. A value of -1 signals that the
!                  corresponding timer is not running
! MXTIMR    [   1] maximum number of simultaneous timings, dimension of
!                  work-arrays LISTTM, TIMERS
! NCUMTM    [ 300] for each section of the application the number of timings
!                  that contributed to time in DCPUTM
! NSECTM    [   1] number of sections that may be defined for an application
! TIMERS    [ 600] start-time of active timers; columns 1,2: cpu-time,
!                  wall-clock time
!
      INTEGER NSECTM, MXTIMR
      PARAMETER (NSECTM=300, MXTIMR=10)
!
      INTEGER          NCUMTM, LASTTM, LISTTM
      DOUBLE PRECISION DCUMTM, TIMERS
!
      COMMON /TIMCPU/ DCUMTM(NSECTM,2), TIMERS(MXTIMR,2),
     &                NCUMTM(NSECTM), LISTTM(MXTIMR), LASTTM
!$OMP THREADPRIVATE(/TIMCPU/)                                             40.31
      SAVE /TIMCPU/
!
