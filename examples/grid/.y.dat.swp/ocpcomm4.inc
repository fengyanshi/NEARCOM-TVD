!              OCEAN PACK COMMON BLOCKS, file 4 of 4
!
!     Contents of this file
!
!     COMMON / REFNRS / *** file unit reference numbers
!     COMMON / TESTDA / *** test parameters
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
!     30.72, Sep 97: Created from the old ocpcomm.ftn file
!     40.23, Feb 03: two parameters added (INAN, RNAN)
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN.
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
! INAN      [   1] integer representing not a number
! RNAN      [   1] real representing not a number
!
      INTEGER INAN
      REAL    RNAN
      PARAMETER (INAN=-1073750760,
     &           RNAN=-1.07374515E+09)
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
!  COMMON / REFNRS / *** file unit reference numbers ***
!
!     Common variables
!
!        [???] default value
!
! EXPORT [   ] not used
! FUNHI  [CAL] highest free unit number
!              =IUNMAX
! FUNLO  [ 21] lowest free unit number
! HIOPEN [   ] highest unit number of an open file
! IMPORT [   ] not used
! INPUTF [  3] unit number for the file with command input ('INPUT')
!              set by file SWANINIT
! ITMOPT [  1] time coding option
!              =1; ????
!              =3; ????
!              set by file SWANINIT
! IUNMIN [  0] minimum unit number
! IUNMAX [ 99] maximum unit number
!              set by file SWANINIT
! PRINTF [  4] unit number for the file with standard output ('PRINT')
!              set by file SWANINIT
! PRTEST [CAL] unit number for the print file containing test output
!              =PRINTF
!              set by file SWANINIT
! SCREEN [  6] unit number for the screen
!              (is for batch-oriented systems equal to PRINTF)
!
      INTEGER             EXPORT,      FUNHI,       FUNLO,       HIOPEN
      INTEGER             IMPORT,      INPUTF,      ITMOPT,      IUNMAX
      INTEGER             IUNMIN,      PRINTF,      PRTEST,      SCREEN
!
      COMMON / REFNRS /   PRINTF,      INPUTF,      IUNMIN,      IUNMAX
      COMMON / REFNRS /   FUNLO,       FUNHI,       SCREEN,      PRTEST
      COMMON / REFNRS /   IMPORT,      EXPORT,      HIOPEN,      ITMOPT
!
!-------------------------------------------------------------------------------
!
!  COMMON / TESTDA / *** test parameters ***
!
!     Common variables
!
!        [???] default value
!
! ITEST  [  0] indicates the amount of test output requested
!              =30; for command TEST
!              =itest; set by command TEST [itest] [itrace]
! ITRACE [  0] a message is printed up to ITRACE times
!              =itrace; set by command TEST [itest] [itrace]
! LEVERR [  0] severity of the errors encountered.
! MAXERR [  1] maximum severity of errors allowed, if larger no computation
!              =1; warnings
!              =2; errors
!              =3; severe errors
!              =4; terminating errors
!              =maxerr; set by command SET ... [maxerr] ...
!
      INTEGER             ITEST,       ITRACE,      LEVERR,      MAXERR
!
! LTRACE [.F.] indicates whether to call STRACE
!              =.T.; when ITRACE>0
!
      LOGICAL             LTRACE
!
      COMMON / TESTDA /   ITEST,       ITRACE,      LTRACE,      LEVERR
      COMMON / TESTDA /   MAXERR
!
