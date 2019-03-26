!              OCEAN PACK COMMON BLOCKS, file 2 of 4
!
!     Contents of this file
!
!     COMMON / NAMES /  *** names and other character strings
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
!     40.13 Nico Booij
!
!  1. Updates
!
!     30.72, Sept 97: Created from the old ocpcomm.ftn file
!     40.00, Mar. 99: DIRCH1 and DIRCH2 added
!     40.03, May  00: VERTXT (version+bugfixes) added
!     40.13, Oct. 01: length of filenames increased
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
!  LENFNM [  80]  length of file names                                    40.03
!
      INTEGER    LENFNM
      PARAMETER  (LENFNM=80)                                              40.13
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
!  COMMON / NAMES / *** names and other character strings ***
!
!     Common variables (These are used in both SWAN and Ocean Pack)
!
!        [??????] default value
!
! DIRCH1 [\     ] directory separation character as appears in input file
! DIRCH2 [\     ] directory separation character replacing DIRCH1         40.00
! FILEA  [      ] not used
! FILEB  [      ] not used
! FILENM [      ] file name of the file currently used for I/O???
! INST   ['Delft University of Technology'] name of the institute
!                 Can be changed in the file SWANINIT
! PROJID ['SWAN'] acronym of the project for which the computation is taking place
!                 ='NAME'; set by command PROJ 'NAME' ...
! PROJNR [CALCUL] =BLANK; run number for the computation
!                 ='NR'; set by command PROJ ... 'NR' ...
! PROJT1 [CALCUL] =BLANK; 1st line of the project title
!                 ='title1'; set by command PROJ ... 'title1' ...
! PROJT2 [CALCUL] =BLANK; 2nd line of the project title
!                 ='title2'; set by command PROJ ... 'title2' ...
! PROJT3 [CALCUL] =BLANK; 3rd line of the project title
!                 ='title3'; set by command PROJ ... 'title3'
! PTITLE [      ] not used
! VERTXT [calcul] program version, character representation               40.03
!
      CHARACTER (LEN=1)      :: DIRCH1,      DIRCH2                             40.00
      CHARACTER (LEN=LENFNM) :: FILEA,       FILEB,       FILENM                40.03
      CHARACTER (LEN=40)     :: INST
      CHARACTER (LEN=16)     :: PROJID
      CHARACTER (LEN=4)      :: PROJNR
      CHARACTER (LEN=72)     :: PROJT1,      PROJT2,      PROJT3
      CHARACTER (LEN=36)     :: PTITLE
      CHARACTER (LEN=20)     :: VERTXT                                          40.03
!
      COMMON / NAMES  /   INST,        PROJID,      PROJNR,      PROJT1
      COMMON / NAMES  /   PROJT2,      PROJT3,      PTITLE,      FILENM
      COMMON / NAMES  /   FILEA,       FILEB,       DIRCH1,      DIRCH2   40.00
      COMMON / NAMES  /   VERTXT                                          40.03
