!              OCEAN PACK COMMON BLOCKS, file 1 of 4
!
!     Contents of this file
!
!     COMMON / LEESDA / *** char. data used by the command reading system
!     COMMON / LEESDN / *** num. data used by the command reading system
!     COMMON / REFTIM / *** Origin for day and time
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
!
!  1. Updates
!
!     30.72, Sept 97: Created from the old ocpcomm.ftn file
!     40.00, Apr. 99: introduced constant LINELN
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
!  LINELN [ 120]  max length of input lines (in command file)                40.00
!
      INTEGER    LINELN
      PARAMETER  (LINELN=120)
!
!  COMMON / LEESDA / *** char. data used by the command reading system ***
!
!     Common variables
!
!        [??????] default value
!
! BLANK  ['    '] blank string
! COMID  [   '$'] character which distinguishes comments in the command input
! ELTEXT [      ] contents of the last string read by reading system
! ELTYPE [      ] type of the element last read by reading system
!                 ='CHAR'; last read data element is the string in ELTEXT
!                 ='EMPT;  empty data field
!                 ='EOF';  end of file has been reached
!                 ='EOR';  end of repeat has been reached
!                 ='ERR';  incorrect data field was encountered
!                 ='INT';  last read data element is the integer in ELINT
!                 ='KEY';  last read data element is the keyword in KEYWRD
!                 ='OTHR'; other
!                 ='REAL'; last read data element is the real in ELREAL
!                 ='USED'; last read data element is processed, new can be read
! KAART  [      ] contents of the input line last read by reading system
! KAR    [      ] character last read by reading system
!                 =COMID; begin or end of comment
!                 =TABC; data element separation mark
!                 =' '; data element separation mark
!                 ='&'; continuation mark
!                 ='('; begin of data group
!                 ='*'; repetition mark
!                 =','; data element seperation mark (only for numbers)
!                 ='/'; end of repetition mark
!                 =':'; assignment mark
!                 =';'; end of record or input line
!                 ='='; assignment mark
!                 ='@'; end of file mark
!                 ='_'; continuation mark
!                 other: letter or digit to be processed by reading system
! KEYWRD [      ] contents of the last keyword read by reading system
! TABC   [CALCUL] =CHAR(9); tabular character
!
      CHARACTER*4         BLANK
      CHARACTER           COMID
      CHARACTER*(LINELN)  ELTEXT
      CHARACTER*4         ELTYPE
      CHARACTER*(LINELN)  KAART
      CHARACTER           KAR
      CHARACTER*8         KEYWRD
      CHARACTER           TABC
!
      COMMON / LEESDA /   ELTYPE,      ELTEXT,      KAART,       KAR
      COMMON / LEESDA /   KEYWRD,      BLANK,       TABC,        COMID
!
!-------------------------------------------------------------------------------
!
!  COMMON / LEESDN / *** num. data used by the command reading system ***
!
!     Common variables
!
!        [??????] default value
!
! ELREAL [      ] last element read from user command, when real or double
!
      DOUBLE PRECISION    ELREAL
!
! ELINT  [      ] last element read from user command, when integer
! KARNR  [      ] position on the input line of character last processed
!                 by the reading system,
!                 =0; no characters read yet
!                 =81; next input line has to be read to the common KAART first
! LENCST [      ] length of the string stored in ELTEXT
!
      INTEGER             ELINT,       KARNR,       LENCST
!
! CHGVAL [      ] whether last read value is different from a given value for
!                 subroutines INREAL, ININTG, INCSTR, INCTIM
!
      LOGICAL             CHGVAL
!
      COMMON / LEESDN /   ELREAL,      ELINT,       KARNR,       CHGVAL
      COMMON / LEESDN /   LENCST
!
!-------------------------------------------------------------------------------
!
!     COMMON / REFTIM /  *** Origin for day and time
!
!     Common variables
!
!           [????] default value
!
! REFDAY    [    ] Day number of the reference day. The first day entered is used
!                  as reference day, the reference time is 0:00 of the reference day.
!
      INTEGER             REFDAY
!
      COMMON / REFTIM /   REFDAY
!
