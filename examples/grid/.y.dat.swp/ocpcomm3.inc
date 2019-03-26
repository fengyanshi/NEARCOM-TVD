!              OCEAN PACK COMMON BLOCKS, file 3 of 4
!
!     Contents of this file
!
!     COMMON / OUTPDA / *** data for output, mainly plotting
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
!        [???] default value
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
!  COMMON / OUTPDA / *** data for output, mainly plotting ***
!
!     Common variables
!
!        [??????] default value
!
! MXQ    [CALCUL] number of grid points of the output frame in X-direction
! MYQ    [CALCUL] number of grid points of the output frame in Y-direction
!
      INTEGER             MXQ,         MYQ
!
! DXQ    [      ] mesh size of the output frame in X-direction
!                 =0.01; if MXQ=1 ????
!                 =XQLEN/(MXQ-1); if MXQ>1
! DYQ    [      ] mesh size of the output frame in Y-direction
!                 =0.01; if MYQ=1 ????
!                 =YQLEN/(MYQ-1); if MYQ>1
! VERNUM [ 30.72] version number of SWAN
!
      REAL                DXQ,         DYQ,         VERNUM
!
      COMMON / OUTPDA /   VERNUM
      COMMON / OUTPDA /   MXQ,         MYQ,         DXQ,         DYQ
