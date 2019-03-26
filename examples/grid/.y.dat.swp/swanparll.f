!
!     SWAN - routines for distributed-memory approach based on MPI
!
!  Contents of this file
!
!     SWINITMPI
!     SWEXITMPI
!     SWSYNC
!     SWSENDNB
!     SWRECVNB
!     SWBROADC
!     SWGATHER
!     SWREDUCE
!     SWSTRIP
!     SWPARTIT
!     SWBLADM
!     SWDECOMP
!     SWEXCHG
!     SWRECVAC
!     SWSENDAC
!     SWCOLOUT
!     SWCOLTAB
!     SWCOLSPC
!     SWCOLBLK
!
!****************************************************************
!
      SUBROUTINE SWINITMPI
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!
!  2. Purpose
!
!     Join parallel application
!
!  3. Method
!
!     Start MPI and initialize common block SWPARALL
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_COMM_RANK    Get rank of processes in MPI communication context
!MPI!     MPI_COMM_SIZE    Get number of processes in MPI communication context
!MPI!     MPI_INIT         Enroll in MPI
!     MSGERR           Writes error message
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Start MPI and initialize common block SWPARALL
!
! 13. Source text
!
      LEVERR = 0
      MAXERR = 1
      ITRACE = 0

!MPI!     --- enroll in MPI
!MPI
!MPI      CALL MPI_INIT ( IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(1)(IF1:IL1)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- initialize common block SWPARALL

      INODE = 0
      NPROC = 1

!MPI!     --- get node number INODE
!MPI
!MPI      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, INODE, IERR )
      INODE = INODE + 1
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(1)(IF1:IL1)//NWL//
!MPI     &            ' Node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
!MPI
!MPI!     --- determine total number of processes
!MPI
!MPI      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, NPROC, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(1)(IF1:IL1)//NWL//
!MPI     &            ' Node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- determine whether this is a parallel run or not

      IF ( NPROC.GT.1 ) THEN
         PARLL = .TRUE.
      ELSE
         PARLL = .FALSE.
      END IF

!MPI!     --- define MPI constants for communication within SWAN
!MPI
!MPI      SWINT  = MPI_INTEGER
!MPI      SWREAL = MPI_REAL
!MPI      SWMAX  = MPI_MAX
!MPI      SWMIN  = MPI_MIN
!MPI      SWSUM  = MPI_SUM

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXITMPI
!
!****************************************************************
!
!MPI      USE MPI
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!
!  2. Purpose
!
!     Exit parallel application
!
!  3. Method
!
!     Wrapper for MPI_FINALIZE
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IERR    :   error value of MPI call
!     PARALMPI:   if true, parallel process is carried out with MPI
!
      INTEGER IERR
      LOGICAL PARALMPI
!
!  8. Subroutines used
!
!MPI!     MPI_ABORT        Abort MPI if severe error occurs
!MPI!     MPI_BARRIER      Blocks until all nodes have called this routine
!MPI!     MPI_INITIALIZED  Indicates whether MPI_Init has been called
!MPI!     MPI_FINALIZE     Cleans up the MPI state and exits
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     if MPI has been initialized
!MPI!        synchronize nodes
!MPI!        if severe error
!MPI!           abort MPI
!MPI!        else
!MPI!           close MPI
!MPI!
! 13. Source text
!
!MPI      CALL MPI_INITIALIZED ( PARALMPI, IERR )
!MPI      IF ( PARALMPI ) THEN
!MPI
!MPI         CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR )
!MPI
!MPI         IF ( LEVERR.GE.4 ) THEN
!MPI
!MPI!        --- in case of a severe error abort all MPI processes
!MPI
!MPI            CALL MPI_ABORT ( MPI_COMM_WORLD, LEVERR, IERR )
!MPI
!MPI         ELSE
!MPI
!MPI!        --- otherwise stop MPI operations on this computer
!MPI
!MPI            CALL MPI_FINALIZE ( IERR )
!MPI
!MPI         END IF
!MPI
!MPI      END IF
!MPI
      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSYNC
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!
!  2. Purpose
!
!     Synchronize nodes
!
!  3. Method
!
!     Wrapper for MPI_BARRIER
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BARRIER      Blocks until all nodes have called this routine
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Blocks until all nodes have called MPI_BARRIER routine.
!     In this way, all nodes are synchronized
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSYNC')

!MPI!     --- blocks until all nodes have called this routine
!MPI
!MPI      CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(1)(IF1:IL1)//NWL//
!MPI     &            ' Node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDNB ( IPTR, ILEN, ITYPE, IDEST, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Data is sent to a neighbour
!
!  3. Method
!
!     Wrapper for MPI_SEND
!
!  4. Argument variables
!
!     IDEST       rank of the destination process
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, IDEST, ITAG
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_SEND         Immediately sends the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Data is sent to a neighbour with command MPI_SEND
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI      CALL MPI_SEND ( IPTR, ILEN, ITYPE, IDEST-1,
!MPI     &                ITAG, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(1)(IF1:IL1)//NWL//
!MPI     &            ' Node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVNB ( IPTR, ILEN, ITYPE, ISOURCE, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Data is received from a neighbour
!
!  3. Method
!
!     Wrapper for MPI_RECV
!
!  4. Argument variables
!
!     ILEN        length of array to be received
!     IPTR        pointer to first element of array to be received
!     ISOURCE     rank of the source process
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, ISOURCE, ITAG
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!MPI!     ISTAT :     MPI status array
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
!MPI      INTEGER      ISTAT(MPI_STATUS_SIZE)
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_RECV         Immediately receives the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
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
!MPI!     Data is received from a neighbour with command MPI_RECV
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI      CALL MPI_RECV ( IPTR, ILEN, ITYPE, ISOURCE-1, ITAG,
!MPI     &                MPI_COMM_WORLD, ISTAT, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(1)(IF1:IL1)//NWL//
!MPI     &            ' Node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBROADC ( IPTR, ILEN, ITYPE )
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Broadcasts data from the master to all other processes
!
!  3. Method
!
!     Wrapper for MPI_BCAST
!
!  4. Argument variables
!
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BCAST        Broadcasts a message from the master
!MPI!                      to all other processes of the group
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     SWTSTA           Start timing for a section of code
!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SNEXTI
!     FLFILE
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Broadcasts data from the master to all other nodes
!MPI!     with command MPI_BCAST
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBROADC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      CALL SWTSTA(201)
!MPI      CALL MPI_BCAST ( IPTR, ILEN, ITYPE, MASTER-1,
!MPI     &                 MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      CALL SWTSTO(201)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWGATHER ( IOPTR, IOLEN, IIPTR, IILEN, ITYPE )
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Gathers different amounts of data from each processor
!     to the master
!
!  3. Method
!
!     Wrapper for MPI_GATHERV
!
!  4. Argument variables
!
!     IILEN       length of input array
!     IIPTR       pointer to first element of input array (local)
!     IOLEN       length of output array
!     IOPTR       pointer to first element of output array (global)
!     ITYPE       type of data
!
      INTEGER IILEN, IIPTR, IOLEN, IOPTR, ITYPE
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     ICOUNT:     array specifying array size of data received
!                 from each processor
!     IDSPLC:     array specifying the starting address of the
!                 incoming data from each processor, relative
!                 to the global array
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER                 I, IENT, IERR, IF, IL
      INTEGER, ALLOCATABLE :: ICOUNT(:), IDSPLC(:)
      CHARACTER*20            INTSTR, CHARS
      CHARACTER*80            MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_GATHER       Gathers data from all nodes to the master
!MPI!     MPI_GATHERV      Gathers different amounts of data from
!MPI!                      all nodes to the master
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLLECT
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!MPI!     gather the array sizes to the master
!MPI!
!     check whether enough space has been allocated
!     for gathered data
!
!     calculate starting address of each local array
!     with respect to the global array
!
!MPI!     gather different amounts of data from each processor
!MPI!     to the master
!MPI!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWGATHER')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF (INODE.EQ.MASTER) THEN
         ALLOCATE(ICOUNT(0:NPROC-1))
         ALLOCATE(IDSPLC(0:NPROC-1))
      END IF

!MPI!     --- gather the array sizes to the master
!MPI
!MPI      CALL MPI_GATHER( IILEN, 1, SWINT, ICOUNT, 1, SWINT,
!MPI     &                 MASTER-1, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- check whether enough space has been allocated
!         for gathered data

      IF (INODE.EQ.MASTER) THEN
         IF ( SUM(ICOUNT).GT.IOLEN ) THEN
            CALL MSGERR(4,
     &                  'Not enough space allocated for gathered data')
            RETURN
         END IF
      END IF

!     --- calculate starting address of each local array
!         with respect to the global array

      IF (INODE.EQ.MASTER) THEN
         IDSPLC(0) = 0
         DO I = 1, NPROC-1
            IDSPLC(I) = ICOUNT(I-1) + IDSPLC(I-1)
         END DO
      END IF

!MPI!     --- gather different amounts of data from each processor
!MPI!         to the master
!MPI
!MPI      CALL MPI_GATHERV( IIPTR, IILEN, ITYPE, IOPTR, ICOUNT, IDSPLC,
!MPI     &                  ITYPE, MASTER-1, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      IF (INODE.EQ.MASTER) DEALLOCATE(ICOUNT,IDSPLC)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCE ( IPTR, ILEN, ITYPE, ITYPRD )
!
!****************************************************************
!
!MPI      USE MPI
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on
!     array (IPTR) of type ITYPE to collect values from
!     all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     ILEN        length of array to be collect
!     IPTR        pointer to first element of array to be collect
!     ITYPE       type of data
!     ITYPRD      type of reduction
!
      INTEGER IPTR, ILEN, ITYPE, ITYPRD
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     IOPTR :     pointer to first element of output array
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL, IOPTR
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_ALLREDUCE    Combines values from all processes and
!MPI!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     SWCOPI           Makes copy of an integer array
!     SWCOPR           Makes copy of a real array
!     SWTSTA           Start timing for a section of code
!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOMP
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
!MPI!     Performs a global reduction of data across all nodes
!MPI!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCE')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      CALL SWTSTA(202)
!MPI      CALL MPI_ALLREDUCE ( IPTR, IOPTR, ILEN, ITYPE,
!MPI     &                     ITYPRD, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error'//NWL//
!MPI     &            ' Return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      IF ( ITYPE.EQ.SWINT ) THEN
         CALL SWCOPI ( IOPTR, IPTR, ILEN )
      ELSE IF ( ITYPE.EQ.SWREAL ) THEN
         CALL SWCOPR ( IOPTR, IPTR, ILEN )
      END IF
      CALL SWTSTO(202)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSTRIP ( IPOWN, IDIR, NPART, IWORK, MXC, MYC )
!
!****************************************************************
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Performs a stripwise partitioning with straight interfaces
!
!  3. Method
!
!     Each active point in a row/column will be assign to a part
!     according to its number and size (stored in IWORK).
!     The remaining points in the row/column will be assign
!     to the same part.
!
!  4. Argument variables
!
!     IDIR        direction of cutting
!                 1 = row
!                 2 = column
!     IPOWN       array giving the subdomain number of each gridpoint
!     IWORK       work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!     NPART       number of parts to be created
!
      INTEGER IDIR, MXC, MYC, NPART
      INTEGER IPOWN(*), IWORK(2,*)
!
!  6. Local variables
!
!     IC    :     index of (IX,IY)-point
!     ICC   :     index of (IX,IY)-point
!     IENT  :     number of entries
!     INCX  :     increment for adressing: 1 for x-dir, MXC for y-dir
!     INCY  :     increment for adressing: MXC for x-dir, 1 for y-dir
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IYY   :     index in y-direction
!     IPART :     a part counter
!     MXCI  :     maximum counter of gridpoints in x/y-direction
!     MYCI  :     maximum counter of gridpoints in y/x-direction
!     NCURPT:     number of currently assigned points to a created part
!     NPREM :     number of remaining points in a row/column
!
      INTEGER IC, ICC, IENT, INCX, INCY, IX, IY, IYY, IPART,
     &        MXCI, MYCI, NCURPT, NPREM
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWPARTIT
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
!     depending on cutting direction, determine indirect addressing
!     create first empty part
!     for all active points do
!         assign this point to the created part
!         if size of created part has been reached
!            determine remaining active points in the current column
!            if no remaining points, create next empty part
!            else remaining points belong to the current part
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSTRIP')

!     --- depending on cutting direction, determine indirect addressing
!         for array IPOWN

      IF ( IDIR.EQ.1 ) THEN
         MXCI = MYC
         MYCI = MXC
         INCX = MXC
         INCY = 1
      ELSE IF ( IDIR.EQ.2 ) THEN
         MXCI = MXC
         MYCI = MYC
         INCX = 1
         INCY = MXC
      END IF

!     --- create first empty part

      IPART  = 1
      NCURPT = 0

!     --- for all active points do

      DO IX = 1, MXCI
         DO IY = 1, MYCI

            IC = IX*INCX + IY*INCY - MXC

            IF ( IPOWN(IC).EQ.1 ) THEN

!              --- assign this point to the created part

               IPOWN(IC) = IWORK(1,IPART)
               NCURPT    = NCURPT + 1

!              --- if size of created part has been reached

               IF ( NCURPT.GE.IWORK(2,IPART) ) THEN

!                 --- determine remaining active points in the
!                     current column

                  NPREM = 0
                  DO IYY = IY+1, MYCI
                     ICC = IX*INCX + IYY*INCY - MXC
                     IF (IPOWN(ICC).EQ.1) NPREM = NPREM +1
                  END DO

                  IF ( NPREM.EQ.0 ) THEN

!                    --- if no remaining points, create next empty part

                     IPART  = IPART + 1
                     NCURPT = 0

                  ELSE

!                    --- else remaining points belong to the current part

                     IWORK(2,IPART  ) = IWORK(2,IPART  ) + NPREM
                     IWORK(2,IPART+1) = IWORK(2,IPART+1) - NPREM

                  END IF

               END IF

            END IF

         END DO
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWPARTIT ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Carries out the partitioning of the SWAN computational grid
!
!  3. Method
!
!     Based on stripwise partitioning
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  6. Local variables
!
!     I     :     loop counter
!     IDIR  :     direction of cutting
!                 1 = row
!                 2 = column
!     IENT  :     number of entries
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IWORK :     work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     NACTP :     total number of active gridpoints
!     NPCUM :     cumulative number of gridpoints
!
      INTEGER I, IDIR, IENT, IX, IY, NACTP, NPCUM
      INTEGER IWORK(2,NPROC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWSTRIP          Performs a stripwise partitioning with straight
!                      interfaces
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWDECOMP
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
!     if not parallel, return
!     determine direction of cutting
!     determine number of active points
!     determine numbers and sizes of parts to be created
!     partition grid
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWPARTIT')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- determine direction of cutting

      IF ( MXC.GT.MYC ) THEN
         IDIR = 2
      ELSE
         IDIR = 1
      END IF

!     --- determine number of active points and
!         set IPOWN to 1 in these points

      NACTP = 0
      DO IX = 1, MXC
         DO IY = 1, MYC
            IF ( KGRPGL(IX,IY).NE.1 ) THEN
               IPOWN(IX,IY) = 1
               NACTP        = NACTP + 1
            END IF
         END DO
      END DO

!     --- determine numbers and sizes of parts to be created

      NPCUM = 0
      DO I = 1, NPROC
         IWORK(1,I) = I
         IWORK(2,I) = (NACTP*I)/NPROC - NPCUM
         NPCUM      = (NACTP*I)/NPROC
      END DO

!     --- partition grid

      CALL SWSTRIP ( IPOWN, IDIR, NPROC, IWORK, MXC, MYC )

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBLADM ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     For the present node, carries out the block administration
!     and determines array bounds with respect to global grid
!
!  3. Method
!
!     Based on domain decomposition, the interface sizes are
!     determined that is needed for the setup of block
!     administration stored as IBLKAD
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  5. Parameter variables
!
!     NWL  [   1] line feed
!DOS!     NWL  [   2] line feed + carriage return
!
      CHARACTER*1 NWL
!DOS      CHARACTER*2 NWL
      PARAMETER (NWL = CHAR(10))
!DOS      PARAMETER (NWL = CHAR(10)//CHAR(13))
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     IC    :     index of (IX,IY)-point
!     ICOFF :     offset of IC-index
!     ICRECV:     array containing positions of unknowns
!                 to be received from neighbour
!     ICSEND:     array containing positions of unknowns
!                 to be sent to neighbour
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     INB   :     neighbour counter
!     ISTART:     startaddress for each size interface in array IBLKAD
!     IX    :     index in x-direction
!     IXOFF :     offset in x-direction
!     IY    :     index in y-direction
!     IYOFF :     offset in y-direction
!     IWORK :     array used to determine interface sizes
!                   IWORK(1,i) = number of the i-th neighbour
!                   IWORK(2,i) = position of the i-th neighbour with
!                                respect to present subdomain
!                                (resp. top, bottom, right, left)
!                   IWORK(3,i) = size of interface to i-th neighbour
!     JOFFS :     offsets at which a point of a neigbhour domain can be found
!     MSGSTR:     string to pass message to call MSGERR
!     MXSIZ :     size of present subdomain in x-direction
!     MYSIZ :     size of present subdomain in y-direction
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!
      INTEGER      I, IC, ICOFF, IDOM, IENT, IF, IL, INB, ISTART,
     &             IX, IXOFF, IY, IYOFF, JOFFS(2,4), MXSIZ, MYSIZ,
     &             NNEIGH, NOVLU
      INTEGER      IWORK(3,NPROC),
     &             ICRECV(NPROC,MAX(MXC,MYC)),
     &             ICSEND(NPROC,MAX(MXC,MYC))
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWDECOMP
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
!     intialize offsets to be used in searching for interfaces
!     determine enclosing box of present subdomain
!     if subdomain appears to be empty
!        give warning and set empty bounding box
!     else
!        extend enclosing box to include halo area
!     localize global bounds in present subdomain
!     determine size of enclosing box
!     determine interface sizes:
!
!        loop over global grid
!           if point belongs to this part
!              for each of the four sizes
!                  if a neighbouring subdomain is found there
!                     find it in the list of neighbours
!                     if not yet in the list, add it
!                     store position of neighbour
!                     update number of overlapping unknowns
!
!     store block administration
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBLADM')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- intialize offsets to be used in searching for interfaces

      JOFFS = RESHAPE((/0,1,0,-1,1,0,-1,0/), (/2,4/))

!     --- determine enclosing box of present subdomain

      IF ( MXC.GT.MYC ) THEN
         MXF = MXC+1
         MXL = 0
      ELSE
         MYF = MYC+1
         MYL = 0
      END IF

      DO IX = 1, MXC
         DO IY = 1, MYC

            IF( IPOWN(IX,IY).EQ.INODE ) THEN

               MXF = MIN(IX,MXF)
               MYF = MIN(IY,MYF)
               MXL = MAX(IX,MXL)
               MYL = MAX(IY,MYL)

            END IF

         END DO
      END DO

!     --- if subdomain appears to be empty

      IF ( MXF.GT.MXL .OR. MYF.GT.MYL ) THEN

!        --- give warning and set empty bounding box

         CHARS = INTSTR(INODE)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'Empty subdomain is detected'//NWL//
     &            ' Node number is '//CHARS(IF:IL)
         CALL MSGERR ( 2, MSGSTR )

         MXF = 1
         MYF = 1
         MXL = 0
         MYL = 0

      ELSE

!        --- extend enclosing box to include halo area

         MXF = MAX(1  ,MXF-IHALOX)
         MYF = MAX(1  ,MYF-IHALOY)
         MXL = MIN(MXC,MXL+IHALOX)
         MYL = MIN(MYC,MYL+IHALOY)

      END IF

!     --- localize global bounds in present subdomain

      IF ( MXCGL.GT.MYCGL ) THEN
         LMXF = MXF.EQ.1     .AND. INODE.EQ.1
         LMXL = MXL.EQ.MXCGL .AND. INODE.EQ.NPROC
         LMYF = MYF.EQ.1
         LMYL = MYL.EQ.MYCGL
      ELSE
         LMXF = MXF.EQ.1
         LMXL = MXL.EQ.MXCGL
         LMYF = MYF.EQ.1     .AND. INODE.EQ.1
         LMYL = MYL.EQ.MYCGL .AND. INODE.EQ.NPROC
      END IF

!     --- determine size of enclosing box

      MXSIZ = MXL - MXF + 1
      MYSIZ = MYL - MYF + 1

      IWORK  = 0
      ICRECV = 0
      ICSEND = 0

!     --- determine interface sizes

      DO IX = 1, MXC
         DO IY = 1, MYC

!           --- if point belongs to this part

            IF ( IPOWN(IX,IY).EQ.INODE ) THEN

!              --- for each of the four sizes

               DO I = 1, 4

                  IXOFF = JOFFS(1,I)
                  IYOFF = JOFFS(2,I)

!                 --- if a neighbouring subdomain is found there

                  IF ( (IX+IXOFF).GT.0.AND.(IX+IXOFF).LE.MXC.AND.
     &                 (IY+IYOFF).GT.0.AND.(IY+IYOFF).LE.MYC ) THEN
                     IF ( IPOWN(IX+IXOFF,IY+IYOFF).NE.0.AND.
     &                    IPOWN(IX+IXOFF,IY+IYOFF).NE.INODE ) THEN
                        IC    = (IY      -MYF)*MXSIZ + (IX      -MXF+1)
                        ICOFF = (IY+IYOFF-MYF)*MXSIZ + (IX+IXOFF-MXF+1)

!                       --- find it in the list of neighbours

                        IDOM = IPOWN(IX+IXOFF,IY+IYOFF)

                        INB = 1
  100                   IF ( INB.LE.NPROC .AND.
     &                       IWORK(1,INB).NE.IDOM .AND.
     &                       IWORK(1,INB).NE.0 )  THEN
                           INB = INB + 1
                           GOTO 100
                        END IF

                        IF ( INB.GT.NPROC ) THEN
                          CALL MSGERR (4,'Found more neighbours than '//
     &                                 'subdomains in the partitioning')
                          RETURN
                        END IF

!                       --- if not yet in the list, add it
                        IF ( IWORK(1,INB).EQ.0 ) IWORK(1,INB) = IDOM
                   
!                    --- store position of neighbour with respect to
!                        present subdomain
                        IWORK(2,INB) = I

!                    --- update number of overlapping unknowns
                        IWORK(3,INB) = IWORK(3,INB) + 1
                        ICSEND(INB,IWORK(3,INB)) = IC
                        ICRECV(INB,IWORK(3,INB)) = ICOFF
                        END IF

                  END IF

               END DO

            END IF

         END DO
      END DO

!     --- store block administration

      NNEIGH    = COUNT(IWORK(1,:)>0)
      IBLKAD(1) = NNEIGH
      ISTART    = 3*NNEIGH+2
      DO INB = 1, NNEIGH
         IBLKAD(3*INB-1) = IWORK(1,INB)
         IBLKAD(3*INB  ) = IWORK(2,INB)
         IBLKAD(3*INB+1) = ISTART
         NOVLU           = IWORK(3,INB)
         IBLKAD(ISTART)  = NOVLU
         DO I = 1, NOVLU
            IBLKAD(ISTART      +I) = ICSEND(INB,I)
            IBLKAD(ISTART+NOVLU+I) = ICRECV(INB,I)
         END DO
         ISTART = ISTART + 2*NOVLU+1
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWDECOMP
!
!****************************************************************
!
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm3.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Carries out domain decomposition meant for
!     distributed-memory approach
!
!  3. Method
!
!     First, carry out the partitioning of the
!     SWAN computational grid and then do the
!     block administration
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!     IX    :     loop counter
!     IY    :     loop counter
!     IPOWN :     array giving the subdomain number of each gridpoint
!
      INTEGER IENT, IX, IY
      INTEGER, ALLOCATABLE :: IPOWN(:,:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWBLADM          Carries out the block administration
!     SWPARTIT         Carries out the partitioning of the SWAN
!                      computational grid
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
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
!     allocate and initialize array for block administration
!     store the original values of MXC, MYC and MCGRD
!     if not parallel, return
!     carry out the partitioning of computational grid
!     carry out the block administration
!     compute MXC, MYC and MCGRD for each subdomain
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWDECOMP')

!     --- allocate and initialize array for block administration

      IF (.NOT.ALLOCATED(IBLKAD)) ALLOCATE(IBLKAD(41+20*MAX(MXC,MYC)))
      IBLKAD = 0

!     --- store the original values of MXC, MYC and MCGRD of
!         global computational grid

      MXCGL   = MXC
      MYCGL   = MYC
      MCGRDGL = MCGRD
      MXF     = 1
      MXL     = MXC
      MYF     = 1
      MYL     = MYC
      LMXF    = MXF.EQ.1
      LMXL    = MXL.EQ.MXCGL
      LMYF    = MYF.EQ.1
      LMYL    = MYL.EQ.MYCGL

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- carry out the partitioning of the SWAN
!         computational grid

      ALLOCATE(IPOWN(MXC,MYC))
      IPOWN = 0
      CALL SWPARTIT( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- carry out the block administration

      CALL SWBLADM( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- compute MXC, MYC and MCGRD for each subdomain

      MXC = MXL - MXF + 1
      MYC = MYL - MYF + 1

      MCGRD = 1
      DO IX = MXF, MXL
         DO IY = MYF, MYL
            IF ( KGRPGL(IX,IY).NE.1 ) MCGRD = MCGRD + 1
         END DO
      END DO

      DEALLOCATE(IPOWN)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXCHG ( FIELD, KGRPNT )
!
!****************************************************************
!
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm3.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Updates geographical field array through exchanging
!     values between neighbouring subdomains
!
!  3. Method
!
!     Made use of MPI by means of SWSENDNB and SWRECVNB
!     and also block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     FIELD       geographical field array for which 'halo' values must
!                 be copied from neighbouring subdomains
!     KGRPNT      indirect addressing for grid points
!
      INTEGER KGRPNT(MXC*MYC)
      REAL    FIELD(MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     ISTART:     pointer in array IBLKAD
!     ITAG  :     message tag for sending and receiving
!     K     :     loop counter
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     WORK  :     work array to store data to be sent to or
!                 received from neighbour
!
      INTEGER IDOM, IENT, INB, ISTART, ITAG, K, NNEIGH, NOVLU
      REAL    WORK(MAX(MXC,MYC))
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!     SWSENDNB         Data is sent to a neighbour
!     SWTSTA           Start timing for a section of code
!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     ---
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
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get subdomain number, pointer and size
!        store data to be sent in array WORK
!        send array WORK
!
!     for all neighbouring subdomains do
!        get subdomain number, pointer and size
!        receive next array and store in WORK
!        store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWEXCHG')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      CALL SWTSTA(203)

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get subdomain number, pointer and size

         IDOM   = IBLKAD(3*INB-1)
         ISTART = IBLKAD(3*INB+1)
         NOVLU  = IBLKAD(ISTART)

!        --- store data to be sent in array WORK

         DO K = 1, NOVLU
            WORK(K) = FIELD(KGRPNT(IBLKAD(ISTART+K)))
         END DO

!        --- send array WORK

         ITAG = 2
         CALL SWSENDNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
         IF (STPNOW()) RETURN

      END DO

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get subdomain number, pointer and size

         IDOM   = IBLKAD(3*INB-1)
         ISTART = IBLKAD(3*INB+1)
         NOVLU  = IBLKAD(ISTART)

!        --- receive next array and store in WORK

         ITAG  = 2
         CALL SWRECVNB ( WORK, NOVLU, SWREAL, IDOM, ITAG )
         IF (STPNOW()) RETURN

!        --- store the received data

         DO K = 1, NOVLU
            FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+K))) = WORK(K)
         END DO

      END DO

      CALL SWTSTO(203)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVAC ( AC2, IS, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm3.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Receives action density from neighbouring subdomains
!     depending on sweep direction
!
!  3. Method
!
!     Use of SWRECVNB and block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     AC2         action density
!     IS          start index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     SWPDIR      sweep direction
!
      INTEGER IS, J, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    AC2(MDC,MSC,MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     IPR   :     array containing positions of neighbours from
!                 which data is to be received
!     ITAG  :     message tag for sending and receiving
!     NNEIGH:     number of neighbouring subdomains
!     WORK  :     work array to store data to received from neighbour
!
      INTEGER IDOM, IENT, INB, IPNB, ITAG, NNEIGH
      INTEGER IPR(2,4)
      REAL    WORK(MDC,MSC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
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
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if position corresponds to sweep selection
!           get subdomain number
!           receive next array and store in WORK
!           store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVAC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IPR = RESHAPE((/2,4,2,3,1,3,1,4/), (/2,4/))

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if position corresponds to sweep selection

         IF ( IPNB.EQ.IPR(1,SWPDIR) .OR. IPNB.EQ.IPR(2,SWPDIR) ) THEN

!           --- get subdomain number

            IDOM   = IBLKAD(3*INB-1)

!           --- receive next array and store in WORK

            ITAG  = 2
            CALL SWRECVNB ( WORK, MDC*MSC, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

!           --- store the received data

            AC2(:,:,KGRPNT(IS,J)) = WORK(:,:)

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDAC ( AC2, IE, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm3.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Sends action density to neighbouring subdomains
!     depending on sweep direction
!
!  3. Method
!
!     Use of SWSENDNB and block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     AC2         action density
!     IE          end index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     SWPDIR      sweep direction
!
      INTEGER IE, J, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    AC2(MDC,MSC,MCGRD)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     IPS   :     array containing positions of neighbours to
!                 which data is to be sent
!     ITAG  :     message tag for sending and receiving
!     NNEIGH:     number of neighbouring subdomains
!     WORK  :     work array to store data to be sent to neighbour
!
      INTEGER IDOM, IENT, INB, IPNB, ITAG, NNEIGH
      INTEGER IPS(2,4)
      REAL    WORK(MDC,MSC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWSENDNB         Data is sent to a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
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
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if position corresponds to sweep selection
!           get subdomain number
!           store data to be sent in array WORK
!           send array WORK
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDAC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IPS = RESHAPE((/1,3,1,4,2,4,2,3/), (/2,4/))

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if position corresponds to sweep selection

         IF ( IPNB.EQ.IPS(1,SWPDIR) .OR. IPNB.EQ.IPS(2,SWPDIR) ) THEN

!           --- get subdomain number

            IDOM   = IBLKAD(3*INB-1)

!           --- store data to be sent in array WORK

            WORK(:,:) = AC2(:,:,KGRPNT(IE,J))

!           --- send array WORK

            ITAG = 2
            CALL SWSENDNB ( WORK, MDC*MSC, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLOUT ( IONOD )
!
!****************************************************************
!
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'timecomm.inc'
      INCLUDE 'ocpcomm4.inc'
	  INCLUDE 'swcomm1.inc'
      INCLUDE 'swcomm3.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May 03: New subroutine
!
!  2. Purpose
!
!     Collects output results
!
!  3. Method
!
!     Read individual process output files containing tables,
!     spectral outputs and block data and write them to
!     generic output files in appropriate manner
!
!  4. Argument variables
!
!     IONOD       array indicating in which subdomain and at what time
!                 output points are located
!
      INTEGER IONOD(MXMIP+1,MXOURQ,0:NPROC-1)
!
!  6. Local variables
!
!     CORQ  :     current item in list of request outputs
!     CUOPS :     current item in list of point sets
!     DIF   :     difference between end and actual times
!     DTTIWR:     to write time string
!     IC    :     loop variable
!     IENT  :     number of entries
!     IRQ   :     request number
!     IT    :     time step counter
!     IT0   :     integer indicating first step of simulation
!     IT1   :     integer indicating last step of simulation
!     IUNIT :     counter for file unit numbers
!     MIP   :     total number of output points
!     MXK   :     number of points in x-direction of output frame
!     MYK   :     number of points in y-direction of output frame
!     OPENED:     logical whether a file is open or not
!     PSTYPE:     type of point set
!     RTYPE :     type of request
!     SNAMPF:     name of plot frame
!     TNEXT :     time of next requested output
!
      INTEGER   IC, IENT, IRQ, IT, IT0, IT1, IUNIT, MIP, MXK, MYK
      REAL      DIF, TNEXT
      LOGICAL   OPENED
      CHARACTER PSTYPE*1, RTYPE*4, SNAMPF*8
      TYPE(ORQDAT), POINTER :: CORQ
      TYPE(OPSDAT), POINTER :: CUOPS
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWCOLBLK         Collects block output
!     SWCOLSPC         Collects spectral output
!     SWCOLTAB         Collects table ouput
!
      LOGICAL   STPNOW
      CHARACTER DTTIWR
!
!  9. Subroutines calling
!
!     SWMAIN
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
!     do for all COMPUTE commands
!        do for all time steps
!           do for all output requests
!              processing of output instructions necessary for collection
!              check time of output action
!              rewrite table output by means of collection of output
!              rewrite spectral output by means of collection of output
!              rewrite block output by means of collection of output
!     close all files and delete process files
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLOUT')

      IF ( NREOQ.EQ.0 ) RETURN

!     --- do for all COMPUTE commands
      DO IC = 1, NCOMPT
         NSTATC = NINT(RCOMPT(IC,1))
         IF ( NSTATC.EQ.1 ) THEN
            IT0 = 0
         ELSE
            IT0 = 1
         END IF
         IT1 = NINT(RCOMPT(IC,2))
         TFINC = RCOMPT(IC,3)
         TINIC = RCOMPT(IC,4)
         DT    = RCOMPT(IC,5)
         TIMCO = TINIC

!        --- do for all time steps
         DO IT = IT0, IT1
            IF (NSTATM.GT.0) CHTIME = DTTIWR(ITMOPT, TIMCO)
!           --- do for all output requests
            CORQ => FORQ
            DO 100 IRQ = 1, NREOQ

!              --- processing of output instructions necessary
!                  for collection

!              --- check time of output action
               DIF = TFINC - TIMCO
               IF ( IT.EQ.IT0 .AND. IC.EQ.1 ) THEN
                  CORQ%OQR(1) = REAL(IONOD(MXMIP+1,IRQ,0))
               END IF
               TNEXT = CORQ%OQR(1)
               IF ( ABS(DIF).LT.0.5*DT .AND. CORQ%OQR(2).LT.0. ) THEN
                  CORQ%OQR(1) = TIMCO
               ELSE IF ( CORQ%OQR(2).GT.0. .AND. TIMCO.GE.TNEXT ) THEN
                  CORQ%OQR(1) = TIMCO + CORQ%OQR(2)
               ELSE
                  GOTO 50
               END IF

               RTYPE  = CORQ%RQTYPE
               SNAMPF = CORQ%PSNAME

               CUOPS => FOPS
               DO
                 IF (CUOPS%PSNAME.EQ.SNAMPF) EXIT
                 IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) GOTO 50
                 CUOPS => CUOPS%NEXTOPS
               END DO
               PSTYPE = CUOPS%PSTYPE

               IF ( PSTYPE.EQ.'F' .OR. PSTYPE.EQ.'H' ) THEN
                  MXK = CUOPS%OPI(1)
                  MYK = CUOPS%OPI(2)
                  MIP = MXK * MYK
               ELSE IF ( PSTYPE.EQ.'C' .OR. PSTYPE.EQ.'P' .OR.
     &                   PSTYPE.EQ.'N' ) THEN
                  MXK = 0
                  MYK = 0
                  MIP = CUOPS%MIP
               END IF

!              --- rewrite table output by means of collection of
!                  output locations

               IF ( RTYPE(1:3).EQ.'TAB' ) THEN
                  CALL SWCOLTAB ( RTYPE, CORQ%OQI, CORQ%IVTYP, MIP, IRQ,
     &                            IONOD )
                  IF (STPNOW()) RETURN
               END IF

!              --- rewrite spectral output by means of collection of
!                  output locations

               IF ( RTYPE(1:2).EQ.'SP' ) THEN
                  CALL SWCOLSPC ( RTYPE, CORQ%OQI, MIP, IRQ, IONOD )
                  IF (STPNOW()) RETURN
               END IF

!              --- rewrite block output by means of collection of process
!                  output data

               IF ( RTYPE(1:3).EQ.'BLK' ) THEN
                  CALL SWCOLBLK ( RTYPE, CORQ%OQI, CORQ%IVTYP, CORQ%FAC,
     &                            SNAMPF, MXK, MYK, IRQ, IONOD )
                  IF (STPNOW()) RETURN
               END IF

  50           CONTINUE
               CORQ => CORQ%NEXTORQ

 100        CONTINUE
            IF ( NSTATC.EQ.1.AND.IT.LT.IT1 ) TIMCO = TIMCO + DT
         END DO

      END DO

!     --- close all files and delete process files

      DO IUNIT = HIOPEN+1, HIOPEN+NREOQ
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE(IUNIT)
      END DO
      DO IUNIT = HIOPEN+NREOQ+1, HIOPEN+NREOQ*(NPROC+1)
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE ( UNIT=IUNIT, STATUS='delete' )
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLTAB ( RTYPE, OQI, IVTYP, MIP, IRQ, IONOD )
!
!****************************************************************
!
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm1.inc'
      INCLUDE 'swcomm3.inc'
      INCLUDE 'swcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May 03: New subroutine
!
!  2. Purpose
!
!     Printing of table output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     IONOD       array indicating in which subdomain and at what time
!                 output points are located
!     IRQ         request number
!     IVTYP       type of variable output
!     MIP         total number of output points
!     OQI         array containing output request data
!     RTYPE       type of output request
!
      INTEGER   IONOD(MXMIP+1,MXOURQ,0:NPROC-1), IRQ, MIP, OQI(4),
     &          IVTYP(OQI(3))
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EMPTY :     logical whether a line is empty or not
!     EXIST :     logical whether a file exist or not
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     JVAR  :     loop counter
!     MSGSTR:     string to pass message to call MSGERR
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     OUTLIN:     output line
!
      INTEGER       IENT, IF, IL, ILPOS, IP, IPROC, IUNIT,
     &              IUT, IVTYPE, JVAR, NLINES, NREF, NVAR
      LOGICAL       EMPTY, EXIST, OPENED
      CHARACTER*80  MSGSTR
      CHARACTER*512 OUTLIN
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
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
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read line-by-line of each process file and write one
!     of them in appropriate manner to generic output file
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLTAB')

      NREF   = OQI(1)
      NVAR   = OQI(3)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of heading

            IF ( RTYPE.NE.'TABD' ) THEN
               IF ( RTYPE.EQ.'TABP' .OR. RTYPE.EQ.'TABI' ) THEN
                  NLINES = NLINES + 7
               ELSE IF ( RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABS' ) THEN
                  NLINES = NLINES + 5
                  IF ( RTYPE.EQ.'TABT' ) THEN
                     NLINES = NLINES + 1
                  ELSE
                     IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2
                     NLINES = NLINES + 2 + MIP
                  END IF
                  DO JVAR = 1, NVAR
                     IVTYPE = IVTYP(JVAR)
                     IF ( OVSVTY(IVTYPE).LE.2 ) THEN
                        NLINES = NLINES + 3
                     ELSE
                        NLINES = NLINES + 6
                    END IF
                  END DO
               END IF
            END IF

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

!     --- read line-by-line of each process file and write one
!         of them in appropriate manner to generic output file

      IF ( NREF.NE.PRINTF ) THEN
         IF ( RTYPE.EQ.'TABS' .AND. NSTATM.EQ.1 ) THEN
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               READ (IUNIT,'(A)') OUTLIN
               CALL TXPBLA(OUTLIN,IF,IL)
               IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO
         END IF
         DO IP = 1, MIP
            EMPTY = .TRUE.
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               READ (IUNIT,'(A)') OUTLIN
               CALL TXPBLA(OUTLIN,IF,IL)
               IF ( EMPTY .AND. IONOD(IP,IRQ,IPROC-1).EQ.IPROC ) THEN
                  WRITE (NREF, '(A)') OUTLIN(1:IL)
                  EMPTY = .FALSE.
               END IF
            END DO
            IF (EMPTY) WRITE (NREF, '(A)') OUTLIN(1:IL)
         END DO
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLSPC ( RTYPE, OQI, MIP, IRQ, IONOD )
!
!****************************************************************
!
      USE OUTP_DATA
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm1.inc'
      INCLUDE 'swcomm3.inc'
      INCLUDE 'swcomm4.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.30: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May 03: New subroutine
!
!  2. Purpose
!
!     Printing of spectral output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     IONOD       array indicating in which subdomain and at what time
!                 output points are located
!     IRQ         request number
!     MIP         total number of output points
!     OQI         array containing output request data
!     RTYPE       type of output request
!
      INTEGER   IONOD(MXMIP+1,MXOURQ,0:NPROC-1), IRQ, MIP, OQI(4)
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EMPTY :     logical whether a line is empty or not
!     EXIST :     logical whether a file exist or not
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IS    :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     MSGSTR:     string to pass message to call MSGERR
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     OPENED:     logical whether a file is open or not
!     OTYPE :     integer indicating dimension of spectrum
!     OUTLIN:     output line
!
      INTEGER       IENT, IF, IL, ILPOS, IP, IPROC, IS,
     &              IUNIT, IUT, NLINES, NREF, OTYPE
      LOGICAL       EMPTY, EXIST, OPENED
      CHARACTER*80  MSGSTR
      CHARACTER (LEN=LENSPO) OUTLIN
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
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
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read line-by-line of each process file and write one
!     of them in appropriate manner to generic output file
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLSPC')

      NREF   = OQI(1)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of first part of heading
           
		    NLINES = NLINES + 5
            IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

!           --- write coordinates of output points to generic
!               output file

            DO IP = 1, MIP
               EMPTY = .TRUE.
               DO IPROC = 1, NPROC
                  IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF ( EMPTY .AND. IONOD(IP,IRQ,IPROC-1).EQ.IPROC ) THEN
                     WRITE (NREF, '(A)') OUTLIN(1:IL)
                     EMPTY = .FALSE.
                  END IF
               END DO
               IF (EMPTY) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO

!           --- count lines of rest of heading and write heading
!               to generic output file

            NLINES = 2 + MSC
            IF (RTYPE(4:4).EQ.'C') THEN
               NLINES = NLINES + 7 + MDC
            ELSE
               NLINES = NLINES + 11
            END IF
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

      IF (RTYPE(4:4).EQ.'C') THEN
         IF (RTYPE.EQ.'SPEC') THEN
            OTYPE = -2
         ELSE
            OTYPE =  2
         END IF
      ELSE
         IF (RTYPE.EQ.'SPE1') THEN
            OTYPE = -1
         ELSE
            OTYPE =  1
         END IF
      END IF

!     --- read line-by-line of each process file and write one
!         of them in appropriate manner to generic output file

      IF ( NSTATM.EQ.1 ) THEN
         DO IPROC = 1, NPROC
            IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
            READ (IUNIT,'(A)') OUTLIN
            CALL TXPBLA(OUTLIN,IF,IL)
            IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
         END DO
      END IF

      DO IP = 1, MIP
         EMPTY = .TRUE.
         DO IPROC = 1, NPROC
            IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
            READ (IUNIT,'(A)') OUTLIN
            CALL TXPBLA(OUTLIN,IF,IL)
            IF ( EMPTY .AND. IONOD(IP,IRQ,IPROC-1).EQ.IPROC ) THEN
               WRITE (NREF, '(A)') OUTLIN(1:IL)
               IF ( OUTLIN(IF:IL).NE.'NODATA' ) THEN
                  IF ( ABS(OTYPE).EQ.1 ) THEN
                     DO IS = 1, MSC
                        READ (IUNIT,'(A)') OUTLIN
                        CALL TXPBLA(OUTLIN,IF,IL)
                        WRITE (NREF, '(A)') OUTLIN(1:IL)
                     END DO
                  ELSE
                     IF ( OUTLIN(IF:IL).NE.'ZERO' ) THEN
                        DO IS = 1, 1+MSC
                           READ (IUNIT,'(A)') OUTLIN
                           CALL TXPBLA(OUTLIN,IF,IL)
                           WRITE (NREF, '(A)') OUTLIN(1:IL)
                        END DO
                     END IF
                  END IF
               END IF
               EMPTY = .FALSE.
            END IF
         END DO
         IF (EMPTY) WRITE (NREF, '(A)') OUTLIN(1:IL)
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLBLK ( RTYPE, OQI, IVTYP, FAC  , PSNAME,
     &                      MXK  , MYK, IRQ  , IONOD )
!
!****************************************************************
!
      USE OUTP_DATA
      USE M_PARALL
!
      IMPLICIT NONE
!
      INCLUDE 'ocpcomm4.inc'
      INCLUDE 'swcomm1.inc'
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
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
!     40.31: Marcel Zijlema
!
!  1. Updates
!
!     40.31, Dec. 03: New subroutine
!
!  2. Purpose
!
!     Writing of block output by means of collecting
!     individual process output files
!
!  4. Argument variables
!
!     FAC         factors of multiplication of block output
!     IONOD       array indicating in which subdomain and at what time
!                 output points are located
!     IRQ         request number
!     IVTYP       type of variable output
!     MXK         number of points in x-direction of output frame
!     MYK         number of points in y-direction of output frame
!     OQI         array containing output request data
!     PSNAME      name of output locations
!     RTYPE       type of output request
!
      INTEGER   IONOD(MXMIP+1,MXOURQ,0:NPROC-1), MXK, MYK, IRQ, OQI(4),
     &          IVTYP(OQI(3))
      REAL      FAC(OQI(3))
      CHARACTER RTYPE*4, PSNAME*8
!
!  6. Local variables
!
!     CTIM  :     string representing date of computation
!     DFAC  :     multiplication factor of block output
!     EXIST :     logical whether a file exist or not
!     FMAX  :     auxiliary real
!     FTIP  :     auxiliary real
!     FTIP1 :     auxiliary real
!     FTIP2 :     auxiliary real
!     IDLA  :     lay-out indicator
!     IF    :     first non-character in string
!     IFAC  :     auxiliary integer
!     IENT  :     number of entries
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPD   :     switch for printing on paper or writing to file
!     IPROC :     loop counter
!     IREC  :     direct access file record counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     MATLAB:     indicates whether binary Matlab files are used
!     MSGSTR:     string to pass message to call MSGERR
!     NAMVAR:     name of MATLAB variable
!     NREF  :     unit reference number
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     OTMP  :     temporary output array
!     VOQ   :     collected output variables
!
      INTEGER      IDLA, IENT, IF, IFAC, IL, ILPOS, IP, IPD, IPROC,
     &             IUNIT, IUT, IXK, IYK, IVTYPE, JVAR, NREF, NVAR
      REAL         DFAC, FMAX, FTIP, FTIP1, FTIP2
      LOGICAL      EXIST, OPENED
      INTEGER, SAVE :: IREC=0
      LOGICAL, SAVE :: MATLAB=.FALSE.
	  CHARACTER*80 MSGSTR
      CHARACTER (LEN=20) :: CTIM
      CHARACTER (LEN=30) :: NAMVAR
      REAL, ALLOCATABLE :: VOQ(:,:), OTMP(:)
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     SBLKPT           Writes block output to an ASCII file
!     STRACE           Tracing routine for debugging
!     SWRMAT           Writes block output to a binary Matlab file
!     TABHED           Prints heading
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
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
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           open individual process output files
!
!     read data of each process file and write it in
!     appropriate manner to generic output file
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLBLK')

      NREF = OQI(1)
      NVAR = OQI(3)
      IDLA = OQI(4)

      IF ( RTYPE.EQ.'BLKP' ) THEN
        IPD = 1
        IF (NREF.EQ.PRINTF) CALL TABHED ('SWAN', PRINTF)
      ELSE IF ( RTYPE.EQ.'BLKD' ) THEN
        IPD = 2
      ELSE
        IPD = 3
      ENDIF

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

            MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.
     &               INDEX (FILENM, '.mat' ).NE.0
            IF (MATLAB) THEN
               CLOSE(NREF)
               OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',
     &              ACCESS='DIRECT', RECL=1)
               IREC = 1
            END IF

!           --- open individual process output files

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
            END DO

         END IF

      END IF

!     --- read data of each process file and write it in
!         appropriate manner to generic output file

      CTIM = CHTIME
      CALL TXPBLA(CTIM,IF,IL)
      CTIM(9:9)='_'
      
	  ALLOCATE(VOQ(MXMIP,2))
      ALLOCATE(OTMP(MXMIP))

      DO JVAR = 1, NVAR

         IVTYPE = IVTYP(JVAR)
         DFAC   = FAC(JVAR)

         VOQ = OVEXCV(IVTYPE)

         IF ( OVSVTY(IVTYPE).LT.3 ) THEN
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               DO IYK = 1, MYK
                  IP = (IYK-1)*MXK
                  READ (IUNIT, FLT_BLKP) (OTMP(IP+IXK), IXK=1,MXK)
               END DO
               DO IP = 1, MXMIP
                  IF (IONOD(IP,IRQ,IPROC-1).EQ.IPROC) VOQ(IP,1)=OTMP(IP)
               END DO
            END DO
         ELSE
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               DO IYK = 1, MYK
                  IP = (IYK-1)*MXK
                  READ (IUNIT, FLT_BLKP) (OTMP(IP+IXK), IXK=1,MXK)
               END DO
               DO IP = 1, MXMIP
                  IF (IONOD(IP,IRQ,IPROC-1).EQ.IPROC) VOQ(IP,1)=OTMP(IP)
               END DO
               DO IYK = 1, MYK
                  IP = (IYK-1)*MXK
                  READ (IUNIT, FLT_BLKP) (OTMP(IP+IXK), IXK=1,MXK)
               END DO
               DO IP = 1, MXMIP
                  IF (IONOD(IP,IRQ,IPROC-1).EQ.IPROC) VOQ(IP,2)=OTMP(IP)
               END DO
            END DO
         END IF

         IF ( IPD.EQ.1 ) THEN
            IF ( DFAC.LE.0. ) THEN
               IF ( OVHEXP(IVTYPE).LT.0.5E10 ) THEN
                  IFAC = INT (10.+LOG10(OVHEXP(IVTYPE))) - 13
               ELSE
                  IF ( OVSVTY(IVTYPE).EQ.1 ) THEN
                     FMAX = 1.E-8
                     DO IP = 1, MXK*MYK
                        FTIP = ABS(VOQ(IP,1))
                        FMAX = MAX (FMAX, FTIP)
                     END DO
                  ELSE IF ( OVSVTY(IVTYPE).EQ.2 ) THEN
                     FMAX = 1000.
                  ELSE IF ( OVSVTY(IVTYPE).EQ.3 ) THEN
                     FMAX = 1.E-8
                     DO IP = 1, MXK*MYK
                        FTIP1 = ABS(VOQ(IP,1))
                        FTIP2 = ABS(VOQ(IP,2))
                        FMAX  = MAX (FMAX, FTIP1, FTIP2)
                     END DO
                  END IF
                  IFAC = INT (10.+LOG10(FMAX)) - 13
               END IF
               DFAC = 10.**IFAC
            END IF
         ELSE
           IF ( DFAC.LE.0. ) DFAC = 1.
         END IF

         IF (OVSVTY(IVTYPE) .LT. 3) THEN
            IF (MATLAB) THEN
               IF (IL.EQ.1) THEN
                  NAMVAR = OVSNAM(IVTYPE)
               ELSE
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_'//CTIM
               END IF
               CALL SWRMAT( MYK, MXK, NAMVAR, VOQ(1,1), NREF,
     &                      IREC, IDLA, OVEXCV(IVTYPE) )
            ELSE
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE), VOQ(1,1) )
            END IF
         ELSE
            IF (MATLAB) THEN
               IF (IL.EQ.1) THEN
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_x'
               ELSE
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_x_'//CTIM
               END IF
               CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,1), NREF, IREC, IDLA, OVEXCV(IVTYPE) )
               IF (IL.EQ.1) THEN
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_y'
               ELSE
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_y_'//CTIM
               END IF
               CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,2), NREF, IREC, IDLA, OVEXCV(IVTYPE) )
            ELSE
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE)//'X-comp',
     &                      VOQ(1,1) )
               CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                      MXK, MYK, IDLA, OVLNAM(IVTYPE)//'Y-comp',
     &                      VOQ(1,2) )
            END IF
         END IF

      END DO

      IF (IPD.EQ.1) WRITE (PRINTF, 111)

      DEALLOCATE(VOQ,OTMP)

  111 FORMAT (///)

      RETURN
      END
