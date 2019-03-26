C******************************************************************************
C PADCIRC RELEASE VERSION 43.03 05/20/2003                                    *
C  last changes in this file VERSION 43.02                                    *
C                                                                             *
C  mod history                                                                *
C  v43.03     - 05/20/03 - rl - from 43.02 - parallel wind stuff (m.brown)    *
C                                          output buffer flush (m.cobb)       *
C                                          3D fixes (k.dresback)              *
C                                          drop MNPROC in fort.15 (t.campbell)*
C                                          various bug fixes in RBCs          *
c                                          ZSURFBUOY/BCPG calc                *
C  v43.02     - 02/06/03 - rl - from 43.01 - code clean up & documentation    *
C  v43.01     - 01/31/03 - jf - from 43.00 - reconcile with v42.07 2D code    *
C                                            changed var. names: NNEIGH->     *
C                                            NNeigh, NEIGH->NeiTab,           *
C                                            NEIGHELE->NeiTabEle,             *
C  v43.00a    - sum  /02 - tc - from 36.01 (3D) & 41.12? (2D), create F90/    *
C                                               parallel unified 2D/3D source *
C  v41.11 - 09/14/01 - RL - eliminated MNWLAT, MNWLON                         *
C  v41.06 - 04/02/01 - RL - eliminated MNWP                                   *
C                                                                             *
C******************************************************************************
C
      MODULE SIZES
      IMPLICIT NONE

C    for curvcirc
      integer, parameter :: nx_max =310
      integer, parameter :: ny_max =450

C
C...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS       
C...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

      INTEGER, PARAMETER :: SZ = 4
      INTEGER, PARAMETER :: NBYTE=4

C...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
C...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
C   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
C   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.
 
      INTEGER, PARAMETER ::  NPREC=7
C
      INTEGER ::  MNPROC,MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,
     &  MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF
C
C     Dimension of vertical FE mesh (To interface 2D & 3D)       - should be moved to global_3dv? (rl 5/2003)
      INTEGER :: MNODES

      LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV
C
C For Definition of Working Directory
C
      INTEGER,SAVE :: MYPROC

      INTEGER,SAVE :: LNAME = 1
      CHARACTER*1,SAVE :: DIRNAME = '.'
      
C---------------------end of data declarations--------------------------------C


      CONTAINS


      SUBROUTINE MAKE_DIRNAME()
      MYPROC=0
      RETURN
      END SUBROUTINE

      END MODULE SIZES
