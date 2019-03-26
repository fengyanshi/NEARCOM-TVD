!  master_time_ctr was modified from subroutine swread in swanpre1.f
! the purpose is to link swan and master program
        subroutine master_time_ctr(COMPUT)
!***********************************************************************

!     Modules

      USE OUTP_DATA                                                       40.13
      USE M_SNL4                                                          40.17
      USE M_GENARR                                                        40.31
      USE M_OBSTA                                                         40.31
      USE M_PARALL                                                        40.31
!   add pass module  - fyshi
      USE PASS
             
      INCLUDE 'timecomm.inc'                                              30.74
      INCLUDE 'ocpcomm1.inc'                                              30.74
!     ocpcomm2.inc is now accessible via USE OUTP_DATA                    40.13
      INCLUDE 'ocpcomm3.inc'                                              30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm1.inc'                                               30.74
      INCLUDE 'swcomm2.inc'                                               30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'                                               30.74

      LOGICAL, SAVE :: RUNMADE = .FALSE.                                  40.03
      CHARACTER PSNAME *8, PNAME *8, COMPUT *(*), PTYPE *1, DTTIWR *18    40.00
! define MasterDT
      REAL  :: MasterDT
!
!  calculate time variables: TINIC, TFINC, TIMCO, RDTIM, MTC, ITMOPT -fyshi
!                       and RUNMADE COMPUT -fyshi

        MasterDT= Master_dt*N_Interval_CallWave

        RUNMADE=.TRUE.
        COMPUT='COMP'
        TINIC=TIMCO
        TFINC=TINIC + MasterDT
        RDTIM=1./DT
        ITMOPT=1
        print*,'TINIC',TINIC,'TFINC',TFINC

        IF (NSTATM.GT.0) CHTIME = DTTIWR(ITMOPT, TIMCO)                   40.00
        NCOMPT = NCOMPT + 1
!  close error message for master program. - fyshi
!        IF (NCOMPT.GT.10) CALL MSGERR (2,
!     &                   'No more than 10 COMPUTE commands are allowed')
!        RCOMPT(NCOMPT,1) = REAL(NSTATC)
!        RCOMPT(NCOMPT,2) = REAL(MTC)
!        RCOMPT(NCOMPT,3) = TFINC
!        RCOMPT(NCOMPT,4) = TINIC
!        RCOMPT(NCOMPT,5) = DT
!     ** Next lines to process the begining and interval times to read
!       the boundary spectra from coarse grid if this is a nested grid.
        IF (NSTATM.EQ.1 .AND. NESRUN .EQ. 1) THEN                              30.00
          CALL DTSTTI(ITMOPT, BEGBOC, TIMARR)
          BEGBOU = DTTIME(TIMARR)
          TIMERB = BEGBOU
        ENDIF
!       set ITERMX equal to MXITST in case of stationary computations
!       and to MXITNS otherwise
        IF (NSTATC.EQ.0) THEN
          ITERMX = MXITST                                                 40.03
        ELSE
          ITERMX = MXITNS                                                 40.03
        ENDIF
        RETURN                                                            30.00
!      ENDIF
      end


