       PROGRAM MK_TIDE
       CHARACTER(LEN=80) :: TIDE_FILE
       INTEGER :: NumConstituent, NumEtaPoint
       INTEGER ::  I,J,K
! wilmington tide
       REAL,DIMENSION(:,:),ALLOCATABLE :: &
                  TidePeriod,TideAmp,TidePha
       REAL,DIMENSION(:),ALLOCATABLE :: &
                  PeriodLoc,AmpLoc,PhaLoc
       REAL,DIMENSION(:),ALLOCATABLE :: TideFac,Tideu0
       REAL :: TideStartDate
       REAL,DIMENSION(:),ALLOCATABLE :: ETA_bd
       INTEGER,DIMENSION(:),ALLOCATABLE :: I_eta_bd,J_eta_bd


       NumConstituent=1
       NumEtaPoint = 220
       ALLOCATE(TidePeriod(NumConstituent,NumEtaPoint), &
                TideAmp(NumConstituent,NumEtaPoint), &
                TidePha(NumConstituent,NumEtaPoint), &
                PeriodLoc(NumConstituent), &
                AmpLoc(NumConstituent), &
                PhaLoc(NumConstituent), &
                TideFac(NumConstituent), &
                Tideu0(NumConstituent), &
                ETA_bd(NumEtaPoint), &
                I_eta_bd(NumEtaPoint), &
                J_eta_bd(NumEtaPoint) &
                           )

       TIDE_FILE='tide_inlet.txt'
       TideStartDate=0.0
       TideFac(1)=1.0
       Tideu0=0.0

! M2(1) at Wilmington NC
       PeriodLoc(1)=12.00
       AmpLoc(1)=1.0
       PhaLoc(1)=0.0
       

       DO I=1,NumEtaPoint
         I_eta_bd(I)=I
         J_eta_bd(I)=1
         DO J=1,NumConstituent
           TidePeriod(J,I)=PeriodLoc(J)
           TideAmp(J,I)=AmpLoc(J)
           TidePha(J,I)=PhaLoc(J)
         ENDDO
       ENDDO

       OPEN(2,FILE=TRIM(TIDE_FILE))
         WRITE(2,*)'tidal boundary conditions'

         WRITE(2,'(f6.1)') TideStartDate
         WRITE(2,'(I3)') NumConstituent

         DO I = 1,NumConstituent
           WRITE(2,'(2F12.3)') TideFac(I),Tideu0(I)
         ENDDO
         WRITE(2,*) NumEtaPoint

         DO I=1,NumEtaPoint
          WRITE(2,'(2I5)')I_eta_bd(I),J_eta_bd(I)
          DO J=1,NumConstituent
            WRITE(2,'(3f12.3)')TidePeriod(J,I),TideAmp(J,I),TidePha(J,I)
          ENDDO
         ENDDO
       CLOSE(2)


       END
