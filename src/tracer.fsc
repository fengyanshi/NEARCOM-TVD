#if defined (TRACKING)
MODULE TRACER
     USE GLOBAL, ONLY : TIME, Jbeg,Jend,Ibeg,Iend,Mloc,Nloc,&
                 Mloc1,Nloc1,MASK,X,Y,U,V,SP,DT,Ubott,Vbott,&
                 Usurf,Vsurf

# if defined (PARALLEL)
     USE MPI
     USE GLOBAL, ONLY : myid,ier,npx,npy,NumberProcessor,px,py,MPI_SP,&
                        n_east,n_west,n_suth,n_nrth
# endif
     IMPLICIT NONE


     INTEGER :: I,J,II,JJ
     REAL(SP),SAVE :: TRACK_TIME_COUNT = 0.0_SP
     REAL(SP),SAVE,DIMENSION(:),ALLOCATABLE :: X_TRACK,Y_TRACK,Sc,S1,S2,S3,&
                   U_TRACER,V_TRACER,X_TRACK_PRE,Y_TRACK_PRE, &
                   Sc_pre,S1_pre,S2_pre,S3_pre
                       
     INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: nx1,ny1,nx2,ny2,nx3,ny3,&
        nx1_pre,ny1_pre,nx2_pre,ny2_pre,nx3_pre,ny3_pre
     INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: LAYER_012
     LOGICAL,SAVE :: FirstCallTrack = .TRUE.
     REAL(SP),SAVE :: TrackerInterval
     REAL(SP),SAVE,DIMENSION(:),ALLOCATABLE :: TRACKER_START
     LOGICAL,SAVE,DIMENSION(:),ALLOCATABLE :: IN_CELL,STOP_SEARCH
     INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: FOUND_IN_DOMAIN
# if defined (PARALLEL)
     INTEGER::myint
# endif
     

     INTEGER,SAVE :: NumTracker=1


     CONTAINS

SUBROUTINE MAIN_TRACK

     IF(.NOT.ALLOCATED(Sc)) THEN
       ALLOCATE(  &
                Sc(NumTracker), &
                S1(NumTracker), &
                S2(NumTracker), &
                S3(NumTracker), &
                nx1(NumTracker), &
                nx2(NumTracker), &
                nx3(NumTracker), &
                ny1(NumTracker), &
                ny2(NumTracker), &
                ny3(NumTracker), &
                Sc_pre(NumTracker), &
                S1_pre(NumTracker), &
                S2_pre(NumTracker), &
                S3_pre(NumTracker), &
                nx1_pre(NumTracker), &
                nx2_pre(NumTracker), &
                nx3_pre(NumTracker), &
                ny1_pre(NumTracker), &
                ny2_pre(NumTracker), &
                ny3_pre(NumTracker), &
                IN_CELL(NumTracker), &
                FOUND_IN_DOMAIN(NumTracker), &
                STOP_SEARCH(NumTracker), &
                X_TRACK_PRE(NumTracker), &
                Y_TRACK_PRE(NumTracker), &
                U_TRACER(NumTracker), &
                V_TRACER(NumTracker) &
               )
        IN_CELL = .FALSE.
        FOUND_IN_DOMAIN = 0
        STOP_SEARCH = .FALSE.
     ENDIF

!   print*,myid,x_track(1),y_track(1),x_track(2),y_track(2)

     IF(FirstCallTrack)THEN
       FirstCallTrack = .FALSE.
       CALL GET_XY_POSITION
     ENDIF

     CALL TRACK_XY
     
END SUBROUTINE MAIN_TRACK

SUBROUTINE GET_XY_POSITION
      REAL(SP) :: x1,x2,x3,y1,y2,y3,&
                  area1,area2,area3
      INTEGER :: I_start,I_end,J_start,J_end,search_count

! ---  areas of the four triangles, area will be negative
!      if an order is clockwise
!       Sc -- triangle 1,2,3
!       S1 -- triangle 2,3,c
!       S2 -- triangle 3,1,c
!       S3 -- triangle 1,2,c

! --- find the triangle includes the points in grid2

  
        DO I=1,NumTracker
          x1=X_TRACK(I)
          y1=Y_TRACK(I)
          
          search_count=0
300       CONTINUE  ! re-search
          search_count=search_count+1

         IF(.NOT.STOP_SEARCH(I)) THEN

          IF(IN_CELL(I))THEN
           I_start=MAX(nx1(I)-2,Ibeg)
           I_end=MIN(nx1(I)+2,Iend)
           J_start=MAX(ny2(I)-2,Jbeg)
           J_end=MIN(ny2(I)+2,Jend)
          ELSEIF(FOUND_IN_DOMAIN(I)<1)THEN
# if defined (PARALLEL)
! use Iend-1 is a wrong idea for parallelization,
! since there's no x or y connected between processor interfaces
           I_start=Ibeg
           if ( n_east .eq. MPI_PROC_NULL ) then
            I_end=Iend-1
           else
            I_end=Iend
           endif
           J_start=Jbeg
           if ( n_nrth .eq. MPI_PROC_NULL ) then
             J_end=Jend-1
           else
             J_end=Jend
           endif
# else
           I_start=Ibeg
           I_end=Iend-1
           J_start=Jbeg
           J_end=Jend-1
# endif

          ELSE  ! for processors don't contain tracers
           I_start=Ibeg
           I_end=Ibeg
           J_start=Jbeg
           J_end=Jbeg 
          ENDIF

          do JJ=J_start,J_end
          do II=I_start,I_end

            x2=X(ii+1,jj)
            y2=Y(ii+1,jj)
            x3=X(ii,jj+1)
            y3=Y(ii,jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=X(ii,jj+1)
            y2=Y(ii,jj+1)
            x3=X(ii,jj)
            y3=Y(ii,jj)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=X(ii,jj)
            y2=Y(ii,jj)
            x3=X(ii+1,jj)
            y3=Y(ii+1,jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              IN_CELL(I) = .TRUE.
              FOUND_IN_DOMAIN(I) = 1
              nx1(I)=ii
              ny1(I)=jj
              nx2(I)=ii+1
              ny2(I)=jj
              nx3(I)=ii
              ny3(I)=jj+1
              S1(I)=area1
              S2(I)=area2
              S3(I)=area3
               
              x1=X(nx1(I),ny1(I))
              y1=Y(nx1(I),ny1(I))
              x2=X(nx2(I),ny2(I))
              y2=Y(nx2(I),ny2(I))
              x3=X(nx3(I),ny3(I))
              y3=Y(nx3(I),ny3(I))
              Sc(I)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 120
            endif

            x2=X(ii+1,jj)
            y2=Y(ii+1,jj)
            x3=X(ii+1,jj+1)
            y3=Y(ii+1,jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=X(ii+1,jj+1)
            y2=Y(ii+1,jj+1)
            x3=X(ii,jj+1)
            y3=Y(ii,jj+1)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=X(ii,jj+1)
            y2=Y(ii,jj+1)
            x3=X(ii+1,jj)
            y3=Y(ii+1,jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              IN_CELL(I) = .TRUE.
              FOUND_IN_DOMAIN(I) = 1
              nx1(I)=ii
              ny1(I)=jj+1
              nx2(I)=ii+1
              ny2(I)=jj
              nx3(I)=ii+1
              ny3(I)=jj+1
              S1(I)=area1
              S2(I)=area2
              S3(I)=area3
               
              x1=X(nx1(I),ny1(I))
              y1=Y(nx1(I),ny1(I))
              x2=X(nx2(I),ny2(I))
              y2=Y(nx2(I),ny2(I))
              x3=X(nx3(I),ny3(I))
              y3=Y(nx3(I),ny3(I))
              Sc(I)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 120
            endif

          enddo
          enddo

            IN_CELL(I) = .FALSE.
            FOUND_IN_DOMAIN(I) = 0

120       continue


# if defined (PARALLEL)
!  gather found_in_domain
     call MPI_ALLREDUCE(FOUND_IN_DOMAIN(I),myint,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier) 
          FOUND_IN_DOMAIN(I) = myint  

          IF(FOUND_IN_DOMAIN(I)<1)THEN
!            WRITE(*,*)'myid',myid,'TRACER',I,'LOST, RE-SEARCH'
            IF(search_count>1) THEN
!              WRITE(*,*)'myid',myid,'TRACER',I,'LOST permanently'
              STOP_SEARCH(I)=.TRUE.
            ELSE
              STOP_SEARCH(I)=.FALSE.
              GOTO 300
            ENDIF
          ENDIF ! couldn't find
   
# else
          IF(FOUND_IN_DOMAIN(I)<1)THEN
!            WRITE(*,*)'TRACER',I,'LOST, RE-SEARCH'
            IF(search_count>1) THEN
!              WRITE(*,*)'TRACER',I,'LOST permanently'
              STOP_SEARCH(I)=.TRUE.
            ELSE
              STOP_SEARCH(I)=.FALSE.
              GOTO 300
            ENDIF
          ENDIF ! couldn't find
# endif

        ENDIF ! end stop search

        ENDDO ! end numtracker

! print*,myid,in_cell(2),Y_TRACK(2),STOP_SEARCH(2),FOUND_IN_DOMAIN(2)
!      CALL SWEXITMPI  
        
END SUBROUTINE GET_XY_POSITION

SUBROUTINE TRACK_XY

       X_TRACK_PRE = X_TRACK
       Y_TRACK_PRE = Y_TRACK
       nx1_pre=nx1
       ny1_pre=ny1
       nx2_pre=nx2
       ny2_pre=ny2
       nx3_pre=nx3
       ny3_pre=ny3
       S1_pre=S1
       S2_pre=S2
       S3_pre=S3
       Sc_pre=Sc

       CALL MOVE_TRACER

# if defined (PARALLEL)
       CALL BROADCAST_XY
# endif

       CALL GET_XY_POSITION

       DO I=1,NumTracker
        IF(IN_CELL(I))THEN
         IF(MASK(nx1(I),ny1(I))+MASK(nx2(I),ny2(I))+MASK(nx2(I),ny2(I))<1)THEN
          nx1(I)=nx1_pre(I)
          ny1(I)=ny1_pre(I)
          nx2(I)=nx2_pre(I)
          ny2(I)=ny2_pre(I)
          nx3(I)=nx3_pre(I)
          ny3(I)=ny3_pre(I)
          S1(I)=S1_pre(I)
          S2(I)=S2_pre(I)
          S3(I)=S3_pre(I)
          Sc(I)=Sc_pre(I)
          X_TRACK(I)=X_TRACK_PRE(I)
          Y_TRACK(I)=Y_TRACK_PRE(I)
         ENDIF
        ENDIF ! end in cell
       ENDDO

END SUBROUTINE TRACK_XY


# if defined (PARALLEL)
SUBROUTINE BROADCAST_XY
     IMPLICIT NONE
     INTEGER :: l
     INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
     REAL(SP),DIMENSION(NumberProcessor) :: xx
     LOGICAL,DIMENSION(NumberProcessor) :: yy
       INTEGER :: IDsend
       REAL(SP),DIMENSION(NumTracker,NumberProcessor) :: X_TRACK_ALL,Y_TRACK_ALL
       LOGICAL,DIMENSION(NumTracker,NumberProcessor) :: IN_CELL_ALL

     call MPI_Gather(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)
     call MPI_Gather(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,ier)

! gather to myid=0
       DO I=1,NumTracker

! x
        call MPI_Gather(X_TRACK(I),1,MPI_SP,&
             xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)

!        if (I.eq.1) call MPI_Barrier(MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
               X_TRACK_ALL(I,l)= xx(l)
           enddo
        endif 

! y       
        call MPI_Gather(Y_TRACK(I),1,MPI_SP,&
             xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
               Y_TRACK_ALL(I,l)= xx(l)
           enddo
        endif 

! in_cell
        call MPI_Gather(IN_CELL(I),1,MPI_LOGICAL,&
             yy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

        if (myid.eq.0) then
           do l=1,px*py
               IN_CELL_ALL(I,l)= yy(l)
           enddo
        endif 

        ENDDO

       DO I=1,NumTracker
         DO l=1,px*py
          IF(IN_CELL_ALL(I,l))THEN
             X_TRACK(I)=X_TRACK_ALL(I,l)
             Y_TRACK(I)=Y_TRACK_ALL(I,l)
          ENDIF
         ENDDO
       ENDDO        

! send x,y to all processors


       DO I=1,NumTracker
! send x
        if (myid.eq.0) then
           do l=1,px*py
              xx(l) =  X_TRACK(I)
           enddo
        endif
        call MPI_Scatter(xx,1,MPI_SP,&
             X_TRACK(I),1,MPI_SP,0,MPI_COMM_WORLD,ier)
! send y
        if (myid.eq.0) then
           do l=1,px*py
              xx(l) =  Y_TRACK(I)
           enddo
        endif
        call MPI_Scatter(xx,1,MPI_SP,&
             Y_TRACK(I),1,MPI_SP,0,MPI_COMM_WORLD,ier)
       ENDDO


END SUBROUTINE BROADCAST_XY
# endif

SUBROUTINE MOVE_TRACER

        DO I=1,NumTracker
           IF(IN_CELL(I).AND.TIME>TRACKER_START(I))THEN            
              IF(Sc(I)==0.0)THEN
                write(*,*) 'Sc is ZERO for TRACER', I
              ELSE
                IF(LAYER_012(I)==1)THEN  ! surface
                  U_TRACER(I)=(S1(I)*Usurf(nx1(I),ny1(I)) &
                        +S2(I)*Usurf(nx2(I),ny2(I)) &
                        +S3(I)*Usurf(nx3(I),ny3(I))) &
                        /Sc(I) 
                  V_TRACER(I)=(S1(I)*Vsurf(nx1(I),ny1(I)) &
                        +S2(I)*Vsurf(nx2(I),ny2(I)) &
                        +S3(I)*Vsurf(nx3(I),ny3(I))) &
                        /Sc(I) 
                ELSEIF(LAYER_012(I)==2)THEN ! bottom
                  U_TRACER(I)=(S1(I)*Ubott(nx1(I),ny1(I)) &
                        +S2(I)*Ubott(nx2(I),ny2(I)) &
                        +S3(I)*Ubott(nx3(I),ny3(I))) &
                        /Sc(I) 
                  V_TRACER(I)=(S1(I)*Vbott(nx1(I),ny1(I)) &
                        +S2(I)*Vbott(nx2(I),ny2(I)) &
                        +S3(I)*Vbott(nx3(I),ny3(I))) &
                        /Sc(I)
                ELSE        ! depth averaged
                  U_TRACER(I)=(S1(I)*U(nx1(I),ny1(I)) &
                        +S2(I)*U(nx2(I),ny2(I)) &
                        +S3(I)*U(nx3(I),ny3(I))) &
                        /Sc(I) 
                  V_TRACER(I)=(S1(I)*V(nx1(I),ny1(I)) &
                        +S2(I)*V(nx2(I),ny2(I)) &
                        +S3(I)*V(nx3(I),ny3(I))) &
                        /Sc(I)
                ENDIF
                X_TRACK(I)=X_TRACK(I)+U_TRACER(I)*DT
                Y_TRACK(I)=Y_TRACK(I)+V_TRACER(I)*DT
              ENDIF 
            ENDIF ! in_cell    
        ENDDO

END SUBROUTINE MOVE_TRACER 
 

END MODULE TRACER
#endif
! end tracking


#if defined (DYE)
SUBROUTINE DYE_MODULE
       USE SEDI
       USE PASS, ONLY : DT_SC
       USE GLOBAL, ONLY : TIME,Jbeg,Jend,Ibeg,Iend,MinDepthFrc,H
       IMPLICIT NONE

      Vol_Dye(75+3:76+3,45+3) = 1.0_SP*H(75+3:76+3,45+3)      

      CALL ADVECTION_DIFFUSION(Vol_Dye,DT_SC)

      DO J=Jbeg,Jend
      DO I=Ibeg,Iend
        Con_Dye(I,J)=Vol_Dye(I,J)/Max(MinDepthFrc,H(I,J))
      ENDDO
      ENDDO

# if defined (PARALLEL)
     CALL PHI_EXCH (Vol_Dye)
# endif

END SUBROUTINE DYE_MODULE
# endif
SUBROUTINE ADVECTION_DIFFUSION(Con,DT)
      USE PARAM
      USE GLOBAL,ONLY : U,V,Mloc,Nloc,H,Mloc1,Nloc1,MASK,Ibeg,Iend,Jbeg,Jend, &
                         Jaco,Ubott,Vbott,tmp4preview,MinDepthFrc, &
                         L11,L12,L22,nu_total 
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mloc,Nloc),INTENT(INOUT) :: Con
      REAL(SP), INTENT(IN) :: DT
      REAL(SP),DIMENSION(Mloc,Nloc) :: Q1,Q2
      REAL(SP),DIMENSION(Mloc1,Nloc) :: Q1L,Q1R,Q1c
      REAL(SP),DIMENSION(Mloc,Nloc1) :: Q2L,Q2R,Q2c
      REAL(SP),DIMENSION(Mloc,Nloc) :: DelxQ1
      REAL(SP),DIMENSION(Mloc,Nloc) :: DelyQ2
      REAL(SP) :: HL11_1,HL11_2,HL12_1,HL12_2,HL22_1,HL22_2, &
                  C_11,C_12,C_22,C_1,C_2,HL11,HL12,HL22

      DO J=1,Nloc
      DO I=1,Mloc
        Q1(I,J)=(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))*Jaco(I,J)*Con(I,J)
        Q2(I,J)=(U(I,J)*L12(I,J)+V(I,J)*L22(I,J))*Jaco(I,J)*Con(I,J)
      ENDDO
      ENDDO

      CALL DelxFun_minmod(1.0_SP,Mloc,Nloc,Q1,DelxQ1)
      CALL DelyFun_minmod(1.0_SP,Mloc,Nloc,Q2,DelyQ2)
     
      CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,1.0_SP,Q1,DelxQ1,Q1L,Q1R)
      CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,1.0_SP,Q2,DelyQ2,Q2L,Q2R)

! 1st order TVD
       DO J=1,Nloc
       DO I=1,Mloc1
        IF((Q1L(I,J)+Q1R(I,J)).GE.ZERO) THEN
          Q1c(I,J)=Q1L(I,J)
        ELSE
          Q1c(I,J)=Q1R(I,J)
        ENDIF
       ENDDO
       ENDDO

       DO J=1,Nloc1
       DO I=1,Mloc
        IF((Q2R(I,J)+Q2L(I,J)).LE.ZERO) THEN
          Q2c(I,J)=Q2R(I,J)
        ELSE
          Q2c(I,J)=Q2L(I,J)
        ENDIF        
       ENDDO
       ENDDO

! deal with masks

      DO J=1,Nloc
      DO I=1,Mloc
         IF(MASK(I,J)<1)THEN
           Q1c(I,J)=ZERO
           Q1c(I+1,J)=ZERO
           Q2c(I,J+1)=ZERO
           Q2c(I,J)=ZERO
         ENDIF
      ENDDO
      ENDDO

      DO J=Jbeg,Jend
      DO I=Ibeg,Iend

       IF(MASK(I,J)>0)THEN
        Con(I,J)=Con(I,J)-1.0_SP/MAX(SMALL,Jaco(I,J))*(Q1c(I+1,J)-Q1c(I,J)  &
                               +Q2c(I,J+1)-Q2c(I,J))*DT                               
       ELSE
        Con(I,J)=ZERO
       ENDIF

      ENDDO
      ENDDO

! diffusion

      DO J=Jbeg,Jend
      DO I=Ibeg,Iend

       IF(MASK(I,J)>0)THEN
         HL11_1=0.5_SP*(H(I+1,J)*nu_total(I+1,J)*L11(I+1,J)-H(I-1,J)*nu_total(I-1,J)*L11(I-1,J))
         HL12_1=0.5_SP*(H(I+1,J)*nu_total(I+1,J)*L12(I+1,J)-H(I-1,J)*nu_total(I-1,J)*L12(I-1,J))
         HL22_1=0.5_SP*(H(I+1,J)*nu_total(I+1,J)*L22(I+1,J)-H(I-1,J)*nu_total(I-1,J)*L22(I-1,J))
         HL11_2=0.5_SP*(H(I,J+1)*nu_total(I,J+1)*L11(I,J+1)-H(I,J-1)*nu_total(I,J-1)*L11(I,J-1))
         HL12_2=0.5_SP*(H(I,J+1)*nu_total(I,J+1)*L12(I,J+1)-H(I,J-1)*nu_total(I,J-1)*L12(I,J-1))
         HL22_2=0.5_SP*(H(I,J+1)*nu_total(I,J+1)*L22(I,J+1)-H(I,J-1)*nu_total(I,J-1)*L22(I,J-1))
         HL11=H(I,J)*L11(I,J)*nu_total(I,J)
         HL12=H(I,J)*L12(I,J)*nu_total(I,J)
         HL22=H(I,J)*L22(I,J)*nu_total(I,J)

         C_1=(Con(I+1,J)-Con(I-1,J))*0.5_SP/MAX(SMALL,H(I,J))
         C_2=(Con(I,J+1)-Con(I,J-1))*0.5_SP/MAX(SMALL,H(I,J))
         C_11=Con(I+1,J)-2.0_SP*Con(I,J)+Con(I-1,J)/MAX(SMALL,H(I,J))
         C_22=Con(I,J+1)-2.0_SP*Con(I,J)+Con(I,J-1)/MAX(SMALL,H(I,J))
         C_12=0.25_SP*(Con(I+1,J+1)-Con(I-1,J+1)  &
            -Con(I+1,J-1)+Con(I-1,J-1))/MAX(SMALL,H(I,J))
         Con(I,J)=Con(I,J) + ( &
           HL11*(HL11*C_11+C_1*HL11_1+HL12*C_12+C_2*HL12_1) &
          +HL12*(HL11*C_12+C_1*HL11_2+HL12*C_22+C_2*HL12_2) &
          +HL12*(HL12*C_11+C_1*HL12_1+HL22*C_12+C_2*HL22_1) &
          +HL22*(HL12*C_12+C_1*HL12_2+HL22*C_22+C_2*HL22_2) )
       ELSE
        Con(I,J)=ZERO
       ENDIF

      ENDDO
      ENDDO


END SUBROUTINE ADVECTION_DIFFUSION
