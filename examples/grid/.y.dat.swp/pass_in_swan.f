
      SUBROUTINE PASS_IN_SWAN (COMPDA,KGRPNT,XCGRID,YCGRID)
 
!
      USE OUTP_DATA                                                       40.31
      USE M_PARALL  
      USE PASS                                                            40.31
!      

      implicit none
      INCLUDE 'timecomm.inc'                                              30.74
      INCLUDE 'ocpcomm4.inc'                                              30.74
      INCLUDE 'swcomm1.inc'                                               30.74
      INCLUDE 'swcomm2.inc'                                               30.74
      INCLUDE 'swcomm3.inc'                                               30.74
      INCLUDE 'swcomm4.inc'

      REAL COMPDA(MCGRD,MCMVAR),XCGRID(MXC,MYC),YCGRID(MXC,MYC)
      integer IX,IY, INDX,KGRPNT(MXC,MYC)
! local variables
      real x(mxc,myc),y(mxc,myc),ang(mxc,myc)

       PRINT*,'PASS_IN_SWAN ...'

        DO IY=1,MYC
        DO IX=1,MXC
         X(IX,IY)=XCGRID(IX,IY)
         Y(IX,IY)=YCGRID(IX,IY)
        ENDDO
        ENDDO

          do IY=1,MYC
          do IX=1,MXC-1
           ang(IX,IY)=atan2(y(IX+1,IY)-y(IX,IY),x(IX+1,IY)-x(IX,IY))
          enddo
          enddo
          do IY=1,MYC
           ang(MXC,IY)=ang(MXC-1,IY)
          enddo

!
! Move values at 'present' time level 2 to 'old' time level 1.
! MCGRD = MXC*MYC+1-#masked cells.
! MXC = # cells x-dir in this tile including halox.
! MYC = # cells y-dir in this tile including haloy.
! COMPDA has only active points + 1.
!
        DO INDX = 2, MCGRD
          COMPDA(INDX,JWLV1)=COMPDA(INDX,JWLV2)
          COMPDA(INDX,JVX1)=COMPDA(INDX,JVX2)
          COMPDA(INDX,JVY1)=COMPDA(INDX,JVY2)
        ENDDO
!
!        DO INDX = 2, MCGRD
!          COMPDA(INDX,JWLV2)=0.0
!          COMPDA(INDX,JVX2)=0.5
!          COMPDA(INDX,JVY2)=0.5
!        ENDDO

         DO IX=1,MXC
         DO IY=1,MYC
          INDX=KGRPNT(IX,IY)
          IF(INDX.GT.1)THEN
         COMPDA(INDX,JWLV2)=Intp_eta_Wave(IX,IY)
          COMPDA(INDX,JVX2)=Intp_U_Wave(IX,IY)*COS(ANG(IX,IY))
     &                     +Intp_V_Wave(IX,IY)*SIN(ANG(IX,IY))
          COMPDA(INDX,JVY2)=Intp_V_Wave(IX,IY)*COS(ANG(IX,IY))
     &                     -Intp_U_Wave(IX,IY)*SIN(ANG(IX,IY))
!          COMPDA(INDX,JVX2)=0.5
!          COMPDA(INDX,JVY2)=0.5
          ENDIF
         ENDDO
         ENDDO
         INDX=KGRPNT(10,10)

!  convert wave force to real x and y directions
         DO IX=1,MXC
         DO IY=1,MYC
          Pass_Wave_Fx(IX,IY)=0.
          Pass_Wave_Fy(IX,IY)=0.
         ENDDO
         ENDDO



 
                   
      RETURN
      END

