
      SUBROUTINE PASS_OUT_SWAN (COMPDA,KGRPNT,XCGRID,YCGRID,
     &           MXK,MYK,VOQR,VOQ)
 
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
      integer IX,IY, INDX,KGRPNT(MXC,MYC),MXK,MYK
      integer VOQR(*)
      real VOQ(MXK*MYK,*)
! local variables
      real x(mxc,myc),y(mxc,myc),ang(mxc,myc)

       PRINT*,'PASS_OUT_SWAN ...'

         if(MXC.ne.MXK.or.MYC.ne.MYK)then
         print*,'MXK=',MXK, '  MYK=',MYK
         print*,'MXC=',MXC, '  MYC',MYC
         print*, 'reset MXK or MYK in INPUT file'
         stop
         endif

c --- peak direction
         if(VOQR(14).gt.0)then
         print*,'passing peak direction'
         DO IY=1,MYK
         DO IX=1,MXK
          INDX=(IY-1)*MXK+IX
          pass_theta(ix,iy)=VOQ(INDX,VOQR(14))
          if(pass_theta(ix,iy).le.-900.)pass_theta(ix,iy)=0.
         ENDDO
         ENDDO
         DO IX=1,MXK
          pass_theta(ix,myk)=pass_theta(ix,myk-1)
         ENDDO
         endif

c --- peak period
         if(VOQR(12).gt.0)then
         print*,'passing peak period'
         DO IX=1,MXK
         DO IY=1,MYK
          INDX=(IY-1)*MXK+IX
          pass_PPER(ix,iy)=VOQ(INDX,VOQR(12))
         ENDDO
         ENDDO
         DO IX=1,MXK
          pass_pper(ix,myk)=pass_pper(ix,myk-1)
         ENDDO
         endif

c --- ubott and Hs
         DO IX=1,MXC
         DO IY=1,MYC
          INDX=KGRPNT(IX,IY)
          IF(INDX.GT.1)THEN
          pass_ubott(ix,iy)=COMPDA(INDX,JUBOT)
          pass_height(ix,iy)=COMPDA(INDX,JHS)
          ENDIF
         ENDDO
         ENDDO

c --- Pass_C

         DO IX=1,MXK
         DO IY=1,MYK
          pass_C(ix,iy)=sqrt(9.81*Depth_Wave(ix,iy))
         ENDDO
         ENDDO

c     ccalculate wave mass flux
       do IY=1,MYC
        do IX=1,MXC
        if(Pass_C(ix,iy).gt.0.1)then
       Pass_MassFlux(ix,iy)=9.8*Pass_Height(ix,iy)
     *        *Pass_Height(ix,iy)/Pass_C(ix,iy)*1./8.
        else
        Pass_MassFlux(ix,iy)=0.
        endif
       Pass_MassFluxU(ix,iy)=Pass_MassFlux(ix,iy)
     &     *cos(Pass_Theta(ix,iy)*3.1415926/180.)
       Pass_MassFluxV(ix,iy)=Pass_MassFlux(ix,iy)
     &     *sin(Pass_Theta(ix,iy)*3.1415926/180.)
        enddo
        enddo

                   
      RETURN
      END

