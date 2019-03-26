c*************************************
        subroutine SediModule()
         
         use PASS

         real por, qx(Nx_Max,Ny_Max), qy(Nx_Max,Ny_Max),
     .        dhdt(Nx_max,Ny_max),ht(Nx_max,Ny_max),dx,dy

	  print*, 'sediment module .......'

         por = 0.35

         dx=X_Sedi(2,1)-X_Sedi(1,1)
         dy=Y_Sedi(1,2)-Y_Sedi(1,1)

          print*,Intp_sedfluxcum_x(10,10),Pass_sedfluxcum_x(10,10)
     &        ,Pass_sedflux_x(10,10)

         do i = 1, Nx_Sedi
          do j = 1, Ny_Sedi
             ht(i,j)=Depth_Sedi(i,j)
             qx(i,j) = Intp_sedfluxcum_x(i,j)/(1.-por)
             qy(i,j) = Intp_sedfluxcum_y(i,j)/(1.-por)      
          end do
         end do

c -- ez scheme for weak interaction (fyshi)
c --- h and dhdt defined at node

        do j=2,Ny_Sedi-1
        do i=2,Nx_Sedi-1
          dhdt(i,j)=(qx(i+1,j)-qx(i-1,j))/dx/2.
     &                 +(qy(i,j+1)-qy(i,j-1))/dy/2.
        enddo
        enddo

        j=1
        do i=2,Nx_Sedi-1
          dhdt(i,j)=(qx(i+1,j)-qx(i-1,j))/dx/2.
     &                 +(qy(i,j+1)-qy(i,j))/dy/1.
        enddo

        j=Ny_Sedi
        do i=2,Nx_Sedi-1
          dhdt(i,j)=(qx(i+1,j)-qx(i-1,j))/dx/2.
     &                 +(qy(i,j)-qy(i,j-1))/dy/1.
        enddo

        i=1
        do j=2,Ny_Sedi-1
          dhdt(i,j)=(qx(i+1,j)-qx(i,j))/dx/1.
     &                 +(qy(i,j+1)-qy(i,j-1))/dy/2.
        enddo

        i=Nx_Sedi
        do j=2,Ny_Sedi-1
          dhdt(i,j)=(qx(i,j)-qx(i-1,j))/dx/1.
     &                 +(qy(i,j+1)-qy(i,j-1))/dy/2.
        enddo

         i=1
         j=1
          dhdt(i,j)=(qx(i+1,j)-qx(i,j))/dx/1.
     &                 +(qy(i,j+1)-qy(i,j))/dy/1.

         i=Nx_Sedi
         j=1
          dhdt(i,j)=(qx(i,j)-qx(i-1,j))/dx/1.
     &                 +(qy(i,j+1)-qy(i,j))/dy/1.

         i=1
         j=Ny_Sedi
          dhdt(i,j)=(qx(i+1,j)-qx(i,j))/dx/1.
     &                 +(qy(i,j)-qy(i,j-1))/dy/1.

         i=Nx_Sedi
         j=Ny_Sedi
          dhdt(i,j)=(qx(i,j)-qx(i-1,j))/dx/1.
     &                 +(qy(i,j)-qy(i,j-1))/dy/1.


c  update depth
         do i = 1, Nx_Sedi
           do j = 1, Ny_Sedi
             ht(i,j)  = ht(i,j) + N_Interval_CallSedi*Master_dt
     &           *dhdt(i,j) !turned on
             Pass_Dupdated(i,j)=ht(i,j)
             Depth_Sedi(i,j)=ht(i,j)
           end do
          end do 

         open(1,file='dhdt.out')
        do i = 1, Nx_Sedi
         write(1,100)(dhdt(i,j), j=1,Ny_Sedi )
        end do
        close(1)

        open(1,file='depth.out')
        do j = 1, Ny_Sedi
         write(1,100)(depth_circ(i,j), i=1,Nx_Sedi )
        end do
        close(1)

  100  format(201(f20.12))


      return
      end


        subroutine Sedflux(ntime)

        use pass

         real por, qx(Nx_Max,0:Ny_Max+1), qy(Nx_Max,0:Ny_Max+1),Wo,
     .        c1,c2,qxum(Nx_Max,Ny_Max),qyum(Nx_Max,Ny_Max),
     .        qxslp(Nx_Max,Ny_Max),qyslp(Nx_Max,Ny_Max),
     .        qxsk(Nx_Max,Ny_Max),qysk(Nx_Max,Ny_Max),fw,pie,
     .        uw(0:1000),uw2, uw2sum,uw3,uw3sum,dhodx(Nx_Max,Ny_Max),
     .        dhody(Nx_Max,Ny_Max),fake(Nx_Max,Ny_Max),Tefuwxsum,
     .        fake1(Nx_Max,Ny_Max),fake2(Nx_Max,Ny_Max),Tefuwysum,
     .        skew(Nx_Max,Ny_Max),uw22(Nx_Max,Ny_Max),magu2,Tefsum,
     .        dhdt(Nx_max,Ny_max),ht(Nx_max,Ny_max),dx,dy,taueff,
     .        qxcr(Nx_Max,Ny_Max),qycr(Nx_Max,Ny_Max),magq,magqcr,
     .        thetacrit,d,s,taucr(nx_max,ny_max),rho,Tef(nx_max,ny_max),
     .        Tefuwx(nx_max,ny_max),Tefuwy,uwsum,uwavg(nx_max,ny_max),
     .        magu2sum,magu2avg(nx_max,ny_max),taumag
        integer ntime

c	write(*,*)'compute sediment fluxes '

         c1 = 0.3
         c2 = 0.3
         por = 0.35
         Wo=0.022
         pie=4.*atan2(1.,1.)
         rho=1015.

c Calculate depth gradients

         do i = 2, Nx_Circ-1
          do j = 2, Ny_Circ-1
             ht(i,j)=Depth_Circ(i,j)
             dx=X_Circ(i-1,j)-X_Circ(i+1,j)
             dy=Y_Circ(i,j+1)-Y_Circ(i+1,j-1)
             if(dx.eq.0..or.dy.eq.0.)then
               write(*,*) 'dx or dy = 0 in Circ sedflux, stop!'
               stop
             endif
             dhodx(i,j)=(ht(i+1,j)-ht(i-1,j))/dx
             dhody(i,j)=(ht(i,j+1)-ht(i,j-1))/dy
          end do
         end do

         do j=2,Ny_Circ-1
             dhodx(1,j)=dhodx(2,j)
             dhodx(Nx_Circ,j)=dhodx(Nx_Circ-1,j)
             dhody(1,j)=dhody(2,j)
             dhody(Nx_Circ,j)=dhody(Nx_Circ-1,j)
         enddo

         do i=1,Nx_Circ
             dhodx(i,1)=dhodx(i,2)
             dhodx(i,Ny_Circ)=dhodx(i,Ny_Circ-1)
             dhody(i,1)=dhody(i,2)
             dhody(i,Ny_Circ)=dhody(i,Ny_Circ-1)
         enddo


c  Calculate critical shear stress
          call taucrit(taucr,rho) 

         do i = 1, Nx_Circ
          do j = 1, Ny_Circ

         magu2sum=0
         uwsum=0
         uw2sum=0
         uw3sum=0
         Tefsum=0
         Tefuwxsum=0
         Tefuwysum=0

c Compute magnitude of u^2 and Tau and then find Tau_eff                  
         do l=1,100
            uw(l)=pass_uw(i,j,l)
            uwsum=uwsum+uw(l)
            uw2sum=uw2sum+uw(l)*uw(l)
            uw3sum=uw3sum+uw(l)*uw(l)*uw(l)
            magu2=uw(l)*uw(l)+Pass_Ub(i,j)*Pass_Ub(i,j)
     .           +2*uw(l)*(Pass_Ub(i,j)*cos(Pass_theta(i,j)*pie/180)
     .           +Pass_Vb(i,j)*sin(Pass_theta(i,j)*pie/180))
     .           +Pass_Vb(i,j)*Pass_Vb(i,j)
            magu2sum=magu2sum+magu2
            taumag=0.5d0*rho*magu2*Pass_fw(i,j)

            taueff=taumag-taucr(i,j)
            if (taueff.lt.0d0) taueff=0.d0
            Tefsum=Tefsum+taueff
           Tefuwxsum=Tefuwxsum+taueff*uw(l)*cos(Pass_theta(i,j)*pie/180)
           Tefuwysum=Tefuwysum+taueff*uw(l)*sin(Pass_theta(i,j)*pie/180)
         end do

c Time averages of Tau_eff and Tau_eff*uw
         magu2avg(i,j)=magu2sum/100
         Tef(i,j)=Tefsum/100
         Tefuwx(i,j)=Tefuwxsum/100
         Tefuwy=Tefuwysum/100
         uw2=uw2sum/100
         uw3=uw3sum/100
         uwavg(i,j)=uwsum/100
c         if(uw3.le.0) uw3=0
         skew(i,j)=uw3/uw2**1.5
 
c Compute cross-shore transport
           qxum(i,j)=2.*c1/9.8/rho*Tef(i,j)*Pass_Ub(i,j)
           qxslp(i,j)=2.*c2*c1*Wo*dhodx(i,j)/9.8/rho*Tef(i,j)
           qxsk(i,j)= 2.*c1/9.8/rho*Tefuwx(i,j)

           qx(i,j) = qxum(i,j)
c            qx(i,j) = 
     .          +qxslp(i,j)
     .            +qxsk(i,j)
           pass_sedflux_x(i,j) = qx(i,j)

c Compute longshore transport
           qyum(i,j)=2.*c1/9.8/rho*Tef(i,j)*Pass_Vb(i,j)
           qyslp(i,j)=2.*c2*c1*Wo*dhody(i,j)/9.8/rho*Tef(i,j)
           qysk(i,j)= 2.*c1/9.8/rho*Tefuwy

           qy(i,j) = qyum(i,j)
c           qy(i,j) =
     .          +qyslp(i,j)
     .            +qysk(i,j)
           pass_sedflux_y(i,j) = qy(i,j)

           end do
         end do

c treat pass_sedfluxcum as a time average flux

 	do j=1,ny_circ
   	  do i=1,nx_circ
           pass_sedfluxcum_x(i,j) = pass_sedfluxcum_x(i,j)
     &           +pass_sedflux_x(i,j)/ntime
           pass_sedfluxcum_y(i,j) = pass_sedfluxcum_y(i,j)
     &           +pass_sedflux_y(i,j)/ntime
           end do
          end do
 
        goto 7432
        open(1,file='qxum.out')
        do i = 1, Nx_Sedi
         write(1,100)(qxum(i,j), j=1,Ny_Sedi )
        end do
        close(1)
           open(1,file='qyum.out')
        do i = 1, Nx_Sedi
         write(1,100)(qyum(i,j), j=1,Ny_Sedi )
        end do
        close(1)
            open(1,file='qxslp.out')
        do i = 1, Nx_Sedi
         write(1,100)(qxslp(i,j), j=1,Ny_Sedi )
        end do
        close(1)
           open(1,file='qyslp.out')
        do i = 1, Nx_Sedi
         write(1,100)(qyslp(i,j), j=1,Ny_Sedi )
        end do
        close(1)
           open(1,file='qxsk.out')
        do i = 1, Nx_Sedi
         write(1,100)(qxsk(i,j), j=1,Ny_Sedi )
c         write(*,*)(qxsk(i,j), j=1,Ny_Sedi )
        end do
        close(1)
           open(1,file='qysk.out')
        do i = 1, Nx_Sedi
         write(1,100)(qysk(i,j), j=1,Ny_Sedi )
        end do
        close(1)
          open(1,file='mqy.out')
        do i = 1, Nx_Sedi
         write(1,100)(qy(i,j), j=1,Ny_Sedi )
        end do
        close(1)
          open(1,file='mqx.out')
        do i = 1, Nx_Sedi
         write(1,100)(qx(i,j), j=1,Ny_Sedi )
        end do
        close(1)

7432    continue


100     format(500f16.8)



         return
         end

c ----------------------------------------------

       subroutine taucrit(taucr,rho) 
!       include 'pass.h'
       use pass

       real taucr(nx_max,ny_max),rho,d,s,thetacrit

        thetacrit=.05
        d=0.0002
        s=2.65
         do i = 1, Nx_Circ
          do j = 1, Ny_Circ
             taucr(i,j)=rho*(s-1)*9.8*d*thetacrit
          end do
        end do

        return
        end
