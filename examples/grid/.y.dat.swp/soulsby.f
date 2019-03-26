        program try_soulsby
        parameter(m=309,n=288)
        real Pass_Ub(m,n),Pass_Vb(m,n),Pass_ubott(m,n),Pass_theta(m,n)
     &       ,Pass_SedFlux_x(m,n),Pass_SedFlux_y(m,n)
     &       ,Pass_MassFluxU(m,n),Pass_MassFluxV(m,n)
     &       ,Pass_U(m,n),Pass_V(m,n)
     &       ,depth(m,n)

        integer nm_first,nm_second,nm_third,nm_fourth,num_file
        character*4 file_name
        character*28 fdir

        fdir='../wave-induced-current_210/'

        print*,'finename=?'
        read(*,*)num_file


        nm_first=mod(num_file/1000,10)
        nm_second=mod(num_file/100,10)
        nm_third=mod(num_file/10,10)
        nm_fourth=mod(num_file,10)


        write(file_name(1:1),'(I1)')nm_first
        write(file_name(2:2),'(I1)')nm_second
        write(file_name(3:3),'(I1)')nm_third
        write(file_name(4:4),'(I1)')nm_fourth


        do j=1,n
        do i=1,m
         Pass_Ub(i,j)=0.0
         Pass_Vb(i,j)=1.
         Pass_ubott(i,j)=0.9
         Pass_theta(i,j)=0.
        enddo
        enddo

      open(2,file=fdir//'wx'//file_name//'.out')
          do j=1,n
            read(2,*)(Pass_MassFluxU(i,j),
     &                        i=1,m)
          enddo
       close(2)

       open(2,file=fdir//'wy'//file_name//'.out')
          do j=1,n
            read(2,*)(pass_massfluxV(i,j),
     &                        i=1,m)
          enddo
        close(2)

      open(2,file=fdir//'cu'//file_name//'.out')
          do j=1,n
            read(2,*)(Pass_U(i,j),
     &                        i=1,m)
          enddo
       close(2)

       open(2,file=fdir//'cv'//file_name//'.out')
          do j=1,n
            read(2,*)(pass_V(i,j),
     &                        i=1,m)
          enddo
        close(2)

      open(2,file=fdir//'ub'//file_name//'.out')
          do j=1,n
            read(2,*)(Pass_Ub(i,j),
     &                        i=1,m)
          enddo
       close(2)

       open(2,file=fdir//'vb'//file_name//'.out')
          do j=1,n
            read(2,*)(pass_Vb(i,j),
     &                        i=1,m)
          enddo
        close(2)

       open(2,file=fdir//'wb'//file_name//'.out')
          do j=1,n
            read(2,*)(pass_ubott(i,j),
     &                        i=1,m)
          enddo
        close(2)


       open(2,file=fdir//'ag'//file_name//'.out')
          do j=1,n
            read(2,*)(pass_theta(i,j),
     &                        i=1,m)
          enddo
        close(2)

      open(2,file=fdir//'dp'//file_name//'.out')
          do j=1,n
            read(2,*)(depth(i,j),
     &                        i=1,m)
          enddo
        close(2)


        do j=1,n
        do i=1,m
          if(depth(i,j).ne.0)then
          Pass_Ub(i,j)=Pass_U(i,j)-Pass_MassFluxU(i,j)/depth(i,j)
          Pass_Vb(i,j)=Pass_V(i,j)-Pass_MassFluxV(i,j)/depth(i,j)
c          Pass_Ub(i,j)=-Pass_MassFluxU(i,j)/depth(i,j)
c          Pass_Vb(i,j)=-Pass_MassFluxV(i,j)/depth(i,j)
          endif
        enddo
        enddo


        call soulsby(Pass_SedFlux_x,Pass_SedFlux_y,m,n,
     &           Pass_Ub,Pass_Vb,Pass_ubott,Pass_theta)


        print*,Pass_SedFlux_x(1,1),Pass_SedFlux_y(1,1)

        open(1,file='flux_x.out')
         do j=1,n
         write(1,111)(pass_sedflux_x(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='flux_y.out')
         do j=1,n
         write(1,111)(pass_sedflux_y(i,j),i=1,m)
         enddo
        close(1)

111     format(500E12.2)

        end

        subroutine soulsby(SedFlux_x,SedFlux_y,m,n,
     &           Ub,Vb,ubott,phi_w)
        implicit none
        integer m,n,i,j
        real Ub(m,n),Vb(m,n),ubott(m,n),phi_w(m,n)
     &       ,SedFlux_x(m,n),SedFlux_y(m,n)
     &       ,phi,phi_c,deg2rad,tau_c,tau_w,tau_m
     &       ,tau_max,Psi_x1,Psi_x2,Psi_x,Psi_y
     &       ,rho,Cd,fw,grav,d,pi,q1,q2
     &       ,th_m,th_w,th_max,th_max1,th_max2,th_cr,delta,A2,s

        pi=4.*atan2(1.,1.)
        deg2rad=pi/180.
        delta=0.2
        rho=1027.
        Cd=0.0023
        fw=0.0316
        grav=9.8
        A2=12
        th_cr=0.05
        s=2.65
        d=0.0002

        do 100 j=1,n
        do 100 i=1,m
         phi_c=atan2(Vb(i,j),Ub(i,j))
         phi=phi_w(i,j)*deg2rad-phi_c
         tau_c=rho*Cd*(Ub(i,j)*Ub(i,j)+Vb(i,j)*Vb(i,j))
         tau_w=0.5*rho*fw*(ubott(i,j)*ubott(i,j))
         tau_m=tau_c*(1.+1.2*(tau_w/(tau_c+tau_w))**3.2)
         tau_max=sqrt((tau_m+tau_w*cos(phi))**2+(tau_w*sin(phi))**2)

         th_m=tau_m/rho/grav/(s-1.)/d
         th_w=tau_w/rho/grav/(s-1.)/d

         th_max1=sqrt((th_m+th_w*(1+delta)*cos(phi))**2+
     &                (th_w*(1+delta)*sin(phi))**2)
         th_max2=sqrt((th_m+th_w*(1-delta)*cos(phi+pi))**2+
     &                (th_w*(1-delta)*sin(phi+pi))**2)
         th_max=max(th_max1,th_max2)

         Psi_x1=A2*sqrt(th_m)*(th_m-th_cr)
         Psi_x2=A2*(0.9534+0.1907*cos(2.*phi))*sqrt(th_w)*th_m
     &         +A2*(0.229*delta*th_w**1.5*cos(phi))
         Psi_x=max(Psi_x1,Psi_x2)

         Psi_y=A2*(0.1907*th_w**2)/(th_w**1.5 +1.5*th_m**1.5)
     &         *(th_m*sin(2.*phi) 
     &         +1.2*delta*th_w*sin(phi))

         if(th_max.le.th_cr)then
          Psi_x=0.
          Psi_y=0.
         endif

         q1=Psi_x*sqrt(grav*(s-1.)*d**3)
         q2=Psi_y*sqrt(grav*(s-1.)*d**3)

         SedFlux_x(i,j)=q1*cos(phi_c)-q2*sin(phi_c)
         SedFlux_y(i,j)=q1*sin(phi_c)+q2*cos(phi_c)

100     continue

        return
        end




