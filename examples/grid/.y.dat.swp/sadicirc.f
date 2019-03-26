c-----------------------------------------------------------
c   This is the Semi-implicit Mode of the circulation module
c used in the nearshore community model (NOPP model). For details,
c refer: Shi F. and Kirby J. T., 2004, Efficient semi-implicit numerical
c schemes for quasi-3D nearshore circulation model with non-orthogonal curvilinear grid,   
c in preparation.
c 02/17/2004
c ------------------------------------------------------------

	Subroutine CircModule()
        USE sadi
        real ddy(m,n),thetay(m,n)
  
! ---  initialize ...

      if (Master_Start.eq.1)then	

	print*,'initialize circulation module ...'

	call init_shorecirc

        call Jacobi

        call alphaetal

        call adu

	call get_gxx_Cf

      else	

	print*,'circulation module ... '

! time integration

	call adaptor_wave_to_circ
        call adaptor_sedi_to_circ
        call initialize_sedflux
c -- time

	kplot=max(int(plotintra/dt),1)

c -- time integration

	ik=0

        do 100 itstep=1,ntime

	  time_circ=Time_Master+(itstep-1)*dt

c --- wind , minute, m/s degree
        
        if(windtime(1).ge.0.)then
	     wind_intv=wind_intv+dt
	    if((windtime(2)-windtime(1))*3600..le.wind_intv)then
	      windtime(1)=windtime(2)
	      windspd(1)=windspd(2)
	      winddir(1)=winddir(2)
           read(7,*,end=981)windtime(2),windspd(2),winddir(2)
            wind_intv=0.
981         if(windtime(2).eq.windtime(1))then
              write(*,*)'no more wind data, use the last read'
	        windtime(2)=windtime(1)+ time_circ*100.
	        wind_intv=0.
	      endif
	    endif
	  endif 
	  
	  call calculate_wind    

        call openb(time_circ)

	  call convertUVZ	

	  call diffusion  

          if(Wave_To_Circ_Radiation)then
          call radiation_stress
          endif

          call solve_zeta_u

          call solve_zeta_v

          call velocity

          call convert_Tau

        if (k_3d.gt.1)then
          call solve_3d
	  endif

        if (disp3d.eq.1)then
          call dispersion
        endif

        if(mod(itstep,kplot).eq.0)then

	print*,'circ-intr-time=',time_circ, 
     &          itstep,ntime,z(int(m/2),int(n/2))

	if(time_circ.ge.plotstart)then
          ik=ik+1
	print*,'plotting...',ik

        call outputcirc(m,n,ik,z)
	call outputcirc(m,n,ik+1000,u)
	call outputcirc(m,n,ik+2000,v)
	call outputcirc(m,n,ik+3000,hh)
	call outputcirc(m,n,ik+4000,pass_sedflux_y)

	endif
        endif

        if(time_circ.ge.time_sta)then
           time_sta=time_sta + sta_interval

	print*,'t=',time_circ, '   output station z u v'

        call output_station(time_circ,n_station,i_sta,j_sta)

        endif

c --- sediment flux

      if(N_Interval_CallSedi.ge.0)then
       call adaptor_circ_to_sedi     
       call Sedflux(ntime)
	endif

100     continue

c -- analysis
        do j=1,n
	do i=1,m
	 if(Ntype(i,j).eq.0)then
c         PREX(i,j)=(z(i,j)-z(i-1,j))/dxi*Yeta(i,j)/sq_g_0(i,j)
c     &            -(z(i,j)-z(i,j-1))/deta*Yxi(i,j)/sq_g_0(i,j)
c         PREY(i,j)=-(z(i,j)-z(i-1,j))/dxi*Xeta(i,j)/sq_g_0(i,j)
c     &            +(z(i,j)-z(i,j-1))/deta*Xxi(i,j)/sq_g_0(i,j)

c         ddy(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
c         thetay(i,j)=atan2(abs(x(i,j)-x(i,j-1)),
c     &                     abs(y(i,j)-y(i,j-1)))
c         PREY(i,j)=(z(i,j)-z(i,j-1))/ddy(i,j)*9.8*
c     &       (depth_circ(i,j)+z(i,j))
c         RADY(i,j)=Intp_Fy_Circ(i,j)*cos(thetay(i,j))
c     &            -Intp_Fx_Circ(i,j)*sin(thetay(i,j))
         endif
	enddo
	enddo

	call adaptor_circ_to_wave

	endif       

c -- time integration over

	return

        end

c ---------------------------------
c       origin
c ---------------------------------
        subroutine init_shorecirc
        USE SADI
	real depmax, spacingmin,disx,disy
        integer kctr,Np
	character*255 Circ_type, Circ_specify_flux,Circ_mask,
     &                Circ_station,Circ_wind

	namelist /fnames/ Circ_type, Circ_specify_flux,Circ_mask,
     &                    Circ_station,Circ_wind
     &           /numerics/ CR,curv_grid,orth_grid,plotintra,
     &                      plotstart,sta_interval
     &           /grid/ k_3d
     &           /boundary/
     &          east_ele2_flx3_grd4_rad5_per6,jstart_E,jend_E,
     &          west_ele2_flx3_grd4_rad5_per6,jstart_W,jend_W,
     &          south_ele2_flx3_grd4_rad5_per6,istart_S,iend_S,
     &          north_ele2_flx3_grd4_rad5_per6,istart_N,iend_N,
     &          east_data1_anly2,west_data1_anly2,
     &          south_data1_anly2,north_data1_anly2
     &                     
     &           /physics/ f_cwc, rmanning, hi, delta,
     &                     c_subgrid, anu_tb_const, anu_t0_const,
     &                     disp3d,
     &                     depth_min,wind_u,wind_v


	open(1,file='curvcircinput.dat')
	  read(1,nml=fnames)
	  read(1,nml=numerics)
          read(1,nml=grid)
	  read(1,nml=boundary)
	  read(1,nml=physics)
	close(1)


	plotstart = 0

        f=3.1416*SIN(3.1416*39.0/180.0)/21600.
c	f=0.


	m=nx_circ
	n=ny_circ

! ----------- depth boundary ...

! boundary
	if (Circ_type .ne. ' ')then
      
        open(1,file=Circ_type)
          do j=n,1,-1
             read(1,*)(Ntype(i,j),i=1,m)
          enddo

        close(1)

	else

	call gridtype

        endif

        call newgridtype


c -- x y coordinates

            do j=1,n
            do i=1,m
	       x(i,j)=X_Circ(i,j)
	       y(i,j)=Y_Circ(i,j)
            enddo
            enddo 

	dx=x(2,1)-x(1,1)
	dy=y(1,2)-y(1,1)

c --- station

         if(Circ_station.ne.' ')then
           open(1,file=Circ_station)
             read(1,*)n_station
             do i=1,n_station
             read(1,*,end=999)i_sta(i),j_sta(i)
             enddo
             go to 998
999          write(*,*)'station number wrong!!!'
             stop
998          continue
           close(1)
           open(83,file='outputz.out')
	     open(84,file='outputuv.out')

         endif

c ---  wind data
         if(Circ_wind.ne.' ')then
           open(7,file=Circ_wind)
             read(7,*) windtime(1),windspd(1),winddir(1)
             read(7,*) windtime(2),windspd(2),winddir(2)
             print*,windtime(1),windtime(2),windspd(1),windspd(2)
	   else
	     windtime(1)=-1.
         endif


c --- tide
         goto 1119
         open(1,file='amp_west.dat')
         do j=1,n
         read(1,200)(amp_west(j,i),i=1,10)
         enddo
         close(1)
         open(1,file='pha_west.dat')
         do j=1,n
         read(1,200)(pha_west(j,i),i=1,10)
         enddo
         close(1)

         open(1,file='amp_south.dat')
         do j=istart_S,iend_S
         read(1,200)(amp_south(j,i),i=1,10)
         enddo
         close(1)
         open(1,file='pha_south.dat')
         do j=istart_S,iend_S
         read(1,200)(pha_south(j,i),i=1,10)
         enddo
         close(1)

         open(1,file='amp_north.dat')
         do j=istart_N,iend_N
         read(1,200)(amp_north(j,i),i=1,10)
         enddo
         close(1)
         open(1,file='pha_north.dat')
         do j=istart_N,iend_N
         read(1,200)(pha_north(j,i),i=1,10)
         enddo
         close(1)

1119     continue
200      format(100f12.6)

c -- depth

	do j=1,n
	do i=1,m	
	depth_c(i,j)=depth_circ(i,j)
	enddo
	enddo

c --- mask
	if (Circ_mask .ne. ' ')then
      
        open(1,file=Circ_mask)
          do j=n,1,-1
             read(1,*)(Mask(i,j),i=1,m)
          enddo
        close(1)
        else
          do j=1,n
          do i=1,m
             Mask(i,j)=1
          enddo
          enddo

        endif

        do j=1,n
          do i=1,m
           if(Mask(i,j).eq.0)then
           depth_c(i,j)=-99.
           endif
          enddo
        enddo
c ---

        call depth_u_v_H

c -- set Z

        do j=1,n
        do i=1,m
          Z(i,j)=0.
          if(H(i,j).lt.0.)Z(i,j)=-H(i,j)-delta
          HH(i,j)=H(i,j)+Z(i,j)
        enddo
        enddo

c --- dt

	depmax=0.
	spacingmin=10000000.
	do j=1,n-1
	do i=1,m-1
	  if(ntype(i,j).eq.0)then
	  if(H(i,j).gt.depmax)depmax=H(i,j)
	  disx=sqrt((x(i+1,j)-x(i,j))**2+(y(i+1,j)-y(i,j))**2)/2.
	  disy=sqrt((x(i,j+1)-x(i,j))**2+(y(i,j+1)-y(i,j))**2)/2.
	  if(disx.lt.spacingmin)spacingmin=disx
	  if(disy.lt.spacingmin)spacingmin=disy
	  if(spacingmin.eq.0)print*,i,j,ntype(i,j)
	  endif
	enddo
	enddo

         if(depmax.eq.0.)then
          write(*,*)'depth input wrong! depmax=0. stop!'
          stop
         endif

	dt=CR*spacingmin/sqrt(9.8*depmax)

         if(dt.eq.0.)then
          write(*,*)'spacingmin=0 or Cr =0, stop!'
          stop
         endif        

	ntime=int(Master_dt*N_interval_CallCirc/dt)

         if(ntime.eq.0.)then
          write(*,*)'dt in circulation module larger than master_dt'
          write(*,*)'stop'
          stop
         endif     

        dt=Master_dt*N_interval_CallCirc/ntime

	print*,'dt=',dt,'    ntime=',ntime


c -- initialize variables ...

        do j=1,n
        do i=1,m
          UU(i,j)=0.
          VV(i,j)=0.
          U(i,j)=0.
          V(i,j)=0.
          Tau_1(i,j)=0.
          Tau_2(i,j)=0.
          diffx(i,j)=0.
          diffy(i,j)=0.
          dispx(i,j)=0.
          dispy(i,j)=0.
        enddo
        enddo

        Tao_ax=0.00
        Tao_ay=0.00


c -- ibeg,iend

	do j=1,n
	ibeg(j)=0
	iend(j)=0
	enddo

	do i=1,m
	jbeg(i)=0
	jend(i)=0
	enddo


      do j=1,n-1

        kctr=0
        do i=1,m-1
          if(kctr.eq.0)then
            Np=Ntype1(i,j)
            if(Np.ge.30.and.Np.le.80)Np=20+mod(Np,10)
            if(Np.eq.14.or.Np.eq.24.or.Np.eq.15.or.Np.eq.25.or.Np.eq.18
     *         .or.Np.eq.28)then
                ibeg(j)=i
                kctr=1
            endif
          endif
        enddo

        kctr=0
        do i=m-1,1,-1
          if(kctr.eq.0)then
            Np=Ntype1(i,j)
            if(Np.ge.30.and.Np.le.80)Np=20+mod(Np,10)
            if(Np.eq.12.or.Np.eq.22.or.Np.eq.16.or.Np.eq.26.or.Np.eq.17
     *         .or.Np.eq.27)then
                iend(j)=i
                kctr=1
            endif
          endif
        enddo

      enddo

c -- ok

c -- jbeg,jend

      do i=1,m-1

        kctr=0
        do j=1,n-1
          if(kctr.eq.0)then
            Np=Ntype1(i,j)
            if(Np.ge.30.and.Np.le.80)Np=20+mod(Np,10)
            if(Np.eq.11.or.Np.eq.21.or.Np.eq.15.or.Np.eq.25.or.Np.eq.16
     *         .or.Np.eq.26)then
                jbeg(i)=j
                kctr=1
            endif
          endif
        enddo

        kctr=0
        do j=n-1,1,-1
          if(kctr.eq.0)then
            Np=Ntype1(i,j)
            if(Np.ge.30.and.Np.le.80)Np=20+mod(Np,10)
            if(Np.eq.13.or.Np.eq.23.or.Np.eq.18.or.Np.eq.28.or.Np.eq.17
     *         .or.Np.eq.27)then
                jend(i)=j
                kctr=1
            endif
          endif
        enddo

      enddo

c -- ok

        return
        end

        subroutine adaptor_sedi_to_circ
!        use sadi
        use sadi 

c -- depth

	do j=1,n
	do i=1,m	
	depth_c(i,j)=depth_circ(i,j)
c --- tmp added to test instability
c          if(depth_c(i,j).lt.depth_min)depth_c(i,j)=depth_min  
	enddo
	enddo

        do j=1,n
          do i=1,m
           if(Mask(i,j).eq.0)then
           depth_c(i,j)=-99.
           endif
          enddo
        enddo
c ---

        call depth_u_v_H

c -- set Z

        do j=1,n
        do i=1,m       
          if(H(i,j).lt.-Z(i,j))Z(i,j)=-H(i,j)-delta
          HH(i,j)=H(i,j)+Z(i,j)
        enddo
        enddo


        return
        end

c ---------------------------------
c       Jacobi
c ---------------------------------
        subroutine Jacobi
!        use sadi
        use sadi
        do j=1,n-1
        do i=1,m-1
          if(Ntype(i,j).lt.80)then
            rJ(i,j)=1./(dxi*deta)
     *              *(0.5*(X(i+1,j)+X(i+1,j+1))-0.5*(X(i,j)+X(i,j+1)))
     *              *(0.5*(Y(i,j+1)+Y(i+1,j+1))-0.5*(Y(i,j)+Y(i+1,j)))
     *             -1./(dxi*deta)
     *              *(0.5*(X(i,j+1)+X(i+1,j+1))-0.5*(X(i,j)+X(i+1,j)))
     *              *(0.5*(Y(i+1,j)+Y(i+1,j+1))-0.5*(Y(i,j)+Y(i,j+1)))
          endif
        enddo
        enddo
        return
        end

c ---------------------------------
c       alpha 
c ---------------------------------
        subroutine Alphaetal
!        use sadi
        use sadi

c -- alpha
        do j=1,n-1
        do i=1,m-1
          if(Ntype(i,j).lt.80)then
            alpha(i,j)=1./(deta**2)
     *            *(0.5*(X(i,j+1)+X(i+1,j+1))-0.5*(X(i,j)+X(i+1,j)))**2
     *                +1./(deta**2)
     *            *(0.5*(Y(i,j+1)+Y(i+1,j+1))-0.5*(Y(i,j)+Y(i+1,j)))**2
          endif
        enddo
        enddo
c -- gamma
        do j=1,n-1
        do i=1,m-1
          if(Ntype(i,j).lt.80)then
            gamma(i,j)=1./(dxi**2)
     *            *(0.5*(X(i+1,j)+X(i+1,j+1))-0.5*(X(i,j)+X(i,j+1)))**2
     *                +1./(dxi**2)
     *            *(0.5*(Y(i+1,j)+Y(i+1,j+1))-0.5*(Y(i,j)+Y(i,j+1)))**2

          endif
        enddo
        enddo
c -- beta
        do j=1,n-1
        do i=1,m-1
	  if(orth_grid.eq.0)then
          if(Ntype(i,j).lt.80)then
            beta(i,j)=1./(dxi*deta)
     *              *(0.5*(X(i+1,j)+X(i+1,j+1))-0.5*(X(i,j)+X(i,j+1)))
     *              *(0.5*(X(i,j+1)+X(i+1,j+1))-0.5*(X(i,j)+X(i+1,j)))
     *               +1./(dxi*deta)
     *              *(0.5*(Y(i+1,j)+Y(i+1,j+1))-0.5*(Y(i,j)+Y(i,j+1)))
     *              *(0.5*(Y(i,j+1)+Y(i+1,j+1))-0.5*(Y(i,j)+Y(i+1,j)))
          endif
	  else
	    beta(i,j)=0.
          endif
        enddo
        enddo
c -- ok
        return
        end

c ---------------------------------
c       adu 
c ---------------------------------
        subroutine adu
!        use sadi
        use sadi
        integer Np,NP1
        real X_xi_eta,Y_x, Y_xi_eta,X_xi,X_eta,Y_eta,X_eta_eta,
     *       Y_eta_eta,X_xi_xi,Y_xi_xi, Y_xi
c -- adu1
        do j=1,n-1
        do i=1,m-1
          Np=Ntype(i,j)
          if(Np.ge.30.and.Np.le.80)Np=20+mod(Np,10)
          if(Np.eq.0.or.Np.eq.12.or.Np.eq.13.or.Np.eq.17
     *      .or.Np.eq.22.or.Np.eq.23.or.Np.eq.27) then
            X_xi_eta=1./(4.*dxi*deta)
     *              *((X(i+1,j+1)-X(i-1,j+1))-(X(i+1,j-1)-X(i-1,j-1)))
            Y_xi=1./(2.*dxi)*(Y(i+1,j)-Y(i-1,j))
            Y_xi_eta=1./(4.*dxi*deta)
     *              *((Y(i+1,j+1)-Y(i-1,j+1))-(Y(i+1,j-1)-Y(i-1,j-1)))
            X_xi=1./(2.*dxi)*(X(i+1,j)-X(i-1,j))
            X_eta=1./(2.*deta)*(X(i,j+1)-X(i,j-1))
            Y_eta=1./(2.*deta)*(Y(i,j+1)-Y(i,j-1))
            X_eta_eta=1./(deta**2)
     *                *(X(i,j+1)-2.*X(i,j)+X(i,j-1))
            Y_eta_eta=1./(deta**2)
     *                *(Y(i,j+1)-2.*Y(i,j)+Y(i,j-1))
            X_xi_xi=1./(dxi**2)
     *             *(X(i+1,j)-2.*X(i,j)+X(i-1,j))
            Y_xi_xi=1./(dxi**2)
     *             *(Y(i+1,j)-2.*Y(i,j)+Y(i-1,j))
          endif

          if(Np.eq.11.or.Np.eq.16
     *       .or.Np.eq.21.or.Np.eq.26)then
            X_xi_eta=1./(2.*dxi*deta)
     *              *((X(i+1,j+1)-X(i-1,j+1))-(X(i+1,j)-X(i-1,j)))
            Y_xi=1./(2.*dxi)*(Y(i+1,j)-Y(i-1,j))
            Y_xi_eta=1./(2.*dxi*deta)
     *              *((Y(i+1,j+1)-Y(i-1,j+1))-(Y(i+1,j)-Y(i-1,j)))
            X_xi=1./(2.*dxi)*(X(i+1,j)-X(i-1,j))
            X_eta=1./(1.*deta)*(X(i,j+1)-X(i,j))
            Y_eta=1./(1.*deta)*(Y(i,j+1)-Y(i,j))
            X_eta_eta=1./(deta**2)
     *                *(X(i,j+2)-2.*X(i,j+1)+X(i,j))
            Y_eta_eta=1./(deta**2)
     *                *(Y(i,j+2)-2.*Y(i,j+1)+Y(i,j))
            X_xi_xi=1./(dxi**2)
     *             *(X(i+1,j)-2.*X(i,j)+X(i-1,j))
            Y_xi_xi=1./(dxi**2)
     *             *(Y(i+1,j)-2.*Y(i,j)+Y(i-1,j))
          endif

          if(Np.eq.14.or.Np.eq.18
     *      .or.Np.eq.24.or.Np.eq.28)then
            X_xi_eta=1./(2.*dxi*deta)
     *              *((X(i+1,j+1)-X(i,j+1))-(X(i+1,j-1)-X(i,j-1)))
            Y_xi=1./(1.*dxi)*(Y(i+1,j)-Y(i,j))
            Y_xi_eta=1./(2.*dxi*deta)
     *              *((Y(i+1,j+1)-Y(i,j+1))-(Y(i+1,j-1)-Y(i,j-1)))
            X_xi=1./(1.*dxi)*(X(i+1,j)-X(i,j))
            X_eta=1./(2.*deta)*(X(i,j+1)-X(i,j-1))
            Y_eta=1./(2.*deta)*(Y(i,j+1)-Y(i,j-1))
            X_eta_eta=1./(deta**2)
     *                *(X(i,j+1)-2.*X(i,j)+X(i,j-1))
            Y_eta_eta=1./(deta**2)
     *                *(Y(i,j+1)-2.*Y(i,j)+Y(i,j-1))
            X_xi_xi=1./(dxi**2)
     *             *(X(i+2,j)-2.*X(i+1,j)+X(i,j))
            Y_xi_xi=1./(dxi**2)
     *             *(Y(i+2,j)-2.*Y(i+1,j)+Y(i,j))
          endif

          if(Np.eq.15.or.Np.eq.25)then
            X_xi_eta=1./(1.*dxi*deta)
     *              *((X(i+1,j+1)-X(i,j+1))-(X(i+1,j)-X(i,j)))
            Y_xi=1./(1.*dxi)*(Y(i+1,j)-Y(i,j))
            Y_xi_eta=1./(1.*dxi*deta)
     *              *((Y(i+1,j+1)-Y(i,j+1))-(Y(i+1,j)-Y(i,j)))
            X_xi=1./(1.*dxi)*(X(i+1,j)-X(i,j))
            X_eta=1./(1.*deta)*(X(i,j+1)-X(i,j))
            Y_eta=1./(1.*deta)*(Y(i,j+1)-Y(i,j))
            X_eta_eta=1./(deta**2)
     *                *(X(i,j+2)-2.*X(i,j+1)+X(i,j))
            Y_eta_eta=1./(deta**2)
     *                *(Y(i,j+2)-2.*Y(i,j+1)+Y(i,j))
            X_xi_xi=1./(dxi**2)
     *             *(X(i+2,j)-2.*X(i+1,j)+X(i,j))
            Y_xi_xi=1./(dxi**2)
     *             *(Y(i+2,j)-2.*Y(i+1,j)+Y(i,j))
          endif

            adu1(i,j)=X_xi_eta*Y_xi-Y_xi_eta*X_xi
            adu2(i,j)=X_xi_eta*Y_eta-Y_xi_eta*X_eta+X_eta_eta*Y_xi
     *                -Y_eta_eta*X_xi
            adu3(i,j)=X_eta_eta*Y_eta-Y_eta_eta*X_eta
            adv1(i,j)=Y_xi_eta*X_eta-X_xi_eta*Y_eta
            adv2(i,j)=Y_xi_xi*X_eta-X_xi_xi*Y_eta+Y_xi_eta*X_xi
     *                -X_xi_eta*Y_xi
            adv3(i,j)=Y_xi_xi*X_xi-X_xi_xi*Y_xi
            Xxi(i,j)=X_xi
            Xeta(i,j)=X_eta
            Yxi(i,j)=Y_xi
            Yeta(i,j)=Y_eta

        enddo
        enddo

c	open(9,file='tmp1.out')
c       do j=1,n
c          write(9,199)(Xeta(i,j),i=1,m)
c	enddo
c	close(9)
c	open(9,file='tmp2.out')
c       do j=1,n
c          write(9,199)(Yeta(i,j),i=1,m)
c	enddo
c	close(9)
c	open(9,file='tmp3.out')
c       do j=1,n
c          write(9,199)(Xxi(i,j),i=1,m)
c	enddo
c	close(9)
c	stop

c199   format(500f16.3)

        return
        end

c ---------------------------------
c       zeta_u
c ---------------------------------
        subroutine solve_zeta_u
!        use sadi
        use sadi
        integer Np,Np1
        real Y_eta,X_eta
        real A(mtrigm),C(mtrigm),D(mtrigm)
        real Ea(mtrigm),Eb(mtrigm),Ztrig(mtrigm)
        real Zorig(mtrim)

        do 100 j=1,n-1
c -- begin row
          if(ibeg(j).eq.iend(j))goto 100

          do 200 i=ibeg(j),iend(j)

            Np=Ntype(i,j)
            Np1=Ntype(i,j)
            if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
c -- u_rJ 
            if(Np.eq.0.or.Np.eq.11.or.Np.eq.21.
     *        or.np.eq.12.or.Np.eq.22.or.Np.eq.16.or.Np.eq.26)then
              u_rj=0.5*(rJ(i-1,j)+rJ(i,j))
              u_H=0.5*(h(i-1,j)+h(i,j)+z(i-1,j)+z(i,j))
              if(u_H.lt.hi)u_H=hi
              u_v=0.25*(vv(i,j)+vv(i,j+1)+vv(i-1,j+1)+vv(i-1,j))
              u_beta=0.5*(beta(i-1,j)+beta(i,j))
              u_alpha=0.5*(alpha(i-1,j)+alpha(i,j))
              u_gamma=0.5*(gamma(i-1,j)+gamma(i,j))
              Y_eta=1./deta*(Y(i,j+1)-Y(i,j))
              X_eta=1./deta*(X(i,j+1)-X(i,j))
            endif

            if(Np.eq.13.or.Np.eq.23.or.Np.eq.17.or.Np.eq.27)then
              u_rj=0.5*(rJ(i-1,j)+rJ(i,j))
              u_H=0.5*(h(i-1,j)+h(i,j)+z(i-1,j)+z(i,j))
              if(u_H.lt.hi)u_H=hi
              u_v=0.25*(vv(i,j)+0.+0.+vv(i-1,j))
              u_beta=0.5*(beta(i-1,j)+beta(i,j))
              u_alpha=0.5*(alpha(i-1,j)+alpha(i,j))
              u_gamma=0.5*(gamma(i-1,j)+gamma(i,j))
              Y_eta=1./deta*(Y(i,j+1)-Y(i,j))
              X_eta=1./deta*(X(i,j+1)-X(i,j))
            endif

            if(Np.eq.14.or.Np.eq.24.or.Np.eq.15.or.np.eq.25)then
              u_rj=rJ(i,j)
              u_H=h(i,j)+z(i,j)
              if(u_H.lt.hi)u_H=hi
              u_v=0.5*(vv(i,j)+vv(i,j+1))
              u_beta=beta(i,j)
              u_alpha=alpha(i,j)
              u_gamma=gamma(i,j)
              Y_eta=1./deta*(Y(i,j+1)-Y(i,j))
              X_eta=1./deta*(X(i,j+1)-X(i,j))
            endif

            if(Np.eq.18.or.Np.eq.28)then
              u_rj=rJ(i,j)
              u_H=h(i,j)+z(i,j)
              if(u_H.lt.hi)u_H=hi
              u_v=0.5*(vv(i,j)+0.)
              u_beta=beta(i,j)
              u_alpha=alpha(i,j)
              u_gamma=gamma(i,j)
              Y_eta=1./deta*(Y(i,j+1)-Y(i,j))
              X_eta=1./deta*(X(i,j+1)-X(i,j))
            endif

c -- u_rj over
c --  \partial u**2/H
            if(Np.eq.0.or.Np.eq.11.or.Np.eq.21.or.Np.eq.13.or.Np.eq.23)
     *        then
              H_1=HH(i,j)
              H_2=HH(i-1,j)
              if(H_1.lt.hi)H_1=hi
              if(H_2.lt.hi)H_2=hi
              UUH_xi=1./(dxi)*(0.25*(UU(i+1,j)+UU(i,j))**2/H_1
     *             -0.25*(UU(i-1,j)+UU(i,j))**2/H_2)
            endif

            if(Np.eq.16.or.Np.eq.17.or.Np.eq.12
     *        .or.Np.eq.26.or.Np.eq.27.or.Np.eq.22)then
              H_1=HH(i,j)
              H_2=HH(i-1,j)
              if(H_1.lt.hi)H_1=hi
              if(H_2.lt.hi)H_2=hi
              UUH_xi=1./(dxi)*(0.25*(0.+UU(i,j))**2/H_1
     *             -0.25*(UU(i-1,j)+UU(i,j))**2/H_2)
            endif

            if(Np.eq.15.or.Np.eq.14.or.Np.eq.18
     *        .or.Np.eq.25.or.Np.eq.24.or.Np.eq.28)then
              H_1=HH(i,j)
              if(H_1.lt.hi)H_1=hi
              UUH_xi=1./(dxi/2.)*(0.25*(UU(i+1,j)+0.)**2/H_1
     *             -0.)
            endif

c ---   periodic
            if(Np1.eq.64)
     *        then
              H_1=HH(i,j)
              if(H_1.lt.hi)H_1=hi
              UUH_xi=1./(dxi)*(0.25*(UU(i+1,j)+UU(i,j))**2/H_1
     *             -0.25*(UU(iend(j),j)+UU(i,j))**2/H_1)
            endif

            if(Np1.eq.62)
     *        then
              H_1=HH(i,j)
              if(H_1.lt.hi)H_1=hi
              UUH_xi=1./(dxi)*(0.25*(UU(ibeg(j),j)+UU(i,j))**2/H_1
     *             -0.25*(UU(i-1,j)+UU(i,j))**2/H_1)
            endif

c --- elevation

            if(Np1.eq.25.or.Np1.eq.26.or.Np1.eq.21
     *        .or.Np1.eq.28.or.Np1.eq.27.or.Np1.eq.23)then
              UVH_xi=0.
            endif

c -- \partial uv/H
            if(Np.eq.0.or.Np.eq.12.or.Np.eq.22)then
              H_1=0.25*(HH(i-1,j+1)+HH(i,j+1)+HH(i-1,j)+HH(i,j))
              H_2=0.25*(HH(i-1,j-1)+HH(i,j-1)+HH(i-1,j)+HH(i,j))
              if(H_1.lt.hi)H_1=hi
              if(H_2.lt.hi)H_2=hi
              UVH_eta=1./(deta)*(
     *                0.5*(UU(i,j+1)+UU(i,j))
     *               *0.5*(VV(i-1,j+1)+VV(i,j+1))
     *               /(H_1)
     *               -0.5*(UU(i,j-1)+UU(i,j))
     *               *0.5*(VV(i-1,j)+VV(i,j))
     *               /(H_2)
     *              )

            endif
            if(Np.eq.11.or.Np.eq.16.or.Np.eq.21.or.Np.eq.26)then
              H_1=0.25*(HH(i-1,j+1)+HH(i,j+1)+HH(i-1,j)+HH(i,j))
              if(H_1.lt.hi)H_1=hi
              UVH_eta=1./(deta)*(
     *                0.5*(UU(i,j+1)+UU(i,j))
     *               *0.5*(VV(i-1,j+1)+VV(i,j+1))
     *               /(H_1)
     *               -0.
     *              )
            endif
            if(Np.eq.13.or.Np.eq.17.or.Np.eq.23.or.Np.eq.27)then
              H_2=0.25*(HH(i-1,j-1)+HH(i,j-1)+HH(i-1,j)+HH(i,j))
              if(H_2.lt.hi)H_2=hi
              UVH_eta=1./(deta)*(
     *                0.
     *               -0.5*(UU(i,j-1)+UU(i,j))
     *               *0.5*(VV(i-1,j)+VV(i,j))
     *               /(H_2)
     *              )
            endif
            if(Np.eq.15.or.Np.eq.18.or.Np.eq.14
     *        .or.Np.eq.25.or.Np.eq.28.or.Np.eq.24)then
              UVH_eta=0.
            endif
            
c --- elevation

            if(Np1.eq.25.or.Np1.eq.26.or.Np1.eq.21
     *        .or.Np1.eq.28.or.Np1.eq.27.or.Np1.eq.23)then
              UVH_eta=0.
            endif

c -- Z_eta
            if(Np.eq.0.or.Np.eq.12.or.Np.eq.22)then
              H_1=Z(i-1,j+1)+H(i-1,j+1)
              H_2=Z(i,j+1)+H(i,j+1)
              H_3=Zorig(i-1)+H(i-1,j-1)
              H_4=Zorig(i)+H(i,j-1)
              if(Ntype(i-1,j+1).eq.9)H_1=-1.*delta
              if(Ntype(i,j+1).eq.9)H_2=-1.*delta
              if(Ntype(i-1,j-1).eq.9)H_3=-1.*delta
              if(Ntype(i,j-1).eq.9)H_4=-1.*delta
              Z_ctr=1.
              if(H_1.gt.0..and.H_2.gt.0.)then
                Z_t=0.5*(Z(i-1,j+1)+Z(i,j+1))
              else
                if(H_1.le.0..and.H_2.gt.0.)Z_t=Z(i,j+1)
                if(H_1.gt.0..and.H_2.le.0.)Z_t=Z(i-1,j+1)
                if(H_1.le.0..and.H_2.le.0.)Z_ctr=0.
              endif

              if(H_3.gt.0..and.H_4.gt.0.)then
                Z_b=0.5*(Zorig(i-1)+Zorig(i))
              else
                if(H_3.le.0..and.H_4.gt.0.)Z_b=Zorig(i)
                if(H_3.gt.0..and.H_4.le.0.)Z_b=Zorig(i-1)
                if(H_3.le.0..and.H_4.le.0.)Z_ctr=0.
              endif

              Z_eta=1./(2.*deta)*(Z_t
     *              -Z_b)*Z_ctr
            endif

            if(Np.eq.11.or.Np.eq.16.or.Np.eq.21.or.Np.eq.26)then
              H_1=Z(i-1,j+1)+H(i-1,j+1)
              H_2=Z(i,j+1)+H(i,j+1)
              H_3=Z(i-1,j)+H(i-1,j)
              H_4=Z(i,j)+H(i,j)
              Z_ctr=1.
              if(H_1.gt.0..and.H_2.gt.0.)then
                Z_t=0.5*(Z(i-1,j+1)+Z(i,j+1))
              else
                if(H_1.le.0..and.H_2.gt.0.)Z_t=Z(i,j+1)
                if(H_1.gt.0..and.H_2.le.0.)Z_t=Z(i-1,j+1)
                if(H_1.le.0..and.H_2.le.0.)Z_ctr=0.
              endif

              if(H_3.gt.0..and.H_4.gt.0.)then
                Z_b=0.5*(Z(i-1,j)+Z(i,j))
              else
                if(H_3.le.0..and.H_4.gt.0.)Z_b=Z(i,j)
                if(H_3.gt.0..and.H_4.le.0.)Z_b=Z(i-1,j)
                if(H_3.le.0..and.H_4.le.0.)Z_ctr=0.
              endif
              Z_eta=1./(1.*deta)*(Z_t
     *              -Z_b)*Z_ctr
            endif

            if(Np.eq.17.or.Np.eq.13.or.Np.eq.27.or.Np.eq.23)then
              H_1=Z(i-1,j)+H(i-1,j)
              H_2=Z(i,j)+H(i,j)
              H_3=Zorig(i-1)+H(i-1,j-1)
              H_4=Zorig(i)+H(i,j-1)
              Z_ctr=1.
              if(H_1.gt.0..and.H_2.gt.0.)then
                Z_t=0.5*(Z(i-1,j)+Z(i,j))
              else
                if(H_1.le.0..and.H_2.gt.0.)Z_t=Z(i,j)
                if(H_1.gt.0..and.H_2.le.0.)Z_t=Z(i-1,j)
                if(H_1.le.0..and.H_2.le.0.)Z_ctr=0.
              endif

              if(H_3.gt.0..and.H_4.gt.0.)then
                Z_b=0.5*(Zorig(i-1)+Zorig(i))
              else
                if(H_3.le.0..and.H_4.gt.0.)Z_b=Zorig(i)
                if(H_3.gt.0..and.H_4.le.0.)Z_b=Zorig(i-1)
                if(H_3.le.0..and.H_4.le.0.)Z_ctr=0.
              endif
              Z_eta=1./(1.*deta)*(Z_t
     *              -Z_b)*Z_ctr
            endif

            if(Np.eq.14.or.Np.eq.24)then
              H_1=Z(i,j+1)+H(i,j+1)
              H_3=Zorig(i)+H(i,j-1)
              Z_ctr=1.
              if(H_1.gt.0.)then
                Z_t=Z(i,j+1)
              else
                Z_ctr=0.
              endif

              if(H_3.gt.0.)then
                Z_b=Zorig(i)
              else
                Z_ctr=0.
              endif
              Z_eta=1./(2.*deta)*(Z_t
     *              -Z_b)*Z_ctr
            endif

            if(Np.eq.15.or.Np.eq.25)then
              H_1=Z(i,j+1)+H(i,j+1)
              H_3=Z(i,j)+H(i,j)
              Z_ctr=1.
              if(H_1.gt.0.)then
                Z_t=Z(i,j+1)
              else
                Z_ctr=0.
              endif

              if(H_3.gt.0.)then
                Z_b=Z(i,j)
              else
                Z_ctr=0.
              endif
              Z_eta=1./(1.*deta)*(Z_t
     *              -Z_b)*Z_ctr
            endif

            if(Np.eq.18.or.Np.eq.28)then
              H_1=Z(i,j)+H(i,j)
              H_3=Zorig(i)+H(i,j-1)
              Z_ctr=1.
              if(H_1.gt.0.)then
                Z_t=Z(i,j)
              else
                Z_ctr=0.
              endif

              if(H_3.gt.0.)then
                Z_b=Zorig(i)
              else
                Z_ctr=0.
              endif
              Z_eta=1./(1.*deta)*(Z_t
     *              -Z_b)*Z_ctr
            endif
c -- get Exxx
            Eadv1=1./u_rj*(UUH_xi+UVH_eta)
            Eadv2=1./(u_rj**2*u_H)*(adu1(i,j)*UU(i,j)+adu2(i,j)*u_v)
            Eadv3=1./(u_rj**2*u_H)*(adu3(i,j))*u_v**2
            Ecor1=-1./u_rj*f*u_beta
            Ecor2=-1./u_rj*f*u_alpha*u_v
            Egrd1=g*u_H/u_rj*u_alpha
            Egrd2=g*u_H/u_rj*u_beta*Z_eta
            Ewind=1./rho*(Yeta(i,j)*Tao_ax-Xeta(i,j)*Tao_ay)
            if(f_cwc.eq.0)then
            f_cwc=9.8*rmanning**2/((u_H/1.)**0.3333)
            endif
            Cd=f_cwc
            Efrc=Cd/(u_rj*u_H**2)*
     *           SQRT(u_gamma*(UU(i,j)+ff_11(i,j)*u_rj)**2+
     *           u_alpha*(u_v+ff_12(i,j)*u_rj)**2
     *       +2.*u_beta*(UU(i,j)+ff_11(i,j)*u_rj)
     *          *(u_v+ff_12(i,j)*u_rj))
            Efrc_corr=-Efrc*(ff_11(i,j)*u_rj)
            Tau_1(i,j)=Efrc
	    Ediff=diffx(i,j)
            Edisp=dispx(i,j)
            Ewave=wavefx(i,j)

c -- balance analysis
            ADVX(i,j)=Eadv1+Eadv2*UU(i,j)+Eadv3
c            ADEY(i,j)=Z_eta*g*u_H/u_rj*u_gamma
            if(Np.eq.0)then
            PREX(i,j)=(z(i,j)-z(i-1,j))/dxi*g*u_H/u_rj*u_alpha
            else
            PREX(i,j)=0.
            endif
            RADX(i,j)=Ewave
            FRCX(i,j)=Efrc*UU(i,j)

c --- specified boundary condition

            if(mod(int(Np1/10),10).ge.2)Ewave=0.


c -- get Ea Eb
c --- wet and dry point

          if(Np.ne.14.and.Np.ne.15.and.Np.ne.18.and.
     *      Np.ne.19.and.Np.lt.80)then
            H_L=HH(i-1,j)
            H_R=HH(i,j)
            Z_=Z(i,j)-Z(i-1,j)

            if((H_R.gt.0..and.H_L.gt.0.).or.
     *         (H_R.gt.0..and.H_L.le.0..and.Z_.gt.delta).or.
     *         (H_R.le.0..and.H_L.gt.0..and.Z_.lt.(-1.*delta)))
     *      then    
              Ea(i)=Egrd1*dt/1./(1.+dt/1.*(Eadv2+Ecor1+Efrc))
              Eb(i)=(UU(i,j)+dt/1.
     *              *(Egrd2+Ewind-Eadv1-Eadv3-Ecor2
     *                +Ediff+Edisp+Ewave+Efrc_corr))
     *                     /(1.+dt/1.*(Eadv2+Ecor1+Efrc))
            else
              Ea(i)=0.
              Eb(i)=0.
            endif
          endif

c --- wall 
            if(Np1.eq.14.or.Np1.eq.15.or.Np1.eq.18)then
              Ea(i)=0.
              Eb(i)=0.
            endif

            if(Np1.eq.12.or.Np1.eq.16.or.Np1.eq.17)then
              Ea(i+1)=0.
              Eb(i+1)=0.
            endif

c --- elevation corners

            if(Np1.eq.25.or.Np1.eq.28)then
              Ea(i)=0.
              Eb(i)=0.
            endif

            if(Np1.eq.26.or.Np1.eq.27)then
              Ea(i+1)=0.
              Eb(i+1)=0.
            endif

c --- flux
c --- east and west

            if(Np1.eq.34.or.Np1.eq.35.or.Np1.eq.38)then
              Ea(i)=0.
              Eb(i)=Uopl(j)
            endif

            if(Np1.eq.32.or.Np1.eq.36.or.Np1.eq.37)then
              Ea(i+1)=0.
              Eb(i+1)=Uopr(j)
            endif

c --- make sure periodic condition
            if(Np1.eq.64.or.Np1.eq.65.or.Np1.eq.68)then
              Ea(i)=Ea(iend(j)+1)
              Eb(i)=Eb(iend(j)+1)
            endif  

            if(Np1.eq.62.or.Np1.eq.66.or.Np1.eq.67)then
              Ea(i+1)=Ea(ibeg(j))
              Eb(i+1)=Eb(ibeg(j))
            endif            

200       continue

 
c -- Ea Eb over
            im=iend(j)-ibeg(j)+1
c -- A B C
            do i=1,im
              Np=Ntype(i+ibeg(j)-1,j)
           if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
              i_tr=i+ibeg(j)-1
              if(Np.eq.0.or.Np.eq.12.or.Np.eq.14)then
                z_rj=rJ(i_tr,j)
                Res=1./z_rj*(Ea(i_tr)+Ea(i_tr+1))
     *              /(dxi**2)+2./dt
                A(i)=-Ea(i_tr)/(z_rj*dxi**2)/Res
                C(i)=-Ea(i_tr+1)/z_rj/(dxi**2)/Res

                if(Ea(i_tr).eq.0..and.Ea(i_tr+1).eq.0.)then
                  D(i)=Z(i_tr,j)
                else
                  D(i)=1./Res*(-1./z_rj*(Eb(i_tr+1)
     *               -Eb(i_tr))/dxi
     *               -1./z_rj*(VV(i_tr,j+1)-VV(i_tr,j))/
     *               deta+2./dt*Z(i_tr,j)
     *                )
                endif

              endif

              if(Np.eq.11.or.Np.eq.15.or.Np.eq.16)then
                z_rj=rJ(i_tr,j)
                Res=1./z_rj*(Ea(i_tr)+Ea(i_tr+1))/(dxi**2)+2./dt
                A(i)=-Ea(i_tr)/(z_rj*dxi**2)/Res
                C(i)=-Ea(i_tr+1)/z_rj/(dxi**2)/Res

                if(Ea(i_tr).eq.0..and.Ea(i_tr+1).eq.0.)then
                  D(i)=Z(i_tr,j)
                else
                  D(i)=1./Res*(-1./z_rj*(Eb(i_tr+1)-Eb(i_tr))/dxi
     *               -1./z_rj*(VV(i_tr,j+1)-0.)/
     *               deta+2./dt*Z(i_tr,j)
     *                )
                endif

              endif

              if(Np.eq.13.or.Np.eq.17.or.Np.eq.18)then
                z_rj=rJ(i_tr,j)
                Res=1./z_rj*(Ea(i_tr)+Ea(i_tr+1))/(dxi**2)+2./dt
                A(i)=-Ea(i_tr)/(z_rj*dxi**2)/Res
                C(i)=-Ea(i_tr+1)/z_rj/(dxi**2)/Res

                if(Ea(i_tr).eq.0..and.Ea(i_tr+1).eq.0.)then
                  D(i)=Z(i_tr,j)
                else
                  D(i)=1./Res*(-1./z_rj*(Eb(i_tr+1)-Eb(i_tr))/dxi
     *               -1./z_rj*(0.-VV(i_tr,j))/
     *               deta+2./dt*Z(i_tr,j)
     *                )
                endif
              endif

              if(Np.ge.80)then
                A(i)=0.
                C(i)=0.
                D(i)=Z(i_tr,j)
              endif

            enddo

c -- A B C ok
c -- remain Z
          do i=ibeg(j),iend(j)
            Zorig(i)=Z(i,j)
          enddo

c -- solve trig
            ileft=Ntype(ibeg(j),j)
            iright=Ntype(iend(j),j)

c -- define btype: 0=wall, flux
c ---              1=gradient, elevation
c                  2=periodic

            if(mod(int(ileft/10),10).eq.1.or.
     &         mod(int(ileft/10),10).eq.3)btype_left=0
            if(mod(int(iright/10),10).eq.1.or.
     &         mod(int(iright/10),10).eq.3)btype_right=0
            if(mod(int(ileft/10),10).eq.2.or.
     &         mod(int(ileft/10),10).eq.4)btype_left=1
            if(mod(int(iright/10),10).eq.2.or.
     &         mod(int(iright/10),10).eq.4)btype_right=1
            if(east_ele2_flx3_grd4_rad5_per6.eq.6)btype_left=2

c -- periodic

            if(btype_left.eq.2)then

             write(*,*)'some problems with a b c'
             stop

             if(a(ibeg(j)).ne.0)then
             mtrig=im
              call trig_periodic(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j)-1,j)=Ztrig(i)
              enddo
              else
              mtrig=im
              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j)-1,j)=Ztrig(i)
              enddo
              endif
            endif

c -- two wall, flux
            if(btype_left.eq.0.and.
     &         btype_right.eq.0)then 
              mtrig=im
              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j)-1,j)=Ztrig(i)
              enddo
            endif
c -- left elevation boundary
            if((ileft.eq.24.or.ileft.eq.25.or.ileft.eq.28).and.
     *         btype_right.eq.0)then
              mtrig=im-1
              C(1)=C(2)
              D(1)=D(2)-A(2)*Zopl(j)
              do i=2,mtrig
                A(i)=A(i+1)
                C(i)=C(i+1)
                D(i)=D(i+1)
              enddo

              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j),j)=Ztrig(i)
              enddo
              Z(ibeg(j),j)=Zopl(j)

            endif
c -- right elevation boundary
            if((iright.eq.22.or.iright.eq.26.or.iright.eq.27).and.
     *         btype_left.eq.0)then
              mtrig=im-1
              D(mtrig)=D(mtrig)-C(mtrig)*Zopr(j)

              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j)-1,j)=Ztrig(i)
              enddo
              Z(iend(j),j)=Zopr(j)

            endif
c -- left and right elevation boundary
            if(mod(int(ileft/10),10).eq.2.and.
     *        mod(int(iright/10),10).eq.2)then
              mtrig=im-2
              C(1)=C(2)
              D(1)=D(2)-A(2)*Zopl(j)
              do i=2,mtrig-1
                A(i)=A(i+1)
                C(i)=C(i+1)
                D(i)=D(i+1)
              enddo
              A(mtrig)=A(mtrig+1)
              D(mtrig)=D(mtrig+1)-C(mtrig+1)*Zopr(j)

              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j),j)=Ztrig(i)
              enddo
              Z(ibeg(j),j)=Zopl(j)
              Z(iend(j),j)=Zopr(j)
            endif

c --  bottom and top elevation boundary
            do i=ibeg(j),iend(j)
              Np=Ntype(i,j)
              if(Np.eq.21.or.Np.eq.25.or.Np.eq.26) Z(i,j)=Zopb(i)
              if(Np.eq.23.or.Np.eq.27.or.Np.eq.28) Z(i,j)=Zopt(i)
            enddo

c --  bottom and top gradient boundary
            do i=ibeg(j),iend(j)
              Np=Ntype(i,j)
              if(Np.eq.41.or.Np.eq.45.or.Np.eq.46) Z(i,j)=Z(i,j+1)
              if(Np.eq.43.or.Np.eq.47.or.Np.eq.48) Z(i,j)=Z(i,j-1)
            enddo

c --- left gradient boundary
            if((ileft.eq.44.or.ileft.eq.45.or.ileft.eq.48).and.
     *         btype_right.eq.0)then
              mtrig=im-1
              C(1)=C(2)/(1+A(2))
              D(1)=D(2)/(1+A(2))
              do i=2,mtrig
                A(i)=A(i+1)
                C(i)=C(i+1)
                D(i)=D(i+1)
              enddo

              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j),j)=Ztrig(i)
              enddo
              Z(ibeg(j),j)=Z(ibeg(j)+1,j)
            endif

c -- right gradient boundary
            if((iright.eq.42.or.iright.eq.46.or.iright.eq.47).and.
     *         btype_left.eq.0)then
              mtrig=im-1
              D(mtrig)=D(mtrig)/(1.+C(mtrig))
              A(mtrig)=A(mtrig)/(1.+C(mtrig))

              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j)-1,j)=Ztrig(i)
              enddo
              Z(iend(j),j)=Z(iend(j)-1,j)

            endif


c -- left and right gradient boundary
            if(mod(int(ileft/10),10).eq.4.and.
     *        mod(int(iright/10),10).eq.4)then
              mtrig=im-2
              C(1)=C(2)/(1+A(2))
              D(1)=D(2)/(1+A(2))
              do i=2,mtrig-1
                A(i)=A(i+1)
                C(i)=C(i+1)
                D(i)=D(i+1)
              enddo

              A(mtrig)=A(mtrig+1)/(1.+C(mtrig+1))
              D(mtrig)=D(mtrig+1)/(1.+C(mtrig)+1)        

              call trig(A,C,D,Ztrig,mtrig)
              do i=1,mtrig
                Z(i+ibeg(j),j)=Ztrig(i)
              enddo
              Z(ibeg(j),j)=Z(ibeg(j)+1,j)
              Z(iend(j),j)=Z(iend(j)-1,j)
            endif

c -- z computed over
c -- computed UU
         do i=ibeg(j),iend(j)
           Np=Ntype(i,j)
           Np1=Ntype(i,j)
           if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
           if(Np.eq.0.or.Np.eq.11.or.Np.eq.16.or.Np.eq.12.or.Np.eq.17
     *     .or.Np.eq.13)then
             UU(i,j)=-Ea(i)/dxi*(Z(i,j)-Z(i-1,j))+Eb(i)
           endif
           if(Np.eq.15.or.Np.eq.18.or.Np.eq.14)then
             UU(i,j)=0.
           endif
           if(Np1.eq.35.or.Np1.eq.38.or.Np1.eq.34)then
              UU(i,j)=Uopl(j)
           endif
           if(Np1.eq.36.or.Np1.eq.37.or.Np1.eq.32)then
              UU(i+1,j)=Uopr(j)
           endif
           if(Np1.eq.65.or.Np1.eq.68.or.Np1.eq.64)then
             UU(i,j)=-Ea(i)/dxi*(Z(i,j)-Z(iend(j),j))+Eb(i)             
           endif
           if(Np1.eq.66.or.Np1.eq.67.or.Np1.eq.62)then
             UU(i+1,j)=-Ea(i+1)/dxi*(Z(ibeg(j),j)-Z(i,j))+Eb(i+1)             
           endif
c           if(Np1.eq.25.or.Np1.eq.26.or.Np1.eq.21)then
c             UU(i,j)=0.
c           endif
c           if(Np1.eq.28.or.Np1.eq.27.or.Np1.eq.23)then
c             UU(i,j)=0.
c           endif

         enddo
c -- OK

100     continue

c --- gradient boundary condition
        do j=1,n-1
        do i=ibeg(j),iend(j)
           Np1=Ntype(i,j)
           if(Np1.eq.45.or.Np1.eq.41.or.Np1.eq.46)then
             UU(i,j)=UU(i,j+1)
           endif
           if(Np1.eq.48.or.Np1.eq.43.or.Np1.eq.47)then
             UU(i,j)=UU(i,j-1)
           endif
        enddo
        enddo

c ---  flux boundary condition
c        do j=1,n-1
c        do i=ibeg(j),iend(j)
c           Np1=Ntype(i,j)
c           if(Np1.eq.35.or.Np.eq.36.or.Np.eq.31)then
c             UU(i,j)=0.
c           endif
c           if(Np1.eq.38.or.Np.eq.37.or.Np.eq.33)then
c             UU(i,j)=0.
c           endif
c        enddo
c        enddo

c -- set dry point and HH

	do j=1,n-1
        do i=ibeg(j),iend(j)
          HH(i,j)=H(i,j)+Z(i,j)	  
	  if(HH(i,j).lt.0.)then
            Z(i,j)=-H(i,j)-delta
	    HH(i,j)=H(i,j)+Z(i,j)
          endif
        enddo
        enddo	

c -- ok

        return
        end

c ---------------------------------
c       trig
c ---------------------------------
        subroutine trig(A,C,D,Z,m)
        real a(m),c(m),d(m),z(m)
        integer i,j

        do i=2,m
          if(a(i).ne.0.)then
            c(i)=c(i)/a(i)/(1./a(i)-c(i-1))
            d(i)=(d(i)/a(i)-d(i-1))/(1./a(i)-c(i-1))
          endif
        enddo

        z(m)=d(m)

        do i=m-1,1,-1
          z(i)=d(i)-c(i)*z(i+1)
        enddo

        return
        end

c ---------------------------------
c       trig periodic
c ---------------------------------
        subroutine trig_periodic(A,C,D,Z,m)
        real a(m),c(m),d(m),z(m)
        real cp(m),dp(m),ep(m),ff,gg
        integer i,j

        if(a(1).ne.0.and.c(m).ne.0)then

        cp(1)=c(1)
        dp(1)=d(1)
        ep(1)=-a(1)
        do i=2,m-1
          cp(i)=c(i)/(1.-a(i)*cp(i-1))
          dp(i)=(d(i)-a(i)*dp(i-1))/(1.-a(i)*cp(i-1))
          ep(i)=-a(i)*ep(i-1)/(1.-a(i)*cp(i-1))
        enddo

        ff=-c(m)/(1-a(m)*cp(m-1)+a(m)*ep(m-1))
        gg=(d(m)-a(m)*dp(m-1))/(1.-a(m)*cp(m-1)+a(m)*ep(m-1))
        
        c(1)=c(1)/(1.+a(1)*ff)
        d(1)=(d(1)-a(1)*gg)/(1.+a(1)*ff)
        a(m)=a(m)/(1.+c(m)/ff)
        d(m)=(d(m-1)+c(m)*gg/ff)/(1.+c(m)/ff)

        endif

        do i=2,m
          if(a(i).ne.0.)then
            c(i)=c(i)/a(i)/(1./a(i)-c(i-1))
            d(i)=(d(i)/a(i)-d(i-1))/(1./a(i)-c(i-1))
          endif
        enddo

        z(m)=d(m)

        do i=m-1,1,-1
          z(i)=d(i)-c(i)*z(i+1)
        enddo


        return
        end

c ---------------------------------
c       zeta_v
c ---------------------------------
        subroutine solve_zeta_v
!        use sadi
        use sadi
        integer Np,Np1
        real Y_xi,X_xi
        real A(mtrigm),C(mtrigm),D(mtrigm)
        real Ea(mtrigm),Eb(mtrigm),Ztrig(mtrigm)
        real Zorig(ntrim)
	
        do 100 i=1,m-1
c -- begin row
          if(jbeg(i).eq.jend(i)) goto 100
          do 200 j=jbeg(i),jend(i)
            Np=Ntype(i,j)
            Np1=Ntype(i,j)
            if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
c -- v_rJ et al

            if(Np.eq.0.or.Np.eq.14.or.Np.eq.24.
     *        or.np.eq.13.or.Np.eq.23.or.Np.eq.18.or.Np.eq.28)then
              v_rj=0.5*(rJ(i,j-1)+rJ(i,j))
              v_H=0.5*(h(i,j-1)+h(i,j)+z(i,j-1)+z(i,j))
              if(v_H.lt.hi) v_H=hi
              v_u=0.25*(uu(i,j)+uu(i,j-1)+uu(i+1,j-1)+uu(i+1,j))
              v_beta=0.5*(beta(i,j-1)+beta(i,j))
              v_alpha=0.5*(alpha(i,j-1)+alpha(i,j))
              v_gamma=0.5*(gamma(i,j-1)+gamma(i,j))
              Y_xi=1./dxi*(Y(i+1,j)-Y(i,j))
              X_xi=1./dxi*(X(i+1,j)-X(i,j))
            endif

            if(Np.eq.12.or.Np.eq.22.or.Np.eq.17.or.Np.eq.27)then
              v_rj=0.5*(rJ(i,j-1)+rJ(i,j))
              v_H=0.5*(h(i,j-1)+h(i,j)+z(i,j-1)+z(i,j))
              if(v_H.lt.hi) v_H=hi
              v_u=0.25*(uu(i,j)+0.+0.+uu(i,j-1))
              v_beta=0.5*(beta(i,j-1)+beta(i,j))
              v_alpha=0.5*(alpha(i,j-1)+alpha(i,j))
              v_gamma=0.5*(gamma(i,j-1)+gamma(i,j))
              Y_xi=1./dxi*(Y(i+1,j)-Y(i,j))
              X_xi=1./dxi*(X(i+1,j)-X(i,j))
            endif

            if(Np.eq.11.or.Np.eq.21.or.Np.eq.15.or.np.eq.25)then
              v_rj=rJ(i,j)
              v_H=h(i,j)+z(i,j)
              if(v_H.lt.hi) v_H=hi
              v_v=0.5*(uu(i,j)+uu(i+1,j))
              v_beta=beta(i,j)
              v_alpha=alpha(i,j)
              v_gamma=gamma(i,j)
              Y_xi=1./dxi*(Y(i+1,j)-Y(i,j))
              X_xi=1./dxi*(X(i+1,j)-X(i,j))
            endif

            if(Np.eq.16.or.Np.eq.26)then
              v_rj=rJ(i,j)
              v_H=h(i,j)+z(i,j)
              if(v_H.lt.hi) v_H=hi
              v_v=0.5*(vv(i,j)+0.)
              v_beta=beta(i,j)
              v_alpha=alpha(i,j)
              v_gamma=gamma(i,j)
              Y_xi=1./dxi*(Y(i+1,j)-Y(i,j))
              X_xi=1./dxi*(X(i+1,j)-X(i,j))
            endif
c -- u_rj over

c --  \partial v**2/H
            if(Np.eq.0.or.Np.eq.14.or.Np.eq.24.or.Np.eq.12.or.Np.eq.22)
     *        then
              H_1=HH(i,j)
              H_2=HH(i,j-1)
              if(H_1.lt.hi) H_1=hi
              if(H_2.lt.hi) H_2=hi
              VVH_eta=1./(deta)*(0.25*(VV(i,j+1)+VV(i,j))**2/H_1
     *             -0.25*(VV(i,j-1)+VV(i,j))**2/H_2)
            endif

            if(Np.eq.18.or.Np.eq.17.or.Np.eq.13
     *        .or.Np.eq.28.or.Np.eq.27.or.Np.eq.23)then
              H_1=HH(i,j)
              H_2=HH(i,j-1)
              if(H_1.lt.hi) H_1=hi
              if(H_2.lt.hi) H_2=hi
              VVH_eta=1./(deta)*(0.25*(0.+VV(i,j))**2/H_1
     *             -0.25*(VV(i,j-1)+VV(i,j))**2/H_2)
            endif

            if(Np.eq.15.or.Np.eq.11.or.Np.eq.16
     *        .or.Np.eq.25.or.Np.eq.21.or.Np.eq.26)then
              H_1=HH(i,j)
              if(H_1.lt.hi) H_1=hi
              VVH_eta=1./(deta/2.)*(0.25*(VV(i,j+1)+0.)**2/H_1
     *             -0.)
            endif

c ---   periodic
            if(Np1.eq.61)
     *        then
              H_1=HH(i,j)
              if(H_1.lt.hi)H_1=hi
              VVH_eta=1./(deta)*(0.25*(VV(i,j+1)+VV(i,j))**2/H_1
     *             -0.25*(VV(i,jend(i))+VV(i,j))**2/H_1)
            endif

            if(Np1.eq.68.or.Np1.eq.63.or.Np1.eq.67)
     *        then
              H_1=HH(i,j)
              if(H_1.lt.hi)H_1=hi
              VVH_eta=1./(deta)*(0.25*(VV(i,jbeg(i))+VV(i,j))**2/H_1
     *             -0.25*(VV(i,j-1)+VV(i,j))**2/H_1)
            endif

c --- elevation

            if(Np1.eq.25.or.Np1.eq.28.or.Np1.eq.24
     *        .or.Np1.eq.26.or.Np1.eq.27.or.Np1.eq.22)then
              UVH_eta=0.
            endif


c -- \partial uv/H
            if(Np.eq.0.or.Np.eq.13.or.Np.eq.23)then
              H_1=0.25*(HH(i+1,j-1)+HH(i+1,j)+HH(i,j-1)+HH(i,j))
              H_2=0.25*(HH(i-1,j-1)+HH(i-1,j)+HH(i,j-1)+HH(i,j))
              if(H_1.lt.hi) H_1=hi
              if(H_2.lt.hi) H_2=hi
              UVH_xi=1./(dxi)*(
     *                0.5*(VV(i+1,j)+VV(i,j))
     *               *0.5*(UU(i+1,j-1)+UU(i+1,j))
     *               /(H_1)
     *               -0.5*(VV(i-1,j)+VV(i,j))
     *               *0.5*(UU(i,j-1)+UU(i,j))
     *               /(H_2)
     *              )
            endif
            if(Np.eq.14.or.Np.eq.18.or.Np.eq.24.or.Np.eq.28)then
              H_1=0.25*(HH(i+1,j-1)+HH(i+1,j)+HH(i,j-1)+HH(i,j))
              if(H_1.lt.hi) H_1=hi
              UVH_xi=1./(dxi)*(
     *                0.5*(VV(i+1,j)+VV(i,j))
     *               *0.5*(UU(i+1,j-1)+UU(i+1,j))
     *               /(H_1)
     *               -0.
     *              )
            endif
            if(Np.eq.12.or.Np.eq.17.or.Np.eq.22.or.Np.eq.27)then
              H_2=0.25*(HH(i-1,j-1)+HH(i-1,j)+HH(i,j-1)+HH(i,j))
              if(H_2.lt.hi) H_2=hi
              UVH_xi=1./(dxi)*(
     *                0.
     *               -0.5*(VV(i-1,j)+VV(i,j))
     *               *0.5*(UU(i,j-1)+UU(i,j))
     *               /(H_2)
     *              )
            endif
            if(Np.eq.15.or.Np.eq.16.or.Np.eq.11
     *        .or.Np.eq.25.or.Np.eq.26.or.Np.eq.21)then
              UVH_xi=0.
            endif

c --- elevation

            if(Np1.eq.25.or.Np1.eq.28.or.Np1.eq.24
     *        .or.Np1.eq.26.or.Np1.eq.27.or.Np1.eq.22)then
              UVH_xi=0.
            endif

c -- Z_xi
            if(Np.eq.0.or.Np.eq.13.or.Np.eq.23)then
              H_1=Z(i+1,j-1)+H(i+1,j-1)
              H_2=Z(i+1,j)+H(i+1,j)
              H_3=Zorig(j-1)+H(i-1,j-1)
              H_4=Zorig(j)+H(i-1,j)
              if(Ntype(i+1,j-1).eq.9)H_1=-1.*delta
              if(Ntype(i+1,j).eq.9)H_2=-1.*delta
              if(Ntype(i-1,j-1).eq.9)H_3=-1.*delta
              if(Ntype(i-1,j).eq.9)H_4=-1.*delta
              Z_ctr=1.
              if(H_1.gt.0..and.H_2.gt.0.)then
                Z_r=0.5*(Z(i+1,j-1)+Z(i+1,j))
              else
                if(H_1.le.0..and.H_2.gt.0.)Z_r=Z(i+1,j)
                if(H_1.gt.0..and.H_2.le.0.)Z_r=Z(i+1,j-1)
                if(H_1.le.0..and.H_2.le.0.)Z_ctr=0.
              endif

              if(H_3.gt.0..and.H_4.gt.0.)then
                Z_l=0.5*(Zorig(j-1)+Zorig(j))
              else
                if(H_3.le.0..and.H_4.gt.0.)Z_l=Zorig(j)
                if(H_3.gt.0..and.H_4.le.0.)Z_l=Zorig(j-1)
                if(H_3.le.0..and.H_4.le.0.)Z_ctr=0.
              endif
              Z_xi=1./(2.*dxi)*(Z_r
     *              -Z_l)*Z_ctr
            endif

            if(Np.eq.14.or.Np.eq.18.or.Np.eq.24.or.Np.eq.28)then
              H_1=Z(i+1,j-1)+H(i+1,j-1)
              H_2=Z(i+1,j)+H(i+1,j)
              H_3=Z(i,j-1)+H(i,j-1)
              H_4=Z(i,j)+H(i,j)
              Z_ctr=1.
              if(H_1.gt.0..and.H_2.gt.0.)then
                Z_r=0.5*(Z(i+1,j-1)+Z(i+1,j))
              else
                if(H_1.le.0..and.H_2.gt.0.)Z_r=Z(i+1,j)
                if(H_1.gt.0..and.H_2.le.0.)Z_r=Z(i+1,j-1)
                if(H_1.le.0..and.H_2.le.0.)Z_ctr=0.
              endif

              if(H_3.gt.0..and.H_4.gt.0.)then
                Z_l=0.5*(Z(i,j-1)+Z(i,j))
              else
                if(H_3.le.0..and.H_4.gt.0.)Z_l=Z(i,j)
                if(H_3.gt.0..and.H_4.le.0.)Z_l=Z(i,j-1)
                if(H_3.le.0..and.H_4.le.0.)Z_ctr=0.
              endif
              Z_xi=1./(1.*dxi)*(Z_r
     *              -Z_l)*Z_ctr
            endif

            if(Np.eq.17.or.Np.eq.12.or.Np.eq.27.or.Np.eq.22)then
              H_1=Z(i,j-1)+H(i,j-1)
              H_2=Z(i,j)+H(i,j)
              H_3=Zorig(j-1)+H(i-1,j-1)
              H_4=Zorig(j)+H(i-1,j)
              Z_ctr=1.
              if(H_1.gt.0..and.H_2.gt.0.)then
                Z_r=0.5*(Z(i,j-1)+Z(i,j))
              else
                if(H_1.le.0..and.H_2.gt.0.)Z_r=Z(i,j)
                if(H_1.gt.0..and.H_2.le.0.)Z_r=Z(i,j-1)
                if(H_1.le.0..and.H_2.le.0.)Z_ctr=0.
              endif

              if(H_3.gt.0..and.H_4.gt.0.)then
                Z_l=0.5*(Zorig(j-1)+Zorig(j))
              else
                if(H_3.le.0..and.H_4.gt.0.)Z_l=Zorig(j)
                if(H_3.gt.0..and.H_4.le.0.)Z_l=Zorig(j-1)
                if(H_3.le.0..and.H_4.le.0.)Z_ctr=0.
              endif
              Z_xi=1./(1.*dxi)*(Z_r
     *              -Z_l)*Z_ctr
            endif

            if(Np.eq.11.or.Np.eq.21)then
              H_1=Z(i+1,j)+H(i+1,j)
              H_3=Zorig(j)+H(i-1,j)
              Z_ctr=1.
              if(H_1.gt.0.)then
                Z_r=Z(i+1,j)
              else
                Z_ctr=0.
              endif

              if(H_3.gt.0.)then
                Z_l=Zorig(j)
              else
                Z_ctr=0.
              endif
              Z_xi=1./(2.*dxi)*(Z_r
     *              -Z_l)*Z_ctr
            endif

            if(Np.eq.15.or.Np.eq.25)then
              H_1=Z(i+1,j)+H(i+1,j)
              H_3=Z(i,j)+H(i,j)
              Z_ctr=1.
              if(H_1.gt.0.)then
                Z_r=Z(i+1,j)
              else
                Z_ctr=0.
              endif

              if(H_3.gt.0.)then
                Z_l=Z(i,j)
              else
                Z_ctr=0.
              endif
              Z_xi=1./(1.*dxi)*(Z_r
     *              -Z_l)*Z_ctr
            endif

            if(Np.eq.16.or.Np.eq.26)then
              H_1=Z(i,j)+H(i,j)
              H_3=Zorig(j)+H(i-1,j)
              Z_ctr=1.
              if(H_1.gt.0.)then
                Z_r=Z(i,j)
              else
                Z_ctr=0.
              endif

              if(H_3.gt.0.)then
                Z_l=Zorig(j)
              else
                Z_ctr=0.
              endif
              Z_xi=1./(1.*dxi)*(Z_r
     *              -Z_l)*Z_ctr
            endif
c -- get Exxx
            Eadv1=1./v_rj*(UVH_xi+VVH_eta)
            Eadv2=1./(v_rj**2*v_H)*(adv1(i,j)*VV(i,j)+adv2(i,j)*v_U)
            Eadv3=1./(v_rj**2*v_H)*(adv3(i,j))*v_U**2
            Ecor1=1./v_rj*f*v_beta
            Ecor2=1./v_rj*f*v_gamma*v_U
            Egrd1=g*v_H/v_rj*v_gamma
            Egrd2=g*v_H/v_rj*v_beta*Z_xi
            Ewind=1./rho*(-Yxi(i,j)*Tao_ax+Xxi(i,j)*Tao_ay)
            if(f_cwc.eq.0.)then
            f_cwc=9.8*rmanning**2/((v_H/1.)**0.3333)
            endif
            Cd=f_cwc
            Efrc=Cd/(v_rj*v_H**2)*
     *           SQRT(v_gamma*(v_U+ff_11(i,j)*v_rj)**2+
     *           v_alpha*(VV(i,j)+ff_12(i,j)*v_rj)**2+
     *       2.*v_beta*(VV(i,j)+ff_12(i,j)*v_rj)*(v_U+ff_11(i,j)*v_rj))
            Efrc_corr=-Efrc*ff_12(i,j)*v_rj
            Tau_2(i,j)=Efrc
            Ediff=diffy(i,j)
            Edisp=dispy(i,j)
            Ewave=wavefy(i,j)
c -- balance analysis
            ADVY(i,j)=Eadv1+Eadv2*VV(i,j)+Eadv3
c            PREX(i,j)=Z_xi*g*v_H/v_rj*v_alpha
            if (Np.eq.0)then
            PREY(i,j)=(z(i,j)-z(i,j-1))/deta*g*v_H/v_rj*v_gamma
            else
            PREY(i,j)=0.
            endif
            RADY(i,j)=Ewave
            FRCY(i,j)=Efrc*VV(i,j)

c --- specified boundary condition

            if(mod(int(Np1/10),10).ge.2)Ewave=0.

c -- get Ea Eb  and wet and dry point

          if(Np.ne.11.and.Np.ne.15.and.Np.ne.16.and.Np.ne.19.
     *     and.Np.lt.80)then
            H_T=HH(i,j)
            H_B=HH(i,j-1)
            Z_=Z(i,j)-Z(i,j-1)
            if((H_T.gt.0..and.H_B.gt.0.).or.
     *         (H_T.gt.0..and.H_B.le.0..and.Z_.gt.delta).or.
     *         (H_T.le.0..and.H_B.gt.0..and.Z_.lt.(-1.*delta)))
     *      then
              Ea(j)=Egrd1*dt/1./(1.+dt/1.*(Eadv2+Ecor1+Efrc))
              Eb(j)=(VV(i,j)+dt/1.
     *              *(Egrd2+Ewind-Eadv1-Eadv3-Ecor2
     *                +Ediff+Edisp+Ewave+Efrc_corr))
     *                     /(1.+dt/1.*(Eadv2+Ecor1+Efrc))
            else
              Ea(j)=0.
              Eb(j)=0.
            endif
          endif
c --- wall
            if(Np1.eq.11.or.Np1.eq.15.or.Np1.eq.16)then
              Ea(j)=0.
              Eb(j)=0.
            endif

            if(Np1.eq.13.or.Np1.eq.18.or.Np1.eq.17)then
              Ea(j+1)=0.
              Eb(j+1)=0.
            endif

c --- elevation corners

            if(Np1.eq.25.or.Np1.eq.26)then
              Ea(j)=0.
              Eb(j)=0.
            endif

            if(Np1.eq.28.or.Np1.eq.27)then
              Ea(j+1)=0.
              Eb(j+1)=0.
            endif

c --- flux
            if(Np1.eq.31.or.Np1.eq.35.or.Np1.eq.36)then
              Ea(j)=0.
              Eb(j)=Vopb(i)
            endif

            if(Np1.eq.33.or.Np1.eq.38.or.Np1.eq.37)then
              Ea(j+1)=0.
              Eb(j+1)=Vopt(i)
            endif

c --- make sure periodic condition
            if(Np1.eq.61.or.Np1.eq.65.or.Np1.eq.66)then
              Ea(j)=Ea(jend(i)+1)
              Eb(j)=Eb(jend(i)+1)
            endif  

            if(Np1.eq.63.or.Np1.eq.68.or.Np1.eq.67)then
              Ea(j+1)=Ea(jbeg(i))
              Eb(j+1)=Eb(jbeg(i))
            endif      

200       continue



c -- Ea Eb over
            jm=jend(i)-jbeg(i)+1
c -- A B C
            do j=1,jm
              Np=Ntype(i,j+jbeg(i)-1)
              if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
              j_tr=j+jbeg(i)-1
              if(Np.eq.0.or.Np.eq.11.or.Np.eq.13)then
                z_rj=rJ(i,j_tr)
                Res=1./z_rj*(Ea(j_tr)+Ea(j_tr+1))
     *              /(deta**2)+2./dt
                A(j)=-Ea(j_tr)/(z_rj*deta**2)/Res
                C(j)=-Ea(j_tr+1)/z_rj/(deta**2)/Res
                if(Ea(j_tr).eq.0..and.Ea(j_tr+1).eq.0.)then
                  D(j)=Z(i,j_tr)
                else
                  D(j)=1./Res*(-1./z_rj*(Eb(j_tr+1)
     *               -Eb(j_tr))/deta
     *               -1./z_rj*(UU(i+1,j_tr)-UU(i,j_tr))/
     *               dxi+2./dt*Z(i,j_tr)
     *                )
                endif
              endif


              if(Np.eq.14.or.Np.eq.15.or.Np.eq.18)then
                z_rj=rJ(i,j_tr)
                Res=1./z_rj*(Ea(j_tr)+Ea(j_tr+1))
     *              /(deta**2)+2./dt
                A(j)=-Ea(j_tr)/(z_rj*deta**2)/Res
                C(j)=-Ea(j_tr+1)/z_rj/(deta**2)/Res

                if(Ea(j_tr).eq.0..and.Ea(j_tr+1).eq.0.)then
                  D(j)=Z(i,j_tr)
                else
                  D(j)=1./Res*(-1./z_rj*(Eb(j_tr+1)
     *               -Eb(j_tr))/deta
     *               -1./z_rj*(UU(i+1,j_tr)-0.)/
     *               dxi+2./dt*Z(i,j_tr)
     *                )
                endif
              endif

              if(Np.eq.12.or.Np.eq.17.or.Np.eq.16)then
                z_rj=rJ(i,j_tr)
                Res=1./z_rj*(Ea(j_tr)+Ea(j_tr+1))
     *              /(deta**2)+2./dt
                A(j)=-Ea(j_tr)/(z_rj*deta**2)/Res
                C(j)=-Ea(j_tr+1)/z_rj/(deta**2)/Res

                if(Ea(j_tr).eq.0..and.Ea(j_tr+1).eq.0.)then
                  D(j)=Z(i,j_tr)
                else
                  D(j)=1./Res*(-1./z_rj*(Eb(j_tr+1)
     *               -Eb(j_tr))/deta
     *               -1./z_rj*(0.-UU(i,j_tr))/
     *               dxi+2./dt*Z(i,j_tr)
     *                )
                endif
              endif

              if(Np.ge.80)then
                A(j)=0.
                C(j)=0.
                D(j)=Z(i,j_tr)
              endif

            enddo
c -- A B C ok
c -- remain Z
          do j=jbeg(i),jend(i)
            Zorig(j)=Z(i,j)
          enddo
c -- solve trig

            jbot=Ntype(i,jbeg(i))
            jtop=Ntype(i,jend(i))

c -- define btype: 0=wall, flux
c ---              1=gradient, elevation

            if(mod(int(jbot/10),10).eq.1.or.
     &         mod(int(jbot/10),10).eq.3)btype_bot=0
            if(mod(int(jtop/10),10).eq.1.or.
     &         mod(int(jtop/10),10).eq.3)btype_top=0
            if(mod(int(jbot/10),10).eq.2.or.
     &         mod(int(jbot/10),10).eq.4)btype_bot=1
            if(mod(int(jtop/10),10).eq.2.or.
     &         mod(int(jtop/10),10).eq.4)btype_top=1
            if(south_ele2_flx3_grd4_rad5_per6.eq.6)btype_bot=2

c -- periodic

            if(btype_bot.eq.2)then

               write(*,*)'some problems with a b c'
               stop

             if(a(jbeg(i)).ne.0)then
             mtrig=im
              call trig_periodic(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i)-1)=Ztrig(j)
              enddo
              else
              mtrig=im
              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i)-1)=Ztrig(j)
              enddo
              endif
            endif

c -- no surface elevation open boundary
            if(btype_bot.eq.0.and.
     &         btype_top.eq.0)then
              mtrig=jm
              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i)-1)=Ztrig(j)
              enddo
            endif
c -- bottom elevation boundary
            if((jbot.eq.21.or.jbot.eq.25.or.jbot.eq.26).and.
     *         btype_top.eq.0)then
              mtrig=jm-1
              C(1)=C(2)
              D(1)=D(2)-A(2)*Zopb(i)
              do j=2,mtrig
                A(j)=A(j+1)
                C(j)=C(j+1)
                D(j)=D(j+1)
              enddo

              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i))=Ztrig(j)
              enddo
              Z(i,jbeg(i))=Zopb(i)

            endif
c -- top elevation boundary
            if((jtop.eq.23.or.jtop.eq.28.or.jtop.eq.27).and.
     *         btype_bot.eq.0)then
              mtrig=jm-1
              D(mtrig)=D(mtrig)-C(mtrig)*Zopt(i)
              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i)-1)=Ztrig(j)
              enddo
              Z(i,jend(i))=Zopt(i)

            endif
c -- bottom and top elevation boundary
            if(mod(int(jtop/10),10).eq.2.and.
     *        mod(int(jbot/10),10).eq.2)then
              mtrig=jm-2
              C(1)=C(2)
              D(1)=D(2)-A(2)*Zopb(i)
              do j=2,mtrig-1
                A(j)=A(j+1)
                C(j)=C(j+1)
                D(j)=D(j+1)
              enddo
              A(mtrig)=A(mtrig+1)
              D(mtrig)=D(mtrig+1)-C(mtrig+1)*Zopt(i)

              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i))=Ztrig(j)
              enddo
              Z(i,jbeg(i))=Zopb(i)
              Z(i,jend(i))=Zopt(i)
            endif
c --  left and right elevation boundary
            do j=jbeg(i),jend(i)
              Np=Ntype(i,j)
              if(Np.eq.24.or.Np.eq.25.or.Np.eq.28) Z(i,j)=Zopl(j)
              if(Np.eq.22.or.Np.eq.27.or.Np.eq.26) Z(i,j)=Zopr(j)
            enddo

c --- left  and right gradient boundary

            do j=jbeg(i),jend(i)
              Np=Ntype(i,j)
              if(Np.eq.44.or.Np.eq.45.or.Np.eq.48) Z(i,j)=Z(i+1,j)
              if(Np.eq.42.or.Np.eq.47.or.Np.eq.46) Z(i,j)=Z(i-1,j)
            enddo

c -- bottom gradient boundary
            if((jbot.eq.41.or.jbot.eq.45.or.jbot.eq.46).and.
     *         btype_top.eq.0)then
              mtrig=jm-1
              C(1)=C(2)/(1.+A(2))
              D(1)=D(2)/(1.+A(2))
              do j=2,mtrig
                A(j)=A(j+1)
                C(j)=C(j+1)
                D(j)=D(j+1)
              enddo

              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i))=Ztrig(j)
              enddo
              Z(i,jbeg(i))=Z(i,jbeg(i)+1)
            endif

c -- top gradient boundary
            if((jtop.eq.43.or.jtop.eq.48.or.jtop.eq.47).and.
     *         btype_bot.eq.0)then
              mtrig=jm-1
              A(mtrig)=A(mtrig)/(1.+C(mtrig))
              D(mtrig)=D(mtrig)/(1.+C(mtrig))
              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i)-1)=Ztrig(j)
              enddo
              Z(i,jend(i))=Z(i,jend(i)-1)
            endif

c -- bottom and top gradient boundary
            if(mod(int(jtop/10),10).eq.4.and.
     *        mod(int(jbot/10),10).eq.4)then
              mtrig=jm-2
              C(1)=C(2)/(1.+A(2))
              D(1)=D(2)/(1.+A(2))
              do j=2,mtrig-1
                A(j)=A(j+1)
                C(j)=C(j+1)
                D(j)=D(j+1)
              enddo

              A(mtrig)=A(mtrig+1)/(1.+C(mtrig+1))
              D(mtrig)=D(mtrig+1)/(1.+C(mtrig+1))

              call trig(A,C,D,Ztrig,mtrig)
              do j=1,mtrig
                Z(i,j+jbeg(i))=Ztrig(j)
              enddo
              Z(i,jbeg(i))=Z(i,jbeg(i)+1)
              Z(i,jend(i))=Z(i,jend(i)-1)
            endif

c -- z computed over
c -- computed VV
         do j=jbeg(i),jend(i)
           Np=Ntype(i,j)
           Np1=Ntype(i,j)
           if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
           if(Np.eq.0.or.Np.eq.14.or.Np.eq.18.or.Np.eq.13.or.Np.eq.17
     *     .or.Np.eq.12)then
             VV(i,j)=-Ea(j)/deta*(Z(i,j)-Z(i,j-1))+Eb(j)
           endif
           if(Np.eq.15.or.Np.eq.16.or.Np.eq.11)then
             VV(i,j)=0.
           endif
           if(Np1.eq.35.or.Np1.eq.36.or.Np1.eq.31)then
              VV(i,j)=Vopb(i)
           endif
           if(Np1.eq.38.or.Np1.eq.37.or.Np1.eq.33)then
              VV(i,j+1)=Vopt(i)
           endif
           if(Np1.eq.65.or.Np1.eq.66.or.Np1.eq.61)then
             VV(i,j)=-Ea(j)/deta*(Z(i,j)-Z(i,jend(i)))+Eb(j)             
           endif
           if(Np1.eq.68.or.Np1.eq.67.or.Np1.eq.63)then
            VV(i,j+1)=-Ea(j+1)/deta*(Z(i,jbeg(i))-Z(i,j))+Eb(j+1)             
           endif
         enddo
c -- OK

100     continue

c --- gradient boundary condition

	do i=1,m-1
        do j=jbeg(i),jend(i)
           Np1=Ntype(i,j)
           if(Np1.eq.45.or.Np1.eq.48.or.Np1.eq.44)then
             VV(i,j)=VV(i+1,j)
           endif
           if(Np1.eq.46.or.Np1.eq.47.or.Np1.eq.42)then
             VV(i,j)=VV(i-1,j)
           endif
        enddo
        enddo      

c -- set dry point and HH

	do i=1,m-1
        do j=jbeg(i),jend(i)
          HH(i,j)=H(i,j)+Z(i,j)	  
	  if(HH(i,j).lt.0.)then
            Z(i,j)=-H(i,j)-delta
	    HH(i,j)=H(i,j)+Z(i,j)
          endif
        enddo
        enddo	

c -- ok

        return
        end

        subroutine velocity
!        use sadi
        use sadi
        integer Np,Np1
c -- U
        do j=1,n-1
        do i=1,m
          Np=Ntype(i,j)
          if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
          if(Np.eq.0.or.Np.eq.11.or.Np.eq.16.or.Np.eq.12)then
          H_=0.5*(HH(i,j)+HH(i-1,j))
          if(H_.lt.hi)H_=hi
          HJ=H_*0.5*(rJ(i,j)+rJ(i-1,j))
          u(i,j)=0.5*(Xxi(i,j)+Xxi(i-1,j))*UU(i,j)/HJ
     *          +0.5*(Xeta(i,j)+Xeta(i-1,j))*0.25
     *          *(VV(i,j)+VV(i,j+1)+VV(i-1,j+1)+VV(i-1,j))/HJ
          endif
          if(Np.eq.13.or.Np.eq.17)then
          H_=0.5*(HH(i,j)+HH(i-1,j))
          if(H_.lt.hi)H_=hi
          HJ=H_*0.5*(rJ(i,j)+rJ(i-1,j))
          u(i,j)=0.5*(Xxi(i,j)+Xxi(i-1,j))*UU(i,j)/HJ
     *          +0.5*(Xeta(i,j)+Xeta(i-1,j))*0.25
     *          *(VV(i,j)+0.+0.+VV(i-1,j))/HJ
          endif
          if(Np.eq.14.or.Np.eq.15)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i,j)
          u(i,j)=Xxi(i,j)*UU(i,j)/HJ
     *          +Xeta(i,j)*0.5
     *          *(VV(i,j)+VV(i,j+1))/HJ
          endif
          if(Np.eq.18)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i,j)
          u(i,j)=Xxi(i,j)*UU(i,j)/HJ
     *          +Xeta(i,j)*0.5
     *          *(VV(i,j)+0.)/HJ
          endif
          if(Np.eq.88.or.Np.eq.83.or.Np.eq.87)then
          H_=HH(i,j-1)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i,j-1)
          u(i,j)=0.
     *          +Xeta(i,j-1)*0.5
     *          *(VV(i,j)+0.)/HJ
          endif

          if(Np.eq.86.or.Np.eq.82)then
          H_=HH(i-1,j)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i-1,j)
          u(i,j)=Xxi(i-1,j)*UU(i,j)/HJ
          endif
        enddo
        enddo

c --- u at c points

        do j=1,n
        do i=1,m

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	  
         
          if(Np.eq.0.or.Np.eq.2.or.Np.eq.7.or.Np.eq.92
     *      .or.Np.eq.3.or.Np.eq.8.or.Np.eq.4.or.Np.eq.9)then   

            uconv(i,j)=0.5*(u(i,j)+u(i,j-1))

          elseif(Np.eq.5.or.Np.eq.1.or.Np.eq.6.or.Np.eq.96)then

            uconv(i,j)=0.5*(3.*u(i,j)-u(i,j+1))

          elseif(Np.eq.98.or.Np.eq.93.or.Np.eq.97)then

            uconv(i,j)=0.5*(3.*u(i,j-1)-u(i,j-2))

          endif

        enddo
        enddo


c -- V
        do j=1,n
        do i=1,m-1
          Np=Ntype(i,j)
          if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
          if(Np.eq.0.or.Np.eq.14.or.Np.eq.18.or.Np.eq.13)then
          H_=0.5*(HH(i,j)+HH(i,j-1))
          if(H_.lt.hi)H_=hi
          HJ=H_*0.5*(rJ(i,j)+rJ(i,j-1))
          v(i,j)=0.5*(Yeta(i,j)+Yeta(i,j-1))*VV(i,j)/HJ
     *          +0.5*(Yxi(i,j)+Yxi(i,j-1))*0.25
     *          *(UU(i,j)+UU(i+1,j)+UU(i+1,j-1)+UU(i,j-1))/HJ
          endif
          if(Np.eq.12.or.Np.eq.17)then
          H_=0.5*(HH(i,j)+HH(i,j-1))
          if(H_.lt.hi)H_=hi
          HJ=H_*0.5*(rJ(i,j)+rJ(i,j-1))
          v(i,j)=0.5*(Yeta(i,j)+Yeta(i,j-1))*VV(i,j)/HJ
     *          +0.5*(Yxi(i,j)+Yxi(i,j-1))*0.25
     *          *(UU(i,j)+0.+0.+UU(i,j-1))/HJ
          endif
          if(Np.eq.11.or.Np.eq.15)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i,j)
          v(i,j)=Yeta(i,j)*VV(i,j)/HJ
     *          +Yxi(i,j)*0.5
     *          *(UU(i,j)+UU(i+1,j))/HJ
          endif
          if(Np.eq.18)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i,j)
          v(i,j)=Yeta(i,j)*VV(i,j)/HJ
     *          +Yxi(i,j)*0.5
     *          *(UU(i,j)+0.)/HJ
          endif
          if(Np.eq.88.or.Np.eq.83)then
          H_=HH(i,j-1)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i,j-1)
          v(i,j)=Yeta(i,j-1)*VV(i,j)/HJ
          endif

          if(Np.eq.87)v(i,j)=0.

          if(Np.eq.86.or.Np.eq.82)then
          H_=HH(i-1,j)
          if(H_.lt.hi)H_=hi
          HJ=H_*rJ(i-1,j)
          v(i,j)=
     *          +Yxi(i-1,j)*0.5
     *          *(UU(i,j)+0.)/HJ
          endif

        enddo
        enddo

c v at c points

        do j=1,n
        do i=1,m

 	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	            

          if(Np.eq.0.or.Np.eq.3.or.Np.eq.7.or.Np.eq.93
     *      .or.Np.eq.2.or.Np.eq.6.or.Np.eq.1)then ! change 19

            vconv(i,j)=0.5*(v(i,j)+v(i-1,j))

	  elseif(Np.eq.9)then

	    vconv(i,j)=0.

          elseif(Np.eq.5.or.Np.eq.4.or.Np.eq.8.or.Np.eq.98)then

            vconv(i,j)=0.5*(3.*v(i,j)-v(i+1,j))

          elseif(Np.eq.96.or.Np.eq.92.or.Np.eq.97)then

            vconv(i,j)=0.5*(3.*v(i-1,j)-v(i-2,j))

          endif

        enddo
        enddo

        return
        end
c ---------------------------------
c     open boundary
c ---------------------------------
        subroutine openb(timecirc)
c        use sadi
        use sadi
        real timecirc,time_real
        real Ht,Gt,T,Ubkg,Vbkg,phase,ele0
        integer kkk

c --- only for sf0 case
      real fac,val,omega,deg2rad
      real pha0(1:ntide),fac0(1:ntide)
      real periodm(1:ntide)
      data pha0/230.8,8.3,223.5,0.0,349.3,352.6,345.3,196.2,201.5,280.8/
      data fac0/0.9650,1.108,1.176,1.000,1.000,0.965,1.176,1.299,1.,1./
      data periodm/12.42,23.934,25.819,12.0,                   
     &   24.066,12.66,26.868,11.97,4382.91,8765.82/
      ele0=0.97

          deg2rad=pi/180.
c       starting date 06/18/2005, dates=31+28+31+30+31+18=169
	time_real=timecirc+169.*24.*3600.


        Ht=0.5
	Gt=90.*3.14/180.
	T=12.42*3600.
        phase=5400./sqrt(9.8*2.5)/T*2.*3.14
c ---  flux m^2/s
        Ubkg=-0.001*sin(2.*3.14159/T*timecirc)
        Vbkg=-0.50*sin(2.*3.14159/T*timecirc)
        
        if(east_data1_anly2.eq.2.or.
     &     west_data1_anly2.eq.2.or.
     &     south_data1_anly2.eq.2.or.
     &     north_data1_anly2.eq.2) then
! --- east elevation 
        if(east_ele2_flx3_grd4_rad5_per6.eq.2)then
	do j=1,n-1
	  Zopr(j)=Ht*sin(2.*3.14159/T*timecirc)
	enddo
	Zopt(m-1)=Ht*sin(2.*3.14159/T*timecirc)
	Zopb(m-1)=Ht*sin(2.*3.14159/T*timecirc)
        endif

! --- west elevation
        if(west_ele2_flx3_grd4_rad5_per6.eq.2)then
c	do j=1,n-1
c	  Zopl(j)=Ht*sin(2.*3.14159/T*timecirc)
c	enddo
c	Zopt(1)=Zopl(n-1)
c	Zopb(1)=Zopl(1)
c        endif

c -- only for sf0 case
         do j=1,n-1
         Zopl(j)=0.
          DO kkk=1,ntide
          fac=fac0(kkk)
          val=amp_west(j,kkk)
          omega=2.0*pi*time_real/(periodm(kkk)*3600.0)                    
     &      +pha0(kkk)*deg2rad
          Zopl(j)= Zopl(j)+                                      
     &       fac*val*COS(omega-pha_west(j,kkk)*deg2rad)
        END DO
c          Zopl(j)=(Zopl(j)+ele0)*TANH(timecirc/3600./24.)    
           Zopl(j)=0.
         enddo
	Zopt(1)=Zopl(n-1)
        Zopb(1)=Zopl(1)
         endif

! --- south elevation
        if(south_ele2_flx3_grd4_rad5_per6.eq.2)then
c	do i=1,m-1
c	  Zopb(i)=Ht*sin(2.*3.14159/T*timecirc)
c	enddo
c	Zopl(1)=Zopb(1)
c	Zopr(1)=Zopb(m-1)
c        endif

c - only for sf0 case
        do i=1,m-1
          Zopb(i)=0.
        DO kkk=1,ntide
          fac=fac0(kkk)
          val=amp_south(i,kkk)
          omega=2.0*pi*time_real/(periodm(kkk)*3600.0)
     &      +pha0(kkk)*deg2rad
          Zopb(i)= Zopb(i)+
     &       fac*val*COS(omega-pha_south(i,kkk)*deg2rad)
        END DO
c          Zopb(i)=(Zopb(i)+ele0)*TANH(timecirc/3600./24.)
           Zopb(i)=0.
        enddo
        Zopl(1)=Zopb(1)
        Zopr(1)=Zopb(m-1)

        endif

! --- north elevation
        if(north_ele2_flx3_grd4_rad5_per6.eq.2)then
c	do i=1,m-1
c	  Zopt(i)=Ht*sin(2.*3.14159/T*timecirc)
c	enddo
c	Zopl(n-1)=Zopt(1)
c	Zopr(n-1)=Zopr(m-1)
c        endif

c -- only for sf0 case

c       do i=istart_N,iend_N
        do i=1,m-1
          Zopt(i)=0.
        DO kkk=1,ntide
          fac=fac0(kkk)
          val=amp_north(i,kkk)
          omega=2.0*pi*time_real/(periodm(kkk)*3600.0)
     &      +pha0(kkk)*deg2rad
          Zopt(i)= Zopt(i)+
     &       fac*val*COS(omega-pha_north(i,kkk)*deg2rad)
        END DO
c          Zopt(i)=(Zopt(i)+ele0)*TANH(timecirc/3600./24.)
           Zopt(i)=0.           
        enddo
        Zopl(n-1)=Zopt(1)
        Zopr(n-1)=Zopt(m-1)


        endif	

! --- east flux

        if(east_ele2_flx3_grd4_rad5_per6.eq.3)then
	do j=1,n-1
	  Uopr(j)=Ubkg*Yeta(m-1,j)
     &         
	enddo
          if(south_ele2_flx3_grd4_rad5_per6.eq.3)then
            Vopb(m-1)=Vbkg*Xxi(m-1,1)
          else
            Vopb(m-1)=0.
          endif
          if(north_ele2_flx3_grd4_rad5_per6.eq.3)then
            Vopt(m-1)=Vbkg*Xxi(m-1,n-1)
          else
            Vopb(m-1)=0.
          endif
        endif

! --- west flux

        if(west_ele2_flx3_grd4_rad5_per6.eq.3)then
	do j=1,n-1
	  Uopl(j)=Ubkg*Yeta(1,j)
	enddo
          if(south_ele2_flx3_grd4_rad5_per6.eq.3)then 
            Vopb(1)=Vbkg*Xxi(1,1)
          else
            Vopb(1)=0.
          endif
          if(north_ele2_flx3_grd4_rad5_per6.eq.3)then 
            Vopt(1)=Vbkg*Xxi(1,n-1)
          else
            Vopt(1)=0.
          endif
        endif

! --- south flux

        if(south_ele2_flx3_grd4_rad5_per6.eq.3)then
	do i=1,m-1
	  Vopb(i)=Vbkg*Xxi(i,1)
	enddo
          if(west_ele2_flx3_grd4_rad5_per6.eq.3)then
           Uopl(1)=Ubkg*Yeta(1,1)
          else
	   Uopl(1)=0.
          endif
          if(east_ele2_flx3_grd4_rad5_per6.eq.3)then
           Uopr(1)=Ubkg*Yeta(m-1,1)
          else
	   Uopr(1)=0.
          endif
        endif

! --- north flux

        if(north_ele2_flx3_grd4_rad5_per6.eq.3)then
	do i=1,m-1
	  Vopt(i)=Vbkg*Xxi(i,n-1)
     &         
	enddo
          if(west_ele2_flx3_grd4_rad5_per6.eq.3)then
            Uopl(n-1)=Ubkg*Yeta(1,n-1)
          else
	    Uopl(n-1)=0.
          endif
          if(east_ele2_flx3_grd4_rad5_per6.eq.3)then
            Uopr(n-1)=Ubkg*Yeta(m-1,n-1)
          else
	    Uopr(n-1)=0.
          endif
        endif

        else

        write(*,*)'Please specify your data format'
        stop

        endif


        return
        end

c ---------------------------------

        subroutine outputcirc(mpr,npr,num_file,varb)
!        use sadi
        use sadi
        integer mpr,npr
        integer nm_first,nm_second,nm_third,nm_fourth,num_file
        real varb(mtrim,ntrim)
        character*4 file_name

        nm_first=mod(num_file/1000,10)
        nm_second=mod(num_file/100,10)
        nm_third=mod(num_file/10,10)
        nm_fourth=mod(num_file,10)


        write(file_name(1:1),'(I1)')nm_first
        write(file_name(2:2),'(I1)')nm_second
        write(file_name(3:3),'(I1)')nm_third
        write(file_name(4:4),'(I1)')nm_fourth

        open(2,file='s'//file_name//'.out')
        do i = 1, mpr 
         write(2,100)(varb(i,j), j=1,npr ) 
        end do 
  
100     format(501(f20.8))
        close(2)
	
	return

        end 

	subroutine interp_from_z_to_c(nx,ny,Ntp,var,var_c)
!	use sadi
        use sadi
	integer Ntp(mtrim,ntrim),nx,ny,Np,Np1
	real var(mtrim,ntrim),var_c(mtrim,ntrim)
        real zz1,zz2


! --- zeta_c

        do j=1,ny
        do i=1,nx
          Np=Ntp(i,j)
            if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)

          if(Np.eq.12.or.Np.eq.17.or.Np.eq.13)then

            var_c(i,j)=0.25*
     *             (var(i,j)+var(i-1,j)+var(i-1,j-1)+var(i,j-1))

	  elseif(Np.eq.0)then
	    Np1=Ntp(i-1,j-1)
	    if(Np1.ge.80)then

            var_c(i,j)=1./3.*
     *             (var(i,j)+var(i-1,j)+var(i,j-1))

	    else

            var_c(i,j)=0.25*
     *             (var(i,j)+var(i-1,j)+var(i-1,j-1)+var(i,j-1))

	    endif

          elseif(Np.eq.11.or.Np.eq.16)then

            var_c(i,j)=0.5*
     *                  (0.5*(3.*var(i-1,j)-var(i-1,j+1))
     *                  +0.5*(3.*var(i,j)-var(i,j+1)))

          elseif(Np.eq.92)then

            var_c(i,j)=0.5*
     *                  (0.5*(3.*var(i-1,j)-var(i-2,j))
     *                  +0.5*(3.*var(i-1,j-1)-var(i-2,j-1)))

          elseif(Np.eq.93)then

            var_c(i,j)=0.5*
     *                  (0.5*(3.*var(i,j-1)-var(i,j-2))
     *                  +0.5*(3.*var(i-1,j-1)-var(i-1,j-2)))

          elseif(Np.eq.14.or.Np.eq.18)then

            var_c(i,j)=0.5*
     *                  (0.5*(3.*var(i,j)-var(i+1,j))
     *                  +0.5*(3.*var(i,j-1)-var(i+1,j-1)))

          elseif(Np.eq.15)then
	
	    zz1=0.5*(3.*var(i,j)-var(i,j+1))
	    zz2=0.5*(3.*var(i+1,j)-var(i+1,j+1))
            var_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.96)then

	    zz1=0.5*(3.*var(i-1,j)-var(i-1,j+1))
	    zz2=0.5*(3.*var(i-2,j)-var(i-2,j+1))
            var_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.97)then

	    zz1=0.5*(3.*var(i-1,j-1)-var(i-1,j-2))
	    zz2=0.5*(3.*var(i-2,j-1)-var(i-2,j-2))
            var_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.98)then

	    zz1=0.5*(3.*var(i,j-1)-var(i,j-2))
	    zz2=0.5*(3.*var(i+1,j-1)-var(i+1,j-2))
            var_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.19)then                                           

!            !zeta_c(i,j)=zeta(i-1,j-1)  !change?

	    var_c(i,j)=0.5*
     *                  (0.5*(3.*var(i-1,j)-var(i-2,j))
     *                  +0.5*(3.*var(i-1,j-1)-var(i-2,j-1)))

          endif

        enddo
        enddo

	return

	end 

! -----------------------
        subroutine get_gxx_Cf
!        use sadi
        use sadi
        integer Np,Np1

!--- initialize g_11 ...

	do j=1,n
	do i=1,m
	sq_g_0(i,j)=1.
	sq_g_0_z(i,j)=1.
	sq_g_0_u(i,j)=1.
	sq_g_0_v(i,j)=1.
	g_0(i,j)=1.
	g_11(i,j)=1.
	g_12(i,j)=0.
	g_22(i,j)=1.
        Cf111(i,j)=0.
        Cf211(i,j)=0.
        Cf212(i,j)=0.
        Cf112(i,j)=0.
        Cf222(i,j)=0.
        Cf122(i,j)=0.
	X_xi_1(i,j)=0.
	X_xi_2(i,j)=0.
	Y_xi_1(i,j)=0.
	Y_xi_2(i,j)=0.
	enddo
	enddo
	        
!--- calculate sq_g_0, g_11, g_12, g_22 
      
        do j=1,n
        do i=1,m

	  if(Ntype1(i,j).gt.80)then
            Np=Ntype1(i,j)
	  else
	    Np=mod(Ntype1(i,j),10)
	  endif

          if(Ntype1(i,j).eq.0.or.Ntype1(i,j).eq.19
     *       .or.Np.eq.2.or.Np.eq.3
     *       .or.Np.eq.7)then

            sq_g_0(i,j)=(x(i+1,j)-x(i-1,j))/2.*(y(i,j+1)-y(i,j-1))/2.
     *               -(x(i,j+1)-x(i,j-1))/2.*(y(i+1,j)-y(i-1,j))/2.

            g_11(i,j)=((x(i+1,j)-x(i-1,j))/2.)**2
     *               +((y(i+1,j)-y(i-1,j))/2.)**2

            g_22(i,j)=((x(i,j+1)-x(i,j-1))/2.)**2
     *               +((y(i,j+1)-y(i,j-1))/2.)**2

            g_12(i,j)=0.5*(x(i+1,j)-x(i-1,j))*0.5*(x(i,j+1)-x(i,j-1))
     *               +0.5*(y(i+1,j)-y(i-1,j))*0.5*(y(i,j+1)-y(i,j-1))

!--- for Ntype1=11,16
          elseif(Np.eq.1.or.Np.eq.6
     *           .or.Ntype1(i,j).eq.91)then

            sq_g_0(i,j)=(x(i+1,j)-x(i-1,j))/2.*(y(i,j+1)-y(i,j))/1.
     *               -(x(i,j+1)-x(i,j))/1.*(y(i+1,j)-y(i-1,j))/2.

            g_11(i,j)=((x(i+1,j)-x(i-1,j))/2.)**2
     *               +((y(i+1,j)-y(i-1,j))/2.)**2

            g_22(i,j)=((x(i,j+1)-x(i,j))/1.)**2
     *               +((y(i,j+1)-y(i,j))/1.)**2

            g_12(i,j)=0.5*(x(i+1,j)-x(i-1,j))*1.0*(x(i,j+1)-x(i,j))
     *               +0.5*(y(i+1,j)-y(i-1,j))*1.0*(y(i,j+1)-y(i,j))

! --- for Ntype1=92
          elseif(Ntype1(i,j).eq.92.or.Ntype1(i,j).eq.82)then

            sq_g_0(i,j)=(x(i,j)-x(i-1,j))/1.*(y(i,j+1)-y(i,j-1))/2.
     *               -(x(i,j+1)-x(i,j-1))/2.*(y(i,j)-y(i-1,j))/1.

            g_11(i,j)=((x(i,j)-x(i-1,j))/1.)**2
     *               +((y(i,j)-y(i-1,j))/1.)**2

            g_22(i,j)=((x(i,j+1)-x(i,j-1))/2.)**2
     *               +((y(i,j+1)-y(i,j-1))/2.)**2

            g_12(i,j)=1.0*(x(i,j)-x(i-1,j))*0.5*(x(i,j+1)-x(i,j-1))
     *               +1.0*(y(i,j)-y(i-1,j))*0.5*(y(i,j+1)-y(i,j-1))

! --- for Ntype1=93
          elseif(Ntype1(i,j).eq.93.or.Ntype1(i,j).eq.83)then

            sq_g_0(i,j)=(x(i+1,j)-x(i-1,j))/2.*(y(i,j)-y(i,j-1))/1.
     *               -(x(i,j)-x(i,j-1))/1.*(y(i+1,j)-y(i-1,j))/2.

            g_11(i,j)=((x(i+1,j)-x(i-1,j))/2.)**2
     *               +((y(i+1,j)-y(i-1,j))/2.)**2

            g_22(i,j)=((x(i,j)-x(i,j-1))/1.)**2
     *               +((y(i,j)-y(i,j-1))/1.)**2

            g_12(i,j)=0.5*(x(i+1,j)-x(i-1,j))*1.0*(x(i,j)-x(i,j-1))
     *               +0.5*(y(i+1,j)-y(i-1,j))*1.0*(y(i,j)-y(i,j-1))

! --- for Ntype1=14,18
          elseif(Np.eq.4.or.Np.eq.8
     *           .or.Ntype1(i,j).eq.94.or.Ntype1(i,j).eq.84)then

            sq_g_0(i,j)=(x(i+1,j)-x(i,j))/1.*(y(i,j+1)-y(i,j-1))/2.
     *               -(x(i,j+1)-x(i,j-1))/2.*(y(i+1,j)-y(i,j))/1.

            g_11(i,j)=((x(i+1,j)-x(i,j))/1.)**2
     *               +((y(i+1,j)-y(i,j))/1.)**2

            g_22(i,j)=((x(i,j+1)-x(i,j-1))/2.)**2
     *               +((y(i,j+1)-y(i,j-1))/2.)**2

            g_12(i,j)=1.0*(x(i+1,j)-x(i,j))*0.5*(x(i,j+1)-x(i,j-1))
     *               +1.0*(y(i+1,j)-y(i,j))*0.5*(y(i,j+1)-y(i,j-1))

! --- for Ntype1=15
          elseif(Np.eq.5.or.Ntype1(i,j).eq.95)then

            sq_g_0(i,j)=(x(i+1,j)-x(i,j))/1.*(y(i,j+1)-y(i,j))/1.
     *               -(x(i,j+1)-x(i,j))/1.*(y(i+1,j)-y(i,j))/1.

            g_11(i,j)=((x(i+1,j)-x(i,j))/1.)**2
     *               +((y(i+1,j)-y(i,j))/1.)**2

            g_22(i,j)=((x(i,j+1)-x(i,j))/1.)**2
     *               +((y(i,j+1)-y(i,j))/1.)**2

            g_12(i,j)=1.0*(x(i+1,j)-x(i,j))*1.0*(x(i,j+1)-x(i,j))
     *               +1.0*(y(i+1,j)-y(i,j))*1.0*(y(i,j+1)-y(i,j))

! --- for Ntype1=96
          elseif(Ntype1(i,j).eq.96.or.Ntype1(i,j).eq.86)then

            sq_g_0(i,j)=(x(i,j)-x(i-1,j))/1.*(y(i,j+1)-y(i,j))/1.
     *               -(x(i,j+1)-x(i,j))/1.*(y(i,j)-y(i-1,j))/1.

            g_11(i,j)=((x(i,j)-x(i-1,j))/1.)**2
     *               +((y(i,j)-y(i-1,j))/1.)**2

            g_22(i,j)=((x(i,j+1)-x(i,j))/1.)**2
     *               +((y(i,j+1)-y(i,j))/1.)**2

            g_12(i,j)=1.0*(x(i,j)-x(i-1,j))*1.0*(x(i,j+1)-x(i,j))
     *               +1.0*(y(i,j)-y(i-1,j))*1.0*(y(i,j+1)-y(i,j))

! --- for Ntype1=97
          elseif(Ntype1(i,j).eq.97.or.Ntype1(i,j).eq.87)then

            sq_g_0(i,j)=(x(i,j)-x(i-1,j))/1.*(y(i,j)-y(i,j-1))/1.
     *               -(x(i,j)-x(i,j-1))/1.*(y(i,j)-y(i-1,j))/1.

            g_11(i,j)=((x(i,j)-x(i-1,j))/1.)**2
     *               +((y(i,j)-y(i-1,j))/1.)**2

            g_22(i,j)=((x(i,j)-x(i,j-1))/1.)**2
     *               +((y(i,j)-y(i,j-1))/1.)**2

            g_12(i,j)=1.0*(x(i,j)-x(i-1,j))*1.0*(x(i,j)-x(i,j-1))
     *               +1.0*(y(i,j)-y(i-1,j))*1.0*(y(i,j)-y(i,j-1))

! --- for Ntype1=98
          elseif(Ntype1(i,j).eq.98.or.Ntype1(i,j).eq.88)then

            sq_g_0(i,j)=(x(i+1,j)-x(i,j))/1.*(y(i,j)-y(i,j-1))/1.
     *               -(x(i,j)-x(i,j-1))/1.*(y(i+1,j)-y(i,j))/1.

            g_11(i,j)=((x(i+1,j)-x(i,j))/1.)**2
     *               +((y(i+1,j)-y(i,j))/1.)**2

            g_22(i,j)=((x(i,j)-x(i,j-1))/1.)**2
     *               +((y(i,j)-y(i,j-1))/1.)**2

            g_12(i,j)=1.0*(x(i+1,j)-x(i,j))*1.0*(x(i,j)-x(i,j-1))
     *               +1.0*(y(i+1,j)-y(i,j))*1.0*(y(i,j)-y(i,j-1))

          endif

        enddo
        enddo

! --- don't consider non-orthogonality

	if(orth_grid.eq.1)then
	do j=1,n
	do i=1,m
	  g_12(i,j)=0.
	enddo
	enddo
	endif
	

! ---  get g_0

	do i=1,m
	do j=1,n
	g_0(i,j) = sq_g_0(i,j)*sq_g_0(i,j)
	enddo
	enddo

! ---   get _u

        do j=1,n-1
        do i=1,m

          Np=Ntype1(i,j)

          if(Np.ne.93.and.Np.ne.98.and.Np.ne.97
     *      .and.Np.ne.83.and.Np.ne.88.and.Np.ne.87)then
            sq_g_0_u(i,j)=0.5*(sq_g_0(i,j+1)+sq_g_0(i,j))
            g_11_u(i,j)=0.5*(g_11(i,j+1)+g_11(i,j))
            g_12_u(i,j)=0.5*(g_12(i,j+1)+g_12(i,j))
            g_22_u(i,j)=0.5*(g_22(i,j+1)+g_22(i,j))
          endif

          if(Np.eq.93.or.Np.eq.98.or.Np.eq.97
     *      .or.Np.eq.83.or.Np.eq.88.or.Np.eq.87)then
            sq_g_0_u(i,j)=sq_g_0(i,j)
            g_11_u(i,j)=g_11(i,j)
            g_12_u(i,j)=g_12(i,j)
            g_22_u(i,j)=g_22(i,j)
          endif

        enddo
        enddo

! ---   get _v

        do j=1,n
        do i=1,m-1

          Np=Ntype1(i,j)

          if(Np.ne.92.and.Np.ne.96.and.Np.ne.97
     *      .and.Np.ne.82.and.Np.ne.86.and.Np.ne.87)then
            sq_g_0_v(i,j)=0.5*(sq_g_0(i+1,j)+sq_g_0(i,j))
            g_11_v(i,j)=0.5*(g_11(i+1,j)+g_11(i,j))
            g_12_v(i,j)=0.5*(g_12(i+1,j)+g_12(i,j))
            g_22_v(i,j)=0.5*(g_22(i+1,j)+g_22(i,j))
          endif

          if(Np.eq.92.or.Np.eq.96.or.Np.eq.97
     *      .or.Np.eq.82.or.Np.eq.86.or.Np.eq.87)then
            sq_g_0_v(i,j)=sq_g_0(i,j)
            g_11_v(i,j)=g_11(i,j)
            g_12_v(i,j)=g_12(i,j)
            g_22_v(i,j)=g_22(i,j)
          endif

        enddo
        enddo

! ---   get sq_g_0_z

        do j=1,n-1
        do i=1,m-1

          sq_g_0_z(i,j)=0.5*(sq_g_0_u(i+1,j)+sq_g_0_u(i,j))
          g_11_z(i,j)=0.5*(g_11_u(i+1,j)+g_11_u(i,j))
          g_22_z(i,j)=0.5*(g_22_u(i+1,j)+g_22_u(i,j))
          g_12_z(i,j)=0.5*(g_12_u(i+1,j)+g_12_u(i,j))

        enddo
        enddo  

	do j=1,n
	do i=1,m
	g_0_z(i,j)=sq_g_0_z(i,j)*sq_g_0_z(i,j)
	enddo
	enddo
        
! ---   Cristoffels were made at + grid point

! --- g_11_xi

	do j=1,n
	do i=1,m
        tr1(i,j)=g_11_v(i,j)
	enddo
	enddo

        call der_xi_at_c_symmetry

        do j=1,n
        do i=1,m
          if(Ntype1(i,j).ne.99) then
          Cf111(i,j)=0.5*g_22(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf211(i,j)=-0.5*g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
	  endif
        enddo
        enddo

! --- g12_xi

	do j=1,n
	do i=1,m
        tr1(i,j)=g_12_v(i,j)
	enddo
	enddo
 
        call der_xi_at_c_symmetry

        do j=1,n
        do i=1,m
	  if(Ntype1(i,j).ne.99) then
          Cf111(i,j)=Cf111(i,j)-g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf211(i,j)=Cf211(i,j)+g_11(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          endif
        enddo
        enddo

! --- g11_et

        do j=1,n
        do i=1,m
        tr1(i,j)=g_11_u(i,j)
	enddo
	enddo

        call der_eta_at_c_symmetry

        do j=1,n
        do i=1,m
	  if(Ntype1(i,j).ne.99) then
          Cf111(i,j)=Cf111(i,j)+0.5*g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf211(i,j)=Cf211(i,j)-0.5*g_11(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf212(i,j)=-0.5*g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf112(i,j)=0.5*g_22(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          endif
        enddo
        enddo

! --- g22_xi

        do j=1,n
        do i=1,m
        tr1(i,j)=g_22_v(i,j)
	enddo
	enddo

        call der_xi_at_c_symmetry

        do j=1,n
        do i=1,m
	  if(Ntype1(i,j).ne.99) then
          Cf212(i,j)=Cf212(i,j)+0.5*g_11(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf112(i,j)=Cf112(i,j)-0.5*g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf222(i,j)=0.5*g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf122(i,j)=-0.5*g_22(i,j)/sq_g_0(i,j)**2*trd1(i,j)
	  endif
        enddo
        enddo

! --- g12_et

        do j=1,n
        do i=1,m
        tr1(i,j)=g_12_u(i,j)
	enddo
	enddo

        call der_eta_at_c_symmetry

        do j=1,n
        do i=1,m
	  if(Ntype1(i,j).ne.99) then
          Cf222(i,j)=Cf222(i,j)-g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf122(i,j)=Cf122(i,j)+g_22(i,j)/sq_g_0(i,j)**2*trd1(i,j)
	  endif
        enddo
        enddo

! --- g22_et

        do j=1,n
        do i=1,m
        tr1(i,j)=g_22_u(i,j)
	enddo
	enddo

        call der_eta_at_c_symmetry

        do j=1,n
        do i=1,m
	  if(Ntype1(i,j).ne.99) then
          Cf222(i,j)=Cf222(i,j)+0.5*g_11(i,j)/sq_g_0(i,j)**2*trd1(i,j)
          Cf122(i,j)=Cf122(i,j)-0.5*g_12(i,j)/sq_g_0(i,j)**2*trd1(i,j)
	  endif
        enddo
        enddo

! --- calculate X_xi_1, X_xi_2, Y_xi_1, Y_xi_2

	do j=1,n
	do i=1,m
          Np=Ntype1(i,j)
	if(Np.le.80.)Np=mod(Np,10)
	if(Ntype1(i,j).eq.0.or.Np.eq.2.or.Np.eq.3
     *    .or.Np.eq.7.or.Np.eq.9)then
	  X_xi_1(i,j)=0.5*(x(i+1,j)-x(i-1,j))
	  X_xi_2(i,j)=0.5*(x(i,j+1)-x(i,j-1))
	  Y_xi_1(i,j)=0.5*(y(i+1,j)-y(i-1,j))
	  Y_xi_2(i,j)=0.5*(y(i,j+1)-y(i,j-1))
	elseif(Np.eq.4.or.Np.eq.8)then
	  X_xi_1(i,j)=x(i+1,j)-x(i,j)
	  X_xi_2(i,j)=0.5*(x(i,j+1)-x(i,j-1))
	  Y_xi_1(i,j)=y(i+1,j)-y(i,j)
	  Y_xi_2(i,j)=0.5*(y(i,j+1)-y(i,j-1))
	elseif(Ntype1(i,j).eq.92.or.Ntype1(i,j).eq.82)then
	  X_xi_1(i,j)=x(i,j)-x(i-1,j)
	  X_xi_2(i,j)=0.5*(x(i,j+1)-x(i,j-1))
	  Y_xi_1(i,j)=y(i,j)-y(i-1,j)
	  Y_xi_2(i,j)=0.5*(y(i,j+1)-y(i,j-1))
	elseif(Np.eq.1.or.Np.eq.6)then
	  X_xi_1(i,j)=0.5*(x(i+1,j)-x(i-1,j))
	  X_xi_2(i,j)=x(i,j+1)-x(i,j)
	  Y_xi_1(i,j)=0.5*(y(i+1,j)-y(i-1,j))
	  Y_xi_2(i,j)=y(i,j+1)-y(i,j)
	elseif(Ntype1(i,j).eq.93.or.Ntype1(i,j).eq.83)then
	  X_xi_1(i,j)=0.5*(x(i+1,j)-x(i-1,j))
	  X_xi_2(i,j)=x(i,j)-x(i,j-1)
	  Y_xi_1(i,j)=0.5*(y(i+1,j)-y(i-1,j))
	  Y_xi_2(i,j)=y(i,j)-y(i,j-1)
	elseif(Np.eq.5)then
	  X_xi_1(i,j)=x(i+1,j)-x(i,j)	  
	  X_xi_2(i,j)=x(i,j+1)-x(i,j)
	  Y_xi_1(i,j)=y(i+1,j)-y(i,j)
	  Y_xi_2(i,j)=y(i,j+1)-y(i,j)
	elseif(Ntype1(i,j).eq.96.or.Ntype1(i,j).eq.86)then
	  X_xi_1(i,j)=x(i,j)-x(i-1,j)
	  X_xi_2(i,j)=x(i,j+1)-x(i,j)
	  Y_xi_1(i,j)=y(i,j)-y(i-1,j)
	  Y_xi_2(i,j)=y(i,j+1)-y(i,j)
	elseif(Ntype1(i,j).eq.97.or.Ntype1(i,j).eq.87)then
	  X_xi_1(i,j)=x(i,j)-x(i-1,j)
	  X_xi_2(i,j)=x(i,j)-x(i,j-1)		
	  Y_xi_1(i,j)=y(i,j)-y(i-1,j)
	  Y_xi_2(i,j)=y(i,j)-y(i,j-1)
	elseif(Ntype1(i,j).eq.98.or.Ntype1(i,j).eq.88)then
	  X_xi_1(i,j)=x(i+1,j)-x(i,j)
	  X_xi_2(i,j)=x(i,j)-x(i,j-1)
	  Y_xi_1(i,j)=y(i+1,j)-y(i,j)
	  Y_xi_2(i,j)=y(i,j)-y(i,j-1)
	endif

	enddo
	enddo

	return
 
        end

! ------------------------------------------------------------------

        subroutine der_xi_at_c_symmetry
!        use sadi
        use sadi
        integer Np

        do j=1,n
        do i=1,m

	if(Ntype1(i,j).lt.80)then
        Np=mod(Ntype1(i,j),10)
	else
	Np=Ntype1(i,j)
	endif

        if(Np.eq.0.or.Np.eq.1.or.Np.eq.6.or.Np.eq.3.or.Np.eq.7
     *     .or.Np.eq.93.or.Np.eq.2)then

          trd1(i,j)=tr1(i,j)-tr1(i-1,j)

        elseif(Np.eq.4.or.Np.eq.5.or.Np.eq.8.or.Np.eq.98.or.Np.eq.88)
     *    then

          trd1(i,j)=0.

        elseif(Np.eq.92.or.Np.eq.96.or.Np.eq.97
     *        .or.Np.eq.82.or.Np.eq.86.or.Np.eq.87)then

          trd1(i,j)=0.

        elseif(Np.eq.9)then

          trd1(i,j)=tr1(i,j)-tr1(i-1,j)

        endif

        enddo
        enddo

	return

        end 

! -------------------------------------------------------------------

       subroutine der_eta_at_c_symmetry
!        use sadi
        use sadi
        integer Np

        do j=1,n
        do i=1,m

	if(Ntype1(i,j).lt.80)then
        Np=mod(Ntype1(i,j),10)
	else
	Np=Ntype1(i,j)
	endif

        if(Np.eq.0.or.Np.eq.4.or.Np.eq.8.or.Np.eq.2.or.Np.eq.7
     *     .or.Np.eq.92..or.Np.eq.82.or.Np.eq.3.or.Np.eq.9)then

          trd1(i,j)=tr1(i,j)-tr1(i,j-1)

        elseif(Np.eq.1.or.Np.eq.5.or.Np.eq.6.or.Np.eq.96.or.Np.eq.86)
     *    then

          trd1(i,j)=0.

        elseif(Np.eq.93.or.Np.eq.98.or.Np.eq.97
     *    .or.Np.eq.83.or.Np.eq.88.or.Np.eq.87)then

          trd1(i,j)=0.

        endif

        enddo
        enddo

	return

        end

! -------------------------------------

	subroutine depth_u_v_H
!	use sadi
        use sadi

	do j=1,n-1
	do i=1,m-1
	  depth(i,j)=0.25*(depth_c(i,j)+depth_c(i+1,j)
     *                    +depth_c(i,j+1)+depth_c(i+1,j+1))
	enddo
	enddo

	
	do j=1,n
	  depth(m,j)=depth(m-1,j)
	enddo
	do i=1,m
	  depth(i,n)=depth(i,n-1)
	enddo
	depth(m,n)=depth(m-1,n-1)

	do j=1,n-1
	do i=1,m-1
	  H(i,j)=depth(i,j)
	enddo
	enddo

	do j=1,n
	do i=1,m-1
	  depth_v(i,j)=0.5*(depth_c(i,j)+depth_c(i+1,j))
	enddo
	enddo

	do j=1,n-1
	do i=1,m
	  depth_u(i,j)=0.5*(depth_c(i,j)+depth_c(i,j+1))
	enddo
	enddo

	return

	end

! -------------------------------------

	subroutine convertUVZ
!	use sadi
        use sadi
        integer Np

	do j=1,n-1
	do i=1,m
	  Q1(i,j)=UU(i,j)/sq_g_0_u(i,j)
          Np=Ntype(i,j)
	  if(Np.eq.24.or.Np.eq.25.or.Np.eq.28)
     *     Q1(i,j)=UU(i+1,j)/sq_g_0_u(i+1,j)
	  if(Np.eq.44.or.Np.eq.45.or.Np.eq.48)
     *     Q1(i,j)=UU(i+1,j)/sq_g_0_u(i+1,j)
	  if(Np.eq.22.or.Np.eq.26.or.Np.eq.27)
     *     Q1(i,j)=UU(i-1,j)/sq_g_0_u(i-1,j)
	  if(Np.eq.82.or.Np.eq.86.or.Np.eq.87)
     *     Q1(i,j)=UU(i-2,j)/sq_g_0_u(i-2,j)
	enddo
	enddo

	do j=1,n
	do i=1,m-1
	  Q2(i,j)=VV(i,j)/sq_g_0_v(i,j)
          Np=Ntype(i,j)
	  if(Np.eq.21.or.Np.eq.25.or.Np.eq.26)
     *     Q2(i,j)=VV(i,j+1)/sq_g_0_v(i,j+1)
	  if(Np.eq.41.or.Np.eq.45.or.Np.eq.46)
     *     Q2(i,j)=VV(i,j+1)/sq_g_0_v(i,j+1)
	  if(Np.eq.23.or.Np.eq.28.or.Np.eq.27)
     *     Q2(i,j)=VV(i,j-1)/sq_g_0_v(i,j-1)
	  if(Np.eq.43.or.Np.eq.48.or.Np.eq.47)
     *     Q2(i,j)=VV(i,j-1)/sq_g_0_v(i,j-1)
	  if(Np.eq.83.or.Np.eq.88.or.Np.eq.87)
     *     Q2(i,j)=VV(i,j-2)/sq_g_0_v(i,j-2)
	enddo
	enddo

	do j=1,n-1
	do i=1,m-1
	  zeta(i,j)=Z(i,j)
	enddo
	enddo

	call calculate_zeta_uvc
	call calculate_Q1_zvc
	call calculate_Q2_zuc

	return
	end

! ----------------------------------------------------------------------------

	subroutine calculate_zeta_uvc
!	use sadi
        use sadi
        integer Np,Np1
        real zz1,zz2

! --- second order averaging

! --- zeta_c
        do j=1,n
        do i=1,m

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=Ntype1(i,j)
	  endif

          if(Np.eq.2.or.Np.eq.7.or.Np.eq.3)then

            zeta_c(i,j)=0.25*
     *             (zeta(i,j)+zeta(i-1,j)+zeta(i-1,j-1)+zeta(i,j-1))

	  elseif(Np.eq.0)then
	    Np1=Ntype1(i-1,j-1)
	    if(Np1.eq.99)then

            zeta_c(i,j)=1./3.*
     *             (zeta(i,j)+zeta(i-1,j)+zeta(i,j-1))

	    else

            zeta_c(i,j)=0.25*
     *             (zeta(i,j)+zeta(i-1,j)+zeta(i-1,j-1)+zeta(i,j-1))

	    endif

          elseif(Np.eq.1.or.Np.eq.6)then

            zeta_c(i,j)=0.5*
     *                  (0.5*(3.*zeta(i-1,j)-zeta(i-1,j+1))
     *                  +0.5*(3.*zeta(i,j)-zeta(i,j+1)))

          elseif(Np.eq.92.or.Np.eq.82)then

            zeta_c(i,j)=0.5*
     *                  (0.5*(3.*zeta(i-1,j)-zeta(i-2,j))
     *                  +0.5*(3.*zeta(i-1,j-1)-zeta(i-2,j-1)))

          elseif(Np.eq.93.or.Np.eq.83)then

            zeta_c(i,j)=0.5*
     *                  (0.5*(3.*zeta(i,j-1)-zeta(i,j-2))
     *                  +0.5*(3.*zeta(i-1,j-1)-zeta(i-1,j-2)))

          elseif(Np.eq.4.or.Np.eq.8)then

            zeta_c(i,j)=0.5*
     *                  (0.5*(3.*zeta(i,j)-zeta(i+1,j))
     *                  +0.5*(3.*zeta(i,j-1)-zeta(i+1,j-1)))

          elseif(Np.eq.5)then
	
	    zz1=0.5*(3.*zeta(i,j)-zeta(i,j+1))
	    zz2=0.5*(3.*zeta(i+1,j)-zeta(i+1,j+1))
            zeta_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.96.or.Np.eq.86)then

	    zz1=0.5*(3.*zeta(i-1,j)-zeta(i-1,j+1))
	    zz2=0.5*(3.*zeta(i-2,j)-zeta(i-2,j+1))
            zeta_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.97.or.Np.eq.87)then

	    zz1=0.5*(3.*zeta(i-1,j-1)-zeta(i-1,j-2))
	    zz2=0.5*(3.*zeta(i-2,j-1)-zeta(i-2,j-2))
            zeta_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.98.or.Np.eq.88)then

	    zz1=0.5*(3.*zeta(i,j-1)-zeta(i,j-2))
	    zz2=0.5*(3.*zeta(i+1,j-1)-zeta(i+1,j-2))
            zeta_c(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.9)then                                           

            !zeta_c(i,j)=zeta(i-1,j-1)  !change?

	    zeta_c(i,j)=0.5*
     *                  (0.5*(3.*zeta(i-1,j)-zeta(i-2,j))
     *                  +0.5*(3.*zeta(i-1,j-1)-zeta(i-2,j-1)))

          endif

        enddo
        enddo

! --- zeta_u

        do j=1,n-1
        do i=1,m

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=Ntype1(i,j)
	  endif	 

          if(Np.ne.4.and.Np.ne.5.and.Np.ne.8
     *      .and.Np.ne.96.and.Np.ne.92.and.Np.ne.86.and.Np.ne.82)then
            zeta_u(i,j)=0.5*(zeta(i-1,j)+zeta(i,j))
          endif

          if(Np.eq.4.or.Np.eq.5.or.Np.eq.8)then
            zeta_u(i,j)=0.5*(3.*zeta(i,j)-zeta(i+1,j))
          endif

          if(Np.eq.96.or.Np.eq.92.or.Np.eq.9
     *     .or.Np.eq.86.or.Np.eq.82)then
            zeta_u(i,j)=0.5*(3.*zeta(i-1,j)-zeta(i-2,j))
          endif

        enddo
        enddo

! --- zeta_v

        do j=1,n
        do i=1,m-1

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=Ntype1(i,j)
	  endif	  

          if(Np.ne.1.and.Np.ne.5.and.Np.ne.6.and.Np.ne.86
     *       .and.Np.ne.99.and.Np.ne.96
     *       .and.Np.ne.93.and.Np.ne.98.and.Np.ne.83.and.Np.ne.88)then
            zeta_v(i,j)=0.5*(zeta(i,j-1)+zeta(i,j))
          endif

          if(Np.eq.1.or.Np.eq.5.or.Np.eq.6)then
            zeta_v(i,j)=0.5*(3.*zeta(i,j)-zeta(i,j+1))
          endif

          if(Np.eq.93.or.Np.eq.98.or.Np.eq.9
     *       .or.Np.eq.83.or.Np.eq.88)then
            zeta_v(i,j)=0.5*(3.*zeta(i,j-1)-zeta(i,j-2))
          endif

        enddo
        enddo

	return

	end 

! ---------------------------------------------------------------------------

	subroutine calculate_Q1_zvc
!	use sadi
        use sadi
        integer Np,Np1
! --- second order averaging

! --- Q1_z

        do j=1,n-1
        do i=1,m-1
          Q1_z(i,j)=0.5*(Q1(i,j)+Q1(i+1,j))
        enddo
        enddo

! --- Q1_c

        do j=1,n
        do i=1,m

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	  
         
          if(Np.eq.0.or.Np.eq.2.or.Np.eq.7.or.Np.eq.92
     *      .or.Np.eq.3.or.Np.eq.8.or.Np.eq.4.or.Np.eq.9)then   

            Q1_c(i,j)=0.5*(Q1(i,j)+Q1(i,j-1))

          elseif(Np.eq.5.or.Np.eq.1.or.Np.eq.6.or.Np.eq.96)then

            Q1_c(i,j)=0.5*(3.*Q1(i,j)-Q1(i,j+1))

          elseif(Np.eq.98.or.Np.eq.93.or.Np.eq.97)then

            Q1_c(i,j)=0.5*(3.*Q1(i,j-1)-Q1(i,j-2))

          endif

        enddo
        enddo

! --- Q1_v the second order (please see previous code for fourth-order)

        do j=1,n
        do i=1,m-1

 	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	          

          if(Np.eq.0.or.Np.eq.4.or.Np.eq.8.or
     *      .Np.eq.3.or.Np.eq.7.or.Np.eq.2) then

             Q1_v(i,j)=0.25*(Q1(i,j)+Q1(i+1,j)+Q1(i+1,j-1)+Q1(i,j-1))

          elseif(Np.eq.5.or.Np.eq.1.or.Np.eq.6)then

             Q1_v(i,j)=0.5*(0.5*(3.*Q1(i,j)-Q1(i,j+1))
     *                    +0.5*(3.*Q1(i+1,j)-Q1(i+1,j+1)))

          elseif(Np.eq.98.or.Np.eq.93.or.Np.eq.9)then

             Q1_v(i,j)=0.5*(0.5*(3.*Q1(i,j-1)-Q1(i,j-2))
     *                    +0.5*(3.*Q1(i+1,j-1)-Q1(i+1,j-2)))

          endif

        enddo
        enddo

	return

	end

! ---------------------------------------------------------------------------

	subroutine calculate_Q2_zuc
!	use sadi
        use sadi
        integer Np,Np1
! --- second order averaging

! --- Q2_z

        do j=1,n-1
        do i=1,m-1
          Q2_z(i,j)=0.5*(Q2(i,j)+Q2(i,j+1))
        enddo
        enddo

! --- Q2_c

        do j=1,n
        do i=1,m

 	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	            

          if(Np.eq.0.or.Np.eq.3.or.Np.eq.7.or.Np.eq.93
     *      .or.Np.eq.2.or.Np.eq.6.or.Np.eq.1)then ! change 19

            Q2_c(i,j)=0.5*(Q2(i,j)+Q2(i-1,j))

	  elseif(Np.eq.9)then

	    Q2_c(i,j)=0.

          elseif(Np.eq.5.or.Np.eq.4.or.Np.eq.8.or.Np.eq.98)then

            Q2_c(i,j)=0.5*(3.*Q2(i,j)-Q2(i+1,j))

          elseif(Np.eq.96.or.Np.eq.92.or.Np.eq.97)then

            Q2_c(i,j)=0.5*(3.*Q2(i-1,j)-Q2(i-2,j))

          endif

        enddo
        enddo

! --- Q2_u the second order (See old code for fourth-order)

        do j=1,n-1
        do i=1,m

 	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	                      

          if(Np.eq.0.or.Np.eq.3.or.Np.eq.7.or
     *      .Np.eq.1.or.Np.eq.6.or.Np.eq.2) then

            Q2_u(i,j)=0.25*(Q2(i,j)+Q2(i,j+1)+Q2(i-1,j+1)+Q2(i-1,j))

          elseif(Np.eq.4.or.Np.eq.5.or.Np.eq.8)then

            Q2_u(i,j)=0.5*(0.5*(3.*Q2(i,j)-Q2(i+1,j))
     *                    +0.5*(3.*Q2(i,j+1)-Q2(i+1,j+1)))

          elseif(Np.eq.92.or.Np.eq.96.or.Np.eq.9)then

            Q2_u(i,j)=0.5*(0.5*(3.*Q2(i-1,j)-Q2(i-2,j))
     *                    +0.5*(3.*Q2(i-1,j+1)-Q2(i-2,j+1)))

         endif

        enddo
        enddo

	return

	end 


! ----------------------------------------------------
! ---------------------------------------------------------------

	subroutine diffusion
!	use sadi
        use sadi
        integer Np,Np1
        real htot
	
c	call eddy_viscosity

	do j=1,n
	do i=1,m
	  anu_t(i,j)=anu_tb_const
	enddo
	enddo

	call subgrid_mixing


! ---  calculate anu at z

	do j=1,n-1
	do i=1,m-1
	  anu_z(i,j)=0.25*(anu(i,j)+anu(i+1,j)+anu(i+1,j+1)+anu(i,j+1))
	enddo
	enddo

! ---  calculate Dup11 at zeta, Dup12 has been done in subgrid_mixing,
!      Dup11 Dup12 and Dup22 at u may be replaced by Dupxx at c and z for
!      saving time

c ---  Dup11 at z point

	do j=1,n
	do i=1,m
	  htot=(depth_u(i,j)+zeta_u(i,j))
	  if(htot.lt.hi)htot=hi
          tr1(i,j)=Q1(i,j)/htot
	enddo
	enddo

        call der_xi_at_z

	do j=1,n
	do i=1,m
           Dup11(i,j)=g_22_z(i,j)/g_0_z(i,j)*trd1(i,j)
	enddo
	enddo

! ---  Dup11 _xi

	do j=1,n-1
	do i=1,m-1
	  htot=(depth(i,j)+zeta(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=sq_g_0_z(i,j)*anu_z(i,j)*htot*Dup11(i,j)
	enddo
	enddo

        call der_xi_at_u_symmetry  

	do j=1,n-1
	do i=1,m
        diffx(i,j)=trd1(i,j)
	enddo
	enddo

! ---  Dup12 _eta

	do j=1,n
	do i=1,m
	  htot=(depth_c(i,j)+zeta_c(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=sq_g_0(i,j)*anu(i,j)*htot*Dup12(i,j)
	enddo
	enddo

        call der_eta_at_u   

	do j=1,n-1
	do i=1,m
        diffx(i,j)=diffx(i,j)+trd1(i,j)
	enddo
	enddo

! --- calculate div T^2

! ---   Dup22 at z, Dup12 done in subgrid_mixing, Dup11_v Dup12_v and Dup22_v
!        can be replaced by Dup11_c, Dup12_c and Dup22_z for saving time

	do j=1,n
	do i=1,m-1
	  htot=(depth_v(i,j)+zeta_v(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=Q2(i,j)/htot
	enddo
	enddo

        call der_eta_at_z

	do j=1,n-1
	do i=1,m-1
        Dup22(i,j)=g_11_z(i,j)/g_0_z(i,j)*trd1(i,j)
	enddo
	enddo

c ---  Dup12 _xi

	do j=1,n
	do i=1,m
	  htot=(depth_c(i,j)+zeta_c(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=sq_g_0(i,j)*anu(i,j)*htot*Dup12(i,j)
	enddo
	enddo

        call der_xi_at_v  

	do j=1,n
	do i=1,m-1
          diffy(i,j)=trd1(i,j)
	enddo
	enddo

c ---   Dup22 _eta

	do j=1,n-1
	do i=1,m-1
	  htot=(depth(i,j)+zeta(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=sq_g_0_z(i,j)*anu_z(i,j)*htot*Dup22(i,j)
	enddo
	enddo

        call der_eta_at_v_symmetry  

	do j=1,n
	do i=1,m-1
          diffy(i,j)=diffy(i,j)+trd1(i,j)
	enddo
	enddo

	return

	end

c -------------------------------------------------

	subroutine subgrid_mixing
!	use sadi
        use sadi
        integer Np,Np1
        real htot

c ---  Dup11 at c point
	
	do j=1,n
	do i=1,m-1
	htot=(depth_v(i,j)+zeta_v(i,j))
	if(htot.lt.hi)htot=hi
        tr1(i,j)=Q1_v(i,j)/htot
	enddo
	enddo

        call der_xi_at_c_anti_symmetry

	do j=1,n
	do i=1,m
         Dup11(i,j)=g_22(i,j)/g_0(i,j)*trd1(i,j)
         Ddn11(i,j)=g_11(i,j)*trd1(i,j)
	enddo
	enddo

c --- Dup12 at c point

	do j=1,n-1
	do i=1,m
	htot=(depth_u(i,j)+zeta_u(i,j))
	if(htot.lt.hi)htot=hi
        tr1(i,j)=Q1(i,j)/htot
	enddo
	enddo

        call der_eta_at_c_symmetry

	do j=1,n
	do i=1,m
        Dup12(i,j)=0.5*g_11(i,j)/g_0(i,j)*trd1(i,j)
        Ddn12(i,j)=0.5*g_11(i,j)*trd1(i,j)
	enddo
	enddo

	do j=1,n
	do i=1,m-1
	htot=(depth_v(i,j)+zeta_v(i,j))
	if(htot.lt.hi)htot=hi
        tr1(i,j)=Q2(i,j)/htot
	enddo
	enddo

        call der_xi_at_c_symmetry

	do j=1,n
	do i=1,m
        Dup12(i,j)=Dup12(i,j)+0.5*g_22(i,j)/g_0(i,j)*trd1(i,j)
        Ddn12(i,j)=Ddn12(i,j)+0.5*g_22(i,j)*trd1(i,j)
	enddo
	enddo

c ---  Dup22 at c

	do j=1,n-1
	do i=1,m
	htot=(depth_u(i,j)+zeta_u(i,j))
	if(htot.lt.hi)htot=hi
        tr1(i,j)=Q2_u(i,j)/htot
	enddo
	enddo

        call der_eta_at_c_anti_symmetry

	do j=1,n
	do i=1,m
        Dup22(i,j)=g_11(i,j)/g_0(i,j)*trd1(i,j)
        Ddn22(i,j)=g_22(i,j)*trd1(i,j)

        enddo
	enddo

! --- deal with various boundaries
        do j=1,n
        do i=1,m
          Np=Ntype(i,j)
          if(Np.lt.80)then 
             Np=mod(Np,10)
          else
             Np=90+mod(Np,10)
          endif
          if(Np.eq.1)then
            Dup11(i,j)=Dup11(i,j+1)
            Dup12(i,j)=Dup12(i,j+1)
            Dup22(i,j)=Dup22(i,j+1)
            Ddn11(i,j)=Ddn11(i,j+1)
            Ddn12(i,j)=Ddn12(i,j+1)
            Ddn22(i,j)=Ddn22(i,j+1)
          endif
          if(Np.eq.3)then
            Dup11(i,j)=Dup11(i,j-1)
            Dup12(i,j)=Dup12(i,j-1)
            Dup22(i,j)=Dup22(i,j-1)
            Ddn11(i,j)=Ddn11(i,j-1)
            Ddn12(i,j)=Ddn12(i,j-1)
            Ddn22(i,j)=Ddn22(i,j-1)
          endif
          if(Np.eq.4)then
            Dup11(i,j)=Dup11(i+1,j)
            Dup12(i,j)=Dup12(i+1,j)
            Dup22(i,j)=Dup22(i+1,j)
            Ddn11(i,j)=Ddn11(i+1,j)
            Ddn12(i,j)=Ddn12(i+1,j)
            Ddn22(i,j)=Ddn22(i+1,j)
          endif
          if(Np.eq.2)then
            Dup11(i,j)=Dup11(i-1,j)
            Dup12(i,j)=Dup12(i-1,j)
            Dup22(i,j)=Dup22(i-1,j)
            Ddn11(i,j)=Ddn11(i-1,j)
            Ddn12(i,j)=Ddn12(i-1,j)
            Ddn22(i,j)=Ddn22(i-1,j)
          endif
          if(Np.eq.5)then
            Dup11(i,j)=Dup11(i+1,j+1)
            Dup12(i,j)=Dup12(i+1,j+1)
            Dup22(i,j)=Dup22(i+1,j+1)
            Ddn11(i,j)=Ddn11(i+1,j+1)
            Ddn12(i,j)=Ddn12(i+1,j+1)
            Ddn22(i,j)=Ddn22(i+1,j+1)
          endif          
          if(Np.eq.7)then
            Dup11(i,j)=Dup11(i-1,j-1)
            Dup12(i,j)=Dup12(i-1,j-1)
            Dup22(i,j)=Dup22(i-1,j-1)
            Ddn11(i,j)=Ddn11(i-1,j-1)
            Ddn12(i,j)=Ddn12(i-1,j-1)
            Ddn22(i,j)=Ddn22(i-1,j-1)
          endif
          if(Np.eq.6)then
            Dup11(i,j)=Dup11(i-1,j+1)
            Dup12(i,j)=Dup12(i-1,j+1)
            Dup22(i,j)=Dup22(i-1,j+1)
            Ddn11(i,j)=Ddn11(i-1,j+1)
            Ddn12(i,j)=Ddn12(i-1,j+1)
            Ddn22(i,j)=Ddn22(i-1,j+1)
          endif          
          if(Np.eq.8)then
            Dup11(i,j)=Dup11(i+1,j-1)
            Dup12(i,j)=Dup12(i+1,j-1)
            Dup22(i,j)=Dup22(i+1,j-1)
            Ddn11(i,j)=Ddn11(i+1,j-1)
            Ddn12(i,j)=Ddn12(i+1,j-1)
            Ddn22(i,j)=Ddn22(i+1,j-1)
          endif
          if(Np.eq.93)then
            Dup11(i,j)=Dup11(i,j-2)
            Dup12(i,j)=Dup12(i,j-2)
            Dup22(i,j)=Dup22(i,j-2)
            Ddn11(i,j)=Ddn11(i,j-2)
            Ddn12(i,j)=Ddn12(i,j-2)
            Ddn22(i,j)=Ddn22(i,j-2)
          endif
          if(Np.eq.98)then
            Dup11(i,j)=Dup11(i+1,j-2)
            Dup12(i,j)=Dup12(i+1,j-2)
            Dup22(i,j)=Dup22(i+1,j-2)
            Ddn11(i,j)=Ddn11(i+1,j-2)
            Ddn12(i,j)=Ddn12(i+1,j-2)
            Ddn22(i,j)=Ddn22(i+1,j-2)
          endif
          if(Np.eq.97)then
            Dup11(i,j)=Dup11(i-1,j-2)
            Dup12(i,j)=Dup12(i-1,j-2)
            Dup22(i,j)=Dup22(i-1,j-2)
            Ddn11(i,j)=Ddn11(i-1,j-2)
            Ddn12(i,j)=Ddn12(i-1,j-2)
            Ddn22(i,j)=Ddn22(i-1,j-2)
          endif
          if(Np.eq.92)then
            Dup11(i,j)=Dup11(i-2,j)
            Dup12(i,j)=Dup12(i-2,j)
            Dup22(i,j)=Dup22(i-2,j)
            Ddn11(i,j)=Ddn11(i-2,j)
            Ddn12(i,j)=Ddn12(i-2,j)
            Ddn22(i,j)=Ddn22(i-2,j)
          endif
          if(Np.eq.96)then
            Dup11(i,j)=Dup11(i-2,j+1)
            Dup12(i,j)=Dup12(i-2,j+1)
            Dup22(i,j)=Dup22(i-2,j+1)
            Ddn11(i,j)=Ddn11(i-2,j+1)
            Ddn12(i,j)=Ddn12(i-2,j+1)
            Ddn22(i,j)=Ddn22(i-2,j+1)
          endif
        enddo
        enddo

! ---   anu_s()

	do j=1,n
	do i=1,m
        anu_s(i,j)=c_subgrid*(sq_g_0(i,j))*
     *              SQRT(ABS(Dup11(i,j)*Ddn11(i,j)
     *                  +2.*Dup12(i,j)*Ddn12(i,j)
     *                  +Dup22(i,j)*Ddn22(i,j)))
         anu(i,j)=anu_t(i,j)+anu_s(i,j)+anu_tb(i,j)
	enddo
	enddo

	return

	end

! --------------------------------------------------------

        subroutine der_xi_at_c_anti_symmetry
!        use sadi
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m

	if(Ntype1(i,j).lt.80)then
        Np=mod(Ntype1(i,j),10)
	else
	Np=Ntype1(i,j)
	endif

	if(Ntype(i,j).lt.80.and.Ntype(i,j).ge.20)then
        Np1=mod(Ntype(i,j),10)
	else
	Np1=Ntype(i,j)
	endif

        if(Np.ne.4.and.Np.ne.5.and.Np.ne.98.and.Np.ne.8
     *     .and.Np.ne.92.and.Np.ne.96.and.Np.ne.97
     *     .and.Np.ne.88.and.Np.ne.82.and.Np.ne.86.and.Np.ne.87)then

          trd1(i,j)=tr1(i,j)-tr1(i-1,j)

          elseif(Np.eq.4.or.Np.eq.5.or.Np.eq.98.or.Np.eq.8)then

	    if(Np1.eq.4.or.Np1.eq.5.or.Np1.eq.8)then
	      trd1(i,j)=0.
	     else
              trd1(i,j)=tr1(i,j)+tr1(i,j)
	     endif

	  elseif(Np.eq.88)then
	    trd1(i,j)=0.

          elseif(Np.eq.92.or.Np.eq.96.or.Np.eq.97)then

            trd1(i,j)=-tr1(i-1,j)-tr1(i-1,j)

          elseif(Np.eq.82.or.Np.eq.86.or.Np.eq.87)then

            trd1(i,j)=0.

        endif

        enddo
        enddo
        
	return

        end

! ------------------------------------------------------------

        subroutine der_eta_at_c_anti_symmetry
!        use sadi
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m

	if(Ntype1(i,j).lt.80)then
	Np=mod(Ntype1(i,j),10)
	else
        Np=Ntype1(i,j)
	endif

	if(Ntype(i,j).lt.80.and.Ntype(i,j).ge.20)then
        Np1=mod(Ntype(i,j),10)
	else
	Np1=Ntype(i,j)
	endif

        if(Np.ne.1.and.Np.ne.5.and.Np.ne.96.and.Np.ne.6
     *     .and.Np.ne.93.and.Np.ne.98.and.Np.ne.97.and.Np.ne.99
     *     .and.Np.ne.86.and.Np.ne.83.and.Np.ne.88.and.Np.ne.87)then

          trd1(i,j)=tr1(i,j)-tr1(i,j-1)

          elseif(Np.eq.1.or.Np.eq.5.or.Np.eq.96.or.Np.eq.6)then

	    if(Np1.eq.1.or.Np1.eq.5.or.Np1.eq.6)then
	     trd1(i,j)=0.
	    else
             trd1(i,j)=tr1(i,j)+tr1(i,j)
	    endif

	  elseif(Np.eq.88)then
	    trd1(i,j)=0.

          elseif(Np.eq.93.or.Np.eq.98.or.Np.eq.97)then
              trd1(i,j)=-tr1(i,j-1)-tr1(i,j-1)

          elseif(Np.eq.83.or.Np.eq.88.or.Np.eq.87)then
              trd1(i,j)=0.

        endif

        enddo
        enddo
        
	return

        end

! ------------------------------------------------------------

        subroutine der_xi_at_z
!        use sadi
         use sadi

        do j=1,n-1
        do i=1,m-1

          trd1(i,j)=tr1(i+1,j)-tr1(i,j)

        enddo
        enddo

	return
        
        end

! ------------------------------------------------------------

        subroutine der_eta_at_z
!        use sadi
        use sadi

        do j=1,n-1
        do i=1,m-1

          trd1(i,j)=tr1(i,j+1)-tr1(i,j)

        enddo
        enddo

	return

        end

! -----------------------------------------------------------

        subroutine der_eta_at_u
!        use sadi
        use sadi
 
        do j=1,n-1
        do i=1,m

          trd1(i,j)=tr1(i,j+1)-tr1(i,j)

        enddo
        enddo

	return

        end

! -----------------------------------------------------------

        subroutine der_xi_at_v
        use sadi

        do j=1,n
        do i=1,m-1

          trd1(i,j)=tr1(i+1,j)-tr1(i,j)

        enddo
        enddo

	return

        end

c -----------------------------------------------------

        subroutine der_xi_at_u_symmetry
        use sadi
        integer Np,Np1

        do j=1,n-1
        do i=1,m

	if(Ntype1(i,j).lt.80)then
	Np=mod(Ntype1(i,j),10)
	else
        Np=90+mod(Ntype1(i,j),10)
	endif	

        if(Np.eq.0.or.Np.eq.1.or.Np.eq.3.or.Np.eq.6.or.Np.eq.7
     *     .or.Np.eq.2)then

          trd1(i,j)=tr1(i,j)-tr1(i-1,j)

          elseif(Np.eq.92.or.Np.eq.96.or.Np.eq.9)then

            trd1(i,j)=0.

          elseif(Np.eq.4.or.Np.eq.5.or.Np.eq.8)then

            trd1(i,j)=0.

        endif

        enddo
        enddo

	return
        end

! -----------------------------------------------------------------

        subroutine der_xi_at_u_anti_symmetry
        use sadi
        integer Np,Np1

        do j=1,n-1
        do i=1,m

	if(Ntype1(i,j).lt.80)then
	Np=mod(Ntype1(i,j),10)
	else
        Np=Ntype1(i,j)
	endif	 

	if(Ntype(i,j).lt.80.and.Ntype(i,j).ge.20)then
        Np1=mod(Ntype(i,j),10)
	else
	Np1=Ntype(i,j)
	endif

        if(Np.eq.0.or.Np.eq.1.or.Np.eq.3.or.Np.eq.6.or.Np.eq.7
     *     .or.Np.eq.2)then

          trd1(i,j)=tr1(i,j)-tr1(i-1,j)

          elseif(Np.eq.92.or.Np.eq.96.or.Np.eq.9)then

            trd1(i,j)=-2.*tr1(i-1,j)

          elseif(Np.eq.82.or.Np.eq.86)then
	
	    trd1(i,j)=0.

          elseif(Np.eq.4.or.Np.eq.5.or.Np.eq.8)then

	    if(Np1.eq.4.or.Np1.eq.5.or.Np1.eq.8)then
	      trd1(i,j)=0.
	    else
              trd1(i,j)=2.*tr1(i,j)
	    endif

        endif

        enddo
        enddo

	return

        end

! -----------------------------------------------------------------

        subroutine der_eta_at_u_symmetry
        use sadi
        integer Np,Np1

        do j=1,n-1
        do i=1,m

          trd1(i,j)=tr1(i,j+1)-tr1(i,j)

        enddo
        enddo

	return
        end
! -----------------------------------------------------------------

        subroutine der_eta_at_u_anti_symmetry
        use sadi
        integer Np,Np1

        do j=1,n-1
        do i=1,m

          trd1(i,j)=tr1(i,j+1)-tr1(i,j)

        enddo
        enddo

	return
        end

! -----------------------------------------------------------------

        subroutine der_eta_at_v_symmetry
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m-1

	if(Ntype1(i,j).lt.80)then
	Np=mod(Ntype1(i,j),10)
	else
        Np=90+mod(Ntype1(i,j),10)
	endif	         

        if(Np.eq.0.or.Np.eq.4.or.Np.eq.2.or.Np.eq.8.or.Np.eq.7
     *     .or.Np.eq.3)then

          trd1(i,j)=tr1(i,j)-tr1(i,j-1)

          elseif(Np.eq.93.or.Np.eq.98.or.Np.eq.9)then

            trd1(i,j)=0.

          elseif(Np.eq.1.or.Np.eq.5.or.Np.eq.6)then

            trd1(i,j)=0.

        endif

        enddo
        enddo

	return
        end 
! -----------------------------------------------------------------

        subroutine der_eta_at_v_anti_symmetry
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m-1

	if(Ntype1(i,j).lt.80)then
	Np=mod(Ntype1(i,j),10)
	else
        Np=Ntype1(i,j)
	endif	

	if(Ntype(i,j).lt.80.and.Ntype(i,j).ge.20)then
        Np1=mod(Ntype(i,j),10)
	else
	Np1=Ntype(i,j)
	endif
  

        if(Np.eq.0.or.Np.eq.4.or.Np.eq.2.or.Np.eq.8.or.Np.eq.7
     *     .or.Np.eq.3)then

          trd1(i,j)=tr1(i,j)-tr1(i,j-1)

          elseif(Np.eq.93.or.Np.eq.98.or.Np.eq.9)then

            trd1(i,j)=-2.*tr1(i,j-1)

          elseif(Np.eq.83.or.Np.eq.88)then

            trd1(i,j)=0.

          elseif(Np.eq.1.or.Np.eq.5.or.Np.eq.6)then

	    if(Np1.eq.1.or.Np1.eq.5.or.Np1.eq.6)then
            trd1(i,j)=0.
	    else
            trd1(i,j)=2.*tr1(i,j)
	    endif

        endif

        enddo
        enddo

	return
        end

! -----------------------------------------------------------------

        subroutine der_xi_at_v_symmetry
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m-1

          trd1(i,j)=tr1(i+1,j)-tr1(i,j)

        enddo
        enddo

	return
        end

! -----------------------------------------------------------------

        subroutine der_xi_at_v_anti_symmetry
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m-1

          trd1(i,j)=tr1(i+1,j)-tr1(i,j)

        enddo
        enddo

	return
        end

! -------------------------------------------------------------------

        subroutine gridtype
        use sadi
        integer Np,Np1

	do j=1,n
	do i=1,m
	  Ntype(i,j)=99
	enddo
	enddo


c --- box
       
	  do j=2,n-1
	   do i=2,m-1
	   Ntype(i,j)=0
           enddo
	  enddo

        do i=1,m-1
          Ntype(i,1)=11
          Ntype(i,n-1)=13
          Ntype(i,n)=93
        enddo

        do j=1,n-1
          Ntype(1,j)=14
          Ntype(m-1,j)=12
          Ntype(m,j)=92
        enddo

        Ntype(1,1)=15
        Ntype(m-1,1)=12
        Ntype(m-1,n-1)=12
        Ntype(1,n-1)=18     

	Ntype(m-1,1)=16
	Ntype(m-1,n-1)=17
        Ntype(m,1)=96
        Ntype(m,n)=97
        Ntype(1,n)=98	

! --- east_elevation boundary ------------------

       if(EAST_ele2_flx3_grd4_rad5_per6.ge.2)then
       Np=EAST_ele2_flx3_grd4_rad5_per6*10
! --- middle
        if(jstart_E.ge.2.and.jend_E.le.n-2)then

        Ntype(m-2,1)=16
        Ntype(m-1,1)=96 
        Ntype(m,1)=99

        do j=2,jstart_E-1
          Ntype(m-2,j)=12
          Ntype(m-1,j)=92
          Ntype(m,j)=99
        enddo
        Ntype(m-1,jstart_E)=Np+6
        Ntype(m,jstart_E)=86
        do j=jstart_E+1,jend_E-1
          Ntype(m-1,j)=Np+2
          Ntype(m,j)=82
        enddo
        Ntype(m-1,jend_E)=Np+7
        Ntype(m,jend_E)=82
        Ntype(m-2,jend_E+1)=12
        Ntype(m-1,jend_E+1)=19
        Ntype(m,jend_E+1)=87
        do j=jend_E+2,n-1
          Ntype(m-2,j)=12
          Ntype(m-1,j)=92
          Ntype(m,j)=99
        enddo	
        Ntype(m-2,n-1)=17
        Ntype(m-2,n)=93
        Ntype(m-1,n)=97
        Ntype(m,n)=99

        endif

! --- bottom - middle
        if(jstart_E.le.1.and.jend_E.le.n-2)then
        do j=2,jend_E-1
          Ntype(m-1,j)=Np+2
          Ntype(m,j)=82
        enddo
	Ntype(m-1,1)=Np+6
        Ntype(m,1)=86
        Ntype(m-1,jend_E)=Np+7
        Ntype(m,jend_E)=82
        Ntype(m-2,jend_E+1)=12
        Ntype(m-1,jend_E+1)=19
        Ntype(m,jend_E+1)=87
        do j=jend_E+2,n-1
          Ntype(m-2,j)=12
          Ntype(m-1,j)=92
          Ntype(m,j)=99
        enddo	
        Ntype(m-2,n-1)=17
        Ntype(m-2,n)=93
        Ntype(m-1,n)=97
        Ntype(m,n)=99
	
        endif

! --- middle - top
        if(jstart_E.ge.2.and.jend_E.ge.n-1)then 
        Ntype(m-2,1)=16
        Ntype(m-1,1)=96 
        Ntype(m,1)=99

        do j=2,jstart_E-1
          Ntype(m-2,j)=12
          Ntype(m-1,j)=92
          Ntype(m,j)=99
        enddo
        Ntype(m-1,jstart_E)=Np+6
        Ntype(m,jstart_E)=86
        do j=jstart_E+1,n-1
          Ntype(m-1,j)=Np+2
          Ntype(m,j)=82
        enddo
	Ntype(m-1,n-1)=Np+7
        Ntype(m,n)=87

        endif

! --- bottom - top
        if(jstart_E.le.1.and.jend_E.ge.n-1)then 
        do j=2,n-1
          Ntype(m-1,j)=Np+2
          Ntype(m,j)=82
        enddo 
	Ntype(m-1,1)=Np+6
        Ntype(m,1)=86     
	Ntype(m-1,n-1)=Np+7
        Ntype(m,n)=87

        endif

        endif

! ---------- west elevation --------------

        if(west_ele2_flx3_grd4_rad5_per6.ge.2)then
           Np=west_ele2_flx3_grd4_rad5_per6*10

! --- middle

         if(jstart_W.ge.2.and.jend_W.le.n-2)then

           Ntype(1,1)=99
           Ntype(2,1)=15
           do j=2,jstart_W-1
             Ntype(1,j)=99
             Ntype(2,j)=14
           enddo
           Ntype(1,jstart_W)=Np+5
           do j=jstart_W+1,jend_W-1
             Ntype(1,j)=Np+4
           enddo
           Ntype(1,jend_W)=Np+8
           do j=jend_W+1,n-1
             Ntype(1,j)=99
             Ntype(2,j)=14
           enddo
           Ntype(2,n-1)=18
           Ntype(1,n)=99      
           Ntype(2,n)=98

           endif
! --- bottom - top
         if(jstart_W.le.1.and.jend_W.ge.n-1)then

           Ntype(1,1)=Np+5
           do j=2,n-2
             Ntype(1,j)=Np+4
           enddo
           Ntype(1,n-1)=Np+8
           Ntype(1,n)=88

         endif
! --- bottom - middle
         if(jstart_W.le.1.and.jend_W.le.n-2)then

           Ntype(1,1)=Np+5
           do j=2,jend_W-1
             Ntype(1,j)=Np+4
           enddo
           Ntype(1,jend_W)=Np+8
           do j=jend_W+1,n-1
             Ntype(1,j)=99
             Ntype(2,j)=14
           enddo
           Ntype(2,n-1)=18
           Ntype(1,n)=99      
           Ntype(2,n)=98

         endif

! --- middle - top
         if(jstart_W.ge.2.and.jend_W.ge.n-1)then
           Ntype(1,1)=99
           Ntype(2,1)=15
           do j=2,jstart_W-1
             Ntype(1,j)=99
             Ntype(2,j)=14
           enddo
           Ntype(1,jstart_W)=Np+5
           do j=jstart_W+1,n-2
             Ntype(1,j)=Np+4
           enddo
           Ntype(1,n-1)=Np+8
           Ntype(1,n)=88

          endif

          endif

! --- south elevation boundary -----------

        if(south_ele2_flx3_grd4_rad5_per6.ge.2)then
          Np=south_ele2_flx3_grd4_rad5_per6*10
! --- middle
          if(istart_S.ge.2.and.iend_S.le.m-2)then

          Ntype(1,1)=99
          Ntype(1,2)=15
          do i=2,istart_S-1
            Ntype(i,1)=99
            Ntype(i,2)=11
          enddo
          Ntype(istart_S,1)=Np+5
          do i=istart_S+1,iend_S-1
            Ntype(i,1)=Np+1
          enddo
          Ntype(iend_S,1)=Np+6
          Ntype(iend_S+1,1)=86
          Ntype(iend_S+1,2)=11
          do i=iend_S+2,m-2
            Ntype(i,1)=99
            Ntype(i,2)=11
          enddo
          Ntype(m-1,1)=99
          Ntype(m-1,2)=16
          Ntype(m,1)=99
          Ntype(m,2)=96
          endif

! ---  left and right
          if(istart_S.le.1.and.iend_S.ge.m-1)then
            Ntype(1,1)=Np+5
            do i=2,m-2
              Ntype(i,1)=Np+1
            enddo
            Ntype(m-1,1)=Np+6
            Ntype(m,1)=86
          endif

! --- left - middle
          if(istart_S.le.1.and.iend_S.le.m-2)then
            Ntype(1,1)=Np+5
            do i=2,iend_S-1
              Ntype(i,1)=Np+1
            enddo
            Ntype(iend_S,1)=Np+6
            Ntype(iend_S+1,1)=86
            Ntype(iend_S+1,2)=11
            do i=iend_S+2,m-2
              Ntype(i,1)=99
              Ntype(i,2)=11
            enddo
            Ntype(m-1,1)=99
            Ntype(m-1,2)=16
            Ntype(m,1)=99
            Ntype(m,2)=96
           endif

! --- middle - right
          if(istart_S.ge.2.and.iend_S.ge.m-1)then

          Ntype(1,1)=99
          Ntype(1,2)=15
          do i=2,istart_S-1
            Ntype(i,1)=99
            Ntype(i,2)=11
          enddo
          Ntype(istart_S,1)=Np+5
          do i=istart_S+1,m-2
            Ntype(i,1)=Np+1
          enddo
          Ntype(m-1,1)=Np+6
          Ntype(m,1)=86
          endif

        endif

! --- north boundary -----------

        if(north_ele2_flx3_grd4_rad5_per6.ge.2)then
         Np=north_ele2_flx3_grd4_rad5_per6*10
! --- middle
          if(istart_N.ge.2.and.iend_N.le.m-2)then
          Ntype(1,n)=99
          Ntype(1,n-1)=98
          Ntype(1,n-2)=18
          do i=2,istart_N-1
            Ntype(i,n)=99
            Ntype(i,n-1)=93
            Ntype(i,n-2)=13      
          enddo
          Ntype(istart_N,n)=88
          Ntype(istart_N,n-1)=Np+8
          do i=istart_N+1,iend_N-1
            Ntype(i,n)=83
            Ntype(i,n-1)=Np+3
          enddo
          Ntype(iend_N,n)=83
          Ntype(iend_N,n-1)=Np+7
          Ntype(iend_N+1,n)=87
          Ntype(iend_N+1,n-1)=19
          Ntype(iend_N+1,n-2)=13
          do i=iend_N+2,m-2
            Ntype(i,n)=99
            Ntype(i,n-1)=93
            Ntype(i,n-2)=13
          enddo  
          Ntype(m-1,n)=99
          Ntype(m-1,n-1)=93
          Ntype(m-1,n-2)=17
          Ntype(m,n)=99
          Ntype(m,n-1)=97            
          endif

! --- left - right
          if(istart_N.le.1.and.iend_N.ge.m-1)then
          Ntype(1,n)=88
          Ntype(1,n-1)=Np+8
          do i=2,m-2
            Ntype(i,n)=83
            Ntype(i,n-1)=Np+3
          enddo
          Ntype(m-1,n)=83
          Ntype(m-1,n-1)=Np+7
          Ntype(m,n)=87
          Ntype(m,n-1)=82 
          endif    

! --- left - middle
          if(istart_N.le.1.and.iend_N.le.m-2)then
          Ntype(1,n)=88
          Ntype(1,n-1)=Np+8
          do i=2,iend_N-1
            Ntype(i,n)=83
            Ntype(i,n-1)=Np+3
          enddo
          Ntype(iend_N,n)=83
          Ntype(iend_N,n-1)=Np+7
          Ntype(iend_N+1,n)=87
          Ntype(iend_N+1,n-1)=19
          Ntype(iend_N+1,n-2)=13
          do i=iend_N+2,m-2
            Ntype(i,n)=99
            Ntype(i,n-1)=93
            Ntype(i,n-2)=13
          enddo  
          Ntype(m-1,n)=99
          Ntype(m-1,n-1)=93
          Ntype(m-1,n-2)=17
          Ntype(m,n)=99
          Ntype(m,n-1)=97            
          endif

! --- middle - right
          if(istart_N.ge.2.and.iend_N.ge.m-1)then
          Ntype(1,n)=99
          Ntype(1,n-1)=98
          Ntype(1,n-2)=18
          do i=2,istart_N-1
            Ntype(i,n)=99
            Ntype(i,n-1)=93
            Ntype(i,n-2)=13      
          enddo
          Ntype(istart_N,n)=88
          Ntype(istart_N,n-1)=Np+8
          do i=istart_N+1,m-2
            Ntype(i,n)=83
            Ntype(i,n-1)=Np+3
          enddo
          Ntype(m-1,n)=83
          Ntype(m-1,n-1)=Np+7
          Ntype(m,n)=87
          Ntype(m,n-1)=82 
          endif

        endif

          return

        end

        subroutine newgridtype
        use sadi

c --- seperate ntype according to boundary conditions

	do j=1,n
	do i=1,m
	  Ntype1(i,j)=Ntype(i,j)
	enddo
	enddo

	do j=1,n-1
	do i=1,m-1
	  if (Ntype(i,j).eq.24.and.Ntype(i+1,j).eq.11)Ntype1(i,j)=15
	  if (Ntype(i,j).eq.24.and.Ntype(i+1,j).eq.13)Ntype1(i,j)=18
	  if (Ntype(i,j).eq.21.and.Ntype(i,j+1).eq.14)Ntype1(i,j)=15
	  if (Ntype(i,j).eq.21.and.Ntype(i,j+1).eq.12)Ntype1(i,j)=16
          if (Ntype(i,j).eq.21.and.Ntype(i,j+1).eq.14)Ntype1(i,j)=15
	enddo
	enddo

	do j=1,n
	do i=2,m
	  if (Ntype(i,j).eq.22.and.Ntype(i-1,j).eq.11)Ntype1(i,j)=16
	  if (Ntype(i,j).eq.22.and.Ntype(i-1,j).eq.13)Ntype1(i,j)=17
	enddo
	enddo

	do j=2,n
	do i=1,m
	  if (Ntype(i,j).eq.23.and.Ntype(i,j-1).eq.14)Ntype1(i,j)=18
	  if (Ntype(i,j).eq.23.and.Ntype(i,j-1).eq.12)Ntype1(i,j)=17
	  if (Ntype(i,j).eq.23.and.Ntype(i,j-1).eq.14)Ntype1(i,j)=18
	enddo
	enddo

	open(1,file='type1.dat')
	do j=n,1,-1
	write(1,991)(Ntype1(i,j),i=1,m)
	enddo
991	format(2000I3)	
	close(1)

	open(1,file='type.dat')
	do j=n,1,-1
	write(1,991)(Ntype(i,j),i=1,m)
	enddo
	close(1)        

        return
        end

	subroutine adaptor_wave_to_circ
	use sadi
!	include 'pass.h'
        integer ny2
	real ttm(ntrim)
	ny2=n/2

! --- pass1

	do j=1,n
	do i=1,m
	S_xx(i,j) = Intp_Sxx_Circ(i,j)/1030.
	S_xy(i,j) = Intp_Sxy_Circ(i,j)/1030.
	S_yy(i,j) = Intp_Syy_Circ(i,j)/1030.
	
	Qwx(i,j) = Intp_MassFluxU_Circ(i,j)
	Qwy(i,j) = Intp_MassFluxV_Circ(i,j)

	dissipation(i,j) = Intp_Diss_Circ(i,j)

	ibrk(i,j)=Intp_ibrk_Circ(i,j)
	
	Angle(i,j) = Intp_Theta_Circ(i,j)*3.1415926/180.

	Height(i,j) = Intp_Height_Circ(i,j)

	u0(i,j) = Intp_ubott_Circ(i,j)

	enddo
	enddo


! --- suggested by Kevin because of refdif output
	do j=1,n
	do i=1,2
	  S_xx(i,j)=S_xx(3,j)
	  S_xy(i,j)=S_xy(3,j)
	  S_yy(i,j)=S_yy(3,j)
	  Qwx(i,j)=Qwx(3,j)
	  Qwy(i,j)=Qwy(3,j)
	  dissipation(i,j)=dissipation(3,j)
	  Angle(i,j)=Angle(3,j)
	  Height(i,j)=Height(3,j)
	  u0(i,j)=u0(3,j)
	enddo
	enddo

! --- for no-wave case, set u0 small value for calculation of bottom friction
	do j=1,n
	do i=1,m
	  if(u0(i,j).eq.0.)u0(i,j)=0.01
	enddo
	enddo

	call convert_Sxx_Sxy_Syy_Qw

        return

	end

! -------------------------------------------------------------------

	subroutine convert_Sxx_Sxy_Syy_Qw
	use sadi

        if(Wave_To_Circ_Radiation)then

	do j=1,n
        do i=1,m
	S_11(i,j)=1./g_0(i,j)*(S_xx(i,j)*Y_xi_2(i,j)*Y_xi_2(i,j)
     *             - 2.*S_xy(i,j)*Y_xi_2(i,j)*X_xi_2(i,j)
     *             + S_yy(i,j)*X_xi_2(i,j)*X_xi_2(i,j))

	S_12(i,j)=1./g_0(i,j)*(-S_xx(i,j)*Y_xi_1(i,j)*Y_xi_2(i,j) 
     *             + S_xy(i,j)*X_xi_1(i,j)*Y_xi_2(i,j)
     *              + S_xy(i,j)*Y_xi_1(i,j)*X_xi_2(i,j) - S_yy(i,j)
     *             *X_xi_1(i,j)*X_xi_2(i,j))

	S_22(i,j)=1./g_0(i,j)*(S_xx(i,j)*Y_xi_1(i,j)*Y_xi_1(i,j) 
     *              - 2.*S_xy(i,j)*X_xi_1(i,j)*Y_xi_1(i,j)
     *              + S_yy(i,j)*X_xi_1(i,j)*X_xi_1(i,j))
        enddo
        enddo

! ---  calculate S_11 and S_22 at z points

        do j=1,n
        do i=1,m
	tmp1(i,j)=S_11(i,j)
	tmp2(i,j)=S_22(i,j)
        enddo
        enddo

	do j=1,n-1
	do i=1,m-1
	S_11(i,j)=0.25*
     *        (tmp1(i,j)+tmp1(i+1,j)+tmp1(i+1,j+1)+tmp1(i,j+1))
	S_22(i,j)=0.25*
     *        (tmp2(i,j)+tmp2(i+1,j)+tmp2(i+1,j+1)+tmp2(i,j+1))
	enddo
	enddo

! --- boundaries
        do j=1,n-1
         S_11(m,j)=S_11(m-1,j)
         S_22(m,j)=S_22(m-1,j)
        enddo
        do i=1,m-1
         S_11(i,n)=S_11(i,n-1)
         S_22(i,n)=S_22(i,n-1)
        enddo

        else
c --- convert Fx and Fy

        do j=1,n
        do i=1,m
c -- set wave forcing to 0 at shoreline
c       if(hh(i,j).le.1.0.and.hh(i,j).ge.0)then
c         Intp_Fx_Circ(i,j)=Intp_Fx_Circ(i,j)*hh(i,j)
c         Intp_Fy_Circ(i,j)=Intp_Fy_Circ(i,j)*hh(i,j)
c       endif

c       if(hh(i,j).lt.0.0)then
c         Intp_Fx_Circ(i,j)=0.
c         Intp_Fy_Circ(i,j)=0.
c       endif

        wavefx(i,j)=-X_xi_2(i,j)*Intp_Fy_Circ(i,j)
     *           +Y_xi_2(i,j)*Intp_Fx_Circ(i,j)
        wavefy(i,j)= X_xi_1(i,j)*Intp_Fy_Circ(i,j)
     *           -Y_xi_1(i,j)*Intp_Fx_Circ(i,j)
        wavefx(i,j)=wavefx(i,j)*9.8
        wavefy(i,j)=wavefy(i,j)*9.8
c        RADX(i,j)=wavefx(i,j)
c        RADY(i,j)=wavefy(i,j)
        enddo
        enddo

        endif


        do j=1,n
        do i=1,m
	Qw1(i,j)=-X_xi_2(i,j)/sq_g_0(i,j)*Qwy(i,j) 
     *           +Y_xi_2(i,j)/sq_g_0(i,j)*Qwx(i,j)
	Qw2(i,j)= X_xi_1(i,j)/sq_g_0(i,j)*Qwy(i,j) 
     *           -Y_xi_1(i,j)/sq_g_0(i,j)*Qwx(i,j)
        enddo
        enddo


        return

	end


! ----------------------------------------------
        subroutine radiation_stress
        use sadi
 
        do j=1,n-1
        do i=1,m-1
          tr1(i,j)=sq_g_0_z(i,j)*S_11(i,j)
        enddo
        enddo

	call der_xi_at_u_symmetry

        do j=1,n-1
        do i=1,m
          wavefx(i,j)=-trd1(i,j)
        enddo
        enddo

        do j=1,n 
        do i=1,m
          tr1(i,j)=sq_g_0(i,j)*S_12(i,j)
        enddo
        enddo

	call der_eta_at_u_symmetry

        do j=1,n-1
        do i=1,m
          wavefx(i,j)=wavefx(i,j)-trd1(i,j)	
        enddo
        enddo

        do j=1,n-1
        do i=1,m-1
          tr1(i,j)=sq_g_0(i,j)*S_12(i,j)
        enddo
        enddo

	call der_xi_at_v_symmetry

        do j=1,n
        do i=1,m-1
          wavefy(i,j)=-trd1(i,j)
        enddo
        enddo

        do j=1,n-1
        do i=1,m-1
          tr1(i,j)=sq_g_0_z(i,j)*S_22(i,j)
        enddo
        enddo

	call der_eta_at_v_symmetry

        do j=1,n
        do i=1,m-1
          wavefy(i,j)=wavefy(i,j)-trd1(i,j)
        enddo
        enddo

        do j=1,n
          wavefy(m,j)=wavefy(m-1,j)
        enddo

        return

        end

! ----------------------------------------------------
	subroutine adaptor_circ_to_wave
	use sadi
	real ttm(mtrim),zz1,zz2,zw(mtrim,ntrim)
        integer ny2,Np,Np1,ny
	ny2=ny/2

! --- pass2

! --  take off dry points

        do j=1,n-1
        do i=1,m-1
          if(HH(i,j).lt.0)then
            zw(i,j)=0.
          else
            zw(i,j)=z(i,j)
          endif
        enddo
        enddo

! --  convert z to c points %%

         do j=1,n
         do i=1,m

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=Ntype1(i,j)
	  endif

          if(Np.eq.2.or.Np.eq.7.or.Np.eq.3)then

            zconv(i,j)=0.25*
     *             (zw(i,j)+zw(i-1,j)+zw(i-1,j-1)+zw(i,j-1))

	  elseif(Np.eq.0)then
	    Np1=Ntype1(i-1,j-1)
	    if(Np1.eq.99)then

            zconv(i,j)=1./3.*
     *             (zw(i,j)+zw(i-1,j)+zw(i,j-1))

	    else

            zconv(i,j)=0.25*
     *             (zw(i,j)+zw(i-1,j)+zw(i-1,j-1)+zw(i,j-1))

	    endif

          elseif(Np.eq.1.or.Np.eq.6)then

            zconv(i,j)=0.5*
     *                  (0.5*(3.*zw(i-1,j)-zw(i-1,j+1))
     *                  +0.5*(3.*zw(i,j)-zw(i,j+1)))

          elseif(Np.eq.92.or.Np.eq.82)then

            zconv(i,j)=0.5*
     *                  (0.5*(3.*zw(i-1,j)-zw(i-2,j))
     *                  +0.5*(3.*zw(i-1,j-1)-zw(i-2,j-1)))

          elseif(Np.eq.93.or.Np.eq.83)then

            zconv(i,j)=0.5*
     *                  (0.5*(3.*zw(i,j-1)-zw(i,j-2))
     *                  +0.5*(3.*zw(i-1,j-1)-zw(i-1,j-2)))

          elseif(Np.eq.4.or.Np.eq.8)then

            zconv(i,j)=0.5*
     *                  (0.5*(3.*zw(i,j)-zw(i+1,j))
     *                  +0.5*(3.*zw(i,j-1)-zw(i+1,j-1)))

          elseif(Np.eq.5)then
	
	    zz1=0.5*(3.*zw(i,j)-zw(i,j+1))
	    zz2=0.5*(3.*zw(i+1,j)-zw(i+1,j+1))
            zconv(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.96.or.Np.eq.86)then

	    zz1=0.5*(3.*zw(i-1,j)-zw(i-1,j+1))
	    zz2=0.5*(3.*zw(i-2,j)-zw(i-2,j+1))
            zconv(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.97.or.Np.eq.87)then

	    zz1=0.5*(3.*zw(i-1,j-1)-zw(i-1,j-2))
	    zz2=0.5*(3.*zw(i-2,j-1)-zw(i-2,j-2))
            zconv(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.98.or.Np.eq.88)then

	    zz1=0.5*(3.*zw(i,j-1)-zw(i,j-2))
	    zz2=0.5*(3.*zw(i+1,j-1)-zw(i+1,j-2))
            zconv(i,j)=0.5*(3.*zz1-zz2)

          elseif(Np.eq.9)then                                           

	    zconv(i,j)=0.5*
     *                  (0.5*(3.*zw(i-1,j)-zw(i-2,j))
     *                  +0.5*(3.*zw(i-1,j-1)-zw(i-2,j-1)))

          endif

        enddo
        enddo
       
! ---

        do j=1,n
        do i=1,m
	  Pass_eta(i,j) = zconv(i,j)
        enddo
        enddo

        do j=1,n
        do i=1,m
	Pass_U(i,j) = uconv(i,j)
	Pass_V(i,j) = vconv(i,j)
        enddo
        enddo

c only for Sandford group

c        do j=1,n
c        do i=1,m
c        Pass_Tau_x(i,j)=-fric_x(i,j)
c        Pass_Tau_y(i,j)=-fric_y(i,j)
c        enddo
c        enddo

        return

	end 

! -------------------------------
	subroutine adaptor_circ_to_sedi
	use sadi

        do j=1,n
        do i=1,m
	Pass_Ub(i,j) = U_3d(i,j,k_3d)
	Pass_Vb(i,j) = V_3d(i,j,k_3d)
        Pass_fw(i,j) = f_cwc
        enddo
        enddo

        return
        end


        subroutine convert_Tau
        use sadi
        integer Np,Np1

        do j=1,n
        do i=1,m
         Tau_1(i,j)=Tau_1(i,j)*(UU(i,j)+ff_11(i,j)*rJ(i,j))*mask(i,j)
         Tau_2(i,j)=Tau_2(i,j)*(VV(i,j)+ff_12(i,j)*rJ(i,j))*mask(i,j)
        enddo
        enddo

c -- Tau_x 
c       Tau_x will be the bottom stress, not Tau_x/H
c ---

        do j=1,n-1
        do i=1,m
          Np=Ntype(i,j)
          if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
          if(Np.eq.0.or.Np.eq.11.or.Np.eq.16.or.Np.eq.12)then
          H_=0.5*(HH(i,j)+HH(i-1,j))
          if(H_.lt.hi)H_=hi
c --- if want Tau_x, replace '1.' with 'H_'
          HJ=1.*0.5*(rJ(i,j)+rJ(i-1,j))
          Tau_x(i,j)=0.5*(Xxi(i,j)+Xxi(i-1,j))*Tau_1(i,j)/HJ
     *          +0.5*(Xeta(i,j)+Xeta(i-1,j))*0.25
     *          *(Tau_2(i,j)+Tau_2(i,j+1)+Tau_2(i-1,j+1)
     *           +Tau_2(i-1,j))/HJ
          endif
          if(Np.eq.13.or.Np.eq.17)then
          H_=0.5*(HH(i,j)+HH(i-1,j))
          if(H_.lt.hi)H_=hi
          HJ=1.*0.5*(rJ(i,j)+rJ(i-1,j))
          Tau_x(i,j)=0.5*(Xxi(i,j)+Xxi(i-1,j))*Tau_1(i,j)/HJ
     *          +0.5*(Xeta(i,j)+Xeta(i-1,j))*0.25
     *          *(Tau_2(i,j)+0.+0.+Tau_2(i-1,j))/HJ
          endif
          if(Np.eq.14.or.Np.eq.15)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i,j)
          Tau_x(i,j)=Xxi(i,j)*Tau_1(i,j)/HJ
     *          +Xeta(i,j)*0.5
     *          *(Tau_2(i,j)+Tau_2(i,j+1))/HJ
          endif
          if(Np.eq.18)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i,j)
          Tau_x(i,j)=Xxi(i,j)*Tau_1(i,j)/HJ
     *          +Xeta(i,j)*0.5
     *          *(Tau_2(i,j)+0.)/HJ
          endif
          if(Np.eq.88.or.Np.eq.83.or.Np.eq.87)then
          H_=HH(i,j-1)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i,j-1)
          Tau_x(i,j)=0.
     *          +Xeta(i,j-1)*0.5
     *          *(Tau_2(i,j)+0.)/HJ
          endif

          if(Np.eq.86.or.Np.eq.82)then
          H_=HH(i-1,j)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i-1,j)
          Tau_x(i,j)=Xxi(i-1,j)*Tau_1(i,j)/HJ
          endif
        enddo
        enddo

c --- Tau_x at c points

        do j=1,n
        do i=1,m

	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	  
         
          if(Np.eq.0.or.Np.eq.2.or.Np.eq.7.or.Np.eq.92
     *      .or.Np.eq.3.or.Np.eq.8.or.Np.eq.4.or.Np.eq.9)then   

            fric_x(i,j)=0.5*(Tau_x(i,j)+Tau_x(i,j-1))

          elseif(Np.eq.5.or.Np.eq.1.or.Np.eq.6.or.Np.eq.96)then

            fric_x(i,j)=0.5*(3.*Tau_x(i,j)-Tau_x(i,j+1))

          elseif(Np.eq.98.or.Np.eq.93.or.Np.eq.97)then

            fric_x(i,j)=0.5*(3.*Tau_x(i,j-1)-Tau_x(i,j-2))

          endif

        enddo
        enddo


c -- Tau_y
        do j=1,n
        do i=1,m-1
          Np=Ntype(i,j)
          if(Np.ge.20.and.Np.le.80)Np=10+mod(Np,10)
          if(Np.eq.0.or.Np.eq.14.or.Np.eq.18.or.Np.eq.13)then
          H_=0.5*(HH(i,j)+HH(i,j-1))
          if(H_.lt.hi)H_=hi
          HJ=1.*0.5*(rJ(i,j)+rJ(i,j-1))
          Tau_y(i,j)=0.5*(Yeta(i,j)+Yeta(i,j-1))*Tau_2(i,j)/HJ
     *          +0.5*(Yxi(i,j)+Yxi(i,j-1))*0.25
     *          *(Tau_1(i,j)+Tau_1(i+1,j)+Tau_1(i+1,j-1)
     *          +Tau_1(i,j-1))/HJ
          endif
          if(Np.eq.12.or.Np.eq.17)then
          H_=0.5*(HH(i,j)+HH(i,j-1))
          if(H_.lt.hi)H_=hi
          HJ=1.*0.5*(rJ(i,j)+rJ(i,j-1))
          Tau_y(i,j)=0.5*(Yeta(i,j)+Yeta(i,j-1))*Tau_2(i,j)/HJ
     *          +0.5*(Yxi(i,j)+Yxi(i,j-1))*0.25
     *          *(Tau_1(i,j)+0.+0.+Tau_1(i,j-1))/HJ
          endif
          if(Np.eq.11.or.Np.eq.15)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i,j)
          Tau_y(i,j)=Yeta(i,j)*Tau_2(i,j)/HJ
     *          +Yxi(i,j)*0.5
     *          *(Tau_1(i,j)+Tau_1(i+1,j))/HJ
          endif
          if(Np.eq.18)then
          H_=HH(i,j)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i,j)
          Tau_y(i,j)=Yeta(i,j)*Tau_2(i,j)/HJ
     *          +Yxi(i,j)*0.5
     *          *(Tau_1(i,j)+0.)/HJ
          endif
          if(Np.eq.88.or.Np.eq.83.or.Np.eq.87)then
          H_=HH(i,j-1)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i,j-1)
          Tau_y(i,j)=Yeta(i,j-1)*Tau_2(i,j)/HJ
          endif

          if(Np.eq.86.or.Np.eq.82)then
          H_=HH(i-1,j)
          if(H_.lt.hi)H_=hi
          HJ=1.*rJ(i-1,j)
          Tau_y(i,j)=
     *          +Yxi(i-1,j)*0.5
     *          *(Tau_1(i,j)+0.)/HJ
          endif

        enddo
        enddo

c Tau_y at c points

        do j=1,n
        do i=1,m

 	  if(Ntype1(i,j).lt.80)then
            Np=mod(Ntype1(i,j),10)
	  else
	    Np=90+mod(Ntype1(i,j),10)
	  endif	            

          if(Np.eq.0.or.Np.eq.3.or.Np.eq.7.or.Np.eq.93
     *      .or.Np.eq.2.or.Np.eq.6.or.Np.eq.1)then ! change 19

            fric_y(i,j)=0.5*(Tau_y(i,j)+Tau_y(i-1,j))

	  elseif(Np.eq.9)then

	    fric_y(i,j)=0.

          elseif(Np.eq.5.or.Np.eq.4.or.Np.eq.8.or.Np.eq.98)then

            fric_y(i,j)=0.5*(3.*Tau_y(i,j)-Tau_y(i+1,j))

          elseif(Np.eq.96.or.Np.eq.92.or.Np.eq.97)then

            fric_y(i,j)=0.5*(3.*Tau_y(i-1,j)-Tau_y(i-2,j))

          endif

        enddo
        enddo

c --- using bottom velocity calculate Tau_x and Tau_y

        goto 812

        do j=1,ny_circ
        do i=1,nx_circ
        Tau_x(i,j)=f_cwc*sqrt(u_3d(i,j,k_3d)**2+v_3d(i,j,k_3d)**2)
     &             *u_3d(i,j,k_3d)
        Tau_y(i,j)=f_cwc*sqrt(u_3d(i,j,k_3d)**2+v_3d(i,j,k_3d)**2)
     &             *v_3d(i,j,k_3d)
        enddo
        enddo

812     continue

        return
        end

c -----------
        subroutine initialize_sedflux
 	use sadi

 	do j=1,ny_circ
   	  do i=1,nx_circ
           pass_sedfluxcum_x(i,j) = 0.
           pass_sedfluxcum_y(i,j) = 0.
           end do
          end do
         
        return
        end

c ------------
         subroutine solve_3d
         use sadi
         real h3d,z3d,anu3d(mtrim,ntrim)
         integer kkk

c --  in Cartesian

          do j=1,n
           do i=1,m
             anu3d(i,j)=anu(i,j)
           enddo
          enddo

c --  cartesian

          do j=1,n
           do i=1,m
             h3d=max(depth_min,H(i,j))

             ff_11_cart(i,j)=-(Tau_x(i,j)-0.)/3./anu3d(i,j)*h3d
     &          -Intp_MassFluxU_Circ(i,j)/h3d
             ee_11_cart(i,j)=(Tau_x(i,j)-0.)/anu3d(i,j)
             dd_11_cart(i,j)=(Tau_x(i,j)-0.)/2./anu3d(i,j)/h3d

             ff_12_cart(i,j)=-(Tau_y(i,j)-0.)/3./anu3d(i,j)*h3d
     &          -Intp_MassFluxV_Circ(i,j)/h3d
             ee_12_cart(i,j)=(Tau_y(i,j)-0.)/anu3d(i,j)
             dd_12_cart(i,j)=(Tau_y(i,j)-0.)/2./anu3d(i,j)/h3d
           enddo
          enddo

c -- boundaries
           j=n
           do i=1,m
             ff_11_cart(i,j)=ff_11_cart(i,j-1)
             ee_11_cart(i,j)=ee_11_cart(i,j-1)
             dd_11_cart(i,j)=dd_11_cart(i,j-1)
          enddo

           i=m
           do j=1,n
             ff_12_cart(i,j)=ff_12_cart(i-1,j)
             ee_12_cart(i,j)=ee_12_cart(i-1,j)
             dd_12_cart(i,j)=dd_12_cart(i-1,j)
           enddo

c --- contravariant

        do j=1,n
        do i=1,m
	ff_11(i,j)=-X_xi_2(i,j)/sq_g_0(i,j)*ff_12_cart(i,j) 
     *           +Y_xi_2(i,j)/sq_g_0(i,j)*ff_11_cart(i,j)
	ff_12(i,j)= X_xi_1(i,j)/sq_g_0(i,j)*ff_12_cart(i,j) 
     *           -Y_xi_1(i,j)/sq_g_0(i,j)*ff_11_cart(i,j)
	ee_11(i,j)=-X_xi_2(i,j)/sq_g_0(i,j)*ee_12_cart(i,j) 
     *           +Y_xi_2(i,j)/sq_g_0(i,j)*ee_11_cart(i,j)
	ee_12(i,j)= X_xi_1(i,j)/sq_g_0(i,j)*ee_12_cart(i,j) 
     *           -Y_xi_1(i,j)/sq_g_0(i,j)*ee_11_cart(i,j)
	dd_11(i,j)=-X_xi_2(i,j)/sq_g_0(i,j)*dd_12_cart(i,j) 
     *           +Y_xi_2(i,j)/sq_g_0(i,j)*dd_11_cart(i,j)
	dd_12(i,j)= X_xi_1(i,j)/sq_g_0(i,j)*dd_12_cart(i,j) 
     *           -Y_xi_1(i,j)/sq_g_0(i,j)*dd_11_cart(i,j)
        enddo
        enddo

c ---     get 3d U and V in Cartesian coordinates
          do j=1,ny_Circ
          do i=1,nx_Circ
            if(H(i,j).gt.0.)then
              do kkk=1,k_3d
              z3d=H(i,j)-H(i,j)*(kkk-1.)/(k_3d-1.)
              U_3d(i,j,kkk)=dd_11_cart(i,j)*z3d*z3d+ee_11_cart(i,j)*z3d
     &                      +ff_11_cart(i,j)
              V_3d(i,j,kkk)=dd_12_cart(i,j)*z3d*z3d+ee_12_cart(i,j)*z3d
     &                      +ff_12_cart(i,j)
              enddo
            else
              do kkk=1,k_3d
              U_3d(i,j,kkk)=0.
              V_3d(i,j,kkk)=0.
              enddo
            endif              
          enddo
          enddo

          do j=1,ny_circ
          do i=1,nx_circ
          do kkk=1,k_3d
            U_3d(i,j,kkk)=uconv(i,j)+U_3d(i,j,kkk)
            V_3d(i,j,kkk)=vconv(i,j)+V_3d(i,j,kkk)            
          enddo
          enddo
          enddo

         return
         end

         subroutine filter5point(vari,mf,nf,numf)
         use sadi
         real tmpf(mtrim,ntrim),vari(mtrim,ntrim)
         integer mf,nf,numf,kf
         
         do kf=1,numf
           do j=1,nf
           do i=1,mf
             tmpf(i,j)=vari(i,j)
           enddo
           enddo
           do j=2,nf-1
           do i=2,mf-1
             vari(i,j)=0.5*tmpf(i,j)+0.125*(tmpf(i-1,j)+tmpf(i+1,j)
     *                   +tmpf(i,j+1)+tmpf(i,j-1))
           enddo
           enddo
           do i=2,mf-1
            vari(i,1)=0.5*tmpf(i,1)+0.25*(tmpf(i-1,1)+tmpf(i+1,1))
            vari(i,nf)=0.5*tmpf(i,nf)+0.25*(tmpf(i-1,nf)+tmpf(i+1,nf))
           enddo
           do j=2,nf-1
            vari(1,j)=0.5*tmpf(1,j)+0.25*(tmpf(1,j-1)+tmpf(1,j+1))     
            vari(mf,j)=0.5*tmpf(mf,j)+0.25*(tmpf(mf,j-1)+tmpf(mf,j+1))     
           enddo

         enddo 

         return
         end

c ---------------------------------------------------

         subroutine dispersion
         use sadi
         real htot

c --- D_ab terms

        do j=1,n
        do i=1,m
        if (H(i,j).gt.0.)then
	D_11(i,j)=1./anu(i,j)*(
     *      1./63.*dd_11(i,j)*dd_11(i,j)*(H(i,j))**6 
     *     +1./36.*(dd_11(i,j)*ee_11(i,j)+dd_11(i,j)*ee_11(i,j))
     *     *(H(i,j))**5
     *     +(1./15.*dd_11(i,j)*ff_11(i,j)+1./15.*dd_11(i,j)*ff_11(i,j)
     *     +1./20.*ee_11(i,j)*ee_11(i,j))
     *       *(H(i,j))**4
     *     +1./8.*(ee_11(i,j)*ff_11(i,j)+ee_11(i,j)*ff_11(i,j))
     *       *(H(i,j))**3
     *     +1./3.*ff_11(i,j)*ff_11(i,j)*(H(i,j))**2
     *         )

	D_12(i,j)=1./anu(i,j)*(
     *      1./63.*dd_11(i,j)*dd_12(i,j)*(H(i,j))**6 
     *     +1./36.*(dd_11(i,j)*ee_12(i,j)+dd_12(i,j)*ee_11(i,j))
     *      *(H(i,j))**5
     *     +(1./15.*dd_11(i,j)*ff_12(i,j)+1./15.*dd_12(i,j)*ff_11(i,j)
     *     +1./20.*ee_11(i,j)*ee_12(i,j))
     *       *(H(i,j))**4
     *     +1./8.*(ee_11(i,j)*ff_12(i,j)+ee_12(i,j)*ff_11(i,j))
     *       *(H(i,j))**3
     *     +1./3.*ff_11(i,j)*ff_12(i,j)*(H(i,j))**2
     *         )

! ---   D_12=D_21

	D_21(i,j)=D_12(i,j)

	D_22(i,j)=1./anu(i,j)*(
     *      1./63.*dd_12(i,j)*dd_12(i,j)*(H(i,j))**6 
     *     +1./36.*(dd_12(i,j)*ee_12(i,j)+dd_12(i,j)*ee_12(i,j))
     *       *(H(i,j))**5
     *     +(1./15.*dd_12(i,j)*ff_12(i,j)+1./15.*dd_12(i,j)*ff_12(i,j)
     *     +1./20.*ee_12(i,j)*ee_12(i,j))
     *       *(H(i,j))**4
     *     +1./8.*(ee_12(i,j)*ff_12(i,j)+ee_12(i,j)*ff_12(i,j))
     *       *(H(i,j))**3
     *     +1./3.*ff_12(i,j)*ff_12(i,j)*(H(i,j))**2
     *         )
        else
         D_11(i,j)=0.
         D_12(i,j)=0.
         D_21(i,j)=0.
         D_22(i,j)=0.
        endif

        enddo
        enddo

c --- have to filter D_ab becasue of high order polynomial

        call filter5point(D_11,m,n,10)
        call filter5point(D_12,m,n,10)
        call filter5point(D_22,m,n,10)
        call filter5point(D_21,m,n,10)

! --- c-grid requires D_ab_c and D_ab_z, so calculate D_ab_z now

	do j=1,n-1
	do i=1,m-1
	  D_11_z(i,j)=0.25
     *      *(D_11(i,j)+D_11(i,j+1)+D_11(i+1,j+1)+D_11(i+1,j))
	  D_12_z(i,j)=0.25
     *      *(D_12(i,j)+D_12(i,j+1)+D_12(i+1,j+1)+D_12(i+1,j))
	  D_22_z(i,j)=0.25
     *      *(D_22(i,j)+D_22(i,j+1)+D_22(i+1,j+1)+D_22(i+1,j))
	enddo
	enddo


c ---  Dup11 at z point i.e. d(Q1/h)/dx

	do j=1,n
	do i=1,m
	  htot=(depth_u(i,j)+zeta_u(i,j))
	  if(htot.lt.hi)htot=hi
          tr1(i,j)=Q1(i,j)/htot
	enddo
	enddo

        call der_xi_at_z

	do j=1,n
	do i=1,m
	  htot=(depth(i,j)+zeta(i,j))
	  if(htot.lt.hi)htot=hi
           Dup11(i,j)=2.*sq_g_0_z(i,j)*htot*D_11_z(i,j)*trd1(i,j)
	enddo
	enddo

! ---  Dup11 _xi

	do j=1,n-1
	do i=1,m-1
        tr1(i,j)=Dup11(i,j)
	enddo
	enddo

        call der_xi_at_u_symmetry  

	do j=1,n-1
	do i=1,m
        dispx(i,j)=trd1(i,j)
	enddo
	enddo

c ---  Dup12 at c point i.e. d(Q1/h)/dy

	do j=1,n-1
	do i=1,m
	htot=(depth_u(i,j)+zeta_u(i,j))
	if(htot.lt.hi)htot=hi
        tr1(i,j)=Q1(i,j)/htot
	enddo
	enddo

        call der_eta_at_c_symmetry

	do j=1,n
	do i=1,m
	htot=(depth_c(i,j)+zeta_c(i,j))
	if(htot.lt.hi)htot=hi
        Dup12(i,j)=sq_g_0(i,j)*D_22(i,j)*htot*trd1(i,j)
	enddo
	enddo

! ---  Dup12 _eta

	do j=1,n
	do i=1,m
        tr1(i,j)=Dup12(i,j)
	enddo
	enddo

        call der_eta_at_u   

	do j=1,n-1
	do i=1,m
        dispx(i,j)=dispx(i,j)+trd1(i,j)
	enddo
	enddo

! --- calculate Dispy

! ---   Dup22 at z i.e. d(Q2/h)/dy

	do j=1,n
	do i=1,m-1
	  htot=(depth_v(i,j)+zeta_v(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=Q2(i,j)/htot
	enddo
	enddo

        call der_eta_at_z

	do j=1,n-1
	do i=1,m-1
	  htot=(depth(i,j)+zeta(i,j))
	  if(htot.lt.hi)htot=hi	
        Dup22(i,j)=2.*D_22(i,j)*sq_g_0_z(i,j)*htot*trd1(i,j)
	enddo
	enddo

c ---   Dup22 _eta

	do j=1,n-1
	do i=1,m-1
	  htot=(depth(i,j)+zeta(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=Dup22(i,j)
	enddo
	enddo

        call der_eta_at_v_symmetry  

	do j=1,n
	do i=1,m-1
          dispy(i,j)=trd1(i,j)
	enddo
	enddo

! ---   Dup12 at c i.e. d(Q2/h)/dx

	do j=1,n-1
	do i=1,m
	  htot=(depth_v(i,j)+zeta_v(i,j))
	  if(htot.lt.hi)htot=hi	
        tr1(i,j)=Q2(i,j)/htot
	enddo
	enddo

        call der_xi_at_c_symmetry

	do j=1,n
	do i=1,m
	  htot=(depth_c(i,j)+zeta_c(i,j))
	  if(htot.lt.hi)htot=hi	
        Dup12(i,j)=D_11(i,j)*sq_g_0(i,j)*htot*trd1(i,j)
	enddo
	enddo

c ---  Dup12 _xi

	do j=1,n
	do i=1,m
        tr1(i,j)=Dup12(i,j)
	enddo
	enddo

        call der_xi_at_v  

	do j=1,n
	do i=1,m-1
          dispy(i,j)=dispy(i,j)+trd1(i,j)
	enddo
	enddo


1000    format(500f12.6)


	return

	end


! ---------------------------- output station

        subroutine output_station(time_c,n_st,i_st,j_st)
         use sadi
         real  time_c
         integer n_st,i_st(nsta_max),j_st(nsta_max)

         write(83,100)time_c, (z(i_st(i),j_st(i)),i=1,n_st)
	   write(84,100)time_c, (u(i_st(i),j_st(i)),i=1,n_st)
     &	   ,(v(i_st(i),j_st(i)),i=1,n_st)


100      format(500f16.4)
        return
        end

! ----------------------------- calculate wind
 
        subroutine calculate_wind
        use sadi
        real wx,wy,wsp,wdir,wpar
      	wpar=0.00063

       
	if(windtime(1).ge.0.)then
         wsp=windspd(1)+(windspd(2)-windspd(1))/
     &	   (windtime(2)-windtime(1))*wind_intv/3600.
	   wdir=winddir(1)+(winddir(2)-winddir(1))/
     &	   (windtime(2)-windtime(1))*wind_intv/3600.
	   wx=-wsp*sin(wdir*3.1415926/180.)
	   wy=-wsp*cos(wdir*3.1415926/180.)
         Tao_ax=wpar*wsp*wx
	   Tao_ay=wpar*wsp*wy
	else 
        Tao_ax=0.
	  Tao_ay=0.
	endif

	  return
	  end
