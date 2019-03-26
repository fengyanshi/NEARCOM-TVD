	parameter(L=308,M=287,ncon=10)
        parameter(L2big=2*L-2,M2big=2*M-2)
        real dr(0:L,0:M),xr(0:L,0:M),yr(0:L,0:M)
        real amp_north(0:L,1:ncon), amp_west(0:M,1:ncon), 
     &       amp_south(0:L,1:ncon)
        real pha_north(0:L,1:ncon), pha_west(0:M,1:ncon), 
     &       pha_south(0:L,1:ncon)
        real amp1(ncon),pha1(ncon)
        real amp2(ncon),pha2(ncon)
        real amp3(ncon),pha3(ncon)
        real amp4(ncon),pha4(ncon)
	real wamp1(ncon),wpha1(ncon)
        real wamp2(ncon),wpha2(ncon)
        real wamp3(ncon),wpha3(ncon)
        real wamp4(ncon),wpha4(ncon)
					
        character*7 AA

c        data amp1/0.535, 0.368,0.227,0.130,0.120,0.119,0.040,0.036/
c        data pha1/191.8, 220.7,206.6,199.9,219.8,165.0,198.1,188.3/
c        data amp2/0.535, 0.368,0.227,0.130,0.120,0.119,0.040,0.036/
c        data pha2/191.8, 220.7,206.6,199.9,219.8,165.0,198.1,188.3/
c        data amp3/0.516,0.366,0.226,0.128,0.120,0.116,0.040,0.036/
c        data pha3/186.3,219.0,204.7,190.8,218.1,160.1,196.1,179.7/
c        data amp4/0.516,0.366,0.226,0.128,0.120,0.116,0.040,0.036/
c        data pha4/186.3,219.0,204.7,190.8,218.1,160.1,196.1,179.7/

c        AA='ut     '
c        write(*,199)AA,50.2,34.8
c        stop

c 1A
        open(1,file='tide_cons98.dat')
        print*,'98 CONSTITUENTS!!!!!!'
        read(1,*)
        read(1,*)
        do i=1,ncon
          read(1,*)AA, amp2(i),pha2(i)
           print*,i,AA,amp2(i),pha2(i)
        enddo

c 1B
        read(1,*)
        read(1,*)
        do i=1,ncon
          read(1,*)AA, amp1(i),pha1(i)
           print*,i,AA,amp1(i),pha1(i)
        enddo


c 2A
        read(1,*)
        read(1,*)
        do i=1,ncon
          read(1,*)AA, amp3(i),pha3(i)
           print*,i,AA,amp3(i),pha3(i)
        enddo

        do i=1,26
          read(1,*)
        enddo

        do i=1,ncon
          read(1,*)AA, amp4(i),pha4(i)
           print*,i,AA,amp4(i),pha4(i)
        enddo
        close(1)

c -- read new tide cons
        open(1,file='tide_cons05.dat')
        print*,'05 CONSTITUENTS!!!!!!'
        read(1,*)
        read(1,*)
        do i=1,ncon
          read(1,*)AA, wamp2(i),wpha2(i)
             print*,i,AA,wamp2(i),wpha2(i)
        enddo
c 1B
        read(1,*)
        read(1,*)
        do i=1,ncon
          read(1,*)AA, wamp1(i),wpha1(i)
             print*,i,AA,wamp1(i),wpha1(i)
        enddo
c 2A
        read(1,*)
        read(1,*)
        do i=1,ncon
          read(1,*)AA, wamp3(i),wpha3(i)
             print*,i,AA,wamp3(i),wpha3(i)
        enddo

            do i=1,26
               read(1,*)
            enddo

            do i=1,ncon
                 read(1,*)AA, wamp4(i),wpha4(i)
            print*,i,AA,wamp4(i),wpha4(i)
            enddo
            close(1)
    
c - read over

c -- adjust
         do i=1,ncon
           pha1(i)=wpha1(i)
           pha2(i)=wpha2(i)
           pha3(i)=wpha3(i)
           pha4(i)=wpha4(i)
	   
	   amp1(i)=wamp1(i)
	   amp2(i)=wamp2(i)
	   amp3(i)=wamp3(i)
	   amp4(i)=wamp4(i)
         enddo

199     format(A3,f5.3, f7.3)


	open(1,file='xymast_09_07.dat')
 	  do j=0,M
	    read(1,*)(xr(i,j),i=0,L)
	  enddo
 	  do j=0,M
	    read(1,*)(yr(i,j),i=0,L)
	  enddo
	close(1)

100	format(500f16.2)

        m_sw=0
        n_sw=0
        m_se=112
        n_se=0
        m_nw=0
        n_nw=287
        m_ne=40
        n_ne=287
        
        do i=0,L
         do k=1,ncon
          amp_north(i,k)=0.
          pha_north(i,k)=0.
          amp_south(i,k)=0.
          pha_south(i,k)=0.
         enddo
        enddo

        do j=0,M
          do k=1,ncon
          amp_west(j,k)=0.
          pha_west(j,k)=0.
          enddo
        enddo

c - west
        do k=1,ncon
         do j=n_sw,n_nw
          dis=sqrt((xr(m_sw,n_sw)-xr(m_nw,n_nw))**2
     &            +(yr(m_sw,n_sw)-yr(m_nw,n_nw))**2)
          r=sqrt((xr(m_sw,j)-xr(m_sw,n_sw))**2
     &            +(yr(m_sw,j)-yr(m_sw,n_sw))**2)
          amp_west(j,k)=amp3(k)+r/dis*(amp2(k)-amp3(k))
          pha_west(j,k)=pha3(k)+r/dis*(pha2(k)-pha3(k))
         enddo
        enddo

c - south
        do k=1,ncon
         do i=m_sw,m_se
          dis=sqrt((xr(m_sw,n_sw)-xr(m_se,n_se))**2
     &            +(yr(m_sw,n_sw)-yr(m_se,n_se))**2)
          r=sqrt((xr(i,n_sw)-xr(m_sw,n_sw))**2
     &            +(yr(i,n_sw)-yr(m_sw,n_sw))**2)
          amp_south(i,k)=amp3(k)+r/dis*(amp4(k)-amp3(k))
          pha_south(i,k)=pha3(k)+r/dis*(pha4(k)-pha3(k))
         enddo
        enddo         
       
        do k=1,ncon
         do i=m_se+1,L
          amp_south(i,k)=amp_south(m_se,k)
          pha_south(i,k)=pha_south(m_se,k)
        enddo
        enddo

c - north
        do k=1,ncon
         do i=m_nw,m_ne
          dis=sqrt((xr(m_nw,n_nw)-xr(m_ne,n_ne))**2
     &            +(yr(m_nw,n_nw)-yr(m_ne,n_ne))**2)
          r=sqrt((xr(i,n_nw)-xr(m_nw,n_nw))**2
     &            +(yr(i,n_nw)-yr(m_nw,n_nw))**2)
          amp_north(i,k)=amp2(k)+r/dis*(amp1(k)-amp2(k))
          pha_north(i,k)=pha2(k)+r/dis*(pha1(k)-pha2(k))
         enddo
        enddo

        do k=1,ncon
         do i=m_ne+1,L
          amp_north(i,k)=amp_north(m_ne,k)
          pha_north(i,k)=pha_north(m_ne,k)
        enddo
        enddo

         open(1,file='amp_west.dat')
         do j=0,M
         write(1,200)(amp_west(j,k),k=1,ncon)
         enddo
         close(1)

         open(1,file='pha_west.dat')
         do j=0,M
         write(1,200)(pha_west(j,k),k=1,ncon)
         enddo
         close(1)

         open(1,file='amp_south.dat')
         do i=m_sw,L
         write(1,200)(amp_south(i,k),k=1,ncon)
         enddo
         close(1)

         open(1,file='pha_south.dat')
         do i=m_sw,L
         write(1,200)(pha_south(i,k),k=1,ncon)
         enddo
         close(1)

         open(1,file='amp_north.dat')
         do i=m_nw,L
         write(1,200)(amp_north(i,k),k=1,ncon)
         enddo
         close(1)

         open(1,file='pha_north.dat')
         do i=m_nw,L
         write(1,200)(pha_north(i,k),k=1,ncon)
         enddo
         close(1)


200      format(100f12.6)

        end

