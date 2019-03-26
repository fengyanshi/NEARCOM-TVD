         parameter(nx=167,ny=434)
         real x(nx,ny),y(nx,ny),depth(nx,ny),depth_sw(nx,ny),
     &        tmp(nx,ny)
         integer ntype(nx,ny)


c	open(2,file='depth_5.dat')
c	do j=1,ny
c	  read(2,*)(depth(i,j),i=1,nx)
c	enddo
110	format(501f12.6)
c	close(2)

        open(1,file='dep_rect.dat')
        do j=1,ny
        read(1,*)(depth(i,j),i=1,nx)
190     format(801f16.8)
        enddo
        close(1)

        do j=1,ny
	do i=1,nx
	 depth(i,j)=depth(i,j)+1.05
         depth_sw(i,j)=depth(i,j)
	enddo
	enddo

	open(2,file='xy_rect.dat')
	do j=1,ny
	  read(2,*)(x(i,j),i=1,nx)
	enddo
	do j=1,ny
	  read(2,*)(y(i,j),i=1,nx)
	enddo

	close(2)
120 	format(501f12.3)

       do k=1,2
       do j=1,ny
         do i=1,nx
          tmp(i,j)=depth(i,j)
         enddo
       enddo

       do j=2,ny-1
       do i=2,nx-1
         depth(i,j)=0.5*tmp(i,j)+0.125*(
     &    tmp(i+1,j)+tmp(i-1,j)+tmp(i,j+1)+tmp(i,j-1))
       enddo
       enddo
       enddo
 
       do j=1,ny
        do i=1,nx
         if(depth_sw(i,j).lt.0.5)depth_sw(i,j)=-99999.0
c         if(depth(i,j).lt.0.1)depth(i,j)=-3.0
c        if(depth(i,j).lt.0.15.and.depth(i,j).ge.0.1)depth(i,j)=0.15
        enddo
        enddo

         nskx=1
         nsky=1

         open(1,file='xxyy0.dat')
         do j=1,ny,nsky
         write(1,100)(x(i,j),i=1,nx,nskx)
         enddo
         do j=1,ny,nsky
         write(1,100)(y(i,j),i=1,nx,nskx)
         enddo
100      format(500f16.2)
         close(1)

         open(1,file='depp0.dat')
         do j=1,ny,nsky
         write(1,100)(depth_sw(i,j),i=1,nx,nskx)
         enddo
         close(1)

         open(1,file='xymast_rect.dat')
         do j=1,ny
	 write(1,100)(x(i,j),i=1,nx)
	 enddo
	 do j=1,ny
	 write(1,100)(y(i,j),i=1,nx)
	 enddo
	 close(1)

	 open(1,file='depth_rect.dat')
	 do j=1,ny
	 write(1,100)(depth(i,j),i=1,nx)
	 enddo
	 close(1)

         end
         
