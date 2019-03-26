         parameter(nx=309,ny=288)
         real x(nx,ny),y(nx,ny),depth(nx,ny)
         integer ntype(nx,ny)


c	open(2,file='depth_5.dat')
c	do j=1,ny
c	  read(2,*)(depth(i,j),i=1,nx)
c	enddo
110	format(501f12.6)
c	close(2)

        open(1,file='depth_09_07_s_2.dat')
        do j=1,ny
        read(1,*)(depth(i,j),i=1,nx)
190     format(801f16.8)
        enddo
        close(1)

	open(2,file='xymast_09_07.dat')
	do j=1,ny
	  read(2,*)(x(i,j),i=1,nx)
	enddo
	do j=1,ny
	  read(2,*)(y(i,j),i=1,nx)
	enddo

	close(2)
120 	format(501f12.3)

	open(2,file='ntype_09_07.dat')
	do j=ny,1,-1
	read(2,300)(ntype(i,j),i=1,nx)
	enddo
300	format(500I3)
	close(2)  

 
       do j=1,ny
        do i=1,nx
         if(depth(i,j).gt.100.)depth(i,j)=100.
         if(depth(i,j).lt.0.1)depth(i,j)=0.1
         if(ntype(i,j).eq.99)depth(i,j)=-99999.0
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
100      format(500f10.2)
         close(1)

         open(1,file='depp0.dat')
         do j=1,ny,nsky
         write(1,100)(depth(i,j),i=1,nx,nskx)
         enddo
         close(1)

         end
         
