        parameter(nx=251,ny=401)
        real x(nx,ny),y(nx,ny),depth(nx,ny),dep_swan(nx,ny)
        

        open(1,file='x.dat')
        do i=1,nx
        read(1,*)x(i,1)
        enddo
        close(1)
        open(1,file='y.dat')
        do j=1,ny
        read(1,*)y(1,j)
        enddo
        close(1)
        do j=1,ny
        do i=1,nx
         x(i,j)=x(i,1)
         y(i,j)=y(1,j)
        enddo
        enddo

        do j=1,ny
        do i=1,nx
         x(i,j)=x(i,j)-3000.
        enddo
        enddo

c parameters can be modified
        slope_land=1./1000.
        slope_beach=1./100.
        slope_shelf=1./1000.
        x_shore=0.
        x_beach_toe=1000.
        depth_limit=100.
c ---------------------------
        depth_toe=(x_beach_toe-x_shore)*slope_beach

        
        do j=1,ny
        do i=1,nx
         if(x(i,j).le.x_shore)then
          depth(i,j)=slope_land*(x(i,j)-x_shore)
         endif
         if(x(i,j).gt.x_shore.and.x(i,j).le.x_beach_toe)then
          depth(i,j)=(x(i,j)-x_shore)*slope_beach
          endif
         if(x(i,j).gt.x_beach_toe)then
          depth(i,j)=depth_toe+(x(i,j)-x_beach_toe-x_shore)*slope_shelf
          endif
         if (depth(i,j).gt.depth_limit)then
           depth(i,j)=depth_limit
         endif
          dep_swan(i,j)=depth(i,j)
          if(dep_swan(i,j).lt.0.)dep_swan(i,j)=-99999.0
        enddo
        enddo

        open(2,file='xymast.dat')
        do j=1,ny
          write(2,120)(x(i,j),i=1,nx)
        enddo
        do j=1,ny
          write(2,120)(y(i,j),i=1,nx)
        enddo
120     format(801f16.6)
        close(2)

        open(2,file='depth.dat')
        do j=1,ny
          write(2,120)(depth(i,j),i=1,nx)
        enddo
        close(2)


        open(2,file='depp0.dat')
         do j=1,ny,2
           write(2,120)(dep_swan(i,j),i=1,nx,2)
         enddo
        close(2)

        open(2,file='xxyy0.dat')
         do j=1,ny,2
          write(2,120)(x(i,j),i=1,nx,2)
         enddo
         do j=1,ny,2
          write(2,120)(y(i,j),i=1,nx,2)
         enddo
        close(2)

         end



