         parameter(m=220,n=100)
         real x(m,n),y(m,n),cori(m,n),dep(m,n),dep_swan(m,n)
         real h(m,n)

         depmax=10.
         depmin=0.002
         slope=1./30.0

         open(1,file='../grid/x.dat')
          do i=1,m
           read(1,*)x(i,1)
          enddo
         close(1)

         open(1,file='../grid/y.dat')
          do j=1,n
           read(1,*)y(1,j)
          enddo
         close(1)

         do j=1,n
           do i=1,m
             x(i,j)=x(i,1)
             y(i,j)=y(1,j)
             cori(i,j)=0.0
           enddo
         enddo

         width_surf=500.0
         width_inlet=1000.0
         rlength_inlet = 500.0
         dep_inlet=5.0
         dep_basin=5.0
         flat=10.0

!   shoal
         a_in=600.0
         a_out=850.0
         b_in=200.0
         b_out=450.0
         height=9.0
         shoal_sitting = 10.0
         pi=3.1415926

         do j=1,n
         do i=1,m
!         dep(i,j)=hconst(i,j)    
         r=sqrt(x(i,j)*x(i,j)+y(i,j)*y(i,j))
         phi=atan2(y(i,j),x(i,j))
         r_out=a_out*a_out*b_out*b_out/
     &         (a_out*a_out*sin(phi)*sin(phi)
     &         +b_out*b_out*cos(phi)*cos(phi))
         r_out=sqrt(r_out)
         r_in=a_in*a_in*b_in*b_in/
     &         (a_in*a_in*sin(phi)*sin(phi)
     &         +b_in*b_in*cos(phi)*cos(phi))
         r_in=sqrt(r_in)

         dis=r_out-r
         rate=r_out-r_in
         h(i,j)=shoal_sitting-height*sin(pi/rate*dis)
!         print*,i,j,r,r_in,r_out,dis,h(i,j)
         if(r.gt.r_out.or.r.lt.r_in) h(i,j)=shoal_sitting
         if(y(i,j).gt.0.)h(i,j)=shoal_sitting
        enddo
        enddo

        open(1,file='tmp.txt')
         do j=1,n
           write(1,100)(h(i,j),i=1,m)
         enddo
        close(1)
        

         do j=1,n
         do i=1,m
          dep(i,j)=flat
! slope
          if(y(i,j).gt.-width_surf)then
            dep(i,j)=flat-flat*(y(i,j)+width_surf)/width_surf
          endif
! channel
          if(abs(x(i,j)).le.500.0)then
            if(dep(i,j).lt.dep_inlet)dep(i,j)=dep_inlet
          endif
! basin
          if(y(i,j).ge.rlength_inlet)then
            dep(i,j)=dep_basin
          endif

         enddo
         enddo


! combine

        do j=1,n
         do i=1,m
           dep(i,j)=min(dep(i,j),h(i,j))
         enddo
        enddo

    
! for swan, make a minimum depth to get equal grid decomposition
       do j=1,n
        do i=1,m
         if(dep(i,j).lt.-2.0)dep(i,j)=-2.0
         dep_swan(i,j)=dep(i,j)
         if(dep_swan(i,j).lt.0.001) dep_swan(i,j)=0.001
        enddo
        enddo

        open(1,file='x_str.txt')
         do j=1,n
           write(1,100)(x(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='y_str.txt')
         do j=1,n
           write(1,100)(y(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='cori_str.txt')
         do j=1,n
           write(1,100)(cori(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='dep_shoal_inlet.txt')
         do j=1,n
           write(1,100)(dep(i,j),i=1,m)
         enddo
        close(1)

! for swan
         open(1,file='xxyy_swan_str.txt')
         do j=1,n
         write(1,100)(x(i,j),i=1,m)
         enddo
         do j=1,n
         write(1,100)(y(i,j),i=1,m)
         enddo
         close(1)

         open(1,file='dep_swan_shoal_inlet.txt')
         do j=1,n
         write(1,100)(dep_swan(i,j),i=1,m)
         enddo
         close(1)

100     format(3000f12.3)

        end
