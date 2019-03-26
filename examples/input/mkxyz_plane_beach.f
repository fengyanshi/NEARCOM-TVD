         parameter(m=220,n=100)
         real x(m,n),y(m,n),cori(m,n),dep(m,n),dep_swan(m,n)

         depmax=10.
         depmin=0.002
         slope=1./30.0
         dx=20.0
         dy=20.0

          do i=1,m
          do j=1,n
             x(i,j)=(i-1)*dx
             y(i,j)=(j-n)*dy
             cori(i,j)=0.0
          enddo
          enddo

         width_surf=500.0
         flat=10.0
         do j=1,n
         do i=1,m
          dep(i,j)=flat
! slope
          if(y(i,j).gt.-width_surf)then
            dep(i,j)=flat-flat*(y(i,j)+width_surf)/width_surf
          endif
         enddo
         enddo


    
! for swan, make a minimum depth to get equal grid decomposition
       do j=1,n
        do i=1,m
         dep_swan(i,j)=dep(i,j)
         if(dep_swan(i,j).lt.0.001) dep_swan(i,j)=0.001
        enddo
        enddo

        open(1,file='x.txt')
         do j=1,n
           write(1,100)(x(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='y.txt')
         do j=1,n
           write(1,100)(y(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='cori.txt')
         do j=1,n
           write(1,100)(cori(i,j),i=1,m)
         enddo
        close(1)

        open(1,file='dep.txt')
         do j=1,n
           write(1,100)(dep(i,j),i=1,m)
         enddo
        close(1)

! for swan
         open(1,file='xxyy_swan.txt')
         do j=1,n
         write(1,100)(x(i,j),i=1,m)
         enddo
         do j=1,n
         write(1,100)(y(i,j),i=1,m)
         enddo
         close(1)

         open(1,file='dep_swan.txt')
         do j=1,n
         write(1,100)(dep_swan(i,j),i=1,m)
         enddo
         close(1)

100     format(3000f12.3)

        end
