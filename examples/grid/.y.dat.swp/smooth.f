         parameter(m=291,n=180)

         real depth(m,n),tmp(m,n)

         open(1,file='depth_09_06.dat')
         do j=1,n
          read(1,*)(depth(i,j),i=1,m)
         enddo
100     format(500f12.3)

         do k=1,4
          do j=1,n
          do i=1,m
           tmp(i,j)=depth(i,j)
          enddo
          enddo

          do j=2,n-1
          do i=2,m-1
           depth(i,j)=0.5*tmp(i,j)+0.125*(tmp(i+1,j)+tmp(i-1,j)
     &                    +tmp(i,j+1)+tmp(i,j-1))
          enddo
          enddo

          enddo

          open(2,file='depth_09_06_s.dat')

          do j=1,n
           write(2,100)(depth(i,j),i=1,m)
          enddo

          end




