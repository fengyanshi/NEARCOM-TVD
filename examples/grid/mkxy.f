       parameter(m=127,n=383)
       real x(m,n),y(m,n)

       open(1,file='x.dat')
         do i=1,m
          read(1,*)xx
          do j=1,n
           x(i,j)=xx
          enddo
         enddo
       close(1)
       open(1,file='y.dat')
         do j=1,n
          read(1,*)yy
          do i=1,m
           y(i,j)=yy
          enddo
         enddo
       close(1)

       end
