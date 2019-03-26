          parameter(m=309,n=288)
          real wind_x(m,n),wind_y(m,n),x(m,n),y(m,n)
     &          ,ang(m,n)
         character*4 file_name
          character*22 folder
          folder='../wind_files_309x288/'
          open(1,file='fwind.dat')
          open(3,file='sfwind_jun_aug.dat')

          open(4,file='xxyy0.dat')
         do j=1,n
         read(4,100)(x(i,j),i=1,m)
         enddo
         do j=1,n
         read(4,100)(y(i,j),i=1,m)
         enddo
100      format(500f10.2)
         close(4)
          
          do j=1,n
          do i=1,m-1
           ang(i,j)=atan2(y(i+1,j)-y(i,j),x(i+1,j)-x(i,j))
          enddo
          enddo
          do j=1,n
           ang(m,j)=ang(m-1,j)
          enddo

c          open(7,file='ang.dat')
c           do j=1,n
c            write(7,100)(ang(i,j),i=1,m)
c           enddo
c          stop
 


          nhour=75*24
          do k=1,nhour
          num_file=k
          nm_first=mod(num_file/1000,10)
          nm_second=mod(num_file/100,10)
          nm_third=mod(num_file/10,10)
          nm_fourth=mod(num_file,10)
   
   
          write(file_name(1:1),'(I1)')nm_first
          write(file_name(2:2),'(I1)')nm_second
          write(file_name(3:3),'(I1)')nm_third
          write(file_name(4:4),'(I1)')nm_fourth

          write(1,111)folder,'wind',file_name,'.dat'
111       format(A22,A4,A4,A4)

c --- write wind files
c          goto 1111
          read(3,*,end=120)time,amp,angle

         do j=1,n
         do i=1,m
          wind_x(i,j)=-amp*sin(angle*3.1415926/180.+ang(i,j))
          wind_y(i,j)=-amp*cos(angle*3.1415926/180.+ang(i,j))
         enddo
         enddo

         open(2,file=folder//'wind'//file_name//'.dat')
         
         do j=1,n
          write(2,*)(wind_x(i,j),i=1,m)
         enddo
         do j=1,n
          write(2,*)(wind_y(i,j),i=1,m)
         enddo

         close(2)

1111     continue
          enddo

120      print*,k-1

          close(1) 

          end





