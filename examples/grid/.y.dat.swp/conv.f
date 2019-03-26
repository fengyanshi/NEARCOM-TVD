         real x
         character*4 time
         character*1000 data

         open(1,file='outputuv.not')
         open(2,file='outputuv.out')

         icount=0
          do kk=1,9999999
          read(1,110,end=200)time          
         if (time.eq.'Time')then
           print*,'Time'
           icount=icount+1
            if(icount.eq.2)then
           print*,'writing'
              do k=1,4
              read(1,120)data
              write(2,120)data
              icount=0
              enddo
             endif
          endif
         enddo
          
200      continue
110      format(A5)
120      format(A1000)
         end
