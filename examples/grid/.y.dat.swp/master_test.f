        Program master
        integer istart,iend,itime


        istart=1
        iend=3

        do itime=istart,iend
        call swan(itime)
        enddo

        end
