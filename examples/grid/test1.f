c	program test1

	include "comdk1.h"

	open(1,file='f.dat')
	open(2,file='u.dat')

	imax=200
	jmax=100

	do j=1,jmax
	yj(j)=(j-1.)*1.
	dely(j)=1.
	enddo

	do i=1,imax
	xi(i)=(i-1.)*1.
	enddo

	waveh=12.5
	flht=40.
	gy=-980.

	do kkkk=10,100
	print*,kkkk
	t=(kkkk-1.)*0.01
	call solitary
	enddo

	end


