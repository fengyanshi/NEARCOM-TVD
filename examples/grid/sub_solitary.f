	program test1

	include "comdk1.h"

c	open(1,file='f.dat')
c	open(2,file='u.dat')


	imax=400
	jmax=80

	do j=1,jmax
	yj(j)=0.5+(j-1.)*1.
        y(j)=1.+(j-1.)*1.
	dely(j)=1.
	enddo

	do i=1,imax
	xi(i)=0.5+(i-1.)*1.
        x(i)=1.+(i-1.)*1.
	enddo

	waveh=21.9
	flht=30.
	gy=-980.
        wavet=0.8

c	do kk=1,1000
c	print*,kk
c	t=(kk-1.)*0.01
	call solitary
c	enddo

	end


c ======================================================================
      subroutine Solitary
c   Purpose -
c     Solitary 3rd order wave theory by Grimshaw(1971).
c     Called when iwave=4 and kl=7
c ======================================================================
c      implicit real*8 (a-h,o-z)
c      common /blk1/ waveh,gy,flht
      include "comdk1.h" 

      waveht=waveh
      g=abs(gy)
      rootgh=sqrt(g*flht)
      es=waveht/flht
      es2=es**2
      es3=es**3
      wvcel=rootgh*sqrt(1.+es-0.05*es2-3./70.*es3)
      wvalpa=sqrt(0.75*es)*(1-5./8.*es+71./128.*es2);

      do i=1,imax
      arg_f=wvalpa*(xi(i)-wvcel*(wavet))/flht          ! <------ 
      s_f=1./cosh(arg_f)
      q_f=tanh(arg_f)
      s2_f=s_f**2
      s4_f=s_f**4
      s6_f=s_f**6
      q2_f=q_f**2

      arg_u=wvalpa*(x(i)-wvcel*(wavet))/flht          ! <------ 
      s_u=1./cosh(arg_u)
      q_u=tanh(arg_u)
      s2_u=s_u**2
      s4_u=s_u**4
      s6_u=s_u**6
      q2_u=q_u**2

  
      eta10=flht*(es*s2_f - 0.75*es2*s2_f*q2_f + es3*
     &            (5./8.*s2_f*q2_f - 101./80.*s4_f*q2_f))+flht

      do 20 j=1,jmax
        ij1=(j-1)*imax+i
        f(ij1)=(eta10-y(j))/dely(j)
        if(f(ij1).le.em6) f(ij1)=0.0
        if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0

        yoh=yj(j)/flht
        yoh_v=y(j)/flht
        z2=yoh**2
        z4=yoh**4
        z2_v=yoh_v**2
        z4_v=yoh_v**2

c -- someting wrong here compared with Grimshaw 1971(shi, 06/10/04)

c        u(ij1)=rootgh*(es*s2_u - es2*(-0.25*s2_u+s4_u
c     &           +z2*(1.5*s2_u-2.25*s4_u)) -
c     &     es3*(19./40.*s2_u + 0.2*s4_u - 1.2*s6_u
c     &           +z2*(-1.5*s2_u - 3.75*s4_u + 
c     &          7.5*s6_u) + z4*(-3./8.*s2_u + 45./16.*s4_u
c     &           - 45./16.*s6_u)))

        u(ij1)=rootgh*(es*s2_u + es2*(-0.75*s2_u+s2_u*q2_u
     &           +z2*(3./4.*s2_u-2.25*s2_u*q2_u)) +
     &     es3*(21./40.*s2_u - s2_u*q2_u - 1.2*s4_u*q2_u
     &           +z2*(-9./4.*s2_u + 3.75*s2_u*q2_u + 
     &          7.5*s4_u*q2_u) + z4*(3./8.*s2_u - 45./16.*s4_u
     &          *q2_u)))

c -- someting wrong here compared with Grimshaw 1971(shi, 06/10/04)
c        v(ij1)=
c     &rootgh*(sqrt(3.*es)*yoh_v*q_f*(es*s2_f - 
c     &     es2*(3./8.*s2_f + 2.*s4_f + z2_v*(0.5*s2_f-1.5*s4_f)) -
c     &     es3*(49./640.*s2_f - 17./20.*s4_f - 18./5.*s6_f+z2_v*(-13./
c     &          16.*s2_f - 25./16.*s4_f + 7.5*s6_f) 
c     &           + z4_v*(-3./40.*s2_f + 
c     &          9./8.*s4_f - 27./16.*s6_f))))

        v(ij1)=
     &rootgh*(sqrt(3.*es)*yoh_v*q_f*(es*s2_f - 
     &     es2*(3./8.*s2_f + 2.*s4_f + z2_v*(0.5*s2_f-1.5*s4_f)) -
     &     es3*(49./640.*s2_f + 17./20.*s4_f + 18./5.*s6_f-z2_v*(-13./
     &          16.*s2_f - 25./16.*s4_f + 7.5*s6_f) 
     &           + z4_v*(-3./40.*s2_f + 
     &          9./8.*s4_f - 27./16.*s6_f))))

        if(f(ij1).le.em6) then
          u(ij1)=0.
          v(ij1)=0.
        end if
 20	  continue
        enddo

	open(1,file='f.dat')
        do j=1,jmax
        write(1,201)(f((j-1)*imax+i),i=1,imax)
        enddo
        close(1)

	open(1,file='u.dat')
        do j=1,jmax
        write(1,201)(u((j-1)*imax+i),i=1,imax)
        enddo
        close(1)

	open(1,file='v.dat')
        do j=1,jmax
        write(1,201)(v((j-1)*imax+i),i=1,imax)
        enddo
        close(1)


201     format(400f16.8)


      return
      end
