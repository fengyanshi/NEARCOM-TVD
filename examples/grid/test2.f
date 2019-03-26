

c ======================================================================
      subroutine Solitary
c   Purpose -
c     Solitary 3rd order wave theory by Grimshaw(1971).
c     Called when iwave=4 and kl=7
c ======================================================================
c      implicit real*8 (a-h,o-z)
c      common /blk1/ waveh,gy,flht
      include "comdk1.h" 

	print*,'test---',waveh,flht,gy

      waveht=waveh
      g=abs(gy)
      rootgh=sqrt(g*flht)
      es=waveht/flht
      es2=es**2
      es3=es**3
      wvcel=rootgh*sqrt(1.+es-0.05*es2-3./70.*es3)
      wvalpa=sqrt(0.75*es)*(1-5./8.*es+71./128.*es2);

      arg=wvalpa*(xi(1)-wvcel*t)/flht          ! <------ 
      s=1./cosh(arg)
      q=tanh(arg)
      s2=s**2
      s4=s**4
      s6=s**6
      q2=q**2
  
      eta10=flht*(es*s2 - 0.75*es2*s2*q2 + es3*
     &            (5./8.*s2*q2 - 101./80.*s4*q2))
	print*,t,flht,eta10

      do 20 j=2,jmax
        ij1=(j-1)*imax+1
        f(ij1)=(eta10-y(j-1))/dely(j)
        if(f(ij1).le.em6) f(ij1)=0.0
        if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0

        yoh=yj(j)/flht
        z2=yoh**2
        z4=yoh**4
        u(ij1)=rootgh*(es*s2 - es2*(-0.25*s2+s4+z2*(1.5*s2-2.25*s4)) -
     &     es3*(19./40.*s2 + 0.2*s4 - 1.2*s6+z2*(-1.5*s2 - 3.75*s4 + 
     &          7.5*s6) + z4*(-3./8.*s2 + 45./16.*s4 - 45./16.*s6)))
        v(ij1)=rootgh*(sqrt(3.*es)*yoh*q*(es*s2 - 
     &     es2*(3./8.*s2 + 2.*s4 + z2*(0.5*s2-1.5*s4)) -
     &     es3*(49./640.*s2 - 17./20.*s4 - 18./5.*s6 + z2*(-13./
     &          16.*s2 - 25./16.*s4 + 7.5*s6) + z4*(-3./40.*s2 + 
     &          9./8.*s4 - 27./16.*s6))))

        if(f(ij1).lt.emf) then
          u(ij1)=0.
          v(ij1)=0.
        end if
 20	  continue

	do j=2,jmax
        ij1=(j-1)*imax+1
	write(1,*)f(ij1)
	write(2,*)u(ij1),v(ij1)
	enddo

	pause

      return
      end
