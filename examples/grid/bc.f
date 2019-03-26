     	subroutine bc
c
c ======================================================================
c
c   Purpose -
c     set boundary conditions
c
c   BC is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE     ACCEL   CONVECT  CONVECTC     SETUP    
c
c
c
c   BC calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
c
c
c   *	added call for wavemaker Airy, Stokes(kl=7), Piston(kl=8)
c   *	Sponge layer and SRC(kr=7)   
c   *	change  the BCs at surface and empty cells 
c
c	Last modify and tested on Oct.11, 2002 Qun Zhao
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c       
      ijl=1
      ijr=imax
      do 100 j=1,jmax
        f(ijl)=f(ijl+1)
        f(ijr)=f(ijr-1)
        p(ijl)=p(ijl+1)
        p(ijr)=p(ijr-1)
        go to (10,20,30,40,30,5,45,48), kl
   5   u(ijl)=uinf(3)
	v(ijl)=v(ijl+1)
        go to 50
   10   u(ijl)=0.0
        v(ijl)=v(ijl+1)
        go to 50
   20   u(ijl)=0.0
        v(ijl)=-v(ijl+1)*delx(1)/delx(2)
        go to 50
   30   if (iter.gt.0) go to 50
        u(ijl)=u(ijl+1)*(x(2)*rx(1)*cyl+1.0-cyl)
        v(ijl)=v(ijl+1)
        go to 50
   40   u(ijl)=u(ijr-2)
        v(ijl)=v(ijr-2)
        f(ijl)=f(ijr-2)
	goto 50
   45   continue
	if (iwave.eq.2) call airy
	if (iwave.eq.3) call stokes
	if (iwave.eq.4) call Solitary
        if (iwave.eq.6.and.ncyc.eq.1.and.j.eq.1) call Solitary_Initial
	if(iwave.eq.5.and.ncyc.ge.1) call  cnoidal	
	goto 50
   48	continue
	if (iwave.eq.2) call wavemk2 
	if (iwave.eq.3) call wavemk3  
   50   go to (60,70,80,90,80,55,92), kr
   55   u(ijr-1)=uinf(4)
        v(ijr)=v(ijr-1)
        go to 95
   60   u(ijr-1)=0.0
        v(ijr)=v(ijr-1)
        go to 95
   70   u(ijr-1)=0.0
        v(ijr)=-v(ijr-1)*delx(imax)/delx(im1)
        go to 95
   80   if (iter.gt.0) go to 95
        u(ijr-1)=u(ijr-2)*(x(ibar)*rx(im1)*cyl+1.0-cyl)
        v(ijr)=v(ijr-1)
        go to 95
   90   u(ijr-1)=u(ijl+1)
        v(ijr-1)=v(ijl+1)
        p(ijr-1)=p(ijl+1)
        f(ijr-1)=f(ijl+1)
        v(ijr)=v(ijl+2)
        f(ijr)=f(ijl+2)
	goto 95
c-------->kr=7 radiation boundary;
 92	 continue
	u(ijr)=1.0/(1./delt+wavec/delx(im1))*
	1(un(ijr)/delt+wavec/delx(im1)*u(ijr-1))
	v(ijr)=1.0/(1./delt+wavec/delx(im1))*
	1(vn(ijr)/delt+wavec/delx(im1)*v(ijr-1))
   95   ijl=ijl+imax
        ijr=ijr+imax
  100 continue
c
      ijb=1
      ijt=imax*jm1+1
      do 200 i=1,imax
        f(ijb)=f(ijb+imax)
        f(ijt)=f(ijt-imax)
        p(ijb)=p(ijb+imax)
        p(ijt)=p(ijt-imax)
        go to (110,120,130,140,130,105), kt
  105   v(ijt-imax)=vinf(2)
        u(ijt)=u(ijt-imax)
        go to 150
  110   v(ijt-imax)=0.0
        u(ijt)=u(ijt-imax)
        go to 150
  120   v(ijt-imax)=0.0
        u(ijt)=-u(ijt-imax)*dely(jmax)/dely(jm1)
        go to 150
  130   if (iter.gt.0) go to 150
        v(ijt-imax)=v(ijt-2*imax)
        u(ijt)=u(ijt-imax)
        go to 150
  140   v(ijt-imax)=v(ijb+imax)
        u(ijt-imax)=u(ijb+imax)
        p(ijt-imax)=p(ijb+imax)
        f(ijt-imax)=f(ijb+imax)
        u(ijt)=u(ijb+2*imax)
        f(ijt)=f(ijb+2*imax)
  150   go to (160,170,180,190,180,155), kb
  155   v(ijb)=vinf(1)
        u(ijb)=u(ijb+imax)
        go to 195
  160   v(ijb)=0.0
        u(ijb)=u(ijb+imax)
        go to 195
  170   v(ijb)=0.0
        u(ijb)=-u(ijb+imax)*dely(1)/dely(2)
        go to 195
  180   if (iter.gt.0) go to 195
        v(ijb)=v(ijb+imax)
        u(ijb)=u(ijb+imax)
        go to 195
  190   v(ijb)=v(ijt-2*imax)
        u(ijb)=u(ijt-2*imax)
        f(ijb)=f(ijt-2*imax)
  195   ijt=ijt+1
        ijb=ijb+1
  200 continue
c
      if (ibcflg.eq.0) go to 9999
c
c.... free surface and sloped boundary conditions
c
      do 421 i=2,im1
        do 420 j=2,jm1
c
          ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax
          imjm=ij-1-imax
          imjp=ij-1+imax
          ipjm=ij-imax+1
          ipjp=ij+imax+1
c
          if (ac(ij).gt.em6) go to 210
c
          sumac=9.*ac(ij)/16. +
     &          3.*(ac(ipj)+ac(imj)+ac(ijp)+ac(ijm))/32. +
     &          1.*(ac(ipjp)+ac(imjp)+ac(ipjm)+ac(imjm))/64.
          sumfac=9.*ac(ij)*f(ij)/16. +
     &           3.*(ac(ipj)*f(ipj)+ac(imj)*f(imj)+
     &               ac(ijp)*f(ijp)+ac(ijm)*f(ijm))/32. +
     &           1.*(ac(ipjp)*f(ipjp)+ac(imjp)*f(imjp)+
     &               ac(ipjm)*f(ipjm)+ac(imjm)*f(imjm))/64.
          f(ij)=sumac*sumfac/(sumac+1.0d-25)**2
c
  210 continue
c
c.......surface cell velocities
c 
	dyb=yj(j)-yj(j-1)
	dxr=xi(i+1)-xi(i)
	dxl=xi(i)-xi(i-1)
	if(nf(ij).eq.0.or.nf(ij).gt.5) goto 420
      nfsb=0	 
      if(f(ipj).lt.emf) nfsb=nfsb+1
      if(f(ijp).lt.emf) nfsb=nfsb+2
      if(f(imj).lt.emf) nfsb=nfsb+4
      if(f(ijm).lt.emf) nfsb=nfsb+8
      if(nfsb.eq.0) goto 420
      if(nfsb.gt.8) goto 240
      goto(250,260,270,280,290,300,310,320), nfsb
 240  nfsb1=nfsb-8
      goto(330,340,350,360,370,380,390), nfsb1
 250	u(ij)=(dyb*u(imj)+delx(i)*u(ijm))/(delx(i)+dyb)
      goto 400
 260  v(ij)=v(ijm)-dely(j)/delx(i)*(u(ij)-u(imj))
      goto 400
 270	u(ij)=(dyb*u(imj)+delx(i)*u(ijm))/(delx(i)+dyb)
      goto 260					
 280	u(imj)=(dyb*u(ij)+delx(i)*u(imjm))/(delx(i)+dyb)
      goto 400
 290   u(imj)=u(imjm)
  	u(ij)=u(ijm)
 300	u(imj)=(dyb*u(ij)+delx(i)*u(imjm))/(delx(i)+dyb)
      goto 260					
 310  u(imj)=u(imjm)
c 310  u(imj)=u(imjm1)  <---------- imjm1 = not assigned
      u(ij)=u(ijm)
      goto 260					
 320  v(ijm)=v(ij)+dely(j)/delx(i)*(u(ij)-u(imj))
      goto 400
 330	u(ij)=u(ijp)
      goto 320
 340  v(ij)=v(imj)
      goto 320
 350  v(ij)=v(imj)
      v(ijm)=v(imjm)
      goto 250
 360	u(imj)=u(imjm)
      goto 320
 370  u(ij)=u(ijp)
      u(imj)=u(imjp)
      goto 320
 380  v(ij)=v(ipj)
      v(ijm)=v(ipjm)
      goto 280
 390	u(ij)=u(ijm)
      v(ijm)=v(ij)
      v(ijp)=v(ij)
 400	continue 

	if(f(ipj).lt.emf.and.f(ipjp).lt.emf.and.f(ijp).gt.emf)
     1  v(ipj)=v(ij)-(xi(i+1)-xi(i))*(u(ijp)-u(ij))/(yj(j+1)-yj(j))
	if(f(ijp).lt.emf.and.f(ipjp).lt.emf.and.f(ipj).gt.emf)
     1  u(ijp)=u(ij)-(yj(j+1)-yj(j))*(v(ipj)-v(ij))/(xi(i+1)-xi(i))
	if(f(imj).lt.emf.and.f(imjp).lt.emf.and.f(ijp).gt.emf)
     1 v(imj)=v(ij)
     2          +(xi(i)-xi(i-1))*(u(imjp)-u(imj))/(yj(j+1)-yj(j))
	if(f(ijm).lt.emf.and.f(ipjm).lt.emf.and.f(ipj).gt.emf)
     1  u(ijm)=u(ij)
     2        +(yj(j)-yj(j-1))*(v(ipjm)-v(ijm))/(xi(i+1)-xi(i))
	 
 420	continue
 421	continue 
c
 9999 continue
c
c.... Special velocity boundary conditions
c
      return
      end



