c=========================================================
	subroutine wavenum(pi,g,h,wavew,k0,xkn1)
c	 Qun Zhao Sep. 22,2002
c=========================================================
	parameter(m=30)
        real*8  pi,g,h,wavew,k0,xkn1(m),x(m),keci
	sigma=wavew*wavew*h/g
c
c	calculate progressive wave number---k0
c	x*tanh(x)=sigma, x=k0*h
c
	if(sigma.ge.1.0) then
	x1=sigma
	else
	x1=sqrt(sigma)
	end if
  10	aa=tanh(x1)*(x1*tanh(x1)-sigma)
	bb=tanh(x1)*tanh(x1)*(1.0-sigma)+sigma
	x2=x1-aa/bb
	error=abs(x1-x2)
	if(error.ge.0.0001) then
	x1=x2
	goto 10
	else
	k0=x2/h
	end if
c
c	calculate envancent wave numbers with  M modes
c	x*tan(x)=-sigma
c	quadatic approximation--Eq.(15),(16)   

	do 20 n=1, m
	beta=(n*pi*pi-4.*n+2.)*pi/2./(8.*n-n*pi*pi-2.)
	gama=((2.*n-1.)*pi+4.*n*beta)/beta/(2.*n-1)/pi
	gama=1./gama
	alfa=(2.*n-1.)*beta/gama
	a=alfa+sigma
	b=sigma*beta+alfa*(gama-pi/2.)
	c=-alfa*gama*pi/2.
	d=b*b-4.*a*c
	x1=(-b+sqrt(d))/2./a
	x2=(-b-sqrt(d))/2./a
	if(x1.gt.0) then
	keci=x1
	else if (x2.gt.0) then
	keci=x2
	end if
	xa=(n-0.5)*pi+keci
c
c	solve equation by Newton iteration--Eq.(17)
c
	x1=xa
 16	t=tan(x1)
	z=(1+t*t)*x1
	x2=(z*x1-sigma)/(z+t)
	error=abs(x2-x1)
	if(error.ge.0.0001) then
	x1=x2
	in=in+1
	goto 16
	else
	x(n)=x2
	xkn1(n)=x(n)/h
	end if
 20	continue

	return

	end



      subroutine airy 
c
c ======================================================================
c
c   Purpose -	 wave maker using analytical solution
c	  	 Airy wave theory. Called when iwave=2 and kl=7
c	 	 Qun Zhao Sep.23,2002
c ======================================================================
c
c##############################################################
        implicit real*8 (a-h,o-z)
c        include "32bit.h"
c##############################################################
c
c############
        include "comdk1.h" 
c############

	data nt/2/

   
	shkd=sinh(wavek*flht)
	chkd=cosh(wavek*flht)
	sh2kd=sinh(2.*wavek*flht)
	wavewt=wavew*t

	waveht=waveh
	if((t/float(nt)/wavet).le.1)
     1 waveht=waveh*(1-exp(-5.0*t/float(nt)/wavet))


	eta10=waveht/2.*sin(wavek*xi(1)-wavewt)+flht
 

	do 20 j=1,jmax

        ij1=(j-1)*imax+1

	f(ij1)=(eta10-y(j)+dely(j))/dely(j)
	if(f(ij1).le.em6) f(ij1)=0.0
	if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0

	u(ij1)=waveht/2.*wavew/shkd*cosh(wavek*yj(j))
     1	*sin(wavek*x(1)-wavewt)

	v(ij1)=waveht/2.*wavew/shkd*sinh(wavek*y(j))
     1	*cos(wavek*xi(1)-wavewt) 

 20	continue

 	return
	end



      subroutine stokes 
c
c ======================================================================
c
c   Purpose -	 wave maker using analytical solution
c	  	 Stokes wave theory. Called when iwave=3 and kl=7
c	 	 Qun Zhao Oct.13,2002 tested.
c ======================================================================
c
c##############################################################
        implicit real*8 (a-h,o-z)
c        include "32bit.h"
c##############################################################
c
c############
        include "comdk1.h" 
c############

	data nt/2/
	shkd=sinh(wavek*flht)
	chkd=cosh(wavek*flht)
	sh2kd=sinh(2.*wavek*flht)
	wavewt=wavew*t

	waveht=waveh
	if((t/float(nt)/wavet).le.1)
     1  waveht=waveh*(1-exp(-5.0*t/float(nt)/wavet))
	 
	eta10=waveht/2.*cos(wavek*xi(1)-wavewt)
     1  +waveht/8.*(pi*waveht/wavel)*chkd/shkd**3
     2  *(2+cosh(2.*wavek*flht))*cos(2.*(wavek*xi(1)-wavewt))
     3  +flht

 

c	do 20 j=1,jmax
	do 20 j=2,jmax

        ij1=(j-1)*imax+1
c	if(j.eq.1) f(ij1)=(eta10-y(j-1))/dely(j)
	f(ij1)=(eta10-y(j-1))/dely(j)
	if(f(ij1).le.em6) f(ij1)=0.0
	if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0
	 

	u(ij1)=pi*waveht/wavet*cosh(wavek*yj(j))/shkd
     1   *cos(wavek*x(1)-wavewt)
     2   +3./4.*pi*waveht/wavet*(pi*waveht/wavel)
     3	 *cosh(2.*wavek*yj(j))/shkd**4
     4   *cos(2.*(wavek*x(1)-wavewt))

	v(ij1)=pi*waveht/wavet*sinh(wavek*y(j))/shkd
     1   *sin(wavek*xi(1)-wavewt)
     2   +3./4.*pi*waveht/wavet*(pi*waveht/wavel)
     3	 *sinh(2.*wavek*y(j))/shkd**4
     4   *sin(2.*(wavek*xi(1)-wavewt))
	
 

        if(f(ij1).lt.emf) then
	u(ij1)=0.
	v(ij1)=0.
	end if
 20	continue

	return
	end



	subroutine  sponge_layer 
c
c ======================================================================
c
c   Purpose - sponge layer for open boundary called when kr=7
c             the SRC is applied at the outer edge of the sponge layer 
c	      Qun Zhao Oct. 13, 2002  
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

	xlr=xl(nkx+1)-xle; ir1=xle/delx(1);
	if(kr.eq.7) then
	do i=ir1+1,im1
	do 111 j=2,jm1
	ij=(j-1)*imax+i
	if(ac(ij).lt.emf) goto 111
	xr=xi(i)-xle
	v(ij)=v(ij)*exp(-rc*xr/xlr)
	if(f(ij).lt.emf) v(ij)=0.
111	continue
	end do
	end if

	return
	end



	 
	subroutine wavemk3 
c
c ======================================================================
c
c   Purpose - piston type abosorbing generating wavemaker, 
c             Stokes wave theory. Called when iwave=3 and kl=8
c	      Qun Zhao Oct. 12, 2002 tested
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
 
	 
	shkd=sinh(wavek*flht)
	chkd=cosh(wavek*flht)
	sh2kd=sinh(2.*wavek*flht)
	www=4.0*shkd*shkd/(2.*wavek*flht+sh2kd)

	if(t.le.wavet) then
	waveht=waveh*t/(wavet)
	wavehtn=waveh*(t-delt)/(wavet)
	else 
	waveht=waveh
	end if
		
c	umk=wavew*waveht/2./www*cos(wavek*x(2)-wavew*t)+
c 	1 waveht**2*wavew/16./flht*(3.*chkd/shkd**3-2./www)
c     2 *cos(2.*(wavek*xi(2)-wavew*t))

c	eta_stokesn=wavehtn/2.*cos(wavek*xi(2)-wavew*(t-delt))
c     1  +wavehtn/8.*(pi*wavehtn/wavel)*chkd/shkd**3*
c     2  (2+cosh(2.*wavek*flht))
c     3  *cos(2.*(wavek*xi(2)-wavew*(t-delt)))+flht


	umk=wavew*waveht/2./www*cos(wavek*x(1)-wavew*t)+
     1  waveht**2*wavew/16./flht*(3.*chkd/shkd**3-2./www)
     2  *cos(2.*(wavek*x(1)-wavew*t))

	eta_stokesn=wavehtn/2.*cos(wavek*x(1)-wavew*(t-delt))
     1  +wavehtn/8.*(pi*wavehtn/wavel)*chkd/shkd**3*
     2  (2+cosh(2.*wavek*flht))
     3  *cos(2.*(wavek*x(1)-wavew*(t-delt)))+flht

c 	deta=etahn(2)-eta_stokesn

	etahm=(etahn(1)+etahn(2))/2.
	deta=etahm-eta_stokesn 
		
	if((t-delt).le.0) then
	deta=0. 
	end if
	um1=deta*wavew/www

	if((t-delt).ge.0) then
	umk=umk-um1
	end if

	do 20 j=2,jm1
	ij1=(j-1)*imax+1
	u(ij1)=umk       
20	continue

	return 
	end







	subroutine wavemk2 
c
c ======================================================================
c
c   Purpose - piston type abosorbing generating wavemaker, 
c             Airy wave theory. Called when iwave=2 and kl=8
c	      Qun Zhao Oct. 13, 2002 tested
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
 
	shkd=sinh(wavek*flht)
	sh2kd=sinh(2.*wavek*flht)
	www=4.0*shkd*shkd/(2.*wavek*flht+sh2kd)

	cn=0.0
	etas2=0.0

	do 10 n=1,mn
	sinknd=sin(xkn1(n)*flht)
	sin2knd=sin(2.*xkn1(n)*flht)
	cn=cn+4.0*sinknd*sinknd/(2.*xkn1(n)*flht+sin2knd)
	etas2=etas2+ 4.0*sinknd*sinknd/(2.*xkn1(n)*flht+sin2knd)
     1             *exp(-xkn1(n)*x(1))
10	continue
	
	waveht=waveh
	if(t.le.wavet)          waveht=waveh*t/(wavet)
	wavehtn=waveh
	if((t-delt).le.wavet)  wavehtn=waveh*(t-delt)/(wavet)
	 	  
        eta0=waveht/2.*cos(wavek*x(1)-wavew*t) 
	umk=eta0*wavew/www

 	eta2=waveht/2.*cos(wavek*x(1)-wavew*(t-delt)) 	 
	etas2n=wavehtn/2.*sin(wavew*(t-delt))/www*etas2	


	etahm=((etahn(1)+etahn(2))/2.-flht)
	deta=etahm-eta2-etas2n
      	
	 if((t-delt).le.0)  deta=0.  
	 um1=deta*wavew/www	
	 umk=umk-um1
c
c	   remend for mass flux: Q=pi*H0/8d
c
c 	if(um.ge.0.) then
c	umk=umk*(1-pi*waveh/8./flht)
c	else if (um.lt.0.) then
c	umk=umk*(1+pi*waveh/8./flht)
c	end if

	 
	do 20 j=1,jmax
	ij1=(j-1)*imax+1
	u(ij1)=umk       
20	continue
 
	return 
	end








	subroutine cnoidal
	 
c ======================================================================
c
c   Purpose -  wave maker using analytical solution 
c             Cnoidal wave theory. Called when iwave=5 and kl=7
c	      Qun Zhao Oct. 13, 2002  tested
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
	 
	external ellf,elle
	 
	data nt/2/
	g=abs(gy)
	wavewt=wavew*t

	waveht=waveh
	if((t/float(nt)/wavet).le.1)
     1 waveht=waveh*(1-exp(-5.0*t/float(nt)/wavet))

	wahd=waveht/flht
	ursell=g*waveht*wavet**2/flht**2
c 
c	solve for kappa1
c
	call sfkappa
	K1=ellf(kappa1)
	E2=elle(kappa1)
	lmr=(1-kappa1*kappa1)/(kappa1*kappa1)
	mu=E2/(kappa1**2*K1)
c
c	calculate wave length  and celerity
c
	call wave_length
	wavec=wavel/wavet
c
c	prepare for eta,u,v
c
	call abpnm
        call waveform

	return
	end



	subroutine waveform
c ======================================================================
c
c   Purpose - calculate water surface, velocity and pressure 
c             
c	      Qun Zhao Oct. 12, 2002  
c ======================================================================

c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############

	pit=pi/4.
	t1=t/wavet
	g=abs(gy)
	emmc=1.-kappa1**2


	xt1=2.*K1*(xi(1)/wavel-t1+pit)
	call sncndn(xt1,emmc,sn,cn,dn)
	et1=a0(0)*cn**(2*0)+a0(1)*cn**(2*1)			      
     1      +a0(2)*cn**(2*2)+a0(3)*cn**(2*3)
	eta10=flht*et1+flht


	do 200 j=1,jmax

 	ij1=(j-1)*imax+1
	f(ij1)=(eta10-y(j-1))/dely(j)
	if(f(ij1).le.em6) f(ij1)=0.0
	if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0


	zdu=yj(j)/flht
	zdv=y(j)/flht

	xtu1=2.*K1*(x(1)/wavel-t/wavet+pit)
	call sncndn(xtu1,emmc,sn,cn,dn)
	u1=0.
	do 50 n=0,3
	do 50 m=0,2
	u1=u1+b0(n,m)*zdu**(2.*m)*cn**(2*n)
 50	continue
	u(ij1)=sqrt(g*flht)*u1

	xtv1=2.*K1*(xi(1)/wavel-t/wavet+pit)
	call sncndn(xtv1,emmc,sn,cn,dn)
	v1=0.	 
	do 55 n=0,3
	do 55 m=0,2
	if(n.ge.1) v1=v1+sn*dn*(4.*K1*flht/wavel)
     1	      *n/(2*m+1)*b0(n,m)*zdv**(2*m+1)*cn**(2*n-1)	
 55	continue
	v(ij1)=sqrt(g*flht)*v1
	
        if(f(ij1).gt.emf) goto 200
	u(ij1)=0.
	v(ij1)=0.

 200	continue

	return
	end




	subroutine abpnm

c ======================================================================
c
c   Purpose -prepare for surface elevation,velocity calculation 
c             
c	      Qun Zhao Oct. 12, 2002  
c ======================================================================

c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############

	e=wahd
	a0(0)=e*(lmr-mu)
	2     +e**2*(-2.*lmr+mu-2.*lmr**2+2.*lmr*mu)/4.
	3	 +e**3*(133.*lmr-16.*mu+399*lmr**2-466*lmr*mu+100*mu**2
     &     +266*lmr**3-466*lmr**2*mu+200*lmr*mu**2)/400.
	a0(1)=e
	1     +e**2*(-3./4.)
	2     +e**3*(50.-lmr-60.*mu)/80.
	a0(2)=e**2*(3./4.)+e**3*(-151.+lmr+60.*mu)/80.
	a0(3)=e**3*(101./80.)

	b0(0,0)=e*(lmr-mu)+e**2*(lmr-mu-2.*lmr**2+2.*mu**2)/4.
	3       +e**3*(-71.*lmr+47.*mu-23.*lmr**2+97.*lmr*mu-50.*mu**2
     &          +153*lmr**3-153.*lmr**2*mu-25.*lmr*mu**2+25.*mu**3)/200.
	b0(0,1)=e**2*(-3.*lmr/4.)+e**3*(6.*lmr+24.*lmr**2-21.*lmr*mu)/8.
	b0(0,2)=e**3*(3.*lmr-3.*lmr**2)/16.
	b0(1,0)=e+e**2*(1.-6.*lmr+2.*mu)/4.
	3 +e**3*(-19.-27.*lmr+10.*mu+101.*lmr**2-100.*lmr*mu+15.*mu**2)/40.
	b0(1,1)=e**2*(-3.+3.*lmr)/2.
	3       +e**3*(6.+36.*lmr-21.*mu-24.*lmr**2+21.*lmr*mu)/4.
	b0(1,2)=e**3*(6.-39.*lmr+6.*lmr**2)/16.
	b0(2,0)=e**2*(-1.)+e**3*(-2.+32.*lmr-15.*mu)/10.
	b0(2,1)=e**2*(9./4.)+e**3*(30.-120.*lmr+63.*mu)/8.
	b0(2,2)=e**3*(-45.+45.*lmr)/16.
	b0(3,0)=e**3*(6./5.)
	b0(3,1)=e**3*(-15./2.)
	b0(3,2)=e**3*(45./16.)

	return
	end


	subroutine wave_length

c ======================================================================
c
c   Purpose - calculate wave length: Eq.(2.122)
c             
c	      Qun Zhao Oct. 12, 2002  
c ======================================================================

c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############

	e=wahd
	g=abs(gy)
	wavel=wavet*sqrt(g*flht)
	1 *(1+e*(1+2.*lmr-3.*mu)/2.
	2 +e**2*(-6.-16.*lmr+5.*mu-16.*lmr**2+10.*lmr*mu+15*mu*2)/40.
	3 +e**3*(150.+1079.*lmr-203*mu+2337*lmr**2-2653*lmr*mu+350.*mu**2
     4 +1558*lmr**3-2653*lmr**2*mu+700.*lmr*mu**2+175.*mu**3)/2800.)
10	continue
	return
	end


	Function func_2120(xkcn)
c ======================================================================
c
c   Purpose - function of modulus kappa1: Eq.(2.120)
c             
c	      Qun Zhao Oct. 12, 2002  
c ======================================================================

c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############
     
	external ellf,elle
	e=wahd
	xK1=ellf(xkcn)
	xE2=elle(xkcn)
	xlmr=(1-xkcn*xkcn)/(xkcn*xkcn)
	xmucn=xE2/(xkcn**2*xK1)

	func_2120=(16.*xkcn**2*xK1**2/3.)/(1.+e*(-1.-2.*xlmr)/4.        
     2 +e**2*(8.+33.*xlmr-10.*xmucn+33.*xlmr**2-20.*xlmr*xmucn)/40.)
     3 -ursell
	return
	end








c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c	Adapted subroutines from Numerical Recipies


	subroutine sfkappa
c ======================================================================
c
c   Purpose - solve for mudulus kappa1 by bisection method
c             	   
c ======================================================================

c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############

      PARAMETER (jmaxk=50,xacc=1.e-8,x1cn=0.0001,x2cn=0.9999)
	external func_2120,ellf,elle

      fmid=func_2120(x2cn)
      f1=func_2120(x1cn)
      if(f1*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f1.lt.0.)then
        rtbis=x1cn
        dx=x2cn-x1cn
      else
        rtbis=x2cn
        dx=x1cn-x2cn
      endif
      do 11 j=1,JMAXk
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func_2120(xmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) goto 12
11    continue
      pause 'too many bisections in rtbis'
12	continue
	kappa1=rtbis
	return
      END


      SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
c ======================================================================
c
c   Purpose - calculate Jacobi function
c             	   
c ======================================================================
      implicit real*8 (a-h,o-z)

      REAL*8 cn,dn,emmc,sn,uu,CA
      PARAMETER (CA=.0003)
      INTEGER i,ii,l
      REAL*8 a,b,c,d,emc,u,em(13),en(13)
      LOGICAL bo
      emc=emmc
      u=uu
      if(emc.ne.0.)then
        bo=(emc.lt.0.)
        if(bo)then
          d=1.-emc
          emc=-emc/d
          d=sqrt(d)
          u=d*u
        endif
        a=1.
        dn=1.
        do 11 i=1,13
          l=i
          em(i)=a
          emc=sqrt(emc)
          en(i)=emc

          c=0.5*(a+emc)
          if(abs(a-emc).le.CA*a)goto 1
          emc=a*emc
          a=c
11      continue
1       u=c*u
        sn=sin(u)
        cn=cos(u)
        if(sn.eq.0.)goto 2
        a=cn/sn
        c=a*c
        do 12 ii=l,1,-1
          b=em(ii)
          a=c*a
          c=dn*c
          dn=(en(ii)+a)/(b+a)
          a=c/b
12      continue
        a=1./sqrt(c**2+1.)
        if(sn.lt.0.)then
          sn=-a
        else
          sn=a
        endif
        cn=c*sn
2       if(bo)then

          a=dn
          dn=cn
          cn=a
          sn=sn/d
        endif
      else
        cn=1./cosh(u)
        dn=cn
        sn=tanh(u)
      endif
      return
      END


	FUNCTION ellf(ak)
c ======================================================================
c
c   Purpose - first kind complete elliptic integral
c             	   
c ======================================================================
      implicit real*8 (a-h,o-z)

      parameter(pi=3.1415927,phi=0.5*pi)
      REAL*8  ellf,ak 
C     USES rf
      REAL*8 s,rf
      s=sin(phi)
      ellf=s*rf(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.d0)
      return
      END


      FUNCTION elle(ak)
c ======================================================================
c
c   Purpose - second kind complete elliptic integral
c             	   
c ======================================================================
      implicit real*8 (a-h,o-z)

	parameter(pi=3.1415927,phi=0.5*pi)
      REAL*8 elle,ak 
C     USE rd,rf
      REAL*8 cc,q,s,rd,rf
      s=sin(phi)
      cc=cos(phi)**2
      q=(1.-s*ak)*(1.+s*ak)
      elle=s*(rf(cc,q,1.d0)-((s*ak)**2)*rd(cc,q,1.d0)/3.)
      return
      END


      FUNCTION rf(x,y,z)
c ======================================================================
c
c   Purpose - Carlson's  elliptic integral  of first kind
c             	   
c ======================================================================
      implicit real*8 (a-h,o-z)

      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,
     *C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)

        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
	end



      FUNCTION rd(x,y,z)

c ======================================================================
c
c   Purpose - Carlson's  elliptic integral  of second kind
c             	   
c ======================================================================

       implicit real*8 (a-h,o-z) 

      REAL*8 rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,
     *C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
      REAL*8 alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     *sqrtz,sum,xt,yt,zt
      if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)

        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=.2*(xt+yt+3.*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*
     *ec+delz*C4*ea)))/(ave*sqrt(ave))

      return

      END 

c ======================================================================
      subroutine Solitary
c   Purpose -
c     Solitary 3rd order wave theory by Grimshaw(1971).
c     Called when iwave=4 and kl=7
c ======================================================================
      implicit real*8 (a-h,o-z)
c      common /blk1/ waveh,gy,flht
      include "comdk1.h" 

      waveht=waveh/2.
      g=abs(gy)
      rootgh=sqrt(g*flht)
      es=waveht/flht
      es2=es**2
      es3=es**3
      wvcel=rootgh*sqrt(1.+es-0.05*es2-3./70.*es3)
      wvalpa=sqrt(0.75*es)*(1-5./8.*es+71./128.*es2);

      arg_f=wvalpa*(xi(1)-wvcel*(t-wavet))/flht          ! <------ 
      s_f=1./cosh(arg_f)
      q_f=tanh(arg_f)
      s2_f=s_f**2
      s4_f=s_f**4
      s6_f=s_f**6
      q2_f=q_f**2

      arg_u=wvalpa*(xi(1)-2.*(x(1)-xi(1))-wvcel*(t-wavet))/flht          ! <------ 
      s_u=1./cosh(arg_u)
      q_u=tanh(arg_u)
      s2_u=s_u**2
      s4_u=s_u**4
      s6_u=s_u**6
      q2_u=q_u**2

  
      eta10=2.*flht*(es*s2_f - 0.75*es2*s2_f*q2_f + es3*
     &            (5./8.*s2_f*q2_f - 101./80.*s4_f*q2_f))+flht


      do 20 j=1,jmax
        ij1=(j-1)*imax+1
        f(ij1)=(eta10-y(j))/dely(j)
        if(f(ij1).le.em6) f(ij1)=0.0
        if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0

        yoh=yj(j)/flht
        yoh_v=y(j)/flht
        z2=yoh**2
        z4=yoh**4
        z2_v=yoh_v**2
        z4_v=yoh_v**2

        u(ij1)=rootgh*(es*s2_u - es2*(-0.25*s2_u+s4_u
     &           +z2*(1.5*s2_u-2.25*s4_u)) -
     &     es3*(19./40.*s2_u + 0.2*s4_u - 1.2*s6_u
     &           +z2*(-1.5*s2_u - 3.75*s4_u + 
     &          7.5*s6_u) + z4*(-3./8.*s2_u + 45./16.*s4_u
     &           - 45./16.*s6_u)))
        u(ij1)=2.0*u(ij1)
        v(ij1)=
     &rootgh*(sqrt(3.*es)*yoh_v*q_f*(es*s2_f - 
     &     es2*(3./8.*s2_f + 2.*s4_f + z2_v*(0.5*s2_f-1.5*s4_f)) -
     &     es3*(49./640.*s2_f - 17./20.*s4_f - 18./5.*s6_f+z2_v*(-13./
     &          16.*s2_f - 25./16.*s4_f + 7.5*s6_f) 
     &           + z4_v*(-3./40.*s2_f + 
     &          9./8.*s4_f - 27./16.*s6_f))))
         v(ij1)=2.0*v(ij1)

        if(f(ij1).lt.emf) then
          u(ij1)=0.
          v(ij1)=0.
        end if
 20	  continue


      return
      end

c ======================================================================
      subroutine Solitary_Initial
c   Purpose -
c     Solitary 3rd order wave theory by Grimshaw(1971).
c     initial case for big wave
c     Called when iwave=6 and kl=7 (iwave=4 for maker)
c ======================================================================
      implicit real*8 (a-h,o-z)
c      common /blk1/ waveh,gy,flht
      include "comdk1.h" 

      waveht=waveh/2.
      g=abs(gy)
      rootgh=sqrt(g*flht)
      es=waveht/flht
      es2=es**2
      es3=es**3
      wvcel=rootgh*sqrt(1.+es-0.05*es2-3./70.*es3)
      wvalpa=sqrt(0.75*es)*(1-5./8.*es+71./128.*es2);

      es_u=waveht/flht
      es2_u=es_u**2
      es3_u=es_u**3
      wvcel_u=rootgh*sqrt(1.+es_u-0.05*es2_u-3./70.*es3_u)
      wvalpa_u=sqrt(0.75*es_u)*(1-5./8.*es_u+71./128.*es2_u);


      do i=1,imax
      arg_f=wvalpa*(xi(i)-wvcel*(wavet))/flht          ! <------ 
      s_f=1./cosh(arg_f)
      q_f=tanh(arg_f)
      s2_f=s_f**2
      s4_f=s_f**4
      s6_f=s_f**6
      q2_f=q_f**2

      arg_u=wvalpa_u*(x(i)-wvcel_u*(wavet))/flht          ! <------ 
      s_u=1./cosh(arg_u)
      q_u=tanh(arg_u)
      s2_u=s_u**2
      s4_u=s_u**4
      s6_u=s_u**6
      q2_u=q_u**2

  
      eta10=2.*flht*(es*s2_f - 0.75*es2*s2_f*q2_f + es3*
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

        u(ij1)=rootgh*(es_u*s2_u - es2_u*(-0.25*s2_u+s4_u
     &           +z2*(1.5*s2_u-2.25*s4_u)) -
     &     es3_u*(19./40.*s2_u + 0.2*s4_u - 1.2*s6_u
     &           +z2*(-1.5*s2_u - 3.75*s4_u + 
     &          7.5*s6_u) + z4*(-3./8.*s2_u + 45./16.*s4_u
     &           - 45./16.*s6_u)))
        u(ij1)=u(ij1)*2.

c -- someting wrong here compared with Grimshaw 1971(shi, 06/10/04)
        v(ij1)=
     &rootgh*(sqrt(3.*es_u)*yoh_v*q_f*(es_u*s2_f - 
     &     es2_u*(3./8.*s2_f + 2.*s4_f + z2_v*(0.5*s2_f-1.5*s4_f)) -
     &   es3_u*(49./640.*s2_f - 17./20.*s4_f - 18./5.*s6_f+z2_v*(-13./
     &          16.*s2_f - 25./16.*s4_f + 7.5*s6_f) 
     &           + z4_v*(-3./40.*s2_f + 
     &          9./8.*s4_f - 27./16.*s6_f))))
        v(ij1)=v(ij1)*2.

         p(ij1)=(eta10-yj(j))*g*rhof

        if(f(ij1).le.em6) then
          u(ij1)=0.
          v(ij1)=0.
          p(ij1)=0.
        end if
 20	  continue
        enddo

      return
      end

c ======================================================================
      subroutine Solitary_Initial_Wei
c   Purpose -
c     Solitary from Wei Boussinesq.
c     initial case for big wave
c     Called when iwave=6 and kl=7 (iwave=4 for maker)
c ======================================================================
      implicit real*8 (a-h,o-z)
c      common /blk1/ waveh,gy,flht
c       parameter(mbar=1200,nbar=192)
      include "comdk1.h" 
        real*8 ui_so(ibar2,jbar2), vi_so(ibar2,jbar2), eti_so(ibar2)       

       call solitary_wei(ui_so,vi_so,eti_so)

      wave_center=4.5*50.

      waveht=waveh/2.
      g=abs(gy)
      rootgh=sqrt(g*flht)
      es=waveht/flht
      es2=es**2
      es3=es**3
      wvcel=rootgh*sqrt(1.+es-0.05*es2-3./70.*es3)
      wvalpa=sqrt(0.75*es)*(1-5./8.*es+71./128.*es2);

      es_u=waveht/flht
      es2_u=es_u**2
      es3_u=es_u**3
      wvcel_u=rootgh*sqrt(1.+es_u-0.05*es2_u-3./70.*es3_u)
      wvalpa_u=sqrt(0.75*es_u)*(1-5./8.*es_u+71./128.*es2_u);


      do i=1,imax
      arg_f=wvalpa*(xi(i)-wave_center)/flht          ! <------ 
      s_f=1./cosh(arg_f)
      q_f=tanh(arg_f)
      s2_f=s_f**2
      s4_f=s_f**4
      s6_f=s_f**6
      q2_f=q_f**2

      arg_u=wvalpa_u*(x(i)-wave_center)/flht          ! <------ 
      s_u=1./cosh(arg_u)
      q_u=tanh(arg_u)
      s2_u=s_u**2
      s4_u=s_u**4
      s6_u=s_u**6
      q2_u=q_u**2

  
c      eta10=2.*flht*(es*s2_f - 0.75*es2*s2_f*q2_f + es3*
c     &            (5./8.*s2_f*q2_f - 101./80.*s4_f*q2_f))+flht

         eta10=eti_so(i)*100.

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

        u(ij1)=rootgh*(es_u*s2_u - es2_u*(-0.25*s2_u+s4_u
     &           +z2*(1.5*s2_u-2.25*s4_u)) -
     &     es3_u*(19./40.*s2_u + 0.2*s4_u - 1.2*s6_u
     &           +z2*(-1.5*s2_u - 3.75*s4_u + 
     &          7.5*s6_u) + z4*(-3./8.*s2_u + 45./16.*s4_u
     &           - 45./16.*s6_u)))
        u(ij1)=u(ij1)*1.6

c -- someting wrong here compared with Grimshaw 1971(shi, 06/10/04)
        v(ij1)=
     &rootgh*(sqrt(3.*es_u)*yoh_v*q_f*(es_u*s2_f - 
     &     es2_u*(3./8.*s2_f + 2.*s4_f + z2_v*(0.5*s2_f-1.5*s4_f)) -
     &   es3_u*(49./640.*s2_f - 17./20.*s4_f - 18./5.*s6_f+z2_v*(-13./
     &          16.*s2_f - 25./16.*s4_f + 7.5*s6_f) 
     &           + z4_v*(-3./40.*s2_f + 
     &          9./8.*s4_f - 27./16.*s6_f))))
        v(ij1)=v(ij1)*1.6

         p(ij1)=(eta10-yj(j))*g*rhof

        if(f(ij1).le.em6) then
          u(ij1)=0.
          v(ij1)=0.
          p(ij1)=0.
        end if
 20	  continue
        enddo

c          do i=1,imax
c          eta10=eti_so(i)*100.

c      do 20 j=1,jmax
c        ij1=(j-1)*imax+i
c        f(ij1)=(eta10-y(j))/dely(j)
c        if(f(ij1).le.em6) f(ij1)=0.0
c        if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0
c
c         u(ij1)=ui_so(i,j)*100.
c         v(ij1)=vi_so(i,j)*100.
c
c         p(ij1)=(eta10-yj(j))*g*rhof
c
c        if(f(ij1).le.em6) then
c          u(ij1)=0.
c          v(ij1)=0.
c          p(ij1)=0.
c        end if
c 20	  continue
c        enddo


        open(1,file='zeta.out')
        do i=1,imax
        write(1,*)i,eti_so(i)
        enddo
        close(1)

        open(1,file='u.out')
        do j=1,jmax
        write(1,100)(u(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        enddo
        close(1)

100 	format(800f16.8)

       return
       end

        subroutine solitary_wei(uii,vii,etii)
        implicit real*8 (a-h,o-z)
        include 'comdk1.h'
c        parameter(mbar=1200,nbar=192) 

        real*8 uii(ibar2,jbar2), vii(ibar2,jbar2), etii(ibar2)
        real*8 ua(ibar2)
        integer ngage,ixg(30),iyg(30),itg(6),ispg(4)
        character*30  f1n, f2n, f3n, f4n, f5n, f6n, f7n

        namelist/data1/ibe,imch,a0,h0,tpd,dx,dy,dt,mx,ny,nt
     &                 ,itbgn,itend,itdel,itscr,itftr, theta
     +                 ,cbkv,delta,slmda,isltb,islte
        namelist/data2/isrc,jsrc,cspg,cspg2,cspg3,ispg,ngage
     &                ,ixg,iyg,itg,cbrk,ck_bt,c_dm
     +                ,cbrk,ck_bt,c_dm,isld,idout,idft

        open (1, file='funwave2d.data')

        read (1, nml = data1)
        read (1, nml = data2)
        close(1)

        g=9.8
        alpha=-0.39

c----------------------------------------------------------------------
c====> Case 3: Solitary Wave   the peak at t=0 is at i=isrc
c----------------------------------------------------------------------

       print *, 'isrc=', isrc

          call sub_sltry(a0, h0, alpha, cph, bue, ae1, ae2, au, C_ph)
          write (*,*) 'cph,bue,ae1,ae2,au = ',cph,bue,ae1,ae2,au

            kk=1.

             do i = 1, imax
                xi_z        = bue*(float(i-isrc)*dx-C_ph*(kk-1.)*dt)
                stmp_z      = 1./(cosh(xi_z))**2
                xi_u   = bue*(float(i-isrc)*dx-dx/2.-C_ph*(kk-1.)*dt)
                stmp_u      = 1./(cosh(xi_u))**2
                etii(i)  = (ae1+ae2*stmp_z)*stmp_z+h0
                ua(i)   = au*stmp_u
             enddo

             do i=2,imax-1
                do j = 1, jmax
                z_so=(j-1.)*dy
                za_so=0.4690*h0  
                uii(i,j)=ua(i)
     &        +1./dx/dx*(ua(i+1)-2.*ua(i)+ua(i-1))
     &        *(0.5*(za_so*za_so-z_so*z_so)+1.*(1.-h0)*(za_so-z_so))
                vii(i,j)=0.
                enddo        
             enddo
             do j=1,jmax
              uii(1,j)=uii(2,j)
              uii(imax,j)=uii(imax-1,j)
              vii(1,j)=vii(2,j)
              vii(imax,j)=vii(imax-1,j)
             enddo
1100       continue

         return
         stop
         end


c-----------------------------------------------------------------------
c        subroutine for obtain solitary wave solution for NBM
c-----------------------------------------------------------------------

         subroutine sub_sltry (a0, h1, alp, r1, bue, ae1, ae2, au, C_ph)

         implicit real*8 (a-h,o-z)
         implicit integer (i-n)
          g  = 9.81

c
c--------coefficients for third order polynomial equations
c               " x**3+p*x**2+q*x+r=0  "

         alp2 = alp + 1./3.
         eps  = a0/h1

         p = -(alp2+2.*alp*(1.0+eps))/(2.0*alp)
         q = eps*alp2/alp
         r = alp2/(2.*alp)

c--------Newton-Rapson's method to solve x ( >1 )for the above equation
         ite = 0
         x=1.2
1        fx  = r+x*(q+x*(p+x))
         fpx = q+x*(2.*p+3.*x)
          x = x-fx/fpx

         ite = ite+1
         if (ite.gt.10) then
            write(*,*) 'no solitary wave solution (check eps = a0/h1)'
            stop
         endif
         if (abs(fx).ge.1e-5) goto 1
         rx = sqrt(x)
         cph = sqrt(g*h1)
         r1 = rx*cph
         C_ph=rx*cph

         write (*,*) rx, r1


c---------coefficients for solitary solutions :
c         "   u =  au/(cosh(bue*xi))**2    "
c         "  et = ae1/(cosh(bue*xi)**2+ae2/(cosh(bue*xi)**4   "

          au =  (x-1.)/(eps*rx)*cph*eps
         bue =  sqrt((x-1.)/(4.*(alp2-alp*x)))/h1
         ae1 =  (x-1.)/(eps*3.*(alp2-alp*x))*a0
         ae2 = -(x-1.)/(2.*eps)*(x-1.)*(2.*alp*x+alp2)
     &                                 /(x*(alp2-alp*x))*a0
         zr  = bue*rx*cph


         return
         end
