
       subroutine stress_vt 
c
c ======================================================================
c
c   Purpose -	Compute turbulence shear stress :  visx_vt,visy_vt 
c                
c	 	 Qun Zhao Nov.8,2002
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
	
	common/stress_matrix/sxx(nxy),sxy(nxy),syy(nxy) 
	real*8 sxx,sxy,syy 

c
c.... get the strain and rotation tensors
c
c     -----------------------------------------------------------------
       call strain(ar,at,ac,x,y,ri,r,un,vn,div,curlu,im1,jm1,nxy,
     &            sxx,syy,sxy)     
c     -----------------------------------------------------------------
c
c.... bottom and top boundaries for sxx, syy, sxy only
c
      ijb=2
      ijt=jm1*imax+2
      deltyb=dely(2)
      deltyt=dely(jm1)
      do 110 i=2,im1
        sxyb=un(ijb-1+imax)/deltyb
        sxyt=-un(ijt-1-imax)/deltyt
        sxx(ijt)=sxx(ijt-imax)
        syy(ijt)=syy(ijt-imax)
        sxy(ijt)=cvmgt(sxyt,zero,kt.eq.2)        
        sxx(ijb)=sxx(ijb+imax)
        syy(ijb)=syy(ijb+imax)
        sxy(ijb+imax)=cvmgt(sxyb,zero,kb.eq.2)        
        ijb=ijb+1
        ijt=ijt+1
  110 continue
c
c.... right and left boundaries for sxx, syy, sxy only
c
      ijl=1
      ijr=imax
      deltxl=delx(2)
      deltxr=delx(im1)
      do 120 j=1,jmax
        if (j.eq.1) then
          ijvl=ijl+1
          ijvr=ijr-1
        else
          ijvl=ijl+1-imax
          ijvr=ijr-1-imax
        endif
        sxyl=vn(ijvl)/deltxl
        sxyr=-vn(ijvr)/deltxr
        sxx(ijr)=sxx(ijr-1)
        syy(ijr)=syy(ijr-1)
        sxy(ijr)=cvmgt(sxyr,zero,kr.eq.2)        
        sxx(ijl)=sxx(ijl+1)
        syy(ijl)=syy(ijl+1)
        sxy(ijl+1)=cvmgt(sxyl,zero,kl.eq.2)        
        ijl=ijl+imax
        ijr=ijr+imax
  120 continue

        do   i=2,im1
        do   j=2,jm1

        ij=(j-1)*imax+i

 	  if(f(ij).lt.emf.or.ac(ij).lt.emf) then
          sxx(ij)=0.d0
          syy(ij)=0.d0
          sxy(ij)=0.d0
	  end if
	
          xnuttij=xnutt(1,ij)     	         
          sxx(ij)=2.*xnuttij*sxx(ij)
          sxy(ij)=2.*xnuttij*sxy(ij)
          syy(ij)=2.*xnuttij*syy(ij)
  
	  end do
	  end do


        do  200 i=2,im1
        do  200  j=2,jm1

        ij=(j-1)*imax+i
          ijm=ij-imax
          ijp=ij+imax
          imj=ij-1
          ipj=ij+1
          ipjp=ipj+imax

          tauxyy=(sxy(ipjp)-sxy(ipj))/dely(j)
          tauxxx=(ri(i+1)*sxx(ipj)-ri(i)*sxx(ij))
     &                /(r(i)*0.5*(delx(i+1)+delx(i)))
          visx_vt(ij)=(tauxxx+tauxyy) 


          tauxyx=(r(i)*sxy(ipjp)-r(i-1)*sxy(ijp))/(ri(i)*delx(i))
          tauyyy=2.*(syy(ijp)-syy(ij))/(dely(j+1)+dely(j))
          visy=(tauyyy+tauxyx)


          visy_vt(ij)=(tauyyy+tauxyx) 
	    

 200	continue

	return
	end



	subroutine input_klmodel(m)
c=======================================================================
c   Purpose - model closures
c  	 Qun Zhao Nov. 15, 2002 
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

	
  	do 100 i=2,im1
	do 100 j=2, jm1	
        ij=(j-1)*imax+i

	prod(m,ij)=0.
	eps(m,ij)=0.
	xnutt(m,ij)=0.				    


        if(f(ij).ge.emf.and.ac(ij).ge.emf) then

	ken(m,ij)=dmax1(ken(m,ij),em6)
	xnutm=c_s*ken(m,ij)**0.5*length(m+1,ij)	

	do k=m+1,max
	xnutt(m,ij)=xnutm+xnutt(k,ij)
	end do
	xnutt(m,ij)=dmax1(xnutt(m,ij),em6)

	eps(m,ij)=c_d*ken(m,ij)**1.5/length(m,ij)

        if(m.eq.1) prod(m,ij)=2.*xnutt(m,ij)*sij(ij)
 
 	if(m.ge.2) prod(m,ij)=eps(m-1,ij)
	if(m.eq.max) eps(m,ij)=prod(m,ij)		  
 	prod(m,ij)=dmax1(prod(m,ij),em6)
 	eps(m,ij)=dmax1(eps(m,ij),em6)

	end if
 100  continue

	 
c
c       boundary condition for prod,eps and xnut
c
        do 200 i=2,im1
        do 200 j=2,jm1

          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax
c
c       free surface no gradient
c
 	
        nff=nf(ij)
        nfr=nf(ipj)
        nfl=nf(imj)
        nft=nf(ijp)
        nfb=nf(ijm)
 
        if(nff.ge.5) then
        if(nft.lt.5) eps(m,ij)=eps(m,ijp)
        if(nfb.lt.5) eps(m,ij)=eps(m,ijm)
        if(nfr.lt.5) eps(m,ij)=eps(m,ipj)
        if(nfl.lt.5) eps(m,ij)=eps(m,imj)

        if(nft.lt.5) prod(m,ij)=prod(m,ijp)
        if(nfb.lt.5) prod(m,ij)=prod(m,ijm)
        if(nfr.lt.5) prod(m,ij)=prod(m,ipj)
        if(nfl.lt.5) prod(m,ij)=prod(m,imj)

        if(nft.lt.5) xnutt(m,ij)=xnutt(m,ijp)
        if(nfb.lt.5) xnutt(m,ij)=xnutt(m,ijm)
        if(nfr.lt.5) xnutt(m,ij)=xnutt(m,ipj)
        if(nfl.lt.5) xnutt(m,ij)=xnutt(m,imj)
	  end if
c
c       wall boundary
c
  
 200  continue


 999	continue

        return
        end



	subroutine klmodel_m(m)
c=======================================================================
c   Purpose - transport equation  for k -upwind  scheme for now
c  	 Qun Zhao Nov. 15, 2002 
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

	do 10 i=2, im1
	do 10 j=2, jm1

          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax

	uij=0.
	vij=0.

	if(f(ij).gt.emf.and.ac(ij).gt.emf) then
	call dkcal(m,ij)

	uij=0.5*(un(ij)+un(imj))
	if(uij.ge.0) fkx=dkl*uij
	if(uij.lt.0) fkx=dkr*uij
 
	vij=0.5*(vn(ij)+vn(ijm))
 	if(vij.ge.0) fky=dkb*vij
 	if(vij.lt.0) fky=dkt*vij
 


	xnuttr=(xnutt(m,ipj)*delx(i)+xnutt(m,ij)*delx(i+1))
     &      /sigmak/(delx(i)+delx(i+1))

	xnuttl=(xnutt(m,imj)*delx(i)+xnutt(m,ij)*delx(i-1))
     &      /sigmak/(delx(i)+delx(i-1))

	viskx=(xnuttr*(ken(m,ipj)-ken(m,ij))/(xi(i+1)-xi(i))
     &       -xnuttl*(ken(m,ij)-ken(m,imj))/(xi(i)-xi(i-1)))/delx(i)

	xnuttt=(xnutt(m,ij)  *dely(j+1)+xnutt(m,ijp)*dely(j))
     &     /sigmak/(dely(j)+dely(j+1))

	xnuttb=(xnutt(m,ij)  *dely(j-1)+xnutt(m,ijm)*dely(j))
     &     /sigmak/(dely(j)+dely(j-1))

	visky=(xnuttt*(ken(m,ijp)-ken(m,ij))/(yj(j+1)-yj(j))
     &       -xnuttr*(ken(m,ij)-ken(m,ijm))/(yj(j)-yj(j-1)))/dely(j)

	ke(m,ij)=ken(m,ij)
	1         +delt*(-fkx-fky+viskx+visky+prod(m,ij)-eps(m,ij))
	
	end if
 10	continue
  

c
c       boundary condition for ke
c
        do 200 i=2,im1
        do 200 j=2,jm1
c
c       free surface no gradient
c
           ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
         
        nff=nf(ij)
        nfr=nf(ipj)
        nfl=nf(imj)
        nft=nf(ijp)
        nfb=nf(ijm)
 
        if(nff.ge.5) then
        if(nft.lt.5) ke(m,ij)=ke(m,ijp)
        if(nfb.lt.5) ke(m,ij)=ke(m,ijm)
        if(nfr.lt.5) ke(m,ij)=ke(m,ipj)
        if(nfl.lt.5) ke(m,ij)=ke(m,imj)
	  end if


	if(f(ij).lt.emf.or.ac(ij).lt.emf) then
	ke(m,ij)=0.d0
	xnutt(m,ij)=0.d0
	end if

c
c       wall boundary
c
 

 200	continue

 999	continue
	return
	end


	subroutine length_scale
c=======================================================================
c   Purpose - define length scale and initiation
c  	 Qun Zhao Nov. 15, 2002 
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
	 
       do i=1,imax
       do j=1,jmax
       ij=(j-1)*imax+i
       do m=1, max 
       length(m,ij)=sqrt(delx(i)*dely(j))/(2**(m-1))       
	ke(m,ij)=em6
	ken(m,ij)=em6
	eps(m,ij)=em6
	prod(m,ij)=em6
	xnutt(m,ij)=em6
	end do
	m=max+1
       length(m,ij)=sqrt(delx(i)*dely(j))/(2**(m-1))    	
       end do
       end do
       
	return
	end			




	subroutine shear

c=======================================================================c
c   Purpose - compute shear strain rate 
c  	 Qun Zhao Nov. 15, 2002 tested O.K.
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
  	do i=2,im1
	do j=2, jm1
	
          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax

	uijp=un(ijp)
	uij=un(ij)
	uijm=un(ijm)

	uimjp=un(imjp)
	uimj=un(imj)
	uimjm=un(imjm)

	vipj=vn(ipj)
	vij=vn(ij)
	vimj=vn(imj)

	vipjm=vn(ipjm)
	vijm=vn(ijm)
	vimjm=vn(imjm)

	if(f(ijp).lt.emf.and.f(ipjp).lt.emf) uijp=0.
	if(f(ij).lt.emf.and.f(ipj).lt.emf)   uij=0.
	if(f(ijm).lt.emf.and.f(ipjm).lt.emf) uijm=0.

	if(f(ijp).lt.emf.and.f(imjp).lt.emf) uimjp=0.
	if(f(ij).lt.emf.and.f(imj).lt.emf) uimj=0.
	if(f(ijm).lt.emf.and.f(imjm).lt.emf) uimjm=0.

	if(f(ipj).lt.emf.and.f(ipjp).lt.emf) vipj=0.
	if(f(ij).lt.emf.and.f(ijp).lt.emf) vij=0.
	if(f(imj).lt.emf.and.f(imjp).lt.emf) vimj=0.

	if(f(ipj).lt.emf.and.f(ipjm).lt.emf) vipjm=0.
	if(f(ij).lt.emf.and.f(ijm).lt.emf) vijm=0.
	if(f(imj).lt.emf.and.f(imjm).lt.emf) vimjm=0.

c ---
	if(ar(ijp).lt.0.or.ar(ipjp).lt.0) uijp=0.
	if(ar(ij).lt.0.or.ar(ipj).lt.0) uij=0.
	if(ar(ijm).lt.0.or.ar(ipjm).lt.0) uijm=0.

	if(ar(ijp).lt.0.or.ar(imjp).lt.0) uimjp=0.
	if(ar(ij).lt.0.or.ar(imj).lt.0) uimj=0.
	if(ar(ijm).lt.0.or.ar(imjm).lt.0) uimjm=0.


	if(at(ipj).lt.0.or.at(ipjp).lt.0) vipj=0.
	if(at(ij).lt.0.or.at(ijp).lt.0) vij=0.
	if(at(imj).lt.0.or.at(imjp).lt.0) vimj=0.


	if(at(ipj).lt.0.or.at(ipjm).lt.0) vipjm=0.
	if(at(ij).lt.0.or.at(ijm).lt.0) vijm=0.
	if(at(imj).lt.0.or.at(imjm).lt.0) vimjm=0.

        urt=(uijp+uij)/2.
        urb=(uij+uijm)/2.
        ult=(uimjp+uimj)/2.
        ulb=(uimj+uimjm)/2.
        ut=(urt+ult)/2.
        ub=(urb+ulb)/2.
	
        ux=(uij-uimj)/delx(i)
	uy=(ut-ub)/dely(j)

        vrt=(vipj+vij)/2.
        vrb=(vipjm+vijm)/2.
        vlt=(vij+vimj)/2.
        vlb=(vijm+vimjm)/2.
	vr=(vrt+vrb)/2.
	vl=(vlt+vlb)/2.

	vx=(vr-vl)/delx(i)
        vy=(vij-vijm)/dely(j)

 
        sij(ij)=(ux*ux+0.5*(uy+vx)*(uy+vx)+vy*vy) 
   
	if(ac(ij).lt.emf.or.f(ij).lt.emf) sij(ij)=0.


	end do
	end do

c	open(39,file='tmpqun.dat')
c	do  j=1,jmax
c	write(39,99) (sij(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
c 99	format(3x,1200(f16.5)) 

	return
	end


	subroutine dkcal(m,ij)
c=======================================================================
c   Purpose - convective term of k -upwind for now
c  	 Qun Zhao Nov. 15, 2002 
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

          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax

        dkl=(ken(m,ij)-ken(m,imj))/(xi(i)-xi(i-1))
	dkr=(ken(m,ipj)-ken(m,ij))/(xi(i+1)-xi(i))
 	dkb=(ken(m,ij)-ken(m,ijm))/(yj(j)-yj(j-1))
 	dkt=(ken(m,ijp)-ken(m,ij))/(yj(j+1)-yj(j))

c	if(ar(ipj).lt.em6) dkr=(-ken(m,imj)-ken(m,ij))*rdx(i+1)
c	if(ar(imj).lt.em6) dkl=(ken(m,ij)+ken(m,ipj))*rdx(i)
c	if(ar(ijp).lt.em6)
c     &    dkt=(-ken(m,ijm)-ken(m,ij))*2.0/(dely(j)+dely(j+1))
c	if(ar(ijm).lt.em6)
c     &    dkb=(ken(m,ij)+ken(m,ijp))*2.0/(dely(j)+dely(j-1))

	if(ar(ipj).lt.em6) dkr=0.
	if(ar(imj).lt.em6) dkl=0.
	if(at(ijp).lt.em6)
     &    dkt=0.
	if(at(ijm).lt.em6)
     &    dkb=0.

	return
	end


	subroutine anothershear	

c=======================================================================c
c   Purpose - compute shear strain rate 
c  	  rate of strain, fyshi 01/29/03
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

       real*8 txx(nxy),txy(nxy),tyy(nxy)
c###########
c
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... get the strain and rotation tensors
c

c     -----------------------------------------------------------------
      call strain(ar,at,ac,x,y,ri,r,un,vn,div,curlu,im1,jm1,nxy,
     &            txx,tyy,txy)
c     -----------------------------------------------------------------
c
c.... bottom and top boundaries
c
      ijb=2
      ijt=jm1*imax+2
      deltyb=dely(2)
      deltyt=dely(jm1)
      do 110 i=2,im1
        sxyb=un(ijb-1+imax)/deltyb
        sxyt=-un(ijt-1-imax)/deltyt
        txx(ijt)=txx(ijt-imax)
        tyy(ijt)=tyy(ijt-imax)
        txy(ijt)=cvmgt(sxyt,zero,kt.eq.2)
        curlu(ijt)=-2.*txy(ijt)
        txx(ijb)=txx(ijb+imax)
        tyy(ijb)=tyy(ijb+imax)
        txy(ijb+imax)=cvmgt(sxyb,zero,kb.eq.2)
        curlu(ijb+imax)=-2.*txy(ijb+imax)
        ijb=ijb+1
        ijt=ijt+1
  110 continue
c
c.... right and left boundaries
c
      ijl=1
      ijr=imax
      deltxl=delx(2)
      deltxr=delx(im1)
      do 120 j=1,jmax
        if (j.eq.1) then
          ijvl=ijl+1
          ijvr=ijr-1
        else
          ijvl=ijl+1-imax
          ijvr=ijr-1-imax
        endif
        sxyl=vn(ijvl)/deltxl
        sxyr=-vn(ijvr)/deltxr
        txx(ijr)=txx(ijr-1)
        tyy(ijr)=tyy(ijr-1)
        txy(ijr)=cvmgt(sxyr,zero,kr.eq.2)
        curlu(ijr)=2.*txy(ijr)
        txx(ijl)=txx(ijl+1)
        tyy(ijl)=tyy(ijl+1)
        txy(ijl+1)=cvmgt(sxyl,zero,kl.eq.2)
        curlu(ijl+1)=2.*txy(ijl+1)
        ijl=ijl+imax
        ijr=ijr+imax
  120 continue
c
      do 100 i=1,imax
        do 100 j=1,jmax
c
          ij=(j-1)*imax+i
	  sij(ij)=txx(ij)*txx(ij)
     &              +tyy(ij)*tyy(ij)+txy(ij)*txy(ij)
c
  100 continue



	open(39,file='tmpsij.dat')
	do  j=1,jmax
	write(39,99) (sij(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
	end do	
	close(39)
 99	format(3x,1200(f16.5)) 



	return
	end