      subroutine turbulence_mut
c
c ======================================================================
c
c   Purpose -
c     compute xmut from the turbulence model
c
c   turbulence_mut is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VTILDE
c
c
c   turbulence_mut calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          STRAIN    ripple
c
c              fyshi 12/12/02
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c###########
       include "comdk1.h" 
       real*8 txx(nxy),txy(nxy),tyy(nxy),tstar
c###########
c
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... get the strain and rotation tensors
c

	if(itur.eq.1)goto 123

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
	  xmusmag(ij)=0.15*0.15*delx(i)*dely(j)
     &              *sqrt(2.*(txx(ij)*txx(ij)
     &              +tyy(ij)*tyy(ij)+txy(ij)*txy(ij)))
c
  100 continue

  123 continue

c --- for turbulence, xmut will be used instead of xmu fyshi 12/12/02

        tstar=1.
	do j=2,jm1
	  do i=2,im1
          ij=(j-1)*imax+i
          if(t.ge.tstar)then
	  xmut(ij)=rhof*(xnu+
     &         turcoef*tanh(3.14/1.*(t-tstar))*(abs(etah(i)-30.)/8.
     &                 +xnutt(1,ij)*0.)
     &                 +xmusmag(ij))

c     &         turcoef*tanh(3.14/1.*(t-tstar))*xnutt(1,ij)*500.
c     &                 +xmusmag(ij))
c          xmut(ij)=xmut(ij)+min(xnutt(1,ij),100)
          else
	  xmut(ij)=rhof*(xnu+xmusmag(ij))
          endif
	  if(i.ge.turl.and.i.le.turr)then
	  xmut(ij)=rhof*(xnu+xmusmag(ij)+turmax)
	  endif
	  enddo
	enddo
	print*,'----',xmut(450),etah(100),xnutt(1,450)

c	open(39,file='tmp.dat')
c	do  j=1,jmax
c	write(39,99) (xmut(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
c 99	format(3x,800(f16.5)) 
c
      return
      end
