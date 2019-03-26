 
      subroutine vofadv
c
c ======================================================================
c
c   Purpose -
c     compute volume of fluid advection fluxes from the newly
c     determined velocity field, updating the volume of fluid
c     function f(i,j)
c
c   VOFADV is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE   SRFEMIN
c
c
c   VOFADV calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple   SETARRY    ripple     SETNF    ripple
c          VOFCOR    ripple    VOFERR    ripple   VOFPACK    ripple
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
c      include "32bit.h"
c##############################################################
c
c############
      include "comdk1.h" 
      include "iccgdk.h" 
c############
c
c      dimension dij(1),umom(1),vmom(1)
      dimension dij(nxy),umom(nxy),vmom(nxy)
      equivalence (dij(1),cg(1)),(umom(1),sc1(1)),
     &            (vmom(1),sc2(1))
      data tiny /1.0d-25/, zero /0.0d0/, tenth /0.10d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... reset the pertinent arrays
c     ---------------------------------
c      call setarry (umom(1),0.0d0,nxy)
c      call setarry (vmom(1),0.0d0,nxy)
      call setarry (umom,0.0d0,nxy)
      call setarry (vmom,0.0d0,nxy)
c     ---------------------------------
c
c.... Compute cell momenta
c
      do 200 j=2,jm1
        do 200 i=2,im1
          ij=(j-1)*imax+i
          ip=i+1
          jp=j+1
          ipj=ij+1
          ijp=ij+imax
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorf=(delx(ip)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(ip))
          rhotf=(dely(jp)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(jp))
          umom(ij)=rhorf*u(ij)
          vmom(ij)=rhotf*v(ij)
  200 continue
c
      if (ncyc.lt.1) go to 100
c
      flgc=0.0
      do 50 j=1,jm1
        do 50 i=1,im1
          ij=(j-1)*imax+i
          vx=u(ij)*delt
          vy=v(ij)*delt
          abvx=abs(vx)
          abvy=abs(vy)
          if (abvx.le.fcvlim*delx(i).and.abvy.le.fcvlim*dely(j))
     &       go to 5
          if (ac(ij).gt.em6) then
            flgc=1.
            write (iotty,150) ncyc,t,delt,i,j,fcvlim,abvx,delx(i),
     &                        abvy,dely(j)
            write (9,150) ncyc,t,delt,i,j,fcvlim,abvx,delx(i),
     &                        abvy,dely(j)
          endif
c
    5     if (mod(ncyc,2).eq.0) go to 30
   10     continue
          if ((ar(ij).lt.em6).and.(mod(ncyc,2).eq.0)) go to 50
          if (ar(ij).lt.em6) go to 30
          ia=i+1
          id=i
          idm=max0(i-1,1)
          ijdm=(j-1)*imax+idm
          ardm=ar(ijdm)
          rb=ar(ij)*r(i)+tiny
          ra=ac(ij+1)*ri(i+1)+tiny
          rd=ac(ij)*ri(i)+tiny
          incf=1
          incu=0
          if (vx.ge.0.0) go to 20
          ia=i
          id=i+1
          idm=min0(i+2,imax)
          ijdm=(j-1)*imax+idm
          ardm=ar(ijdm-1)
          ra=ac(ij)*ri(i)+tiny
          rd=ac(ij+1)*ri(i+1)+tiny
          incf=-1
          incu=-1
c
   20     continue
          ija=(j-1)*imax+ia
          ijd=(j-1)*imax+id
          ijdu=(j-1)*imax+id+incu
          iad=ia
          if (nf(ijd).eq.3.or.nf(ijd).eq.4) iad=id
          if (fn(ija).lt.emf.or.fn(ijdm).lt.emf) iad=ia
          ijad=(j-1)*imax+iad
          fdm=dmax1(fn(ijdm),fn(ijd),tenth)
          if(ardm.lt.em6) fdm=1.0d0
          fx1=fn(ijad)*abs(vx)+dmax1((fdm-fn(ijad))*abs(vx)-
     &         (fdm-fn(ijd))*delx(id),zero)
          fx=dmin1(fx1,fn(ijd)*delx(id)*rd/rb)
          delfd=fx*rdx(id)*(rb/rd)
          epsd=f(ijd)-delfd
          delfd=cvmgt(delfd+epsd,delfd,epsd.lt.zero)
          delfa=rdx(ia)*rd*delfd/(rdx(id)*ra)
          f(ijd)=f(ijd)-delfd
          f(ija)=f(ija)+delfa
          if (ac(ijd).lt.em6) f(ijd)=fn(ijd)
          if (ac(ija).lt.em6) f(ija)=fn(ija)
          delma=delfa*rhof*cvol(ijd)
          delmd=delfd*rhof*cvol(ijd)
          umom(ija)=umom(ija)+delma*u(ijdu)
          umom(ijd)=umom(ijd)-delmd*u(ijdu)
          if (mod(ncyc,2).eq.0) go to 50
c
   30     if (at(ij).lt.em6) go to 45
          ja=j+1
          jd=j
          jdm=max0(j-1,1)
          ijdm=(jdm-1)*imax+i
          atdm=at(ijdm)
          rb=at(ij)+tiny
          ra=ac(ij+imax)+tiny
          rd=ac(ij)+tiny
          incf=1
          incv=0
          if (vy.ge.0.0) go to 40
          ja=j
          jd=j+1
          jdm=min0(j+2,jmax)
          ijdm=(jdm-1)*imax+i
          atdm=at(ijdm-imax)
          ra=ac(ij)+tiny
          rd=ac(ij+imax)+tiny
          incf=-1
          incv=-1
c
   40     continue
          ija=(ja-1)*imax+i
          ijd=(jd-1)*imax+i
          ijdv=(jd-1+incv)*imax+i
          jad=ja
          if (nf(ijd).eq.1.or.nf(ijd).eq.2) jad=jd
          if (fn(ija).lt.emf.or.fn(ijdm).lt.emf) jad=ja
          ijad=(jad-1)*imax+i
          fdm=dmax1(fn(ijdm),fn(ijd),tenth)
          if(atdm.lt.em6) fdm=1.0d0
          fy1=fn(ijad)*abs(vy)+dmax1((fdm-fn(ijad))*abs(vy)-
     &         (fdm-fn(ijd))*dely(jd),zero)
          fy=dmin1(fy1,fn(ijd)*dely(jd)*rd/rb)
          delfd=fy*rdy(jd)*(rb/rd)
          epsd=f(ijd)-delfd
          delfd=cvmgt(delfd+epsd,delfd,epsd.lt.zero)
          delfa=rdy(ja)*rd*delfd/(rdy(jd)*ra)
          f(ijd)=f(ijd)-delfd
          f(ija)=f(ija)+delfa
          if (ac(ijd).lt.em6) f(ijd)=fn(ijd)
          if (ac(ija).lt.em6) f(ija)=fn(ija)
          delma=delfa*rhof*cvol(ijd)
          delmd=delfd*rhof*cvol(ijd)
          vmom(ija)=vmom(ija)+delma*v(ijdv)
          vmom(ijd)=vmom(ijd)-delmd*v(ijdv)
   45     if (mod(ncyc,2).eq.0) go to 10
   50 continue
   35 continue
c
c.... defoam it by packing vertically
c
      if (npack.ne.0)
c       -------------------------------------------
     &  call vofpack(2,im1,2,jbar,imax,at,nf,f)
c       -------------------------------------------
c
c.... divergence correction
c
      if (idiv.ne.0)
c       -------------------------------------------
     &  call vofcor(2,im1,2,jm1,imax,nf,fn,f,
     &              ac,at,ar,rdx,rdy,r,ri,u,v,
     &              dij,delt)
c       -------------------------------------------
c
  100 continue
c
c.... make sure the new vof function is
c     above the floor and below the ceiling
c
c     --------------------------------------------------
      call voferr(2,im1,2,jm1,imax,t,ncyc,f,delx,
     &            dely,ri,ac,ar,at,vchgt)
c     --------------------------------------------------
c
c.... Special boundary conditions for f
c
c.... Update boundary conditions
c
      ibcflg=1
c     -------------
      call setnf
      call bc
c     -------------
      ibcflg=0
c
 9999 return
c
  150 format (1x,12hVOFADV ERROR/1x,5hncyc=,i7,1x,2ht=,1pe14.6,1x,
     1         5hdelt=,e12.4,1x,4hi,j=,2i4,1x,7hfcvlim=,e11.3/3x,
     2         5habvx=,e12.4,1x,5hdelx=,e12.4,1x,5habvy=,e12.4,1x,
     3         5hdely=,e12.4)
      end

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C		SUBROUTINES
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 
      subroutine vofpack(i1,i2,j1,j2,nqj,at,nf,f)
c
c ======================================================================
c
c   Purpose -
c     pack the vof function vertically downward
c
c   VOFPACK is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VOFADV
c
c
c   VOFPACK calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
      dimension at(1),nf(1),f(1)
      data emf /1.0d-06/, emf1 /0.999999d0/, one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 60 i=i1,i2
        do 60 j=j1,j2
          ij=(j-1)*nqj+i
          ijp=ij+nqj
          if (at(ij).lt.emf) go to 60
          if (nf(ij).ne.0) go to 60
          if (f(ij).ge.emf1) go to 60
          fadd=dmin1(one-f(ij),f(ijp))
          f(ij)=f(ij)+fadd
          f(ijp)=f(ijp)-fadd
   60 continue
c
      return
      end



 
      subroutine vofcor(i1,i2,j1,j2,nqj,nf,fn,f,ac,at,ar,
     &                  rdx,rdy,r,ri,u,v,div,dt)
c
c ======================================================================
c
c   Purpose -
c     correct the vof function for velocity divergence errors
c
c   VOFCOR is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VOFADV
c
c
c   VOFCOR calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
      dimension at(1),nf(1),f(1),fn(1),ac(1),ar(1),
     &          rdx(1),rdy(1),r(1),ri(1),u(1),v(1),
     &          div(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 80 j=j1,j2
        do 80 i=i1,i2
          ij=(j-1)*nqj+i
          imj=ij-1
          ijm=ij-nqj
          if (ac(ij).le.0.) go to 80
          div(ij)=(rdx(i)*(ar(ij)*r(i)*u(ij)-
     1             ar(imj)*r(i-1)*u(imj))+rdy(j)*
     2             (at(ij)*ri(i)*v(ij)-at(ijm)*ri(i)*
     3             v(ijm)))/(ac(ij)*ri(i))
          f(ij)=f(ij)+dt*fn(ij)*div(ij)
   80 continue
c
      return
      end


 
      subroutine voferr(i1,i2,j1,j2,nqj,t,ncyc,f,dx,dy,
     &                  ri,ac,ar,at,vchgt)
c
c ======================================================================
c
c   Purpose -
c     Correct the VOF function if it lies below the floor
c     (emf) or above the ceiling (emf1)
c
c   VOFERR is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VOFADV
c
c
c   VOFERR calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
      dimension f(1),dx(1),dy(1),ri(1),ac(1),
     &          ar(1),at(1)
      data emf /1.0d-06/,emf1 /0.999999d0/,tpi /6.283185307d0/,
     &     em6 /1.0d-06/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 140 j=j1,j2
        do 140 i=i1,i2
          ij=(j-1)*nqj+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+nqj
          ijm=ij-nqj
          if (ac(ij).lt.em6) go to 140
          if ((f(ij).eq.0.0).or.(f(ij).eq.1.0d0)
     &      .or.(f(ij).gt.emf.and.f(ij).lt.emf1)) go to 140
          vchg=0.0
          if (f(ij).ge.emf1) go to 110
          vchg=-f(ij)
          if (abs(vchg).gt.100.*emf)
     &             write (9,160) ncyc,t,i,j,f(ij)
          f(ij)=0.0
          go to 120
c
  110     continue
          vchg=1.0d0-f(ij)
          if (f(ij).gt.1.0d0+100.*emf)
     &             write (9,170) ncyc,t,i,j,-vchg
          f(ij)=1.0d0
c
  120     vchgt=vchgt+vchg*dx(i)*dy(j)*ac(ij)*ri(i)*tpi
  140 continue
c
      return
  160 format (" VOF too small on cycle ",i5,", time ",
     &        1pe13.5," : F(",i2,",",i2,") = ",1pe14.7)
  170 format (" VOF too large on cycle ",i5,", time ",
     &        1pe13.5," : F(",i2,",",i2,") - 1.0 = ",1pe13.7)
      end
