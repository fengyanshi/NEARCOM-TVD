      subroutine vtilde
c
c ======================================================================
c
c   Purpose -
c     compute the incremental velocity change due to
c     gravitational acceleration
c
c   VTILDE is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE   SRFEMIN
c
c   VTILDE calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple   CONVECT    ripple  CONVECTC    ripple
c          STRESS    ripple
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
       include "iccgdk.h" 
c############
c
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 1 j=1,jmax
        do 1 i=1,imax
          ij=(j-1)*imax+i
          u(ij)=un(ij)
          v(ij)=vn(ij)
    1 continue
c
c.... get the stress tensor
c
c     --------------
      call turbulence_mut
c     --------------

c     ------------
      call stress
c     ------------
c
c.... get the turbulent stress tensor
c

        do   j=2,jm1
        do   i=2,im1
          ij=(j-1)*imax+i
	    visx_vt(ij)=0.d0
	    visy_vt(ij)=0.d0
	  end do
	  end do

	if (itur.eq.1)then

c     ----------------
      call stress_vt 
c     ----------------

	endif

      if (conserve) then
c       --------------
        call convectc
c       --------------
      else
c       --------------
        call convect
c       --------------
      endif
c


      visx=0.0
      visy=0.0
      tiny=1.0d-25

c ---    add porous - fyshi 10/18/02
      do 20 j=2,jm1
        do 20 i=2,im1
c
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
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          rhox=rhorc+tiny
          rhoy=rhotc+tiny
c
c....     check for surrounding fluid
          rhobar=rhox
          if (ar(ij).lt.em6) go to 10
c
          if (xnu.eq.0.0) go to 5
c....     compute the x-direction viscous stress
c
          tauxyy=(tauxy(ipjp)-tauxy(ipj))/dely(j)
          tauxxx=(ri(i+1)*tauxx(ipj)-ri(i)*tauxx(ij))
     &                /(r(i)*0.5*(delx(i+1)+delx(i)))
          visx=(tauxxx+tauxyy)/rhof+visx_vt(ij)
c
    5     continue
c....     explicitly update the x-velocity with the
c         x-direction convective flux, viscous stress,
c         body force, and pressure force
          rhox=rhobar*(delx(i+1)+delx(i))
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pr(ii,ij).eq.1)then
          u(ij)=u(ij)+delt*visx*1./(1.+ca(ii))
     &       +delt*gx*porosity(ii)/(1.+ca(ii))
     &       -delt*porosity(ii)/(1.+ca(ii))
     &       *(grav*ap(ii,ij)*u(ij)+grav*bp(ii,ij)
     &        *sqrt(u(ij)*u(ij)+v(ij)*v(ij))*u(ij))   
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          u(ij)=u(ij)+delt*visx+delt*gx
          endif
   10     u(ij)=cvmgt(zero,u(ij),(rhorc.lt.frsurf).or.
     &               (ar(ij).lt.em6))
c
c....     reset y-velocity and check for surrounding fluid
          rhobar=rhoy
          if (at(ij).lt.em6) go to 25
c
          if (xnu.eq.0.0) go to 15
c....     compute the y-direction viscous stress
c
          tauxyx=(r(i)*tauxy(ipjp)-r(i-1)*tauxy(ijp))/(ri(i)*delx(i))
          tauyyy=2.*(tauyy(ijp)-tauyy(ij))/(dely(j+1)+dely(j))
          visy=(tauyyy+tauxyx)/rhof+visy_vt(ij)
c
   15     continue
c....     explicitly update the y-velocity with the
c         y-direction convective flux, viscous stress,
c         body force, and pressure force
          rhoy=rhobar*(dely(j+1)+dely(j))
c porous
          nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pt(ii,ij).eq.1)then
          v(ij)=v(ij)+delt*visy*1./(1.+ca(ii))
     &       +delt*gy*porosity(ii)/(1.+ca(ii))
     &       -delt*porosity(ii)/(1.+ca(ii))
     &       *(grav*ap(ii,ij)*v(ij)+grav*bp(ii,ij)
     &        *sqrt(u(ij)*u(ij)+v(ij)*v(ij))*v(ij))   
          nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          v(ij)=v(ij)+delt*visy+delt*gy
	  endif

   25     v(ij)=cvmgt(zero,v(ij),(rhotc.lt.frsurf).or.
     &               (at(ij).lt.em6))
c
   20 continue
c
c     ---------
      call bc
c     ---------
c
      return
      end


 

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C		SUBROUTINES
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



 
      subroutine stress
c
c ======================================================================
c
c   Purpose -
c     compute the viscous stress tensor
c
c   STRESS is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VTILDE
c
c
c   STRESS calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          STRAIN    ripple
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c###########
       include "comdk1.h" 
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
     &            tauxx,tauyy,tauxy)
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
        tauxx(ijt)=tauxx(ijt-imax)
        tauyy(ijt)=tauyy(ijt-imax)
        tauxy(ijt)=cvmgt(sxyt,zero,kt.eq.2)
        curlu(ijt)=-2.*tauxy(ijt)
        tauxx(ijb)=tauxx(ijb+imax)
        tauyy(ijb)=tauyy(ijb+imax)
        tauxy(ijb+imax)=cvmgt(sxyb,zero,kb.eq.2)
        curlu(ijb+imax)=-2.*tauxy(ijb+imax)
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
        tauxx(ijr)=tauxx(ijr-1)
        tauyy(ijr)=tauyy(ijr-1)
        tauxy(ijr)=cvmgt(sxyr,zero,kr.eq.2)
        curlu(ijr)=2.*tauxy(ijr)
        tauxx(ijl)=tauxx(ijl+1)
        tauyy(ijl)=tauyy(ijl+1)
        tauxy(ijl+1)=cvmgt(sxyl,zero,kl.eq.2)
        curlu(ijl+1)=2.*tauxy(ijl+1)
        ijl=ijl+imax
        ijr=ijr+imax
  120 continue
c
      do 100 i=1,imax
        do 100 j=1,jmax
c
          ij=(j-1)*imax+i
          tauxx(ij)=2.*xmut(ij)*tauxx(ij)
          tauyy(ij)=2.*xmut(ij)*tauyy(ij)
          tauxy(ij)=2.*xmut(ij)*tauxy(ij)
c
  100 continue
c
      return
      end


 
      subroutine strain(ar,at,ac,x,y,ri,r,u,v,div,curlu,im1,jm1,
     &                  nxy,sxx,syy,sxy)
c
c ======================================================================
c
c   Purpose -
c     compute the strain rate tensor
c
c   STRAIN is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         IMPLCTP    STRESS
c
c
c   STRAIN calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c         SETARRY    ripple
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c      include "32bit.h"
c##############################################################
c
      dimension ar(1),at(1),ac(1),x(1),y(1),r(1),
     &          ri(1),u(1),v(1),div(1),curlu(1),
     &          sxx(1),sxy(1),syy(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... reset the pertinent arrays
c
c     --------------------------------
      call setarry (div(1),0.0d0,nxy)
      call setarry (curlu(1),0.0d0,nxy)
      call setarry (sxx(1),0.0d0,nxy)
      call setarry (syy(1),0.0d0,nxy)
      call setarry (sxy(1),0.0d0,nxy)
c     --------------------------------
c
      nxp=im1+1
      do 100 i=2,im1
c
        delx=x(i)-x(i-1)
        delxp=x(i+1)-x(i)
c
        do 100 j=2,jm1
c
          ij=(j-1)*nxp+i
          if (ac(ij).lt.1.0d-06) go to 100
          ipj=ij+1
          imj=ij-1
          ijp=ij+nxp
          ijm=ij-nxp
          ipjp=ijp+1
          ipjm=ijm+1
          imjp=ijp-1
          imjm=ijm-1
c
          dely=y(j)-y(j-1)
          delyp=y(j+1)-y(j)
c
          drutdxc=(ar(ij)*r(i)*u(ij)-ar(imj)*r(i-1)*u(imj))/delx
          dudxc=(u(ij)-u(imj))/delx
          dvtdyc=(at(ij)*v(ij)-at(ijm)*v(ijm))/dely
          dvdyc=(v(ij)-v(ijm))/dely
          dudytr=(u(ijp)-u(ij))/(0.5*(dely+delyp))
          dvdxtr=(v(ipj)-v(ij))/(0.5*(delx+delxp))
c
          sxx(ij)=dudxc
          syy(ij)=dvdyc
          sxy(ipjp)=0.5*(dudytr+dvdxtr)
c
          div(ij)=(drutdxc/ri(i)+dvtdyc)/ac(ij)
          curlu(ipjp)=dvdxtr-dudytr
c
  100 continue
c
      return
      end


 
      subroutine convectc
c
c ======================================================================
c
c   Purpose -
c     compute the momentum-conservative transport terms needed for
c     changing the velocities from a Lagrangian to an
c     Eulerian frame
c
c   CONVECTC is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VTILDE
c
c
c   CONVECTC calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple     DVCAL    ripple   VELGRAD    ripple
c
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
       include "iccgdk.h" 
c############
c
      dimension ufc(nxy),vfc(nxy),rhotf(nxy),rhorf(nxy),
     &          uvx(nxy),vvx(nxy),util(nxy),vtil(nxy),
     &          rgux(nxy),rgvy(nxy),rguy(nxy),rgvx(nxy)
      equivalence (ufc(1),a(1)),(vfc(1),b(1)),
     &            (rhotf(1),c(1)),(rhorf(1),d(1)),
     &            (rgux(1),e(1)),(rgvy(1),af(1)),
     &            (uvx(1),ag(1)),(vvx(1),bg(1)),
     &            (util(1),ah(1)),(vtil(1),bh(1)),
     &            (rguy(1),bf(1)),(rgvx(1),cf(1))
c
      data tiny /1.0d-25/, zero /0.0d0/, one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      dtadv=0.5*delt
      nsubcyc=0
c
c.... Get the flux velocities and and face densities
c
      nsurf=0
      do 5 j=1,jmax
        do 5 i=1,imax
          ip=i+1
          if (i.eq.imax) ip=i
          im=i-1
          if (i.eq.1) im=i
          jp=j+1
          if (j.eq.jmax) jp=j
          jm=j-1
          if (j.eq.1) jm=j
          ij=(j-1)*imax+i
          ijfr(ij)=0
          ipj=ij+1
          if (i.eq.imax) ipj=ij
          ijp=ij+imax
          if (j.eq.jmax) ijp=ij
          imj=ij-1
          if (i.eq.1) imj=ij
          ijm=ij-imax
          if (j.eq.1) ijm=ij
          imjm=imj-imax
          if (j.eq.1) imjm=imj
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorf(ij)=(delx(ip)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(ip))
          rhotf(ij)=(dely(jp)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(jp))
          ufc(ij)=0.5*(u(ij)+u(imj))
          vfc(ij)=0.5*(v(ij)+v(ijm))
          vvx(ij)=(delx(i)*v(imjm)+delx(im)*v(ijm))/(delx(i)+delx(im))
          uvx(ij)=(dely(j)*u(imjm)+dely(jm)*u(imj))/(dely(j)+dely(jm))
          util(ij)=u(ij)
          vtil(ij)=v(ij)
          srfij=(1.0-f(ij))*(f(ij)-0.0)
          srfipj=(1.0-f(ipj))*(f(ipj)-0.0)
          srfimj=(1.0-f(imj))*(f(imj)-0.0)
          srfijp=(1.0-f(ijp))*(f(ijp)-0.0)
          srfijm=(1.0-f(ijm))*(f(ijm)-0.0)
          srf = srfij+srfipj+srfimj+srfijp+srfijm
          if (srf.gt.1.0d-12 .and. i.ne.1
     &        .and. i.ne.imax .and. j.ne.1
     &        .and. j.ne.jmax) then
            nsurf=nsurf+1
            ijfr(nsurf)=ij
          endif
    5 continue
c
c     --------------------------------------------------------
      call velgrad(rhorf,rhotf,util,vtil,rgux,rguy,rgvx,rgvy)
c     --------------------------------------------------------
c
c.... X sweep
c
    1 do 20 j=2,jm1
        do 20 i=2,im1
c
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
          volrf=r(i)*dely(j)*0.5*(delx(i)+delx(i+1))
          voltf=ri(i)*delx(i)*0.5*(dely(j)+dely(j+1))
c
          rdelx=1.0/(0.5*(delx(i)+delx(i+1)))
          deltx=1./rdelx
          rdelxl=1.0/(0.5*(delx(i)+delx(i-1)))
          deltxl=1./rdelxl
          rdely=1.0/(0.5*(dely(j)+dely(j+1)))
          if (ar(ij).lt.em6) go to 10
c
c....     X-momentum
c
          sgurr=(1.-sign(one,ufc(ipj)))/2.
          sgurl=(1.+sign(one,ufc(ipj)))/2.
          pmr=sign(one,ufc(ipj))
          flxr=1.-abs(ufc(ipj))*dtadv/delx(i+1)
          sgulr=(1.-sign(one,ufc(ij)))/2.
          sgull=(1.+sign(one,ufc(ij)))/2.
          pml=sign(one,ufc(ij))
          flxl=1.-abs(ufc(ij))*dtadv/delx(i)
          urr=0.5*opalp*util(ipj)+0.5*omalp*util(ij)
     &         +vanleer*pmr*0.5*delx(i+1)*rgux(ipj)*flxr
          url=0.5*opalp*util(ij)+0.5*omalp*util(ipj)
     &          +vanleer*pmr*0.5*delx(i+1)*rgux(ij)*flxr
          donorr=sgurr*urr+sgurl*url
          ulr=0.5*opalp*util(ij)+0.5*omalp*util(imj)
     &          +vanleer*pml*0.5*delx(i)*rgux(ij)*flxl
          ull=0.5*opalp*util(imj)+0.5*omalp*util(ij)
     &          +vanleer*pml*0.5*delx(i)*rgux(imj)*flxl
          donorl=sgulr*ulr+sgull*ull
          fuxr=ac(ipj)*ri(i+1)*ufc(ipj)*donorr
          fuxl=ac(ij)*ri(i)*ufc(ij)*donorl
c
          fux=dely(j)*(fuxl-fuxr)
c
          volrf=ar(ij)*volrf+tiny
          util(ij)=util(ij)+dtadv*fux/volrf
   10     util(ij)=cvmgt(zero,util(ij),(rhorf(ij).lt.frsurf).or.
     &               (ar(ij).lt.em6))
c
          if (at(ij).lt.em6) go to 25
c
c....     Y-momentum
c
          aravg=(dely(j+1)*ar(ij)+dely(j)*ar(ijp))/
     &                (dely(j)+dely(j+1))
          alavg=(dely(j+1)*ar(imj)+dely(j)*ar(imjp))/
     &                (dely(j)+dely(j+1))
          sgvrr=(1.-sign(one,uvx(ipjp)))/2.
          sgvrl=(1.+sign(one,uvx(ipjp)))/2.
          pmr=sign(one,uvx(ipjp))
          flxr=1.-abs(uvx(ipjp))*dtadv*rdelx
          sgvlr=(1.-sign(one,uvx(ijp)))/2.
          sgvll=(1.+sign(one,uvx(ijp)))/2.
          pml=sign(one,uvx(ijp))
          flxl=1.-abs(uvx(ijp))*dtadv*rdelxl
          vrr=0.5*opalp*vtil(ipj)+0.5*omalp*vtil(ij)
     &          +vanleer*pmr*0.5*deltx*rgvx(ipj)*flxr
          vrl=0.5*opalp*vtil(ij)+0.5*omalp*vtil(ipj)
     &          +vanleer*pmr*0.5*deltx*rgvx(ij)*flxr
          donorr=sgvrr*vrr+sgvrl*vrl
          vlr=0.5*opalp*vtil(ij)+0.5*omalp*vtil(imj)
     &          +vanleer*pml*0.5*deltxl*rgvx(ij)*flxl
          vll=0.5*opalp*vtil(imj)+0.5*omalp*vtil(ij)
     &          +vanleer*pml*0.5*deltxl*rgvx(imj)*flxl
          donorl=sgvlr*vlr+sgvll*vll
          fvxr=aravg*uvx(ipjp)*r(i)*donorr
          fvxl=alavg*uvx(ijp)*r(i-1)*donorl
c
          fvx=(fvxl-fvxr)/rdely
c
          voltf=at(ij)*voltf+tiny
          vtil(ij)=vtil(ij)+dtadv*fvx/voltf
   25     vtil(ij)=cvmgt(zero,vtil(ij),(rhotf(ij).lt.frsurf).or.
     &               (at(ij).lt.em6))
c
   20 continue
c
      nsubcyc=nsubcyc+1
      if (nsubcyc.eq.4) go to 3
c     -------------------------------------------------------
      call velgrad(rhorf,rhotf,util,vtil,rgux,rguy,rgvx,rgvy)
c     -------------------------------------------------------
c
c.... Y sweep
c
    2 do 40 j=2,jm1
        do 40 i=2,im1
c
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
          volrf=r(i)*dely(j)*0.5*(delx(i)+delx(i+1))
          voltf=ri(i)*delx(i)*0.5*(dely(j)+dely(j+1))
c
          rdelx=1.0/(0.5*(delx(i)+delx(i+1)))
          rdely=1.0/(0.5*(dely(j)+dely(j+1)))
          delty=1./rdely
          rdelyb=1.0/(0.5*(dely(j)+dely(j-1)))
          deltyb=1./rdely
          if (ar(ij).lt.em6) go to 15
c
c....     X-momentum
c
          atavg=(delx(i)*at(ipj)+delx(i+1)*at(ij))/
     &                   (delx(i)+delx(i+1))
          abavg=(delx(i)*at(ipjm)+delx(i+1)*at(ijm))/
     &                   (delx(i)+delx(i+1))
          sgvbt=(1.-sign(one,vvx(ipj)))/2.
          sgvbb=(1.+sign(one,vvx(ipj)))/2.
          pmb=sign(one,vvx(ipj))
          flxb=1.-abs(vvx(ipj))*dtadv*rdelyb
          sgvtt=(1.-sign(one,vvx(ipjp)))/2.
          sgvtb=(1.+sign(one,vvx(ipjp)))/2.
          pmt=sign(one,vvx(ipjp))
          flxt=1.-abs(vvx(ipjp))*dtadv*rdely
          ubb=0.5*opalp*util(ijm)+0.5*omalp*util(ij)
     &          +vanleer*pmb*0.5*deltyb*rguy(ijm)*flxb
          ubt=0.5*opalp*util(ij)+0.5*omalp*util(ijm)
     &          +vanleer*pmb*0.5*deltyb*rguy(ij)*flxb
          donorb=sgvbb*ubb+sgvbt*ubt
          utb=0.5*opalp*util(ij)+0.5*omalp*util(ijp)
     &          +vanleer*pmt*0.5*delty*rguy(ij)*flxt
          utt=0.5*opalp*util(ijp)+0.5*omalp*util(ij)
     &          +vanleer*pmt*0.5*delty*rguy(ijp)*flxt
          donort=sgvtb*utb+sgvtt*utt
          fuyb=abavg*vvx(ipj)*donorb
          fuyt=atavg*vvx(ipjp)*donort
c
          fuy=(fuyb-fuyt)*r(i)/rdelx
c
          volrf=ar(ij)*volrf+tiny
          util(ij)=util(ij)+dtadv*fuy/volrf
   15     util(ij)=cvmgt(zero,util(ij),(rhorf(ij).lt.frsurf).or.
     &               (ar(ij).lt.em6))
c
          if (at(ij).lt.em6) go to 35
c
c....     Y-momentum
c
          sgvtt=(1.-sign(one,vfc(ijp)))/2.
          sgvtb=(1.+sign(one,vfc(ijp)))/2.
          pmt=sign(one,vfc(ijp))
          flxt=1.-abs(vfc(ijp))*dtadv/dely(j+1)
          sgvbt=(1.-sign(one,vfc(ij)))/2.
          sgvbb=(1.+sign(one,vfc(ij)))/2.
          pmb=sign(one,vfc(ij))
          flxb=1.-abs(vfc(ij))*dtadv/dely(j)
          vtt=0.5*opalp*vtil(ijp)+0.5*omalp*vtil(ij)
     &          +vanleer*pmt*0.5*dely(j+1)*rgvy(ijp)*flxt
          vtb=0.5*opalp*vtil(ij)+0.5*omalp*vtil(ijp)
     &          +vanleer*pmt*0.5*dely(j+1)*rgvy(ij)*flxt
          donort=sgvtt*vtt+sgvtb*vtb
          vbt=0.5*opalp*vtil(ij)+0.5*omalp*vtil(ijm)
     &          +vanleer*pmb*0.5*dely(j)*sgvbt*rgvy(ij)*flxb
          vbb=0.5*opalp*vtil(ijm)+0.5*omalp*vtil(ij)
     &          +vanleer*pmb*0.5*dely(j)*sgvbb*rgvy(ijm)*flxb
          donorb=sgvbt*vbt+sgvbb*vbb
          fvyt=ac(ijp)*vfc(ijp)*donort
          fvyb=ac(ij)*vfc(ij)*donorb
c
          fvy=(fvyb-fvyt)*ri(i)*delx(i)
c
          voltf=at(ij)*voltf+tiny
          vtil(ij)=vtil(ij)+dtadv*fvy/voltf
   35     vtil(ij)=cvmgt(zero,vtil(ij),(rhotf(ij).lt.frsurf).or.
     &               (at(ij).lt.em6))
c
   40 continue
c
c     --------------------------------------------------------
      call velgrad(rhorf,rhotf,util,vtil,rgux,rguy,rgvx,rgvy)
c     --------------------------------------------------------
      nsubcyc=nsubcyc+1
      if (nsubcyc.eq.3) go to 1
      go to 2
c
    3 do 60 j=2,jm1
        do 60 i=2,im1
          ij=(j-1)*imax+i
          u(ij)=util(ij)
          v(ij)=vtil(ij)
   60 continue
      if (nsurf.eq.0) go to 45
c
      alphc=1.0d0

c ---     add porous media -- fyshi 10/24

      do 70 nc=1,nsurf
c
          ij=ijfr(nc)
          j=ij/imax + 1
          i=ij-(j-1)*imax
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax
c
c....     get the velocity fluxes
c         -----------------
          call dvcal(un,vn)
c         -----------------
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          aut=(delx(i+1)*at(ij)+delx(i)*at(ipj))/(delx(i)+delx(i+1))
          aub=(delx(i+1)*at(ijm)+delx(i)*at(ipjm))/(delx(i)+delx(i+1))
          avr=(dely(j+1)*ar(ij)+dely(j)*ar(ijp))/(dely(j)+dely(j+1))
          avl=(dely(j+1)*ar(imj)+dely(j)*ar(imjp))/(dely(j)+dely(j+1))
c
          dudr=ac(ipj)*dudr
          dudl=ac(ij)*dudl
          dudt=aut*dudt
          dudb=aub*dudb
          dvdt=ac(ijp)*dvdt
          dvdb=ac(ij)*dvdb
          dvdr=avr*dvdr
          dvdl=avl*dvdl
c
c....     reset x-velocity and check for surrounding fluid
          rdelx=1.0/(delx(i)+delx(i+1))
          rdely=1.0/(dely(j)+dely(j+1))
          if (ar(ij).lt.em6) go to 50
c
c....     compute the x-direction convective flux
          sgu=sign(one,un(ij))
          rdxa=delx(i)+delx(i+1)+alphc*sgu*(delx(i+1)-delx(i))
          rdxa=1.0/rdxa
          fux=rdxa*un(ij)*(delx(i)*dudr+delx(i+1)*dudl+
     &         alphc*sgu*(delx(i+1)*dudl-delx(i)*dudr))
          vbt=(delx(i)*vn(ipj)+delx(i+1)*vn(ij))*rdelx
          vbb=(delx(i)*vn(ipjm)+delx(i+1)*vn(ijm))*rdelx
          vav=0.5*(vbt+vbb)
          dyt=0.5*(dely(j)+dely(j+1))
          dyb=0.5*(dely(j-1)+dely(j))
          sgv=sign(one,vav)
          dya=dyt+dyb+alphc*sgv*(dyt-dyb)
          fuy=(vav/dya)*(dyb*dudt+dyt*dudb+alphc*sgv*
     &         (dyt*dudb-dyb*dudt))
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype	
	  if(ninter.eq.0.and.pr(ii,ij).eq.1)then
          u(ij)=un(ij)-delt*(fux+fuy)/(ar(ij)+tiny)
     &    *1./(1.+ca(ii))/porosity(ii)	  
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          u(ij)=un(ij)-delt*(fux+fuy)/(ar(ij)+tiny)
	  endif

   50     u(ij)=cvmgt(zero,u(ij),(rhorc.lt.frsurf).or.
     &               (ar(ij).lt.em6))
c
c....     reset y-velocity and check for surrounding fluid
          if (at(ij).lt.em6) go to 55
c
c....     compute the y-direction convective flux
          ubr=(dely(j+1)*un(ij)+dely(j)*un(ijp))*rdely
          ubl=(dely(j+1)*un(imj)+dely(j)*un(imjp))*rdely
          uav=0.5*(ubr+ubl)
          dxr=0.5*(delx(i)+delx(i+1))
          dxl=0.5*(delx(i)+delx(i-1))
          sgu=sign(one,uav)
          dxa=dxr+dxl+alphc*sgu*(dxr-dxl)
          fvx=(uav/dxa)*(dxl*dvdr+dxr*dvdl+alphc*sgu*
     &         (dxr*dvdl-dxl*dvdr))
          sgv=sign(one,vn(ij))
          dya=dely(j+1)+dely(j)+alphc*sgv*(dely(j+1)-dely(j))
          fvy=(vn(ij)/dya)*(dely(j)*dvdt+dely(j+1)*dvdb+
     &          alphc*sgv*(dely(j+1)*dvdb-dely(j)*dvdt))
c porous
          nswi=0
	  do ii=1,nportype
	  if(pt(ii,ij).eq.1.)then
          v(ij)=vn(ij)-delt*(fvx+fvy)/(at(ij)+tiny)
     &         *1./(1.+ca(ii))/porosity(ii)	  
	  nswi=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          v(ij)=vn(ij)-delt*(fvx+fvy)/(at(ij)+tiny)
	  endif
   55     v(ij)=cvmgt(zero,v(ij),(rhotc.lt.frsurf).or.
     &               (at(ij).lt.em6))
c
   70 continue
c
c.... update the boundary conditions
c
   45 ibcflg=1
c     ---------
      call bc
c     ---------
      ibcflg=0
c
      return
      end

      subroutine velgrad(rhorf,rhotf,util,vtil,rgux,rguy,rgvx,rgvy)
c
c ======================================================================
c
c   Purpose -
c     Compute velocity gradients du/dx, du/dy, dv/dx, and dv/dy
c     needed for 2nd-order Van Leer-limited advection.
c     Apply boundary conditions on the intermediate tilde velocities.
c
c   VELGRAD is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c        CONVECTC
c
c
c   VELGRAD calls the following subroutines and functions -
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
c############
       include "comdk1.h" 
       include "iccgdk.h" 
c############
c
      dimension util(1),vtil(1),rgux(1),rguy(1),
     &          rgvx(1),rgvy(1),rhorf(1),rhotf(1)
      data tiny /1.0d-25/, zero /0.0d0/, one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... Apply boundary conditions on the tilde velocities
c
      ijl=1
      ijr=imax
      do 100 j=1,jmax
        go to (10,20,30,40,30,5), kl
    5   util(ijl)=uinf(3)
        vtil(ijl)=vtil(ijl+1)
        go to 50
   10   util(ijl)=0.0
        vtil(ijl)=vtil(ijl+1)
        go to 50
   20   util(ijl)=0.0
        vtil(ijl)=-vtil(ijl+1)*delx(1)/delx(2)
        go to 50
   30   if (iter.gt.0) go to 50
        util(ijl)=util(ijl+1)*(x(2)*rx(1)*cyl+1.0-cyl)
        vtil(ijl)=vtil(ijl+1)
        go to 50
   40   util(ijl)=util(ijr-2)
        vtil(ijl)=vtil(ijr-2)
   50   go to (60,70,80,90,80,55), kr
   55   util(ijr-1)=uinf(4)
        vtil(ijr)=vtil(ijr-1)
        go to 95
   60   util(ijr-1)=0.0
        vtil(ijr)=vtil(ijr-1)
        go to 95
   70   util(ijr-1)=0.0
        vtil(ijr)=-vtil(ijr-1)*delx(imax)/delx(im1)
        go to 95
   80   if (iter.gt.0) go to 95
        util(ijr-1)=util(ijr-2)*(x(ibar)*rx(im1)*cyl+1.0-cyl)
        vtil(ijr)=vtil(ijr-1)
        go to 95
   90   util(ijr-1)=util(ijl+1)
        vtil(ijr-1)=vtil(ijl+1)
        vtil(ijr)=vtil(ijl+2)
   95   ijr=ijr+imax
        ijl=ijl+imax
  100 continue
c
      ijb=1
      ijt=imax*jm1+1
      do 200 i=1,imax
        go to (110,120,130,140,130,105), kt
  105   vtil(ijt-imax)=vinf(2)
        util(ijt)=util(ijt-imax)
        go to 150
  110   vtil(ijt-imax)=0.0
        util(ijt)=util(ijt-imax)
        go to 150
  120   vtil(ijt-imax)=0.0
        util(ijt)=-util(ijt-imax)*dely(jmax)/dely(jm1)
        go to 150
  130   if (iter.gt.0) go to 150
        vtil(ijt-imax)=vtil(ijt-2*imax)
        util(ijt)=util(ijt-imax)
        go to 150
  140   vtil(ijt-imax)=vtil(ijb+imax)
        util(ijt-imax)=util(ijb+imax)
        util(ijt)=util(ijb+2*imax)
  150   go to (160,170,180,190,180,155), kb
  155   vtil(ijb)=vinf(1)
        util(ijb)=util(ijb+imax)
        go to 195
  160   vtil(ijb)=0.0
        util(ijb)=util(ijb+imax)
        go to 195
  170   vtil(ijb)=0.0
        util(ijb)=-util(ijb+imax)*dely(1)/dely(2)
        go to 195
  180   if (iter.gt.0) go to 195
        vtil(ijb)=vtil(ijb+imax)
        util(ijb)=util(ijb+imax)
        go to 195
  190   vtil(ijb)=vtil(ijt-2*imax)
        util(ijb)=util(ijt-2*imax)
  195   ijt=ijt+1
        ijb=ijb+1
  200 continue
c
      if (vanleer.eq.0.0) go to 9999
c
      do 6 j=1,jmax
        do 6 i=1,imax
          ip=i+1
          if (i.eq.imax) ip=i
          im=i-1
          if (i.eq.1) im=i
          jp=j+1
          if (j.eq.jmax) jp=j
          jm=j-1
          if (j.eq.1) jm=j
          ij=(j-1)*imax+i
          ipj=ij+1
          if (i.eq.imax) ipj=ij
          ijp=ij+imax
          if (j.eq.jmax) ijp=ij
          imj=ij-1
          if (i.eq.1) imj=ij
          ijm=ij-imax
          if (j.eq.1) ijm=ij
          rgux(ij)=0.0
          rguy(ij)=0.0
          rgvx(ij)=0.0
          rgvy(ij)=0.0
c
          if (rhorf(ij)/rhof .lt. 0.90) go to 7
c
          rgux(ij)=(util(ipj)-util(imj))/(delx(i)+delx(ip))
          unmin=dmin1(util(imj),util(ipj))
          unmax=dmax1(util(imj),util(ipj))
          ugl=util(ij)-0.5*delx(i)*rgux(ij)
          ugr=util(ij)+0.5*delx(ip)*rgux(ij)
          ugmin=dmin1(ugl,ugr)
          ugmax=dmax1(ugl,ugr)
          denom1=dmax1(tiny,ugmax-util(ij))
          denom2=dmax1(tiny,util(ij)-ugmin)
          alph1=(unmax-dmin1(util(ij),unmax))/denom1
          alph2=(dmax1(util(ij),unmin)-unmin)/denom2
          alphvl=dmax1(zero,dmin1(alph1,alph2,one))
          rgux(ij)=alphvl*rgux(ij)
c
          delty=dely(j)+0.5*(dely(jp)+dely(jm))
          rguy(ij)=(util(ijp)-util(ijm))/delty
          unmin=dmin1(util(ijm),util(ijp))
          unmax=dmax1(util(ijm),util(ijp))
          ugb=util(ij)-0.5*dely(j)*rguy(ij)
          ugt=util(ij)+0.5*dely(j)*rguy(ij)
          ugmin=dmin1(ugb,ugt)
          ugmax=dmax1(ugb,ugt)
          denom1=dmax1(tiny,ugmax-util(ij))
          denom2=dmax1(tiny,util(ij)-ugmin)
          alph1=(unmax-dmin1(util(ij),unmax))/denom1
          alph2=(dmax1(util(ij),unmin)-unmin)/denom2
          alphvl=dmax1(zero,dmin1(alph1,alph2,one))
          rguy(ij)=alphvl*rguy(ij)
c
    7     if (rhotf(ij)/rhof .lt. 0.90) go to 6
c
          rgvy(ij)=(vtil(ijp)-vtil(ijm))/(dely(j)+dely(jp))
          vnmin=dmin1(vtil(ijm),vtil(ijp))
          vnmax=dmax1(vtil(ijm),vtil(ijp))
          vgb=vtil(ij)-0.5*dely(j)*rgvy(ij)
          vgt=vtil(ij)+0.5*dely(jp)*rgvy(ij)
          vgmin=dmin1(vgb,vgt)
          vgmax=dmax1(vgb,vgt)
          denom1=dmax1(tiny,vgmax-vtil(ij))
          denom2=dmax1(tiny,vtil(ij)-vgmin)
          alph1=(vnmax-dmin1(vtil(ij),vnmax))/denom1
          alph2=(dmax1(vtil(ij),vnmin)-vnmin)/denom2
          alphvl=dmax1(zero,dmin1(alph1,alph2,one))
          rgvy(ij)=alphvl*rgvy(ij)
c
          deltx=delx(i)+0.5*(delx(ip)+delx(im))
          rgvx(ij)=(vtil(ipj)-vtil(imj))/deltx
          vnmin=dmin1(vtil(imj),vtil(ipj))
          vnmax=dmax1(vtil(imj),vtil(ipj))
          vgl=vtil(ij)-0.5*delx(i)*rgvx(ij)
          vgr=vtil(ij)+0.5*delx(i)*rgvx(ij)
          vgmin=dmin1(vgl,vgr)
          vgmax=dmax1(vgl,vgr)
          denom1=dmax1(tiny,vgmax-vtil(ij))
          denom2=dmax1(tiny,vtil(ij)-vgmin)
          alph1=(vnmax-dmin1(vtil(ij),vnmax))/denom1
          alph2=(dmax1(vtil(ij),vnmin)-vnmin)/denom2
          alphvl=dmax1(zero,dmin1(alph1,alph2,one))
          rgvx(ij)=alphvl*rgvx(ij)
c
    6   continue
c
 9999 return
      end


 
      subroutine dvcal (uc,vc)
c
c ======================================================================
c
c   Purpose -
c     compute the x- and y-velocity fluxes for cell ij
c
c   DVCAL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         CONVECT  CONVECTC
c
c
c   DVCAL calls the following subroutines and functions -
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
c############
       include "comdk1.h" 
c############
c
c.... provide free-slip-like (islip=0) or standard free-slip (islip=1)
c     boundary conditions for all obstacle surface flow cells
c
      data islip / 0 /
      dimension uc(1),vc(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c                                dvdt
c                                 x
c
c
c
c                               vc(ij)
c                   dvdl x--------*--------x dudt,dvdr
c                        |                 |
c                        |                 |
c                        |                 |
c                        |    dudl,dvdb    |      dudr
c                 uc(imj)*        x        *uc(ij)  x
c                        |                 |
c                        |                 |
c                        |                 |
c                        |                 |
c                        x--------*--------x dudb
c                               vc(ijm)
c
      ij=(j-1)*imax+i
      ipj=ij+1
      ijp=ij+imax
      ipjp=ipj+imax
      ipjm=ipj-imax
      ijm=ij-imax
      imj=ij-1
      imjp=imj+imax
      imjm=imj-imax
c
      dudr=(uc(ipj)-uc(ij))*rdx(i+1)
      dudl=(uc(ij)-uc(imj))*rdx(i)
      dudt=(uc(ijp)-uc(ij))*2.0/(dely(j)+dely(j+1))
      dudb=(uc(ij)-uc(ijm))*2.0/(dely(j)+dely(j-1))
      dvdr=(vc(ipj)-vc(ij))*2.0/(delx(i)+delx(i+1))
      dvdl=(vc(ij)-vc(imj))*2.0/(delx(i)+delx(i-1))
      dvdt=(vc(ijp)-vc(ij))*rdy(j+1)
      dvdb=(vc(ij)-vc(ijm))*rdy(j)
c
      if (islip.ne.0) go to 20
      if (ar(ipj).lt.em6) dudr=0.0
      if (ar(imj).lt.em6) dudl=0.0
      if (ar(ijp).lt.em6) dudt=0.0
      if (ar(ijm).lt.em6) dudb=0.0
      if (at(ipj).lt.em6) dvdr=0.0
      if (at(imj).lt.em6) dvdl=0.0
      if (at(ijp).lt.em6) dvdt=0.0
      if (at(ijm).lt.em6) dvdb=0.0
      go to 9999
c
   20 if (islip.ne.1) go to 9999
      if (ar(ijp).lt.em6) dudt=0.0
      if (ar(ijm).lt.em6) dudb=0.0
      if (at(ipj).lt.em6) dvdr=0.0
      if (at(imj).lt.em6) dvdl=0.0
c
 9999 return
      end


 
      subroutine convect
c
c ======================================================================
c
c   Purpose -
c     compute the convective transport terms needed for
c     changing the velocities from a Lagrangian to an
c     Eulerian frame
c
c   CONVECT is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VTILDE
c
c
c   CONVECT calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple     DVCAL    ripple
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
        include "comdk1.h" 
        include "iccgdk.h" 
c############
c
      dimension uc(1),vc(1)
      equivalence (uc(1),a(1)),(vc(1),b(1))
c
      data tiny /1.0d-25/, zero /0.0d0/, one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... save the Lagrangian velocities
c
      do 5 ij=1,nxy
        uc(ij)=un(ij)
        vc(ij)=vn(ij)
    5 continue
c
c.... compute the convective transport terms and
c     transform the Lagrangian velocities to Eulerian
c
c ---  add porous media fyshi 10/24/02
      do 20 j=2,jm1
        do 20 i=2,im1
c
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
c....     get the velocity fluxes
c         -----------------
          call dvcal(uc,vc)
c         -----------------
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          aut=(delx(i+1)*at(ij)+delx(i)*at(ipj))/(delx(i)+delx(i+1))
          aub=(delx(i+1)*at(ijm)+delx(i)*at(ipjm))/(delx(i)+delx(i+1))
          avr=(dely(j+1)*ar(ij)+dely(j)*ar(ijp))/(dely(j)+dely(j+1))
          avl=(dely(j+1)*ar(imj)+dely(j)*ar(imjp))/(dely(j)+dely(j+1))
c
          dudr=ac(ipj)*dudr
          dudl=ac(ij)*dudl
          dudt=aut*dudt
          dudb=aub*dudb
          dvdt=ac(ijp)*dvdt
          dvdb=ac(ij)*dvdb
          dvdr=avr*dvdr
          dvdl=avl*dvdl
c
c....     reset x-velocity and check for surrounding fluid
          rdelx=1.0/(delx(i)+delx(i+1))
          rdely=1.0/(dely(j)+dely(j+1))
          if (ar(ij).lt.em6) go to 10
c
c....     compute the x-direction convective flux
          sgu=sign(one,uc(ij))
          rdxa=delx(i)+delx(i+1)+alpha*sgu*(delx(i+1)-delx(i))
          rdxa=1.0/rdxa
          fux=rdxa*uc(ij)*(delx(i)*dudr+delx(i+1)*dudl+
     &         alpha*sgu*(delx(i+1)*dudl-delx(i)*dudr))
          vbt=(delx(i)*vc(ipj)+delx(i+1)*vc(ij))*rdelx
          vbb=(delx(i)*vc(ipjm)+delx(i+1)*vc(ijm))*rdelx
          vav=0.5*(vbt+vbb)
          dyt=0.5*(dely(j)+dely(j+1))
          dyb=0.5*(dely(j-1)+dely(j))
          sgv=sign(one,vav)
          dya=dyt+dyb+alpha*sgv*(dyt-dyb)
          fuy=(vav/dya)*(dyb*dudt+dyt*dudb+alpha*sgv*
     &         (dyt*dudb-dyb*dudt))
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pr(ii,ij).eq.1.)then
          u(ij)=u(ij)-delt*(fux+fuy)/(ar(ij)+tiny)
     &         *1./(1.+ca(ii))/porosity(ii)
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          u(ij)=u(ij)-delt*(fux+fuy)/(ar(ij)+tiny)
	  endif
   10     u(ij)=cvmgt(zero,u(ij),(rhorc.lt.frsurf).or.
     &               (ar(ij).lt.em6))
c
c....     reset y-velocity and check for surrounding fluid
          if (at(ij).lt.em6) go to 25
c
c....     compute the y-direction convective flux
          ubr=(dely(j+1)*uc(ij)+dely(j)*uc(ijp))*rdely
          ubl=(dely(j+1)*uc(imj)+dely(j)*uc(imjp))*rdely
          uav=0.5*(ubr+ubl)
          dxr=0.5*(delx(i)+delx(i+1))
          dxl=0.5*(delx(i)+delx(i-1))
          sgu=sign(one,uav)
          dxa=dxr+dxl+alpha*sgu*(dxr-dxl)
          fvx=(uav/dxa)*(dxl*dvdr+dxr*dvdl+alpha*sgu*
     &         (dxr*dvdl-dxl*dvdr))
          sgv=sign(one,vc(ij))
          dya=dely(j+1)+dely(j)+alpha*sgv*(dely(j+1)-dely(j))
          fvy=(vc(ij)/dya)*(dely(j)*dvdt+dely(j+1)*dvdb+
     &          alpha*sgv*(dely(j+1)*dvdb-dely(j)*dvdt))
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pt(ii,ij).eq.1.)then
          v(ij)=v(ij)-delt*(fvx+fvy)/(at(ij)+tiny)
     &         *1./(1.+ca(ii))/porosity(ii)
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          v(ij)=v(ij)-delt*(fvx+fvy)/(at(ij)+tiny)
	  endif

   25     v(ij)=cvmgt(zero,v(ij),(rhotc.lt.frsurf).or.
     &               (at(ij).lt.em6))
c
   20 continue
c
c.... update the boundary conditions
c
      ibcflg=1
c     ---------
   45 call bc
c     ---------
      ibcflg=0
c
      return
      end


