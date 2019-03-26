      subroutine tension
c
c ======================================================================
c
c   Purpose -
c     Compute the volume forces due to surface tension
c
c   TENSION is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE   SRFEMIN
c
c
c   TENSION calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c         BDYCELL    ripple  CRVATURE    ripple   SETARRY    ripple
c           ISMAX   cftmath   MOLLIFY    ripple
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
      dimension ro(1),gradftx(1),gradfty(1)
      equivalence (a(1),ro(1)),(af(1),gradftx(1)),(bf(1),gradfty(1))
      data tiny /1.0d-25/, zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... zero-out the relevant arrays
c
c     ----------------------------------
      call setarry (ro(1),0.0d0,nxy)
      call setarry (gradrox(1),0.0d0,nxy)
      call setarry (gradroy(1),0.0d0,nxy)
      call setarry (tensx(1),0.0d0,nxy)
      call setarry (tensy(1),0.0d0,nxy)
      call setarry (fsv(1),0.0d0,nxy)
      call setarry (ftilde(1),0.0d0,nxy)
      call setarry (kappa(1),0.0d0,nxy)
c     ----------------------------------
c
c.... Compute cell densities and copy the
c     raw color function, f(ij), into the
c     smoothed color, ftilde(ij)
c
      do 50 j=2,jm1
        do 50 i=2,im1
          ij=(j-1)*imax+i
          ro(ij)=f(ij)*rhof
          ftilde(ij)=ro(ij)
   50 continue
c
c.... Reflect ro(ij) and ftilde(ij) in ghost cells
c
c     -------------------------------------------------------
      call bdycell(im1,jm1,1,1,1,1,0.,0.,0.,0.,pbc,ro)
      call bdycell(im1,jm1,1,1,1,1,0.,0.,0.,0.,pbc,ftilde)
c     -------------------------------------------------------
c
c.... Optionally smooth the color function
c
c     ----------------------------------------------------------
      if (smooth)
     &   call mollify(2,im1,2,jm1,1,imax,nsmooth,ro,ftilde,af)
c     ----------------------------------------------------------
c
c.... Compute vertex-centered normal vectors
c
      do 35 j=2,jm1+1
        do 35 i=2,im1+1
          ij=(j-1)*imax+i
          imjm=ij-1-imax
          ijm=ij-imax
          imj=ij-1
          im=i-1
          jm=j-1
          deltay=0.5*(dely(j)+dely(jm))
          deltax=0.5*(delx(i)+delx(im))
          rhot=(delx(i)*ro(imj)+delx(im)*ro(ij))/(delx(i)+delx(im))
          trhot=(delx(i)*ftilde(imj)+delx(im)*ftilde(ij))/
     &            (delx(i)+delx(im))
          rhob=(delx(i)*ro(imjm)+delx(im)*ro(ijm))/(delx(i)+delx(im))
          trhob=(delx(i)*ftilde(imjm)+delx(im)*ftilde(ijm))/
     &            (delx(i)+delx(im))
          rhol=(dely(j)*ro(imjm)+dely(jm)*ro(imj))/(dely(j)+dely(jm))
          trhol=(dely(j)*ftilde(imjm)+dely(jm)*ftilde(imj))/
     &            (dely(j)+dely(jm))
          rhor=(dely(j)*ro(ijm)+dely(jm)*ro(ij))/(dely(j)+dely(jm))
          trhor=(dely(j)*ftilde(ijm)+dely(jm)*ftilde(ij))/
     &            (dely(j)+dely(jm))
          gradrox(ij)=(rhor-rhol)/deltax
          gradftx(ij)=(trhor-trhol)/deltax
          gradroy(ij)=(rhot-rhob)/deltay
          gradfty(ij)=(trhot-trhob)/deltay
   35 continue
c
c.... Impose wall adhesion via an equilibrium
c     contact angle in all obstacle cells
c
      do 2350 j=2,jmax
        do 2350 i=2,imax
          ij=(j-1)*imax+i
c
	  do ii=1,nobshp
          if (ijobs(ii,ij).eq.0) go to 2350
          if (ijobs(ii,ij).eq.1) thetaeq=cangleb
          if (ijobs(ii,ij).eq.2) thetaeq=canglet
          if (ijobs(ii,ij).eq.3) thetaeq=canglel
          if (ijobs(ii,ij).eq.4) thetaeq=cangler
	  enddo
c
          grdmgro=sqrt(gradrox(ij)**2 + gradroy(ij)**2) + tiny
          grdmgft=sqrt(gradftx(ij)**2 + gradfty(ij)**2) + tiny
          tx=-gradnwy(ij)
          ty=gradnwx(ij)
          csthtar=(tx*gradrox(ij)+ty*gradroy(ij))/grdmgro
          csthtaf=(tx*gradftx(ij)+ty*gradfty(ij))/grdmgft
          txr=cvmgt(-tx,tx,csthtar.lt.zero)
          tyr=cvmgt(-ty,ty,csthtar.lt.zero)
          txf=cvmgt(-tx,tx,csthtaf.lt.zero)
          tyf=cvmgt(-ty,ty,csthtaf.lt.zero)
c
          grdtwrx=grdmgro*txr
          grdtwry=grdmgro*tyr
          grdnwrx=grdmgro*gradnwx(ij)
          grdnwry=grdmgro*gradnwy(ij)
          grdtwfx=grdmgft*txf
          grdtwfy=grdmgft*tyf
          grdnwfx=grdmgft*gradnwx(ij)
          grdnwfy=grdmgft*gradnwy(ij)
c
          gradrox(ij)=grdnwrx*cos(thetaeq)+grdtwrx*sin(thetaeq)
          gradroy(ij)=grdnwry*cos(thetaeq)+grdtwry*sin(thetaeq)
          gradftx(ij)=grdnwfx*cos(thetaeq)+grdtwfx*sin(thetaeq)
          gradfty(ij)=grdnwfy*cos(thetaeq)+grdtwfy*sin(thetaeq)
c
 2350 continue
c
c.... Compute the mean curvature, kappa, using
c     the smoothed color function, ftilde(ij)
c     ----------------------------------------------------------
      call crvature(2,im1,2,jm1,1,imax,ftilde,gradftx,gradfty,
     &              delx,dely,r,ri,kappa,tiny)
c     ----------------------------------------------------------
c
      do 110 j=2,jm1
        do 100 i=2,im1
          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ipjp=ij+imax+1
          ijp=ij+imax
          ijm=ij-imax
          im=i-1
          ip=i+1
          jm=j-1
          jp=j+1
c
          avnx=0.25*(gradrox(ipj)+gradrox(ipjp)
     &            +gradrox(ijp)+gradrox(ij))
          avny=0.25*(gradroy(ipj)+gradroy(ipjp)
     &            +gradroy(ijp)+gradroy(ij))
c
          tensin=sigma*kappa(ij)/rhof
          tensx(ij)=avnx*tensin
          tensx(ij)=cvmgt(zero,tensx(ij),ac(ij).lt.em6)
          tensy(ij)=avny*tensin
          tensy(ij)=cvmgt(zero,tensy(ij),ac(ij).lt.em6)
          fsv(ij)=sqrt(tensx(ij)**2 + tensy(ij)**2)
c
  100   continue
  110 continue
c
      ijmax=ismax(nxy,fsv,1)
      fmax=fsv(ijmax)
      ffloor=0.001*fmax
      do 200 j=2,jm1
        do 200 i=2,im1
          ij = (j-1)*imax + i
          kappa(ij)=cvmgt(kappa(ij),zero,fsv(ij).gt.ffloor)
  200 continue
c
 9999 return
      end


 
      subroutine srffrce
c
c ======================================================================
c
c   Purpose -
c     compute fluid acceleration due to surface forces
c
c   SRFFRCE is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE   SRFEMIN
c
c
c   SRFFRCE calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple
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
      data tiny /1.0d-25/, zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 20 j=2,jm1
        do 20 i=2,im1
          ij=(j-1)*imax+i
          ipj=ij+1
          ijp=ij+imax
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          rhox=rhorc+tiny
          rhoy=rhotc+tiny
c
          fxvij=tensx(ij)
          fyvij=tensy(ij)
          fxvipj=tensx(ipj)
          fyvijp=tensy(ijp)
c
          fxrc=(delx(i)*fxvipj+delx(i+1)*fxvij)/(delx(i)+delx(i+1))
          if (gfnctn) fxrc=2.*rhox*fxrc/rhof
          fytc=(dely(j)*fyvijp+dely(j+1)*fyvij)/(dely(j)+dely(j+1))
          if (gfnctn) fytc=2.*rhoy*fytc/rhof
 
          if (ar(ij).lt.1.0d-06) go to 10
          u(ij)=u(ij)+delt*fxrc/rhox
   10     u(ij)=cvmgt(zero,u(ij),(rhorc.lt.frctn*rhof).or.
     &                     (ar(ij).lt.1.0d-06))
c
          if (at(ij).lt.1.0d-06) go to 15
          v(ij)=v(ij)+delt*fytc/rhoy
   15     v(ij)=cvmgt(zero,v(ij),(rhotc.lt.frctn*rhof).or.
     &                (at(ij).lt.1.0d-06))
c
   20 continue
c
c     ----------
      call bc
c     ----------
c
      return
      end
