

 
      subroutine aset(ii)
c
c ======================================================================
c
c   Purpose -
c     Initialize any interior obstacles with the function f(x,y);
c     for f(x,y) < 0, the point (x,y) is either open or closed
c     to the flow for ioh=0 or 1, respectively
c
c   ASET is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   ASET calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c             OBS    ripple   SETARRY    ripple       FXY    ripple
c
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
      dimension iflg(5),dis(4),xm(5),ym(5),
     &          grdnwcx(nxy),grdnwcy(nxy),tsin(nxy),sine(nxy)
      equivalence (b(1),grdnwcx(1)),(c(1),grdnwcy(1)),
     &            (d(1),tsin(1)),(e(1),sine(1))
c
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      xmode=2.*pi/(x(im1)-x(1))
      ymode=2.*pi/(y(jm1)-y(1))
      rad2o2=sqrt(2.0d0)/2.0d0
c
c     -------------------------------------
      call setarry (tsin(1),-500.0d0,nxy)
      call setarry (sine(1),0.0d0,nxy)
      call setarry (a(1),0.0d0,nxy)
      call setarry (grdnwcx(1),0.0d0,nxy)
      call setarry (grdnwcy(1),0.0d0,nxy)
      call setarry (gradnwx(1),0.0d0,nxy)
      call setarry (gradnwy(1),0.0d0,nxy)
c     -------------------------------------
c
      if (nobs(ii).le.0) go to 240
c
      do 230 k=1,nobs(ii)
        do 221 j=1,jmax
          do 220 i=1,imax
            ij=(j-1)*imax+i
            rdxdy=1.0/(delx(i)*dely(j))
            do 60 m=1,4
              go to (10,20,30,40), m
   10         x1=x(i)
		    if(j.ne.1)y1=y(j-1)
		    if(j.eq.1)y1=y(j)-dely(j)
c              y1=cvmgt(y(j)-dely(j),y(j-1),j.eq.1)
              dis(1)=dely(j)
              go to 50
   20         y1=y(j)
              x1=x(i)
              dis(2)=delx(i)
              go to 50
   30		    if(i.ne.1)x1=x(i-1)
		    if(i.eq.1)x1=x(i)-delx(i)
c  30         x1=cvmgt(x(i)-delx(i),x(i-1),i.eq.1)
              y1=y(j)
              dis(3)=dely(j)
              go to 50

   40		    if(j.ne.1)y1=y(j-1)
		    if(j.eq.1)y1=y(j)-dely(j)
   		    if(i.ne.1)x1=x(i-1)
		    if(i.eq.1)x1=x(i)-delx(i)
c  40         y1=cvmgt(y(j)-dely(j),y(j-1),j.eq.1)
c              x1=cvmgt(x(i)-delx(i),x(i-1),i.eq.1)
              dis(4)=delx(i)
   50         iflg(m)=0
              fconic=fxy(x1,y1,oa1(ii,k),oa2(ii,k),ob1(ii,k),
     &           ob2(ii,k),oc1(ii,k),
     &           oc2(ii,k),od1(ii,k),nxo(ii,k)*xmode,od2(ii,k),
     &           mxo(ii,k)*xmode,
     &           oe1(ii,k),nyo(ii,k)*ymode,oe2(ii,k),myo(ii,k)*ymode)
              if (abs(fconic).le.1.0d-12) fconic=0.0
              if (fconic.le.0.0) iflg(m)=1
              xm(m)=x1
              ym(m)=y1
   60       continue
c
            iflg(5)=iflg(1)
            xm(5)=xm(1)
            ym(5)=ym(1)
            iflgs=0
c
            do 70 m=1,4
              iflgs=iflgs+iflg(m)
   70       continue
c
            brij=0.0
            btij=0.0
            if (iflgs.eq.0) go to 220
            if (iflgs.lt.4) go to 80
            bij=1.0d0
            brij=1.0d0
            btij=1.0d0
            go to 200
   80       if (iflg(1).eq.1.and.iflg(2).eq.1) brij=1.0d0
            if (iflg(2).eq.1.and.iflg(3).eq.1) btij=1.0d0
c
            do 160 m=1,4
              if (iflg(m).eq.iflg(m+1)) go to 160
              x1=xm(m)
              y1=ym(m)
              x2=xm(m+1)
              y2=ym(m+1)
              if (iflg(m).eq.0) go to 90
              x2=xm(m)
              y2=ym(m)
              x1=xm(m+1)
              y1=ym(m+1)
   90         epsif=0.001d0*(abs(x2-x1)+abs(y2-y1))
              smn=0.0
              fmn=fxy(x2,y2,oa1(ii,k),oa2(ii,k),ob1(ii,k),
     &           ob2(ii,k),oc1(ii,k),
     &           oc2(ii,k),od1(ii,k),nxo(ii,k)*xmode,od2(ii,k),
     &           mxo(ii,k)*xmode,
     &           oe1(ii,k),nyo(ii,k)*ymode,oe2(ii,k),myo(ii,k)*ymode)
              smx=1.0d0
              fmx=fxy(x1,y1,oa1(ii,k),oa2(ii,k),ob1(ii,k),
     &           ob2(ii,k),oc1(ii,k),
     &           oc2(ii,k),od1(ii,k),nxo(ii,k)*xmode,od2(ii,k),
     &           mxo(ii,k)*xmode,
     &           oe1(ii,k),nyo(ii,k)*ymode,oe2(ii,k),myo(ii,k)*ymode)
              s=0.5
  100         xt=s*x1+(1.0-s)*x2
              yt=s*y1+(1.0-s)*y2
              fs=fxy(xt,yt,oa1(ii,k),oa2(ii,k),ob1(ii,k),
     &           ob2(ii,k),oc1(ii,k),
     &           oc2(ii,k),od1(ii,k),nxo(ii,k)*xmode,
     &           od2(ii,k),mxo(ii,k)*xmode,
     &           oe1(ii,k),nyo(ii,k)*ymode,oe2(ii,k),myo(ii,k)*ymode)
              if (abs(fs).lt.epsif) go to 130
              if (fs.ge.0.0) go to 110
              fden=abs(fs-fmn)+1.0d-10
              se=s-fs*(s-smn)/fden
              if (se.gt.smx) se=smx
              fmn=fs
              smn=s
              go to 120
  110         fden=abs(fmx-fs)+1.0d-10
              se=s-fs*(smx-s)/fden
              if (se.lt.smn) se=smn
              fmx=fs
              smx=s
  120         si=s-fs*(smx-smn)/(fmx-fmn)
              s=0.5*(se+si)
              go to 100
  130         dis(m)=sqrt((xt-x2)**2+(yt-y2)**2)
              go to (140,150,160,160), m
  140         brij=dis(1)/dely(j)
              go to 160
  150         btij=dis(2)/delx(i)
  160       continue
c
            m=0
            bij=0.0
  170       continue
            m=m+1
            if (m.eq.5) go to 190
            if (iflg(m).eq.0) go to 170
            mp1=m+1
            if (mp1.eq.5) mp1=1
            mm1=m-1
            if (mm1.eq.0) mm1=4
            bij=bij+dis(m)*dis(mm1)
            if (iflg(mp1).eq.1) go to 180
            dis2=dis(m)
  180       continue
            if (iflg(mm1).eq.1) go to 170
            dis1=dis(mm1)
            go to 170
  190       continue
            if (iflgs.eq.3) bij=bij-dis1*dis2
            bij=0.5*bij*rdxdy
            if (bij.gt.1.0) bij=1.0d0
  200       continue
            if (ioh(ii,k).eq.0) go to 210
            bij=-bij
            brij=-brij
            btij=-btij
  210       ac(ij)=ac(ij)+bij
            if(ac(ij).gt.0.9999d0) ac(ij)=1.0d0
            if(ac(ij).lt.0.0001d0) ac(ij)=0.0
            ar(ij)=ar(ij)+brij
            if(ar(ij).gt.0.9999d0) ar(ij)=1.0d0
            if(ar(ij).lt.0.0001d0) ar(ij)=0.0
            at(ij)=at(ij)+btij
            if(at(ij).gt.0.9999d0) at(ij)=1.0d0
            if(at(ij).lt.0.0001d0) at(ij)=0.0
  220     continue
  221   continue
  230 continue
c
c.... get the physical location of the obstacle boundary
c
c     ---------
      call obs(ii)
c     ---------
c
c.... compute obstacle wall outward normals
c
  240 continue
c
c.... The csin array contains the angle (in degrees)
c     the wall makes with the +y-axis according to
c     the following convention:
 
c        wall on the bottom:  135 < |theta| < 180 (default -90.0)
c        wall on the top:  0 < |theta| < 45 (default 90.0)
c        wall on the left:  45 < |theta| < 135 (default -180.0)
c        wall on the right:  45 < |theta| < 135 (default 0.0)
c
c.... Also compute the cell-centered
c     unit outward normal to the wall
c
      do 400 j=2,jm1
        do 400 i=2,im1
          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          imjp=imj+imax
          ipjm=ipj-imax
          imjm=imj-imax
          if ((i.eq.2).and.(ijobs(ii,ij).eq.0)) then
            grdnwx=-1.0d0
            grdnwy=0.0
            tsin(ij)=0.0
            go to 390
          endif
          if ((i.eq.im1).and.(ijobs(ii,ij).eq.0)) then
            grdnwx=1.0d0
            grdnwy=0.0
            tsin(ij)=0.0
            go to 390
          endif
          if ((j.eq.2).and.(ijobs(ii,ij).eq.0)) then
            grdnwy=-1.0d0
            grdnwx=0.0
            tsin(ij)=-1.0d0
            go to 390
          endif
          if ((j.eq.jm1).and.(ijobs(ii,ij).eq.0)) then
            grdnwy=1.0d0
            grdnwx=0.0
            tsin(ij)=1.0d0
            go to 390
          endif
          if (ijobs(ii,ij).ne.0) then
            deltx=0.5*(delx(i+1)+delx(i-1))+delx(i)
            delty=0.5*(dely(j+1)+dely(j-1))+dely(j)
            gradacxt=(ac(ipjp)-ac(imjp))/deltx
            gradacxm=(ac(ipj)-ac(imj))/deltx
            gradacxb=(ac(ipjm)-ac(imjm))/deltx
            gradacyr=(ac(ipjp)-ac(ipjm))/delty
            gradacym=(ac(ijp)-ac(ijm))/delty
            gradacyl=(ac(imjp)-ac(imjm))/delty
            gradacy=(gradacyl+2.*gradacym+gradacyr)/4.
            gradacx=(gradacxt+2.*gradacxm+gradacxb)/4.
            grdacmg=sqrt(gradacx**2+gradacy**2)+1.0d-10
            x1o=xobs(ii,ijobs(ii,ij))
            x2o=xobs(ii,ijobs(ii,ij)+1)
            y1o=yobs(ii,ijobs(ii,ij))
            y2o=yobs(ii,ijobs(ii,ij)+1)
            xveco=x2o-x1o
            yveco=y2o-y1o
            vecmag=sqrt(yveco**2+xveco**2)
            grdnwx=-yveco/vecmag
            grdnwy=xveco/vecmag
            costhta=(grdnwx*gradacx+grdnwy*gradacy)/grdacmg
            grdnwx=cvmgt(-grdnwx,grdnwx,costhta.gt.0.0)
            grdnwy=cvmgt(-grdnwy,grdnwy,costhta.gt.0.0)
            tsin(ij)=grdnwy
            sine(ij)=grdnwx
          endif
  390     grdnwcx(ij)=grdnwx
          grdnwcy(ij)=grdnwy
  400 continue
c
      do 235 ij=1,nxy
        ijobs(ii,ij)=0
  235 continue
c
c.... Convert the cell-centered unit
c     wall normals to vertex-centered
c
      do 410 j=2,jm1
        do 410 i=2,im1
        ij=(j-1)*imax+i
        ipj=ij+1
        ijp=ij+imax
        ipjp=ij+1+imax
        thta=abs(tsin(ij))
        cthta=tsin(ij)
        if (thta.eq.500.0) go to 410
c....   ``Right'' or ``Left'' wall
        if ((cthta.le.rad2o2).and.(cthta.gt.-rad2o2)) then
          if (grdnwcx(ij).gt.0.0) then
            ij1=ipj
            ij2=ipjp
            ijobs(ii,ij1)=4
            ijobs(ii,ij2)=4
          else
            ij1=ijp
            ij2=ij
            ijobs(ii,ij1)=3
            ijobs(ii,ij2)=3
          endif
          go to 5
        endif
c....   ``Bottom'' wall
        if ((cthta.le.-rad2o2).and.(cthta.ge.-1.0d0)) then
          ij1=ij
          ij2=ipj
          ijobs(ii,ij1)=1
          ijobs(ii,ij2)=1
          go to 5
        endif
c....   ``Top'' wall
        if ((cthta.le.1.0d0).and.(cthta.gt.rad2o2)) then
          ij1=ipjp
          ij2=ijp
          ijobs(ii,ij1)=2
          ijobs(ii,ij2)=2
          go to 5
        endif
    5   gradnwx(ij1)=(a(ij1)*gradnwx(ij1)+grdnwcx(ij))/(a(ij1)+1.)
        gradnwy(ij1)=(a(ij1)*gradnwy(ij1)+grdnwcy(ij))/(a(ij1)+1.)
        a(ij1)=a(ij1)+1.0d0
        gradnwx(ij2)=(a(ij2)*gradnwx(ij2)+grdnwcx(ij))/(a(ij2)+1.)
        gradnwy(ij2)=(a(ij2)*gradnwy(ij2)+grdnwcy(ij))/(a(ij2)+1.)
        a(ij2)=a(ij2)+1.0d0
  410 continue
      do 420 j=2,jmax
        do 420 i=2,imax
          ij=(j-1)*imax+i
          gradmag=sqrt(gradnwx(ij)**2 + gradnwy(ij)**2)+
     &            1.0d-25
          gradnwx(ij)=gradnwx(ij)/gradmag
          gradnwy(ij)=gradnwy(ij)/gradmag
  420 continue
c.... Zero out the corners
      gradnwx(imax+2)=0.0
      gradnwx(2*imax)=0.0
      gradnwx(jm1*imax+2)=0.0
      gradnwx(jmax*imax)=0.0
      gradnwy(imax+2)=0.0
      gradnwy(2*imax)=0.0
      gradnwy(jm1*imax+2)=0.0
      gradnwy(jmax*imax)=0.0
      ijobs(ii,imax+2)=0.0
      ijobs(ii,2*imax)=0.0
      ijobs(ii,jm1*imax+2)=0.0
      ijobs(ii,jmax*imax)=0.0
c
      ijl=1
      ijr=imax
      do 280 j=1,jmax
        ar(ijl)=cvmgt(zero,ar(ijl),kl.le.2)
        at(ijl)=cvmgt(zero,at(ijl),kl.le.2)
        ac(ijl)=cvmgt(em10,ac(ijl),kl.le.2)
        ar(ijr-1)=cvmgt(zero,ar(ijr-1),kr.le.2)
        ar(ijr)=cvmgt(zero,ar(ijr),kr.le.2)
        at(ijr)=cvmgt(zero,at(ijr),kr.le.2)
        ac(ijr)=cvmgt(em10,ac(ijr),kr.le.2)
        ijr=ijr+imax
        ijl=ijl+imax
  280 continue
c
      ijb=1
      ijt=jm1*imax+1
      do 300 i=1,imax
        at(ijb)=cvmgt(zero,at(ijb),kb.le.2)
        ar(ijb)=cvmgt(zero,ar(ijb),kb.le.2)
        ac(ijb)=cvmgt(em10,ac(ijb),kb.le.2)
        at(ijt-1)=cvmgt(zero,at(ijt-1),kt.le.2)
        at(ijt)=cvmgt(zero,at(ijt),kt.le.2)
        ar(ijt)=cvmgt(zero,ar(ijt),kt.le.2)
        ac(ijt)=cvmgt(em10,ac(ijt),kt.le.2)
        ijb=ijb+1
        ijt=ijt+1
  300 continue
c
      do 311 j=2,jm1
        do 310 i=2,im1
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          if (ac(ij).gt.em6) go to 310
          ar(ij)=0.0
          ar(imj)=0.0
          at(ij)=0.0
          at(ijm)=0.0
  310   continue
  311 continue
c
c.... set ar and at values for inflow and outflow boundary segments
c     here as update modification depending on application
c
      return
      end


      subroutine obs(ii)
c
c ======================================================================
c
c   Purpose -
c     find the physical location of any obstacles
c
c   OBS is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET
c
c
c   OBS calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
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
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      npt=1
      do 171 i=2,im1
        atr=1.0-em6
        atl=1.0-em6
        atcc=1.0-em6
        do 170 j=2,jm1
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          if (ac(ij).lt.em6) go to 170
          afr=1.0d0
          aft=1.0d0
          afl=1.0d0
          afb=1.0d0
          if (ar(ij).lt.atr) afr=ar(ij)/atr
          if (at(ij).lt.atcc) aft=at(ij)/atcc
          if (ar(imj).lt.atl) afl=ar(imj)/atl
          if (at(ijm).lt.atcc) afb=at(ijm)/atcc
          go to 7
    7     continue
          if (ac(ij).ge.atcc) go to 120
          if ((aft+afb).lt.em6.or.(afl+afr).lt.em6) go to 170
          m=1
          amn=afb+afr
          if ((afr+aft).gt.amn) go to 10
          m=2
          amn=afr+aft
   10     if ((aft+afl).gt.amn) go to 20
          m=3
          amn=aft+afl
   20     if ((afl+afb).gt.amn) go to 30
          m=4
   30     go to (40,60,80,100), m
   40     x1=x(i-1)+aft*delx(i)
          y1=y(j)
          if (aft.lt.1.0d0) go to 50
          y1=y1-afr*dely(j)
   50     x2=x(i-1)
          y2=y(j)-afl*dely(j)
          if (afl.lt.1.0d0) go to 160
          x2=x2+afb*delx(i)
          go to 160
   60     x1=x(i-1)
          y1=y(j-1)+afl*dely(j)
          if (afl.lt.1.0d0) go to 70
          x1=x1+aft*delx(i)
   70     x2=x(i-1)+afb*delx(i)
          y2=y(j-1)
          if (afb.lt.1.0d0) go to 160
          y2=y2+afr*dely(j)
          go to 160
   80     x1=x(i)-afb*delx(i)
          y1=y(j-1)
          if (afb.lt.1.0d0) go to 90
          y1=y1+afl*dely(j)
   90     x2=x(i)
          y2=y(j-1)+afr*dely(j)
          if (afr.lt.1.0d0) go to 160
          x2=x2-aft*delx(i)
          go to 160
  100     x1=x(i)
          y1=y(j)-afr*dely(j)
          if (afr.lt.1.0d0) go to 110
          x1=x1-afb*delx(i)
  110     x2=x(i)-aft*delx(i)
          y2=y(j)
          if (aft.lt.1.0d0) go to 160
          y2=y2-afl*dely(j)
          go to 160
  120     if (afr.gt.em6) go to 130
          x1=x(i)
          y1=y(j-1)
          x2=x1
          y2=y(j)
          xobs(ii,npt)=x1
          xobs(ii,npt+1)=x2
          yobs(ii,npt)=y1
          yobs(ii,npt+1)=y2
          ijobs(ii,ij)=npt
          npt=npt+2
  130     if (aft.gt.em6) go to 140
          x1=x(i-1)
          y1=y(j)
          x2=x(i)
          y2=y1
          xobs(ii,npt)=x1
          xobs(ii,npt+1)=x2
          yobs(ii,npt)=y1
          yobs(ii,npt+1)=y2
          ijobs(ii,ij)=npt
          npt=npt+2
  140     if (afl.gt.em6) go to 150
          x1=x(i-1)
          y1=y(j)
          x2=x1
          y2=y(j-1)
          xobs(ii,npt)=x1
          xobs(ii,npt+1)=x2
          yobs(ii,npt)=y1
          yobs(ii,npt+1)=y2
          ijobs(ii,ij)=npt
          npt=npt+2
  150     if (afb.gt.em6) go to 170
          x1=x(i-1)
          y1=y(j-1)
          x2=x(i)
          y2=y1
          xobs(ii,npt)=x1
          xobs(ii,npt+1)=x2
          yobs(ii,npt)=y1
          yobs(ii,npt+1)=y2
          ijobs(ii,ij)=npt
          npt=npt+2
  160     xobs(ii,npt)=x1
          xobs(ii,npt+1)=x2
          yobs(ii,npt)=y1
          yobs(ii,npt+1)=y2
          ijobs(ii,ij)=npt
          npt=npt+2
  170   continue
  171 continue
      nobspt=npt-1
c
      return
      end

