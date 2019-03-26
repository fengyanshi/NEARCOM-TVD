
 
      subroutine initvoff
c
c ======================================================================
c
c   Purpose -
c     initialize any free surfaces (via the volume of fluid function)
c     with the function f(x,y); for f(x,y) < 0, the point (x,y) is
c     either within the fluid or beyond the free surface for ifh=0
c     or 1, respectively
c
c   INITVOFF is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   INITVOFF calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c             FXY    ripple
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
       dimension iflg(5), dis(4), xm(5), ym(5)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      xmode=2.*pi/(x(im1)-x(1))
      ymode=2.*pi/(y(jm1)-y(1))
c
      do 230 k=1,nfrsrf
        do 221 j=2,jm1
          do 220 i=2,im1
            ij=(j-1)*imax+i
            rdxdy=1.0/(delx(i)*dely(j))
            do 60 m=1,4
              go to (10,20,30,40), m
   10         x1=x(i)
              y1=y(j-1)
              dis(1)=dely(j)
              go to 50
   20         y1=y(j)
              x1=x(i)
              dis(2)=delx(i)
              go to 50
   30         x1=x(i-1)
              y1=y(j)
              dis(3)=dely(j)
              go to 50
   40         y1=y(j-1)
              x1=x(i-1)
              dis(4)=delx(i)
   50         iflg(m)=0
c             ---------------------------------------------------------
              fconic=fxy(x1,y1,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                   fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                   fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             ---------------------------------------------------------
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
   90         epsif=0.001*(abs(x2-x1)+abs(y2-y1))
              smn=0.0
c             -------------------------------------------------------
              fmn=fxy(x2,y2,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             -------------------------------------------------------
              smx=1.0d0
c             -------------------------------------------------------
              fmx=fxy(x1,y1,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             -------------------------------------------------------
              s=0.5
  100         xt=s*x1+(1.0-s)*x2
              yt=s*y1+(1.0-s)*y2
c             -------------------------------------------------------
              fs=fxy(xt,yt,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             -------------------------------------------------------
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
            if (bij.gt.1.0d0) bij=1.0d0
  200       continue
            if (ifh(k).eq.0) go to 210
            bij=-bij
            brij=-brij
            btij=-btij
  210       f(ij)=f(ij)+bij
            if(f(ij).gt.0.9999) f(ij)=1.0d0
            if(f(ij).lt.0.0001) f(ij)=0.0
            fn(ij)=f(ij)
  220     continue
  221   continue
  230 continue
c
      return
      end

 
      subroutine initreg
c
c ======================================================================
c
c   Purpose -
c     initialize region quantities
c
c   INITREG is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   INITREG calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c         SETARRY    ripple
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
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... default initial velocities are zero
c
c     -------------------------------
      call setarry (u(1),0.0d0,nxy)
      call setarry (v(1),0.0d0,nxy)
c     -------------------------------
c
c.... set initial velocity field into u and v arrays
c
      do 100 i=1,imax
        do 100 j=1,jmax
          ij=(j-1)*imax+i
          ijp=ij+imax
          ipj=ij+1
          if ((at(ij).gt.em6).and.((f(ij).gt.emf).or.
     &          (f(ijp).gt.emf))) v(ij)=vi
          if ((ar(ij).gt.em6).and.((f(ij).gt.emf).or.
     &          (f(ipj).gt.emf))) u(ij)=ui
          un(ij)=u(ij)
          vn(ij)=v(ij)
  100 continue
c
      return
      end



 

