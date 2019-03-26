      subroutine setup
c
c ======================================================================
c
c   Purpose -
c     do the problem setup
c
c   SETUP is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   SETUP calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c         INITREG    ripple  INITVOFF    ripple      ASET    ripple
c              BC    ripple    RINPUT    ripple   SETARRY    ripple
c           SETNF    ripple   SRFEMIN    ripple     TAPIN    ripple
c           EQUIB    ripple   MESHSET    ripple
c
c ======================================================================
c
c##############################################################
        implicit real*8 (a-h,o-z)
c	include "32bit.h"
c##############################################################
c
c############
	include "comdk1.h"    
	include "iccgdk.h" 
c############
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

c
c.... read the input data

c
c     ------------------
      call rinput
c     ------------------

c
c.... wave conditions
c

c     ------------------
      call wave_condition
c     ------------------


c.... If the restart flag is greater than zero,
c     read a restart dump and skip the rest of the setup routine.
c
      if (nrestart.le.0) go to 40
c     -----------
c      call tapin		 
c     -----------
      go to 9998
c
c.... generate the computational mesh
c
  40  continue
c     -------------
c   40 call meshset  ! move to rinput
c     -------------

c.... set default cell array values
c     ----------------------------------------------
      call setarry (ac(1),1.0d0,imax*jmax)
      call setarry (ar(1),1.0d0,imax*jmax)
      call setarry (at(1),1.0d0,imax*jmax)
      call setarry (ac0(1),1.0d0,imax*jmax)
      call setarry (ar0(1),1.0d0,imax*jmax)
      call setarry (at0(1),1.0d0,imax*jmax)

      call setarry (f(imax+1),1.0d0,imax*(jmax-1))
c     ----------------------------------------------
c
c.... generate any interior obstacles or special boundary conditions
c.... set the corresponding values of the flow arrays ar,ac,at
c
c     ----------

      do ii=1,nobshp
      call aset(ii)
      do 20 j=2,jm1
        do 20 i=2,im1
          ij=(j-1)*imax+i	
          ar(ij)=cvmgt(zero,ar(ij),(ar(ij).le.em6).or.
     &               (ar0(ij).le.em6))
	  ar0(ij)=ar(ij)
          at(ij)=cvmgt(zero,at(ij),(at(ij).le.em6).or.
     &               (at0(ij).le.em6))
	  at0(ij)=at(ij)
          ac(ij)=cvmgt(zero,ac(ij),(ac(ij).le.em6).or.
     &               (ac0(ij).le.em6))
	  ac0(ij)=ac(ij)
   20  continue
      enddo



c	open(39,file='ar.dat')
c	do  j=1,jmax
c	write(39,99) (ar(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
c	open(39,file='at.dat')
c	do  j=1,jmax
c	write(39,99) (at(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
c	open(39,file='ac.dat')
c	do  j=1,jmax
c	write(39,99) (ac(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
 99	format(300(f3.0)) 
c	stop

c     ----------

c --- generate porous medium -- fyshi 10/18/02

      do 2343 itype=1,nportype
	do ii=1,imax*jmax
          pc(itype,ii)=1.0d0
          pr(itype,ii)=1.0d0
          pt(itype,ii)=1.0d0
          pbdy(itype,ii)=0.0d0
	enddo
      call porousset(itype)
2343  continue

c
c.... bond number
c
      bond=1.d+100
      if (sigma.gt.0.) bond=rhof*gy*x(im1)**2/sigma
c
c.... determine whether the free-surface is
c     initialized with free-surface functions
c     or a constant fluid height (flht) above
c     the x-axis
c
      if (nfrsrf.gt.0) go to 90
      nzone=nequib
c
c.... determine the equilibrium free-surface relative
c     to a free-surface at a constant height y=flht
c
c     --------------------------------------------------------
c      if (iequib.gt.0)
c     &          call equib (a,af,nzone,bond,dangle,cyl,iotty)  ! qun
c     --------------------------------------------------------

c
c.... Compute the vof functions associated with
c     the equilibrium free-suface configuration
c
      sflht=flht
      do 120 i=1,imax
        do 110 j=2,jmax
          jj=j
          if (.not.(upright)) jj=jmax+2-j
          ij=(jj-1)*imax+i
          f(ij)=1.0d0
          if (iequib.le.0) go to 100
          ldck=float(nequib-1)*xi(i)/x(im1)+1.000001d0
          ldck=min0(nequib,ldck)
          ldck=max0(1,ldck)
          flht=sflht+a(ldck)*x(im1)
          if (.not.(upright)) flht=y(jm1)-flht
  100     continue
          if (flht.gt.y(jj-1).and.flht.lt.y(jj)) then
            if (upright) then
              f(ij)=rdy(jj)*(flht-y(jj-1))
            else
              f(ij)=rdy(jj)*(y(jj)-flht)
            endif
          endif
          if (upright) then
            if (y(jj-1).ge.flht) f(ij)=0.0
          else
            if (y(jj).le.flht) f(ij)=0.0
          endif
  110   continue
        ijb=imax+i
        f(ijb-imax)=f(ijb)
  120 continue
      flht=sflht
      go to 95
c
c.... compute initial vof functions characterized
c     by the user-input free-surface functions
c
c     --------------
   90 call initvoff
c     --------------

c
c.... calculate dtvis and dtsft
c
   95 dlt=1.0d+10
      dtvis=1.0d+10
      dtsft=1.0d+10
      ijsft=0
      ijvis=0
      if ((xnu.eq.0.0).and.(isurf10.eq.0)) go to 140
      do 130 i=2,im1
        do 130 j=2,jm1
          ij=(j-1)*imax+i
          sqj=sqrt(jcb(ij))
          dxsq=delx(i)**2
          dysq=dely(j)**2
          rdsq=dxsq*dysq/(dxsq+dysq)
          rdsq=rdsq/(3.0*xnu+1.0d-60)
          dtvis=dmin1(dtvis,rdsq)
          dlt=dmin1(sqj,dlt)
          if (dlt.eq.sqj) ijsft=ij
          if (dtvis.eq.rdsq) ijvis=ij
  130 continue
      jvis=ijvis/imax + 1
      ivis=ijvis-(jvis-1)*imax
      jsft=ijsft/imax + 1
      isft=ijsft-(jsft-1)*imax
      sigx=sigma
      if (sigx.eq.0.0) sigx=em10
      dtsft=sqrt(rhof*dlt**3/(4.0*pi*sigx))
      if (isurf10.eq.0) dtsft=1.0d+10
      if (xnu.eq.0.0) dtvis=1.0d+10

 140	continue
c
c     ------------------------
c  140 if (emin) call srfemin    ! qun
c     ------------------------

c
c.... initialize region quantities
c
c     -------------
      call initreg
c     -------------
c
c.... update boundary conditions
c
c     ----------------
      call bc
c     ----------------
c
c.... Set the nf flag
c
c     ----------------
      call setnf
c     ----------------
c
c.... optional graphics dump
c
 9998   continue
c	if (dump) then
c        write (8) prbname
c        write (8) jnm,dat,tim,ochn
c        write (8) nxy,nvar,ibar,imax,jbar,jmax,im1,jm1
c        write (8) ixsymplt,iysymplt,dxmin,dymin,scale,vmxfrctn
c        write (8) (x(i),i=1,imax),(y(i),i=1,jmax),
c     &            (xi(i),i=1,imax),(yj(i),i=1,jmax),
c     &            (xf(i),i=1,5),(yf(i),i=1,5)
c      endif


c
 9999 return
  300 format (1x,2hk=,1pe12.4,2x,3hxi=,e12.4,2x,4hper=,e12.4)
  330 format (a80)
  340 format (10x,i10)
  350 format (10x,e20.6)
  360 format (3x,7hndump =,i10/4x,6hqvol =,1pe20.6)
  370 format (10x,7i10)
  380 format (1x,4hcon=,1pe10.3,1x,11hand fcvlim=,e10.3,1x,
     &         41hare incompatible. setting fcvlim=1.3*con.)
  390 format (1x,13hbond number =,1pe12.4)
  400 format (/1x,8hnequib =,i5,1x,29his too big for the dimensions,2i5)
  410 format (1x,17hcutting nequib to,i5,15h and continuing/)
  420 format (1h1)
  440 format (2x,i3,3x,i3,6(3x,1pe12.5))
      end




C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C			SUBROUTINES
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


 
      subroutine rinput
c
c ======================================================================
c
c   Purpose -
c     read the input data; initialize selected and default variables
c
c   RINPUT is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP     TAPIN
c
c
c   RINPUT calls the following subroutines and functions -
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
c############
      include "comdk1.h" 
c############
c
c.... Namelist input
c
      namelist /numparam/ delt,twfin,prtdt,pltdt,
     1                    alpha,kl,kr,kt,kb,autot,
     2                    npack,con,dmpdt,
     3                    dtmax,idiv,
     4                    erriccg,fcvlim,nrestart,
     5                    gfnctn,frctn,conserve,
     6                    cray,itmxiccg,sym
      namelist /fldparam/ xnu,icyl,gx,gy,ui,vi,psat,sigma,
     1                    isurf10,rhof,uinf,cangler,
     2                    canglel,canglet,cangleb,vinf,pbc
      namelist /mesh/ nkx,xl,xc,nxl,nxr,dxmn,nky,yl,yc,nyl,nyr,dymn
      namelist /obstcl/ nobs,oa2,oa1,ob2,ob1,oc2,oc1,ioh,od1,od2,
     &                  oe1,oe2,nxo,mxo,nyo,myo,nobshp
      namelist /freesurf/ nfrsrf,fa2,fa1,fb2,fb1,fc2,fc1,ifh,fd1,fd2,
     &                    fe1,fe2,nxf,mxf,nyf,myf,
     &                    nsmooth,iequib,flht,upright,
     &                    smooth,emin,deltemin,itmxemin,itskpemn
      namelist /graphics/ plots,dump,iout,scale,vmxfrctn,unfrmmsh,
     &                    ixsymplt,iysymplt,ixskip,iyskip,icolor
      namelist /wavecase/ iwave,waveh,wavet,rc,xle

c --- the porous parameters and location
      namelist /porousparam/ porosity,gammaa,grav,alphap,betap,
     &         Ttypical,Dp50,xnumo
      namelist /porous/ nportype,
     &                  nporous,pa2,pa1,pb2,pb1,pc2,pc1,iph,pd1,pd2,
     &                  pe1,pe2,nxp,mxp,nyp,myp
      namelist /smag/ turl, turr, turcoef, turmax
	namelist /turb/ itur
c
      data rsflg /0.0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... Branch around defaults if restarting
c
      if (rsflg.ne.0.0) go to 60
      rsflg = 1.0d0
c
c.... Set defaults for input variables
c
c.... Namelist /numparam/
c
      con = 0.30d0
      dtmax = 0.001d0
      alpha = 1.0d0
      autot = 1.0d0
      erriccg = 1.0d-08
      fcvlim = 0.39d0
      frctn = 5.0d-07
      npack = 0
      nrestart = -1
      idiv = 1
      kl = 1
      kr = 1
      kt = 1
      kb = 1
      itmxiccg = 1000
      gfnctn = .true.
      sym = .true.
      conserve = .true.
      cray = .false.
c
c.... Namelist /fldparam/
c
       sigma = 0.0d0
       rhof = 1.0d0
       psat = 0.0d0
       xnu = 0.0d0
       gx = 0.0d0
       gy = 0.0d0
       ui = 0.0d0
       vi = 0.0d0
       cangler = 90.0d0
       canglel = 90.0d0
       canglet = 90.0d0
       cangleb = 90.0d0
       isurf10 = 0
       icyl = 0
c
c.... Namelist /obstcl/
c
      nobshp = 0
c
c.... Namelist /freesurf/
c
      flht = 0.0d0
      nfrsrf = 0
      nsmooth = 1
      upright = .true.
      smooth = .false.
      emin = .false.
      itmxemin = 500
      itskpemn = 50
      iequib = 0
c
c.... Namelist /graphics/
c
      scale = 2.0d0
      vmxfrctn = 0.001d0
      plots = .false.
      dump = .false.
      unfrmmsh = .false.
      ixsymplt = 0
      iysymplt = 0
      ixskip = 1
      iyskip = 1
      icolor = 0
c ---                initial porous medium - fyshi 10/18/02
c.... Namelist /porous/
c
      do ii=1,nportype
      nporous(ii) = 0
      enddo
c
c.... Initialize constants
c
      emf = 1.0d-06
      em6 = 1.0d-06
      em10 = 1.0d-10
      ep10 = 1.0d+10
      tquit = 60.d0
      em6p1 = 1.000001d0
      pi = 3.14159265359d0
      rpd = 0.0174532925d0
      tbeg = 0.d0
      tpi = 6.283185307d0
      em61 = 0.999999d0
c
c.... Initialize selected variables
c
      ndump = 0
      t = 0.d0
      flgc = 0.d0
      twdmp = 0.d0
      twplt = 0.d0
      twprt = 0.d0
      vchgt = 0.d0
      pbcl = 0.d0
      pbcr = 0.d0
      pbct = 0.d0
      pbcb = 0.d0
      iter = 0
      ncyc = 0
      nflgc = 0
      nocon = 0
      nvar = 18
      ibcflg = 0
      vanleer=0.0
c
c.... set default array values

c
c     ---------------------------------------
      call setarry (uinf(1),0.0d0,4)
      call setarry (vinf(1),0.0d0,4)
      call setarry (pbc(1),0.0d0,4)

	do ii=1,nshp
	do iii=1,nobs(nshp)
	  oa1(ii,iii)=0.0d0
	  oa2(ii,iii)=0.0d0
	  ob1(ii,iii)=0.0d0
	  ob2(ii,iii)=0.0d0
	  oc1(ii,iii)=0.0d0
	  oc2(ii,iii)=0.0d0
	  od1(ii,iii)=0.0d0
	  od2(ii,iii)=0.0d0
	  oe1(ii,iii)=0.0d0
	  oe2(ii,iii)=0.0d0
	  nxo(ii,iii)=0.0d0
	  mxo(ii,iii)=0.0d0
	  nyo(ii,iii)=0.0d0
	  myo(ii,iii)=0.0d0
	enddo
	enddo

c      call setarry (oa1(1),0.0d0,nobd)
c      call setarry (oa2(1),0.0d0,nobd)
c      call setarry (ob1(1),0.0d0,nobd)
c      call setarry (ob2(1),0.0d0,nobd)
c      call setarry (oc1(1),0.0d0,nobd)
c      call setarry (oc2(1),0.0d0,nobd)
c      call setarry (od1(1),0.0d0,nobd)
c      call setarry (od2(1),0.0d0,nobd)
c      call setarry (oe1(1),0.0d0,nobd)
c      call setarry (oe2(1),0.0d0,nobd)
c      call setarry (nxo(1),0.0d0,nobd)
c      call setarry (mxo(1),0.0d0,nobd)
c      call setarry (nyo(1),0.0d0,nobd)
c      call setarry (myo(1),0.0d0,nobd)

      call setarry (fa1(1),0.0d0,nfrsrfd)
      call setarry (fa2(1),0.0d0,nfrsrfd)
      call setarry (fb1(1),0.0d0,nfrsrfd)
      call setarry (fb2(1),0.0d0,nfrsrfd)
      call setarry (fc1(1),0.0d0,nfrsrfd)
      call setarry (fc2(1),0.0d0,nfrsrfd)
      call setarry (fd1(1),0.0d0,nfrsrfd)
      call setarry (fd2(1),0.0d0,nfrsrfd)
      call setarry (fe1(1),0.0d0,nfrsrfd)
      call setarry (fe2(1),0.0d0,nfrsrfd)
      call setarry (nxf(1),0.0d0,nfrsrfd)
      call setarry (mxf(1),0.0d0,nfrsrfd)
      call setarry (nyf(1),0.0d0,nfrsrfd)
      call setarry (myf(1),0.0d0,nfrsrfd)
c     ---------------------------------------
c
c.... write out run time info
c
   60	continue

c	read(5,400) prbname
c      write(9,25)
c      write(13,25)
c      write(9,50) dat,tim,ochn
c      write(13,50) dat,tim,ochn
c      write(9,100) prbname
c      write(13,100) prbname
c      write(9,150)
c
c...  read initial input data; write it out to paper and film
c


      read (1,numparam)
      if (ndump.ne.0) then
        write(13,175)
      endif
      write(13,numparam)
      read (1,fldparam)
      write(13,fldparam)
      read (1,mesh)
      write(13,mesh)
      read (1,obstcl)
      write(13,obstcl)
      read (1,freesurf)
      write(13,freesurf) 
      read (1,wavecase)
      write(13,wavecase)

c --- read porous medium -- fyshi 10/18/02
      read (1,porousparam)
      write(13,porousparam)
      read (1,porous)
      write(13,porous)

c --- read turbulence -- fyshi 12/10/02
      read (1,smag)
      write(13,smag)

	read (1,turb)
	write(13,turb)

c --- porous parameter
      do ii=1,nportype
        ca(ii)=gammaa(ii)*(1.-porosity(ii))/porosity(ii)
      enddo

c
c.... Get problem name and print it out
c      if (cray) then
c        iotty=59
c       -----------------------------
c        open (unit=iotty,file="tty")
c       -----------------------------
c      else
c        iotty=6
c      endif
c      jnm=" RIPPLE "
c      write (iotty,500) prbname
c
c.... Compute constant terms and initialize necessary variables
c
      pi=acos(-1.0d0)
      if (icyl.eq.0) tpi=1.0d0
      cyl=float(icyl)
      omcyl=1.-cyl
      emf1=1.d00-emf
c.... contact angle
      dangle=cangler
      if (cyl.eq.1.0d0) canglel=90.0
      canglel=canglel*pi/180.
      cangler=cangler*pi/180.
      canglet=canglet*pi/180.
      cangleb=cangleb*pi/180.
c.... flux limiter
      if (con*1.3.le.fcvlim) go to 10
      write (13,200) con,fcvlim
      write (iotty,200) con,fcvlim
      fcvlim=1.3*con
   10 continue
c.... minimum time step
      dtend=0.002*delt
      if (dmpdt.eq.0.0) dmpdt=twfin
      twdmp=dmpdt
      if (deltemin.eq.0.0) deltemin=delt
      frsurf=frctn*rhof
      if (alpha.gt.1.0d0) then
        vanleer=1.0d0
        omalp=0.0
        opalp=2.0
      else
        omalp=1.-alpha
        opalp=1.+alpha
      endif
      if (kl.eq.5) pbcl=1.0d0
      if (kr.eq.5) pbcr=1.0d0
      if (kt.eq.5) pbct=1.0d0
      if (kb.eq.5) pbcb=1.0d0
      if ((kb.eq.5).or.(kt.eq.5).or.
     &     (kl.eq.5).or.(kr.eq.5)) sym=.false.

c ---- meshset is moved to here fyshi 12/12/02
c      from setup

       call meshset	

      xmu=xnu*rhof

	do j=2,jm1
	  do i=2,im1
          ij=(j-1)*imax+i
	    xmusmag(ij)=0.
	  enddo
	enddo



      if ((alpha.gt.1.0d0).and.(.not.(conserve))) alpha=1.0d0
c
      return
   25 format(32x," RIPPLE ")
   50 format(15x,"Date: ",a8,2x,"Time: ",a8,2x,"Machine: ",a8)
  100 format (/,5x,a80,/)
  150 format(/,15x,"* * * * * * * * * * * * * * * * * * * *",
     &       /,15x,"           ERROR DIAGNOSTICS           ",
     &       /,15x,"* * * * * * * * * * * * * * * * * * * *",/)
  175 format(/,15x,"* * * * * * * * * * * * * * * * * * * *",
     &       /,15x,"               RESTART                 ",
     &       /,15x,"* * * * * * * * * * * * * * * * * * * *",/)
  200 format (1x,4hcon=,1pe10.3,1x,11hand fcvlim=,e10.3,1x,
     &         41hare incompatible. setting fcvlim=1.3*con.)
  300 format (1x,13hbond number =,1pe12.4)
  400 format(a80)
  500 format (" RIPPLE:  ",a80,/," Processing input data . . .")
      end


 
      subroutine setarry (a,c,n)
c
c ======================================================================
c
c   Purpose -
c     sets the contents of array a of length n equal to c
c
c   SETARRY is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET      GEOM   INITREG   PLOTOUT    RINPUT
c           SETUP    STRAIN   TENSION    VOFADV    
c
c
c   SETARRY calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
c	include "32bit.h"
c##############################################################
c
      dimension a(n)
c
c  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 10 i=1,n
        a(i)=c
   10 continue
c
      return
      end


 
      subroutine meshset
c
c ======================================================================
c
c   Purpose -
c     mesh generator
c
c   MESHSET is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   MESHSET calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c           ISMIN   cftmath      KILL    ripple      GEOM    ripple
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
      include "comdk1.h"  	! qun
c############
c
      data one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      i=1
      j=1
      x(1)=xl(1)
      y(1)=yl(1)
c
      do 40 k=1,nkx
        if (nxl(k).eq.0) go to 20
        dxml=(xc(k)-xl(k))/nxl(k)
        nt=nxl(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dxmn1=dmin1(dxmn(k),dxml)
        cmc=(xc(k)-xl(k)-tn*dxmn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=xc(k)-xl(k)-cmc
        do 10 l=1,nt
          i=i+1
          rln=(float(l)-tn)/tn
   10   x(i)=xc(k)+bmc*rln-cmc*rln*rln
   20   if (nxr(k).eq.0) go to 40
        dxmr=(xl(k+1)-xc(k))/nxr(k)
        nt=nxr(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dxmn1=dmin1(dxmn(k),dxmr)
        cmc=(xl(k+1)-xc(k)-tn*dxmn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=xl(k+1)-xc(k)-cmc
        do 30 l=1,nt
          i=i+1
          rln=float(l)/tn
   30   x(i)=xc(k)+bmc*rln+cmc*rln*rln
   40 continue
c
      if (kr.ne.4) go to 50
      i=i+1
      x(i)=x(i-1)+x(2)-x(1)
   50 continue
c
      do 90 k=1,nky
        if (nyl(k).eq.0) go to 70
        dyml=(yc(k)-yl(k))/nyl(k)
        nt=nyl(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dymn1=dmin1(dymn(k),dyml)
        cmc=(yc(k)-yl(k)-tn*dymn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=yc(k)-yl(k)-cmc
        do 60 l=1,nt
          j=j+1
          rln=(float(l)-tn)/tn
   60   y(j)=yc(k)+bmc*rln-cmc*rln*rln
   70   if (nyr(k).eq.0) go to 90
        dymr=(yl(k+1)-yc(k))/nyr(k)
        nt=nyr(k)
        tn=nt
        tn=dmax1(tn,one+em6)
        dymn1=dmin1(dymn(k),dymr)
        cmc=(yl(k+1)-yc(k)-tn*dymn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=yl(k+1)-yc(k)-cmc
        do 80 l=1,nt
          j=j+1
          rln=float(l)/tn
   80   y(j)=yc(k)+bmc*rln+cmc*rln*rln
   90 continue
c
      if (kt.ne.4) go to 100
      j=j+1
      y(j)=y(j-1)+y(2)-y(1)
  100 continue
      numx=i
      numy=j
      numxm1=numx-1
      numym1=numy-1
      numxp1=numx+1
      numyp1=numy+1
      ibar=numx-1
      jbar=numy-1
      imax=ibar+2
      jmax=jbar+2
      im1=imax-1
      jm1=jmax-1
c
c.... calculate values needed for variable mesh
c
      do 120 i=1,numx
        if (x(i).eq.0.0) go to 110
        rx(i)=1.0/x(i)
        go to 120
  110   rx(i)=0.0
  120 continue
c
      do 130 i=2,numx
        xi(i)=0.5*(x(i-1)+x(i))
        delx(i)=x(i)-x(i-1)
  130 rdx(i)=1.0/delx(i)
c
      delx(1)=delx(2)
      xi(1)=xi(2)-delx(2)
      rdx(1)=1.0/delx(1)
      delxa=delx(numx)
      if (kr.eq.4) delxa=delx(3)
      delx(numxp1)=delxa
      xi(numxp1)=xi(numx)+delxa
      x(numxp1)=xi(numxp1)+0.5*delx(numxp1)
      rdx(numxp1)=1.0/delx(numxp1)
c
      do 140 i=2,numy
        yj(i)=0.5*(y(i-1)+y(i))
        dely(i)=y(i)-y(i-1)
        rdy(i)=1.0/dely(i)
  140 continue
c
      dely(1)=dely(2)
      rdy(1)=1.0/dely(1)
      yj(1)=yj(2)-dely(2)
      delya=dely(numy)
      if (kt.eq.4) delya=dely(3)
      dely(numyp1)=delya
      yj(numyp1)=yj(numy)+delya
      y(numyp1)=yj(numyp1)+0.5*dely(numyp1)
      rdy(numyp1)=1.0/dely(numyp1)
c
c.... set r and ri array for plane or cylindrical geometry
c
      do 145 i=1,imax
        r(i)=x(i)
        ri(i)=xi(i)
        if(icyl.eq.1) go to 145
        r(i)=1.0d0
        ri(i)=1.0d0
  145 continue
c
c.... set constant terms for plotting
c
      xmin=x(1)
      xmax=x(im1)
      ymin=y(1)
      ymax=y(jm1)
      xb(1)=xmin
      xb(2)=xmax
      xb(3)=xmax
      xb(4)=xmin
      xb(5)=xmin
      yb(1)=ymin
      yb(2)=ymin
      yb(3)=ymax
      yb(4)=ymax
      yb(5)=ymin
      if (iysymplt.gt.0) xmin=-xmax
      if (ixsymplt.gt.0) ymin=-ymax
      xf(1)=xmin
      xf(2)=xmax
      xf(3)=xmax
      xf(4)=xmin
      xf(5)=xmin
      yf(1)=ymin
      yf(2)=ymin
      yf(3)=ymax
      yf(4)=ymax
      yf(5)=ymin
      idxmin=ismin(im1-1,delx(2),1)
      dxmin=delx(idxmin)
      idymin=ismin(jm1-1,dely(2),1)
      dymin=dely(idymin)
c
c.... get the metric for the mesh
c
c     ----------
      call geom
c     ----------

      open(23,file='x.dat')
      do i=1,numxp1
         write(23,*)x(i)
      enddo
      close(23)

      open(23,file='y.dat')
      do i=1,numyp1
         write(23,*)y(i)
      enddo
      close(23)
      stop

c
c.... write out mesh information
c
      write (13,210)
      do 190 i=1,numxp1
        write (13,220) i,x(i),i,rx(i),i,delx(i),i,rdx(i),i,xi(i)
  190 continue
c
      write (13,210)
      do 200 i=1,numyp1
        write (13,230) i,y(i),i,dely(i),i,rdy(i),i,yj(i)
  200 continue
c
c.... test array size
c
      if (imax.le.ibar2.and.jmax.le.jbar2) go to 9999
c
c                    * * * * error section * * * *
c
      write (9,240)
c
c     ---------------------------------
c      call kill (iotty,ncyc,"MESHSET")   ! qun
	stop				  ! qun
c     ---------------------------------
c
 9999 return
c
  210 format (1h1)
  220 format (1x,2hx(,i3,2h)=,1pe12.5,2x,3hrx(,i3,2h)=,1pe12.5,2x,
     &        5hdelx(,i3,2h)=,1pe12.5,1x,4hrdx(,i3,2h)=,1pe12.5,2x,
     &        3hxi(,i3,2h)=,1pe12.5)
  230 format (1x,2hy(,i3,2h)=,1pe12.5,3x,5hdely(,i3,2h)=,1pe12.5,3x,
     &        4hrdy(,i3,2h)=,1pe12.5,3x,3hyj(,i3,2h)=,1pe12.5)
  240 format (45h mesh size inconsistent with array dimensions)
      end


 
      subroutine geom
c
c ======================================================================
c
c   Purpose -
c     calculate the metric for the mesh
c
c   GEOM is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         MESHSET
c
c
c   GEOM calls the following subroutines and functions -
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
c############
      include "comdk1.h"   ! qun
c############
c
      real*8 jacob
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... zero-out the volume arrays
c
c     ---------------------------------
      call setarry (cvol(1),0.0d0,nxy)
c     ---------------------------------
c
c.... Compute the metric coefficients used
c     in the implicit pressure solution
c
      do 100 j=2,jm1
        js=j-1
        do 100 i=2,im1
          is=i-1
          ij=(j-1)*imax+i
          x1=0.5*(x(is+1)+x(is+1)-x(is)-x(is))
          x2=0.5*(x(is+1)+x(is)-x(is)-x(is+1))
          y1=0.5*(y(js)+y(js+1)-y(js)-y(js+1))
          y2=0.5*(y(js+1)+y(js+1)-y(js)-y(js))
          jacob=x1*y2-x2*y1
c
          jcb(ij)=jacob
          alp(ij)=(x2**2 + y2**2)/jacob
          gam(ij)=(x1**2 + y1**2)/jacob
          cvol(ij)=ri(i)*delx(i)*dely(j)
c
  100 continue
c
      jb=1
      jt=jmax
      ijb=2
      ijt=jm1*imax+2
      do 101 i=2,im1
        alp(ijt)=alp(ijt-imax)
        gam(ijt)=gam(ijt-imax)
        jcb(ijt)=jcb(ijt-imax)
        cvol(ijt)=cvol(ijt-imax)
        alp(ijb)=alp(ijb+imax)
        gam(ijb)=gam(ijb+imax)
        jcb(ijb)=jcb(ijb+imax)
        cvol(ijb)=cvol(ijb+imax)
        ijb=ijb+1
        ijt=ijt+1
  101 continue
c
      il=1
      ir=imax
      ijl=1
      ijr=imax
      do 102 j=1,jmax
        alp(ijr)=alp(ijr-1)
        gam(ijr)=gam(ijr-1)
        jcb(ijr)=jcb(ijr-1)
        cvol(ijr)=cvol(ijr-1)
        alp(ijl)=alp(ijl+1)
        gam(ijl)=gam(ijl+1)
        jcb(ijl)=jcb(ijl+1)
        cvol(ijl)=cvol(ijl+1)
        ijl=ijl+imax
        ijr=ijr+imax
  102 continue
c
      return
      end


