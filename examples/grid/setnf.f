 
      subroutine setnf
c
c ======================================================================
c
c   Purpose -
c     determine the surface cell type nf(i,j)
c
c   SETNF is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP    VOFADV
c
c
c   SETNF calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c      
c
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
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... Set the default value of nf
c
       do 100 ij=1,nxy
        nf(ij)=0
  100  continue       
c
      do 751 i=2,im1
        do 750 j=2,jm1
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
c....     if cell is obstacle cell skip to end of loops
          if (ac(ij).lt.em6) go to 750
c
c....     declare empty cell to be a void cell
          if (f(ij).lt.emf) nf(ij)=6
c
c....     if cell is empty (or full but psat = 0.0) skip to end of loops;
c....     cells will have default values
          if (f(ij).lt.emf.or.(f(ij).gt.emf1.and.psat.eq.0.0))
     &       go to 750
c
c....     four tests to see whether one of the four neighbor cells is both
c         empty and open to flow from (i,j) cell; if so, enter main do
c         loops thru 190
          if (f(ipj).lt.emf.and.ar(ij).gt.em6) go to 190
          if (f(imj).lt.emf.and.ar(imj).gt.em6) go to 190
          if (f(ijp).lt.emf.and.at(ij).gt.em6) go to 190
          if (f(ijm).lt.emf.and.at(ijm).gt.em6) go to 190
c
c....     cell is not a surface cell, obstacle cell,empty cell or partic-
c         ular type of full cell. if it satisfies pressure test, set nf=5
c         for isolated cell and skip  to end of loops; otherwise cell is
c         fluid cell and we skip to end of loops with default values
          if (p(ij).le.psat*em6p1.and.psat.gt.0.0) nf(ij)=5
          go to 750
c
c....     we now enter calculational parts of main do loops
  190     continue
c
c....     calculate the partial derivatives of f
c
c....     distances from midpoint of cell to midpoint of neighbor cells
c         distance to right and left neighbors
          dxr=0.5*(delx(i)+delx(i+1))
          dxl=0.5*(delx(i)+delx(i-1))
c
c....     distance to top and bottom neighbors
          dyt=0.5*(dely(j)+dely(j+1))
          dyb=0.5*(dely(j)+dely(j-1))
c
c....     denominators for finite difference formulas for partial
c         derivatives in x and y directions
          rxden=1.0/(dxr*dxl*(dxr+dxl))
          ryden=1.0/(dyt*dyb*(dyt+dyb))
c
c....     fofm (and fofp) indicate whether cells with lesser (greater) in-
c         dices contribute to average fluid heights in three cell array
c         fofm=1.0 when cell contributes;=0.0 otherwise
c         index is j for vertical heights; i for horizontal heights
c         obstacle cell does not contribute if no fluid in neighbor cell of
c         three cell array
          fofm=1.0d0
          if (ac(ipjm).le.em6.and.f(ijm).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(ipjp).le.em6.and.f(ijp).lt.emf) fofp=0.0
c
c....     y fluid height in cells to right = avfr
c....     y heights measured from floor of ( j - 1 ) cells
          avfr=(1.0+ac(ipjm)*(f(ipjm)-1.0))*fofm*dely(j-1)+
     1         (1.0+ac(ipj)*(f(ipj)-1.0))*dely(j)+
     2         (1.0+ac(ipjp)*(f(ipjp)-1.0))*fofp*dely(j+1)
          fofm=1.0d0
          if (ac(imjm).le.em6.and.f(ijm).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(imjp).le.em6.and.f(ijp).lt.emf) fofp=0.0
c
c....     y fluid height in cells to left = avfl
          avfl=(1.0+ac(imjm)*(f(imjm)-1.0))*fofm*dely(j-1)+
     1         (1.0+ac(imj)*(f(imj)-1.0))*dely(j)+
     2         (1.0+ac(imjp)*(f(imjp)-1.0))*fofp*dely(j+1)
          fofm=1.0d0
          if (ac(ijm).le.em6.and.f(imjm)+f(ipjm).lt.emf)
     &       fofm=0.0
          fofp=1.0d0
          if (ac(ijp).le.em6.and.f(imjp)+f(ipjp).lt.emf)
     &       fofp=0.0
c
c....     y fluid height in central cells = avfcx
          avfcx=(1.0+ac(ijm)*(f(ijm)-1.0))*fofm*dely(j-1)+
     1          (1.0+ac(ij)*(f(ij)-1.0))*dely(j)+
     2          (1.0+ac(ijp)*(f(ijp)-1.0))*fofp*dely(j+1)
          fofm=1.0d0
          if (ac(imjp).le.em6.and.f(imj).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(ipjp).le.em6.and.f(ipj).lt.emf) fofp=0.0
c
c....     x fluid width in cells above = avft
c         x fluid width measured from floor of (i-1) cells
          avft=(1.0+ac(imjp)*(f(imjp)-1.0))*fofm*delx(i-1)+
     1         (1.0+ac(ijp)*(f(ijp)-1.0))*delx(i)+
     2         (1.0+ac(ipjp)*(f(ipjp)-1.0))*fofp*delx(i+1)
          fofm=1.0d0
          if (ac(imjm).le.em6.and.f(imj).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(ipjm).le.em6.and.f(ipj).lt.emf) fofp=0.0
c
c....     x fluid width in cells below = avfb
          avfb=(1.0+ac(imjm)*(f(imjm)-1.0))*fofm*delx(i-1)+
     1         (1.0+ac(ijm)*(f(ijm)-1.0))*delx(i)+
     2         (1.0+ac(ipjm)*(f(ipjm)-1.0))*fofp*delx(i+1)
          fofm=1.0d0
          if (ac(imj).le.em6.and.f(imjm)+f(imjp).lt.emf)
     &       fofm=0.0
          fofp=1.0d0
          if (ac(ipj).le.em6.and.f(ipjp)+f(ipjm).lt.emf)
     &    fofp=0.0
c
c....     x fluid width in central cells = avfcy
          avfcy=(1.0+ac(imj)*(f(imj)-1.0))*fofm*delx(i-1)+
     1          (1.0+ac(ij)*(f(ij)-1.0))*delx(i)+
     2          (1.0+ac(ipj)*(f(ipj)-1.0))*fofp*delx(i+1)
c
c....     avfl set by convention in first column of cells in x and y directions
c....     obstacles are placed on floor from which distances are measured;
c....     fluid is then above obstacles
c....     if nf = 2 or 4 floor is at top of cell; at i+1 or j+1 respectively
          if(i.eq.2) avfl=avfcx
c
c....     three point finite difference formulas for surface slopes with
c         variable mesh sizes; formulas exact for quadratics
c         slope dh/dx for almost horizontal fluid = pfx
          pfx=rxden*((avfr-avfcx)*dxl**2+(avfcx-avfl)*dxr**2)
c
c....     slope dw/dy for almost vertical fluid = pfy
          pfy=ryden*((avft-avfcy)*dyb**2+(avfcy-avfb)*dyt**2)
c
c....     pf = sum of squares of tangents; used as flag to differentiate
c         surface cells from isolated cells (nf=5)
          pf=pfx**2+pfy**2
c
c....     if pf very small, cell is declared isolated (instead of surface)
c         and we continue; otherwise go to 660
          if (pf.gt.em10) go to 660
c
c....     set nf(i,j) and p(i,j) for isolated cell; bypass determination
c         of pressure interpolation cell,calculation of surface pressure
  655     nf(ij)=5
c
c....     having set nf for the isolated cell,
c         we now skip to end of main loops
          go to 750
  660     continue
c
c....     for surface cells we pick up calculations of main loops; having
c         determined slopes and some auxiliary quantities,
c         we now get nf
c....     in order to set flags to be used later, we sum the f's in cols.
c         to right and left and in rows above and below the (i,j) cell
c
c         sfim = sum of f's in (i-1) col.
c         sfip = sum of f's in (i+1) col.
c         sfjp = sum of f's in (j+1) row
c         sfjm = sum of f's in (j-1) row
c         sfic = sum of f's in i col.
c         sfjc = sum of f's in j row
c
          sfim=f(imjp)+f(imj)+f(imjm)
          sfic=f(ijp)+f(ij)+f(ijm)
          sfip=f(ipjp)+f(ipj)+f(ipjm)
          sfjp=f(ipjp)+f(ijp)+f(imjp)
          sfjc=f(ipj)+f(ij)+f(imj)
          sfjm=f(ipjm)+f(ijm)+f(imjm)
c
c....     if there is little fluid in three by three array of cells: set
c         cell pressure as for an isolated cell and go to end of main loop
c         flags are initially set = 0
          if(sfim+sfic+sfip+sfjm+sfjc+sfjp.lt.0.10) go to 655
c
c....     if any cell face is completely closed to flow: skip row and
c         column tests
          if(ar(ij).lt.em6.or.ar(imj).lt.em6.or.at(ij).lt.em6.
     &         or.at(ijm).lt.em6) go to 670
          iflgx=0
          jflgy=0
c
c....     if either col. to left or col. to right is empty or is full:
c         iflgx = 1
          if (sfim.lt.emf.or.sfim.gt.3.0-emf) iflgx=1
          if (sfip.lt.emf.or.sfip.gt.3.0-emf) iflgx=1
c
c....     if either row above or row below is empty or is full: jflgy = 1
          if (sfjp.lt.emf.or.sfjp.gt.3.0-emf) jflgy=1
          if (sfjm.lt.emf.or.sfjm.gt.3.0-emf) jflgy=1
c
c....     if both flags = 1: continue execution at 670 without intervention
          if (iflgx.eq.1.and.jflgy.eq.1) go to 670
c
c....     if exactly one flag = 1: change the corresponding slope
          if (iflgx.eq.1) pfx=1.0d10*pfx
          if (jflgy.eq.1) pfy=1.0d10*pfy
  670     continue
c
c....     we have concluded slope increases
c
c....     determine the pressure interpolation cell nf
c
c         algorithm guarantees that a neighboring fluid cell (one or more
c         always exist) always lies at floor of the surface cell. the fluid
c         cell at the floor is used as the interpolation neighbor cell in
c         presit or prescr
c
c....     get absolute value of slopes; minimum will determine whether
c         surface has near horizontal or near vertical orientation
          abpfx=abs(pfx)
          abpfy=abs(pfy)
c
c....     set default values of the indices of the interpolation neighbor
c         cell (l,m): (l,m) = (i,j)
          l=i
          m=j
c
c....     if surface is more nearly horizontal: go to 680. otherwise
c         surface is more nearly vertical and we set nf = 2 or 1
          if (abpfy.ge.abpfx) go to 680
c
c....     set minimum slope, l index, and nf consistently with nf = 2
          pfmn=pfy
          nf(ij)=2
          l=i+1
c
c....     compute length intervals needed for evaluating the pressure
c         interpolation
          dmx=delx(i)
          dmin=0.5*(dmx+delx(i+1))
          dnbr=delx(i+1)
c
c....     Sign of larger absolute slope (pfx) determines nf value.
c         if pfx > 0.0 then nf = 2 and we skip to end of nf routine at 690
c         otherwise nf = 1 and we go on, repeating immediately previous
c         steps appropriately for nf = 1
          if (pfx.gt.0.0) go to 690
          nf(ij)=1
          pfmn=-pfy
          l=i-1
          dmx=delx(i)
          dmin=0.5*(dmx+delx(i-1))
          dnbr=delx(i-1)
c
c....     having completed calculations for nf = 1, we skip to end of nf
c         routine holding nf = 1 data
          go to 690
  680     continue
c
c....     we are now in surface more nearly horizontal case.
c         we start with nf = 4
c         set minimum slope, l index and nf consistently with nf = 4
          pfmn=-pfx
          nf(ij)=4
          m=j+1
          dmx=dely(j)
          dmin=0.5*(dmx+dely(j+1))
          dnbr=dely(j+1)
c
c....     Sign of larger slope (pfy) determines nf value.
c         if pfy > 0.0 then nf = 4 and we skip to end of nf routine
c         otherwise nf = 3 and we go on, repeating immediately previous
c         steps appropriately for nf = 3
          if (pfy.gt.0.0) go to 690
          nf(ij)=3
          pfmn=pfx
          m=j-1
          dmx=dely(j)
          dmin=0.5*(dmx+dely(j-1))
          dnbr=dely(j-1)
  690     continue
c
  750   continue
  751 continue
c
      return
      end


      subroutine surface
c
c ======================================================================
c
c   Purpose - integrate surface elevation
c
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
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
	if(ncyc.eq.0) then
	do i=1,imax
	etah(i)=flht
	end do
	goto 99
	end if
	do 200 i=1,imax
	do 100 j=jm1,2,-1
      ij=(j-1)*imax+i
      ijm=ij-imax
	if(f(ij).lt.0.5.and.f(ijm).gt.0.5) then
	iy=j-1
	goto 150
	end if
100	continue
      ij2=imax+i
	if (f(ij2).lt.emf) etah(i)=0.
      ijm1=(jm1-1)*imax+i
	if (f(ijm1).ge.0.5) etah(i)=y(jm1)+0.5*dely(jm1)*f(ijm1)
	goto 200
150	continue
	ijiy=(iy-1)*imax+i
	ijiyp=ijiy+imax 
	yeta=dely(iy)*(f(ijiy)-0.5)+dely(iy+1)*f(ijiyp)
	etah(i)=yj(iy)+dmax1(0.d0,yeta)
200	continue
99	continue

	do i=1,imax
	etahn(i)=etah(i)
	end do

	return
	end



 
	subroutine wave_condition 
c
c ======================================================================
c
c   Purpose - calculate wave condition
c	      Qun Zhao Oct. 12, 2002 
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

	wavew=2.*pi/wavet
	if(iwave.eq.2.or.iwave.eq.3) then
        call  wavenum(pi,-gy,flht,wavew,wavek,xkn1)
	wavel=2.*pi/wavek
	wavec=wavel/wavet
	write(*,*) '	L=',wavel,'		C=',wavec, '	k=',wavek
	write(*,*) '  D/L=',flht/wavel, '	H/D=',waveh/flht
	end if

	return 
	end


