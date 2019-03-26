 
      subroutine implctp
c
c ======================================================================
c
c   Purpose -
c     calculate the implicit pressure
c
c   IMPLCTP is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE   SRFEMIN
c
c
c   IMPLCTP calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c           ACCEL    ripple    ILUCGJ      iccg   BDYCELL    ripple
c           CGITJ      iccg    STRAIN    ripple
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
       include "iccgdk.h" 
c###########
c
      data tiny /1.0d-25/, zerod /0.0d0/, oned /1.0d0/
      data isw /0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... get the strain and rotation tensors
c
c     -----------------------------------------------------------------
      call strain(ar,at,ac,x,y,ri,r,u,v,div,curlu,im1,jm1,nxy,
     &            tauxx,tauyy,tauxy)
c     -----------------------------------------------------------------

c ---   set boundary condition of div for porous media -fyshi 10/24	

c
c.... reflect pressures to ghost cells
c
c     --------------------------------------------------------
      call bdycell(im1,jm1,1,1,1,1,pbcl,pbcr,pbct,pbcb,pbc,p)
c     --------------------------------------------------------
c
c.... load the iccg arrays, which are the diagonals of
c     the matrix to be inverted.  The arrays represent
c     the coefficients for cell ij and each of its 8
c     neighbors, as shown below
c
c                       ----------------------------
c                       |        |        |        |
c                       |        |        |        |
c                       |  imjp  |  ijp   |  ipjp  |
c                       |        |        |        |
c                       |  (C)   |  (D)   |  (E)   |
c                       ----------------------------
c                       |        |        |        |
c                       |        |        |        |
c                       |  imj   |   ij   |  ipj   |
c                       |        |        |        |
c                       |  (BF)  |  (A)   |  (B)   |
c                       ----------------------------
c                       |        |        |        |
c                       |        |        |        |
c                       | imjm   |  ijm   |  ipjm  |
c                       |        |        |        |
c                       |  (EF)  |  (DF)  |  (CF)  |
c                       ----------------------------
c
c ... divs at porous media interface are taken out fyshi 10/31/02
c
c
      do 751 ii=1,nportype	
      do 750 i=2,im1
        do 750 j=2,jm1
          ij=(j-1)*imax+i
         if(pbdy(ii,ij).eq.1.)div(ij)=div(ij+imax)	
         if(pbdy(ii,ij).eq.3.)then
	   if(j.ne.2) div(ij)=div(ij+imax)	   
c	   div(ij+imax)=div(ij+2*imax)
	 endif
         if(pbdy(ii,ij).eq.2.)div(ij)=div(ij-1)
         if(pbdy(ii,ij).eq.4.)div(ij)=div(ij+1)
         if(pbdy(ii,ij).eq.5.)div(ij)=div(ij+1)
         if(pbdy(ii,ij).eq.6.)div(ij)=div(ij-1)
         if(pbdy(ii,ij).eq.7.)div(ij)=div(ij-imax-1)
         if(pbdy(ii,ij).eq.8.)div(ij)=div(ij-imax+1)		 		 
750   	 continue
751	 continue      

c --- check div
c	open(39,file='tmp.dat')
c	do  j=1,jmax
c	write(39,99) (div(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
c	open(39,file='tmp1.dat')
c	do  j=1,jmax
c	write(39,99) (pbdy(1,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
c	end do	
c	close(39)
c 99	format(3x,800(f5.1)) 
c	pause

      nm=1
c      coef=4.*xmu/3. ! can not find usage - fyshi

c --- add porous                                 - fyshi 10/18/02
c

      do 4001 j=2,jm1
        do 4000 i=2,im1
          ij=(j-1)*imax+i
          ipj=ij+1
          ijp=ij+imax
          ipjp=ipj+imax
          imjp=ijp-1
          imj=ij-1
          imjm=imj-imax
          ijm=ij-imax
          ipjm=ipj-imax
          if ((i.eq.2).and.(cyl.eq.1.0d0)) imjm=ijm
          if ((i.eq.2).and.(cyl.eq.1.0d0)) imj=ij
          if ((i.eq.2).and.(cyl.eq.1.0d0)) imjp=ijp
          p(ij)=pn(ij)-psat
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhoimj=f(imj)*rhof
          rhoijm=f(ijm)*rhof
          zero=cvmgt(zerod,oned,
     &                (rhoij.lt.frsurf).or.(ac(ij).lt.emf))
c
          rl=omcyl+cyl*0.5*(x(i-1)+x(i-1))
          rhol=(delx(i-1)*rhoij+delx(i)*rhoimj)/(delx(i)+delx(i-1))
          rrhol=1./(rhol+tiny)
          rr=omcyl+cyl*0.5*(x(i)+x(i))
          rhor=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rrhor=1./(rhor+tiny)
          rt=omcyl+cyl*0.5*(x(i-1)+x(i))
          rhot=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          rrhot=1./(rhot+tiny)
          rb=omcyl+cyl*0.5*(x(i-1)+x(i))
          rhob=(dely(j-1)*rhoij+dely(j)*rhoijm)/(dely(j)+dely(j-1))
          rrhob=1./(rhob+tiny)
c
          alfr=0.5*ar(ij)*rr*rrhor*(alp(ipj)+alp(ij))
          alfl=0.5*ar(imj)*rl*rrhol*(alp(imj)+alp(ij))
          betr=0.0
          betl=0.0
          bett=0.0
          betb=0.0
          gamt=0.5*at(ij)*rt*rrhot*(gam(ijp)+gam(ij))
          gamb=0.5*at(ijm)*rb*rrhob*(gam(ijm)+gam(ij))
c
          a(nm)=(alfr+alfl+gamb+gamt)
          b(nm)=-(alfr+0.25*(betb-bett))
          bf(nm)=-(alfl+0.25*(bett-betb))
          c(nm)=-0.25*(betl+bett)
          cf(nm)=-0.25*(betr+betb)
          d(nm)=-(gamt+0.25*(betl-betr))
          df(nm)=-(gamb+0.25*(betr-betl))
          e(nm)=0.25*(betr+bett)
          ef(nm)=0.25*(betl+betb)
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pc(ii,ij).eq.1)then
          srce(nm)=-ac(ij)*cvol(ij)*div(ij)/delt
     &             *(1.+ca(ii))/porosity(ii)
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          srce(nm)=-ac(ij)*cvol(ij)*div(ij)/delt
	  endif

          soln(nm)=p(ij)
c
          a(nm)=cvmgt(oned/tiny,a(nm),ac(ij).lt.emf)
          b(nm)=zero*b(nm)
          bf(nm)=zero*bf(nm)
          c(nm)=zero*c(nm)
          cf(nm)=zero*cf(nm)
          d(nm)=zero*d(nm)
          df(nm)=zero*df(nm)
          e(nm)=zero*e(nm)
          ef(nm)=zero*ef(nm)
          srce(nm)=zero*srce(nm)
          soln(nm)=zero*soln(nm)
c
          nm=nm+1
c
 4000   continue
 4001 continue
c
c.... Boundary conditions
c
c.... Internal obstacle boundaries
c     (enforce Neumann conditions)
c
      do 400 j=2,jm1
        do 400 i=2,im1
c
          ij=(j-1)*imax+i
          if (ac(ij).lt.em6) go to 400
          nm=(j-2)*ibar+i-1
          ipj=ij+1
          ijp=ij+imax
          ijm=ij-imax
          imj=ij-1
c
          if (ac(ijp).lt.em6.and.j.lt.jm1) then
            a(nm)=a(nm)+d(nm)
            d(nm)=0.0
          endif
          if (ac(ijm).lt.em6.and.j.gt.2) then
            a(nm)=a(nm)+df(nm)
            df(nm)=0.0
          endif
          if (ac(imj).lt.em6.and.i.gt.2) then
            a(nm)=a(nm)+bf(nm)
            bf(nm)=0.0
          endif
          if (ac(ipj).lt.em6.and.i.lt.im1) then
            a(nm)=a(nm)+b(nm)
            b(nm)=0.0
          endif
c
  400 continue
c
c.... Left boundary
c
      if (kl.eq.4) go to 115
      wsl=0.0
      if (kl.eq.5) wsl=1.0d0
      imj=imax+1
      imjm=imj-imax
      imjp=imj+imax
      nm=1
      do 110 j=2,jm1
        a(nm)=a(nm)+(1.-wsl)*bf(nm)
        d(nm)=d(nm)+(1-wsl)*c(nm)
        df(nm)=df(nm)+(1-wsl)*ef(nm)
        srce(nm)=srce(nm)
     &       -wsl*(c(nm)*p(imjp)+bf(nm)*p(imj)
     &            +ef(nm)*p(imjm))
        bf(nm)=0.0
        ef(nm)=0.0
        c(nm)=0.0
        imj=imj+imax
        imjm=imjm+imax
        imjp=imjp+imax
        nm=nm+im1-1
  110 continue
c
c.... Top boundary
c
  115 if (kt.eq.4) go to 125
      wst=0.0
      if (kt.eq.5) wst=1.0d0
      ij=jm1*imax+2
      imj=ij-1
      ipj=ij+1
      nm=(jm1-2)*(im1-1)+1
      do 120 i=2,im1
        a(nm)=a(nm)+(1.-wst)*d(nm)
        bf(nm)=bf(nm)+(1-wst)*c(nm)
        b(nm)=b(nm)+(1-wst)*e(nm)
        srce(nm)=srce(nm)-
     &           wst*(c(nm)*p(imj)+d(nm)*p(ij)
     &           +e(nm)*p(ipj))
        c(nm)=0.0
        d(nm)=0.0
        e(nm)=0.0
        nm=nm+1
        imj=imj+1
        ij=ij+1
        ipj=ipj+1
  120 continue
c
c.... Right boundary
c
  125 if (kr.eq.4) go to 135
      wsr=0.0
      if (kr.eq.5) wsr=1.0d0
      nm=im1-1
      ij=2*imax
      ijp=ij+imax
      ijm=ij-imax
      do 130 j=2,jm1
        a(nm)=a(nm)+(1.-wsr)*b(nm)
        d(nm)=d(nm)+(1-wsr)*e(nm)
        df(nm)=df(nm)+(1-wsr)*cf(nm)
        srce(nm)=srce(nm)-wsr*(b(nm)*p(ij)+e(nm)*p(ijp)
     &           +cf(nm)*p(ijm))
        b(nm)=0.0
        e(nm)=0.0
        cf(nm)=0.0
        ij=ij+imax
        ijm=ijm+imax
        ijp=ijp+imax
        nm=nm+im1-1
  130 continue
c
c.... Bottom boundary
c
  135 if (kb.eq.4) go to 145
      wsb=0.0
      if (kb.eq.5) wsb=1.0d0
      ijm=2
      imjm=1
      ipjm=3
      nm=1
      do 140 i=2,im1
        a(nm)=a(nm)+(1.-wsb)*df(nm)
        bf(nm)=bf(nm)+(1-wsb)*ef(nm)
        b(nm)=b(nm)+(1-wsb)*cf(nm)
        srce(nm)=srce(nm)-wsb*
     &           (ef(nm)*p(imjm)+df(nm)*p(ijm)
     &          +cf(nm)*p(ipjm))
        ef(nm)=0.0
        df(nm)=0.0
        cf(nm)=0.0
        nm=nm+1
        ijm=ijm+1
        imjm=imjm+1
        ipjm=ipjm+1
  140 continue
c
  145 continue
      itmax=itmxiccg
      its=itmxiccg
      error=erriccg
      nxm=ibar
      nym=jbar
c
      if (sym) then

c       ------------------------------------------------------------
        call cgitj(a,b,c,d,e,af,bf,cf,df,ef,nxm,nym,soln,srce,residu,
     &              sc1,sc2,its,error,xratio,yratio,xnorm,ynorm)
c       ------------------------------------------------------------
      else
c       ------------------------------------------------------------
        call ilucgj(ef,df,cf,bf,a,b,c,d,e,eh,dh,ch,bh,ag,bg,cg,dg,
     &              eg,nxm,nym,soln,srce,residu,sc1,sc2,its,error,
c     &              itmax,xratio,yratio,xnorm,ynorm)
c -- fyshi, there is no itmax in ilucgj subroutine, remove it now

     &              xratio,yratio,xnorm,ynorm)
c       ------------------------------------------------------------
      endif
c
      iter=its
      if (iter .gt. itmax) then
        write(iotty,150) ncyc,t,its
        write(13,150) ncyc,t,its
        nocon=nocon+1
      endif
c
      nm=1
      do 200 j=2,jm1
        do 200 i=2,im1
          ij=(j-1)*imax+i
          p(ij)=soln(nm)+psat
          nm=nm+1
  200 continue
c
c.... Set pressures in the ghost cells
c
c     --------------------------------------------------------
      call bdycell(im1,jm1,1,1,1,1,pbcl,pbcr,pbct,pbcb,pbc,p)
c     --------------------------------------------------------
c
      if (nobshp.eq.0) go to 9999
c
c.... Enforce Neumann conditions across internal
c     obstacle boundaries by quadratically
c     interpolating an AC-weighted pressure
c
      do 300 j=2,jm1
        do 300 i=2,im1
          ij=(j-1)*imax+i
          if (ac(ij).gt.em6) go to 300
c
          ipj=ij+1
          ijp=ij+imax
          ijm=ij-imax
          imj=ij-1
          imjm=ij-1-imax
          imjp=ij-1+imax
          ipjm=ij-imax+1
          ipjp=ij+imax+1
c
          sumac=9.*ac(ij)/16. +
     &          3.*(ac(ipj)+ac(imj)+ac(ijp)+ac(ijm))/32. +
     &          1.*(ac(ipjp)+ac(imjp)+ac(ipjm)+ac(imjm))/64.
          sumpac=9.*ac(ij)*p(ij)/16. +
     &           3.*(ac(ipj)*p(ipj)+ac(imj)*p(imj)+
     &               ac(ijp)*p(ijp)+ac(ijm)*p(ijm))/32. +
     &           1.*(ac(ipjp)*p(ipjp)+ac(imjp)*p(imjp)+
     &               ac(ipjm)*p(ipjm)+ac(imjm)*p(imjm))/64.
          p(ij)=sumac*sumpac/(sumac+tiny)**2
c
  300 continue
c     --------------------------------------------------------
      call bdycell(im1,jm1,1,1,1,1,pbcl,pbcr,pbct,pbcb,pbc,p)
c     --------------------------------------------------------
c
c     -----------
 9999 call accel
c     -----------
c
      return
  150 format("IMPLCTP pressure solution for cycle ",i5," , time ",
     &        1pe12.5,/,"did not converge after ",i5," iterations")
      end



C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C		SUBROUTINES
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


 
      subroutine bdycell(nx,ny,wl,wr,wt,wb,bcl,bcr,bct,bcb,
     &                    bound,p)
c
c ======================================================================
c
c   Purpose -
c     fills the ghost cells with the boundary quantities
c
c   BDYCELL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         IMPLCTP   TENSION
c
c
c   BDYCELL calls the following subroutines and functions -
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
      dimension p(1),bound(4)
      integer wl,wr,wt,wb
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      nxp=nx+1
      nyp=ny+1
      ijt=ny*nxp+2-nxp
      ijb=2+nxp
      ijl=2
      ijr=nx
 
c.... bottom boundary
c
      ij=2
      do 300 i=2,nx
        p(ij)=bcb*bound(1)+(1.-bcb)*p(ij+nx+1)
        ij=ij+1
        ijt=ijt+1
  300 continue
c
c.... top boundary
c
      ij=ny*(nx+1)+2
      do 400 i=2,nx
        p(ij)=bct*bound(2)+(1.-bct)*p(ij-nx-1)
        ij=ij+1
        ijb=ijb+1
  400 continue
c
c.... left boundary
c
      ij=1
      do 100 j=1,ny+1
        p(ij)=bcl*bound(3)+(1.-bcl)*p(ij+1)
        ij=ij+nx+1
        ijr=ijr+nx+1
  100 continue
c
c.... right boundary
c
      ij=nx+1
      do 200 j=1,ny+1
        p(ij)=bcr*bound(4)+(1.-bcr)*p(ij-1)
        ij=ij+nx+1
        ijl=ijl+nx+1
  200 continue
c
      return
      end


      subroutine cgitj(a,b,c,d,e,af,bf,cf,df,ef,
     1     n,m,x,y,r,up,aux,istop,ercg,
     2     xratio,ratio,xnorm,ynorm)
c
c***purpose
c
c      given a 9-point matrix multiply routine called matmul9
c      and an approximate factorization a = l*d*lt, this routine
c      provides a conjugate gradient acceleration to iterate
c      to an approximate solution to a*x=y.
c***description
c
c       cgitj requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. five diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the ldlt factors are stored similarly in five
c      lower diagonal arrays.
c
c     on entry
c
c        a,b,c,d,e   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the upper diagonals of a as one sweeps a left to right.
c
c        el,dl,cl,bl,af   are 1-d arrays of length
c            at least nxm. the arrays contain from left to right
c            the lower diagonals of l,d as they are swept left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is an estimate to the solution vector of size nxm.
c
c        y     is an r.h.s. vector of dimenison at least nxm.
c
c        r,up,aux  are work variables at least nxm in length.
c               the final residual is left array r.
c
c        istop  is the maximum number of iterations allowed.
c               upon exit it contains the actual number of
c               iterations required to converge to accuracy
c               specified.
c
c        ercg   is the error tolerance for convegence
c
c     on return
c
c        x     is the solution vector of dimension at least
c              nxm containing the solution vector.
c
c        xratio   is the ratio of dxnorm to xnorm.
c
c        ratio    is the ratio of rnorm to ynorm
c
c        xnorm    is the l2 norm of the solution vector x.
c
c        ynorm    is the l2 norm of the r.h.s.
c
c-------------------------------------------------------
c***routines called  ldlt9f matmul9 (ldlt)(JTK,6/17)
c-------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 e(m*n),d(m*n),c(m*n),b(m*n),a(m*n),ef(m*n),
     &       df(m*n),cf(m*n),
     &       bf(m*n),af(m*n),x(m*n),y(m*n),r(m*n),up(m*n),aux(m*n),
     &	   el(m*n),dl(m*n),cl(m*n),bl(m*n)
c  el,dl,cl,bl is added by Seo
      mn = m*n
      e3=ercg
c

c --- the linux f77 compiler can not make b(0) with the defination b(m*n)
c      so use nfun instead of 0 

	nfun=0

ccc copy a matrix transposed into factor matrix
c
      do  3 i=1,n+1
        ef(i) = 0.
        df(i) = 0.
        cf(i) = 0.
    3 continue
      do 5 i=1,mn-n-1
        ef(i+n+1) = e(i)
        cf(i+n) = c(i+1)
    5 continue
      do 10 i= 1,mn-n
        df(i+n) = d(i)
   10 continue
      do 15 i=1,mn-1
        bf(i+1) = b(i)
   15 continue
      do 20 i=1,mn
        af(i) = a(i)
   20 continue
	
	do i=1,mn-1
		bl(i+1) = b(i)
	enddo
	do i=1,mn-n+1
		cl(i+n-1) = c(i)
	enddo
	do i=1,mn-n
		dl(i+n) = d(i)
	enddo
	do i=1,mn-n-1
		el(i+n+1) = e(i)
	enddo
c
ccc periodic block copy
c
      do 22 i=2,n
        ef(i) = e(mn-n+i-1)
   22 continue
      do 25 i=1,n
        df(i) = d(mn-n+i)
   25 continue
      do 28 i=1,n-1
        cf(i) = c(mn-n+i+1)
   28 continue

c-------------------------------------------------------
      call ldlt9f(n,m,ef,df,cf,bf,af)
c-------------------------------------------------------
c
c initialize the residual r

c-------------------------------------------------------
c      call matmul9(e(-n),d(1-n),c(2-n),b(nfun),a,b,c,d,e,n,m,x,r)
      call matmul9(el,dl,cl,bl,a,b,c,d,e,n,m,x,r)
c-------------------------------------------------------

      do 30 i=1,mn
        r(i) = y(i)-r(i)
   30 continue
      if (dasum(mn,r,1) .eq. 0.)     then
        i = 0
        go to 200
      end if
      ynorm = dasum(mn,y,1)
c
ccc solve l * d * l(trans)*aux = r0
c
c-------------------------------------------------------
      call ldlt(n,m,ef,df,cf,bf,af,bf(2),cf(n),df(n+1),ef(n+2),
     1          up,r)
c-------------------------------------------------------
c
ccc dot r with aux
c
      rdot = ddot(mn,r,1,up,1)
c
c put up into aux (p0)
c
      do 40 i4=1,mn
 40     aux(i4) = up(i4)
c p is stored in aux and mp or (ldlt)-1(r) are in up.
c
ccc begin main loop
c
      do 100 i = 1,istop
c
c---- find m * p and store it in up
c

c-------------------------------------------------------
c     call matmul9(e(-n),d(1-n),c(2-n),b(nfun),a,b,c,d,e,n,m,aux,up)
      call matmul9(el,dl,cl,bl,a,b,c,d,e,n,m,aux,up)
c-------------------------------------------------------

c
ccc compute alpha(i)
c
      alpha = ddot(mn,up,1,aux,1)
      alpha = rdot/alpha
c
c
ccc compute x(i+1)
c
      do 80 i8 = 1,mn
 80     x(i8) = x(i8) + alpha*aux(i8)
c
ccc compute r(i+1)
c
      do 90 i9=1,mn
        r(i9) = r(i9) - alpha*up(i9)
 90   continue
c
ccc termination test
ccc
ccc stop whenever norm(aux) / norm (x) .lt. e3
c
      rnorm = dasum(mn,r,1)
      dxnorm=dasum(mn,aux,1)*dabs(alpha)
      xnorm=dasum(mn,x,1)
c***
ccc output ratio
c
      ratio = rnorm/ynorm
      xratio=dxnorm/xnorm
c
c***
      if ( xratio .le. e3.and.ratio.le.e3 ) go to 200
c
ccc compute beta(i)
c
c-------------------------------------------------------
      call ldlt(n,m,ef,df,cf,bf,af,bf(2),cf(n),df(n+1),ef(n+2),
     1          up,r)
c-------------------------------------------------------
      alpha = ddot(mn,r,1,up,1)
      beta = alpha/rdot
      rdot = alpha
c get new p (stored in aux)
      do 95  i11=1,mn
 95     aux(i11)=beta*aux(i11)+up(i11)
  100 continue
 102  format (1h0,'cg terminates because of too many iterations. ', /,
     x         'xratio = ',e12.4, 'ratio = ',e12.4)
  200 istop=i

      return
      end


      subroutine ldlt9f(n,m,e,d,c,b,a)
c
c***purpose
c
c      factor (incomplete) a=l*d*lt, where the matrix a is
c      a 5-diagonal matrix associated with a 9-point difference
c      s.p.d, operator.
c
c***description
c
c       ldlt9f requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. five diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the lower diagonals of a as one sweeps a left to right.
c
c        n    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c     on return
c
c        e,d,c,b,a    contain the factors
c            ldlt. the reciprocal of d is stored in array a.
c
c***routines called   none
c
c*** first executable statement   ldlt9f
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m)
c
ccc   form a = l*v
c
      do 100 j=1,m
c
ccc factor central block
c
        a(1,j) = 1./a(1,j)
        do 5 i=2,n
          a(i,j) = 1./(a(i,j)-b(i,j)*(a(i-1,j)*b(i,j)))
    5   continue
        if (j .ne. m)    then
c
ccc   form lower diagonals
c
        do 10 i=2,n
          d(i,j+1) = d(i,j+1)-(a(i-1,j)*e(i,j+1))*b(i,j)
   10   continue
        do 15 i=1,n-1
          c(i,j+1) = c(i,j+1)-(a(i,j)*d(i,j+1))*b(i+1,j)
   15   continue
c
ccc factor lower periodic diagonals
c
        if (j .eq. 1)    then
        do 11 i=2,n
          d(i,1) = d(i,1)-(a(i-1,1)*e(i,1))*b(i,1)
   11   continue
        do 16 i=1,n-1
          c(i,1) = c(i,1)-(a(i,1)*d(i,1))*b(i+1,1)
   16   continue
        end if
c
ccc   form central diagonals
c
        do 20 i=2,n
          b(i,j+1) = b(i,j+1)-(a(i-1,j)*e(i,j+1))*d(i-1,j+1)
     1            -(a(i,j)*d(i,j+1))*c(i-1,j+1)
   20   continue
          a(1,j+1) = a(1,j+1)-(a(1,j)*d(1,j+1))*d(1,j+1)
     1             -(a(2,j)*c(1,j+1))*c(1,j+1)
        do 30 i=2,n-1
          a(i,j+1) = a(i,j+1)-(a(i-1,j)*e(i,j+1))*e(i,j+1)
     1            -(a(i,j)*d(i,j+1))*d(i,j+1)
     2            -(a(i+1,j)*c(i,j+1))*c(i,j+1)
   30   continue
          a(n,j+1) = a(n,j+1)-(a(n-1,j)*e(n,j+1))*e(n,j+1)
     1            -(a(n,j)*d(n,j+1))*d(n,j+1)
c
ccc factor central periodic block
c
        if (j .eq. 1)    then
        do 21 i=2,n
          b(i,m) = b(i,m)-(a(i-1,1)*e(i,1))*d(i-1,1)
     1            -(a(i,1)*d(i,1))*c(i-1,1)
   21   continue
          a(1,m) = a(1,m)-(a(1,1)*d(1,1))*d(1,1)
     1             -(a(2,1)*c(1,1))*c(1,1)
        do 31 i=2,n-1
          a(i,m) = a(i,m)-(a(i-1,1)*e(i,1))*e(i,1)
     1            -(a(i,1)*d(i,1))*d(i,1)
     2            -(a(i+1,1)*c(i,1))*c(i,1)
   31   continue
          a(n,m) = a(n,m)-(a(n-1,1)*e(n,1))*e(n,1)
     1            -(a(n,1)*d(n,1))*d(n,1)
        end if
c
ccc scale diagonals
c
          do 40 i=2,n
          e(i,j+1) = (a(i-1,j)*e(i,j+1))
   40   continue
        do 50 i=1,n
          d(i,j+1) = (a(i,j)*d(i,j+1))
   50   continue
        do 60 i=1,n-1
          c(i,j+1) = (a(i+1,j)*c(i,j+1))
   60   continue
c
ccc scale periodic block
c
        if (j .eq. 1)    then
          do 41 i=2,n
          e(i,1) = (a(i-1,1)*e(i,1))
   41   continue
        do 51 i=1,n
          d(i,1) = (a(i,1)*d(i,1))
   51   continue
        do 61 i=1,n-1
          c(i,1) = (a(i+1,1)*c(i,1))
   61   continue
        end if
        end if
        do 70 i=2,n
          b(i,j) = a(i-1,j)*b(i,j)
   70 continue
  100 continue
      return
      end



      subroutine ldlt(n,m,e,d,c,b,a,bu,cu,du,eu,x,y)
c
c***purpose
c
c      solve (incomplete) ax=(l*d*u)x=y, where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator. factors are generated lu9p.
c
c***description
c
c       ldlt requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        y     is an input vector of dimenisonat least nxm.
c
c     on return
c
c        x     is the output vector of dimension at least
c              nxm containing the solution vector.
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m),
     1          bu(n,m),cu(n,m),du(n,m),eu(n,m),x(n,m),y(n,m)
c*** first executable statement   ldlt
c
ccc   copy y into x
c
      do 1 j=1,m
        do 1 i=1,n
          x(i,j) = y(i,j)
    1 continue
c
ccc   forward solve
c
      do 100 j=1,m
        if (j .ne. 1)    then
          do 15 i=1,n-1
            x(i,j) = x(i,j)-c(i,j)*x(i+1,j-1)
   15     continue
          do 20 i=1,n
            x(i,j) = x(i,j)-d(i,j)*x(i,j-1)
   20     continue
          do 25 i=2,n
            x(i,j) = x(i,j)-e(i,j)*x(i-1,j-1)
   25     continue
c
ccc vertical periodic equations
c
        if (j .eq. m)    then
          do 16 i=1,n-1
            x(i,m) = x(i,m)-c(i,1)*x(i+1,1)
   16     continue
          do 21 i=1,n
            x(i,m) = x(i,m)-d(i,1)*x(i,1)
   21     continue
          do 26 i=2,n
            x(i,m) = x(i,m)-e(i,1)*x(i-1,1)
   26     continue
        end if
        end if
c
ccc   van der worst approximation
c
cdir$ ivdep
        do 5 i=n,2,-1
          x(i,j) = x(i,j)-b(i,j)*x(i-1,j)
    5   continue
cdir$ ivdep
        do 10 i=n,3,-1
          x(i,j) = x(i,j)+b(i-1,j)*b(i,j)*x(i-2,j)
   10   continue
  100 continue
c
ccc   diagonal solve
c
      do 30 j=1,m
        do 30 i=1,n
          x(i,j) = a(i,j)*x(i,j)
   30 continue
c
ccc  backward solve
c
      do 200 j=m,1,-1
        if (j .ne. m)    then
          do 45 i=2,n
            x(i,j) = x(i,j)-cu(i,j)*x(i-1,j+1)
   45     continue
          do 50 i=1,n
            x(i,j) = x(i,j)-du(i,j)*x(i,j+1)
   50     continue
          do 55 i=1,n-1
            x(i,j) = x(i,j)-eu(i,j)*x(i+1,j+1)
   55     continue
c
ccc vertical periodic equations
c
        if (j .eq. 1)   then
          do 46 i=2,n
            x(i,1) = x(i,1)-c(i-1,1)*x(i-1,m)
   46     continue
          do 51 i=1,n
            x(i,1) = x(i,1)-d(i,1)*x(i,m)
   51     continue
          do 56 i=1,n-1
            x(i,1) = x(i,1)-e(i+1,1)*x(i+1,m)
   56     continue
        end if
        end if
c
ccc van der worst approximation
c
cdir$ ivdep
        do 35 i=1,n-1
          x(i,j) = x(i,j)-bu(i,j)*x(i+1,j)
   35   continue
cdir$ ivdep
        do 40 i=1,n-2
          x(i,j) = x(i,j)+bu(i+1,j)*bu(i,j)*x(i+2,j)
   40   continue
  200 continue
      return
      end


      subroutine ilucgj(e,d,c,b,a,bu,cu,du,eu,
     1      el,dl,cl,bl,af,bf,cf,df,ef,
     2      n,m,x,y,r,up,aux,istop,ercg,
     3      xratio,ratio,xnorm,ynorm)
c
c***purpose
c
c      given a 9-point matrix multiply routine called matmul9
c      and an approximate factorization a = l*d*u, this routine
c      provides a conjugate gradient acceleration to iterate
c      to an approximate solution to a*x=y.
c
c***description
c
c       ilucgj requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the ldu factors are stored similarly in nine
c       diagonal arrays.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        el,dl,cl,bl,af,bf,cf,df,ef   are 1-d arrays of length
c            at least nxm. the arrays contain from left to right
c            the diagonals of l,d,u as they are swept left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is an estimate to the solution vector of size nxm.
c
c        y     is an r.h.s. vector of dimenison at least nxm.
c
c        r,up,aux  are work variables at least nxm in length.
c               the final residual is left array r.
c
c        istop  is the maximum number of iterations allowed.
c               upon exit it contains the actual number of
c               iterations required to converge to accuracy
c               specified.
c
c        ercg  is the error tolerance for convegence
c
c     on return
c
c        x     is the output vector of dimension at least
c              nxm containing the solution vector.
c
c        xratio   is the ratio of dxnorm to xnorm.
c
c        ratio    is the ratio of rnorm to ynorm
c
c        xnorm    is the l2 norm of the solution vector x.
c
c        ynorm    is the l2 norm of the r.h.s.
c
c-------------------------------------------------------
c***routines called   lu9p matmul9 ldlt rypax ymax9p
c-------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      real*8 e(m*n),d(m*n),c(m*n),b(m*n),a(m*n),bu(m*n),
     1	     cu(m*n),du(m*n),eu(m*n),el(m*n),dl(m*n),cl(m*n),bl(m*n),
     2       af(m*n),bf(m*n),cf(m*n),df(m*n),ef(m*n),
     3       x(m*n),y(m*n),r(m*n),up(m*n),aux(m*n)
      mn = m*n
      e3=ercg
c
ccc copy a into factor storage b

c
      do 10 i=1,mn
        el(i) = e(i)
        dl(i) = d(i)
        cl(i) = c(i)
        bl(i) = b(i)
        af(i) = a(i)
        bf(i) = bu(i)
        cf(i) = cu(i)
        df(i) = du(i)
        ef(i) = eu(i)
   10 continue
c
ccc factor b into incomplete ldu
c
c-------------------------------------------------------
      call lu9p(n,m,el,dl,cl,bl,af,bf,cf,df,ef)
c-------------------------------------------------------

c initialize the residual r
      ynorm = dnrm2(mn,y,1)
c
ccc compute r0
c
c-------------------------------------------------------
      call matmul9(e(1),d(1),c(1),b(1),a(1),
     1             bu(1),cu(1),du(1),eu(1),n,m,x,r)
c-------------------------------------------------------
      do 30 i=1,mn
        r(i) = y(i)-r(i)
   30 continue
      if (dnrm2(mn,r,1) .eq. 0.)     return
c
ccc solve l * l(trans)*aux = r0
c
c-------------------------------------------------------
      call slve9(n,m,el(1),dl(1),cl(1),bl(1),af(1),
     1         bl(2),cl(n),dl(n+1),el(n+2),aux,r)
c-------------------------------------------------------
c
ccc dot r with ldlt[-1]r
c
      rdot = ddot(mn,r,1,aux,1)
c
ccc compute q(0)
c
	nfun=0
c-------------------------------------------------------
      call rypax(eu(-n),du(1-n),cu(2-n),bu(nfun),a(1),
     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,0.d0)
c-----------------------------------------------------
c
ccc compute p0
c

c-------------------------------------------------------
      call slve9(n,m,ef(-n),df(1-n),cf(2-n),bf(nfun),af(1),
     1             bf(1),cf(1),df(1),ef(1),aux,up)
c-------------------------------------------------------

c
ccc begin main loop
c
      do 100 i100 = 1,istop
c
ccc compute alpha
c
      alpha = ddot(mn,up,1,aux,1)
      alpha = rdot/alpha
c
ccc multiply aux by alpha
c
      do 70 i7 = 1,mn
 70   aux(i7) = alpha*aux(i7)
c
ccc compute new x
c
      do 80 i8 = 1,mn
 80   x(i8) = x(i8) + aux(i8)
c
ccc r = r-ax
c
c-------------------------------------------------------
      call ymax9p(e(1),d(1),c(1),b(1),a(1),
     1             bu(1),cu(1),du(1),eu(1),n,m,aux,r)
c-------------------------------------------------------
c
ccc termination test
ccc
ccc stop whenever norm(aux) / norm (x) .lt. e3
c
      rnorm = dnrm2(mn,r,1)
      dxnorm = dnrm2(mn,aux,1)
      xnorm = dnrm2(mn,x,1)
c***
ccc output ratio
c
      ratio = rnorm/ynorm
      xratio=dxnorm/xnorm
c
c***
      if ( xratio .le. e3.and.ratio.le.e3 ) go to 200
c
ccc compute beta
c
c-------------------------------------------------------
      call slve9(n,m,el(1),dl(1),cl(1),bl(1),af(1),
     1         bl(2),cl(n),dl(n+1),el(n+2),aux,r)
c-------------------------------------------------------
      alpha = ddot(mn,r,1,aux,1)
      beta = alpha/rdot
      rdot = alpha
c
ccc compute q
c

c-------------------------------------------------------
      call rypax(eu(-n),du(1-n),cu(2-n),bu(nfun),a(1),
     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,beta)
c-------------------------------------------------------

c
ccc compute p
c

c-------------------------------------------------------
      call slve9(n,m,ef(-n),df(1-n),cf(2-n),bf(nfun),af(1),
     1             bf(1),cf(1),df(1),ef(1),aux,up)
c-------------------------------------------------------

  100 continue

 102  format (1h0,'cg terminates because of too many iterations. ', /,
     x         'xratio = ',e12.4, 'ratio = ',e12.4)
  200 istop=i100
      return
      end



C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C		SUBROUTINES
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


      subroutine lu9p(n,m,e,d,c,b,a,bu,cu,du,eu)
c
c***purpose
c
c      factor (incomplete) a=l*d*u, where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c
c       lu9p requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        n    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c     on return
c
c        e,d,c,b,a,bu,cu,du,eu    contain the factors
c            ldu. the reciprocal of d is stored in array a.
c
c***routines called   none
c
c*** first executable statement   lu9p
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m),
     1        bu(n,m),cu(n,m),du(n,m),eu(n,m)
c
ccc   form a = l*v
c
      do 100 j=1,m
c
ccc factor central block
c
        a(1,j) = 1./a(1,j)
        do 5 i=2,n
          b(i,j) = a(i-1,j)*b(i,j)
          a(i,j) = 1./(a(i,j)-b(i,j)*bu(i-1,j))
    5   continue
        if (j .ne. m)    then
c
ccc   form upper diagonals
c
          do 10 i=2,n
            du(i,j) = du(i,j)-b(i,j)*eu(i-1,j)
   10     continue
          do 15 i=2,n
            cu(i,j) = cu(i,j)-b(i,j)*du(i-1,j)
   15     continue
c
ccc   form lower diagonals
c
          do 20 i=2,n
           e(i,j+1) = a(i-1,j)*e(i,j+1)
   20     continue
          d(1,j+1) = a(1,j)*d(1,j+1)
          do 25 i=2,n
            d(i,j+1) = a(i,j)*(d(i,j+1)-e(i,j+1)*bu(i-1,j))
   25     continue
          do 30 i=1,n-1
            c(i,j+1) = a(i+1,j)*(c(i,j+1)-d(i,j+1)*bu(i,j))
   30     continue
c
ccc   form central diagonals
c
          do 35 i=2,n
            b(i,j+1) = b(i,j+1)-e(i,j+1)*du(i-1,j)-d(i,j+1)*cu(i,j)
   35     continue
          do 40 i=1,n-1
            bu(i,j+1) = bu(i,j+1)-d(i,j+1)*eu(i,j)-c(i,j+1)*du(i+1,j)
   40     continue
          a(1,j+1) = a(1,j+1)-d(1,j+1)*du(1,j)-c(1,j+1)*cu(2,j)
          do 45 i=2,n-1
            a(i,j+1) = a(i,j+1)-e(i,j+1)*eu(i-1,j)
     1                -d(i,j+1)*du(i,j)-c(i,j+1)*cu(i+1,j)
   45     continue
          a(n,j+1) = a(n,j+1)-e(n,j+1)*eu(n-1,j)-d(n,j+1)*du(n,j)
        end if
  100 continue
      do 110 j=1,m
        do 110 i=1,n
          bu(i,j) = a(i,j)*bu(i,j)
  110 continue
c
ccc   form u = d[-1]*v
c
      do 120 j=1,m-1
        do 120 i=1,n
          cu(i,j) = a(i,j)*cu(i,j)
          du(i,j) = a(i,j)*du(i,j)
          eu(i,j) = a(i,j)*eu(i,j)
  120 continue
      return
      end


      subroutine matmul9(e,d,c,b,a,bu,cu,du,eu,nh,nv,x,y)
c
c***purpose
c
c      form the matrix product a*x=y,where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c
c       matmul9 requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nhxnv. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        nh    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        nv    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is the input vector of dimenison at least nhxnv.
c
c     on return
c
c        y     is the output vector of dimension at least
c              nhxnv containing a*x.
c
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      real*8 e(*),d(*),c(*),b(*),a(*),bu(*),cu(*),du(*),eu(*),
     &     x(*),y(*)

c      print*,'b(0)=',b(0),'e(-n)',e(-n),'m,n,=',m,n
c	open(2,file='out1.txt')
c	do k=-3200,3200
c	write(2,*)k,b(k),c(k),d(k),e(k)
c	enddo
c	close(2)
c	stop
c
c*** first executable statement   matmul9
      n = nh*nv
      nhm1 = nh-1
      nhp1 = nh+1
      i = 1
      y(i)=a(i)*x(i)+bu(i)*x(i+1)
     1      +du(i)*x(i+nh)+eu(i)*x(i+nhp1)
      do 10 i=2,nh
      y(i)=b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     1        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1)
     2        )))))
   10 continue
      i = nhp1
      y(i)=d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nhm1)+du(i)*x(i+nh)+eu(i)*x(i+nhp1)
      do 20 i=nh+2,n-nhp1
      y(i)=e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1)
     3        ))))))))
   20 continue
      i=n-nh
      y(i)=e(i)*x(i-nhp1)+d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nh-1)+du(i)*x(i+nh)
      do 30 i=n-nhm1,n-1
      y(i)=e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        )))))
   30 continue
      i = n
      y(i)=e(i)*x(i-nhp1)+d(i)*x(i-nh)
     1      +b(i)*x(i-1)+a(i)*x(i)
c
ccc periodic equations
c
      l = n-nh
      k = l-1
      m = l+1
      i = 1
      y(i) = y(i)+du(i+l)*x(i+l)+eu(i+l)*x(i+m)
      y(i+l) = y(i+l)+du(i+l)*x(i)+cu(i+m)*x(i+1)
      do 40 i=2,nhm1
        y(i) = y(i)+cu(i+l)*x(i+k)+du(i+l)*x(i+l)
     1        +eu(i+l)*x(i+m)
        y(i+l) = y(i+l)+eu(i+k)*x(i-1)+du(i+l)*x(i)
     1        +cu(i+m)*x(i+1)
   40 continue
      i = nh
      y(i) = y(i)+cu(i+l)*x(i+k)+du(i+l)*x(i+l)
      y(i+l) = y(i+l)+eu(i+k)*x(i-1)+du(i+l)*x(i)
      return
      end



       subroutine slve9(n,m,e,d,c,b,a,bu,cu,du,eu,x,y)
c
c***purpose
c
c      solve (incomplete) ax=(l*lt)x=y, where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator. factors are generated lu9p.
c
c***description
c
c       ldlt requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        y     is an input vector of dimenison at least nxm.
c
c     on return
c
c        x     is the output vector of dimension at least
c              nxm containing the solution vector.
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m),
     1          bu(n,m),cu(n,m),du(n,m),eu(n,m),x(n,m),y(n,m)
c
ccc   copy y into x
c
      do 1 j=1,m
        do 1 i=1,n
          x(i,j) = y(i,j)
    1 continue
c
ccc   forward solve
c
      do 100 j=1,m
        if (j .ne. 1)    then
        x(1,j) = x(1,j)-c(1,j)*x(2,j-1)-d(1,j)*x(1,j-1)
cdir$ ivdep
        do 20 i=2,n-1
          x(i,j) = x(i,j)-c(i,j)*x(i+1,j-1)
     1           -d(i,j)*x(i,j-1)-e(i,j)*x(i-1,j-1)
   20   continue
        x(n,j) = x(n,j)-d(n,j)*x(n,j-1)-e(n,j)*x(n-1,j-1)
        end if
c
ccc   van der worst approximation
c
cdir$ ivdep
        do 5 i=n,2,-1
          x(i,j) = x(i,j)-b(i,j)*x(i-1,j)
    5   continue
cdir$ ivdep
        do 10 i=n,3,-1
          x(i,j) = x(i,j)+b(i-1,j)*b(i,j)*x(i-2,j)
   10   continue
  100 continue
c
ccc diagonal solve
c
      do 150 j=1,m
        do 150 i=1,n
          x(i,j) = a(i,j)*x(i,j)
  150 continue
c
ccc  backward solve
c
      do 200 j=m,1,-1
        if (j .ne. m)    then
        x(n,j) = x(n,j)-cu(n,j)*x(n-1,j+1)-du(n,j)*x(n,j+1)
cdir$ ivdep
        do 50 i=2,n-1
          x(i,j) = x(i,j)-cu(i,j)*x(i-1,j+1)
     1           -du(i,j)*x(i,j+1)-eu(i,j)*x(i+1,j+1)
   50   continue
        x(1,j) = x(1,j)-du(1,j)*x(1,j+1)-eu(1,j)*x(2,j+1)
        end if
c
ccc van der worst approximation
c
cdir$ ivdep
        do 35 i=1,n-1
          x(i,j) = x(i,j)-bu(i,j)*x(i+1,j)
   35   continue
cdir$ ivdep
        do 40 i=1,n-2
          x(i,j) = x(i,j)+bu(i+1,j)*bu(i,j)*x(i+2,j)
   40   continue
  200 continue
      return
      end


       subroutine rypax(e,d,c,b,a,bu,cu,du,eu,nh,nv,x,y,r)
c
c***purpose
c
c      form the matrix product y=r*y+a[t]*x,where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c
c       rypax requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nhxnv. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        nh    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        nv    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is an input vector of dimenison at least nhxnv.
c
c        y     is an input vector of dimension at least nhxnv.
c
c        r     is a scalar multiplier of the input vector y.
c
c     on return
c
c        y     is the output vector of dimension at least
c              nhxnv containing r*y+a[t]*x.
c
c***routines called   none
c
      implicit real*8 (a-h,o-z) 
      real*8 e(*),d(*),c(*),b(*),a(*),bu(*),cu(*),du(*),eu(*),
     &       x(*),y(*)
c
c*** first executable statement   rypax
      n = nh*nv   
      nhm1 = nh-1
      nhp1 = nh+1
      i = 1
      y(i) = r*y(i)+(a(i)*x(i)+bu(i)*x(i+1)
     1      +du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 10 i=2,nh
      y(i) = r*y(i)+(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     1        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     2        )))))
   10 continue
      i = nhp1
      y(i) = r*y(i)+(d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nhm1)+du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 20 i=nh+2,n-nhp1
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     3        ))))))))
   20 continue
      i=n-nh
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nh-1)+du(i)*x(i+nh))
      do 30 i=n-nhm1,n-1
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1))
     2        )))))
   30 continue
      i = n	
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+d(i)*x(i-nh)
     1      +b(i)*x(i-1)+a(i)*x(i))
      return
      end



       subroutine ymax9p(e,d,c,b,a,bu,cu,du,eu,nh,nv,x,y)
c
c***purpose
c
c      form the matrix product y=y-a*x,where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c      
c       ymax9p requires that the coefficients of the ith 
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c   
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nhxnv. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c       
c              a is viewed as a block tridiagonal system of equations.
c      
c        nv    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c    
c        x     is an input vector of dimenison at least nhxnv.
c            
c        y     is an input vector of dimension at least nhxnv.
c
c     on return
c              
c        y     is the output vector of dimension at least
c              nhxnv containing y-a*x.
c
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      real*8 e(*),d(*),c(*),b(*),a(*),bu(*),cu(*),du(*),eu(*),
     &       x(*),y(*)
c
c*** first executable statement   ymax9p
      n = nh*nv
      nhm1 = nh-1
      nhp1 = nh+1
      i = 1
      y(i) = y(i)-(a(i)*x(i)+bu(i)*x(i+1)
     1      +du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 10 i=2,nh
      y(i) = y(i)-(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     1        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     2        )))))
   10 continue
      i = nhp1
      y(i) = y(i)-(d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nhm1)+du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 20 i=nh+2,n-nhp1
      y(i) = y(i)-(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     3        ))))))))
   20 continue
      i=n-nh
      y(i) = y(i)-(e(i)*x(i-nhp1)+d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nh-1)+du(i)*x(i+nh))
      do 30 i=n-nhm1,n-1
      y(i) = y(i)-(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1))
     2        )))))
   30 continue
      i = n
      y(i) = y(i)-(e(i)*x(i-nhp1)+d(i)*x(i-nh)
     1      +b(i)*x(i-1)+a(i)*x(i))
      return
      end

 
      subroutine accel
c
c ======================================================================
c
c   Purpose -
c     compute Lagrangian velocities caused by the implicit
c     pressure and explicit viscous forces
c
c   ACCEL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         IMPLCTP
c
c
c   ACCEL calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple
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
      data tiny /1.0d-25/, zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c ---   plus porous - fyshi 10/18/02

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
          rhox=rhobar*(delx(i+1)+delx(i))
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pr(ii,ij).eq.1)then
          u(ij)=u(ij)+delt*(p(ij)-p(ipj))*2.0/rhox
     *    *porosity(ii)/(1+ca(ii))
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          u(ij)=u(ij)+delt*(p(ij)-p(ipj))*2.0/rhox
	  endif

   10     u(ij)=cvmgt(zero,u(ij),(rhorc.lt.frsurf).or.
     &               (ar(ij).lt.em6))
c
c....     reset y-velocity and check for surrounding fluid
          rhobar=rhoy
          if (at(ij).lt.em6) go to 25
c
          rhoy=rhobar*(dely(j+1)+dely(j))
c porous
	  nswi=0
	  ninter=0
	  do ii=1,nportype
	  if(ninter.eq.0.and.pt(ii,ij).eq.1)then
          v(ij)=v(ij)+delt*(p(ij)-p(ijp))*2.0/rhoy
     *    *porosity(ii)/(1+ca(ii))
	  nswi=1
	  ninter=1
	  endif
	  enddo
	  if(nswi.eq.0)then
          v(ij)=v(ij)+delt*(p(ij)-p(ijp))*2.0/rhoy
	  endif
   25     v(ij)=cvmgt(zero,v(ij),(rhotc.lt.frsurf).or.
     &               (at(ij).lt.em6))
c
   20 continue
c
c.... update the boundary conditions
c
c     ---------
      call bc
c     ---------
c
      return
      end
