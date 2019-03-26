
C**********************************************************************
C**                 FUNCTION DNRM2                                   **
C**********************************************************************
      real*8 function dnrm2(n,dx,incx)
      save cutlo, cuthi, zero, one
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      integer next
      real*8 dx(1), cutlo, cuthi, hitest, sum, xmax, zero, one
      data zero, one /0.0d0, 1.0d0/
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
      if(n .gt. 0) go to 10
      dnrm2 = zero
      go to 300
   10 assign 30 to next
      sum = zero
      nn = n * incx
C
C.....BEGIN MAIN LOOP
      i = 1
   20 go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
C
C.....PHASE 1:  SUM IS ZERO
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
C
C.....PREPARE FOR PHASE 2
      assign 70 to next
      go to 105
C
C.....PREPARE FOR PHASE 4
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
C
C.....PHASE 2:  SUM IS SMALL
C     SCALE TO AVOID DESTRUCTIVE UNDERFLOW
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
C
C.....COMMON CODE FOR PHASES 2 AND 4
C     IN PHASE 4 SUM IS LARGE SO SCALE TO AVOID OVERFLOW
  110 if( dabs(dx(i)) .le. xmax ) go to 115
      sum = one + sum * (xmax / dx(i))**2
      xmax = dabs(dx(i))
      go to 200
C
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
C
C.....PREPARE FOR PHASE 3
   75 sum = (sum * xmax) * xmax
C
C.....FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
   85 hitest = cuthi/float( n )
C
C.....PHASE 3:  SUM IS MID-RANGE THEREFORE NO SCALING
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95 sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
C
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
C
C.....END OF MAIN LOOP
C     COMPUTE SQUARE ROOT AND ADJUST FOR SCALING
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
C2345678901234567890123456789012345678901234567890123456789012345678901
C**********************************************************************
C**  Modification history: Created by C.L. Lawson (jpl), R.J. Hanson **
C**                        (snla), D.R. Kincaid (u. of texas), F.T.  **
C**                        Krogh (jpl) OCT79.  Revised DEC86         **
C**                        Adapted from NASA-VOF2D by MCW, JAN90     **
C**                                                                  **
C**    Purpose: Euclidean length (l2 norm) of d.p. vector.  See      **
C**             NASA-VOF2D for further details.                      **
C**                                                                  **
C**    Called from:        ILUCGJ                                    **
C**    Calls to:                                                     **
C**    External fctns:                                               **
C**                                                                  **
C**********************************************************************
      end

       real*8 function cvmgt(x1,x2,x3)
c
c ======================================================================
c
c   purpose -
c     Cray vectorization function:
c     Return x1 if x3 is true, otherwise return x2
c
c   cvmgt is called by -
c
c   cvmgt calls the following subroutines and functions -
c ======================================================================
c
      implicit real*8 (a-h,o-z)
      logical x3
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      if (x3) then
	cvmgt=x1
      else
	cvmgt=x2
      endif
c
      return
      end

       real*8 function dasum(n,dx,incx)                       
c
c***purpose  sum of magnitudes of d.p. vector components              
c***description                                                      
c                                                                   
c                b l a s  subprogram                               
c    description of parameters                                    
c                                                                
c     --input--                                                 
c        n  number of elements in input vector(s)              
c       dx  double precision vector with n elements           
c     incx  storage spacing between elements of dx           
c                                                           
c     --output--                                           
c    dasum  double precision result (zero if n .le. 0)    
c                                                        
c     returns sum of magnitudes of double precision dx. 
c     dasum = sum from 0 to n-1 of dabs(dx(1+i*incx))  
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  dasum                                                 
c                                                                       
      real*8 dx(1)                                         
c***first executable statement  dasum                               
      dasum = 0.d0                                                 
      if(n.le.0)return                                            
      if(incx.eq.1)goto 20                                       
c                                                               
c        code for increments not equal to 1.                   
c                                                             
      ns = n*incx                                            
          do 10 i=1,ns,incx                                 
          dasum = dasum + dabs(dx(i))                      
   10     continue                                        
      return                                             
c                                                       
c        code for increments equal to 1.               
c                                                     
c                                                    
c        clean-up loop so remaining vector length is a multiple of 6.   
c                                                                      
   20 m = mod(n,6)                                                    
      if( m .eq. 0 ) go to 40                                          
      do 30 i = 1,m                                                   
         dasum = dasum + dabs(dx(i))                                 
   30 continue                                                      
      if( n .lt. 6 ) return                                        
   40 mp1 = m + 1                                                 
      do 50 i = mp1,n,6                                          
         dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))  
     1   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))             
   50 continue                                                      
      return                                                       
      end                                        


C**********************************************************************
C**                  FUNCTION DDOT                                   **
C**********************************************************************
      real*8 function ddot(n,dx,incx,dy,incy)
      real*8 dx(1),dy(1)
      ddot = 0.d0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
C
C.....CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      return
C
C.....CODE FOR BOTH INCREMENTS EQUAL TO 1
C     CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         ddot = ddot + dx(i)*dy(i)
   30 continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
         ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     1   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      return
C
C.....CODE FOR POSITIVE EQUAL INCREMENTS .NE.1
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          ddot = ddot + dx(i)*dy(i)
   70     continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      return
C2345678901234567890123456789012345678901234567890123456789012345678901
C**********************************************************************
C**  Modification history: Created by C.L. Lawson (jpl), R.J. Hanson **
C**                        (snla), D.R. Kinkaid (u. of texas), F.T.  **
C**                        Krogh (jpl) OCT79. Revised DEC86          **
C**                        Adapted from NASA-VOF2D by MCW, JAN90     **
C**                                                                  **
C**   Purpose:  d.p. inner product of d.p. vectors. See NASA-VOF2D   **
C**             for further details.                                 **
C**                                                                  **
C**   Called from:         CGITJ,ILUCGJ                              **
C**   Calls to:                                                      **
C**   External fctns:                                                **
C**                                                                  **
C**********************************************************************
      end
  

 
      function fxy(x,y,a1,a2,b1,b2,c1,c2,d1,
     &                              kxc,d2,kxs,e1,kyc,e2,kys)
c
c ======================================================================
c
c   Purpose -
c     solve the function f(x,y) for a point (x,y), where:
c     f(x,y) = a1*x + a2*x**2 + b1*y + b2*y**2 + c1 + c2*x*y
c            +d1*cos(kxc*x)+d2*sin(kxs*x)+e1*cos(kyc*y)+e2*sin(kys*y)
c
c   FXY is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET  INITVOFF
c
c
c   FXY calls the following subroutines and functions -
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
      real*8 kxc,kxs,kyc,kys
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      fxy=a2*x*x+a1*x+b2*y*y+b1*y+c2*x*y+c1
     &      +d1*cos(kxc*x)+d2*sin(kxs*x)
     &      +e1*cos(kyc*y)+e2*sin(kys*y)
c
      return
      end

 
      subroutine mollify(i1,i2,j1,j2,nqi,nqj,nsmooth,raw,smooth,tmp)
c
c ======================================================================
c
c   Purpose -
c     Mollify (smooth) cell-based data
c
c   MOLLIFY is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         TENSION
c
c
c   MOLLIFY calls the following subroutines and functions -
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
      dimension raw(1),smooth(1),tmp(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 10 j=j1-1,j2+1
        do 10 i=i1-1,i2+1
          ij=(j-1)*nqj + (i-1)*nqi + 1
          tmp(ij)=raw(ij)
   10 continue
c
      do 1 ns=1,nsmooth
c
        do 20 j=j1,j2
          do 20 i=i1,i2
            ij=(j-1)*nqj + (i-1)*nqi + 1
            ipj=ij+nqi
            ijp=ij+nqj
            ijm=ij-nqj
            imj=ij-nqi
            imjm=ij-nqi-nqj
            imjp=ij-nqi+nqj
            ipjm=ij-nqj+nqi
            ipjp=ij+nqj+nqi
c
c....       Quadratic
            smooth(ij)=9.*tmp(ij)/16. +
     &               3.*(tmp(ipj)+tmp(imj)+tmp(ijp)+tmp(ijm))/32. +
     &               1.*(tmp(ipjp)+tmp(imjp)+tmp(ipjm)+tmp(imjm))/64.
c
   20   continue
c
        ijb=2
        ijt=nqj*j2+2
        do 30 i=i1,i2
          smooth(ijb)=smooth(ijb+nqj)
          smooth(ijt)=smooth(ijt-nqj)
          ijb=ijb+1
          ijt=ijt+1
   30   continue
        ijr=nqj
        ijl=1
        do 40 j=j1-1,j2+1
          smooth(ijr)=smooth(ijr-nqi)
          smooth(ijl)=smooth(ijl+nqi)
          ijr=ijr+nqj
          ijl=ijl+nqj
   40   continue
c
        do 50 j=j1-1,j2+1
          do 50 i=i1-1,i2+1
            ij=(j-1)*nqj + (i-1)*nqi + 1
            tmp(ij)=smooth(ij)
   50   continue
c
    1 continue
c
      return
      end



C******************************************************
C        Function ismax
C  The function should be in cftmath lib. But we could not find it
C  in gnu.
C  Description:
C  Index of the real vector element with maximum value 
C  Fengyan 08/14/02
C******************************************************
	integer function ismax(n,x,incx)
	integer n,incx
	real*8  x(1000000)

	rmax=-3.E+38

	if ((1+(n-1)*incx).gt.1000000) then
	print*, 'n exceed 1000000 in ismax function, stop'
	stop
	else
	do k=1,n,incx
	  if(x(k).gt.rmax)then
	    rmax=x(k)
            ismax=k
	  endif
	enddo

	endif

	end

C******************************************************
C        Function ismin
C  The function should be in cftmath lib. But we could not find it
C  in gnu.
C  Description:
C  Index of the real vector element with minimum value 
C  Fengyan 08/14/02
C******************************************************
	integer function ismin(n,x,incx)
	integer n,incx
	real*8  x(1000000)

	rmin=3.E+38

	if ((1+(n-1)*incx).gt.1000000) then
	print*, 'n exceed 1000000 in ismin function, stop'
	stop
	else
	do k=1,n,incx
	  if(x(k).lt.rmin)then
	    rmin=x(k)
            ismin=k
	  endif
	enddo

	endif

	end


 
      subroutine crvature(i1,i2,j1,j2,nqi,nqj,ctilde,gradcx,gradcy,
     &                    delx,dely,r,ri,kappa,tiny)
c
c ======================================================================
c
c   Purpose -
c     Compute the curvature kappa(ij) from the color function ctilde(ij)
c
c   CRVATURE is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         TENSION
c
c
c   CRVATURE calls the following subroutines and functions -
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
      dimension ctilde(1),gradcx(1),gradcy(1),delx(1),
     &          dely(1),r(1),ri(1),kappa(1)
      real*8 kappa
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 110 j=j1,j2
        do 100 i=i1,i2
          ij=(j-1)*nqj+(i-1)*nqi+1
          ipj=ij+nqi
          imj=ij-nqi
          ipjp=ij+nqj+nqi
          ijp=ij+nqj
          ijm=ij-nqj
          im=i-1
          ip=i+1
          jm=j-1
          jp=j+1
c
          avnx=0.25*(gradcx(ipj)+gradcx(ipjp)
     &            +gradcx(ijp)+gradcx(ij))
          avny=0.25*(gradcy(ipj)+gradcy(ipjp)
     &            +gradcy(ijp)+gradcy(ij))
          alphc=sqrt(avnx**2 + avny**2) + tiny
          ralph=1./(alphc+tiny)
c
          dnxdx=(gradcx(ipj)+gradcx(ipjp)-gradcx(ijp)-gradcx(ij))
     &            /(2.*delx(i))
          dnxdy=(gradcx(ipjp)+gradcx(ijp)-gradcx(ipj)-gradcx(ij))
     &            /(2.*dely(j))
          dnydx=(gradcy(ipj)+gradcy(ipjp)-gradcy(ijp)-gradcy(ij))
     &            /(2.*delx(i))
          dnydy=(gradcy(ipjp)+gradcy(ijp)-gradcy(ipj)-gradcy(ij))
     &            /(2.*dely(j))
          gradxa=avnx*ralph*(avnx*ralph*dnxdx+avny*ralph*dnydx)
          gradya=avny*ralph*(avnx*ralph*dnxdy+avny*ralph*dnydy)
          divn=(r(i)*(gradcx(ipj)+gradcx(ipjp))
     &          -r(i-1)*(gradcx(ijp)+gradcx(ij)))/(2.*ri(i)*delx(i))
     &         +(gradcy(ijp)+gradcy(ipjp)-gradcy(ipj)-gradcy(ij))
     &           /(2.*dely(j))
          divnht=-ralph*(divn-(gradxa+gradya))
c
          kappa(ij)=divnht
c
  100   continue
  110 continue
c
      return
      end
